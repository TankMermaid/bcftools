# bcftools SNP Calling 全流程梳理

> 对象:当前仓库 `bcftools`(git HEAD `905ea784`)核心代码
> 范围:`bcftools mpileup` + `bcftools call` 构成的经典 SNP/indel calling 流水线
> 本文行号引用基于当前工作区源码,供代码定位用

---

## 1. 总览:SNP calling 的两段式流水线

bcftools 的变异检测不是单进程一步完成,而是经典的两段式:

```
                    ┌──────────────────────────────────────────────────────────┐
                    │ 阶段一: bcftools mpileup (mpileup.c 为主控)              │
  BAM/CRAM ───────► │  reads 过滤 → pileup → BAQ 重比对 → 位点似然(PL)          │
  FASTA(参考)       │  → 等位基因统计/偏差检测 → 输出 BCF(含 PL/QS/I16/AD...)  │
                    └───────────────────────┬──────────────────────────────────┘
                                             ▼
                    ┌──────────────────────────────────────────────────────────┐
                    │ 阶段二: bcftools call (vcfcall.c 为主控)                  │
  BCF ────────────► │  读入 PL → 基因型/等位基因频率统计推断 → GT/GQ/QUAL 调用  │
                    │  (-c consensus / -m multiallelic / -C 约束 / gVCF)       │
                    └──────────────────────────────────────────────────────────┘
```

设计要点:**阶段一只做"测序数据 → 基因型似然(PL)"的纯统计,不做任何生物学决策**;所有先验、频率估计、变异/基因型判定都推迟到阶段二。因此同一份 BCF 可以被 `-c`、`-m`、`-C`、`-G`、`-F` 等不同策略重放,而不必重读 BAM。

---

## 2. 模块划分(8 个模块)

| 模块 | 文件 | 职责 | 输入 → 输出 |
|---|---|---|---|
| M1 reads 预处理与 pileup | `mpileup.c` | 逐 read 过滤、覆盖度堆叠、按样本分组 | BAM → 每个位点的 `bam_pileup1_t` 列 |
| M2 BAQ 重比对 | `mpileup.c`(`mplp_realn`) | 局部重比对,压低 indel 邻域的系统性碱基质量错误 | 重叠 indel 的 reads → 更新后 BQ |
| M3 位点/基因型似然与等位基因统计 | `bam2bcf.c` + `htslib/errmod` | 每条 read 的误差模型、PL 计算、等位基因排序、覆盖度/偏差注释 | pileup 列 → `bcf_callret1_t` → `bcf_call_t` → BCF 记录 |
| M4 indel 候选检测与打分 | `bam2bcf_indel.c` / `bam2bcf_edlib.c` / `bam2bcf_iaux.c` | 三类 indel 实现:经典 consensus 比对、edlib 长读、indels-2.0 | 含 indel 的列 → indel 型 PL/AD/IDV/IMF |
| M5 call 主控与模式调度 | `vcfcall.c` | 参数解析、ploidy、过滤、分派 ccall/mcall/qcall、输出 | BCF → 调用结果 |
| M6 consensus caller(-c) | `ccall.c` + `prob1.c` + `em.c` | EM 频率估计 + AFS(等位基因频率谱)贝叶斯推断 + 基因型 MAP | BCF → GT/GQ/QUAL + AF1/G3/HWE/FQ 等 |
| M7 multiallelic caller(-m) | `mcall.c` | 质量加权等位基因频率先验、等位基因子集枚举、基因型 MAP | BCF → GT/GQ/GP + AC/AN + QUAL |
| M8 gVCF 输出 | `gvcf.c` | 连续纯参考位点按 DP 档折叠成 block | ref 位点流 → END/MIN_DP 块记录 |

依赖关系:`M1→M2→M3`(阶段一,每个位点一次);`M4` 在 M3 之外对同一列做 indel 通道;`M5` 读入 BCF 后按模式驱动 `M6`/`M7`;`M8` 挂在阶段一(`mpileup -g`)或阶段二(`call --gvcf`)的输出端。

---

## 3. 模块详解

### M1 — reads 预处理与 pileup(`mpileup.c`)

**技术路线**(按每个位点的处理链):

1. **逐 read 过滤**(`mplp_func`,L195):跳过未比对、不满足 FLAG 过滤(`-f/-F`)、不在 `-r/-R` region、样本不在 `-s/-S` 列表的 read;`--illumina-1.3+` 时把 Q>31 的碱基质量截断到 31(L230);按 `-C` 阈值用 `sam_cap_mapq` 做 mapQ 校正(L288);剔除 `min_mq` 以下的 read。
2. **pileup 构建**(`pileup_constructor`,L305):htslib `bam_mplp_auto` 同步推进多文件、按坐标堆叠 reads;每个 read 缓存一个 `plp_cd_t`(bam2bcf.h L101),记录 **样本号、是否含 soft-clip、是否含 indel、是否已重比对、NM 值**,避免在每个位点重复计算。
3. **按样本分组**(`group_smpl`,L361):同一 BAM 中的不同样本(或 `--samples-file` 合并的样本)被拆到独立 pileup 列,之后每个样本单独算 GL。
4. **重叠配对处理**:默认开启 `MPLP_SMART_OVERLAPS`(L939,`bam_mplp_init_overlaps`),双端 read 重叠区只取质量更高的一端,避免 PCR 重复计数。

**默认参数**(`main_mpileup`,L1373 起):`min_baseQ=1`、`max_baseQ=60`、`capQ_thres=0`、`max_depth=250`、`max_indel_depth=250`、`max_read_len=500`、`indel_win_size=110`、flag 默认 `MPLP_NO_ORPHAN|MPLP_REALN|MPLP_REALN_PARTIAL|MPLP_SMART_OVERLAPS`。

---

### M2 — BAQ 重比对(`mplp_realn`,L420)

**动机**:indel 附近比对错误是 SNP 假阳性的主要来源;BAQ(Base Alignment Quality)通过局部重比对把"可能是错配"的碱基质量压低。

**技术路线**(启发式 + 质量重估):

1. **只在"值得重比对"的位点触发**:扫描该列,若存在 indel(`p->indel` 或 CIGAR 含 I/D/N)或有 soft-clip 的 read 比例过高,才进入重比对(L429-443);`-D` 关闭 partial 模式时条件更严格(`has_indel==0` 或"单一种类 indel 且占比 <10%"直接跳过,L445-449)。
2. **对每个 read 做成本裁剪**:已被重比对过的 read 跳过(L468);过长 read(`> max_read_len`)跳过;对短读,若 indel 位于 read 中部(`lm/rm` 两侧匹配段都足够长,阈值 `REALN_DIST = 40+10*(nt<40)+10*(nt<20)`,L487)则跳过——深覆盖下容忍"偷懒",浅覆盖下从严;长读(>500bp)计算比对带宽成本,过大则跳过(L537)。
3. **质量重估**:调用 htslib `sam_prob_realn(b, ref, ref_len, k)`(L548)对 read 重新比对并输出新碱基质量(`-E` 用更激进的 `k=7`,默认 `k=3`);结果写回 BAM 记录,后续 `bcf_call_glfgen` 直接使用新质量。

注意:`mplp_realn` 会修改 `bam1_t`(追加 ZQ tag),因此必须在 pileup 建好哈希之前预留空间(mpileup.c L248-286 的"fudge"技巧)。

---

### M3 — 位点/基因型似然与等位基因统计(`bam2bcf.c`)

这是阶段一的核心统计模块,三个核心函数形成严格流水线:

```
bcf_call_glfgen(单样本列) ──► bcf_callret1_t(每样本)
        │ 每条 read → 25 项基因型似然 p[25] + 注释 anno[16]
        ▼
bcf_call_combine(跨样本)  ──► bcf_call_t(全位点)
        │ 等位基因排序 + PL 数组 + 偏差统计
        ▼
bcf_call2bcf ──► bcf1_t(BCF 记录)
```

#### 3.1 单样本:位点似然生成(`bcf_call_glfgen`,L250)

输入:该样本该位点的 pileup 列 + 参考碱基(4-bit)。输出:`bcf_callret1_t`,其中 `p[25]` 是 5 个等位基因(ACGTN/indel 类型)两两组合的 25 项 phred 基因型似然。

**每条 read 的质量链(SNP 通道,L422-463)**:
1. 取碱基 `b = seq_nt16_int[...]`;比对 N 时按参考碱基补齐(L424-425);
2. `baseQ = min(自身Q, 左邻Q+delta_baseQ, 右邻Q+delta_baseQ)`——**相邻碱基质量差阈值**,压制 homopolymer 处的质量虚高(L430-435);
3. 过滤 `min_baseQ`,截断 `max_baseQ`;
4. 最终进入似然计算的质量 `q = min(baseQ, seqQ=99, mapQ_capped)`(L459-461):mapQ 先按 `capQ` 截断,`mapQ=0` 的 read 强制按 BQ=4 处理(L452-463),再统一 `q∈[4,63]`;
5. 编码进 `bca->bases[]`:`q<<5 | is_rev<<4 | b`(L464),之后一次性交给误差模型。

**误差模型 `errmod_cal`(L570)**:`bcf_callaux_t->e = errmod_init(1-theta)`(L55)由 htslib 提供。它按每条 read 的错误概率 `ε=10^(-Q/10)` 建立"真实等位基因 → 观测碱基"的 5×5 转移矩阵,对所有 read 求积,得到每个基因型(等位基因对)的似然;`theta`(默认 0.83)调节同源/杂合基因型的先验权重。`p[25]` 中索引为 `g = a*5 + b`(L1031 附近)。

**注释统计(L480-489)**:`anno[16]` 是 16 项累加量,布局(bam2bcf.h L150-162):
```
0-3   深度(ref-fwd, ref-rev, alt-fwd, alt-rev)
4-7   baseQ、baseQ²  (ref / alt)
8-11  mapQ、mapQ²    (ref / alt)
12-15 到 read 端的 minDist、minDist² (ref / alt)
```
另有 `QS`(质量累加)、`QM`(平均质量)、`ADF/ADR`、MQ0 计数、`SCR`(soft-clip read 数),以及供偏差检验的**直方图**:碱基位置 `ref_pos/alt_pos`(按 read 内相对位置缩放到 100 桶)、`ref_mq/alt_mq`、`ref_bq/alt_bq`、`fwd_mqs/rev_mqs`、soft-clip 长度 `ref_scl/alt_scl`、错配数 `ref_nm/alt_nm`(L491-537)。

indel 通道(L312-421):质量来自 `p->aux`(mpileup 在 indel 打分时写入的 indelQ),并在高深度/低支持时做启发式降权(如 `seqQ = min(seqQ, seqQ_offset-5*min(20,_n))`,L397)。

#### 3.2 跨样本:合并与等位基因排序(`bcf_call_combine`,L959)

1. **质量分数 QS**:每个样本内把各等位基因的 `QS` 按样本覆盖度归一化后跨样本求和(L968-976),得到每个等位基因的全群体质量分;
2. **等位基因排序**:参考碱基固定为 `a[0]`,其余按 QS 降序进入 `a[1..]`;QS=0 的等位基因被剔除;对 SNP,若存在"本次未观测"的碱基,则追加一个占位等位基因 `a[unseen] = <*>`(L990-1007)——这样 PL 数组在变异/非变异位点间维度一致,`call` 阶段可以据此知道"这个 ALT 其实没看到";
3. **PL 数组生成(L1024-1049)**:对每个样本,从 25 项 `p[]` 里按等位基因排序抽出 `n_alleles*(n_alleles+1)/2` 项基因型似然,**减去该样本最小似然再取整(cap 255)**——即 VCF 的 phred-scaled PL;
4. **注释对齐**:`DP4`、`ADF/ADR`、`QS`、`QM` 都按新等位基因顺序重排(L1050-1134);
5. **偏差统计(仅当存在真 ALT,L1150-1197)**:
   - `FS`:正反链 2×2 表的 **Fisher exact test**(`kt_fisher_exact`);
   - `SGB`:Segregation-based bias,按泊松模型比较"变异真实/假阳性"两种情形下各样本变异 read 数分布的对数似然(`calc_SegBias`,L895);
   - `RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ/NMBZ`:ref 与 alt 两组在位置、mapQ、baseQ、链、soft-clip、NM 直方图上的 **Mann-Whitney U 检验 Z 分数**(`calc_mwu_biasZ`,L817;经典 `mann_whitney_1947` 在 L695);
   - `VDB`:变异 read 在 read 内位置的分布随机性(`calc_vdb`,L600,对 100bp read 拟合参数)。

#### 3.3 BCF 记录生成(`bcf_call2bcf`,L1200)

组装 INFO(`DP/ADF/ADR/AD/DPR/SCR/I16/QS/DP4`(旧格式)/`VDB/SGB/NM/RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ/NMBZ/FS/MQ0F`)与 FORMAT(`PL/DP/DV/SP/AD/ADF/ADR/SCR/QS`),indel 额外写 `INDEL/IDV/IMF`。**QUAL 在此阶段恒为 0**,留给 call 阶段计算。

---

### M4 — indel 候选检测与打分(3 个实现)

`mpileup_reg` 的 indel 通道(L591-613):同一 pileup 列,先跑 SNP 通道,再在 `total_depth < max_indel_depth` 时对含 indel 的位点做第二遍处理。三个实现由 `--indels-2.0` / `--indels-cns`(edlib)选择:

| 实现 | 入口 | 技术路线 |
|---|---|---|
| 经典 v1 | `bcf_call_gap_prep`(bam2bcf_indel.c L698) | ① 找 indel 类型集 `bcf_cgp_find_types`;② 以参考为模板做"多数一致性"读段 `bcf_cgp_ref_sample`,把 ≥70% read 支持的错配掩成 N(消除 SNP 干扰,L736-744);③ 对每种 indel 类型构造参考片段与 read 片段,用 `bcf_cgp_align_score`(banded 动态规划,`openQ=40/extQ=20`)打分;④ `bcf_cgp_compute_indelQ` 把"该 read 支持 indel i vs 支持参考"的似然比转成 indelQ,写入 `p->aux`;⑤ 同源重复区(est_indelreg)对 indelQ 降权 |
| edlib 长读 | `bcf_edlib_gap_prep`(bam2bcf_edlib.c L1319) | 同上思路,但用 edlib 做全局/局部比对(`bcf_edlib_realign`),插入一致性序列由 `bcf_cgp_calc_ins_cons` 按频率构造;对 >500bp 长读、CCS 数据有专门带宽/质量裁剪 |
| indels-2.0 | `bcf_iaux_gap_prep`(bam2bcf_iaux.c L701) | 按样本先构造一致性序列 `iaux_set_consensus`,再对每条 read 用序列上下文模型打分 `iaux_align_read`/`iaux_score_reads`,输出更细的 indel 类型与质量 |

三类实现都通过 `p->aux` 把 indelQ 传回 `bcf_call_glfgen`(M3 的 indel 通道),与 SNP 共用同一套 PL/合并/注释机制。

---

### M5 — call 主控(`vcfcall.c`)

**入口 `main_vcfcall`(L1003)**:解析参数(L1025-1161),关键默认值:
`theta=1.1e-3`(`-P` 先验)、`pref=0.5`(`-p` 变异判定阈值)、`min_perm_p=0.01`、`min_lrt=1`、trio 突变率 `Pm_SNPs=1-1e-8`、`Pm_ins=Pm_del=1-1e-9`。

**主循环(L1198-1261)**:
1. 过滤:`-V snps/indels`、`REF` 首碱基为 N 的位点(默认 `-N` 行为);
2. 识别"unseen"等位基因(`<*>` 或旧 `X`,L1210-1219),判定 `is_ref`;
3. `bcf_unpack` + 按 sex/contig 设 ploidy(`set_ploidy`);
4. 分派:`-m` → `mcall`(M7),否则 `-c` → `ccall`(M6);`-C alleles`/`-C trio` 通过 `call->flag` 进入对应实现;
5. 输出:普通 VCF、`-v` 只出变异位点、`--gvcf` 走 M8。

---

### M6 — consensus caller(`-c`):`ccall.c` + `prob1.c` + `em.c`

这是经典贝叶斯调用模型(源自 Li 2011, `doi:10.1093/bioinformatics/btr509` 的框架)。

**流水线**(`ccall`,L314-338):

```
PL ──► set_pdg3: PL→P(D|g)(仅前 2 个等位基因,3 种基因型)
    ──► bcf_em1: EM 估计等位基因频率 + 基因型频率 + HWE/LRT
    ──► bcf_p1_cal: AFS 贝叶斯推断(后验频率谱、p_ref/p_var、AC、CI)
    ──► bcf_p1_call_gt: 逐样本 HWE 先验 + MAP 基因型
    ──► update_bcf1: QUAL/FQ/PV4/AC1/DP4/MQ + 等位基因裁剪
```

**统计架构**:

1. **P(D|g) 提取**(`set_pdg3`,L90):`pdg = [P(D|hom-alt), P(D|het), P(D|hom-ref)]`,注意索引按 VCF 顺序反转(ccall 内部用"先 hom-alt"布局)。
2. **EM 频率估计**(`bcf_em1`,em.c L167):
   - 初值 `est_freq`:按每个样本的最大似然基因型粗计数(L44);
   - `freqml`:EM 迭代(`freq_iter`,L87)最多 10 次,不收敛则切换 **Brent 一维最小化**(`kmin_brent`,L109-121)求 MLE 频率;
   - `g3_iter`(L124)迭代 3 类基因型频率(偏离 HWE),用于 **HWE 检验**:LRT 统计量 `2·log(P(D|ĝ)/P(D|HWE))` 经 `kf_gammaq` 转 p 值(L192-197);
   - 若 `-1` 指定了样本分组,输出两组频率 `x[5],x[6]` 与 1/2 自由度 LRT p 值(L198-221)。
3. **AFS 贝叶斯推断**(`bcf_p1_cal`,prob1.c L453):
   - **核心 DP**(`mc_cal_y_core`,L212):令 `z[k] = P(观测数据 | 群体中恰有 k 条 ALT 染色体)`,逐样本用多项式系数动态规划累积(二倍体转移系数 `(M0-k+1)(M0-k+2)、k(M0-k+2)、k(k-1)`,对应 0/0、0/1、1/1,L237-240),数值上每步归一化并取对数(L241-243),支持混合倍性;
   - **后验 AFS**(`mc_cal_afs`,L423):`afs1[k] = phi[k]·z[k]/Σ`,其中 `phi` 是先验(默认 `MC_PTYPE_FULL`:`phi[i]=theta/(M-i)`,即每个额外 ALT 等位基因乘一次 θ,L47-62;indel 用 `phi_indel = phi·0.15`,L39-45);
   - 输出:`f_exp`(期望频率)、`p_ref`/`p_var`(折叠/非折叠的变异后验概率)、`ac`(最大后验 ALT 染色体数)、95% 等尾可信区间 `cil/cih`(L486-498);
   - `do_contrast` 时对两组样本做超几何加权的频率对比检验 `contrast2`(L345,`-U` 旧功能)。
4. **变异判定与 QUAL**(`update_bcf1`,L139-233):`is_var = p_ref < pref`(`-p`,默认 0.5);`QUAL = -4.343·log(min(p_ref,p_var))`(L230);`FQ` 带符号的折叠后验评分(L202);`PV4` 由 `test16` 给出 4 个偏差 p 值——正反链 Fisher + baseQ/mapQ/tail-distance 三个 t 检验(L115-130)。
5. **基因型调用**(`bcf_p1_call_gt`,prob1.c L181):对每个样本 `P(g) ∝ P(D|g)·HWE(g|f_exp)`,取 MAP;非变异位点强制 0/0(L202);GQ = phred(1-max posterior)。

---

### M7 — multiallelic caller(`-m`):`mcall.c`

目标:一次考虑所有等位基因(含多等位基因),适用于稀有变异/群体数据。

**统计架构**:

1. **频率先验 = 质量加权 QS**(`mcall`,L1422-1527):从 INFO/QS 读每等位基因质量分并归一化为频率 `f_x`(L1445,注释明确写 "f_x = QS/N in Eq. 1 in call-m math notes");`-G` 分组模式改为用 FORMAT/AD 或 QS 按组聚合(L1458-1496);`-F AN,AC` 参考 panel 先验时把 panel AF 与样本 QS 加权混合(L1499-1519)。
2. **等位基因子集选择**(`mcall_find_best_alleles`,L594-713):枚举 1、2、3 个等位基因的所有组合,对每组计算全体样本的对数似然:
   - 纯合 `i/i`:`Σ_s log P(D_s|i/i)`;杂合组合 `i/j` 用频率加权混合 `f_i²·P(i/i)+f_j²·P(j/j)+2f_i f_j·P(i/j)`(L632-651);三等位同理(L656-701);
   - **每个参与的非参考等位基因加一次先验惩罚 `theta`**(L616,`call->theta` 在 `mcall_init` 里乘了 Watterson 因子 `aM=Σ1/i` 后取 log,L397-410);
   - 取似然最大的组合 `als_new`,并保证 REF 在场(L1556);
3. **基因型调用**(`mcall_call_genotypes`,L748-889):对每个样本,`GP(g) ∝ P(D|g)·(f_i f_j)`,取 MAP 得 GT,`GQ = -4.343·log(1-max/sum)`(L880);计算 AC/AN(L842-843,1644-1647);
4. **QUAL**:变异位点 `QUAL = -4.343·(ref_lk - logsumexp(lk_sum, ref_lk))`,即"仅 REF 模型 vs 最优组合模型"的似然比(L1546,1631);纯参考位点用 `lk_sum` 侧(L1636-1641);
5. **后处理**:裁剪未入选等位基因(PL 重排 `mcall_trim_and_update_PLs`,L1161;Number=R 的 INFO 重排 `mcall_trim_and_update_numberR`)、`-A` 保留全部 ALT、`-C alleles` 强制给定等位基因集(`mcall_constrain_alleles`,L1274)、`-C trio` 孟德尔约束(代码中已禁用,L1608,公式见 L892-910 注释)。

---

### M8 — gVCF 输出(`gvcf.c`)

`gvcf_write`(L88):把连续的纯参考位点(`is_ref`)按 `--gvcf DP 档位`聚合成块:

- 每位的 per-sample DP 最小值决定档位 `dp_range`(`-g 2,5,10,...` 的档位表,L114-117);
- 块内记录 `END`、`MIN_DP`、聚合后的 `QS/PL/DP`(各样本取块内最小值,L192-194);
- 遇到变异位点、档位变化、染色体切换或坐标不连续时先 flush 前块(L130-166)。

---

## 4. 关键统计架构速查表

| 统计量 | 含义 | 计算位置 | 依赖模型 |
|---|---|---|---|
| PL | 基因型 phred 似然 | `bcf_call_glfgen` L570 + `bcf_call_combine` L1024 | htslib `errmod`(5×5 误差矩阵) |
| QS | 等位基因质量分(归一化) | `bcf_call_combine` L968 | 碱基质量累加 |
| I16 | 16 项深度/质量/距离累加量 | `bcf_call_glfgen` L480 | — |
| DP4/AD/ADF/ADR | 正反链×ref/alt 深度 | `bcf_call_combine` L1050 | — |
| RPBZ/MQBZ/BQBZ/MQSBZ/SCBZ/NMBZ | 偏差 Z 分数 | `calc_mwu_biasZ` bam2bcf.c L817 | Mann-Whitney U |
| FS | 链偏差 p 值 | `bcf_call_combine` L1156 | Fisher exact |
| VDB | 变异位点分布偏差 | `calc_vdb` L600 | 拟合参数 |
| SGB | 分离偏差对数似然 | `calc_SegBias` L895 | 泊松混合 |
| AF1/G3/HWE/AF2/LRT | 频率与 HWE/分组检验 | `bcf_em1` em.c L167 | EM + Brent |
| AFS/AC/CI/FQ | 后验频率谱/等位基因计数/可信区间 | `bcf_p1_cal` prob1.c L453 | DP 累积 `mc_cal_y_core` |
| QUAL(ccall) | -4.343·log(变异后验) | `update_bcf1` L230 | p_ref |
| QUAL(mcall) | 频率模型似然比 | `mcall` L1546 | logsumexp |
| GQ | phred(1-MAP) | `bcf_p1_call_gt` L205 / `mcall_call_genotypes` L880 | MAP |
| indelQ | indel 支持度质量 | `bcf_cgp_compute_indelQ`(bam2bcf_indel.c L597) | banded DP 比对 |

---

## 5. 主调用链(文件粒度)

```
main.c ──► main_mpileup (mpileup.c L1373)
             └─ mpileup() (L618) / mpileup_reg() (L555)
                  ├─ mplp_realn()                  [M2]
                  ├─ group_smpl()
                  ├─ bcf_call_glfgen()  ──► errmod_cal()      [M3]
                  ├─ bcf_call_combine()                          [M3]
                  ├─ bcf_call2bcf()                              [M3]
                  └─ bcf_call_gap_prep() / bcf_edlib_gap_prep()
                     / bcf_iaux_gap_prep()                       [M4]

main.c ──► main_vcfcall (vcfcall.c L1003)
             └─ 主循环 (L1198)
                  ├─ ccall() (ccall.c L314)
                  │    ├─ bcf_em1() (em.c L167)                  [M6]
                  │    ├─ bcf_p1_cal() (prob1.c L453)            [M6]
                  │    └─ update_bcf1() (ccall.c L139)           [M6]
                  ├─ mcall() (mcall.c L1422)                     [M7]
                  │    ├─ mcall_find_best_alleles() (L594)
                  │    └─ mcall_call_genotypes() (L748)
                  └─ gvcf_write() (gvcf.c L88)                   [M8]
```

---

## 6. 相关文件清单

| 文件 | 角色 |
|---|---|
| `mpileup.c` | 阶段一主控、reads 过滤、pileup、BAQ 重比对 |
| `bam2bcf.c` / `bam2bcf.h` | 位点似然、等位基因统计、偏差检验、BCF 组装;核心数据结构 `bcf_callaux_t`/`bcf_callret1_t`/`bcf_call_t` |
| `bam2bcf_indel.c` | 经典 indel 打分(consensus + banded DP) |
| `bam2bcf_edlib.c` | edlib 长读 indel 打分 |
| `bam2bcf_iaux.c` | indels-2.0 indel 打分 |
| `vcfcall.c` | 阶段二主控、参数、ploidy、模式调度 |
| `call.h` | call 阶段共享数据结构 `call_t`、`smpl_grp_t`、`family_t` |
| `ccall.c` | consensus caller 主流程与输出 |
| `prob1.c` / `prob1.h` | AFS 贝叶斯推断、基因型调用、先验 |
| `em.c` | EM 频率/基因型频率估计、LRT 检验 |
| `mcall.c` | multiallelic caller 主流程 |
| `gvcf.c` / `gvcf.h` | gVCF 块折叠 |
| `ploidy.c` / `ploidy.h` | 倍性解析 |
| `htslib/errmod.c`(外部依赖) | 5×5 测序误差模型,`errmod_cal` 生成基因型似然 |
| `HMM.c`(不相关) | 仅用于 `color-chrs` 插件(CNV),与 SNP calling 无关 |

---

## 7. 备注

- 阶段一 QUAL=0、不做变异判定,保证"统计与决策解耦";阶段二可在同一 BCF 上重放不同先验/策略。
- `-c` 与 `-m` 的哲学差异:ccall 假设"每样本最多一个 ALT",通过群体频率谱全局推断(群体先验);mcall 显式枚举等位基因组合,频率先验来自 QS/panel,更面向多等位与稀有变异。
- 涉及的关键文献(代码注释引用):Li 2011 *Bioinformatics*(`doi:10.1093/bioinformatics/btr509`,AFS/贝叶斯框架)、VDB 来源 PMID:22524474。
