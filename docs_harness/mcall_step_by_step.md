# `bcftools call -m` 逐步详解与调试输出指南

> 配套文档:`snp_calling_architecture.md`(总体架构与 8 模块划分)。本文把 **M6(mcall)** 拆到"每一步",并在代码中埋了分级调试输出,方便复现中间结果。
>
> 代码基准:git HEAD `905ea784`;行号以本仓库 `mcall.c` 当前版本为准(插入调试代码后行号已偏移,下文行号为**当前**文件行号,可用 `grep -n "step=" mcall.c` 复查)。

---

## 1. 总览:`call -m` 的处理流水线

```
bcftools mpileup -Ob -o x.bcf ref.fa reads.bam     # 阶段一:纯统计(产生 PL/QS/AD/I16)
bcftools call -m -Ov x.bcf > x.vcf                 # 阶段二:贝叶斯推断(本文范围)
```

`-m`(multiallelic calling model)是默认的 call 模型,核心特点:

- 位点水平:在**候选等位基因集合**上做最大似然选择(最多三等位组合),决定"这个位点是不是变异、留下哪些 ALT";
- 样本水平:在入选等位基因上做 MAP 基因型调用,输出 GT/GQ/GP;
- 统计源:完全依赖 mpileup 阶段一产出的 `FORMAT/PL`(基因型似然)与 `INFO/QS`(等位基因质量和,充当 AF 先验)。

主循环在 `vcfcall.c main_vcfcall`(L1198 起),每个位点调 `mcall()`(`mcall.c` L1518)。

---

## 2. 逐步分解(共 4 个阶段、16 步)

### 阶段 0:初始化(仅一次,进程启动时)

**S0 `mcall_init`**(`mcall.c` L422)

| 子项 | 位置 | 内容 |
|---|---|---|
| 样本分组 | `init_sample_groups` L237 | 默认全部样本一组;`-G` 文件或 `-G -`(每样本一组)分组 |
| PL→P 查表 | `call_init_pl2p` L82 | `pl2p[PL] = 10^(-PL/10)`,PL 上限 255 |
| 数组预分配 | L436-440 | `als_map`、`pl_map`、`gts`(按二倍体) |
| trio(-T) | L441-448 | `CALL_CONSTR_TRIO` 时初始化家族结构与 CGT/UGT 头 |
| `-C` 模式 | L449 | `vcmp_init` 等位基因比较器 |
| 头声明 | L451-463 | 追加 GT/GQ/GP/AC/AN/DP4/MQ/PV4 头行 |
| **theta 先验** | L465-486 | Watterson 校正:`aM = Σ_{i=2}^{2N} 1/i`,再 `theta = log(theta_raw * aM)`;`theta_raw` 默认 `1.1e-3`(`vcfcall.c` L1147) |

### 阶段 1:主循环预处理(每个位点,`vcfcall.c` L1207-1237)

**S1 记录过滤与符号等位基因识别**

1. `-v/--variants-only`(CF_INDEL_ONLY)、`-i`(CF_NO_INDEL)、`-a`(CF_ACGT_ONLY)过滤(L1209-1211);
2. 识别 unseen 等位基因:扫描 ALT 中的 `X`/`<X>`/`<*>`(L1213-1222),记入 `args.aux.unseen`;
3. 判定纯参考位点 `is_ref`(L1224):`n_allele==1` 或只有 REF+`<*>` 两等位;
4. `bcf_unpack(BCF_UN_ALL)` + `set_ploidy`(L1226-1227);
5. 调 `mcall(&args.aux, bcf_rec)`(L1245)→ 返回 `-2` 跳过、`0` 非变异、`>0` 变异;
6. 输出:`--variants-only` 过滤、`gvcf_write` 折叠、`bcf_write1`(L1247-1255)。

### 阶段 2:`mcall()` 位点级推断(`mcall.c` L1518)

**S2-1 强制等位基因 `-C`**(L1531-1537,`mcall_constrain_alleles` L1370)

仅当 `-C alleles` 启用:把 targets 中的等位基因映射进记录,重排 PL/QS/Number=R 字段;targets 里有但 VCF 没有的等位基因用 unseen 的似然估算。之后 `nals_ori` 取**约束后**的等位基因数。

**S2-2 读取并校验 PL**(L1542-1545)

```c
call->nPLs = bcf_get_format_int32(hdr, rec, "PL", &PLs, &mPLs);
if ( nPLs != nsmpl*nals*(nals+1)/2 && nPLs != nsmpl*nals ) error(...);
```

合法形态:全体二倍体 `nsmpl·C(nals,2)` 或全体单倍体 `nsmpl·nals`。

**S2-3 PL → P(D|G) 概率**(L1547-1549,`set_pdg` L520)

- 每个基因型 `g`:`pdg[g] = 10^(-PL[g]/10)`,再除以组内和归一化(Σ pdg = 1);
- missing 处理:未观察等位基因(`unseen`)缺失的基因型用 AX/XX 的似然补齐;全缺失的样本 `pdg≡0`(`FLAT_PDG_FOR_MISSING=0` 时);
- 注意:pdg 是**似然**(P(D|G)),不是后验。

**S2-4 组装 AF 先验 qsum(逐等位基因,充当 f_x = QS/N)**(L1552-1660)

| 情况 | 位置 | 公式 |
|---|---|---|
| 单组(默认) | L1554-1606 | `qsum[x] = INFO/QS[x]`(mpileup 的等位基因质量分总和);不足 `nals_ori` 的补 0 |
| 多组 `-G` | L1607-1651 | 按 `-G [AD|QS:]` 指定的 FORMAT 标签聚合:`qsum[x] += AD[x]/ΣAD`(每样本先归一) |
| panel `-F AN,AC` | L1653-1683 | `qsum'[x] = (qsum[x] + 0.5·AC[x]) / (nsmpl + 0.5·AN)` |
| **归一化** | L1686-1699 | 每组 `qsum[x] /= Σqsum`,使先验和为 1 |

然后删除 INFO/QS(L1701)。

**S2-5 等位基因集合选择**(L1708-1741,`mcall_find_best_alleles` L656)

每组独立枚举(至多三等位组合),打分并挑最大似然组合:

- **1 个等位基因**(L660-686):只考虑纯合 `ia/ia`,`lk_tot = Σ_s log P(D_s | ia/ia)`;`ia≠0` 时加先验 `theta`;`ia==0` 的结果记作 `ref_lk`;
- **2 个等位基因**(L688-727):`fa,f b = qsum 归一`;每样本 `P(D|{ia,ib}) = fa²·P(ia/ia) + fb²·P(ib/ib) + 2fa·fb·P(ia/ib)`;非参考等位每个加 `theta`;
- **3 个等位基因**(L729-783):同上推广到 6 个基因型项;
- 累积:`max_lk`(最优组合)、`lk_sum`(`logsumexp2` 累加所有组合,用于归一化)、`ref_lk`;
- 组间合并:`call->als_new |= grp->als`;`max_qual = max over groups of -4.343·(ref_lk - logsumexp2(lk_sum, ref_lk))`(L1732)。

**S2-6 REF 强制在场与变异判定**(L1743-1748)

```c
if ( !(als_new&1) ) als_new |= 1;
is_variant = (als_new==1) ? 0 : 1;
if ( CALL_VARONLY && !is_variant ) return 0;
```

**S2-7 确定输出等位基因集 nals_new 与映射表**(L1752-1767)

- unseen 的处理(`CALL_KEEP_UNSEEN`)、`-k/--keep-alts`(CALL_KEEPALT)强制保留全部 ALT;
- `init_allele_trimming_maps`(L615):`als_map[old]=new`、`pl_map[new]=old`。

**S2-8 基因型调用(三选一分支)**(L1772-1830)

| 分支 | 条件 | 动作 |
|---|---|---|
| A | `als_new==1`(只有 REF) | `mcall_set_ref_genotypes` 全部 0/0(或深度 0 → ./.),**删除 PL** |
| B | 非变异(因 `-A/--keep-alts` 强制) | `mcall_set_ref_genotypes` + 保留并裁剪 PL |
| C | 变异 | `mcall_call_genotypes` 做 MAP 调用(见 S3),写 GP/GQ,裁剪 PL |

**S2-9 裁剪 Number=R 的 INFO/FORMAT 字段**(L1832-1833,`mcall_trim_and_update_numberR` L1265)

被丢弃等位基因对应位置的 AD/QS 等按 `als_map` 重排(或删除到 1 个)。

**S2-10 设置 QUAL**(L1835-1845)

```c
if ( nAC )            rec->qual = max_qual;                     // 变异位点质量
else if ( lk_sum != -HUGE_VAL ) rec->qual = -4.343*(lk_sum - logsumexp2(lk_sum, ref_lk));  // 参考位点
else if ( call->ac[0] ) rec->qual = -4.343*theta;              // 无高质量读数支持
else                  bcf_float_set_missing(rec->qual);
```

**S2-11 AC/AN**(L1847-1861)`:AC[x] = Σ_s GT 中等位基因 x 的拷贝数`;`AN = Σ AC`。

**S2-12 写回记录**(L1863-1880):裁剪 alleles、写 GT、DP4/MQ/PV4(`test16` 来自 `call.c`),移除 I16。

### 阶段 3:样本级基因型调用(`mcall_call_genotypes` L828)

**S3-1 MAP 基因型扫描**(L843-938,逐样本)

- 只扫描入选等位基因(`grp->als` 位掩码);
- 纯合 `ia/ia`:`lk = pdg[iaa] · qsum[ia]²`(单倍体 `·qsum[ia]`);
- 杂合 `ia/ib`:`lk = 2·pdg[iab] · qsum[ia]·qsum[ib]`;
- `best_lk` 最大的基因型 → GT(`USE_PRIOR_FOR_GTS=0`,GT 不用 theta 先验);
- 无深度(全 0 pdg)→ `./.`、`gps[0]=-1`。

**S3-2 GP 归一化与 GQ**(L939-985)

- `GP[g] = lk[g]/Σ lk`(只对 nmax 个有效基因型);
- `GQ = -4.34294·log(1 - max/Σ)`,上限 `INT8_MAX`。

---

## 3. 关键数学速查

| 量 | 定义 |
|---|---|
| `pdg[g]` | `P(D_s|G=g) ∝ 10^(-PL/10)`,每样本归一 |
| `qsum[x]` | AF 先验,mpileup QS 归一(或 -G/-F 混合),Σ=1 |
| 组合似然(2 等位) | `P(D|{a,b}) = Π_s [fa²P(ia/ia)+fb²P(ib/ib)+2fa·fb·P(ia/ib)]` |
| 先验惩罚 | 组合中每个非 REF 等位基因乘 `exp(theta)`(log 域加 `theta`) |
| 位点 QUAL | `-4.343·(ref_lk − logsumexp2(lk_sum, ref_lk))` |
| GQ | `-4.343·log(1 − max_g GP[g])` |
| theta | `log(theta_raw · Σ_{i=2}^{2N} 1/i)`,默认 `theta_raw=1.1e-3` |

---

## 4. 调试输出(本次新增)

### 4.1 开关

环境变量 **`BCFTOOLS_DEBUG_MCALL`**,三档,默认关闭:

| 值 | 级别 | 内容 |
|---|---|---|
| 0 / 未设置 | 关 | 零输出、零行为变化(默认) |
| 1 | 站点级 | 每站点:begin、qsum 先验、best_alleles、site、alleles、genotypes 分支、qual、ac_an、final |
| 2 | 样本级 | 加上:每样本 PL、pdg、gt_sample(GT+GP_raw)、gq_sample(GQ+GP) |
| 3 | 组合级 | 加上:find_best_alleles 中每个候选等位基因组合的 lk_tot/fa/fb/fc |

所有输出走 **stderr**,不污染 stdout 的 VCF/BCF 流。测试时:

```bash
BCFTOOLS_DEBUG_MCALL=1 bcftools call -m -Ov x.bcf > x.vcf 2> mcall.log
BCFTOOLS_DEBUG_MCALL=2 bcftools call -m -Ov x.bcf 2> mcall.log
BCFTOOLS_DEBUG_MCALL=3 bcftools call -m -Ov x.bcf 2> mcall.log   # 输出量大,慎用
```

### 4.2 输出行格式

统一前缀 `[MCALL DBG] site=CHR:POS step=...`,key=value 空格分隔,一行一个事件;数组值逗号分隔:

```
[MCALL DBG] site=chr1:12003 step=begin nals_ori=3 nsmpl=1 nsmpl_grp=1 theta=-6.80914 unseen=2 ploidy=all-diploid flag=0x100
[MCALL DBG] site=chr1:12003 step=PL sample=0 vals=0,13,37,44,60,83
[MCALL DBG] site=chr1:12003 step=pdg sample=0 vals=1,0.0501,0.002,0.001,0.0001,5e-09
[MCALL DBG] site=chr1:12003 step=qsum_raw group=0 vals=0.9942,0.0058,0
[MCALL DBG] site=chr1:12003 step=qsum group=0 nsmpl=1 vals=0.9942,0.0058,0
[MCALL DBG] site=chr1:12003 step=best_alleles group=0 als=0x3 nals=2 max_lk=-3.14159 ref_lk=-1.2 lk_sum=-2.5
[MCALL DBG] site=chr1:12003 step=site als_new=0x3 max_qual=24.5
[MCALL DBG] site=chr1:12003 step=alleles is_variant=1 als_new=0x3 nals_new=2 als_map=0,1,-1
[MCALL DBG] site=chr1:12003 step=genotypes mode=variant
[MCALL DBG] site=chr1:12003 step=qual nAC=1 qual=24.5 lk_sum=-2.5 ref_lk=-1.2
[MCALL DBG] site=chr1:12003 step=ac_an ac=1,1,0 an=2
[MCALL DBG] site=chr1:12003 step=final nals=2 alleles=A,C qual=24.5 ret=2
```

### 4.3 调试点清单(与步骤对应)

| 调试点 | step 名 | 级别 | 位置(mcall.c) | 说明 |
|---|---|---|---|---|
| 初始化 | `init` | 1 | L478 | theta_raw、aM、校正后 theta、n、分组数 |
| 位点入口 | `begin` | 1 | L1527 | nals_ori、nsmpl、分组数、theta、unseen、flag |
| `-C` 后 | `constrain_alleles` | 1 | L1537 | 约束后 nals、unseen(仅 -C) |
| 读 PL | `PL` | 2 | L1558 | 每样本原始 PL 数组(二倍体顺序) |
| set_pdg | `pdg` | 2 | L1577 | 每样本归一化后的 P(D|G) |
| 单组 QS | `qsum_raw` | 1 | L1599 | INFO/QS 原始值 |
| 多组聚合 | `qsum_grouped` | 2 | L1647 | -G 模式聚合结果 |
| panel 混合 | `qsum_panel` | 2 | L1681 | -F 混合后(仅 -F) |
| 归一化 | `qsum` | 1 | L1699 | **最终 AF 先验**(Σ=1) |
| 组合枚举 | `combo` | 3 | L683/723/775 | 每个候选组合:combo 位掩码、fa/fb/fc、lk_tot |
| 组最优 | `best_alleles` | 1 | L1724 | 每组:als 位掩码、max_lk、ref_lk、lk_sum |
| 位点判定 | `site` | 1 | L1739 | 合并后 als_new、max_qual |
| 输出等位集 | `alleles` | 1 | L1765 | is_variant、nals_new、als_map |
| 分支 | `genotypes` | 1 | L1773/1779/1787 | mode=ref_only / ref_with_PL / variant |
| 样本 GT | `gt_sample` | 2 | L927 | 每样本 GT、best_lk、GP_raw(未归一) |
| 样本 GQ | `gq_sample` | 2 | L980 | 每样本 GQ、归一化 GP |
| 位点 QUAL | `qual` | 1 | L1843 | nAC、qual、lk_sum、ref_lk |
| AC/AN | `ac_an` | 1 | L1853 | ac 数组、an |
| 结束 | `final` | 1 | L1893 | 输出等位基因串、qual、返回码 |

字段说明:

- `als`/`als_new`:位掩码,bit i 表示第 i 个等位基因入选(bit 0 = REF)。`0x3` = REF+第一个 ALT;
- `als_map`:`旧等位基因下标 → 新下标`,`-1` = 被丢弃;
- `combo`:该候选组合的位掩码;
- `GP_raw` 与 `GP`:`GP_raw` 是未归一化的 `lk`(可含 0 占位),`GP` 是输出到 FORMAT/GP 的归一化后验(截断到 nmax 个);
- `gt`:调用结果,`0/1` 形式,`.` 表示缺失;`best_lk` 是赢得 MAP 的原始 lk;
- `ploidy`:`all-diploid` 表示未用 `-P`;per-sample 时每样本 GT 行带 `ploidy=N`。

### 4.4 中间结果复现示例

按站点抽取某个位点的完整中间状态:

```bash
BCFTOOLS_DEBUG_MCALL=2 bcftools call -m -Ov x.bcf 2> mcall.log
grep "chr1:12003" mcall.log
```

验证"qsum 归一化"这一步(Σ 应 ≈ 1):

```bash
grep "step=qsum " mcall.log | head -1 | \
  awk -F'vals=' '{split($2,a,","); s=0; for(i in a) s+=a[i]; print "sum="s}'
```

验证 GT 与 AC 一致性(ac_an 行的 ac 应为各样本 GT 等位基因计数之和):

```bash
grep "step=ac_an" mcall.log | head -1
grep "step=gt_sample" mcall.log
```

对比不同先验对选择的影响(同一 BCF 重放):

```bash
bcftools call -m -Ov x.bcf 2>&1 >/dev/null | grep "step=best_alleles"   # 默认 theta
bcftools call -m -P 1e-2 -Ov x.bcf 2>&1 >/dev/null | grep "step=best_alleles"
```

---

## 5. 实现说明(本次代码改动)

- 全部改动位于 `mcall.c`,约 12 处 `fprintf(stderr, ...)`,全部包在 `if ( mcall_dbg_level()>=N )` 中;
- `mcall_dbg_level()`(L46)惰性读取一次环境变量,静态缓存;
- 辅助:`MCALL_DBG_PFX`(L56)打印站点前缀;`mcall_dbg_vals_i32/f/d` 打印数组;`mcall_dbg_gt`(L77)把编码 GT 转成 `0/1` 文本;
- 唯一非打印改动:`mcall()` 入口 `call->rec = rec`(L1521),供无 rec 参数的内部函数打印站点名(`call->rec` 原本未被 mcall 路径使用);
- **默认关闭,不改变任何既有输出与行为**;开启时仅向 stderr 追加文本,stdout 的 VCF/BCF 字节流不变(可 diff 验证)。

### 5.1 回归验证方法

```bash
# 1) 默认行为不变(无调试输出)
bcftools call -m -Ov x.bcf > a.vcf 2> err_default.log
[ -s err_default.log ] && echo "unexpected stderr" || echo "OK: stderr empty"

# 2) 开调试后 stdout 字节一致
BCFTOOLS_DEBUG_MCALL=3 bcftools call -m -Ov x.bcf > b.vcf 2> /dev/null
cmp a.vcf b.vcf && echo "OK: stdout identical"
```
