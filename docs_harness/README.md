# docs_harness — bcftools SNP Calling 代码梳理

本目录是 `bcftools` 核心代码中 **SNP calling(含 indel)流水线**的结构化梳理文档,基于当前仓库源码(`git HEAD 905ea784`)逐函数通读后整理。

## 内容

- [`snp_calling_architecture.md`](./snp_calling_architecture.md) — 主文档:
  - 两段式流水线总览(`bcftools mpileup` → `bcftools call`)
  - 8 个模块划分(M1–M8)及其输入/输出
  - 每个关键模块的技术路线(含代码位置引用)
  - 统计架构速查表(统计量 ↔ 计算公式 ↔ 代码位置)
  - 主调用链与文件清单
- [`mcall_step_by_step.md`](./mcall_step_by_step.md) — `bcftools call -m` 逐步详解:
  - 16 步全流程分解(初始化 → 位点推断 → 基因型调用 → 输出)
  - 关键数学公式(组合似然、QUAL、GQ、theta 先验)
  - **分级调试输出** `BCFTOOLS_DEBUG_MCALL=1|2|3` 的开关、调试点清单与中间结果复现示例

## 一句话概括

- **阶段一 `mpileup`**:读 BAM → 过滤/BAQ 重比对 → 用 htslib `errmod` 误差模型把每条 read 转成基因型似然 PL → 等位基因排序与偏差统计 → 输出"纯统计"BCF(QUAL=0)。
- **阶段二 `call`**:读 PL → `-c` 走 EM + AFS 贝叶斯(`prob1.c`/`em.c`),`-m` 走质量加权频率 + 等位基因子集枚举(`mcall.c`)→ 输出 GT/GQ/QUAL/gVCF。
