# RiboGRAM 核苷酸依赖图谱：首轮结果

## 一句话结论

仿照 Relome 的逐碱基饱和突变框架，RiboGRAM 在 3,868 个 Rfam 家族代表上恢复了显著的 RNA 结构依赖：距离匹配 AUC 均值为 **0.834**，**97.6%** 的家族高于随机水平；这些结果证明模型内部存在可重复、家族特异的有向依赖架构，但尚不能把未注释边解释为新的物理接触。

## 1. 问题与分析边界

参考《The relational architecture of the human genome》所采用的 nucleotide-dependency analysis，我们没有直接扫描人类基因组，而是扫描 RiboGRAM 的自然适用域——非编码 RNA。首轮目标是先回答三个基础问题：

1. 逐碱基替换产生的模型依赖能否恢复 Rfam 共识配对？
2. 依赖是否具有方向性、互惠性和家族特异性？
3. 超出 Rfam 配对的强边是否高于保持序列组成的随机背景？

本轮选择 Rfam 15.1 中至少具有一个可映射共识配对、长度为 20–1,024 nt 的一个无缺口代表序列。共覆盖 3,868 个家族、450,883 nt；其中 Gene/ncRNA 3,299 个、Cis-regulator 552 个、Intron 17 个。

## 2. 饱和突变定义

对参考 RNA 中每个源位置 \(i\)，依次替换为其余三种碱基。对每个目标位置 \(j\) 和四种目标碱基，计算替换前后 RiboGRAM 预测的 log2 odds 变化。定向依赖分数定义为：

\[
D_{i\rightarrow j}=\max_{a'\neq a_i,\,b\in\{A,C,G,U\}}
\left|\operatorname{logit}_2 P(b_j\mid x_{i\leftarrow a'})-
\operatorname{logit}_2 P(b_j\mid x)\right|.
\]

因此 \(D_{i\rightarrow j}\) 与 \(D_{j\rightarrow i}\) 可以不同。对角线自依赖被置零。完整扫描包含 1,356,517 条参考或突变输入，分成四个计算 shard，在 4 × NVIDIA A800 80GB 上约 17.5 分钟完成，失败数为 0。

## 3. 结构恢复结果

- 以 Rfam 共识配对为正例、其余位置对为负例，家族 ROC-AUC 均值为 **0.8383**，中位数为 **0.8896**。
- 为排除“近距离位置天然分数更高”的影响，在相同序列距离分箱内比较，距离匹配 AUC 均值仍为 **0.8341**。
- **3,775/3,868（97.6%）** 的家族距离匹配 AUC 高于 0.5；2,702 个家族高于 0.75。
- 已知配对的依赖强度相对非配对的中位富集为 **4.49 倍**。
- 对每条已知配对，在相同距离的非配对背景中计算百分位；家族层面的中位结果显示，**46.4%** 的已知配对进入同距离背景的前 1%。

这些结果支持一个明确结论：RiboGRAM 的 substitution-dependency map 不是简单的距离效应，而是能够恢复大量已知 RNA 结构关系。

## 4. 家族与长度分层

距离匹配 AUC 的主要分层结果为：Gene/ncRNA **0.852**、Intron **0.781**、Cis-regulator **0.731**。典型亚型中，miRNA 为 **0.926**、CD-box snoRNA 为 **0.884**，而 IRES 为 **0.591**。

长度效应明显：20–80 nt、81–160 nt、161–320 nt、321–512 nt、513–1,024 nt 的距离匹配 AUC 均值依次为 **0.869、0.853、0.704、0.652、0.626**。长度与距离匹配 AUC 的 Spearman 相关为 **−0.317**。

因此不能只展示 tRNA、miRNA 等短 RNA 的成功案例。更可靠的表述是：RiboGRAM 对短至中等长度、结构约束较强的 RNA 恢复最好；长 RNA 的依赖分辨率明显下降，可能同时受到模型有效上下文、复杂结构、代表序列选择和参考标注稀疏性的影响。

## 5. 方向性与互惠性

- 已知配对两个方向的中位不对称度为 **0.127**。
- 已知配对的中位互惠度为 **0.775**。
- q99.9 强边从低坐标指向高坐标的比例为 **0.483**，没有明显的全局 5′→3′ 单向偏置。

这表明 RiboGRAM 的边不是完全对称的结构分数，但许多已知配对在两个方向上都能得到响应。更值得研究的是家族内部的 source hub、target hub 和局部非对称模块，而不是宣称存在统一的方向规律。

## 6. 二核苷酸保持的 null 校准

为避免碱基组成、GC、长度和局部序列距离造成伪富集，构建了 136 条精确保留二核苷酸计数的随机 RNA，共评估 7,822,452 条定向 null 边。阈值按 RNA 长度、GC 含量和序列距离匹配。

- 在真实 RNA 中，超过匹配 q99.9 阈值的位置对占 **2.04%**；null 中为 **0.180%**，富集 **11.32 倍**。
- 真实 q99.9 强边的互惠率为 **22.7%**；null 为 **11.9%**。
- 真实 q99.9 强边中 canonical/wobble 碱基组合占 **42.6%**；null 为 **38.3%**。
- q99.9 阈值召回 **39.7%** 的 Rfam 共识配对；q99 阈值召回 **51.9%**。

全体 38,211,716 个候选位置对中共有 780,489 个 q99.9 强边，其中 46,537 个为 Rfam 共识配对，733,952 个不在当前 Rfam 配对集合中。后者只能称为**候选模型依赖**，不能称为新 RNA 分子相互作用。

## 7. 当前 Atlas 实现

Rfam Atlas 的序列详情页现在包含三个彼此区分的分析层：

1. masked-nucleotide contextual recovery；
2. zero-shot APC-corrected attention pairing；
3. directed saturation-mutagenesis dependency map。

第三层采用与 Relome 一致的弧图语义：蓝色表示 5′→3′ 依赖，绿色表示 3′→5′ 依赖，透明度编码依赖强度，虚线表示 Rfam 共识配对。网页只加载稀疏、null 校准后的强边，每家族最多展示 250 条弧线。

共生成 3,868 个压缩家族数据包、698,660 条展示弧线，总大小约 10.7 MB。由于主 Atlas 仍使用 Rfam 15.0，当前有 3,811 个代表序列可精确链接；另有 49 个 15.1 新家族和 8 个代表序列发生版本变化。

## 8. 可支持与不可支持的结论

当前数据支持：

> RiboGRAM substitution-dependency maps recover Rfam structural relationships and reveal family-specific directed model architectures.

当前数据不支持：

- “发现了 733,952 个新的 RNA 物理相互作用”；
- “依赖边等同于真实碱基配对概率”；
- “所有 RNA 家族都具有相同质量的恢复能力”；
- “模型存在统一的 5′→3′ 因果方向”。

## 9. 下一阶段

下一阶段应从“边恢复”进入“RNA relational modules”：在家族内识别 block、anti-diagonal、offset-diagonal、source hub、target hub 等重复几何结构；在同一家族的多个 Seed 序列中检验模块保守性；再与 Rfam 共识结构、实验 probing、功能位点和已知致病变异交叉验证。只有经过独立数据支持的模块，才适合升级为生物学机制假说。

## 10. 代码与结果位置

- 远程计算目录：`/home/wanghui/ribogram/riborel`
- 本地分析代码：`C:\Users\A\Documents\ChatGPT\RiboScan\research\riborel`
- 本地结果汇总：`C:\Users\A\Documents\ChatGPT\RiboScan\research\riborel\remote_results`
- Atlas 数据包：`C:\Users\A\Documents\ChatGPT\RiboScan\data\dependencies`

参考：Krause et al., *The relational architecture of the human genome*, bioRxiv (2026), DOI 10.64898/2026.10.15.759350；Tomaz da Silva et al., *Nucleotide dependency analysis of genomic language models detects functional elements*, Nature Genetics 57, 2589–2602 (2025)。
