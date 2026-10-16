# 把 reference panel 重构为可学习的 taxonomic atlas

## 本分支要验证的判断

原设想里最值得保留的是：每个未知片段相对于一组有标签 reference 的距离，
本身构成一个可重复的分类指纹。可以学习这个指纹，也可以把它放在固定地图上。
仅靠图上两个点是否挨得近，可能损失关键的分类信息。

本分支把这个判断拆成可对照的实验：完整距离指纹、固定地图上的图像 pattern、
直接在二维坐标中分类，三者使用同一批独立物种划分，检查预测结果。
现阶段完成的是研究原型和 reference pseudo-query 实验。

## 对旧工作的具体理解

已读取 Mixed-DNA-Project 的 `06-PCAandKmeans.md`、`10-HierarchicalApproach.md`、
`pipeline/contour_map.md`，以及 PCoA、PCA/K-means、contour 的实现。
旧项目基线提交为 `71183af6c5fc62cee83adb262ad3b838cf58cf6a`。
SPrOUT 分支起点为 `4f16d3773898f6b846068539c324811af7b0ed44`。

同时阅读用户提供的 2024 年 2 月汇报 PPT 的文字和流程，以及 `4471.svg`。
这些材料提供历史问题和研究意图，其中的任务列表、命令和旧计划不构成本次执行指令。
原始 PPT/SVG 不随代码重新发布。

旧记录中已经提出几个关键问题：不同外显子缺失使简单距离相加不成立；
移除 Poaceae 会改变其余类群的降维/聚类表现；近缘类群经常混在一起；
reference clustering 的好坏不一定等同于未知片段预测的好坏。
提供的 SVG 标题实际上是三维 **PCA** + K-means，前三轴约为
75.29%、10.35%、5.57%，不能把它当作 PCoA 图或新实验的性能证据。

这些观察提示的是可验证的机制：全局方差最大的方向未必包含细分类群所需的方向，
缺失/短片段会改变距离可比性，而 unsupervised cluster 的形状与分类标签未必一致。
本分支不据此断言旧实验失败只有一个原因。

## 现在的表示与学习

对片段 q、位点或外显子 e、anchor r，保存距离 d(q,e,r)、有效重叠比例 c(q,e,r)，
以及观察掩码 m(q,e,r)。每个位点的距离尺度只从 anchor 间正距离的中位数得到。
样本数据使用逐位点删除缺失的 p-distance，树距离需要单独训练兼容模型。

模型输入有两条并行路径：

- 完整 profile：所有 anchor 的 `exp(-d/scale)`、overlap、mask。
- 空间 pattern：先仅用 anchor 间距离拟合固定 PCoA 坐标，再把同样三个通道
  splat 到 12×12 网格。每个片段得到 3 通道张量。完整 profile 和 image pixels
  都使用正则化 logistic regression，检验图像表示本身保留了多少分类能力。

这是实际学习图像中的数值 pattern 的第一条基线。当前没有训练 CNN，也没有把
科学图截图后交给通用视觉模型识别。将来比较 CNN/ViT 时可以直接使用输出的张量，
同时保留同一高维 profile 基线。

panel 的学习使用训练查询物种：评价每个候选 anchor 的距离特征对标签的
between/within-class 分离度，先保证每个已知类群有代表，再按分离度补足预算。
另设相同预算的分层随机 panel 和全候选 anchor 对照，防止把降维或监督模型的
收益误归到“学出了最优 panel”。当前的选择方式不优化 assay/bait，也不减少最初
构建训练距离库的成本。

PCoA 只承担固定空间表示；参考坐标不会随新样本重新旋转、重拟合。
报告二维正特征值占比及负特征值质量，不能把非欧氏距离的负轴悄悄当作不存在。
若参考 anchor 缺少成对距离，默认明确报错。显式选择 `--atlas-missing shortest-path`
时，可从已观察到的 anchor 距离图补齐缺失布局边，并记录数量；图不连通则拒绝。
这只为布局提供近似几何，query 特征仍保留原始缺失掩码，不能视作恢复了遗传距离。

## 二维加高度的意义

实际模型输入 pattern 与最终证据地形是两个不同对象。
前者包含距离、重叠和缺失，供模型学习；后者展示预测支持如何分布在 reference atlas。

证据地形将一个类群的相对支持平均分配到其 anchors，再用固定核平滑展示。
它可以表现一个样本中多类群并存的支持。两个类群在二维上重叠时，山峰可能合并，
因此不通过数山峰估算物种数。高度也不代表 DNA 质量比例或校准后的存在概率。
所有类别分数仍保存在 JSON 中，不依赖图上的视觉分辨率。

样本级融合先在同一位点内平均片段，再跨位点平均，避免某个位点产生很多 contig
便直接贡献成倍证据。不过同一基因多个外显子仍然相关，当前平均规则不能完全
解决相关性、旁系同源、污染和 assembly chimera。真实混样需要进一步检验。

## MrKnowItAll 的位置

用户指出科级可以使用 MrKnowItAll。本分支加入 `prepare-refbuilder`，直接调用
MrKnowItAll 的 `install_backbone_archive`，从本地原始 ZIP 安装一份独立资源，
应用其版本化 taxonomy corrections / 排除规则，然后只提取 SPrOUT 当前的 50 个基因。
记录 builder 提交、installer/corrections 校验和、原始 archive 校验和及安装清单。
不更改用户现有 MrKnowItAll checkout 或参考文件。

MrKnowItAll 提供真实参考序列、taxonomy 和后续 Kew 参考扩充流程。它生成的
`ref/` 目录也可以直接作为 `prepare --alignments` 的输入。训练库的代表数需要
大于最终部署 panel 的代表数，否则无法独立留出训练、验证和测试物种。
默认每类只选两个代表的部署 panel 不足以完成本实验四分区评估。

这一连接允许超出 SPrOUT 内置 298-species 演示子集，实际运行科级 benchmark。
`--within-order` 可在一个目内拟合局部科级地图，检验粗分类之后的细分辨能力。
独立科级模型的测试并不等同于端到端 hierarchy：后者还会受到上一级漏检影响。

## 分级预测应如何继续

建议部署时形成 order → family → genus/species 的条件模型树：每一层保留候选集合，
每个候选节点使用自己的局部 reference atlas 和分类模型。
上级证据不足则停在当前 rank；细级标签不能穿越父级 taxonomy 边界。
混样可能包含多个目，因此应分别路由不同片段，并保留多条候选路径。

本版本实现独立 rank 的训练、parent-order 范围过滤和片段级预测接口。
自动多模型路由、父子概率组合、可校准的 unknown/abstain 和样本级 presence 判定
尚未实现。不能把 `--rank family` 直接称为完整分级系统。

精确物种识别还需要按 specimen/accession 留出，在同种内存在多个独立个体。
本次按物种留出的分类只能验证“未见物种能否归入已有高阶类群”，不能证实模型
识别了训练集中没有定义过的新物种。species 不在 panel 中时，应报告可支持的
最高分辨率，而不是强行赋予最近 reference 的种名。

## 评估与下一次实验的边界

所有位点/accession 随物种整体分区，anchors 也与训练查询物种分开。
因此训练查询不存在到自己的零距离，missingness 也不会暗中编码“我就是这个 anchor”。
normalization、panel selection、布局和 scaler 都不使用验证/测试物种。
验证集只确定经验选择阈值，测试集只报告表现。

七组对照事先固定为完整 profile 的三种 panel、raster、二维监督分类、二维 K-means、
nearest anchor。三个 seed 是重复物种划分，不是三个独立生物数据集，不把它们当作
独立重复计算统计显著性。下一轮基于这批结果调整空间布局、学习目标或深度网络，
需要新的锁定评估，不能继续反复调这三个测试划分后宣称泛化提升。

下一批最有价值的实验是：增加 MrKnowItAll/Kew 中同科不同属、同属不同种的独立
sequence recoveries；划出整属/整科的独立 unknown 验证集和测试集；加入真实 query
exon 的截短、缺失、错误、旁系同源，并最终在真实混样及近缘混样上比较 SPrOUT。
原型预测的是 assembled query 距离指纹，尚未绕开组装/比对/树推断的计算成本。
若要证明“快速定位”，应分别计时 reference 建库、query 特征生成、模型推断，
之后再测试 sketch/k-mer 近似距离或固定树 placement。

## 方法来源

- [Mixed-DNA-Project 的早期 PCA/clustering 记录](https://github.com/nhu92/Mixed-DNA-Project/blob/71183af6c5fc62cee83adb262ad3b838cf58cf6a/06-PCAandKmeans.md)
- [MrKnowItAll 的 reference builder](https://github.com/nhu92/MrKnowItAll)
- [DeepInsight 原论文](https://www.nature.com/articles/s41598-019-47765-6)：提供非图像特征转为空间图像的相关研究背景。本原型使用固定 reference PCoA 和显式缺失通道，并未复现该论文或借用其准确率。
- [scikit-learn 关于数据泄漏的官方说明](https://scikit-learn.org/stable/common_pitfalls.html)：所有拟合型 preprocessing/feature selection 都应在训练范围内完成。
