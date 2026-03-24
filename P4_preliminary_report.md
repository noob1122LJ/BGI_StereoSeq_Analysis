# P4 样本（Stereo-seq：cellbin / bin20 / bin50）初步分析报告

> 数据来源：`P4.ipynb` 中对 3 个 h5ad 结果文件的读取与初步 QC/可视化检查。

## 1. 数据与对象完整性检查

从 Notebook 记录看，P4 目录下包含：

- `Y40001A513.cellbin_1.0.adjusted.h5ad`
- `Y40001A513.bin20_1.0.h5ad`
- `Y40001A513.bin50_1.0.h5ad`
- `LC_4B.raw.gef`
- `ssDNA/`

三个 h5ad 文件均存在并可成功读取。

## 2. 三种分辨率数据规模对比（关键事实）

- **cellbin**：`216,796 × 49,350`（obs × vars）
- **bin20**：`249,163 × 49,344`
- **bin50**：Notebook 输出显示已成功读取并进入后续分析流程（同样包含 leiden / spatial / QC 指标）

解释：
- 在该样本中，`bin20` 的 spot/bin 数量高于 `cellbin`，符合“规则网格分箱通常覆盖更连续”的常见现象。
- `cellbin` 保留了细胞分割相关信息（如 `dnbCount`、`area`、`cell_border`），适合细胞层面的下游分析。

## 3. 预处理状态判读

矩阵值检查显示：三个对象的 `X` 非整数、最小值 > 0（非零项），且存在 `layers['log1p']`。
这强烈提示当前 `X` 不是原始 UMI counts，而是经过标准化/变换（很可能 log1p 或等价变换）后的表达矩阵。

对后续分析的影响：
- 差异分析、伪 bulk、部分空间模型（尤其假设原始计数分布的方法）不应直接把当前 `X` 当 raw counts 使用。
- 若要做严格统计建模，优先确认 raw counts 是否在 `layers['counts']`、`raw` 或原始 GEF/矩阵可回溯。

## 4. 空间聚类与 QC 可视化（初步结论）

Notebook 已完成：

- `leiden` 在 spatial 坐标上的可视化（cellbin/bin20/bin50）
- `total_counts`、`n_genes_by_counts`、`pct_counts_mt` 的空间分布可视化

初步可认为：
- 三种分辨率对象均具备可用的空间结构与聚类标签，流程可继续。
- 已有 QC 指标可用于组织边缘低质量区域、线粒体偏高区域、低复杂度区域识别。

> 当前 notebook 输出中未直接给出阈值统计（如过滤前后数量、异常区域面积占比、各 cluster 的 QC 中位数），因此本结论为“可继续分析”的初筛结论。

## 5. 风险点与注意事项

1. **表达矩阵是否为 raw counts 未完全闭环确认**  
   已知 `X` 近似为变换后值；仍需明确 raw 层位置。

2. **cellbin 与 bin20/bin50 的结论不可直接等价**  
   分辨率改变会影响邻域结构、聚类粒度和 marker 稳定性。

3. **mt 比例阈值不能照搬单细胞 RNA 的固定阈值**  
   空间数据与组织类型差异较大，应以样本内分布和空间位置联合判断。

## 6. 建议的后续分析路线（建议按优先级执行）

### A. 数据层与可重复性校验（最高优先级）

- 明确每个对象：
  - raw counts 存储位置（`raw` / `layers` / 外部文件）
  - 当前 `X` 的变换方式（是否 log1p、是否归一化到固定总量）
- 固化版本与参数：记录 stereopy/scanpy 版本、关键参数（neighbors、resolution、HVG 策略）。

### B. 统一 QC 报告（建议形成表格）

对 cellbin/bin20/bin50 分别输出：

- `n_obs`、`n_vars`
- `total_counts` / `n_genes_by_counts` / `pct_counts_mt` 的分位数（P1/P5/P25/P50/P75/P95/P99）
- 过滤阈值方案与过滤后保留比例
- 空间异常区（高 mt、低 counts）占比

目标：让三种分辨率结果可横向比较，而不是只看图。

### C. 聚类稳健性与跨分辨率映射

- 在三个对象上做参数敏感性：
  - Leiden resolution 梯度（如 0.2, 0.4, 0.6, 1.0）
  - 邻居数（15/30/50）
- 进行跨分辨率 label transfer / 最近邻映射：
  - 评估 cellbin↔bin20↔bin50 的 cluster 一致性（ARI/NMI/重叠率）

目标：明确“生物学结构”与“分辨率伪影”的边界。

### D. 组织学联合验证（强烈建议）

- 将 cluster 边界叠加到 ssDNA 图像。
- 检查高 mt 区域是否与组织边缘、破损区域重合。
- marker 基因空间分布是否与形态学区域一致。

### E. 生物学解释与报告交付

- 每个主 cluster 提供：
  - top markers（含 logFC、表达占比）
  - 空间位置描述
  - 候选细胞类型注释（参考已知 marker）
- 输出最小可交付：
  - 1 份 QC 总表
  - 1 份跨分辨率一致性表
  - 1 份 marker+空间图谱 PDF

## 7. 你现在可以直接执行的“下一步清单”

1. 先做 raw layer 核对并统一分析输入（raw vs log1p）。
2. 生成三种分辨率的 QC 分位数表 + 过滤建议阈值。
3. 补做跨分辨率 cluster 一致性评估（至少 ARI/NMI + 混淆矩阵）。
4. 用 ssDNA 做空间形态学交叉验证。
5. 再进入 marker / 注释 / 通路等生物学解释阶段。

---

如果你愿意，我下一步可以直接给你一套“可运行的标准化分析脚本模板”（按 cellbin/bin20/bin50 批量出 QC 表、聚类稳健性和跨分辨率一致性指标）。
