# %% [markdown]
"""
BGI Stereo-seq 多样本整合（Stereopy 版本）
=====================================

你已经完成前 4 个 cell 后，继续执行本脚本。

设计原则：
1) **严格使用 Stereopy API**（`stereo` 命名空间）。
2) 先单样本 QC / 过滤，再做多样本整合。
3) 所有关键阈值参数集中在一个配置区，方便你根据图形结果快速调参。

> 说明：不同版本 Stereopy 在少量函数签名上可能存在差异，
> 若你的本地版本参数名略有不同，请以 `help(对应函数)` 为准。
"""

# %%
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import stereo as st
from stereo.core.ms_data import MSData


# %%
@dataclass
class QCConfig:
    # 细胞/spot 过滤
    min_gene_by_cell: int = 200
    max_gene_by_cell: int = 8000
    min_count_by_cell: int = 500
    max_pct_counts_mt: float = 15.0

    # 基因过滤
    min_cell_by_gene: int = 3

    # 标准化与降维
    hvg_n_top_genes: int = 3000
    pca_n_pcs: int = 50
    neighbors_n_neighbors: int = 15
    neighbors_n_pcs: int = 30
    leiden_resolution: float = 0.6


CFG = QCConfig()


# %%
# 0) 输入：请确保你已有每个样本的 StereoExpData 列表
#    变量名默认为 `samples`，元素为单样本 StereoExpData。
assert "samples" in globals(), "未找到 `samples` 变量，请先准备好多样本 StereoExpData 列表。"
assert isinstance(samples, Iterable), "`samples` 必须是可迭代对象（list/tuple）。"
samples = list(samples)
assert len(samples) >= 2, "至少需要 2 个样本才能进行整合。"


# %%
# 1) 单样本 QC 与过滤
#    根据你前面绘图（基因数分布、UMI 分布、线粒体占比）调节 CFG 参数。
for i, data in enumerate(samples):
    if data.sn is None or str(data.sn).strip() == "":
        data.sn = f"sample_{i}"

    # 1.1 计算 QC 指标
    data.tl.cal_qc()

    # 1.2 细胞过滤（按项目需要可放宽或收紧）
    data.tl.filter_cells(
        min_gene=CFG.min_gene_by_cell,
        max_gene=CFG.max_gene_by_cell,
        min_count=CFG.min_count_by_cell,
        pct_counts_mt=CFG.max_pct_counts_mt,
        inplace=True,
    )

    # 1.3 基因过滤
    data.tl.filter_genes(min_cell=CFG.min_cell_by_gene, inplace=True)

    # 保存 raw，便于后续回溯
    data.tl.raw_checkpoint()

    print(f"[{data.sn}] cells={data.n_cells}, genes={data.n_genes}")


# %%
# 2) 构建 MSData（多样本容器）
ms_data = MSData(samples)
print(ms_data)


# %%
# 3) 归一化、HVG、降维、聚类（Stereopy 多样本流程）
# 注意：部分版本中参数名可能是 n_neighbors / n_pcs / pca_res_key 等，
# 若你本地版本不同，请使用 help(ms_data.tl.<method>) 对照。

ms_data.tl.normalize_total(target_sum=1e4)
ms_data.tl.log1p()

ms_data.tl.highly_variable_genes(
    n_top_genes=CFG.hvg_n_top_genes,
    batch_key="sn",  # 按样本名分层挑选 HVG
)

ms_data.tl.scale(max_value=10)
ms_data.tl.pca(n_pcs=CFG.pca_n_pcs)

# 多样本整合：优先 Harmony
# 如果你的版本函数名不是 harmony_integrate，常见为 harmony / integrate / batch_correct。
ms_data.tl.harmony_integrate(batch_key="sn", basis="X_pca")

ms_data.tl.neighbors(
    n_neighbors=CFG.neighbors_n_neighbors,
    n_pcs=CFG.neighbors_n_pcs,
    use_rep="X_pca_harmony",
)
ms_data.tl.umap(min_dist=0.3)
ms_data.tl.leiden(resolution=CFG.leiden_resolution, key_added="leiden_r06")


# %%
# 4) 可视化检查整合效果（按样本、聚类上色）
ms_data.pl.umap(color=["sn", "leiden_r06"], wspace=0.35)


# %%
# 5) Marker 基因与导出
ms_data.tl.rank_genes_groups(groupby="leiden_r06", method="wilcoxon")
ms_data.pl.rank_genes_groups(n_genes=15, sharey=False)

# 将 top20 marker 导出为表格
rgg = ms_data.tl.result["rank_genes_groups"]
clusters = rgg["names"].dtype.names
marker_df = pd.concat(
    [
        pd.DataFrame(
            {
                "cluster": c,
                "gene": rgg["names"][c][:20],
                "logfoldchanges": rgg["logfoldchanges"][c][:20],
                "pvals_adj": rgg["pvals_adj"][c][:20],
            }
        )
        for c in clusters
    ],
    ignore_index=True,
)
marker_df.to_csv("markers_leiden_r06_top20.csv", index=False)


# %%
# 6) 保存整合结果
out_dir = Path("results")
out_dir.mkdir(parents=True, exist_ok=True)

ms_data.write_h5ms(str(out_dir / "integrated_samples.h5ms"))
print("已保存:")
print(f"- {out_dir / 'integrated_samples.h5ms'}")
print("- markers_leiden_r06_top20.csv")
