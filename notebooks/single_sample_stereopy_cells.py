# %% [markdown]
"""
# 单样本 Stereo-seq（Stereopy 官方流程风格）

这个脚本按 Jupyter cell 组织，覆盖：
- 读取数据
- 线粒体基因识别与 QC（重点解决 `pct_counts_mt=0`）
- 细胞/基因过滤
- 标准化、HVG、PCA、邻居图、UMAP、Leiden
- marker 导出与结果保存

> 你可以直接把每个 `# %%` 当作 notebook cell 运行。
"""

# %%
from __future__ import annotations

from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd
import stereo as st




def _decode_if_bytes(v):
    if isinstance(v, (bytes, bytearray)):
        return v.decode("utf-8", errors="ignore")
    return str(v)


def extract_gene_symbols_from_gef(gef_path: str) -> tuple[pd.Index | None, str | None, str | None]:
    """Try reading gene symbol directly from common GEF HDF5 gene tables.

    Returns
    -------
    (symbols, table_path, column_name)
    """
    candidate_tables = [
        "/geneExp/bin1/gene",
        "/geneExp/bin1/geneExp",
        "/geneExp/gene",
        "/gene/gene",
        "/gene",
    ]
    candidate_symbol_fields = ["gene_name", "geneName", "symbol", "geneSymbol", "gene"]

    try:
        with h5py.File(gef_path, "r") as f:
            for table_path in candidate_tables:
                if table_path not in f:
                    continue
                ds = f[table_path]
                if not hasattr(ds, "dtype") or ds.dtype.names is None:
                    continue

                names = set(ds.dtype.names)
                symbol_col = next((c for c in candidate_symbol_fields if c in names), None)
                if symbol_col is None:
                    continue

                vals = [_decode_if_bytes(x) for x in ds[symbol_col]]
                vals = pd.Index(vals)
                if len(vals) > 0 and vals.str.len().mean() > 0:
                    return vals, table_path, symbol_col
    except Exception as e:
        print(f"读取 GEF 基因注释失败（将回退默认读取）: {e}")

    return None, None, None


def is_ensg_like(index: pd.Index) -> bool:
    if len(index) == 0:
        return False
    probe = index[: min(500, len(index))].astype(str)
    frac = np.mean([x.startswith(("ENSG", "ENSMUSG", "ENSGALG", "ENSDARG")) for x in probe])
    return frac > 0.5

# %% [markdown]
"""
## 1) 参数区（先按默认跑，再按图调整）
"""

# %%
INPUT_GEF = "path/to/your_sample.gef"
BIN_TYPE = "bins"
BIN_SIZE = 100
SAMPLE_NAME = "sample_01"

MIN_GENE = 200
MAX_GENE = 8000
MIN_COUNT = 500
MAX_PCT_MT = 20.0
MIN_CELL_BY_GENE = 3

N_TOP_HVG = 3000
N_PCS = 50
N_NEIGHBORS = 15
N_PCS_NEIGHBOR = 30
LEIDEN_RES = 0.6

OUT_DIR = Path("results_single")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# %% [markdown]
"""
## 2) 读取单样本数据

如果你是 GEM 格式，请替换为对应读取 API（如 `st.io.read_gem`）。
"""

# %%
data = st.io.read_gef(file_path=INPUT_GEF, bin_type=BIN_TYPE, bin_size=BIN_SIZE)
data.sn = SAMPLE_NAME

# 若默认读到的是 ENSG，尝试直接从 GEF 注释表提取 gene symbol 并替换
if hasattr(data, "gene_names"):
    current_genes = pd.Index([str(x) for x in data.gene_names])
else:
    current_genes = pd.Index([])

if is_ensg_like(current_genes):
    gene_symbols, hit_table, hit_col = extract_gene_symbols_from_gef(INPUT_GEF)
    if gene_symbols is not None and len(gene_symbols) == len(current_genes):
        data.gene_names = np.asarray(gene_symbols)
        print("✅ 已从 GEF 注释表直接替换为 gene symbol。")
        print(f"   命中 table path: {hit_table}")
        print(f"   命中 column name: {hit_col}")
    else:
        print("⚠️ 未能从 GEF 直接提取可对齐的 gene symbol，将继续使用默认基因名。")

print(data)


# %% [markdown]
"""
## 2.1) 诊断：打印 GEF 实际命中的 table path + column name

这个小 cell 可以单独运行，快速确认你当前 GEF 用的是哪张基因注释表、哪一列作为 symbol。
"""

# %%
_symbols, _table, _col = extract_gene_symbols_from_gef(INPUT_GEF)
if _table is None:
    print("未命中候选注释表/字段（可能该 GEF 结构不同，需扩展候选路径）。")
else:
    print(f"GEF 命中 table path: {_table}")
    print(f"GEF 命中 column name: {_col}")
    print(f"提取到基因数: {len(_symbols)}")

# %% [markdown]
"""
## 3) 解释为什么你会出现 `线粒体比例 = 0`
常见原因：
1. **基因名不匹配**：人常见 `MT-`，鼠常见 `mt-`；某些参考基因组是 `Mt-` 或 Ensembl ID（没有 `MT-` 前缀）。
2. **未显式标注 mt 基因**：某些版本 `cal_qc()` 不会自动识别全部命名风格。
3. **基因符号列读取错位**：你看的并不是 symbol，而是 feature id。

下面先显式识别 mt 基因，再计算一个 `pct_counts_mt_manual` 进行交叉验证。
"""

# %%
# 尝试从 StereoExpData 提取基因名
if hasattr(data, "gene_names"):
    gene_names = pd.Index([str(x) for x in data.gene_names])
elif hasattr(data, "genes") and hasattr(data.genes, "index"):
    gene_names = pd.Index([str(x) for x in data.genes.index])
else:
    raise ValueError("无法在 data 中找到基因名，请检查 StereoExpData 结构。")

# 常见线粒体基因命名规则（人/鼠/大小写混用）
mt_regex = re.compile(r"^(MT-|mt-|Mt-)")
mt_mask = gene_names.to_series().str.match(mt_regex, na=False).values

print(f"检测到 mt 基因数量: {int(mt_mask.sum())}")
if mt_mask.sum() == 0:
    print("⚠️ 未匹配到任何 mt 前缀基因。你很可能使用的是 Ensembl ID 或其他命名体系。")
    print("   建议：先做 ID->symbol 映射，再计算线粒体占比。")


# %%
# 先跑官方常规 QC
# 你的版本若支持 mt 相关参数（例如 mt_prefix / mito_genes），建议显式传入。
data.tl.cal_qc()

# 手动计算 mt 占比，避免因命名识别失败导致全 0
X = data.exp_matrix
if hasattr(X, "toarray"):
    total_counts = np.asarray(X.sum(axis=1)).reshape(-1)
    mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).reshape(-1) if mt_mask.sum() > 0 else np.zeros(X.shape[0])
else:
    total_counts = X.sum(axis=1)
    mt_counts = X[:, mt_mask].sum(axis=1) if mt_mask.sum() > 0 else np.zeros(X.shape[0])

pct_mt_manual = np.divide(mt_counts, np.maximum(total_counts, 1), dtype=float) * 100

# 写入 cells 元数据，供过滤与可视化
if hasattr(data, "cells") and isinstance(data.cells, pd.DataFrame):
    data.cells["pct_counts_mt_manual"] = pct_mt_manual
else:
    print("⚠️ data.cells 不是 DataFrame，无法直接写入 pct_counts_mt_manual，请按你版本对象结构调整。")

print(pd.Series(pct_mt_manual).describe(percentiles=[0.5, 0.9, 0.95, 0.99]))


# %% [markdown]
"""
## 4) 可视化 QC 指标并定阈值
你可以重点看：
- n_genes 分布是否有低质量长尾
- total_counts 是否有极低 UMI 与潜在 doublet 高尾
- pct_counts_mt_manual 的 95/99 分位
"""

# %%
# 以下绘图函数名在不同版本可能略有差异，请以 help(data.plt) 为准。
# 常见是 violin/scatter 的组合。
if hasattr(data, "plt") and hasattr(data.plt, "violin"):
    data.plt.violin(keys=["n_genes_by_counts", "total_counts", "pct_counts_mt_manual"], group_by=None)


# %% [markdown]
"""
## 5) 过滤（建议优先用手动 mt 比例字段）
"""

# %%
# 先按基因数和 counts 过滤
data.tl.filter_cells(
    min_gene=MIN_GENE,
    max_gene=MAX_GENE,
    min_count=MIN_COUNT,
    inplace=True,
)

# 再按手动 mt 比例过滤（避免官方自动 mt=0 的问题）
if hasattr(data, "cells") and "pct_counts_mt_manual" in data.cells.columns:
    keep_mask = data.cells["pct_counts_mt_manual"].values <= MAX_PCT_MT
    data = data.sub_by_index(cell_index=np.where(keep_mask)[0])

# 基因过滤
data.tl.filter_genes(min_cell=MIN_CELL_BY_GENE, inplace=True)

# 保存 raw
data.tl.raw_checkpoint()
print(f"After filtering: cells={data.n_cells}, genes={data.n_genes}")


# %% [markdown]
"""
## 6) 标准化 -> HVG -> PCA -> 邻居图 -> UMAP -> Leiden
"""

# %%
data.tl.normalize_total(target_sum=1e4)
data.tl.log1p()
data.tl.highly_variable_genes(n_top_genes=N_TOP_HVG)
data.tl.scale(max_value=10)
data.tl.pca(n_pcs=N_PCS)
data.tl.neighbors(n_neighbors=N_NEIGHBORS, n_pcs=N_PCS_NEIGHBOR)
data.tl.umap(min_dist=0.3)
data.tl.leiden(resolution=LEIDEN_RES, key_added="leiden_r06")

if hasattr(data, "plt") and hasattr(data.plt, "umap"):
    data.plt.umap(color=["leiden_r06"])


# %% [markdown]
"""
## 7) marker 基因与保存
"""

# %%
data.tl.rank_genes_groups(groupby="leiden_r06", method="wilcoxon")
if hasattr(data, "plt") and hasattr(data.plt, "rank_genes_groups"):
    data.plt.rank_genes_groups(n_genes=15, sharey=False)

rgg = data.tl.result["rank_genes_groups"]
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
marker_df.to_csv(OUT_DIR / "markers_leiden_r06_top20.csv", index=False)

data.write_h5ad(str(OUT_DIR / "single_sample_processed.h5ad"))
print(f"Saved to: {OUT_DIR}")
