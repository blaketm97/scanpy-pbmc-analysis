"""Figures for QC, clustering and annotation. All saved as PNG to results/figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Colorblind-checked categorical order; identity is also shown with direct labels on every plot.
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
               "#008300", "#6250d6", "#e34948"]
NEUTRAL = "#8c8b85"
PASS_COLOR, FAIL_COLOR = "#2a78d6", "#eb6834"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": False, "grid.color": GRID,
    "legend.frameon": False,
})


def cell_type_colors(cell_types: list[str]) -> dict:
    """Fixed color per cell type (in the order given); any beyond eight get a neutral gray."""
    return {ct: (CATEGORICAL[i] if i < len(CATEGORICAL) else NEUTRAL) for i, ct in enumerate(cell_types)}


def qc_distributions(adata, cfg, out: Path) -> None:
    metrics = [
        ("n_genes_by_counts", "Genes detected per cell", [cfg.min_genes_per_cell, cfg.max_genes_per_cell]),
        ("total_counts", "UMI counts per cell", []),
        ("pct_counts_mt", "% mitochondrial reads", [cfg.max_pct_mito]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    rng = np.random.default_rng(0)
    for ax, (col, title, cuts) in zip(axes, metrics):
        vals = adata.obs[col].to_numpy()
        fail = adata.obs["fail_fixed_qc"].to_numpy()
        parts = ax.violinplot(vals, showextrema=False, widths=0.8)
        for b in parts["bodies"]:
            b.set_facecolor(GRID); b.set_edgecolor(INK_2); b.set_alpha(1)
        jitter = rng.uniform(-0.25, 0.25, len(vals)) + 1
        ax.scatter(jitter[~fail], vals[~fail], s=2, color=PASS_COLOR, alpha=0.5, linewidths=0)
        ax.scatter(jitter[fail], vals[fail], s=8, color=FAIL_COLOR, linewidths=0)
        for c in cuts:
            ax.axhline(c, color=INK_2, lw=1, ls="--")
            ax.text(0.58, c, f"{c:g} ", va="bottom", ha="left", fontsize=8, color=INK_2)
        ax.set_title(title); ax.set_xticks([])
    axes[0].scatter([], [], color=PASS_COLOR, s=12, label="Kept")
    axes[0].scatter([], [], color=FAIL_COLOR, s=12, label="Removed by QC")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle(f"Per-cell quality metrics before filtering (n = {adata.n_obs:,} cells)", x=0.01, y=1.06, ha="left")
    fig.savefig(out); plt.close(fig)


def qc_scatter(adata, cfg, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    fail = adata.obs["fail_fixed_qc"].to_numpy()
    for ax, (x, y, xl, yl) in zip(axes, [
        ("total_counts", "n_genes_by_counts", "UMI counts per cell", "Genes detected"),
        ("total_counts", "pct_counts_mt", "UMI counts per cell", "% mitochondrial reads"),
    ]):
        ax.scatter(adata.obs.loc[~fail, x], adata.obs.loc[~fail, y], s=4, color=PASS_COLOR, alpha=0.5, linewidths=0, label="Kept")
        ax.scatter(adata.obs.loc[fail, x], adata.obs.loc[fail, y], s=10, color=FAIL_COLOR, linewidths=0, label="Removed")
        ax.set_xlabel(xl); ax.set_ylabel(yl)
    axes[0].axhline(cfg.max_genes_per_cell, color=INK_2, lw=1, ls="--")
    axes[0].axhline(cfg.min_genes_per_cell, color=INK_2, lw=1, ls="--")
    axes[1].axhline(cfg.max_pct_mito, color=INK_2, lw=1, ls="--")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("QC thresholds against the joint distribution", x=0.01, y=1.03, ha="left")
    fig.savefig(out); plt.close(fig)


def pca_variance(adata, cfg, out: Path) -> None:
    ratio = adata.uns["pca"]["variance_ratio"] * 100
    fig, ax = plt.subplots(figsize=(6, 3.4))
    pcs = np.arange(1, len(ratio) + 1)
    ax.plot(pcs, ratio, color=PASS_COLOR, lw=2, marker="o", ms=3.5)
    ax.axvline(cfg.n_pcs_neighbors, color=INK_2, lw=1, ls="--")
    ax.text(cfg.n_pcs_neighbors + 0.5, ratio.max() * 0.9, f"{cfg.n_pcs_neighbors} PCs used", fontsize=8, color=INK_2)
    ax.set_xlabel("Principal component"); ax.set_ylabel("% variance explained")
    ax.set_title("PCA elbow plot", loc="left")
    fig.savefig(out); plt.close(fig)


def stability(stab: pd.DataFrame, chosen: float, out: Path) -> None:
    """Two panels on separate axes: reproducibility (ARI) and number of clusters, by resolution."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    axes[0].plot(stab["resolution"], stab["mean_pairwise_ARI"], color=PASS_COLOR, lw=2, marker="o", ms=5)
    axes[0].fill_between(stab["resolution"], stab["min_pairwise_ARI"], stab["mean_pairwise_ARI"],
                         color=PASS_COLOR, alpha=0.15, lw=0)
    axes[0].set_ylim(min(0.6, stab["min_pairwise_ARI"].min() - 0.05), 1.01)
    axes[0].set_ylabel("Agreement across 5 seeds (ARI)")
    axes[0].set_title("Clustering reproducibility", loc="left")
    axes[1].plot(stab["resolution"], stab["n_clusters_median"], color=PASS_COLOR, lw=2, marker="o", ms=5)
    axes[1].set_ylabel("Number of clusters"); axes[1].set_title("Clusters found", loc="left")
    for ax in axes:
        ax.axvline(chosen, color=INK_2, lw=1, ls="--")
        ax.set_xlabel("Leiden resolution")
    axes[1].text(chosen, axes[1].get_ylim()[1], f" chosen: {chosen}", fontsize=8, color=INK_2, va="top")
    fig.savefig(out); plt.close(fig)


def umap_cell_types(adata, colors: dict, out: Path) -> None:
    """UMAP colored by cell type; each cluster's ID is written at its center, and the legend maps IDs to types."""
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    xy = adata.obsm["X_umap"]
    for ct, col in colors.items():
        m = (adata.obs["cell_type"] == ct).to_numpy()
        clusters = sorted(adata.obs.loc[m, "leiden"].unique(), key=int)
        ax.scatter(xy[m, 0], xy[m, 1], s=5, color=col, linewidths=0, alpha=0.85,
                   label=f"{ct}  (cl. {', '.join(clusters)}; {m.sum():,} cells)")
    for cl in adata.obs["leiden"].cat.categories:
        m = (adata.obs["leiden"] == cl).to_numpy()
        cx, cy = np.median(xy[m], axis=0)
        ax.text(cx, cy, cl, ha="center", va="center", fontsize=9, fontweight="bold", color=INK,
                bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec=INK_2, lw=0.6, alpha=0.9))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), markerscale=3, fontsize=9,
              title="Cell type (Leiden cluster)", title_fontsize=9, alignment="left")
    ax.set_title(f"PBMC 3k: {adata.n_obs:,} cells, {adata.obs['leiden'].nunique()} clusters, "
                 f"{adata.obs['cell_type'].nunique()} cell types", loc="left")
    fig.savefig(out); plt.close(fig)


def marker_dotplot(adata, cfg, order: list[str], out: Path) -> None:
    genes, used = {}, set()
    for ct in order:  # up to three markers per type, never repeating a gene across types
        picks = [g for g in cfg.marker_sets[ct] if g in adata.raw.var_names and g not in used][:3]
        used.update(picks)
        genes[ct] = picks
    dp = sc.pl.dotplot(adata, genes, groupby="cell_type", categories_order=order, use_raw=True,
                       standard_scale="var", cmap="Blues", return_fig=True)
    dp.savefig(out); plt.close("all")


def marker_umaps(adata, genes: list[str], out: Path) -> None:
    fig = sc.pl.umap(adata, color=genes, use_raw=True, cmap="Blues", ncols=4, frameon=False,
                     size=8, vmax="p99", return_fig=True, show=False)
    fig.savefig(out); plt.close(fig)


def composition(adata, colors: dict, out: Path) -> None:
    counts = adata.obs["cell_type"].value_counts()
    order = [ct for ct in colors if ct in counts.index]
    counts = counts[order][::-1]
    pct = counts / counts.sum() * 100
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.barh(counts.index, counts.values, color=[colors[c] for c in counts.index], height=0.7)
    for i, (n, p) in enumerate(zip(counts.values, pct.values)):
        ax.text(n + counts.max() * 0.01, i, f"{n:,} ({p:.1f}%)", va="center", fontsize=8, color=INK_2)
    ax.set_xlim(0, counts.max() * 1.22)
    ax.set_xlabel("Cells"); ax.set_title("Cell type composition", loc="left")
    fig.savefig(out); plt.close(fig)
