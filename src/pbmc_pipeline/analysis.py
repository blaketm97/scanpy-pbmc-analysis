"""Normalization, feature selection, embedding, clustering, cluster stability and marker genes."""

import itertools

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from .config import Config


def normalize_and_select(adata: sc.AnnData, cfg: Config) -> sc.AnnData:
    """Library-size normalize, log-transform, and pick highly variable genes."""
    sc.pp.normalize_total(adata, target_sum=cfg.target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.n_top_genes, flavor="seurat")
    adata.raw = adata  # keep all log-normalized genes for marker testing and plots
    return adata


def embed(adata: sc.AnnData, cfg: Config) -> sc.AnnData:
    """Scale HVGs, run PCA, build the kNN graph and compute UMAP."""
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=cfg.n_pcs, random_state=cfg.random_state)
    sc.pp.neighbors(
        adata, n_neighbors=cfg.n_neighbors, n_pcs=cfg.n_pcs_neighbors, random_state=cfg.random_state
    )
    sc.tl.umap(adata, random_state=cfg.random_state)
    return adata


def leiden(adata: sc.AnnData, resolution: float, seed: int, key: str) -> None:
    sc.tl.leiden(
        adata, resolution=resolution, random_state=seed, key_added=key,
        flavor="igraph", n_iterations=2, directed=False,
    )


def cluster_with_stability(adata: sc.AnnData, cfg: Config) -> pd.DataFrame:
    """Cluster at several resolutions and measure how much labels change across random seeds.

    Mean pairwise adjusted Rand index (ARI) across seeds near 1.0 means the clustering
    is reproducible at that resolution rather than an artifact of one random start.
    """
    rows = []
    for res in cfg.resolutions_to_test:
        labelings = []
        for seed in cfg.stability_seeds:
            key = f"_leiden_{res}_{seed}"
            leiden(adata, res, seed, key)
            labelings.append(adata.obs[key].to_numpy())
            del adata.obs[key]
        aris = [adjusted_rand_score(a, b) for a, b in itertools.combinations(labelings, 2)]
        rows.append({
            "resolution": res,
            "n_clusters_median": int(np.median([len(set(l)) for l in labelings])),
            "mean_pairwise_ARI": round(float(np.mean(aris)), 3),
            "min_pairwise_ARI": round(float(np.min(aris)), 3),
        })
    leiden(adata, cfg.leiden_resolution, cfg.random_state, "leiden")
    return pd.DataFrame(rows)


def find_markers(adata: sc.AnnData, cfg: Config, groupby: str = "leiden") -> pd.DataFrame:
    """Wilcoxon rank-sum test of each cluster vs. all others on log-normalized expression."""
    sc.tl.rank_genes_groups(adata, groupby, method="wilcoxon", use_raw=True, pts=True)
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0)]
    top = df.groupby("group", observed=True).head(cfg.n_marker_genes).copy()
    top = top.rename(columns={"group": groupby, "names": "gene", "pct_nz_group": "pct_in_cluster",
                              "pct_nz_reference": "pct_in_other_cells"})
    cols = [groupby, "gene", "logfoldchanges", "pvals_adj", "pct_in_cluster", "pct_in_other_cells"]
    return top[cols].round({"logfoldchanges": 2, "pct_in_cluster": 3, "pct_in_other_cells": 3})
