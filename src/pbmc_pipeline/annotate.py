"""Marker-based cell type annotation with a recorded confidence for every cluster."""

import pandas as pd
import scanpy as sc

from .config import Config


def score_marker_sets(adata: sc.AnnData, cfg: Config) -> list[str]:
    """Add one per-cell score column per cell type (mean marker expression minus a matched control set)."""
    score_cols = []
    for cell_type, genes in cfg.marker_sets.items():
        present = [g for g in genes if g in adata.raw.var_names]
        col = "score_" + cell_type.replace(" ", "_")
        sc.tl.score_genes(adata, present, score_name=col, use_raw=True, random_state=cfg.random_state)
        score_cols.append(col)
    return score_cols


def annotate_clusters(adata: sc.AnnData, cfg: Config, groupby: str = "leiden") -> pd.DataFrame:
    """Label each cluster with the cell type whose marker score is highest on average.

    Confidence is "high" when the winning score beats the runner-up by at least
    cfg.min_annotation_margin, otherwise "review" so a person checks the markers.
    """
    score_cols = score_marker_sets(adata, cfg)
    means = adata.obs.groupby(groupby, observed=True)[score_cols].mean()
    means.columns = [c.removeprefix("score_").replace("_", " ") for c in score_cols]

    rows = []
    for cluster, scores in means.iterrows():
        ranked = scores.sort_values(ascending=False)
        margin = ranked.iloc[0] - ranked.iloc[1]
        rows.append({
            groupby: cluster,
            "n_cells": int((adata.obs[groupby] == cluster).sum()),
            "cell_type": ranked.index[0],
            "top_score": round(ranked.iloc[0], 3),
            "runner_up": ranked.index[1],
            "runner_up_score": round(ranked.iloc[1], 3),
            "margin": round(margin, 3),
            "confidence": "high" if margin >= cfg.min_annotation_margin else "review",
        })
    table = pd.DataFrame(rows)
    mapping = dict(zip(table[groupby], table["cell_type"]))
    adata.obs["cell_type"] = adata.obs[groupby].map(mapping).astype("category")
    return table


def review_evidence(adata: sc.AnnData, table: pd.DataFrame, cfg: Config, groupby: str = "leiden") -> pd.DataFrame:
    """For each low-confidence cluster, compare it to the rest of the data on the marker genes of
    the two competing cell types, so a reviewer can confirm or overturn the call."""
    rows = []
    expr = adata.raw.to_adata()
    for _, r in table[table["confidence"] == "review"].iterrows():
        in_cl = (adata.obs[groupby] == r[groupby]).to_numpy()
        for candidate in (r["cell_type"], r["runner_up"]):
            for gene in cfg.marker_sets[candidate]:
                if gene not in expr.var_names:
                    continue
                x = expr[:, gene].X
                x = x.toarray().ravel() if hasattr(x, "toarray") else x.ravel()
                rows.append({
                    groupby: r[groupby], "call": r["cell_type"], "candidate": candidate, "gene": gene,
                    "mean_expr_in_cluster": round(float(x[in_cl].mean()), 2),
                    "mean_expr_elsewhere": round(float(x[~in_cl].mean()), 2),
                    "pct_cells_expressing": round(float((x[in_cl] > 0).mean() * 100), 1),
                })
    return pd.DataFrame(rows)
