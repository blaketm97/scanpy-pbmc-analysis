"""Run the full PBMC 3k analysis end to end.

Usage:
    python run_pipeline.py                 # default resolution from config
    python run_pipeline.py --resolution 0.8
"""

import argparse
import json
import platform
import sys
import warnings
from importlib.metadata import version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
warnings.filterwarnings("ignore", category=FutureWarning)

import scanpy as sc  # noqa: E402

from pbmc_pipeline import analysis, annotate, plots, qc  # noqa: E402
from pbmc_pipeline.config import Config  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resolution", type=float, default=None, help="Leiden resolution (default from config)")
    args = parser.parse_args()

    cfg = Config()
    if args.resolution is not None:
        cfg.leiden_resolution = args.resolution
    for d in (cfg.figures_dir, cfg.tables_dir, cfg.processed_h5ad.parent):
        d.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    fig, tab = cfg.figures_dir, cfg.tables_dir

    print("1/6  Loading counts and computing QC metrics")
    adata = qc.add_qc_metrics(qc.load_counts(cfg))
    adata.obs = adata.obs.join(qc.flag_cells(adata, cfg))
    plots.qc_distributions(adata, cfg, fig / "01_qc_distributions.png")
    plots.qc_scatter(adata, cfg, fig / "02_qc_scatter.png")
    adata.obs = adata.obs.drop(columns=["too_few_genes", "too_many_genes", "high_mito", "fail_fixed_qc", "mad_outlier"])

    print("2/6  Filtering cells and genes, removing doublets")
    adata, qc_log = qc.run_qc(adata, cfg)
    qc_log.to_csv(tab / "qc_summary.csv", index=False)

    print("3/6  Normalizing, selecting variable genes, PCA, UMAP")
    adata = analysis.normalize_and_select(adata, cfg)
    adata = analysis.embed(adata, cfg)
    plots.pca_variance(adata, cfg, fig / "03_pca_variance.png")

    print("4/6  Clustering and testing stability across resolutions and seeds")
    stab = analysis.cluster_with_stability(adata, cfg)
    stab.to_csv(tab / "cluster_stability.csv", index=False)
    plots.stability(stab, cfg.leiden_resolution, fig / "04_cluster_stability.png")

    print("5/6  Finding marker genes and annotating cell types")
    markers = analysis.find_markers(adata, cfg)
    markers.to_csv(tab / "cluster_markers_top10.csv", index=False)
    ann = annotate.annotate_clusters(adata, cfg)
    ann.to_csv(tab / "cluster_annotation.csv", index=False)
    annotate.review_evidence(adata, ann, cfg).to_csv(tab / "annotation_review_evidence.csv", index=False)

    order = [ct for ct in cfg.marker_sets if ct in set(adata.obs["cell_type"])]
    colors = plots.cell_type_colors(order)
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(order)
    counts = adata.obs["cell_type"].value_counts()[order]
    (counts.rename("n_cells").to_frame()
        .assign(pct_of_cells=lambda d: (d["n_cells"] / d["n_cells"].sum() * 100).round(1))
        .rename_axis("cell_type").to_csv(tab / "cell_type_counts.csv"))

    print("6/6  Drawing figures and saving results")
    plots.umap_cell_types(adata, colors, fig / "05_umap_cell_types.png")
    plots.marker_dotplot(adata, cfg, order, fig / "06_marker_dotplot.png")
    plots.marker_umaps(adata, ["IL7R", "CCR7", "S100A4", "CD8A", "NKG7", "MS4A1", "LYZ", "FCGR3A",
                               "CD14", "FCER1A", "PPBP", "GNLY"], fig / "07_marker_umaps.png")
    plots.composition(adata, colors, fig / "08_cell_type_composition.png")

    adata.uns["cell_type_colors"] = [colors[c] for c in order]
    adata.write_h5ad(cfg.processed_h5ad, compression="gzip")

    summary = {
        "cells_retained": int(adata.n_obs),
        "genes_retained": int(adata.raw.n_vars),
        "highly_variable_genes": int(adata.n_vars),
        "leiden_resolution": cfg.leiden_resolution,
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "n_cell_types": int(adata.obs["cell_type"].nunique()),
        "clusters_needing_review": ann.loc[ann["confidence"] == "review", "leiden"].tolist(),
        "versions": {p: version(p) for p in ["scanpy", "anndata", "numpy", "pandas", "scikit-learn",
                                             "leidenalg", "igraph", "matplotlib"]},
        "python": platform.python_version(),
    }
    (tab / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "versions"}, indent=2))


if __name__ == "__main__":
    main()
