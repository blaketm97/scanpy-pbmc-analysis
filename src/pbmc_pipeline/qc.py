"""Loading the 10x count matrix and cell/gene quality control."""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

from .config import Config


def load_counts(cfg: Config) -> sc.AnnData:
    """Read the raw 10x matrix; fall back to Scanpy's built-in download of the same dataset."""
    if (cfg.raw_dir / "matrix.mtx").exists():
        adata = sc.read_10x_mtx(cfg.raw_dir, var_names="gene_symbols", cache=False)
    else:
        adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    return adata


def add_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """Annotate mitochondrial / ribosomal / hemoglobin genes and compute per-cell QC metrics."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.match(r"^HB(?!P)")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True
    )
    return adata


def mad_outlier(values: pd.Series, n_mads: float) -> pd.Series:
    """Flag values more than n_mads median absolute deviations from the median."""
    med = np.median(values)
    mad = median_abs_deviation(values)
    return (values < med - n_mads * mad) | (values > med + n_mads * mad)


def flag_cells(adata: sc.AnnData, cfg: Config) -> pd.DataFrame:
    """Return one boolean column per QC rule so every removal has a recorded reason."""
    obs = adata.obs
    flags = pd.DataFrame(index=obs.index)
    flags["too_few_genes"] = obs["n_genes_by_counts"] < cfg.min_genes_per_cell
    flags["too_many_genes"] = obs["n_genes_by_counts"] > cfg.max_genes_per_cell
    flags["high_mito"] = obs["pct_counts_mt"] > cfg.max_pct_mito
    flags["fail_fixed_qc"] = flags[["too_few_genes", "too_many_genes", "high_mito"]].any(axis=1)
    # Data-driven check: would a MAD rule have caught different cells?
    flags["mad_outlier"] = (
        mad_outlier(obs["log1p_total_counts"], cfg.mad_cutoff)
        | mad_outlier(obs["log1p_n_genes_by_counts"], cfg.mad_cutoff)
    )
    return flags


def detect_doublets(adata: sc.AnnData, cfg: Config) -> sc.AnnData:
    """Score likely doublets with Scrublet on raw counts."""
    sc.pp.scrublet(
        adata,
        expected_doublet_rate=cfg.expected_doublet_rate,
        random_state=cfg.random_state,
        verbose=False,
    )
    return adata


def run_qc(adata: sc.AnnData, cfg: Config) -> tuple[sc.AnnData, pd.DataFrame]:
    """Apply fixed-threshold QC, gene filtering and doublet removal; return filtered data and a log."""
    log = [("Cells loaded", adata.n_obs), ("Genes loaded", adata.n_vars)]

    flags = flag_cells(adata, cfg)
    adata.obs = adata.obs.join(flags)
    log += [
        ("Cells with < %d genes" % cfg.min_genes_per_cell, int(flags["too_few_genes"].sum())),
        ("Cells with > %d genes" % cfg.max_genes_per_cell, int(flags["too_many_genes"].sum())),
        ("Cells with > %.0f%% mitochondrial reads" % cfg.max_pct_mito, int(flags["high_mito"].sum())),
        ("Cells removed by fixed QC (any rule)", int(flags["fail_fixed_qc"].sum())),
        ("Cells flagged by %.0f-MAD rule (for comparison)" % cfg.mad_cutoff, int(flags["mad_outlier"].sum())),
    ]

    adata = adata[~adata.obs["fail_fixed_qc"]].copy()
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=cfg.min_cells_per_gene)
    log.append(("Genes removed (in < %d cells)" % cfg.min_cells_per_gene, n_genes_before - adata.n_vars))

    adata = detect_doublets(adata, cfg)
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    log.append(("Predicted doublets removed (Scrublet)", n_doublets))
    adata = adata[~adata.obs["predicted_doublet"]].copy()

    log += [("Cells retained", adata.n_obs), ("Genes retained", adata.n_vars)]
    return adata, pd.DataFrame(log, columns=["step", "count"])
