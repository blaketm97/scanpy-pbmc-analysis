"""Fast checks on small synthetic data: QC rules flag the right cells, and annotation
picks the cell type whose markers a cluster actually expresses."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pbmc_pipeline import annotate, qc  # noqa: E402
from pbmc_pipeline.config import Config  # noqa: E402


@pytest.fixture
def cfg():
    c = Config()
    c.min_genes_per_cell, c.max_genes_per_cell, c.max_pct_mito = 5, 40, 10.0
    return c


def make_adata(counts: np.ndarray, genes: list[str]) -> sc.AnnData:
    a = sc.AnnData(sparse.csr_matrix(counts.astype(np.float32)))
    a.var_names = genes
    a.obs_names = [f"cell{i}" for i in range(counts.shape[0])]
    return a


def test_qc_flags_each_rule(cfg):
    genes = ["MT-CO1"] + [f"G{i}" for i in range(59)]
    counts = np.zeros((4, 60))
    counts[0, 1:21] = 5                 # healthy: 20 genes, 0% mito
    counts[1, 1:3] = 5                  # too few genes
    counts[2, 1:51] = 1                 # too many genes
    counts[3, 1:21] = 5; counts[3, 0] = 50  # 33% mito
    a = qc.add_qc_metrics(make_adata(counts, genes))
    flags = qc.flag_cells(a, cfg)
    assert flags["fail_fixed_qc"].tolist() == [False, True, True, True]
    assert flags.loc["cell1", "too_few_genes"]
    assert flags.loc["cell2", "too_many_genes"]
    assert flags.loc["cell3", "high_mito"]


def test_mad_outlier_flags_extreme_value():
    s = pd.Series([10, 11, 9, 10, 12, 10, 100])
    assert qc.mad_outlier(s, 5).tolist() == [False] * 6 + [True]


def test_annotation_matches_expressed_markers(cfg):
    cfg.marker_sets = {"B": ["MS4A1", "CD79A"], "NK": ["GNLY", "NKG7"]}
    rng = np.random.default_rng(0)
    background = [f"G{i}" for i in range(200)]
    genes = ["MS4A1", "CD79A", "GNLY", "NKG7"] + background
    counts = rng.poisson(1, size=(100, len(genes))).astype(float)
    counts[:50, 0:2] += 30   # first half expresses B markers
    counts[50:, 2:4] += 30   # second half expresses NK markers
    a = make_adata(counts, genes)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    a.obs["leiden"] = pd.Categorical(["0"] * 50 + ["1"] * 50)
    table = annotate.annotate_clusters(a, cfg).set_index("leiden")
    assert table.loc["0", "cell_type"] == "B"
    assert table.loc["1", "cell_type"] == "NK"
    assert (table["confidence"] == "high").all()
