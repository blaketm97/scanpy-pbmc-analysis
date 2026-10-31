"""Analysis parameters in one place, so every threshold is visible and easy to change."""

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


@dataclass
class Config:
    # Paths
    raw_dir: Path = PROJECT_ROOT / "data" / "raw" / "filtered_gene_bc_matrices" / "hg19"
    processed_h5ad: Path = PROJECT_ROOT / "data" / "processed" / "pbmc3k_processed.h5ad"
    figures_dir: Path = PROJECT_ROOT / "results" / "figures"
    tables_dir: Path = PROJECT_ROOT / "results" / "tables"

    # Quality control (standard thresholds for 10x v1 PBMC data)
    min_genes_per_cell: int = 200      # drop empty droplets / debris
    max_genes_per_cell: int = 2500     # drop likely multiplets
    max_pct_mito: float = 5.0          # drop dying / damaged cells
    min_cells_per_gene: int = 3        # drop genes seen in almost no cells
    mad_cutoff: float = 5.0            # data-driven outlier flag, reported alongside fixed thresholds

    # Doublet detection (Scrublet)
    expected_doublet_rate: float = 0.022   # 10x rule of thumb: ~0.8% per 1,000 cells recovered (~2.7k here)

    # Normalization / feature selection
    target_sum: float = 1e4
    n_top_genes: int = 2000

    # Dimensionality reduction and clustering
    n_pcs: int = 50
    n_pcs_neighbors: int = 30
    n_neighbors: int = 15
    leiden_resolution: float = 1.0   # chosen from results/figures/04_cluster_stability.png; see README
    resolutions_to_test: tuple = (0.2, 0.4, 0.6, 0.8, 1.0, 1.4)
    stability_seeds: tuple = (0, 1, 2, 3, 4)
    random_state: int = 0

    # Marker detection
    n_marker_genes: int = 10
    min_annotation_margin: float = 0.25  # score gap between top two cell types needed for a "high" confidence call

    # Canonical PBMC marker genes used to label clusters (from published PBMC atlases / 10x guidance)
    marker_sets: dict = field(default_factory=lambda: {
        "Naive CD4 T": ["IL7R", "CCR7", "LEF1", "TCF7", "CD4"],
        "Memory CD4 T": ["IL7R", "S100A4", "AQP3", "CD40LG", "CD4"],
        "CD8 T": ["CD8A", "CD8B", "CCL5", "GZMK"],
        "NK": ["GNLY", "NKG7", "GZMB", "PRF1", "KLRD1"],
        "B": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "CD14 Monocyte": ["CD14", "LYZ", "S100A8", "S100A9"],
        "FCGR3A Monocyte": ["FCGR3A", "MS4A7", "LST1", "IFITM3"],
        "Dendritic": ["FCER1A", "CST3", "CLEC10A", "CD1C"],
        "Platelet": ["PPBP", "PF4", "GP9"],
    })
