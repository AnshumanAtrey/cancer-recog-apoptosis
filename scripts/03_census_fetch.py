#!/usr/bin/env python3
"""
Step 2a — CELLxGENE Census data layer (EXPLORE + FETCH).

Step 2 goal: a ranked shortlist of receptors enriched on CANCER cells vs matched
NORMAL cells, that participate in cell-cell communication — feeds Step 3 (specificity
audit) and Step 4 (reward: bind cancer target, avoid normal homolog).

This script is the DATA LAYER only. Analysis (LIANA+ communication inference) is
scripts/04, designed AFTER we see what `explore` reports — do not assume disease
strings or cell-type labels; look at the data first.

WHY LIANA+ (not CellChat-in-R), decided in scripts/04: LIANA+ natively reimplements
the CellChat method in Python (+ CellPhoneDB/NATMI/Connectome/log2FC + a consensus
rank-aggregate), in the scanpy/AnnData ecosystem — faithful to the plan's CellChat
intent, multi-method robust, and no R/Bioconductor on Colab.

TWO MODES (the notebook runs EXPLORE first, then FETCH):
  explore  — read obs metadata only (cheap). For each candidate tissue, print the
             available `disease` values + cell counts, the top `cell_type` values,
             and whether 'malignant cell' / 'epithelial cell' are annotated. This
             tells us the EXACT filter strings and whether malignant cells are
             labelled (else Step 2b needs CNV inference).
  fetch    — for each target, pull a capped, subsampled AnnData of tumour cells and
             of normal cells, set var_names to gene symbols, save .h5ad. Idempotent.

REQUIREMENTS:
  pip install cellxgene-census scanpy
  CPU is fine; high-RAM Colab helps for large pulls. No GPU needed.

USAGE:
  python scripts/03_census_fetch.py explore     # do this FIRST
  python scripts/03_census_fetch.py fetch       # after confirming disease strings
  python scripts/03_census_fetch.py             # explore then fetch
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "cellxgene"

# Pin the Census release for reproducibility (a dated LTS snapshot).
CENSUS_VERSION = "2024-07-01"

# tissue_general values to explore (CL/UBERON-standardised in Census).
CANDIDATE_TISSUES = ["lung", "breast", "colon"]

# Fetch targets. `cancer_disease` strings are BEST-GUESSES — the explore step prints
# the exact values; fix these if a target reports 0 cells. normal uses disease=='normal'.
TARGETS = [
    {"label": "lung",  "tissue_general": "lung",   "cancer_disease": "lung adenocarcinoma"},
    {"label": "breast", "tissue_general": "breast", "cancer_disease": "breast carcinoma"},
    {"label": "colon", "tissue_general": "colon",  "cancer_disease": "colorectal cancer"},
]

# Memory guard: random-subsample each condition to at most this many cells.
MAX_CELLS_PER_CONDITION = 20000
# Only well-behaved single-cell data.
BASE_FILTER = "is_primary_data == True"
SUBSAMPLE_SEED = 20260529

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", stream=sys.stdout)
log = logging.getLogger("step2a")


def open_human(census):
    return census["census_data"]["homo_sapiens"]


def read_obs(human, value_filter: str, columns: list[str]):
    """Read obs metadata only (no expression) as a pandas DataFrame."""
    return human.obs.read(value_filter=value_filter, column_names=columns).concat().to_pandas()


# ---------------- EXPLORE ----------------
def explore(census) -> None:
    human = open_human(census)
    log.info("Census version=%s", CENSUS_VERSION)
    for tissue in CANDIDATE_TISSUES:
        log.info("=" * 60)
        log.info("[explore] tissue_general == '%s'", tissue)
        try:
            df = read_obs(human,
                          f"tissue_general == '{tissue}' and {BASE_FILTER}",
                          ["disease", "cell_type", "assay"])
        except Exception as e:
            log.error("[explore] obs read failed for %s: %s: %s", tissue, type(e).__name__, e)
            continue
        log.info("[%s] total primary cells = %d", tissue, len(df))

        dz = df["disease"].value_counts()
        log.info("[%s] disease values (top 12):", tissue)
        for name, n in dz.head(12).items():
            log.info("      %8d  %s", n, name)

        ct = df["cell_type"].value_counts()
        log.info("[%s] cell_type values (top 12):", tissue)
        for name, n in ct.head(12).items():
            log.info("      %8d  %s", n, name)

        for marker in ("malignant cell", "neoplastic cell", "epithelial cell"):
            present = marker in ct.index
            log.info("[%s] cell_type '%s' annotated? %s%s", tissue, marker, present,
                     f" ({ct[marker]} cells)" if present else "")


# ---------------- FETCH ----------------
def _subsampled_coords(human, value_filter: str, cap: int):
    """Return soma_joinids for a filter, randomly capped to `cap` (or None = all)."""
    import numpy as np
    ids = human.obs.read(value_filter=value_filter, column_names=["soma_joinid"]).concat() \
        .to_pandas()["soma_joinid"].to_numpy()
    n = len(ids)
    if n == 0:
        return None, 0
    if n <= cap:
        return ids, n
    rng = np.random.default_rng(SUBSAMPLE_SEED)
    return np.sort(rng.choice(ids, size=cap, replace=False)), n


def fetch_condition(census, label: str, condition: str, value_filter: str, cap: int) -> bool:
    """Pull a capped AnnData for one (target, condition) and save .h5ad. Idempotent."""
    import cellxgene_census
    out = DATA_DIR / f"{label}__{condition}.h5ad"
    if out.exists():
        log.info("[%s/%s] SKIP — %s exists", label, condition, out.name)
        return True

    human = open_human(census)
    coords, total = _subsampled_coords(human, value_filter, cap)
    if coords is None or total == 0:
        log.warning("[%s/%s] 0 cells for filter: %s  → fix disease string (see explore)",
                    label, condition, value_filter)
        return False
    log.info("[%s/%s] %d cells match; pulling %d (cap=%d)", label, condition, total, len(coords), cap)

    adata = cellxgene_census.get_anndata(
        census=census, organism="Homo sapiens",
        obs_coords=[int(i) for i in coords],
        obs_column_names=["disease", "tissue_general", "cell_type", "assay", "donor_id", "dataset_id"],
    )
    # LIANA needs gene SYMBOLS as var_names; Census var has feature_name (symbol).
    if "feature_name" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names
        adata.var_names = adata.var["feature_name"].astype(str)
        adata.var_names_make_unique()
    adata.obs["condition"] = condition
    adata.obs["target_label"] = label

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    log.info("[%s/%s] saved → %s  shape=%s  cell_types=%d",
             label, condition, out.relative_to(PROJECT_ROOT), adata.shape, adata.obs["cell_type"].nunique())
    return True


def fetch(census) -> int:
    ok = 0
    for t in TARGETS:
        label, tissue, cancer = t["label"], t["tissue_general"], t["cancer_disease"]
        log.info("=" * 60)
        log.info("[fetch] target=%s tissue=%s cancer_disease='%s'", label, tissue, cancer)
        tumour_ok = fetch_condition(
            census, label, "tumour",
            f"tissue_general == '{tissue}' and disease == '{cancer}' and {BASE_FILTER}",
            MAX_CELLS_PER_CONDITION)
        normal_ok = fetch_condition(
            census, label, "normal",
            f"tissue_general == '{tissue}' and disease == 'normal' and {BASE_FILTER}",
            MAX_CELLS_PER_CONDITION)
        if tumour_ok and normal_ok:
            ok += 1
    log.info("=" * 60)
    log.info("[fetch] %d/%d targets fully fetched into %s", ok, len(TARGETS), DATA_DIR.relative_to(PROJECT_ROOT))
    return ok


def main() -> int:
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode not in ("explore", "fetch", "all"):
        log.error("usage: python scripts/03_census_fetch.py [explore|fetch|all]"); return 2
    try:
        import cellxgene_census
    except ImportError:
        log.error("missing dependency — run: pip install -q cellxgene-census scanpy"); return 2

    log.info("cancer-recon-apoptosis — Step 2a — CELLxGENE Census data layer (mode=%s)", mode)
    try:
        census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
    except Exception as e:
        log.error("could not open Census %s: %s: %s", CENSUS_VERSION, type(e).__name__, e); return 3
    try:
        if mode in ("explore", "all"):
            explore(census)
        if mode in ("fetch", "all"):
            fetch(census)
    finally:
        census.close()
    log.info("✓ done (mode=%s)", mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
