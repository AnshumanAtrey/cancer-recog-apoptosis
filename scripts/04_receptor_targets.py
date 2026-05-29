#!/usr/bin/env python3
"""
Step 2b — cancer-enriched RECEPTOR shortlist (differential expression).

Input : data/cellxgene/<tissue>__{tumour,normal}.h5ad  (from scripts/03)
Output: data/cellxgene/targets_shortlist.csv  + per-tissue CSVs
        ranked receptors over-expressed on CANCER cells vs matched NORMAL epithelium.

This is the load-bearing Step-2 deliverable (feeds Step 3 specificity audit + Step 4
reward). LIANA+ cell-cell *communication* annotation is a second pass (scripts/05);
here we use LIANA only for its curated receptor gene list.

ADAPTIVE CANCER-CELL SELECTION (the explore step showed malignant labels differ by
tissue): prefer explicit 'malignant cell' / 'neoplastic cell'; else fall back to
epithelial-lineage cell types (valid — these are carcinomas, malignant = aberrant
epithelium). Normal population = epithelial-lineage cells in the normal sample, so
the contrast is epithelium-vs-epithelium (controls for the heterogeneous normal pool).

METHOD: normalize (CP10k + log1p) → Wilcoxon rank_genes_groups(cancer vs normal) on the
receptor gene set → per tissue rank by log2FC (padj<0.05, expressed in >=10% of cancer
cells) → aggregate across tissues (a good PAN-CANCER target is enriched in several).

USAGE:  python scripts/04_receptor_targets.py
REQS :  scanpy, liana (CPU; no GPU)
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "cellxgene"

TISSUES = ["lung", "breast", "colon"]

EXPLICIT_CANCER_LABELS = ["malignant cell", "neoplastic cell"]
# NOTE: keep keywords lineage-specific. Do NOT use broad tissue words like "mammary"
# — that wrongly matches "fibroblast of mammary gland" (stroma). True epithelial types
# are caught by "epithel"/"luminal"/"basal cell" (incl. myoepithelial → has "epithel").
EPITHELIAL_KEYWORDS = [
    "epithel", "pneumocyte", "enterocyte", "colonocyte", "luminal", "basal cell",
    "goblet", "club cell", "secretory", "ductal", "acinar", "keratinocyte",
    "tuft", "paneth", "enteroendocrine", "ionocyte", "hillock", "serous",
]
MIN_CELLS = 50                 # need at least this many cells in each group
MIN_PCT_CANCER = 0.10          # receptor must be expressed in >=10% of cancer cells
LOG2FC_MIN = 0.5               # enrichment threshold for "enriched in this tissue"
PADJ_MAX = 0.05

# Step-1 reference receptors — report their ranks as a sanity anchor.
REFERENCE_RECEPTORS = {"TNFRSF10B": "DR5 (Step-1 target)", "TNFRSF10A": "DR4 (paralog)"}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", stream=sys.stdout)
log = logging.getLogger("step2b")


# ---------- pure helpers (unit-tested locally) ----------
def is_epithelial(cell_type: str) -> bool:
    c = cell_type.lower()
    return any(k in c for k in EPITHELIAL_KEYWORDS)


def select_cancer_celltypes(celltypes) -> tuple[list[str], str]:
    """Return (selected cell types, mode) for the cancer population in a tumour sample."""
    present = set(celltypes)
    explicit = [c for c in EXPLICIT_CANCER_LABELS if c in present]
    if explicit:
        return explicit, "explicit-malignant"
    epi = sorted({c for c in present if is_epithelial(c)})
    if epi:
        return epi, "epithelial-fallback"
    return [], "none"


def aggregate_across_tissues(per_tissue: dict) -> "object":
    """per_tissue: {tissue: DataFrame[receptor, log2fc, padj, pct_cancer, pct_normal]}.

    Returns a combined DataFrame ranked by pan-cancer enrichment:
    n_tissues_enriched (desc), then mean_log2fc (desc).
    """
    import pandas as pd
    frames = []
    for tissue, df in per_tissue.items():
        d = df.copy()
        d["tissue"] = tissue
        d["enriched"] = (d["log2fc"] >= LOG2FC_MIN) & (d["padj"] < PADJ_MAX) & (d["pct_cancer"] >= MIN_PCT_CANCER)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    alld = pd.concat(frames, ignore_index=True)
    g = alld.groupby("receptor")
    summary = pd.DataFrame({
        "n_tissues_enriched": g["enriched"].sum().astype(int),
        "n_tissues_tested": g.size().astype(int),
        "mean_log2fc": g["log2fc"].mean(),
        "max_log2fc": g["log2fc"].max(),
        "mean_pct_cancer": g["pct_cancer"].mean(),
        "min_padj": g["padj"].min(),
    }).reset_index()
    summary = summary.sort_values(["n_tissues_enriched", "mean_log2fc"], ascending=[False, False])
    return summary.reset_index(drop=True)


# ---------- receptor gene set ----------
def receptor_genes() -> set[str]:
    """Unique receptor subunit symbols from LIANA's curated consensus resource."""
    try:
        import liana
        res = liana.resource.select_resource("consensus")
        recs = set()
        for r in res["receptor"].astype(str):
            recs.update(r.split("_"))   # split complexes e.g. TGFBR1_TGFBR2
        recs.discard("")
        log.info("receptor set from LIANA consensus: %d unique subunits", len(recs))
        return recs
    except Exception as e:
        log.warning("LIANA resource unavailable (%s) — falling back to a small curated set", e)
        return {
            "TNFRSF10A", "TNFRSF10B", "FAS", "TNFRSF1A", "EGFR", "ERBB2", "ERBB3", "MET",
            "ALK", "KDR", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2", "IGF1R", "NOTCH1", "NOTCH2",
            "CD19", "MS4A1", "TNFRSF8", "TNFRSF17", "MSLN", "FOLR1", "EPCAM", "PROM1",
        }


# ---------- per-tissue differential ----------
def differential_for_tissue(tissue: str, receptors: set[str]):
    import numpy as np, scanpy as sc, pandas as pd
    tum_p = DATA_DIR / f"{tissue}__tumour.h5ad"
    nor_p = DATA_DIR / f"{tissue}__normal.h5ad"
    if not tum_p.exists() or not nor_p.exists():
        log.warning("[%s] missing h5ad (tumour=%s normal=%s) — skip", tissue, tum_p.exists(), nor_p.exists())
        return None

    tum = sc.read_h5ad(tum_p)
    nor = sc.read_h5ad(nor_p)
    log.info("[%s] tumour=%s normal=%s", tissue, tum.shape, nor.shape)

    # ---- characterise: cell-type composition ----
    log.info("[%s] tumour cell_type (top 10):", tissue)
    for ct, n in tum.obs["cell_type"].value_counts().head(10).items():
        log.info("        %6d  %s", n, ct)

    cancer_cts, mode = select_cancer_celltypes(tum.obs["cell_type"].unique())
    log.info("[%s] cancer-cell selection mode=%s → %s", tissue, mode,
             cancer_cts if len(cancer_cts) <= 8 else f"{cancer_cts[:8]}… ({len(cancer_cts)})")
    cancer_mask = tum.obs["cell_type"].isin(cancer_cts).to_numpy()
    normal_epi_cts = sorted({c for c in nor.obs["cell_type"].unique() if is_epithelial(c)})
    normal_mask = nor.obs["cell_type"].isin(normal_epi_cts).to_numpy()
    log.info("[%s] cancer cells=%d  normal-epithelial cells=%d", tissue, int(cancer_mask.sum()), int(normal_mask.sum()))
    if cancer_mask.sum() < MIN_CELLS or normal_mask.sum() < MIN_CELLS:
        log.warning("[%s] too few cells in a group (cancer=%d normal=%d) — skip",
                    tissue, int(cancer_mask.sum()), int(normal_mask.sum()))
        return None

    # ---- combine, normalise, Wilcoxon cancer vs normal ----
    cc = tum[cancer_mask].copy();  cc.obs["grp"] = "cancer"
    nn = nor[normal_mask].copy();  nn.obs["grp"] = "normal"
    common = cc.var_names.intersection(nn.var_names)
    comb = sc.concat([cc[:, common], nn[:, common]], join="inner")
    comb.obs["grp"] = comb.obs["grp"].astype("category")
    sc.pp.normalize_total(comb, target_sum=1e4)
    sc.pp.log1p(comb)

    sc.tl.rank_genes_groups(comb, "grp", groups=["cancer"], reference="normal",
                            method="wilcoxon", pts=True)
    res = sc.get.rank_genes_groups_df(comb, group="cancer")   # names, logfoldchanges, pvals_adj, pct_nz_group/reference
    res = res.rename(columns={"names": "receptor", "logfoldchanges": "log2fc", "pvals_adj": "padj"})
    pct_g = "pct_nz_group" if "pct_nz_group" in res.columns else None
    pct_r = "pct_nz_reference" if "pct_nz_reference" in res.columns else None
    res["pct_cancer"] = res[pct_g] if pct_g else np.nan
    res["pct_normal"] = res[pct_r] if pct_r else np.nan

    present_recs = [g for g in receptors if g in set(res["receptor"])]
    out = res[res["receptor"].isin(present_recs)][["receptor", "log2fc", "padj", "pct_cancer", "pct_normal"]].copy()
    out = out.sort_values("log2fc", ascending=False).reset_index(drop=True)
    log.info("[%s] receptors tested=%d  enriched(log2fc>=%.1f,padj<%.2f,pct>=%.0f%%)=%d",
             tissue, len(out), LOG2FC_MIN, PADJ_MAX, MIN_PCT_CANCER * 100,
             int(((out.log2fc >= LOG2FC_MIN) & (out.padj < PADJ_MAX) & (out.pct_cancer >= MIN_PCT_CANCER)).sum()))
    # reference receptor sanity
    for sym, desc in REFERENCE_RECEPTORS.items():
        row = out[out["receptor"] == sym]
        if len(row):
            r = row.iloc[0]
            log.info("[%s] ref %s (%s): log2fc=%.2f padj=%.1e pct_cancer=%.2f",
                     tissue, sym, desc, r["log2fc"], r["padj"], r["pct_cancer"])
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(DATA_DIR / f"targets_{tissue}.csv", index=False)
    return out


def main() -> int:
    log.info("cancer-recon-apoptosis — Step 2b — cancer-enriched receptor shortlist")
    try:
        import scanpy  # noqa: F401
        import pandas as pd  # noqa: F401
    except ImportError as e:
        log.error("missing dependency: %s (pip install scanpy liana)", e); return 2
    if not DATA_DIR.exists():
        log.error("no data at %s — run scripts/03 fetch first", DATA_DIR); return 3

    receptors = receptor_genes()
    per_tissue = {}
    for tissue in TISSUES:
        log.info("=" * 64)
        try:
            df = differential_for_tissue(tissue, receptors)
        except Exception as e:
            log.error("[%s] differential failed: %s: %s", tissue, type(e).__name__, e)
            df = None
        if df is not None and len(df):
            per_tissue[tissue] = df

    if not per_tissue:
        log.error("no tissue produced results — check inputs"); return 1

    summary = aggregate_across_tissues(per_tissue)
    out_path = DATA_DIR / "targets_shortlist.csv"
    summary.to_csv(out_path, index=False)

    log.info("=" * 64)
    log.info("PAN-CANCER RECEPTOR SHORTLIST (top 25 by tissues-enriched, then mean log2FC)")
    log.info("%-12s %6s %6s %9s %9s %9s", "receptor", "nEnr", "nTest", "meanL2FC", "maxL2FC", "meanPct")
    for _, r in summary.head(25).iterrows():
        log.info("%-12s %6d %6d %9.2f %9.2f %9.2f", r["receptor"], r["n_tissues_enriched"],
                 r["n_tissues_tested"], r["mean_log2fc"], r["max_log2fc"], r["mean_pct_cancer"])
    log.info("shortlist saved → %s (%d receptors)", out_path.relative_to(PROJECT_ROOT), len(summary))
    log.info("✓ done — feeds Step 3 (specificity audit). Communication annotation = scripts/05 (LIANA).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
