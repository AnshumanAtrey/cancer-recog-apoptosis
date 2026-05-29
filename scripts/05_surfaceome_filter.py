#!/usr/bin/env python3
"""
Step 2b.2 — surfaceome filter: keep only bona-fide cell-surface receptors.

The receptor differential (scripts/04) ranks by cancer-vs-normal enrichment, but
LIANA's receptor list includes proteins that are NOT cell-surface with an ectodomain
(e.g. FLOT1 flotillin, RACK1 scaffold, MAGED1) — you cannot design an extracellular
ligand against those. This script intersects the shortlist with the human surfaceome
so the output is an ACTIONABLE ligand-design target list.

SURFACEOME SOURCE (primary → fallback):
  1. OmniPath intercell locational annotation 'plasma_membrane_transmembrane'
     (transmembrane PM proteins → have an extracellular domain). OmniPath incorporates
     the Bausch-Fluck SURFY surfaceome (PNAS 2018) + others. Already a LIANA dependency.
  2. Curated static list of well-characterised druggable surface receptors (offline).

OUTPUT: data/cellxgene/targets_surface_shortlist.csv  — surface receptors only,
ranked (tissues-enriched desc, mean log2FC desc) + the actionable set
(surface AND enriched in >=2 tissues).

USAGE:  python scripts/05_surfaceome_filter.py
REQS :  pandas, omnipath (CPU; no GPU)
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "cellxgene"
SHORTLIST = DATA_DIR / "targets_shortlist.csv"
SURFACE_GENES_CACHE = DATA_DIR / "surfaceome_genes.txt"
OUT = DATA_DIR / "targets_surface_shortlist.csv"

# Actionable-target thresholds (mirror scripts/04).
MIN_TISSUES_ENRICHED = 2

# Offline fallback — well-characterised druggable cell-surface receptors (symbols).
CURATED_SURFACE = {
    "EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "ALK", "ROS1", "RET", "KDR", "FLT1", "FLT4",
    "PDGFRA", "PDGFRB", "KIT", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "IGF1R", "INSR", "NTRK1",
    "NTRK2", "NTRK3", "AXL", "MERTK", "EPHA2", "EPHA3", "EPHB2", "EPHB4", "NOTCH1", "NOTCH2",
    "NOTCH3", "NOTCH4", "TNFRSF10A", "TNFRSF10B", "FAS", "TNFRSF1A", "TNFRSF1B", "TNFRSF8",
    "TNFRSF17", "CD19", "MS4A1", "CD22", "CD33", "CD38", "CD74", "MSLN", "FOLR1", "FOLH1",
    "EPCAM", "PROM1", "CEACAM5", "MUC1", "ERBB2", "ADGRG1", "ADGRE5", "SDC1", "SDC2", "SDC3",
    "SDC4", "ITGB4", "ITGA3", "ITGAV", "ITGB1", "CDH1", "CDH3", "MCAM", "ALCAM", "NCAM1",
    "L1CAM", "PTPRF", "PTK7", "ROR1", "ROR2", "LRP5", "LRP6", "TFRC", "SLC39A6", "TSPAN8",
    "CLDN6", "CLDN18", "GPC3", "DLL3", "NECTIN4", "TACSTD2", "MET", "LGR5",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", stream=sys.stdout)
log = logging.getLogger("step2b2")


def surface_genes_from_omnipath() -> set[str]:
    """Genes annotated by OmniPath as plasma-membrane transmembrane (have an ectodomain)."""
    import omnipath
    ic = omnipath.requests.Intercell.get()           # intercell annotation table
    cols = list(ic.columns)
    log.info("omnipath intercell rows=%d cols=%s", len(ic), cols)
    cat_col = "category" if "category" in cols else ("parent" if "parent" in cols else None)
    gene_col = "genesymbol" if "genesymbol" in cols else ("genesymbols" if "genesymbols" in cols else None)
    if cat_col is None or gene_col is None:
        raise RuntimeError(f"unexpected intercell schema: {cols}")
    cats = set(ic[cat_col].dropna().unique())
    # categories that mean "transmembrane / on the cell surface"
    want = {c for c in cats if any(k in str(c).lower()
            for k in ("plasma_membrane_transmembrane", "cell_surface", "surfaceome"))}
    log.info("surface-indicating categories matched: %s", sorted(want) or "NONE")
    if not want:
        raise RuntimeError(f"no surface categories found in {sorted(cats)[:20]}…")
    genes = set(ic[ic[cat_col].isin(want)][gene_col].dropna().astype(str))
    genes = {g for g in genes if g and g != "nan"}
    log.info("omnipath surface gene set: %d genes", len(genes))
    return genes


def get_surface_genes() -> tuple[set[str], str]:
    # cached file (reproducible / offline re-runs)
    if SURFACE_GENES_CACHE.exists():
        genes = {ln.strip() for ln in SURFACE_GENES_CACHE.read_text().splitlines() if ln.strip()}
        if genes:
            log.info("surface gene set from cache %s: %d genes", SURFACE_GENES_CACHE.name, len(genes))
            return genes, "cache"
    try:
        genes = surface_genes_from_omnipath()
        if genes:
            SURFACE_GENES_CACHE.parent.mkdir(parents=True, exist_ok=True)
            SURFACE_GENES_CACHE.write_text("\n".join(sorted(genes)) + "\n")
            log.info("cached surface gene set → %s", SURFACE_GENES_CACHE.relative_to(PROJECT_ROOT))
            return genes, "omnipath"
    except Exception as e:
        log.warning("omnipath surfaceome unavailable (%s: %s) — using curated fallback",
                    type(e).__name__, e)
    return set(CURATED_SURFACE), "curated-fallback"


def main() -> int:
    log.info("cancer-recon-apoptosis — Step 2b.2 — surfaceome filter")
    try:
        import pandas as pd
    except ImportError as e:
        log.error("missing dependency: %s", e); return 2
    if not SHORTLIST.exists():
        log.error("no shortlist at %s — run scripts/04 first", SHORTLIST); return 3

    df = pd.read_csv(SHORTLIST)
    log.info("loaded shortlist: %d receptors", len(df))
    surface, src = get_surface_genes()
    log.info("surface gene set source=%s size=%d", src, len(surface))

    df["is_surface"] = df["receptor"].isin(surface)
    n_surface = int(df["is_surface"].sum())
    log.info("of %d receptors, %d are cell-surface (%.0f%%)", len(df), n_surface, 100 * n_surface / max(len(df), 1))

    surf = df[df["is_surface"]].copy()
    surf = surf.sort_values(["n_tissues_enriched", "mean_log2fc"], ascending=[False, False]).reset_index(drop=True)
    surf.to_csv(OUT, index=False)

    # how many of the original top-25 were surface (shows the filter's effect)
    top25 = df.head(25)
    log.info("of the original top-25, %d/%d are cell-surface", int(top25["is_surface"].sum()), len(top25))

    # actionable set
    action = surf[(surf["n_tissues_enriched"] >= MIN_TISSUES_ENRICHED) & (surf["mean_log2fc"] >= 0.5)]
    log.info("=" * 64)
    log.info("ACTIONABLE SURFACE TARGETS (surface ∧ enriched in >=%d tissues ∧ meanL2FC>=0.5): %d",
             MIN_TISSUES_ENRICHED, len(action))
    log.info("%-12s %5s %9s %9s %9s", "receptor", "nEnr", "meanL2FC", "maxL2FC", "meanPct")
    for _, r in action.head(30).iterrows():
        log.info("%-12s %5d %9.2f %9.2f %9.2f", r["receptor"], r["n_tissues_enriched"],
                 r["mean_log2fc"], r["max_log2fc"], r["mean_pct_cancer"])

    # sanity: DR5/DR4 should be flagged surface (they are TM receptors)
    for sym in ("TNFRSF10B", "TNFRSF10A"):
        row = df[df["receptor"] == sym]
        if len(row):
            log.info("sanity %s is_surface=%s (expect True)", sym, bool(row.iloc[0]["is_surface"]))

    log.info("surface shortlist saved → %s (%d receptors)", OUT.relative_to(PROJECT_ROOT), len(surf))
    log.info("✓ done — pick the target thesis from the actionable surface list.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
