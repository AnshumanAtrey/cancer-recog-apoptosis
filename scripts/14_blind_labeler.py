#!/usr/bin/env python3
"""
RUNG 2 / Stage-pre — BLIND geometry labeler (adversary C1 fix).

Assigns the GEOMETRY inputs (arm-reach r_arm, epitope proximity, geometry-match) that the clustering
sim (scripts/13) and the ranker (scripts/12) consume — from FORMAT + epitope_crd (STRUCTURAL facts)
ONLY. It is *endpoint-blind by construction*: it loads dr5_agonist_ladder.csv and DROPS `potency_rank`
before computing anything (asserted below). No human in the loop; deterministic rules.

WHY (the honest core): the four bivalent IgG1s have equal valency and hinge, so blind structure has no
basis to separate them -> their geometry inputs come out (near-)identical BY DESIGN. If the within-valency
contrast can't be resolved from blind structure, RUNG 2 must report a valency lookup table, not a
geometry predictor. The one feature that *could* separate them (the literature crosslink-independence
flag for tigatuzumab) is recorded but EXCLUDED from the blind vector (construct-level leakage), and the
epitope feature is gated by the scrambled-epitope negative control in scripts/12.

USAGE:  python scripts/14_blind_labeler.py
OUT  :  data/dr5_agonists/blind_geometry_labels.csv  (+ refreshes label_provenance footer)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LADDER = PROJECT_ROOT / "data" / "dr5_agonists" / "dr5_agonist_ladder.csv"
OUT = PROJECT_ROOT / "data" / "dr5_agonists" / "blind_geometry_labels.csv"

# --- structural priors: arm-reach / construct SPAN (nm) from CONSTRUCT FORMAT only (endpoint-blind) ---
# A multivalent scaffold's physical span grows with its size/arm-count: a single domain spans little, a
# decavalent IgM spans tens of nm. Span sets how many receptors the construct can simultaneously engage.
# nanobody ~4nm (monovalent); IgG1 Fab-Fab ~13nm; native TRAIL trimer ~9nm; tetravalent sdAb-Fc ~14nm;
# hexavalent TRAIL-Fc ~18nm; IgM decavalent pentamer-of-dimers ~22nm. (Structural facts, not potency.)
R_ARM_NM = {
    "nanobody": 4.0,
    "IgG1": 13.0,
    "trimer": 9.0,
    "tetravalent-sdAb-Fc": 14.0,
    "hexavalent-TRAIL-Fc": 18.0,
    "IgM-decavalent": 22.0,
}
IDEAL_PITCH_NM = 6.0   # DR5 dimer-of-trimers hexagonal pitch (the lattice that propagates)
PITCH_WIDTH = 4.0      # gaussian tolerance for geometry-match

# epitope proximity from epitope_crd (STRUCTURAL binding-site fact): membrane-distal CRD1 -> proximal CRD3.
# Proximal binding geometrically favours productive lattice contacts. Literature-derived -> gated by the
# scrambled-epitope negative control in scripts/12.
EPITOPE_PROX = {
    "CRD1": 0.0, "CRD2": 1.0, "CRD2-CRD3": 1.5, "CRD3": 2.0, "CRD3-proximal": 2.3,
    "native-trimeric": 1.5, "native-TRAIL": 1.5,
}


def geometry_match(r_arm_nm: float) -> float:
    """g in [0,1]: how compatible the arm-reach is with the ~6nm dimer-of-trimers pitch."""
    import math
    return float(math.exp(-((r_arm_nm - IDEAL_PITCH_NM) / PITCH_WIDTH) ** 2))


def main() -> int:
    raw = pd.read_csv(LADDER, comment="#")
    # --- ENFORCE BLINDNESS: drop the endpoint before computing any feature ---
    assert "potency_rank" in raw.columns, "expected potency_rank column in Table A"
    blind = raw.drop(columns=["potency_rank", "crosslink_dependence", "assay_cell_line", "source"])
    #            ^ potency_rank = the y; crosslink_dependence = literature leakage flag (excluded)

    rows = []
    for _, r in blind.iterrows():
        fmt = str(r["format"]).strip()
        r_arm = R_ARM_NM.get(fmt)
        if r_arm is None:
            print(f"[blind-labeler] WARNING unknown format '{fmt}' for {r['name']} — defaulting r_arm=10nm")
            r_arm = 10.0
        epi = EPITOPE_PROX.get(str(r["epitope_crd"]).strip(), 1.0)
        rows.append({
            "name": r["name"],
            "format": fmt,
            "valency": int(r["valency"]),
            "r_arm_nm": r_arm,
            "epitope_proximity": epi,
            "geometry_match": round(geometry_match(r_arm), 4),
        })
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    # honest check: are the bivalent IgG1s tied on blind geometry (excluding the leakage flag)?
    biv = out[out["valency"] == 2]
    arm_tied = biv["r_arm_nm"].nunique() == 1
    gm_tied = biv["geometry_match"].nunique() == 1
    epi_unique = biv["epitope_proximity"].nunique()
    print("[blind-labeler] wrote", OUT.relative_to(PROJECT_ROOT))
    print(out.to_string(index=False))
    print(f"[blind-labeler] bivalent IgG1s: r_arm tied={arm_tied}, geometry_match tied={gm_tied}, "
          f"epitope distinct values={epi_unique}/4")
    print("[blind-labeler] -> blind structure separates the bivalents ONLY via epitope_proximity "
          "(literature-derived); the scrambled-epitope negative control in scripts/12 is the decisive gate.")
    print("[blind-labeler] potency_rank was DROPPED before any feature was computed (blind by construction).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
