# RUNG-30 — de novo binder to KRAS-G12D/HLA-A\*03:01 (the read-the-mutation target; complements PIK3CA)

RUNG-28's position triage (`../rung28_neoag_screen/position_triage.json`) found the KRAS-G12 family as the
best **read-the-mutation** candidate — the opposite mechanism to PIK3CA-E545K's presentation flip. Two
complementary shots on goal.

## Why KRAS-G12D is the right read-the-mutation target (and why it should work where IDH1 didn't)
- **Strongly presented, mutation up-facing.** `VVVGADGVGK` on A\*03:01 (71 nM, pres 0.87; A\*11:01 47 nM); the
  mutation sits at **p6 — the central bulge**, the classic TCR/binder-exposed position (not an anchor).
- **Presence-vs-absence chemistry.** G12 is **Glycine** — WT has *no sidechain* at p6. G12D adds an **Asp**.
  A binder pocket built around the Asp finds **empty space** (and loses a salt bridge) on WT-Gly → real
  discrimination. This is the *strongest* possible handle — fundamentally unlike IDH1's His↔Arg (swap one
  bulky/H-bonding residue for another, which a single pocket tolerates → RUNG-26 NULL) or BRAF's buried anchor.
- **Both MUT and WT present (~equal)** → no presentation flip → selectivity *must* come from the binder, so
  the **hotspot goes ON the mutation (p6)** — and here that can actually work.
- KRAS-G12D is the **#1 oncogenic driver** (PDAC, CRC, LUAD); A\*03:01 common (A\*11:01, more South-Asian,
  to be added from IMGT).

## Pipeline (proven RUNG-26/29 machinery)
1. **Fold** MUT (`VVVGADGVGK`) + WT (`VVVGAGGVGK`) pMHC on Protenix, 2 chains (groove_A0301 + peptide),
   seed-matched, 5 samples. Download both.
2. **Pre-crop + exposure audit** (`scripts/58`): confirm p6-Asp is up-facing/solvent-exposed (the same SASA +
   within-vs-between check that *closed* BRAF — here we expect it to PASS, central-bulge).
3. **PXDesign Extended**, target = cropped MUT pMHC, **hotspot = p6** (the Asp), full batch.
4. **Score** MUT vs WT (AF2 `binder_specificity` repointed + Protenix). WIN = WT clearly worse (pae↑≥3 /
   iptm↓≥0.15) — the binder grips the Asp and finds nothing on Gly. If positive design + hotspot is marginal,
   apply the **RUNG-26f negative-design** method (now built + validated).

## Honest framing
- This is the read-the-mutation route's *best* shot, but it's still unproven — IDH1/BRAF taught us not to
  assume. The exposure audit (step 2) is the go/no-go before PXDesign quota, exactly as for BRAF.
- In-silico fold/affinity ≠ measured; immunopeptidomic presentation of `VVVGADGVGK` (well-documented KRAS
  neoantigen) + binder SPR/cellular = wet-lab residual.
- Runs in parallel with PIK3CA-E545K (presentation flip) — two mechanisms, two targets, decided by data.
