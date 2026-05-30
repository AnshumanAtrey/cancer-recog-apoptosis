# RUNG 2 — label provenance & blinding statement

**Why this file exists (adversary C1 fix).** RUNG 2 asks whether clustering *geometry* adds
predictive signal *beyond valency*. The within-valency contrast (the four bivalent IgG1s) is the
only place valency is mathematically blind, so any "geometry signal" there must not be a label
written while looking at the answer key. This file records exactly where every geometry feature
value comes from, and certifies the labeler was blind to the potency rank.

## Blinding certification

`scripts/14_blind_labeler.py` is a **deterministic rule script**. It reads
`dr5_agonist_ladder.csv` and **drops the `potency_rank` column before computing any feature**
(asserted in code). It assigns geometry inputs from **format + epitope_crd (structural facts) only**:

| input | source | endpoint-blind? |
|-------|--------|-----------------|
| `valency` | construct format (# DR5-binding arms) | yes — structural |
| `r_arm_nm` | format → known arm-reach (nanobody ~4 nm, IgG1 Fab–Fab ~13 nm, TRAIL trimer ~6 nm pitch, IgM pentamer ~15 nm, sdAb-Fc ~8 nm) | yes — structural |
| `epitope_proximity` | `epitope_crd` column (CRD1 distal → CRD3 membrane-proximal) | structural, BUT see leakage note |
| `geometry_match` | computed: r_arm vs ~6 nm dimer-of-trimers pitch | yes — derived |

## Leakage note (honest, load-bearing)

- The `epitope_crd` values and the `crosslink_dependence` flags in Table A are drawn from the
  **same primary literature** as the potency ranks (construct-level overlap). Pre-registration and
  blinding cannot fully separate them. **Therefore:** the literature `crosslink_dependence` flag is
  recorded but **EXCLUDED from the default blind feature vector** (`use_literature_flag=False` in the
  manifest). With it excluded, the four bivalent IgG1s have **near-identical blind geometry by design**,
  because blind structure has no basis to separate four IgG1s of equal valency and hinge. That is the
  honest core: *if the within-valency contrast cannot be resolved from blind structure, RUNG 2 must
  report a valency lookup table, not a geometry predictor.*
- The `epitope_proximity` feature IS allowed (it is a structural binding-site fact), but because it
  too is literature-derived, the **scrambled-epitope negative-control arm** in `scripts/12` is the
  decisive gate: if scrambling the epitope labels reproduces any within-valency ordering, the geometry
  headline is VOID.

## Boltz-2 status on this receptor (oracle ceiling)

The geometry inputs are **literature/format-derived, not Boltz-2-derived**, in this committed build,
because `runs/step1_boltz/interface_metrics_colab_transcript.json` shows Boltz-2 cannot discriminate a
binder from a non-binder on this exact DR5 ectodomain (a scrambled decoy beat real TRAIL on ipTM; lysozyme
tied on pDockQ). The Colab notebook carries the **multi-seed Boltz-2 resolution probe** (seed-to-seed spread
of arm spacing/angle on the three bivalents); the manifest records that if within-molecule spread ≥
between-molecule difference, the same-valency inputs are declared oracle-noise-dominated and the
within-valency test is label-driven, not physics-driven.

_Labeler: deterministic script (no human in the loop), blind-to-endpoint by construction. Committed before the rank was used in any fit._
