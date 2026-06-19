# RUNG-30c — KRAS-G12D/A\*11:01 FREE pMHC fold → go/no-go RESOLVED: the G12D Asp IS up-facing (GO)

The gate from RUNG-30b: the crystal (7OW6) showed the G12D Asp only ~5% exposed, but that's a **TCR-bound**
conformation (induced fit). This folds the **free** pMHC (no TCR) to measure the Asp exposure a de novo binder
would actually meet. Protenix Add Prediction, A\*11:01 groove (extracted from 7OW6) + peptide, seed-matched,
template-on; downloaded 2026-06-19. Raw archived (`folds/{MUT,WT}/`, 5 CIFs + confidences + inputs each).

## Result — the crystal "buried" was an artifact; free Asp is EXPOSED
| | p6 G12D Asp SASA | exposure | iptm |
|---|---|---|---|
| **Free fold (this run)** | **58.4 Å²** | **30% — UP-facing** | 0.969 |
| TCR-bound crystal 7OW6 | 8.8 Å² | 5% — buried | (2.64 Å) |

In the free pMHC, p6-Asp is the most-exposed peptide residue (30%); p5-Ala 22%, p8-Val 19%. So the screen's
"p6 central-bulge up-facing" call was **right** — the crystal just had the TCR pushing the Asp down. (Bonus:
template-on did NOT bias the fold toward the crystal's buried Asp → the free conformation is genuine.)

## Verdict — KRAS-G12D/A\*11:01 is a real READ-THE-MUTATION target (GO)
- Gold-standard neoantigen (R32b: MS-eluted + x-ray + T-cell-positive on A\*11:01).
- **Gly→Asp = presence-vs-absence** (WT Gly has no sidechain) = the strongest discrimination chemistry, and the
  Asp is **accessible (30%)** → a de novo binder CAN grip it (unlike IDH1 His↔Arg or buried BRAF).
- WT IS presented (R32) → the binder MUST discriminate G12D from G12 — but with the Asp up-facing and Gly's
  empty space on WT, discrimination-by-construction is achievable.

## Design target staged (`../staging/`)
From the FREE fold (Asp up-facing), cropped to peptide + 10 Å groove, **hotspot = B6 (the G12D Asp)**. This is
the read-the-mutation analogue of the PIK3CA staging — but here we DO pin the mutation (it's up-facing). MUT-vs-WT
scoring will be the make-or-break (WT is real self, R32). Next: PXDesign Extended on the cropped target, hotspot
B6, max batch (after the PIK3CA design queue clears — webserver quota).

*Honest residual: free-fold conformation is a prediction (iptm 0.97); the real presented conformation + binder
affinity/specificity (SPR/cellular) remain the wet-lab line.*
