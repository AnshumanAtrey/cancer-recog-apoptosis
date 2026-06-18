# RUNG-29 — PIK3CA-E545K / HLA-A\*03:01 pMHC folds (raw, landed from Protenix webserver 2026-06-19)

The design + scoring targets for the screen-chosen presentation-flip binder target (RUNG-28/29). Protenix
webserver Add Prediction, 2 protein chains each, **seed-matched**, 5 samples; downloaded 2026-06-19 03:07.
Archived bit-for-bit (CIFs + per-sample confidences + inputs.json); bulky re-derivable MSAs excluded.

## Register verified — only p11 differs (E545K)
- `MUT/` — groove HLA-A\*03:01 (182 aa) + peptide **`STRDPLSEITK`** (E545K → **Lys at p11**, the basic-C-term anchor).
- `WT/`  — same groove + **`STRDPLSEITE`** (germline **Glu at p11**).
- Chain A (groove) identical MUT/WT; chain B differs **only at p11** (K↔E) = the E545K substitution.

## Fold confidence (Protenix iptm, sample_0)
| | iptm | ptm |
|---|---|---|
| MUT | 0.954 | 0.971 |
| WT  | 0.968 | 0.971 |

Both pMHCs fold confidently — **and WT folds slightly higher**, which is the expected reminder that *fold
confidence ≠ presentation*. Protenix places either peptide in the groove regardless of binding affinity; our
selectivity rests on the **MHCflurry presentation flip** (MUT 51 nM / WT 20 µM, ~397×, RUNG-28), not on the
fold. sample_0 = canonical structure for design/scoring.

## Next (built on this commit)
1. Pre-crop MUT sample_0 (peptide + 10 Å groove) → PXDesign upload target. **NO hotspot on the mutation**
   (p11 is a buried C-terminal anchor) — design a STRONG binder to the up-facing surface; selectivity is
   supplied by presentation, not by the binder reading the mutation.
2. Stage MUT/WT sample_0 + meta.json (for AF2/Protenix MUT-vs-WT scoring as a sanity, though selectivity does
   not depend on the binder discriminating).
3. PXDesign Extended + ODesign cross-check → strong MUT-pMHC binder = the artifact.
