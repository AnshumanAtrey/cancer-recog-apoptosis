# RUNG-29 — de novo binder to PIK3CA-E545K/HLA-A\*03:01 (the RUNG-28 screen winner; the *right* target)

**Why this, not IDH1/BRAF:** RUNG-28 (MHCflurry, 16 drivers × 13 HLA-I) found PIK3CA-E545K as a
**presentation-flip** target — MUT `STRDPLSEITK` binds A\*03:01 at **51 nM** (A\*11:01 42 nM) while WT
`STRDPLSEITE` is **20 µM** (A\*11:01 13 µM) = **300-400× flip**. E545K is the **C-terminal anchor** (p11):
mutant **Lys** is the perfect A\*03/A\*11 basic-C-term anchor; germline **Glu** is anchor-dead. **The mutation
creates the anchor that displays the peptide** → WT pMHC is essentially absent from normal cells.

**The key consequence:** a binder to the MUT pMHC is **inherently mutation-specific by presentation** — it does
NOT need to read the mutation. This sidesteps the wall that sank IDH1 (His↔Arg subtle, wrapped) and BRAF
(V600E buried + weakly presented). PXDesign already proved (RUNG-26c/d) it builds strong pMHC binders on a ≤5%
target — it only lacked discrimination, which here is supplied for free by the MHC.

## Pipeline (reuses the proven RUNG-26 machinery; one critical change)
1. **Fold MUT pMHC** — Protenix webserver Add Prediction, 2 protein chains: groove = `groove_A0301`
   (`pik3ca_e545k_inputs.json`) + peptide `STRDPLSEITK`. (Also fold WT `STRDPLSEITE` for the record.) Seed-match.
2. **Pre-crop** MUT to peptide + 10 Å groove (`scripts/58` machinery) → PXDesign upload target.
3. **PXDesign Extended** on MUT-pMHC — **NO hotspot on the mutation** (p11 is a buried C-terminal anchor, not
   surface). Free footprint / hotspot on an up-facing residue → a STRONG binder. Full batch.
4. **Confirm** the top binder folds confidently vs the MUT pMHC (AF2 + Protenix). Optionally fold vs WT pMHC —
   but selectivity does NOT depend on the binder discriminating; it depends on WT not being presented. So the
   bar here is simply: a strong, confident MUT-pMHC binder. **No negative design needed.**
5. Repeat on **A\*11:01** (second common allele; even stronger flip) for coverage.

## What changed vs the IDH1/BRAF runs
- Target chosen by **screen, not prestige** (presentation flip + strong MUT presentation).
- Discrimination mechanism = **differential presentation**, not binder-reads-mutation → the design problem
  collapses to "make a strong pMHC binder," which is the thing PXDesign is good at.

## Honest residuals (carried from RUNG-28)
- C-terminus must be proteasomally generated at residue 545 (MHCflurry processing supportive; immunopeptidomics
  is the real proof).
- WT 13-20 µM is *weak*, not zero → residual WT cross-presentation, mitigated by the ~300× margin.
- In-silico fold/affinity ≠ measured affinity/specificity (SPR/cellular) — the standard wet-lab residual.
