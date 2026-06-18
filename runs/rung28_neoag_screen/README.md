# RUNG-28 — neoantigen target-selection screen → PIK3CA-E545K is the binder target (the step we skipped)

We hand-picked IDH1-R132H then BRAF-V600E by *driver prestige* and both failed for **target-geometry**
reasons, not method failure:
- **IDH1-R132H** — His↔Arg subtle + binder wraps the conserved core (RUNG-26c/d NULL across AF2-IG, ColabDesign-AF2, Protenix; negative design RUNG-26f marginal).
- **BRAF-V600E** — the mutation is a **buried p3 anchor** (2% solvent-exposed) with no reproducible surface change, AND the peptide barely presents on A\*01:01 (MHCflurry MUT 3971 nM / pres 0.025). (RUNG-26e audit.)

**The criterion we never checked:** a de-novo binder can be mutation-specific only if (A) the MUT peptide is
**strongly presented**, AND either (B1) the mutated residue is **solvent-exposed/up-facing** (binder reads it)
or (B2) the mutation **flips presentation** (MUT presented, WT not → selectivity is built into the MHC, the
binder only needs to bind MUT). This screen does the cheap half — (A) + (B2) — with **MHCflurry** over the
16-hotspot driver panel × 13 common HLA-I alleles (CDS-derived mutant peptides, `scripts/59`).

## Result — 3 presentation-FLIP targets (MUT presented, WT ≥5× weaker); PIK3CA-E545K is the standout
| driver | allele | MUT peptide | MUT affinity / pres | WT affinity / pres | flip |
|---|---|---|---|---|---|
| **PIK3CA E545K** | **HLA-A\*03:01** | `STRDPLSEITK` | **51 nM / 0.83** | 20331 nM / 0.006 | **397×** |
| **PIK3CA E545K** | **HLA-A\*11:01** | `STRDPLSEITK` | **42 nM / 0.85** | 13227 nM / 0.009 | **318×** |
| TP53 R248W | HLA-B\*57:01 | `SSCMGGMNW` | 39 nM | 7194 nM | 185× |

37 (driver,allele) pairs present the MUT strongly (≤500 nM & pres ≥0.5); most (KRAS-G12 on A\*03/A\*11, NRAS-Q61
on A\*01:01, IDH1-R132H on B\*35:01) present WT *equally* → those need the **read-the-mutation** route (exposure
check). Only 3 FLIP.

## Why PIK3CA-E545K / A\*03:01 (+ A\*11:01) is THE target — audited mechanism
E545K sits at the **C-terminal anchor** of `STRDPLSEITK` (MUT **K** / WT **E**, position 11). A\*03:01 and
A\*11:01 are the **basic-C-terminus** alleles (K/R anchor) → mutant **Lys** = perfect anchor (42-51 nM);
germline **Glu** (acidic) = anchor-dead (13-20 µM, unpresented). **The mutation creates the anchor that
displays it.** So WT pMHC is essentially absent from normal cells → a binder to the MUT pMHC is *inherently*
mutation-specific. This **sidesteps the wall that killed IDH1/BRAF** (binder needn't read a subtle/buried
residue); and PXDesign already proved it builds strong pMHC binders (RUNG-26c/d) — it only lacked
discrimination, which here is supplied by presentation. Coverage: A\*03:01 (~European-common) + A\*11:01
(~Asian-common); PIK3CA E545K is a top breast/CRC/endometrial driver.

## Honest residuals
- **Processing:** requires the C-terminus to be generated exactly at residue 545 (proteasomal cut after Lys).
  MHCflurry's processing score supports it, but mass-spec immunopeptidomics is the real proof.
- **Flip margin:** WT is 13-20 µM = *weak*, not strictly zero → residual WT cross-presentation risk, mitigated
  by the ~300× margin.
- **In-silico:** MHCflurry = a presentation predictor; binder affinity/specificity (SPR/cellular) and the fold
  remain the wet-lab residual.
- Read-the-mutation targets (KRAS-G12 etc.) are NOT closed — they need the exposure check (fold + SASA); logged
  for a later rung, not dropped.

## Next (decided, in motion)
Fold PIK3CA-E545K MUT pMHC (A\*03:01 groove + `STRDPLSEITK`) on the Protenix webserver → PXDesign Extended →
strong MUT binder (no negative design needed; presentation supplies selectivity) → score; A\*11:01 as the
second allele. `scripts/59` + `screen.json` hold the full ranked table.
