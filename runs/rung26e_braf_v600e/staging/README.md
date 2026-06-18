# RUNG-26e staging + AUDIT — BRAF-V600E/A\*01:01 CLOSED as a binder target (mutation buried + weakly presented)

Staged the committed BRAF folds into binder-pipeline artifacts (`scripts/58`) AND audited whether V→E is
actually a viable de-novo-binder discrimination target *before* spending PXDesign quota (rules 1 + 5). It is not.

## Artifacts (staging/)
- `braf_mut_pmhc.pdb` / `braf_wt_pmhc.pdb` — full pMHC (chains A=groove, B=peptide), sample_0 = scoring targets.
- `braf_mut_cropped.pdb` — peptide + 10 Å groove (PXDesign upload target, had it been viable).
- `meta.json` — mut/wt paths + hotspot=3. `readability.json` — the p3 exposure check.

## AUDIT — why BRAF is closed (two independent failures)
1. **Mutation is BURIED, not readable.** p3 (V600E) solvent exposure = **2%** in both MUT-Glu (5.5 Å²) and
   WT-Val (2.4 Å²); the up-facing residues (p4-Lys 67-92 Å², p6-Arg 120-136, p10-Ser 72-74) are **conserved**
   MUT=WT. A binder above the peptide reads conserved surface → can't discriminate (same failure as IDH1).
2. **No reproducible surface signal from the buried mutation.** sample_0 showed a p7/p8 bulge difference, but
   across **5 samples/side** the within-WT spread (up to 3.9 Å at p8) swamps the between-state difference
   (ratio 0.4-1.4 = noise). The bulge is intrinsically disordered; V600E does not deterministically reshape it
   (p3 itself: between 0.49 Å vs within 0.24/0.31, ratio 1.8 = marginal).
3. **Weakly presented anyway.** MHCflurry A\*01:01: MUT 3971 nM / pres 0.025, WT 20328 nM / pres 0.006. The
   V600E does create differential presentation (~5×, P3-acidic-anchor direction confirmed) but **MUT itself is
   a poor binder** — consistent with the known immunology that BRAF-V600E presents poorly on class-I HLA.

→ BRAF-V600E/A\*01:01 fails both the read-the-mutation route (buried) and the differential-presentation route
(MUT too weak). My earlier "V→E gives a binder a handle" rationale assumed p3 was exposed; the fold says it
isn't. Quota saved. **Replaced by the RUNG-28 screen winner PIK3CA-E545K/A\*03:01** (300-400× presentation flip,
MUT 42-51 nM). The BRAF folds remain archived in `../folds/` as the run record.
