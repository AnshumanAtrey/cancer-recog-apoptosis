# RUNG-30b — KRAS-G12D / HLA-A\*11:01 design target from a REAL CRYSTAL (PDB 7OW6) + an exposure caveat

The follow-up flagged in RUNG-32b: KRAS-G12D `VVVGADGVGK` has x-ray structures in IEDB, so design against a
real experimental pMHC instead of a fold. `scripts/66` stages it; the grab succeeded **and surfaced a finding.**

## What we got (a genuine upgrade)
- **PDB 7OW6 (2.64 Å)** — real HLA-**A\*11:01** + `VVVGADGVGK` (the *better*-validated KRAS allele: MS-eluted +
  T-cell-positive, R32b; the one UniProt lacked a clean sequence for). TCR stripped → groove (chain A) + peptide
  (chain C) → cropped to peptide + 10 Å groove. Higher quality than a fold, on the right allele.
- **Bonus:** extracted the **verified A\*11:01 α1α2 groove** from the crystal (7 polymorphisms vs A\*03:01:
  8F>Y, 89A>D, 104S>P, 151E>A, 155L>Q, 160D>E, 162T>R) — now we can fold the free A\*11:01 pMHC.
- Crystal choice: 7OW6 over 7OW4 (1.81 Å, free pMHC) because **7OW4 has p5–p6 incl. the G12D Asp DISORDERED**
  (unmodeled) — unusable; and over 7PB2 (3.41 Å, lower res). 7OW6 is the only good full-peptide structure.

## The finding (rule-5): the G12D Asp is NOT cleanly up-facing
SASA on the isolated pMHC (TCR stripped): **p6 G12D Asp = only ~5% exposed (8.8 Å²) = buried**; the most-exposed
peptide residue is p8-Val (23%, "mid"). This **contradicts the screen's position-heuristic** ("p6 central bulge
= up-facing"), which was a coarse proxy. (First version computed SASA on the full crystal *with the TCR* → ~0%
everywhere = artifact; corrected to isolated-pMHC.)

**But the conformation is ambiguous:** the only crystals are TCR-**bound** (7OW6/7PB2 = induced-fit conformation)
or free-but-Asp-disordered (7OW4) — so the *free*-pMHC Asp exposure is genuinely unknown. And TCRs **do** read
KRAS-G12D (gold-standard, R32b), so it IS recognizable — likely via induced fit + the buried charge — but a de
novo binder gripping a partially-buried charged residue is **harder than the "easy read-the-mutation" framing
assumed** (it echoes the IDH1/BRAF wall). Honest downgrade of KRAS's tractability — not a closure.

## Next (decided)
1. **Fold the FREE A\*11:01 pMHC** (Protenix, the extracted groove + `VVVGADGVGK`) → measure the Asp exposure in
   the free conformation vs the crystal. This resolves whether KRAS is a real read-the-mutation shot or another
   buried-mutation wall — **do this BEFORE spending PXDesign quota on KRAS.**
2. If the free Asp is accessible → PXDesign on the cropped target, hotspot **B6** (the Asp), max batch; KRAS WT
   is presented (R32) so MUT-vs-WT scoring is the make-or-break.
3. If still buried → KRAS de novo binder is bounded (like IDH1/BRAF); KRAS stays covered by the internal CRISPR
   key, and PIK3CA (presentation flip, no exposure requirement) becomes the cleaner external-binder shot.

*Files: `kras_g12d_A1101_mut_pmhc.pdb`, `kras_g12d_A1101_mut_cropped.pdb`, `meta.json`. Source 7OW6 gitignored.*
