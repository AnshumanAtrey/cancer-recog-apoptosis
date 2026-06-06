# RUNG 11 — Public-neoantigen Addressability × Discriminability (the SECOND recognition axis)

**The pivot.** RUNG 5→10b exhausted ONE recognition axis — *"which genes are ON/OFF"* (the EXPRESSION axis
on the scRNA atlas) — and proved it bounded: surface logic gates can't be worst-donor-safe, no surface
NOT-blocker even exists, and the genetic HLA-LOH NOT-gate only reaches ~14–28% of patients. This rung opens
the axis AlphaFold/ESM actually live on: **"is this protein MUTATED?"** A somatic mutation is, *by
construction*, absent from every healthy cell — the tumour-**EXCLUSIVE** signal the expression axis could
never produce. Concretely: **public neoantigens** from recurrent driver hotspots (KRAS G12D/V/C, TP53 R175H,
PIK3CA H1047R/E545K, BRAF V600E, EGFR L858R, IDH1 R132H, CTNNB1, NRAS Q61) on common HLA-I alleles.

## How to run (laptop, no GPU, ~3 min)
```
pip install mhcflurry pandas numpy matplotlib scipy
mhcflurry-downloads fetch models_class1_presentation
python scripts/33_neoantigen_addressability.py            # -> JSON + figure
python scripts/33_neoantigen_addressability.py selftest   # offline logic checks (15/15, no deps)
```
Protein windows come from UniProt canonical sequences (pinned in `data/refs/uniprot_*.fasta`, validated
`wt == sequence[pos-1]`). 17 hotspots × 20 HLA-I alleles → 1,134 peptides → 22,680 (peptide,allele) scores.

## Two questions, answered empirically
1. **Addressability** — fraction of patients (per cancer) carrying a driver hotspot whose MUTANT peptide is
   presented (MHCflurry %rank ≤ 2) on an allele they carry = Σ P(mutation)·P(carries restricting allele),
   grouped so same-codon variants are mutually exclusive and genes independent.
2. **Discriminability (safety)** — does the WILD-TYPE counterpart ALSO get presented? Tiers:
   - **clean** — WT not presented → WT pMHC absent from healthy cells → cross-reactivity structurally
     impossible at the pMHC level (the safest handle).
   - **anchor** — mutation at an MHC anchor (P2/C-terminus) → MHC-level discrimination assists.
   - **tcr_dependent** — both presented → safety rests entirely on the TCR (the MAGE-A3 failure class →
     RUNG-12).

## Result (validated; oracle 3/3)
The pipeline recovers all three clinically-proven public-neoantigen handles as presented mutants
(KRAS-G12D/C\*08:02, KRAS-G12V/A\*11:01, TP53-R175H/A\*02:01) — built-in validation, not opinion.

- **The sequence axis decisively OUT-REACHES the expression axis.** BROAD addressability (all presented
  handles): **PDAC 73%, GLIOMA 57%, MELANOMA 43%, CRC 36%, NSCLC 29%** — all ≥ the HLA-LOH ceiling (14–28%),
  on a tumour-EXCLUSIVE signal. *The positive we were hunting.*
- **But most of that reach is `tcr_dependent`** (104 of 175 presented handles): the WT pMHC is also on
  healthy cells, so safety depends on the TCR telling two near-identical surfaces apart. The
  clinically-proven KRAS handles sit in this tier — proof that `tcr_dependent` ≠ impossible, but it ≠ free.
- **A structurally-SAFE core exists and is rankable.** SAFE (clean) addressability is reported as a
  threshold-honest range strict..lenient (strict = WT clearly off MHC, %rank > 4): **MELANOMA ~20% (robust,
  NRAS Q61R/K on A\*03:01/A\*11:01/A\*68:01, WT ranks 12–13)**, CRC ~8%, BRCA ~7%, PDAC ~5%.
- **8 GOLD handles** (strong mutant binder + WT clearly off MHC) — *every one is an anchor mutation*
  (mechanistically exact): PIK3CA E545K/A\*11:01·A\*03:01, EGFR L858R/A\*68:01, CTNNB1 S37F/C\*07,
  PIK3CA H1047R/C\*06:02. These need NO TCR finesse — the most defensible targets.

## Honest ceiling — five irreducible caveats (also in the JSON)
1. **PREDICTED, not measured** — MHCflurry presentation ≠ surface immunopeptidome; every hit is a hypothesis.
2. **Processing modelled coarsely** (proteasome/TAP).
3. **TCR-existence residual** — a presented mutant pMHC is necessary, not sufficient; a discriminating TCR
   must exist/be engineerable. That is RUNG-12 (AlphaFold-Multimer / ESM) + wet-lab. The `tcr_dependent` tier
   marks exactly where this residual dominates.
4. **Frequencies are population-dependent literature point estimates** joining independent datasets (not a
   same-patient cohort). East/South-Asian populations raise A\*11:01/A\*24:02 → KRAS coverage HIGHER there.
5. **Class-I / CD8 only.**

## Where this points
The `tcr_dependent` reach (incl. the proven KRAS handles) is the prize, gated by one question: **is the
mutant pMHC surface structurally distinguishable from wild-type?** → **RUNG-12**: AlphaFold-Multimer / ESM
model the top `tcr_dependent` + GOLD handles and quantify mut-vs-WT pMHC discriminability — certifying which
of the 73%/57%/43% BROAD reach is actually targetable. GPU earns its place there, exactly as it did for the
surface sweep.

## Provenance
`scripts/33_neoantigen_addressability.py` (selftest 15/15). MHCflurry 2.0
`Class1PresentationPredictor`. UniProt canonical sequences (`data/refs/uniprot_*.fasta`). Outputs:
`rung11_neoantigen_addressability.json`, `rung11_neoantigen.png`.
