# RUNG-26e — BRAF-V600E / HLA-A\*01:01 pMHC folds (raw, landed from the Protenix webserver 2026-06-18)

The MUT + WT pMHC structures that the binder-design + specificity-scoring pipeline needs. Folded on the
Protenix webserver (Add Prediction, 2 protein chains each, **seed-matched 19931**, 5 samples), downloaded
2026-06-18 11:06. Archived here **bit-for-bit** (CIFs + per-sample summary confidences + inputs.json);
the bulky re-derivable MSAs (~3 MB/side) were intentionally not committed.

## What's here
- `MUT/` — groove (HLA-A\*01:01, 182 aa) + peptide **`ATEKSRWSGSH`** (V600E → **Glu at p3**). 5 CIF samples + confidences + inputs.json.
- `WT/`  — same groove + peptide **`ATVKSRWSGSH`** (germline **Val at p3**). 5 CIF samples + confidences + inputs.json.

## Register verified (the only difference is p3: E↔V)
Confirmed from both `inputs.json` and the folded CIF chains: chain A groove identical between MUT/WT;
chain B peptides differ **only at position 3** (MUT Glu / WT Val) = the V600E substitution. Seed identical
(19931) → MUT and WT are directly comparable for two-state specificity scoring.

## Fold confidence (Protenix iptm/ptm, top sample = sample_0)
| | sample_0 iptm | ptm | ranking |
|---|---|---|---|
| MUT | 0.9404 | 0.9676 | 0.9458 |
| WT  | 0.9421 | 0.9666 | 0.9470 |

Both pMHCs fold with high confidence (iptm ~0.94) — the peptide sits confidently in the groove in both
states. `sample_0` is top-ranked for both; use it as the canonical MUT/WT structure.

## What these unblock (next, built ON this committed data — not before it)
1. Pre-crop MUT `sample_0` (peptide + 10 Å groove) → PXDesign Extended upload target.
2. Stage MUT/WT `sample_0` + a `meta.json` (mut_pdb/wt_pdb/hotspot=B3) for MUT-vs-WT scoring on Drive.
3. PXDesign Extended on BRAF-MUT, hotspot **B3** → BRAF binder designs.
4. Score MUT vs WT (AF2 `binder_specificity` repointed + Protenix); if positive design doesn't discriminate,
   run the RUNG-26f negative-design notebook repointed to BRAF on the BRAF backbones.

**Why BRAF (vs IDH1):** V→E is small-hydrophobic→charged-carboxylate — a binder optimized for the Glu
leaves a buried unsatisfied charge on WT-Val → a real discrimination handle that positive design *can*
exploit, unlike IDH1's His↔Arg (RUNG-26c/d NULL). Honest: high fold confidence ≠ binder specificity;
that's still the open test.
