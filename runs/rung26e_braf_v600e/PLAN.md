# RUNG-26e — de novo mutation-specific binder for BRAF-V600E/HLA-A*01:01 (the tractable discrimination target)

**Why BRAF-V600E, not IDH1-R132H:** RUNG-26c/d proved de novo binders to IDH1-R132H bind the pMHC but can't discriminate His↔Arg (both bulky → a pocket fitting one tolerates the other; the binder wraps all 10 conserved residues and the mutation is energetically marginal). **BRAF V600E is V→E** — small hydrophobic → charged carboxylate. A binder optimized for the Glu (e.g. a salt-bridging Lys/Arg) leaves a **buried unsatisfied charge** when the residue reverts to WT-Val → real destabilization → a genuine discrimination handle positive design *can* exploit. RUNG-27a: BRAF-V600E has a real catalogued epitope but **no natural class-I TCR** → de novo is the route. Melanoma ~50% driver.

**Target (`braf_v600e_inputs.json`):** HLA-A\*01:01 (same groove we already have) + peptide `ATEKSRWSGSH` (mut) / `ATVKSRWSGSH` (wt). Mutation at **p3** (V→E). Hotspot = **B3**.

## Pipeline (same one proven in RUNG-26c/d)
1. **Fold** BRAF-V600E pMHC AND WT pMHC — Protenix webserver "Add Prediction", 2 chains each:
   - chain1 (MHC) = the groove in `braf_v600e_inputs.json`
   - chain2 (peptide) = `ATEKSRWSGSH` (MUT job) / `ATVKSRWSGSH` (WT job)
   → download both CIFs (these are also the MUT/WT scoring targets).
2. **Pre-crop** the MUT pMHC to peptide + 10 Å groove (script mirrors the IDH1 pre-crop), upload as the design target.
3. **PXDesign Extended**, target = cropped BRAF-MUT pMHC, **hotspot = B3** (the Glu), binder length ~80–100, full batch.
4. **Score** the Protenix-passers vs MUT and WT pMHC (AF2 `binder_specificity` notebook repointed here, + Protenix webserver MUT/WT), compute Δ. WIN = WT clearly worse (the V→E unsatisfied-charge penalty shows up).

## In parallel — the principled method (for the hard substitutions too)
Build **negative design**: ProteinMPNN with an off-target penalty (Baker pMHC-paper approach — score P(target peptide seq | complex) > P(off-target), penalize designs that bind WT), or two-state AF2 (maximize MUT iptm, minimize WT iptm). This is the only method that can discriminate His↔Arg-class swaps; positive design + hotspot structurally cannot (see RUNG-26d autopsy).

## Honest framing
The binder ROUTE is proven (RUNG-26c/d: strong dual-validated pMHC binders on a ≤5% target). The open frontier is SPECIFICITY of a single substitution. BRAF-V600E is the chemically-tractable shot at the first mutation-specific de novo binder; negative design is the general fix. IDH1-R132H stays covered by the internal CRISPR key.
