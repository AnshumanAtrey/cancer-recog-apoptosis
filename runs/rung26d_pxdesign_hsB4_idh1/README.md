# RUNG-26d — PXDesign Extended, hotspot FORCED on the mutated His (B4) for IDH1-R132H/A*01:01

The discrimination-forced redesign that RUNG-26c's specificity NULL pointed to. RUNG-26c's rank_1 bound the pMHC but read the **conserved** p5–p7 core (MUT≈WT, Δiptm 0.0002) → not mutation-specific. Fix: pin the design hotspot on **chain B residue 4** (the R132H His) so every binder *must* contact the mutation.

**Run:** PXDesign **Extended Mode** (protenix-server.com), target = `pmhc_mut_cropped.pdb` (peptide + 10 Å groove, pre-cropped because the in-browser crop is a manual viewer step), hotspot **B4**, binder length 100. Extended gives **both AF2-IG and Protenix** filter scores (vs Preview's AF2-IG only).

## Result — 65 designs; 1 passes BOTH filters

| filter | passers / 65 |
|---|---|
| AF2-IG-success | 8 |
| Protenix-success | 13 |
| **BOTH (AF2-IG ∧ Protenix)** | **1 — rank_1** |

**rank_1** (100 aa): af2_ipAE **5.73**, af2_iptm 0.83, **ptx_iptm 0.939**, **ptx_iptm_binder 0.894**, bound↔unbound RMSD **0.35** (pre-organized). Dual-model–validated binder — stronger than RUNG-26c's rank_1 (AF2-IG only). Sequence in `summary.csv`; structure in `top_designs/rank_1.cif` (rank_2, rank_3 = next-best Protenix passers).

**Difficulty (`difficulty_gauge.png`): still ≤5% (hardest tier)** on both AF2-IG and Protenix passing rates even with the hotspot set — so 1 dual-passer out of 65 is a genuine result, not noise.

## OPEN TEST — specificity (the whole point), NOT yet done
The B4 hotspot **biases** every design to contact the mutated His, but passing the binder filters ≠ discriminating MUT from WT — a binder can contact B4 *and* the conserved core and still bind both. **Pending (next, GPU):** score rank_1 (and the other dual/Protenix passers) vs **MUT and WT** pMHC on AF2 (`binder_specificity_rank1_colab`, repoint to this summary.csv) **and** Protenix webserver (rank_1 + MUT/WT, seed-matched), compute Δiptm / Δ-ipAE.
- **WIN** = WT clearly worse than MUT (the B4 hotspot worked → first mutation-specific de novo binder for a no-natural-TCR neoantigen).
- **NULL** = MUT≈WT again → the binder reads B4 *plus* conserved residues; tighten further (hotspot B4 only + crop tighter, or shorter binder, or penalize conserved-core contact).

## Ceiling
In-silico; AF2-IG (permissive, initial-guess) + Protenix (stricter) confidence = structural plausibility, NOT measured affinity or specificity; hotspot biases contact, doesn't prove discrimination; mut-vs-WT + proteome off-target + expression/SPR = wet-lab residual. A prioritized candidate pending the specificity gate.

*Source: `summary.csv` (65 designs, AF2-IG + Protenix scores), `task_info.json`, `difficulty_gauge.png`, `top_designs/` (rank_1–3 .cif). Full 8.7 MB tarball reproducible from the webserver job.*
