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

## SPECIFICITY VERDICT (2026-06-17) — NULL again: binds, doesn't discriminate (forcing the hotspot on the mutation did NOT solve it)

AF2 cross-check (`specificity_af2.json`) of the 13 Protenix-passers vs MUT & WT pMHC:
- **2/13 are AF2 binders** (MUT pae ≤10): **rank_1** (MUT pae **5.12**, plddt 93, ptx_iptm 0.939) and rank_11 (MUT 8.56).
- **0/13 mutant-specific.** Both binders bind **WT as well or slightly better**: rank_1 DISC = **−0.87** (WT 4.25 < MUT 5.12), rank_11 DISC −0.53. Non-binders' DISC is noise (incl. rank_4 which prefers WT, −12.7).

**This is a robust negative, not an AF2 blind spot:** AF2 *confidently* binds rank_1 to BOTH MUT (5.12) and WT (4.25) — it's not "can't resolve the difference," it's "clearly binds both." And it converges with RUNG-26c (Protenix MUT≈WT, Δiptm 0.0002). So across **two design strategies** (free footprint vs hotspot-forced-on-B4) and **three models** (PXDesign AF2-IG, ColabDesign AF2, Protenix): de novo binders to IDH1-R132H/A*01:01 **bind the pMHC but do not discriminate the His↔Arg substitution.**

**Why:** His→Arg are both bulky; a pocket that fits one tends to tolerate the other. **Positive design (PXDesign + a hotspot) only makes the binder *contact* the residue — it does not make the contact *His-specific*.** Discrimination needs **explicit negative design** (design vs MUT *and against* WT — e.g. Baker's ProteinMPNN off-target-penalty pipeline, or two-state AF2), which neither PXDesign nor a hotspot provides. At generate-and-filter scale (65 designs → 2 binders → 0 specific), the discriminating fraction for this subtle swap is ~0.

**Consequence for the arc:** the EXTERNAL key (de novo binder) is **bounded for subtle substitutions** like R132H. The route that *does* discriminate IDH1-R132H is the INTERNAL key — RUNG-27c's allele-specific **CRISPR/DNA** sensor (no His↔Arg ambiguity at the DNA level). This sharpens, not breaks, the thesis: subtle-substitution neoantigens → internal mutation-sensing; the external binder needs negative design or a chemically-bigger substitution (e.g. BRAF V600E, V→E).

**NEXT options (Anshuman's call):** (a) explicit **negative design** (ProteinMPNN off-target / two-state) — the principled fix, more engineering; (b) re-target **BRAF-V600E** (V→E = far bigger chemical change → easier to discriminate; RUNG-27a: real epitope, no natural TCR); (c) accept the honest bound and lean on the internal CRISPR key for IDH1-R132H.

## Ceiling
In-silico; AF2-IG (permissive, initial-guess) + Protenix (stricter) confidence = structural plausibility, NOT measured affinity or specificity; hotspot biases contact, doesn't prove discrimination; mut-vs-WT + proteome off-target + expression/SPR = wet-lab residual. A prioritized candidate pending the specificity gate.

*Source: `summary.csv` (65 designs, AF2-IG + Protenix scores), `task_info.json`, `difficulty_gauge.png`, `top_designs/` (rank_1–3 .cif). Full 8.7 MB tarball reproducible from the webserver job.*
