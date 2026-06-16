# RUNG-26b — RFdiffusion + ProteinMPNN binder attempt for IDH1-R132H/HLA-A*01:01 (provisional NEGATIVE)

De novo binder to the same clean neoantigen pMHC as RUNG-26c, via the field-standard RFdiffusion (hotspot on the up-facing mutated peptide residue) → ProteinMPNN → AF2 scoring (`scripts/52`, `notebooks/binder_rfdiff_design_colab.ipynb`). Runs on a free Colab T4 (16 GB) — RFdiffusion+ProteinMPNN fit there (unlike PXDesign/BindCraft which need ≥32 GB).

## Two batches

**v1 (the `design_*/metrics.json` in this directory) — DIVERSITY BUG, kept for the record.**
73 designs but **1 unique sequence** (all pae_interaction 26.3, plddt 94) — caused by calling RFdiffusion with `num_designs=1` inside a loop (fixed seed → identical backbone). Not a real batch; archived only as the bug record. `rung26_binder_design.json` is this buggy batch's summary.

**v2 (`colab_scoring_20260616_8of50/`) — diversity fixed, real NEGATIVE.**
Re-run with one `num_designs=50` call → **50 unique** backbones/sequences (RFdiffusion ran on the T4 ~195 s/backbone, ~60× faster than CPU — confirmed GPU). Scoring (AF2/ColabDesign) reached **8/50 before the session timebox** — all 8 are clear non-binders: interface-pAE **25–27** (bar ≤10), nowhere near the threshold. **0 binders.**

## Honest read
0 binders, but this is the **expected field rate**, not "impossible": PXDesign rates this target ≤5% passing (RUNG-26c difficulty gauge), and at 50 designs (1 ProteinMPNN seq/backbone) the expected number of passers is ~0–2. The batch is too small, and the ProteinMPNN sequences here are notably low-complexity (poly-Ala-ish) → weak interfaces. **The same target was cracked by PXDesign (RUNG-26c, rank_1 iptm 0.81 / ipAE 5.75)** — so the route works; RFdiffusion at this scale/sequence-design just didn't find it.

## Ceiling
In-silico; AF2 interface-pAE = standard de novo filter (Bennett 2023), confidence not affinity; specificity (mut-vs-WT) not even reached (no binder passed step 1); v2 scoring truncated at 8/50 by free-Colab timebox.
