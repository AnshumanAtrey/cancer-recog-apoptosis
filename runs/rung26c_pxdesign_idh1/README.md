# RUNG-26c — PXDesign de novo binder for IDH1-R132H/HLA-A*01:01 (the first CREDIBLE pMHC binder of the arc)

**Tool:** PXDesign (ByteDance Seed; bioRxiv 2025.08.15.670450) via the free Protenix webserver `protenix-server.com`, **Preview mode**, job `f64b2075`. PXDesign-d (diffusion) backbones → ProteinMPNN sequences → AF2-IG confidence filtering.

**Target:** IDH1-R132H neoantigen on HLA-A*01:01 — the **clean** handle RUNG-12 certified per-cell-safe and RUNG-27a showed has **no natural class-I TCR** (so de novo design is the *only* recognition route, and exactly the neoantigen-with-no-experimental-structure case the Baker *Science* 2025 pMHC paper validated). pMHC = groove `SHSMRYFFT…LQRTD` + mutant peptide `IIGHHAYGDQY` (mutated His at p4).

## Result — 5 designs, 3 pass AF2-IG-easy; rank_1 is a genuine candidate

| rank | iptm | interface-pAE | plddt | monomer-plddt | bound↔unbound RMSD | verdict |
|---|---|---|---|---|---|---|
| **1** | **0.81** | **5.75** | 0.93 | 0.98 | **0.25** | passes AF2-IG + strict — **credible binder** |
| 2 | 0.74 | 9.27 | 0.94 | 0.96 | 8.37 | passes easy but floppy (large conf. change) |
| 3 | 0.69 | 10.18 | 0.91 | 0.94 | 13.27 | passes easy but very floppy |
| 4 | 0.56 | 12.49 | 0.89 | — | 1.69 | fails |
| 5 | 0.43 | 15.34 | 0.84 | — | 1.61 | fails |

rank_1 binder (80 aa): `SAEEELLAAEARASELEVRVRRLALEQGDEEALRRLDDIGTETRERLNAARAAGASTEERLAIVREALARLEALLAEVEA`.
interface-pAE 5.75 is comfortably inside the ≤10 binder bar; for contrast AfDesign (RUNG-26) gave 0/6 and RFdiffusion (RUNG-26b) 0/50 at pae 25–27. **This is the best result of the entire binder arc.**

## Contact analysis (rank_1.cif, parsed locally — see `rank1_contact_analysis.json`)
The binder **arches over the pMHC** (Baker-style): 16 MHC contacts + **13 peptide contacts**, and it **does touch the mutated His at p4 (closest 2.81 Å)**. BUT its footprint **centers on the conserved C-terminal peptide residues** (p7 Tyr = 9 contacts; p5–p6), with the mutation at the footprint *edge* (1 contact).

## Difficulty (see `difficulty_gauge.png`)
PXDesign rates this target at the **≤5% AF2-IG-easy passing rate — its hardest tier** (with IL17A/TNF/SARS-RBD). So a credible design from a Preview batch is a real result, but most attempts on this target fail → a full (non-Preview) campaign needs many designs.

## HONEST STATUS — specificity is UNTESTED here (the make-or-break)
PXDesign does **no negative/specificity design**. rank_1 binds the *mutant* pMHC well, but our thesis requires **binds MUT, not WT**. Because the interface centers on residues identical between MUT (`IIGHHAYGDQY`) and WT (`IIGRHAYGDQY`), discrimination is **plausible but not proven**.
- **Pending test (submitted 2026-06-17):** fold rank_1 vs WT pMHC on the Protenix webserver (jobs `rank1_MUT` / `rank1_WT`, identical settings, seed-matched) and compare iptm/interface-pAE.
- **WIN** = WT clearly worse (iptm ↓≥0.15 or ipAE ↑≥3). **NULL** = WT≈MUT → re-run PXDesign with the hotspot pinned on the mutated His (p4) + a full batch.
- **Honest prior (from the contact map):** discrimination likely weak → expect to need the hotspot-on-p4 redesign.

## Ceiling
In-silico, Preview batch; iptm/interface-pAE = AF2 structural confidence, NOT measured affinity; mut-vs-WT specificity + proteome-wide off-target + expression/SPR = the wet-lab residual. A prioritised candidate, not a validated binder.

*Source archive: `design_outputs/` (PXDesign output, bit-for-bit), `difficulty_gauge.png`, `rank1_contact_analysis.json`.*
