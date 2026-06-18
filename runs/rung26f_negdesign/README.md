# RUNG-26f — ProteinMPNN two-state **negative design**: the principled fix RUNG-26d's autopsy demanded

> **FINAL VERDICT (2026-06-18, GPU AF2 confirm, `af2_confirm_idh1.json`): NULL — 0/12 mutant-specific.**
> Ran the full pipeline on a clean T4: ProteinMPNN two-state generation on the 3 PXDesign dual-passer
> backbones → 12 top-dscore candidates → AF2 two-state confirm. **All 12 are no-bind (MUT pae 22-25 ≫ 10)
> and every discrimination is ≤ 0** (several bind WT *better*, e.g. −4.06). So across **positive design
> (NULL ×3 models)** AND the **principled negative design (NULL)**, IDH1-R132H's His↔Arg is **intractable
> for a de novo binder** — the substitution is too subtle/buried to grip. IDH1-R132H stays covered by the
> **internal CRISPR key** (RUNG-27d). This is exactly why the external-binder effort pivoted to a
> **presentation-flip** target (PIK3CA-E545K, RUNG-28/29) where the binder needn't discriminate at all.
> The method itself is validated and reusable (the scorer doesn't hallucinate discrimination; see below).

RUNG-26c/d proved de novo binders to IDH1-R132H/HLA-A\*01:01 **bind but don't discriminate** His↔Arg
(NULL across AF2-IG, ColabDesign-AF2, Protenix). The autopsy showed *why*: positive design + a hotspot
makes the binder **contact** the mutation but cannot make binding **depend** on it. The only method that
can is **negative design** — score/optimize the binder against the MUT pMHC **and against the WT pMHC**,
keeping only sequences the model finds much less likely on WT. This rung builds + validates that method on
M2/CPU and ships the GPU notebook that runs it at scale and AF2-confirms it.

## Mechanic (`scripts/57_negdesign_proteinmpnn.py`)
Same backbone for both states (the whole point). The peptide's mutated position is `HIS` in the MUT
target / `ARG` in the WT target (IDH1 R132H: WT=Arg → mut=His, at p4). ProteinMPNN uses backbone atoms for
geometry and the **fixed (visible) target identities** as decoding context, so the binder chain's
conditional NLL changes with that one residue. Design chain C (binder), fix A (groove) + B (peptide):
- `score = mean NLL over binder residues` (lower = ProteinMPNN finds the binder more likely in that context)
- **`dscore = NLL_wt − NLL_mut`  → positive means the binder prefers the MUT context = discriminating.**

## Result 1 — VALIDATION GATE (`validate.json`): the scorer does NOT hallucinate discrimination
Ran the two-state scorer (16 decoding orders) on the existing **non-specific** PXDesign binders (rank_1/2/3,
already NULL under three models). Gate = |dscore| ≈ 0.

| design | NLL_mut | NLL_wt | **dscore (WT−MUT)** | noise sd |
|---|---|---|---|---|
| rank_1 | 1.084 | 1.079 | **−0.005** | 0.031 |
| rank_2 | 1.151 | 1.154 | **+0.003** | 0.031 |
| rank_3 | 1.185 | 1.181 | **−0.004** | 0.048 |

Every dscore is **smaller than its own order-noise sd** → no discrimination invented where none exists.
This is a **4th independent angle** convergently confirming the RUNG-26d NULL, and it **calibrates the noise
floor (~0.03–0.05 NLL units)** that a real negative-design dscore must clear.

## Result 2 — GENERATION path validated (`generate_demo.json`): works; IDH1 His↔Arg is hard even here
ProteinMPNN-sampled new binder sequences conditioned on the MUT context, two-state-scored each with
**paired** decoding orders (cancels most order-noise in the difference), ranked by dscore:

| sample size | best dscore | mean | interpretation |
|---|---|---|---|
| N=8  (temp 0.20) | **+0.0006** | −0.0048 | ≈ 0 |
| N=24 (temp 0.25) | **+0.0235** | −0.0024 | weak, right at the noise floor |

The two-state tail **grows with oversampling** but stays **marginal** on IDH1 — consistent with the
autopsy (His and Arg both bulky, both H-bond → the swap is energetically marginal). Sampled sequences are
clean binder-like helices, so the *mechanic* is sound; the *target* is intrinsically resistant.
**MPNN dscore is a cheap screen, not proof** — the real gate is AF2 (below).

## The GPU run (`notebooks/negdesign_two_state_colab.ipynb`)
Reuses the **M2-validated** `scripts/57` functions verbatim (so the GPU run can't drift from what we gated):
1. **Generate** two-state-selected binders on every available backbone — the 3 repo PXDesign complexes
   (always present) + the RFdiffusion backbones on Drive if found — GPU-fast sampling, ranked by dscore.
2. **AF2-confirm** the top candidates: fold each vs MUT and WT pMHC with the same `score()` used in all
   prior RUNG-26 runs. **Gate:** `mutant_specific = pae_mut ≤ 10 AND binder_plddt ≥ 80 AND (pae_wt − pae_mut) ≥ 3`.
3. **Repoints to BRAF-V600E** (V→E, the chemically tractable case) by editing the `CFG` block once the BRAF
   MUT/WT folds land — a buried unsatisfied charge on WT-Val gives negative design a real handle, unlike His↔Arg.

## Honest framing & residual
- This rung delivers a **built + validated method**, not a confirmed specific binder. The IDH1 generation
  signal is marginal and **unconfirmed by AF2** — no specificity is claimed.
- A positive MPNN dscore is **necessary, not sufficient**; only the AF2 two-state gate (and ideally a
  seed-matched Protenix webserver MUT/WT cross-check) supports a discrimination claim.
- In-silico likelihood/structure ≠ measured affinity/specificity; proteome off-target, expression, and
  SPR/cellular assays remain the wet-lab residual.
- IDH1-R132H stays covered by the **internal CRISPR key** (RUNG-27d: 7/7 wobble drivers DNA-addressable);
  the external specific-binder artifact is being pursued where it's tractable (BRAF).

*Files: `validate.json`, `generate_demo.json`, `notebooks/negdesign_two_state_colab.ipynb`,
`scripts/57_negdesign_proteinmpnn.py` (selftest + `--generate`), `scripts/_build_negdesign_nb.py` (notebook builder).
ProteinMPNN cloned to `.tools/` (gitignored, not an artifact).*
