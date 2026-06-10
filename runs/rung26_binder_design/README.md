# RUNG-26 — de novo mutant-specific binder, target IDH1-R132H / HLA-A*01:01

**Provisional negative (honest scope stated up front).** AfDesign (AF2-hallucination, ColabDesign v1.1.1)
was run on the clean, *exposed* (RUNG-12 E=0.7) neoantigen pMHC IDH1-R132H/A\*01:01, hotspot forced on the
mutated His. Folded with Boltz (proven RUNG-20), designed + scored on a free Colab T4.

## Result: 0 / 6 credible binders → 0 mutant-specific

| design | mut pae_interaction *(bar ≤10)* | binder pLDDT *(bar ≥80)* | binder? |
|---|---|---|---|
| design_0002 | 14.18 | 68.2 | no (closest) |
| design_0005 | 15.08 | 57.4 | no |
| design_0001 | 17.07 | 51.2 | no |
| design_0000 | 22.35 | 48.3 | no |
| design_0003 | 22.86 | 48.5 | no |
| design_0004 | 25.67 | 57.4 | no |

Filters (Bennett 2023): `pae_interaction ≤ 10` AND `binder_plddt ≥ 80`; mutant-specific needs `pae_WT − pae_MUT ≥ 3`.
None crossed the binder bar, so none were WT-folded (discrimination N/A).

## What this is — and what it is NOT (rule 5)

**This is a PROVISIONAL bound, not a verdict.** Three reasons it is thin:
1. **N = 6.** The 4-hour batch was killed by the free-Colab session at design 6 (it is resumable — re-running
   Cell 5 accumulates more). Six samples of a stochastic generator cannot prove "no binder exists."
2. **Single binder length (70 aa).** The spec range is 60–100; only 70 was sampled.
3. **Single tool.** AfDesign-hallucination is the *weakest* de novo method for pMHC — a tiny, mostly-MHC-surface
   target. The field-proven approach for hard flat targets is **RFdiffusion** partial-diffusion onto the hotspot
   → ProteinMPNN → AF2 filter (Bennett/Watson), which AfDesign does not replicate.

**The honest read:** *vanilla AfDesign cannot make a credible binder to this pMHC* — "wrong tool," **not** "binding
impossible." The trend (best pae 14.2, ~40% over bar; pLDDT plateauing 48–68) is not approaching the threshold,
consistent with the literature that de novo pMHC binders need RFdiffusion + heavy oversampling + wet-lab screening
(<1% in-silico success even there).

## Why this target matters (RUNG-27a)
Reality-grounding (VDJdb+IEDB) shows IDH1-R132H has **no natural class-I TCR on record** → de novo design is the
*only* route to a recognition molecule here, which is exactly why it was chosen. The failure is the method, not the
target's worth.

## Next levers (in EV order)
1. **BRAF-V600E** (parallel, account #2, `notebooks/binder_design_braf_colab.ipynb`) — different groove, bigger
   V→E change; controlled test of whether *any* clean pMHC is AfDesign-designable.
2. **RFdiffusion** (RUNG-26b, GPU) — the field lever for hard pMHC. Heavier; still a low-yield + wet-lab frontier.
3. **The orthogonal routes don't depend on cracking this:** KRAS-G12D already has a natural TCR (RUNG-27a), and
   the autonomous mutation-sensing circuit (RUNG-27b) opens the MHC-free route entirely.

*Run log: `runs/logs/rung26_binder_design_20260609T225744Z_1c2874a.log`. Result JSON: `rung26_binder_design.json`.*
