# RUNG 2 — DR5 clustering-geometry consistency check (methodology + honest result)

**One line:** RUNG 2 asks whether the *geometry* of how a DR5 agonist clusters the receptor predicts
its potency *beyond what valency (arm-count) already tells you* — and it is wired so the default,
most-likely, fully-acceptable answer is the honest **"no, it's a valency lookup table."** It does
**not** prove agonism. (Pre-registration: `runs/rung2_clustering/manifest_PREREG.yaml`.)

This rung was hardened by an adversarial design review **before** any code was written. The three
attacks — and the fixes — are baked into the pipeline.

---

## The result (what actually came out)

| finding | value | meaning |
|---|---|---|
| **Spearman(potency_rank, valency)** | **0.94** | valency alone explains ~88% of the real ladder — *that's the law* |
| max achievable global Δρ | 0.06 | the most geometry could possibly add globally |
| permutation-null 95th pct of Δρ | 0.67 | what noise produces at n=9 → **the global test is unreachable** (reported sign-only) |
| **PRIMARY within-valency** (4 bivalents) | 5/6 concordant, **exact p = 0.167** | directionally consistent but **NOT significant** (best case possible is p=0.042) |
| **NEGATIVE CONTROL** | scrambled features reach the same in **18%** of permutations | the apparent within-valency signal is **reproduced by noise → headline VOID** |
| **VERDICT** | **VALENCY LOOKUP TABLE** | geometry did not demonstrably add trustworthy signal beyond valency |
| **POSITIVE side-result** | Spearman(sim firing, potency) = **0.80** | the percolation sim *reproduces* the ladder — it **explains** the valency law mechanistically |

**Plain version:** across the real DR5 drug ladder, how many arms a drug has already predicts ~90% of
how well it works. Once you know that, the *shape* of the clustering adds nothing we can trust at this
sample size — and the one apparent flicker of "shape signal" is reproduced when we scramble the labels,
which proves it was noise. That is an honest negative, and it is the expected outcome. The one thing
that *is* defensible and positive: a from-scratch percolation simulation reproduces the valency→potency
ladder (0.80), i.e. it gives a physical *reason* why valency works (bigger arm-count → bigger receptor
clusters → crosses the percolation threshold that nucleates the death complex).

---

## How it's built (file by file)

1. **`data/dr5_agonists/dr5_agonist_ladder.csv`** — Table A, the frozen n=9 potency answer key
   (TAS266 split out to Table B because its rank was a *toxicity* phenotype, not potency).
2. **`scripts/14_blind_labeler.py`** — assigns geometry inputs (arm-reach, epitope) from **format/
   structure only**, after *dropping* the potency column (blind by construction; `label_provenance.md`).
3. **`scripts/13_clustering_sim.py`** — a real 2D percolation/union-find clustering simulation: a
   valency-N agonist crosslinks DR5 receptors; components ≥ the DISC-nucleation size are firing-competent.
   Higher valency → bigger clusters → percolation. Outputs the geometry features `g1…g6`.
4. **`scripts/12_clustering_rank.py`** — the honest machinery: leakage-free leave-one-out (scaling +
   isotonic refit *inside every fold*), the printed ceiling, the within-valency exact-p PRIMARY test,
   the scrambled-feature negative control, the collinearity guard, and the verdict logic that **defaults
   to the honest failure report**.

Reproduce: `python scripts/14_blind_labeler.py && python scripts/13_clustering_sim.py && python scripts/12_clustering_rank.py`
(CPU, ~12 s). Colab mirror: `notebooks/rung2_clustering_colab.ipynb` (adds the Boltz-2 resolution probe).

---

## The three adversarial fixes (why you can trust the negative)

1. **"It's rigged to fail" (true, and now stated up front).** Because valency explains ~88%, the global
   "does geometry beat valency" test is mathematically unreachable at n=9. We *print the ceiling before
   fitting* and demote Δρ to a sign, moving the real test to the only place valency is blind — the four
   same-valency bivalent antibodies.
2. **"The geometry comes from an unreliable oracle" (true).** Our own committed evidence shows Boltz-2
   can't tell a binder from a fake on this exact receptor. So geometry inputs are literature/format-derived
   here, the within-valency inputs are treated as oracle-noise-dominated by default, and the Colab notebook
   carries the multi-seed Boltz-2 probe to confirm it.
3. **"At n=9 it could be a fluke / leakage" (true).** Fixed with: leakage-free LOO, a **scrambled-feature
   negative control** as the decisive gate (it *did* reproduce the signal → headline void), exact
   small-sample p with the honest power floor printed, no bootstrap, and an explicit blinding statement.

**HARD RULE (asserted in code):** the RUNG-2 clustering score is **never** multiplied or combined with the
RUNG-1 apoptosis-kinetics gate-strength into a single "efficacy" number. They are two separate axes with
separate ceilings and no measured map between them — two coordinates, never a product.

---

## Your wet-lab directive, answered honestly (per assay)

You asked: whenever a wet lab is needed, find a laptop/Colab alternative — but never fake it. Here is the
honest per-assay map. **Bottom line first: no in-silico proxy here substitutes for any wet-lab assay. The
sim changes only *which constructs you prioritize*, never *whether* the irreducible experiment must run,
and never *how many* experiments are needed.** Cloud-lab cost/access figures below are **UNVERIFIED** —
confirm before relying on them.

| wet-lab assay | best laptop/Colab proxy | faithfulness | irreducible residual | cheapest real path (student) |
|---|---|---|---|---|
| **Crosslink-dependence** (±anti-Fc / FcγR) — the axis RUNG 2 targets | clustering sim scores geometric *capability* only | **low–moderate (capability only)** | membrane reality (lateral diffusion, density, glycocalyx, actin corral) decides if lattice geometry actually fires; shares the caspase-8 residual | caspase/viability ±anti-Fc at a low-cost India CRO (e.g. Bioneeds) — *cost UNVERIFIED* |
| **Caspase-8 activation** (Caspase-Glo 8) — THE crux | **none** (capability hypothesis only) | **low** | the cooperative DISC enzymatic threshold — literally the definition of agonism; unskippable | Emerald Cloud Lab / Strateos credits, or India CRO — *access UNVERIFIED* |
| Caspase-3/7 + apoptosis commitment | **RUNG 1 (done)** PySB/EARM | **high *only given a measured caspase-8 input*** | the per-construct caspase-8 input itself (the crux) is unmeasured; never combine with the sim | iXCells / Bioneeds fee-for-service |
| Cell viability / killing (CellTiter-Glo EC50/Emax) | **none** for the number (rank-order hypothesis only) | low | emergent EC50/Emax (receptor density, internalization, FLIP, heterogeneity) | CellTiter-Glo dose-response at Strateos / ECL / India CRO |
| DISC co-IP (FADD + caspase-8) | AlphaFold/Boltz static plausibility | **low** (builds interfaces even for non-binders — proven by our lysozyme control) | actual cellular recruitment/stoichiometry vs cFLIP | university proteomics core (IIT/IISc/NCBS/inStem) with a sponsoring PI |
| Super-res cluster imaging (dSTORM) | sim lattice as an *expected picture* | **low** (hypothesis-grade) | absolute live-membrane nanocluster size; size ≠ function | university SMLM core via a collaborator |

Front-end for all paid assays: community-bio / iGEM / biofoundry bench access (~\$100/mo, *UNVERIFIED*) to
express/purify the construct, then hand off. The honest chain that must run regardless of the sim:
**Caspase-Glo 8 firing → CellTiter-Glo EC50/Emax on DR5+ lines → dSTORM/co-IP for mechanism.**

---

## What RUNG 2 does NOT claim

- It is **not** a geometry predictor (the default verdict is a valency lookup table).
- It does **not** prove any construct fires caspase-8 (the agonism crux is wet-lab).
- It makes **no** generalization claim — every molecule is in the calibration set (zero out-of-class
  validation), so this is retrospective description of an n=9 panel only.
- The cross-assay nature of the rank (different cell lines/conditions) is the dominant uncertainty; the
  tigatuzumab crosslink-independence label is **contested**, and there is no head-to-head assay to settle it.

**The single most likely and fully acceptable outcome — realized here — is the honest failure report.**
That is the result surviving scrutiny instead of overclaiming. See `EVIDENCE_AND_HANDOFF.md` for the
project-wide evidence ceiling and the wet-lab handoff.
