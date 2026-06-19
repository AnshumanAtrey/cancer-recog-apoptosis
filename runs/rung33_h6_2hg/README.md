# RUNG-33 (H6) — 2-HG oncometabolite recognition RECOVERS IDH1-R132H (the target the binder couldn't)

We **closed** IDH1-R132H as a de novo *binder* target — His↔Arg on the pMHC is undiscriminable (RUNG-26 NULL
across positive + negative design). H6 asks: the mutant IDH1 enzyme floods the cell with the oncometabolite
**D-2-hydroxyglutarate (D-2-HG)** at ~mM (vs ~µM normal) — can that be an *orthogonal, intracellular*
recognition handle for the same target? The crux is the same discrimination problem as the binder route, but
for a small molecule: **tell D-2-HG from its abundant near-twin α-ketoglutarate** (α-KG, the precursor, in all
cells; D-2-HG = α-KG with the C2 ketone reduced to a hydroxyl + a new stereocenter). `scripts/67`.

## Model
A 2-HG-gated recognizer (competitive occupancy; a *specific* sensor treats α-KG as a weak **non-agonist**):
`fire = ([2HG]/K) / (1 + [2HG]/K + [αKG]/(S·K))`, gate CLEAN = fire_mutant ≥ 0.5 AND fire_normal ≤ 0.05.
Concentrations (mM, literature): D-2-HG mutant **10** (5–35, Dang 2009 Nature) / normal **0.005** (~5 µM);
α-KG **0.5** (both).

## Result — FEASIBLE, and EASIER than the binder route
Clean gates span a **wide region** of (affinity × selectivity) — **even at ~1× binding selectivity** (48/64 grid
points). The non-obvious finding:

> **The discrimination is driven by the ~2000× 2-HG concentration differential, NOT by binding selectivity.**
> In normal cells α-KG (0.5 mM) out-competes the trace 2-HG (5 µM) at the sensor but **doesn't fire**; in mutant
> cells 2-HG (10 mM) dominates and fires. The cell's own metabolism **pre-amplifies** the signal.

So unlike the de novo binder (His↔Arg = two equal-abundance bulky residues → intractable), the 2-HG route gets
discrimination **for free** from the concentration gap. The real **requirement** is *chemical specificity*:
α-KG / succinate / L-2-HG must be **non-agonists** (not trigger the gate), not merely weaker binders. That is
**precedented** — D2HGDH is catalytically 2-HG-specific, and engineered 2-HG biosensors exist.

Easiest operating point: K_2HG ≈ 0.001 mM, S ≈ 2× → fire_mut 0.98, fire_norm 0.02.

## What this means for the project
IDH1-R132H now has **two viable recognition routes despite the binder route failing**:
1. **Internal CRISPR key** (DNA-level, RUNG-27d — already built + off-target-clean).
2. **2-HG oncometabolite gate** (this rung — feasible, concentration-differential-driven).
The de novo binder was the *only* route that failed for IDH1; the target was never actually lost.

## Honest residuals (this is feasibility, not a designed sensor — cf. RUNG-25 for RNA)
- **α-KG non-agonist assumption is load-bearing.** A false-agonist α-KG would demand high binding selectivity.
- Needs **D- vs L-2-HG chirality** + rejection of succinate/glutarate/malate (other dicarboxylates) — the full
  specificity panel, not just α-KG.
- ~2% normal-cell leak at the single-gate level → **AND-gate with a 2nd signal** (the project's logic) drives it
  lower.
- α-KG may *drop* in IDH1-mut (consumed to make 2-HG) → would *help* discrimination; modeled conservatively equal.
- Next (like RUNG-25 → ODesign): actual **sensor-molecule design** (ODesign ligand model / aptamer) + the
  cross-metabolite specificity panel + in-cell validation. Recovers *recognition*; coupling to apoptosis = the
  internal-key circuit.

*Result: `feasibility.json`. Updates HYPOTHESIS_LEDGER H6 (OPEN → feasibility CONFIRMED).*
