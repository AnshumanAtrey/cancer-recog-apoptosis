#!/usr/bin/env python3
"""
RUNG-34 — KILL-COUPLING: turn the (validated, leaky) recognizers into an autonomous self-destruct, and
quantify whether it kills mutant cells while sparing normal. This is the "recognition -> death" step that
makes the internal key an actual root-kill, not just a detector.

ARCHITECTURE (modular self-destruct, MHC-free):
  recognizer x N  ->  AND-gate  ->  apoptosis effector
  - recognizers (grounded, per driver):
      * toehold RNA-sensor  (non-wobble missense)         -- RUNG-25
      * allele-specific CRISPR (wobble G>A drivers)        -- RUNG-27d (deep-seed, genome-wide off-target-clean R31b)
      * DhdR D-2-HG gate (IDH1 oncometabolite)             -- RUNG-33 (detuned ~0.1-1mM)
  - AND-gate: a SPLIT death effector (each half gated by one recognizer; both halves needed to reconstitute),
    OR serial toeholds, OR dual-guide-required transcription -> effector level ~ product of recognizer activities.
  - death effector: iCasp9 (rimiducid-inducible caspase-9 -- the CLINICAL CAR-T safety switch -> engineered
    apoptosis-on-demand is proven) / constitutively-active caspase-3 / BAX. Triggers the cell's OWN apoptosis.

WHY AND, not single: each recognizer leaks (mis-fires in some normal cells). A single death switch kills that
leak fraction -> at ~1e11 normal cells, even 1% leak = 1e9 dead normal cells = lethal. AND-gating two
INDEPENDENT recognizers makes the joint false-positive ~ L1*L2, and the cooperative apoptosis threshold
(bistable MOMP/caspase switch, eARM RUNG-1) filters any uniform sub-threshold leak. (RUNG-22 logic, applied to
KILLING instead of escape.)

The model propagates measured per-recognizer discrimination -> AND effector -> P(death), and reports
tumour-kill, normal-cell false-death, and therapeutic index vs (N inputs, per-recognizer leak).
"""
import os, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung34_kill_coupling")

# grounded recognizer discrimination: (mutant ON-activity, normal-cell leak)
RECOGNIZERS = {
    "toehold_RNA":   (0.90, 0.05),   # RUNG-25 ON/OFF ~10-20x
    "allele_CRISPR": (0.90, 0.05),   # RUNG-27d deep-seed single mismatch; WT residual ~few %
    "DhdR_2HG":      (0.98, 0.02),   # RUNG-33 detuned occupancy gate (fire_mut~0.98 / fire_norm~0.02)
}


def p_death(effector, K=0.30, n=5.0):
    """Apoptosis commitment: bistable, cooperative caspase/MOMP threshold (eARM, RUNG-1).
    K = commitment threshold (frac), n = Hill cooperativity (caspase amplification)."""
    e = np.asarray(effector, float)
    return e**n / (K**n + e**n)


def and_effector(activities):
    """Split / serial AND: effector reconstituted only if all recognizers active -> product."""
    out = 1.0
    for a in activities:
        out = out * a
    return out


def selftest():
    assert p_death(1.0) > 0.99 and p_death(0.0) < 1e-6
    assert abs(and_effector([0.9, 0.9]) - 0.81) < 1e-9
    assert p_death(0.81) > p_death(0.0025)   # mutant >> normal(N=2,L=0.05)
    print("[selftest] death threshold + AND-gate math OK")


def main():
    os.makedirs(OUT, exist_ok=True)
    selftest()
    print("\nKILL-COUPLING: AND-gate(N) of leaky recognizers -> apoptosis (K=0.30, n=5)")
    print("per-recognizer: mutant ON ~0.9, normal leak L. tumour-kill | normal false-death | therapeutic index\n")
    rows = []
    MUT_ON = 0.90
    for N in (1, 2, 3):
        for L in (0.01, 0.02, 0.05, 0.10, 0.20):
            mut_eff = and_effector([MUT_ON] * N)
            norm_eff = and_effector([L] * N)
            kill = float(p_death(mut_eff))
            fdeath = float(p_death(norm_eff))
            ti = kill / fdeath if fdeath > 0 else float("inf")
            rows.append({"N": N, "leak": L, "tumour_kill": round(kill, 4),
                         "normal_false_death": fdeath, "therapeutic_index": ti})
            safe = fdeath < 1e-4
            print(f"  N={N} L={L:.2f}: kill={kill:.3f} | normal_death={fdeath:.2e} | TI={ti:.1e}"
                  f"{'  <- SAFE (normal<1e-4)' if safe else ''}")
    # design requirement: smallest N giving normal_death<1e-4 at a realistic leak (5%)
    safe_at_5 = [r for r in rows if r["leak"] == 0.05 and r["normal_false_death"] < 1e-4]
    minN = min((r["N"] for r in safe_at_5), default=None)
    # single-cell kill at the safe design (the trade-off) + bystander rescue
    design = next((r for r in rows if r["N"] == minN and r["leak"] == 0.05), None)
    print(f"\nDESIGN: min inputs for normal_death<1e-4 at 5% leak = N={minN} "
          f"(single-cell kill {design['tumour_kill'] if design else '?'}).")
    verdict = (
        f"FEASIBLE root-kill. A {minN}-input AND-gated self-destruct on 5%-leak recognizers gives normal-cell "
        f"false-death ~{design['normal_false_death']:.0e} (therapeutic index ~{design['therapeutic_index']:.0e}) "
        f"-- killing the cell via its OWN apoptosis only when {minN} independent mutation/metabolite signals "
        f"co-occur. The kill/specificity TRADE-OFF: single-cell tumour-kill drops to ~{design['tumour_kill']:.2f} "
        f"at N={minN} (a cell escapes if not all N fire) -- but the resistance-agnostic BYSTANDER DEATH WAVE "
        f"(RUNG-13/21, cures across 4-13% escapees) clears those, so partial single-cell kill is sufficient at "
        f"the tumour level. Death module precedent: iCasp9 (clinical CAR-T safety switch) = engineered "
        f"apoptosis-on-demand is proven.")
    print(f"\nVERDICT: {verdict}")
    json.dump({"tag": "rung34_kill_coupling",
               "architecture": "recognizer(toehold/allele-CRISPR/DhdR) x N -> split/serial AND -> iCasp9/caspase-3/BAX apoptosis",
               "apoptosis_model": {"threshold_K": 0.30, "hill_n": 5.0, "basis": "bistable MOMP/caspase switch, eARM RUNG-1"},
               "recognizers_grounded": RECOGNIZERS,
               "min_inputs_for_safe_gate_at_5pct_leak": minN,
               "design_point": design, "sweep": rows, "verdict": verdict,
               "residuals": "Assumes recognizer leaks are INDEPENDENT -- correlated failure (a stressed normal "
                            "cell mis-firing several recognizers; shared delivery/cell-state) degrades the AND "
                            "multiplication and is the key risk to test. Also: the death module itself must be "
                            "non-leaky (basal iCasp9 must not trigger); delivery of the multi-component circuit; "
                            "kinetics (does effector cross threshold before the cell adapts?) -- a full eARM ODE "
                            "of the coupled circuit + in-cell validation are the next/wet steps. This quantifies "
                            "the SPECIFICITY logic of the kill-coupling, not a built circuit."},
              open(os.path.join(OUT, "kill_coupling.json"), "w"), indent=2)
    print(f"[saved] {OUT}/kill_coupling.json")


if __name__ == "__main__":
    main()
