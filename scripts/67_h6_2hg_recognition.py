#!/usr/bin/env python3
"""
RUNG-33 (H6) — 2-HG oncometabolite recognition gate for IDH1-R132H: feasibility + the selectivity it requires.

WHY: de novo binders CANNOT discriminate IDH1-R132H His<->Arg on the pMHC (RUNG-26 NULL) -> we CLOSED the
external-binder route for IDH1. But mutant IDH1 floods the cell with the oncometabolite D-2-hydroxyglutarate
(D-2-HG) at ~mM, vs ~uM in normal cells. So 2-HG is an ORTHOGONAL, intracellular recognition handle for the
same target -- IF a sensor can tell D-2-HG from its abundant near-twin alpha-ketoglutarate (alpha-KG, the
precursor, present in ALL cells; D-2-HG = alpha-KG with the C2 ketone reduced to a hydroxyl + a new stereocenter).

This is the SAME discrimination problem as the binder route (tell two similar molecules apart) but for a small
metabolite, where the CONCENTRATION DIFFERENTIAL helps. We model a 2-HG-gated recognizer and ask: what sensor
affinity (K_2HG) + selectivity (S = K_alphaKG / K_2HG) makes the gate fire in mutant cells and stay OFF in
normal? Like RUNG-25 (RNA-sensor thermodynamic feasibility), this is a feasibility/requirement calc, not a
designed molecule -- the honest first rung of H6.

Model (competitive occupancy; a SPECIFIC sensor treats alpha-KG as a weak competitor, NOT a false agonist):
  fire(cell) = ([2HG]/K_2HG) / (1 + [2HG]/K_2HG + [aKG]/K_aKG),   K_aKG = S * K_2HG
  CLEAN GATE = fire_mutant >= 0.5 (ON) AND fire_normal <= 0.05 (<5% leak in normal tissue).
"""
import os, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung33_h6_2hg")

# concentrations in mM (intracellular, literature-grounded)
HG_MUT = 10.0     # D-2-HG in IDH1-mutant tumours: 5-35 mM (Dang 2009 Nature; ~10-30 mM) -> 10 conservative
HG_NORM = 0.005   # D-2-HG in normal/IDH1-WT: low uM (~1-5 uM) -> 5 uM
AKG = 0.5         # alpha-KG: ~0.1-1 mM in most cells, roughly constant -> 0.5 mM (both states)
FIRE_ON, LEAK_MAX = 0.5, 0.05


def fire(K_hg, S, hg, akg=AKG):
    """Fraction of sensor productively bound by 2-HG; alpha-KG competes (non-agonist)."""
    K_akg = S * K_hg
    n = hg / K_hg
    return n / (1.0 + n + akg / K_akg)


def selftest():
    # no alpha-KG, saturating 2-HG -> ~full occupancy; trace 2-HG -> ~0
    assert fire(0.1, 10, 100.0, 0.0) > 0.99
    assert fire(0.1, 10, 1e-6, 0.0) < 0.01
    # alpha-KG with S=1 (no selectivity) at high conc competes down the signal
    assert fire(0.1, 1, 0.5, 10.0) < fire(0.1, 1, 0.5, 0.0)
    print("[selftest] occupancy model OK")


def main():
    os.makedirs(OUT, exist_ok=True)
    selftest()
    K_grid = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]     # mM, sensor 2-HG affinity
    S_grid = [1, 2, 3, 5, 10, 30, 100, 300]                    # selectivity (2-HG over alpha-KG)
    print(f"\n2-HG: mutant {HG_MUT} mM / normal {HG_NORM*1000:.0f} uM | alpha-KG {AKG} mM (both)")
    print("CLEAN gate (fire_mut>=0.5 AND fire_norm<=0.05) over (K_2HG mM x selectivity S):\n")
    header = "  K_2HG\\S  " + "".join(f"{s:>6}" for s in S_grid)
    print(header)
    feasible = []
    for K in K_grid:
        row = f"  {K:>6.3f}  "
        for S in S_grid:
            fm, fn = fire(K, S, HG_MUT), fire(K, S, HG_NORM)
            ok = fm >= FIRE_ON and fn <= LEAK_MAX
            row += f"{'  GATE' if ok else '    . '}"
            if ok:
                feasible.append({"K_2HG_mM": K, "selectivity": S, "fire_mut": round(fm, 3), "fire_norm": round(fn, 4)})
        print(row)
    # minimum selectivity that yields ANY clean gate, and the easiest operating point
    min_S = min((f["selectivity"] for f in feasible), default=None)
    best = max(feasible, key=lambda f: f["fire_mut"] - f["fire_norm"]) if feasible else None
    print(f"\nfeasible (K,S) gates: {len(feasible)} | min selectivity for a clean gate: {min_S}x")
    if best:
        print(f"easiest operating point: K_2HG={best['K_2HG_mM']} mM, S={best['selectivity']}x "
              f"-> fire_mut={best['fire_mut']}, fire_norm={best['fire_norm']}")
    verdict = (f"FEASIBLE, and EASIER than the binder route. Clean gates span a wide (affinity x selectivity) "
               f"region -- even at LOW binding selectivity (down to ~{min_S}x). KEY INSIGHT: the discrimination "
               f"is driven by the ~2000x 2-HG CONCENTRATION DIFFERENTIAL (mutant ~mM vs normal ~uM), NOT by "
               f"binding selectivity -- in normal cells alpha-KG out-competes the trace 2-HG but does not FIRE; in "
               f"mutant cells 2-HG dominates. So the cell's own metabolism pre-amplifies the signal. The real "
               f"REQUIREMENT is CHEMICAL SPECIFICITY (alpha-KG / succinate / L-2-HG must be NON-agonists, not just "
               f"weaker binders) -- which is precedented (D2HGDH is catalytically 2-HG-specific; engineered 2-HG "
               f"biosensors exist). Unlike the de novo binder (His<->Arg = two equal-abundance bulky residues, "
               f"intractable, RUNG-26), 2-HG gets discrimination FOR FREE from the concentration gap. -> 2-HG "
               f"RECOVERS IDH1-R132H recognition.") if feasible else "NOT feasible at modeled concentrations."
    print(f"\nVERDICT: {verdict}")
    json.dump({"tag": "rung33_h6_2hg_recognition_feasibility",
               "concentrations_mM": {"2HG_mutant": HG_MUT, "2HG_normal": HG_NORM, "alphaKG": AKG},
               "concentration_sources": "D-2-HG IDH1-mut 5-35 mM (Dang 2009 Nature); normal low-uM; alpha-KG ~0.1-1 mM.",
               "gate_criteria": {"fire_on": FIRE_ON, "leak_max": LEAK_MAX},
               "min_selectivity_for_clean_gate": min_S, "n_feasible_points": len(feasible),
               "easiest_operating_point": best, "feasible_points": feasible,
               "verdict": verdict,
               "residuals": "Feasibility/occupancy model, not a designed sensor (cf. RUNG-25 for RNA). Assumes "
                            "alpha-KG is a competitive NON-agonist (a SPECIFIC sensor); a false-agonist alpha-KG "
                            "would need much higher selectivity. Also needs D- vs L-2-HG chirality + succinate/"
                            "glutarate cross-talk + actual sensor design (ODesign ligand model / aptamer) + in-cell "
                            "validation. Recovers IDH1 RECOGNITION; coupling to apoptosis is the internal-key circuit."},
              open(os.path.join(OUT, "feasibility.json"), "w"), indent=2)
    print(f"[saved] {OUT}/feasibility.json")


if __name__ == "__main__":
    main()
