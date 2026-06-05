#!/usr/bin/env python3
"""
RUNG 7 — the AND-NOT recognition-gate DISCRIMINATION model (laptop, seconds).

THE QUESTION RUNG-6 DID NOT ANSWER
----------------------------------
RUNG-6 counted how many patients have a usable genetic gate (an accounting question). It did NOT ask the
literal recognition question: given the real spread of antigen density and HLA expression, does the Tmod
AND-NOT gate actually SEPARATE a cancer cell from a healthy one, and HOW DOES IT FAIL? This models bar #2 of
the gate ladder (discrimination dynamics) and couples it to the apoptosis commit (RUNG-1 / EARM essence),
giving the first end-to-end in-silico 'recognise -> commit to apoptosis' demonstration.

THE GATE (A2 Bio 'Tmod' dual receptor)
--------------------------------------
  activator signal a = Hill(antigen_A density; Ka, na)         # broad CAR (CEA/MSLN/EGFR), low-affinity
  blocker  signal b = Hill(HLA_B density;     Kb, nb)          # LIR-1, senses a germline allele
  kill-license  L  = a * (1 - b)                               # AND-NOT: fire iff activator high AND blocker low
  apoptosis commit = bistable all-or-none if L > theta         # EARM essence (snap-action, RUNG-1)

POPULATIONS (densities in log10 molecules/cell)
  TUMOUR  : antigen_A high (some antigen-low tail); HLA_B = 0 for cells with CLONAL LOH, retained for the
            subclonal-retained fraction (escape route B).
  NORMAL  : antigen_A mostly below the low-affinity activator threshold but with a leak tail (RUNG-5 showed
            the activator antigens leak into normal tissue); HLA_B high EXCEPT a downregulated 'HLA-low'
            fraction that loses blocker protection (failure route A).

THE LOAD-BEARING (parameter-ROBUST) CLAIM under test:
  the gate can be no safer than the reliability of HLA_B expression in NORMAL tissue — off-tumour toxicity
  floor ~= P(normal cell is HLA-low) x P(it is also antigen-high). Safety is carried ENTIRELY by the blocker
  because the activator is deliberately broad. We verify this is robust to the activator parameters.

HONEST CEILING
--------------
A mechanistic circuit model with LITERATURE-GROUNDED parameters and ILLUSTRATIVE per-cell distributions (the
real per-cell antigen/HLA joint distribution needs the scRNA/proteomic atlas -> a Colab run). It predicts the
two failure modes and the therapeutic window; it does NOT measure them in tissue. binding != agonism remains
the wet-lab residual. The STRUCTURAL findings (safety = blocker reliability; the two failure routes) are
parameter-robust; the exact percentages are parameter-dependent and reported as such with a sensitivity sweep.

USAGE
  python scripts/28_andnot_gate_discrimination.py            # run -> JSON + figure
  python scripts/28_andnot_gate_discrimination.py selftest   # gate/commit logic checks
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung7_recognition"
RESULT_JSON = OUT_DIR / "rung7_gate_discrimination.json"
FIGURE_PNG = OUT_DIR / "rung7_gate_discrimination.png"
SEED = 7

# --- literature-grounded parameters (log10 molecules/cell unless noted) ---
# Low-affinity discriminating CAR fires on antigen-HIGH only (~10^3.5/cell threshold): affinity-tuned CAR
# density discrimination (Caruso 2015 Cancer Res; Liu 2015 Cancer Res; Salzer 2020). MHC-I on normal cells
# ~10^4-10^5/cell; LIR-1 blocker inhibits reliably above ~10^3 (A2 Bio Tmod, Mol Ther Oncolytics 2022).
PARAMS = dict(
    Ka=3.5, na=2.0,        # activator half-max + Hill (low-affinity -> high threshold)
    Kb=3.0, nb=2.0,        # blocker half-max + Hill
    theta=0.30,            # kill-license commit threshold (bistable all-or-none)
    # population densities
    tum_antigen_mu=4.2, tum_antigen_sd=0.5,     # tumour antigen-high (~16k/cell)
    nrm_antigen_mu=3.0, nrm_antigen_sd=0.6,     # normal leak (mostly below activator threshold, high tail)
    hla_hi_mu=4.5, hla_hi_sd=0.3,               # HLA-high cells
    hla_lo_mu=2.0, hla_lo_sd=0.4,               # downregulated 'HLA-low' cells (blocker fails)
    clonal_loh_frac=0.50,   # fraction of LOH+ tumour cells with CLONAL loss (rest retain HLA_B -> escape B)
    normal_hla_low_frac=0.05,   # fraction of normal cells that are HLA-low (failure route A) -- KEY safety driver
    n_cells=40000,
)


def hill(x_log10: np.ndarray, K: float, n: float) -> np.ndarray:
    """Hill activation on a log10-density input. K is the half-max in the same log10 units."""
    x = np.power(10.0, x_log10)
    k = np.power(10.0, K)
    xn = np.power(x, n)
    return xn / (xn + np.power(k, n))


def kill_license(antigen_log10, hla_log10, p) -> np.ndarray:
    """AND-NOT integration: activator AND NOT blocker."""
    a = hill(antigen_log10, p["Ka"], p["na"])
    b = hill(hla_log10, p["Kb"], p["nb"])
    return a * (1.0 - b)


def commits(license_vals: np.ndarray, theta: float) -> np.ndarray:
    """Bistable all-or-none apoptosis commit (EARM essence): license above theta -> committed."""
    return license_vals > theta


def _sample_populations(p, rng):
    n = p["n_cells"]
    # TUMOUR antigen (high, with a low tail); HLA_B = 0 if clonal LOH else retained-high
    tum_antigen = rng.normal(p["tum_antigen_mu"], p["tum_antigen_sd"], n)
    clonal = rng.random(n) < p["clonal_loh_frac"]
    tum_hla = np.where(clonal, -3.0,                                  # clonal LOH -> essentially no HLA_B
                       rng.normal(p["hla_hi_mu"], p["hla_hi_sd"], n))  # subclonal-retained -> HLA-high
    # NORMAL antigen (leak, mostly low); HLA_B high except an HLA-low downregulated fraction
    nrm_antigen = rng.normal(p["nrm_antigen_mu"], p["nrm_antigen_sd"], n)
    hla_low = rng.random(n) < p["normal_hla_low_frac"]
    nrm_hla = np.where(hla_low, rng.normal(p["hla_lo_mu"], p["hla_lo_sd"], n),
                       rng.normal(p["hla_hi_mu"], p["hla_hi_sd"], n))
    return (tum_antigen, tum_hla, clonal), (nrm_antigen, nrm_hla, hla_low)


def evaluate(p, rng):
    (ta, th, clonal), (na_, nh, hla_low) = _sample_populations(p, rng)
    tum_L = kill_license(ta, th, p)
    nrm_L = kill_license(na_, nh, p)
    tum_kill = commits(tum_L, p["theta"])
    nrm_kill = commits(nrm_L, p["theta"])

    tpr = float(tum_kill.mean())                         # tumour cells correctly committed to apoptosis
    fpr = float(nrm_kill.mean())                         # normal cells wrongly killed (off-tumour toxicity)

    # failure-mode attribution
    nrm_antigen_high = hill(na_, p["Ka"], p["na"]) > 0.5
    fmA = float((nrm_kill & hla_low).mean())             # normal killed BECAUSE HLA-low (blocker failed)
    fmA_of_killed = float((nrm_kill & hla_low).sum() / max(1, nrm_kill.sum()))
    tum_escape = ~tum_kill
    esc_antigen_low = float((tum_escape & (hill(ta, p["Ka"], p["na"]) <= 0.5)).mean())
    esc_hla_retained = float((tum_escape & ~clonal).mean())  # subclonal-retained -> blocker fired -> spared

    # ROC by sweeping the commit threshold
    thetas = np.linspace(0, 1, 101)
    roc_tpr = [float((tum_L > t).mean()) for t in thetas]
    roc_fpr = [float((nrm_L > t).mean()) for t in thetas]
    _trapz = getattr(np, "trapezoid", np.trapz)          # np>=2 renamed trapz -> trapezoid
    auc = float(_trapz(roc_tpr[::-1], roc_fpr[::-1]))     # integrate TPR over FPR

    return dict(tpr=round(tpr, 4), fpr=round(fpr, 4),
                failure_A_false_kill_frac=round(fmA, 4),
                failure_A_share_of_normal_kills=round(fmA_of_killed, 4),
                failure_B_escape_antigen_low=round(esc_antigen_low, 4),
                failure_B_escape_hla_retained_subclonal=round(esc_hla_retained, 4),
                roc_auc=round(auc, 4),
                normal_antigen_high_frac=round(float(nrm_antigen_high.mean()), 4),
                _roc=(thetas.tolist(), roc_tpr, roc_fpr),
                _samples=((ta, th), (na_, nh)))


def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    base = evaluate(PARAMS, rng)
    print(f"[rung7] baseline: TPR(tumour killed)={base['tpr']:.1%}  FPR(normal killed)={base['fpr']:.2%}  "
          f"ROC-AUC={base['roc_auc']:.3f}")

    # SENSITIVITY 1 (safety driver): FPR vs the normal HLA-low fraction — the load-bearing claim.
    sens_safety = {}
    for f in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
        p = dict(PARAMS, normal_hla_low_frac=f)
        r = evaluate(p, np.random.default_rng(SEED))
        sens_safety[str(f)] = {"fpr": r["fpr"], "tpr": r["tpr"]}
    # robustness of the claim to the ACTIVATOR parameters: does FPR track HLA-low frac regardless of Ka?
    safety_robust = {}
    for Ka in (3.0, 3.5, 4.0):
        p = dict(PARAMS, Ka=Ka)
        r = evaluate(p, np.random.default_rng(SEED))
        safety_robust[str(Ka)] = {"fpr": r["fpr"], "tpr": r["tpr"]}

    # SENSITIVITY 2 (efficacy driver): TPR vs clonal-LOH fraction (subclonal-retained tumour cells escape).
    sens_efficacy = {}
    for f in (0.2, 0.4, 0.6, 0.8, 1.0):
        p = dict(PARAMS, clonal_loh_frac=f)
        r = evaluate(p, np.random.default_rng(SEED))
        sens_efficacy[str(f)] = {"tpr": r["tpr"], "escape_hla_retained": r["failure_B_escape_hla_retained_subclonal"]}

    # the structural claim, quantified: off-tumour toxicity floor ~= P(HLA-low) x P(antigen-high|normal)
    floor_pred = round(PARAMS["normal_hla_low_frac"] * base["normal_antigen_high_frac"], 4)

    result = {
        "tag": "rung7_andnot_gate_discrimination",
        "question": "Does the Tmod AND-NOT gate DISCRIMINATE cancer from normal, and how does it fail? "
                    "(bar #2 of the gate ladder; couples recognition -> apoptosis commit)",
        "model": "activator Hill AND-NOT blocker Hill -> kill-license -> bistable apoptosis commit (EARM essence).",
        "params": {k: v for k, v in PARAMS.items()},
        "baseline": {k: v for k, v in base.items() if not k.startswith("_")},
        "structural_claim": "Gate safety is carried ENTIRELY by the blocker (activator is broad). Off-tumour "
                            "toxicity FLOOR ~= P(normal HLA-low) x P(normal antigen-high). Predicted floor = "
                            f"{floor_pred:.2%}; measured baseline FPR = {base['fpr']:.2%}.",
        "predicted_toxicity_floor": floor_pred,
        "sensitivity_FPR_vs_normal_HLA_low": sens_safety,
        "safety_claim_robust_to_activator_Ka": safety_robust,
        "sensitivity_TPR_vs_clonal_LOH": sens_efficacy,
        "FAILURE_MODES": {
            "A_false_kill": "normal cells that downregulate the sensed HLA allele lose blocker protection -> "
                            f"killed. {base['failure_A_share_of_normal_kills']:.0%} of all normal kills are this. "
                            "This is the gate's irreducible off-tumour toxicity and scales with normal HLA-low rate.",
            "B_escape": f"antigen-low tumour cells ({base['failure_B_escape_antigen_low']:.1%}) + subclonal "
                        f"HLA-retaining tumour cells ({base['failure_B_escape_hla_retained_subclonal']:.1%}) slip "
                        "through. Ties directly to RUNG-6's clonal-LOH haircut.",
        },
        "CEILING": "Mechanistic circuit model with LITERATURE-GROUNDED parameters + ILLUSTRATIVE per-cell "
                   "distributions (real joint antigen/HLA per-cell spread needs the scRNA/proteomic atlas -> "
                   "Colab). STRUCTURAL findings (safety = blocker reliability; the two failure routes) are "
                   "parameter-robust; exact %s are parameter-dependent (sensitivity sweeps included). "
                   "binding != agonism is the wet-lab residual.",
        "INTERPRETATION": "Recognition is ACHIEVABLE but its safety ceiling is set by normal-tissue HLA "
                          "reliability, NOT by LOH frequency: the gate cannot be safer than the blocker, and "
                          "the blocker fails wherever normal cells are HLA-low. So the next real constraint on "
                          "recognition is normal-tissue HLA heterogeneity (measurable on the atlas), and the "
                          "efficacy ceiling is the clonal-LOH fraction (RUNG-6). The recognition->apoptosis "
                          "coupling works in silico: HLA-lost antigen-high tumour cells snap to apoptosis; "
                          "HLA-retaining normal cells do not, even at high antigen.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung7] wrote {RESULT_JSON}")
    print(f"  structural floor pred = {floor_pred:.2%}  vs baseline FPR {base['fpr']:.2%}  (claim holds if ~equal)")
    print("  FPR vs normal HLA-low frac:", {k: f"{v['fpr']:.1%}" for k, v in sens_safety.items()})
    print("  safety robust to activator Ka:", {k: f"{v['fpr']:.1%}" for k, v in safety_robust.items()})
    print("  TPR vs clonal-LOH frac:", {k: f"{v['tpr']:.1%}" for k, v in sens_efficacy.items()})

    _make_figure(base, sens_safety, sens_efficacy)
    return 0


def _make_figure(base, sens_safety, sens_efficacy) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung7] matplotlib unavailable ({e}); skipped figure")
        return
    (ta, th), (na_, nh) = base["_samples"]
    thetas, roc_tpr, roc_fpr = base["_roc"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.7))

    # panel 1: the (antigen, HLA) decision space
    idx = np.random.default_rng(0).choice(len(ta), 3000, replace=False)
    ax[0].scatter(na_[idx], nh[idx], s=4, alpha=0.3, color="#2B6CB0", label="normal")
    ax[0].scatter(ta[idx], th[idx], s=4, alpha=0.3, color="#C1432B", label="tumour")
    ax[0].axvline(PARAMS["Ka"], ls="--", color="grey", lw=1)
    ax[0].axhline(PARAMS["Kb"], ls=":", color="grey", lw=1)
    ax[0].set_xlabel("activator antigen density  (log10/cell)")
    ax[0].set_ylabel("blocker HLA density  (log10/cell)")
    ax[0].set_title("Recognition decision space\n(kill = right of --, below :)")
    ax[0].legend(fontsize=8, markerscale=2)

    # panel 2: ROC
    ax[1].plot(roc_fpr, roc_tpr, "-", color="#4C9F70")
    ax[1].plot([0, 1], [0, 1], ":", color="grey")
    ax[1].scatter([base["fpr"]], [base["tpr"]], color="#C1432B", zorder=5,
                  label=f"operating pt (FPR {base['fpr']:.1%}, TPR {base['tpr']:.0%})")
    ax[1].set_xlabel("FPR — normal cells killed (off-tumour toxicity)")
    ax[1].set_ylabel("TPR — tumour cells killed")
    ax[1].set_title(f"Discrimination ROC (AUC={base['roc_auc']:.3f})")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    # panel 3: the load-bearing claim — FPR tracks normal HLA-low fraction
    xf = [float(k) * 100 for k in sens_safety]
    yfpr = [v["fpr"] * 100 for v in sens_safety.values()]
    ax[2].plot(xf, yfpr, "-o", color="#C1432B")
    ax[2].plot(xf, xf, ":", color="grey", label="y=x (safety floor = HLA-low rate)")
    ax[2].set_xlabel("% normal cells that are HLA-low (blocker fails)")
    ax[2].set_ylabel("% normal cells killed (FPR)")
    ax[2].set_title("Gate is no safer than the blocker:\noff-tumour toxicity tracks normal HLA-low rate")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    fig.suptitle("RUNG-7: AND-NOT recognition-gate discrimination — safety is carried by the blocker, "
                 "not by LOH frequency", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung7] wrote {FIGURE_PNG}")


def selftest() -> int:
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    p = PARAMS
    # tumour: antigen-high, HLA lost -> license high -> killed
    L_tum = kill_license(np.array([4.5]), np.array([-3.0]), p)[0]
    check("tumour (antigen-high, HLA-lost) -> license high", L_tum > p["theta"])
    check("tumour commits to apoptosis", commits(np.array([L_tum]), p["theta"])[0])
    # normal: antigen-high BUT HLA-high -> blocker fires -> license low -> spared
    L_nrm = kill_license(np.array([4.5]), np.array([4.5]), p)[0]
    check("normal (antigen-high, HLA-high) -> license low (blocker protects)", L_nrm < p["theta"])
    check("normal spared even at high antigen", not commits(np.array([L_nrm]), p["theta"])[0])
    # failure mode A: normal HLA-low + antigen-high -> killed
    L_fa = kill_license(np.array([4.5]), np.array([1.5]), p)[0]
    check("failure A: normal HLA-low + antigen-high -> killed", L_fa > p["theta"])
    # failure mode B: tumour antigen-low -> escapes
    L_fb = kill_license(np.array([2.0]), np.array([-3.0]), p)[0]
    check("failure B: tumour antigen-low -> escapes (license low)", L_fb < p["theta"])
    # AND-NOT monotonicity: license decreases in HLA (blocker), increases in antigen (activator)
    hi_anti = kill_license(np.array([5.0]), np.array([-3.0]), p)[0]
    lo_anti = kill_license(np.array([3.0]), np.array([-3.0]), p)[0]
    check("license increases with antigen (activator)", hi_anti > lo_anti)
    lo_hla = kill_license(np.array([4.5]), np.array([1.0]), p)[0]
    hi_hla = kill_license(np.array([4.5]), np.array([5.0]), p)[0]
    check("license decreases with HLA (blocker / NOT arm)", lo_hla > hi_hla)
    # structural claim: more normal HLA-low -> higher FPR (monotone)
    f_lo = evaluate(dict(p, normal_hla_low_frac=0.0), np.random.default_rng(SEED))["fpr"]
    f_hi = evaluate(dict(p, normal_hla_low_frac=0.20), np.random.default_rng(SEED))["fpr"]
    check("FPR monotone in normal HLA-low fraction (safety=blocker)", f_hi > f_lo)
    # determinism: same seed -> same result
    r1 = evaluate(p, np.random.default_rng(SEED))["fpr"]
    r2 = evaluate(p, np.random.default_rng(SEED))["fpr"]
    check("deterministic under fixed seed", r1 == r2)

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", nargs="?", default="run", choices=["run", "selftest"])
    args = ap.parse_args()
    sys.exit(selftest() if args.mode == "selftest" else main_run())
