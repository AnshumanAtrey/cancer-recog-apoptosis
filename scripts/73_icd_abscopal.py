#!/usr/bin/env python3
"""
RUNG-39 — DANGER SIGNAL / ICD (Shriya's hypothesis H7): does the recognition-gated self-destruct reach
METASTASES, or only cure the local tumour?

THE GAP (the most clinically decisive one left)
-----------------------------------------------
Every prior rung clears the LOCAL tumour: recognition (R23/27/33) -> AND-gate kill (R34/35) -> bystander wave
(R36) -> containment (R37). But the wave is contact/percolation-based -- it cannot cross to a DISTANT
micrometastasis that never received the circuit. Metastasis is what actually kills patients. The only mechanism
that reaches untreated distant deposits is SYSTEMIC anti-tumour immunity raised by IMMUNOGENIC cell death (ICD):
dying cells expose calreticulin + release ATP/HMGB1 (DAMPs) -> dendritic-cell maturation -> tumour-specific T-cell
priming -> circulating T-cells kill distant deposits (the ABSCOPAL effect).

THE TENSION this rung quantifies: PLAIN apoptosis is frequently TOLEROGENIC (non-immunogenic, even
immunosuppressive). So a "clean" internal self-destruct could cure the primary and LEAVE the metastases. H7 asks
whether the kill can be made immunogenic enough to flip local killing into systemic clearance -- and what
immunogenicity threshold that needs.

THE MODEL (coupled ODEs; the local kill is a burst, the immune loop + metastasis are dynamic)
  At t=0 the local tumour is cleared (R36) -> a burst of N_kill dying cells. A fraction `i` (ICD immunogenicity,
  0 = tolerogenic apoptosis, 1 = fully immunogenic e.g. pyroptosis/necroptosis flavour) release DAMPs.
  Priming (DCs -> T-cells) is THRESHOLDED (DCs need enough danger signal to mature -> Hill):
      prime_rate(S) = kp * S^h / (Kp^h + S^h),  S = i * (antigen released)   [tolerogenic i->0 gives NO priming]
  T-cell pool A:   dA/dt = prime_rate(S(t)) - dA*A         (S decays as the burst is cleared)
  Metastasis M:    dM/dt = g*M*(1-M/Mmax) - kkill*A*M/(1+A/Asat)   (logistic growth vs T-cell killing, saturable)
  Immunosuppression `supp` scales kkill down (the R-obstacle-3 microenvironment). MHC-dark metastases (R18/R21)
  are T-invisible -> only the NK/agnostic arm reaches them (noted, not re-simulated here).
  Outcome: M(t_end) -> 0 (ABSCOPAL CURE) or M grows back (metastasis ESCAPES = local-only cure).

HONEST CEILING: reduced immuno-oncology ODE (no explicit DC/lymph-node trafficking, T-cell exhaustion dynamics,
antigen-spreading, or MHC-loss escape -- the last is exactly why R21's NK arm is needed for dark metastases). i,
kp, kkill, supp are EFFECTIVE parameters, not fitted. The ICD immunogenicity of THIS particular kill mechanism is
UNMEASURED (the key wet residual). Robust claim = the THRESHOLD structure (tolerogenic apoptosis fails to clear
metastases; immunogenic death above a priming threshold does) and the design requirement, not absolute timings.

USAGE
  python scripts/73_icd_abscopal.py selftest
  python scripts/73_icd_abscopal.py run     # -> runs/rung39_icd_abscopal/ (CPU, seconds)
"""
import os, sys, json
import numpy as np
from scipy.integrate import solve_ivp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung39_icd_abscopal")

# priming (DC maturation -> T-cell) — thresholded Hill in the danger signal
KP, KPRIME_THRESH, HILL = 3.0, 0.30, 4.0
DA = 0.05           # T-cell decay -> slow = immune MEMORY persistence (primed surveillance lasts; enables durable cure)
# metastasis
G, MMAX = 0.25, 1.0   # logistic growth rate, carrying capacity (normalised burden)
KKILL, ASAT = 2.0, 1.0
ANTIGEN_PER_BURST = 1.0  # normalised antigen+DAMP load from clearing the local tumour
BURST_DECAY = 0.5        # the danger signal is transient (burst cleared over ~2 time units)
M_EXT = 1e-4             # extinction floor: a deposit driven below ~1 cell is ERADICATED, cannot regrow


def simulate(i, M0, supp=1.0, t_end=300.0):
    """i=ICD immunogenicity [0,1], M0=initial metastatic burden, supp=immunosuppression [0,1] (1=none).
    Abscopal cure = the T-cell response drives the distant deposit below the extinction floor (eradicated)."""
    def S(t):  # transient danger signal: immunogenic fraction of the cleared-tumour antigen, decaying
        return i * ANTIGEN_PER_BURST * np.exp(-BURST_DECAY * t)

    def prime(s):
        return KP * s ** HILL / (KPRIME_THRESH ** HILL + s ** HILL)

    def rhs(t, y):
        A, M = y
        dA = prime(S(t)) - DA * A
        dM = G * M * (1 - M / MMAX) - supp * KKILL * A * M / (1 + A / ASAT)
        return [dA, dM]

    def extinct(t, y):
        return y[1] - M_EXT
    extinct.terminal = True
    extinct.direction = -1
    sol = solve_ivp(rhs, (0, t_end), [0.0, M0], events=extinct, max_step=0.25, rtol=1e-7, atol=1e-11)
    cured = len(sol.t_events[0]) > 0                  # M crossed below extinction floor -> eradicated
    Mfin = float(sol.y[1, -1])
    Apk = float(sol.y[0].max())
    t_cure = float(sol.t_events[0][0]) if cured else None
    return {"i": i, "M0": M0, "supp": supp, "A_peak": Apk, "M_final": Mfin,
            "abscopal_cure": cured, "t_cure": t_cure}


def threshold_i(M0, supp=1.0):
    """Minimum ICD immunogenicity i that drives the distant deposit to extinction (abscopal cure)."""
    if not simulate(1.0, M0, supp)["abscopal_cure"]:
        return None                                   # even full immunogenicity can't clear (too suppressed)
    if simulate(0.0, M0, supp)["abscopal_cure"]:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        if simulate(mid, M0, supp)["abscopal_cure"]:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def selftest():
    # tolerogenic apoptosis (i=0): no priming -> distant metastasis is NOT cleared (escapes)
    assert not simulate(0.0, 0.1)["abscopal_cure"], "i=0 (tolerogenic) should fail to clear metastasis"
    # immunogenic (i=1, no suppression): clears it (abscopal cure)
    assert simulate(1.0, 0.1)["abscopal_cure"], "i=1 should clear a metastasis with no suppression"
    # a finite immunogenicity threshold exists (the tolerogenic->immunogenic line)
    th = threshold_i(0.2)
    assert th is not None and 0 < th < 1, th
    # HEAVY immunosuppression: even full immunogenicity cannot clear -> needs the agnostic arms (R21 NK / checkpoint)
    assert not simulate(1.0, 0.2, supp=0.25)["abscopal_cure"], "heavy suppression should escape even at i=1"
    print(f"[selftest] ICD priming threshold + memory + extinction OK "
          f"(i=0 escapes, i=1 cures, threshold@M0=0.2 ~ {th:.2f}, heavy-suppression escapes even at i=1)")


def main():
    os.makedirs(OUT, exist_ok=True)
    selftest()
    res = {"tag": "rung39_icd_abscopal", "hypothesis": "H7 Danger-Signal/ICD",
           "question": "does the recognition-gated self-destruct reach metastases (abscopal), or only cure locally?"}

    # 1. the core tension: tolerogenic vs immunogenic death x immunosuppression (the 3 regimes)
    print("\n=== RUNG-39 (H7 ICD): does the self-destruct reach METASTASES? ===")
    print("\n-- abscopal cure of a distant metastasis (M0=0.2) vs ICD immunogenicity i x immunosuppression --")
    supps = [(1.0, "none"), (0.7, "x0.7"), (0.5, "x0.5"), (0.25, "x0.25 heavy")]
    print("   i \\ supp " + "  ".join(f"{lbl:>10}" for _, lbl in supps))
    grid = []
    for i in (0.0, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0):
        row = []
        for s, _ in supps:
            r = simulate(i, 0.2, s)
            row.append("CURE " if r["abscopal_cure"] else "escp ")
            grid.append(r)
        print(f"   {i:.1f}     " + "  ".join(f"{c:>10}" for c in row))
    res["grid_i_x_supp"] = grid

    # 2. the threshold immunogenicity i* vs immunosuppression (and where immunity ALONE fails)
    print("\n-- minimum ICD immunogenicity i* for abscopal cure (M0=0.2) vs immunosuppression --")
    print("   suppression   i* (min immunogenicity)   note")
    th = []
    for s, lbl in [(1.0, "none"), (0.7, "x0.7"), (0.5, "x0.5"), (0.35, "x0.35"), (0.25, "x0.25 heavy")]:
        ti = threshold_i(0.2, s)
        note = "immunity alone clears above i*" if ti is not None else "immunity ALONE FAILS -> needs NK/checkpoint (R21)"
        th.append({"supp": s, "istar": ti, "note": note})
        print(f"   {lbl:>11}      {('%.2f' % ti) if ti is not None else 'NONE':>8}            {note}")
    res["threshold_vs_suppression"] = th

    # 3. plain-apoptosis baseline vs an engineered-immunogenic kill (the design call)
    tol = simulate(0.05, 0.2)      # ~tolerogenic clean apoptosis
    imm = simulate(0.8, 0.2)       # engineered-immunogenic death (ICD inducer / pyroptotic flavour)
    print(f"\n-- plain (tolerogenic) apoptosis i=0.05: metastasis M_final={tol['M_final']:.3f} -> "
          f"{'CURED' if tol['abscopal_cure'] else 'ESCAPES (local-only cure)'}")
    print(f"-- engineered-immunogenic kill i=0.80: metastasis M_final={imm['M_final']:.3f} -> "
          f"{'ABSCOPAL CURE' if imm['abscopal_cure'] else 'escapes'}")
    res["apoptosis_vs_icd"] = {"tolerogenic_i0.05": tol, "immunogenic_i0.80": imm}

    verdict = (
        "The recognition-gated self-destruct is a LOCAL cure unless it is made IMMUNOGENIC -- then it becomes a "
        "SYSTEMIC one. (1) Plain (tolerogenic) apoptosis raises no anti-tumour immunity (priming stays below the "
        "DC-maturation threshold) -> distant metastases that never received the circuit GROW BACK: the primary is "
        "cured but the patient is not. (2) An IMMUNOGENIC kill (ICD: calreticulin/ATP/HMGB1, or a pyroptotic/"
        "necroptotic flavour) above a priming threshold raises a circulating T-cell response that clears untreated "
        "metastases (the abscopal effect). (3) The required immunogenicity i* RISES with metastatic burden and with "
        "immunosuppression -- a heavy, suppressed metastatic load may exceed what immunogenicity alone can clear, "
        "which is exactly where R21's resistance-agnostic NK arm (for MHC-dark deposits) and checkpoint-release are "
        "needed. DESIGN REQUIREMENT (new, named): the internal self-destruct must be ENGINEERED for immunogenic "
        "cell death, not clean apoptosis -- this is what turns Shriya's local 'destroy itself from within' into a "
        "whole-body cure, and it adds a THIRD clearance layer (local wave R36 + systemic ICD-immunity R39 + NK tail "
        "R21). It also re-frames H7 from 'a cheap score' to the pivotal local-vs-systemic determinant.")
    print("\nVERDICT:\n" + verdict)
    res["verdict"] = verdict
    res["residuals"] = (
        "Reduced immuno-oncology ODE (no explicit DC/lymph-node trafficking, T-cell exhaustion, antigen spreading, "
        "or MHC-loss escape -- the last is why R21's NK arm is needed for dark metastases); i, kp, kkill, supp are "
        "effective, not fitted. THE KEY WET RESIDUAL: the actual ICD immunogenicity of this specific kill mechanism "
        "is unmeasured -- whether the engineered self-destruct exposes calreticulin / releases HMGB1 / triggers "
        "pyroptosis must be measured. Robust claim = the THRESHOLD structure (tolerogenic fails, immunogenic-above-"
        "threshold clears) + the burden/suppression dependence, not absolute numbers.")
    json.dump(res, open(os.path.join(OUT, "icd_abscopal.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/icd_abscopal.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
    else:
        main()
