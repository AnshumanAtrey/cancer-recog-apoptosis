#!/usr/bin/env python3
"""
RUNG-35 — KILL-CIRCUIT KINETICS + LEAK-CORRELATION SENSITIVITY.

This rung stress-tests the TWO things RUNG-34 (kill-coupling) left as assumptions:
  (1) KINETICS  — RUNG-34 used a STATIC Hill p_death(effector). Real apoptosis is an irreversible BISTABLE
                  snap (MOMP, eARM RUNG-1): a RACE between caspase accumulation and signal clearance. Does the
                  effector cross the commitment threshold before a transient leak clears?
  (2) LEAK MODEL + CORRELATION — RUNG-34's headline 4e-11 (N=2, 5% leak) assumed leak is a *uniform sub-threshold
                  LEVEL* L that the steep Hill crushes (L^N -> through Hill -> ~0). That is OPTIMISTIC. Off-target
                  CRISPR cuts / toehold mis-triggers are ALL-OR-NONE per cell: leak = a FRACTION of normal cells
                  that FULLY mis-fire. In those cells the effector is FULL -> Hill does NOT save you -> false-death
                  ~ P(all N fire) ~ L^N (independent) ... up to ~L (correlated). 8 orders worse than RUNG-34.

So this script does NOT assume RUNG-34's number is right. It re-derives the honest false-death under the
conservative probabilistic-leak model with a correlation knob rho (one-factor Gaussian copula = the "correlated
defaults" / shared-stress math), and then shows what the irreversible KINETICS win back (a temporal-coincidence
filter), yielding the REAL engineering requirements: leak amplitude L, leak DURATION (transience), correlation rho,
and number of inputs N.

Outputs runs/rung35_kill_kinetics/kill_kinetics.json. Pure numpy/scipy; runs on the M2 in seconds.
"""
import os, json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung35_kill_kinetics")

# ----------------------------------------------------------------------------------------------------------------
# PART A — bistable apoptosis ODE (reduced eARM caspase switch; the irreversible MOMP snap of RUNG-1)
# ----------------------------------------------------------------------------------------------------------------
# dC/dt = u + Vmax * C^h/(Kf^h + C^h) - kdeg * C
#   C    = active effector caspase (caspase-3), normalised death output
#   u    = AND-gated recognizer drive (production of active initiator). RUNG-34's `and_effector` feeds in here.
#   Hill term = caspase-3 -> caspase-9 -> caspase-3 positive feedback (ultrasensitive, h=4) + IAP-degradation loop
#   kdeg = caspase/IAP turnover. Time is in units of 1/kdeg (~0.5-2 h biologically; stated, not fitted).
# At u=0 the system is BISTABLE (survival C~0 AND death C~0.87 both stable) -> once it snaps high it STAYS high
# even when the input is removed = MOMP irreversibility (RUNG-1).
H, KF, VMAX, KDEG = 4.0, 0.5, 1.0, 1.0


def dCdt(C, u):
    return u + VMAX * C**H / (KF**H + C**H) - KDEG * C


def fixed_points(u):
    """Return sorted real fixed points of dC/dt=0 on C in [0, 2] (low-stable, unstable-threshold, high-stable)."""
    grid = np.linspace(0.0, 2.0, 4000)
    f = dCdt(grid, u)
    roots = []
    for i in range(len(grid) - 1):
        if f[i] == 0.0:
            roots.append(grid[i])
        elif f[i] * f[i + 1] < 0:
            roots.append(brentq(lambda c: dCdt(c, u), grid[i], grid[i + 1]))
    # dedupe
    out = []
    for r in roots:
        if not any(abs(r - o) < 1e-6 for o in out):
            out.append(r)
    return sorted(out)


def u_saddle_node():
    """Smallest input u at which the low (survival) state DISAPPEARS -> any u>=u_crit forces commitment.
    This is the DYNAMIC input threshold, to be compared with RUNG-34's static Hill K=0.30."""
    def n_low(u):
        fps = fixed_points(u)
        # count stable fixed points below the unstable threshold (i.e. a surviving low state exists)
        return len(fps)
    # scan u upward until only one fixed point (high) remains
    us = np.linspace(0.0, 0.6, 1201)
    last_multi = 0.0
    for u in us:
        if len(fixed_points(u)) >= 3:
            last_multi = u
        elif len(fixed_points(u)) == 1 and last_multi > 0:
            return 0.5 * (last_multi + u)
    return float("nan")


def commits(u_func, t_max=200.0, C0=0.0, commit_level=0.5):
    """Integrate the caspase ODE under a (possibly time-varying) input u_func(t); return (committed, t_commit, Cfin).
    Committed = effector caspase crosses into the death basin (C>commit_level at t_max, irreversibly)."""
    def rhs(t, y):
        return [dCdt(y[0], u_func(t))]

    def ev(t, y):
        return y[0] - commit_level
    ev.terminal = False
    ev.direction = 1
    sol = solve_ivp(rhs, (0.0, t_max), [C0], events=ev, max_step=0.5, rtol=1e-7, atol=1e-9)
    Cfin = float(sol.y[0, -1])
    committed = Cfin > commit_level
    t_commit = float(sol.t_events[0][0]) if (committed and len(sol.t_events[0]) > 0) else None
    return committed, t_commit, Cfin


def min_pulse_duration(u_amp, t_max=400.0):
    """For a supra-threshold pulse of amplitude u_amp, the MINIMUM duration it must persist to commit the cell.
    This is the kinetic temporal-coincidence filter: a brief leak (even above u_crit) does NOT kill."""
    if not commits(lambda t: u_amp, t_max=t_max)[0]:
        return float("inf")  # even sustained doesn't commit

    def f(T):
        return 1.0 if commits(lambda t: u_amp if t < T else 0.0, t_max=t_max)[0] else -1.0
    lo, hi = 0.01, t_max
    if f(hi) < 0:
        return float("inf")
    if f(lo) > 0:
        return lo
    # bisection on the step-duration
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------------------------------------------------
# PART B — leak as PROBABILISTIC per-cell firing + correlation (one-factor Gaussian copula = Vasicek)
# ----------------------------------------------------------------------------------------------------------------
def joint_fire(L, N, rho):
    """P(all N recognizers mis-fire in the SAME normal cell), each with marginal leak L, pairwise correlation rho
    via a shared latent stress factor Z:  X_i = sqrt(rho) Z + sqrt(1-rho) eps_i,  fire if X_i > t,  t=Phi^-1(1-L).
    rho=0 -> L^N (RUNG-34's independence). rho=1 -> L (all fire together = catastrophic)."""
    if N == 1:
        return L
    if rho <= 1e-12:
        return L**N
    if rho >= 1 - 1e-12:
        return L
    t = norm.ppf(1.0 - L)
    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)

    def integrand(z):
        p = norm.cdf((a * z - t) / b)   # P(fire | shared stress z)
        return norm.pdf(z) * p**N
    val, _ = quad(integrand, -8.0, 8.0, limit=200)
    return float(val)


def required_N(L, rho, target):
    """Smallest N such that P(all N fire) <= target false-death (conservative: every joint-fire kills)."""
    for N in range(1, 13):
        if joint_fire(L, N, rho) <= target:
            return N
    return None


# ----------------------------------------------------------------------------------------------------------------
def selftest():
    # bistability at u=0: three fixed points (survival, threshold, death)
    fp0 = fixed_points(0.0)
    assert len(fp0) == 3, f"expected bistable u=0, got {fp0}"
    assert fp0[0] < 0.05 and fp0[-1] > 0.7, fp0
    # irreversibility: from the HIGH state with NO input, the cell STAYS dead
    assert commits(lambda t: 0.0, C0=fp0[-1])[0], "death state must be stable at u=0 (MOMP irreversible)"
    # survival: from C=0 with no input, stays alive
    assert not commits(lambda t: 0.0, C0=0.0)[0]
    # dynamic threshold exists and is finite
    uc = u_saddle_node()
    assert 0.0 < uc < 0.6, uc
    # copula limits + monotonicity
    assert abs(joint_fire(0.05, 2, 0.0) - 0.05**2) < 1e-6
    assert abs(joint_fire(0.05, 2, 1.0) - 0.05) < 1e-6
    seq = [joint_fire(0.05, 2, r) for r in (0.0, 0.3, 0.6, 0.9)]
    assert all(seq[i] <= seq[i + 1] + 1e-9 for i in range(len(seq) - 1)), seq
    print("[selftest] bistable ODE (3 FPs, irreversible) + Gaussian-copula limits/monotonicity OK")
    return fp0, uc


def main():
    os.makedirs(OUT, exist_ok=True)
    fp0, uc = selftest()
    res = {"tag": "rung35_kill_kinetics"}

    # ---- PART A: kinetics --------------------------------------------------------------------------------------
    print("\n=== PART A — bistable kill-circuit kinetics (reduced eARM, irreversible MOMP) ===")
    print(f"u=0 fixed points (survival / threshold / death): "
          f"{fp0[0]:.3f} / {fp0[1]:.3f} / {fp0[2]:.3f}  -> BISTABLE at zero input = irreversible snap")
    print(f"DYNAMIC input threshold u_crit (saddle-node, low state vanishes): {uc:.3f}")
    print(f"(RUNG-34 used a STATIC Hill K=0.30 on the *effector level*; here the threshold is on the *drive* u.)")

    # mutant: AND-gate of strong recognizers (RUNG-34 mut_eff ~0.9^2=0.81) sustained -> commits, fast
    u_mut = 0.81
    c_m, t_m, cf_m = commits(lambda t: u_mut)
    # normal sustained sub-threshold leak (graded interpretation, RUNG-34 norm_eff ~0.05^2=0.0025)
    u_norm_graded = 0.0025
    c_ng, _, cf_ng = commits(lambda t: u_norm_graded)
    print(f"\nMUTANT  (sustained drive u={u_mut}):  committed={c_m}  t_commit~{t_m:.1f} (1/kdeg units)  Cfin={cf_m:.3f}")
    print(f"NORMAL  (sustained graded leak u={u_norm_graded}): committed={c_ng}  Cfin={cf_ng:.4f}  "
          f"-> sub-threshold leak NEVER commits (kinetic filter holds)")

    # the temporal-coincidence filter: a normal cell that FULLY mis-fires (u above u_crit) but only TRANSIENTLY
    tmin_full = min_pulse_duration(u_mut)               # how long a full mis-fire must persist to kill
    tmin_atcrit = min_pulse_duration(uc * 1.5)          # 1.5x threshold pulse
    print(f"\nKINETIC FILTER (minimum sustained duration of a supra-threshold mis-fire to actually commit):")
    print(f"  full mis-fire (u={u_mut}): T_min ~ {tmin_full:.2f} (1/kdeg)")
    print(f"  1.5x-threshold pulse (u={uc*1.5:.2f}): T_min ~ {tmin_atcrit:.2f} (1/kdeg)")
    print(f"  -> a leak event SHORTER than T_min does NOT kill, even at full amplitude. Kinetics add a")
    print(f"     TEMPORAL-coincidence filter on top of the AND logic.")
    res["partA_kinetics"] = {
        "u0_fixed_points": fp0, "u_crit_saddle_node": uc,
        "mutant": {"u": u_mut, "committed": c_m, "t_commit_invkdeg": t_m, "Cfin": cf_m},
        "normal_sustained_graded_leak": {"u": u_norm_graded, "committed": c_ng, "Cfin": cf_ng},
        "Tmin_full_misfire_invkdeg": tmin_full, "Tmin_1p5x_threshold_invkdeg": tmin_atcrit,
        "note": "Bistable irreversible switch. Graded sub-threshold leak never commits; a full mis-fire commits "
                "only if SUSTAINED past T_min -> kinetics filter brief leaks. Time in 1/kdeg (~0.5-2h)."}

    # ---- PART B: the honest leak re-derivation (probabilistic + correlation) -----------------------------------
    print("\n=== PART B — HONEST false-death: probabilistic all-or-none leak + correlation rho ===")
    print("RUNG-34 assumed leak = uniform sub-threshold LEVEL -> Hill crushes it -> 4e-11 (OPTIMISTIC).")
    print("Conservative reality: leak = FRACTION of normal cells that FULLY mis-fire -> false-death ~ P(all N fire).")
    print("\n  N  L     rho=0(indep)  rho=0.3   rho=0.6   rho=1(full)   [RUNG-34 graded-Hill for ref]")
    from math import isclose
    sweepB = []
    for N in (2, 3, 4):
        for L in (0.02, 0.05, 0.10):
            j0 = joint_fire(L, N, 0.0)
            j3 = joint_fire(L, N, 0.3)
            j6 = joint_fire(L, N, 0.6)
            j1 = joint_fire(L, N, 1.0)
            # RUNG-34's graded-through-Hill number, for the apples-to-apples contrast
            graded = float((L**N)**5 / (0.30**5 + (L**N)**5))
            sweepB.append({"N": N, "L": L, "rho0": j0, "rho0.3": j3, "rho0.6": j6, "rho1": j1,
                           "rung34_graded_hill": graded})
            print(f"  {N}  {L:.2f}  {j0:.2e}   {j3:.2e}  {j6:.2e}  {j1:.2e}     {graded:.2e}")
    res["partB_leak_correlation"] = sweepB

    # required N to hit safety targets, vs correlation
    print("\n  Required # of AND inputs N for normal false-death <= target (5% per-recognizer leak):")
    print("  target      rho=0   rho=0.3  rho=0.6  rho=0.9")
    reqN = []
    for target in (1e-4, 1e-6, 1e-9):
        row = {"target": target}
        line = f"  {target:.0e}    "
        for rho in (0.0, 0.3, 0.6, 0.9):
            n = required_N(0.05, rho, target)
            row[f"rho{rho}"] = n
            line += f"  {str(n):>5}"
        reqN.append(row)
        print(line)
    res["partB_required_N"] = reqN

    # ---- PART C: synthesis — kinetics wins back what the probabilistic model lost ------------------------------
    print("\n=== PART C — synthesis: do the kinetics rescue the AND-gate? ===")
    # A mis-fire is a TRANSIENT event of duration T_leak (in 1/kdeg). For the coincident AND of N full mis-fires
    # to actually COMMIT the cell, the overlap must persist >= the commitment time T_commit (~T_min at full drive,
    # since coincident full mis-fires reconstitute the full effector). A mis-fire event is therefore "lethal" only
    # if it lasts long enough for commitment to complete within it:
    #     rescue-reduced marginal leak   L_eff = L * max(0, 1 - T_commit / T_leak)
    # T_leak <= T_commit  -> L_eff = 0  (brief leaks self-clear; FULL kinetic rescue)
    # T_leak -> infinity   -> L_eff -> L (sustained leaks; recovers Part B, NO rescue)
    # then the coincident false-death is the copula joint of the REDUCED marginal:  joint_fire(L_eff, N, rho).
    T_commit = tmin_full
    Lmarg = 0.05
    synth = []
    print(f"  (coincident full mis-fires commit only if overlap >= T_commit~{T_commit:.2f}/kdeg; "
          f"L_eff = L*(1 - T_commit/T_leak))")
    print("  T_leak/T_commit   L_eff    N=2 eff-false-death(rho=0 / 0.3)    N=3 (rho=0 / 0.3)")
    for ratio in (0.5, 1.0, 2.0, 5.0, 20.0):
        T_leak = ratio * T_commit
        L_eff = Lmarg * max(0.0, 1.0 - T_commit / T_leak)
        n2_0, n2_3 = joint_fire(L_eff, 2, 0.0), joint_fire(L_eff, 2, 0.3)
        n3_0, n3_3 = joint_fire(L_eff, 3, 0.0), joint_fire(L_eff, 3, 0.3)
        synth.append({"T_leak_over_T_commit": ratio, "L_eff": L_eff,
                      "N2_rho0": n2_0, "N2_rho0.3": n2_3, "N3_rho0": n3_0, "N3_rho0.3": n3_3})
        print(f"     {ratio:>5.1f}        {L_eff:.4f}    {n2_0:.2e} / {n2_3:.2e}        {n3_0:.2e} / {n3_3:.2e}")
    res["partC_synthesis"] = {"T_commit_invkdeg": T_commit, "L_marginal": Lmarg,
                              "model": "L_eff = L*(1 - T_commit/T_leak); false-death = joint_fire(L_eff,N,rho)",
                              "transience_table": synth}

    verdict = (
        "RUNG-34's headline therapeutic index (~1e10, false-death 4e-11 at N=2) was OPTIMISTIC: it assumed each "
        "recognizer's leak is a uniform sub-threshold LEVEL that the bistable Hill switch filters away. Under the "
        "conservative, mechanistically-honest model where leak = a FRACTION of normal cells that FULLY mis-fire "
        "(all-or-none, as off-target CRISPR cuts / toehold mis-triggers actually are), a 2-input AND-gate at 5% "
        "per-recognizer leak gives false-death ~ 2.5e-3 when leaks are INDEPENDENT (rho=0) and rises toward ~5e-2 "
        "as leaks become CORRELATED (rho->1) -- i.e. the AND-gate alone is NOT safe (need ~1e-4 to 1e-9). "
        "HOWEVER the irreversible bistable KINETICS supply a second, orthogonal filter the static model ignored: "
        "a mis-fire only commits the cell if it is SUSTAINED past T_min and the N mis-fires are TEMPORALLY "
        "coincident. If recognizer leaks are made TRANSIENT (fast-clearing) and INDEPENDENT, this temporal-AND "
        "drives effective false-death back down by several more orders. The corrected design requirement is "
        "therefore explicit and testable: (i) per-recognizer leak amplitude L, (ii) leak DURATION << commitment "
        "time, (iii) leak CORRELATION rho kept low (independent delivery/cell-state triggers), (iv) N>=3 if rho "
        "cannot be driven near zero. The single most important wet measurement remains the recognizer-leak "
        "CORRELATION in stressed normal cells -- it sets which row of Part B you live on.")
    print("\nVERDICT:\n" + verdict)
    res["verdict"] = verdict
    res["residuals"] = (
        "Dimensionless time (1/kdeg ~ 0.5-2h, stated not fitted); reduced 1-variable caspase switch (captures the "
        "irreversible bistable snap of eARM/RUNG-1 but not the full tBid/Bax/Apaf cascade); the transience model in "
        "Part C is a monotone surrogate for the sustained-coincidence probability, not a fitted leak-kinetics. This "
        "quantifies the SPECIFICITY logic + its dynamic robustness, not a built circuit. Wet steps: measure L, the "
        "leak time-constant, and rho across recognizers in primary normal cells.")
    json.dump(res, open(os.path.join(OUT, "kill_kinetics.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/kill_kinetics.json")


if __name__ == "__main__":
    main()
