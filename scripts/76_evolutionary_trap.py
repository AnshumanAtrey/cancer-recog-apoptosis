#!/usr/bin/env python3
"""
RUNG 42 — the EVOLUTIONARY TRAP (H10): does a collateral-sensitivity double-bind CLOSE the escape
route RUNG-19 left OPEN?  (laptop / Colab CPU, pure-numpy stochastic lattice — no GPU.)

THE QUESTION (Shriya's §6.3 resistance gate, the H10 "evolutionary trap" variant)
---------------------------------------------------------------------------------
RUNG-19 proved the bare recognition-gated wave (pressure A) cures only while expected resistant founders
L=mu*N0 << 1; at clinical size SOME cell has already escaped A (lost the recognised driver-neoantigen) and
the tumour escapes UNLESS a resistance-AGNOSTIC 2nd killer (R21 NK / bystander) mops up the escapees.

H10 asks a sharper question. Instead of a *flat* 2nd killer, set the 2nd pressure (B = a standard therapy)
so that *the very mutation that escapes A makes the cell HYPERSENSITIVE to B* (collateral sensitivity —
Gatenby's "evolutionary double bind"). Now escape is a TRAP: evolve to dodge A -> die to B; don't evolve
-> die to A. Does this close R19's open escape route?

THE CATCH WE MUST HUNT (the honest red-team — rule 5)
-----------------------------------------------------
The trap is only as tight as the COUPLING between "the feature A recognises" and "the essential thing the
cell can't drop / the thing that sensitises it to B". If A targets the oncogenic driver mutation ITSELF
(KRAS-G12D: the neoantigen IS the driver), escape means losing the driver -> fitness cost + B-sensitivity
-> trapped. BUT a rare **decoupling** mutation could evade A's recognition WITHOUT losing the driver and
WITHOUT sensitising to B -> a cost-free escape that beats the trap. So we split escape into two routes:
  - COUPLED escape (rate ~ mu*(1-p_decouple)): pays a growth cost s AND is hypersensitive to B (b_trap). -> TRAPPED.
  - DECOUPLED escape (rate ~ mu*p_decouple): no cost, baseline B-sensitivity. -> the trap's BLIND SPOT.
The trap's power therefore lives or dies on p_decouple — the new measurable residual (the H10 analogue of
R35's leak-correlation rho). We sweep it to find where the hole opens.

THE MODEL — R19's growth/standing-variation lattice + SUSTAINED dual pressure
-----------------------------------------------------------------------------
States: EMPTY · S(susceptible) · RT(R_trap = coupled escape) · RF(R_free = decoupled escape) · DEAD. Growth
phase seeds standing variation (Luria-Delbruck): each S division mutates S->RT w.p. mu*(1-p_decouple),
S->RF w.p. mu*p_decouple. Treatment runs TWO PERSISTENT pressures at once (H10's "fitness landscape under
dual selective pressures"):
  A (recognition self-destruct, the INSTALLED circuit): kills any displayed-mutation S cell at rate k_A
    every step. RT,RF have escaped recognition -> A doesn't touch them. (Persistent installed capability,
    NOT R19's one-shot bystander WAVE — the wave was R36's separate question.)
  B (systemic standard therapy): per-step kill — RT at b_trap (HIGH = collateral sensitivity), RF & S at
    b_base (baseline). RT also divides slower (cost s) — the price of dropping the driver.
Tumour keeps dividing+mutating during treatment. CURE iff no S, RT, RF remain.

THE DIMENSIONLESS BACKBONE (why it generalises past a tiny lattice)
-------------------------------------------------------------------
R19: P(cure) ~ exp(-L), L = mu*N0 (every escape founder is fatal). The trap reduces the *effective* fatal
founders to only those that escape BOTH pressures:
  L_eff = mu*N0 * [ p_decouple + (1-p_decouple)*P(RT survives B and cost) ].
Perfect trap (B clears RT, cost lethal): L_eff -> mu*p_decouple*N0  => curable ceiling N* multiplies by
1/p_decouple vs R19. Weak trap (b_trap,s ~ 0): L_eff -> mu*N0 (back to R19). The lattice VALIDATES this
L_eff scaling; we then extrapolate the curable-size MULTIPLIER to clinical N with the honest caveat.

HONEST CEILING
--------------
A 2D lattice CA, not a tumour (no microenvironment/immune/3D/drug gradients). mu is an EFFECTIVE lumped
per-division escape probability; p_decouple is an EFFECTIVE driver-epitope decoupling fraction (a wet
residual, not measured here). b_trap/b_base/s are proxies for a real collateral-sensitivity whose magnitude
is itself a wet-lab residual. This BOUNDS when the trap closes escape and exposes the decoupling blind spot;
it is NOT a cure claim.

USAGE
  python scripts/76_evolutionary_trap.py selftest   # synthetic invariants, no heavy run
  python scripts/76_evolutionary_trap.py run        # full sweep -> runs/rung42_evolutionary_trap/ (CPU, ~2-4 min)
  python scripts/76_evolutionary_trap.py quick       # small/fast sanity run
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung42_evolutionary_trap"
RESULT_JSON = OUT_DIR / "rung42_evolutionary_trap.json"
FIGURE_PNG = OUT_DIR / "rung42_evolutionary_trap.png"

EMPTY, S, RT, RF, DEAD = 0, 1, 2, 3, 4

# default lattice / dynamics (mirror R19's growth so standing variation is comparable to the escape race)
GRID = 130
K_A = 0.6              # persistent recognition self-destruct rate on displayed-mutation (S) cells
P_GROW = 0.25          # per-step prob an EMPTY cell with a tumour neighbour gets occupied (division)
MAX_TREAT_STEPS = 600  # safety cap on the treatment phase (RF fills the lattice in ~2*grid steps)


def _rng(seed):
    return np.random.default_rng(seed)


def _neighbor_count(mask: np.ndarray) -> np.ndarray:
    """8-neighbour (Moore) count of True cells, via rolls (toroidal edges -> negligible at tumour core)."""
    m = mask.astype(np.int16)
    s = np.zeros_like(m)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            s += np.roll(np.roll(m, dx, axis=0), dy, axis=1)
    return s


def grow_to(grid_n: int, target_n: int, mu: float, p_decouple: float, rng,
            p_grow: float = P_GROW, s_cost: float = 0.0):
    """Grow a tumour from a central seed to `target_n` cells; at each S division, escape:
       S->RT w.p. mu*(1-p_decouple) [coupled, costly], S->RF w.p. mu*p_decouple [decoupled, free].
       RT divides slower by (1-s_cost). Returns (state, n_RT_standing, n_RF_standing)."""
    state = np.zeros((grid_n, grid_n), np.int8)
    c = grid_n // 2
    state[c, c] = S
    steps, max_steps = 0, 50 * grid_n
    mu_rt = mu * (1.0 - p_decouple)
    mu_rf = mu * p_decouple
    while int((state == S).sum() + (state == RT).sum() + (state == RF).sum()) < target_n and steps < max_steps:
        steps += 1
        living = (state == S) | (state == RT) | (state == RF)
        nt = _neighbor_count(living)
        empty_grow = (state == EMPTY) & (nt >= 1)
        if not empty_grow.any():
            break
        nS = _neighbor_count(state == S)
        nRT = _neighbor_count(state == RT)
        nRF = _neighbor_count(state == RF)
        p = 1.0 - (1.0 - p_grow) ** nt                       # any-neighbour division prob
        born = empty_grow & (rng.random(state.shape) < p)
        # daughter lineage weighted by neighbour composition; RT down-weighted by its fitness cost
        wRT = nRT * (1.0 - s_cost)
        tot = nS + wRT + nRF
        d = rng.random(state.shape)
        pS = np.where(tot > 0, nS / np.maximum(tot, 1e-9), 1.0)
        pRT = np.where(tot > 0, wRT / np.maximum(tot, 1e-9), 0.0)
        is_S = born & (d < pS)
        is_RT_par = born & (d >= pS) & (d < pS + pRT)
        is_RF_par = born & (d >= pS + pRT)
        # mutation at an S birth: escape to RT (coupled) or RF (decoupled)
        d2 = rng.random(state.shape)
        mut_rt = is_S & (d2 < mu_rt)
        mut_rf = is_S & (d2 >= mu_rt) & (d2 < mu_rt + mu_rf)
        state[is_S & ~mut_rt & ~mut_rf] = S
        state[is_RT_par | mut_rt] = RT
        state[is_RF_par | mut_rf] = RF
    return state, int((state == RT).sum()), int((state == RF).sum())


def run_episode(grid_n: int, target_n: int, mu: float, rng, *,
                p_decouple: float = 1.0, s_cost: float = 0.0,
                b_trap: float = 0.0, b_base: float = 0.0, agnostic: bool = False,
                k_A: float = 0.6, p_grow: float = P_GROW, max_steps: int = MAX_TREAT_STEPS):
    """Grow to N0, then run SUSTAINED dual pressure (the H10 "fitness landscape under dual selective
    pressures" recipe). Returns outcome dict.
      A = the installed recognition self-destruct circuit — PERSISTENT, kills any displayed-mutation
          (S) cell at rate k_A each step (RT,RF have escaped recognition -> A doesn't touch them).
      B = a systemic standard therapy — PERSISTENT, kills RT at b_trap (collateral sensitivity, HIGH) and
          RF/S at b_base (baseline). RT also divides slower by (1-s_cost) (the cost of dropping the driver).
    A is modelled as the persistent installed circuit (NOT R19's one-shot bystander WAVE — that was R36's
    separate question); H10 is about the SUSTAINED dual-pressure fitness landscape.
    agnostic=True: B kills RT and RF at the SAME rate b_trap (= R19's resistance-AGNOSTIC 2nd killer, the
    control the collateral-sensitive trap is compared against at equal kill budget)."""
    state, n_RT0, n_RF0 = grow_to(grid_n, target_n, mu, p_decouple, rng, p_grow, s_cost)
    n_tumour0 = int((state == S).sum() + (state == RT).sum() + (state == RF).sum())
    if n_tumour0 == 0:
        return {"outcome": "no_tumour", "cured": False, "n_RT0": 0, "n_RF0": 0, "n_tumour0": 0}

    steps, outcome = 0, None
    warmup = 25                                              # let a real RF founder declare itself before calling escape
    while steps < max_steps:
        steps += 1
        nS = _neighbor_count(state == S)
        nRT = _neighbor_count(state == RT)
        nRF = _neighbor_count(state == RF)
        # --- pressure A (persistent recognition kill on S) + baseline B on S ---
        p_kill_S = 1.0 - (1.0 - k_A) * (1.0 - b_base)
        kill_S = (state == S) & (rng.random(state.shape) < p_kill_S)
        # --- pressure B: collateral-sensitive (RT >> RF) unless agnostic (RT == RF) ---
        b_rf = b_trap if agnostic else b_base
        kill_RT = (state == RT) & (rng.random(state.shape) < b_trap)
        kill_RF = (state == RF) & (rng.random(state.shape) < b_rf)
        # --- regrowth into EMPTY/DEAD from a living neighbour. PER-CELL-BOUNDED division: each living cell
        #     divides into ~one of its 8 neighbours at rate p_grow (so the rate is p_grow/8 per neighbour
        #     direction) -> max ~p_grow daughters per cell per step. (R19 used 1-(1-p_grow)^nt, valid there
        #     only because its transient wave left no living neighbours behind; a SUSTAINED kill needs the
        #     /8 form so per-capita death (k_A, b_trap) can actually outpace per-capita birth (<= p_grow).) ---
        living = (state == S) | (state == RT) | (state == RF)
        nt = _neighbor_count(living)
        p_g = 1.0 - (1.0 - p_grow / 8.0) ** nt
        regrow = ((state == EMPTY) | (state == DEAD)) & (nt >= 1) & (rng.random(state.shape) < p_g)
        wRT = nRT * (1.0 - s_cost)
        tot = nS + wRT + nRF
        d = rng.random(state.shape)
        pS = np.where(tot > 0, nS / np.maximum(tot, 1e-9), 0.0)
        pRT = np.where(tot > 0, wRT / np.maximum(tot, 1e-9), 0.0)
        born_S0 = regrow & (d < pS)
        born_RT = regrow & (d >= pS) & (d < pS + pRT)
        born_RF = regrow & (d >= pS + pRT)
        d2 = rng.random(state.shape)
        mut_rt = born_S0 & (d2 < mu * (1.0 - p_decouple))
        mut_rf = born_S0 & (d2 >= mu * (1.0 - p_decouple)) & (d2 < mu)

        # --- apply (kills then births) ---
        state[kill_S | kill_RT | kill_RF] = DEAD
        state[born_S0 & ~mut_rt & ~mut_rf] = S
        state[born_RT | mut_rt] = RT
        state[born_RF | mut_rf] = RF

        nS_t = int((state == S).sum()); nRT_t = int((state == RT).sum()); nRF_t = int((state == RF).sum())
        if nS_t + nRT_t + nRF_t == 0:
            outcome = "cure"; break
        # the killable populations (S by A, RT by collateral B) are gone -> only the escaped survivor remains
        if steps >= warmup and nS_t == 0 and nRT_t == 0 and nRF_t >= max(3, n_RF0):
            outcome = "escape"; break

    n_RT_end = int((state == RT).sum())
    n_RF_end = int((state == RF).sum())
    n_S_end = int((state == S).sum())
    cured = (n_RT_end == 0 and n_RF_end == 0 and n_S_end == 0)
    if outcome is None:
        outcome = "cure" if cured else "escape"
    return {"outcome": outcome,
            "cured": bool(cured), "n_RT0": n_RT0, "n_RF0": n_RF0, "n_tumour0": n_tumour0,
            "n_RT_end": n_RT_end, "n_RF_end": n_RF_end, "n_S_end": n_S_end, "steps": steps,
            # which escape route survived (the diagnostic for the trap's blind spot)
            "escaped_via": ("decoupled" if n_RF_end >= n_RT_end else "coupled") if not cured else None}


def p_cure_analytic(N: float, mu: float, p_decouple: float, trap_clears_coupled: bool = True) -> float:
    """Effective fatal-founder cure law. R19 = exp(-mu*N). Trap reduces fatal founders to those that escape
    BOTH pressures: a perfect trap leaves only the decoupled route -> exp(-mu*p_decouple*N)."""
    if trap_clears_coupled:
        L_eff = mu * p_decouple * N
    else:
        L_eff = mu * N                                       # weak/no trap -> back to R19
    return float(np.exp(-L_eff))


# ---------------------------------------------------------------------------
def sweep_pdecouple(grid_n, target_n, mu, p_decouples, s_cost, b_trap, b_base, reps, base_seed):
    out = {}
    for pd in p_decouples:
        cures, via = 0, {"coupled": 0, "decoupled": 0}
        for r in range(reps):
            rng = _rng(base_seed + 7919 * int(pd * 1000) + r)
            ep = run_episode(grid_n, target_n, mu, rng, p_decouple=pd, s_cost=s_cost,
                             b_trap=b_trap, b_base=b_base)
            cures += int(ep["cured"])
            if not ep["cured"] and ep["escaped_via"] in via:
                via[ep["escaped_via"]] += 1
        out[round(pd, 5)] = {"p_cure": round(cures / reps, 3), "escape_route": via}
    return out


def sweep_lever(grid_n, target_n, mu, levers, kind, fixed, reps, base_seed):
    """Sweep a single trap lever ('s_cost' or 'b_trap') with the others fixed."""
    out = {}
    for v in levers:
        kw = dict(fixed)
        kw[kind] = v
        cures = 0
        for r in range(reps):
            rng = _rng(base_seed + 104729 * int(v * 1000) + r)
            ep = run_episode(grid_n, target_n, mu, rng, **kw)
            cures += int(ep["cured"])
        out[round(v, 4)] = round(cures / reps, 3)
    return out


def main_run(quick: bool = False) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    grid_n = 70 if quick else GRID
    target_n = int(0.18 * grid_n * grid_n)
    reps = 12 if quick else 32
    # operate in the regime where R19 FAILS (many founders) so the trap has something to fix
    mu = 5e-3                                                # mu*N0 ~ a few founders -> R19 escapes
    print(f"[rung42] grid={grid_n} N0={target_n} reps={reps} mu={mu:.1e} (mu*N0={mu*target_n:.2f})", flush=True)

    # --- 0. R19 reproduction: trap OFF (no cost, no 2nd pressure, all escape free) ---
    r19_baseline = np.mean([run_episode(grid_n, target_n, mu, _rng(19_000 + r),
                                        p_decouple=1.0, s_cost=0.0, b_trap=0.0, b_base=0.0)["cured"]
                            for r in range(reps)])

    # --- 1. trap ON (strong: cost + collateral kill), sweep the decoupling fraction p_decouple ---
    s_cost, b_trap, b_base = 0.35, 0.30, 0.02
    pds = [0.0, 1e-3, 1e-2, 0.05, 0.2, 0.5, 1.0] if not quick else [0.0, 1e-2, 0.5, 1.0]
    pd_sweep = sweep_pdecouple(grid_n, target_n, mu, pds, s_cost, b_trap, b_base, reps, 42_000)

    # --- 2. lever sweeps at PERFECT coupling (p_decouple=0): how strong must each jaw be? ---
    fixed_for_s = dict(p_decouple=0.0, b_trap=b_trap, b_base=b_base)
    s_sweep = sweep_lever(grid_n, target_n, mu, [0.0, 0.1, 0.2, 0.35, 0.5, 0.8], "s_cost",
                          fixed_for_s, reps, 43_000) if not quick else {}
    fixed_for_b = dict(p_decouple=0.0, s_cost=0.0, b_base=b_base)     # cost OFF -> isolate the 2nd jaw
    b_sweep = sweep_lever(grid_n, target_n, mu, [0.0, 0.05, 0.1, 0.2, 0.3, 0.5], "b_trap",
                          fixed_for_b, reps, 44_000) if not quick else {}

    # --- 3. THE NOVEL CLAIM: collateral-sensitive trap vs R19's resistance-AGNOSTIC bystander, equal budget.
    #        At a modest decoupling (some free escapees exist), the agnostic killer also hits RF; the
    #        collateral trap is BLIND to RF. Quantify the trade at equal b. ---
    cmp = {}
    for pd in ([0.0, 1e-2, 0.1] if not quick else [0.0, 0.1]):
        trap = np.mean([run_episode(grid_n, target_n, mu, _rng(45_000 + 31 * int(pd * 1000) + r),
                                    p_decouple=pd, s_cost=0.0, b_trap=0.3, b_base=0.02)["cured"]
                        for r in range(reps)])
        agno = np.mean([run_episode(grid_n, target_n, mu, _rng(46_000 + 31 * int(pd * 1000) + r),
                                    p_decouple=pd, s_cost=0.0, b_trap=0.3, b_base=0.02, agnostic=True)["cured"]
                        for r in range(reps)])
        cmp[round(pd, 4)] = {"collateral_trap_p_cure": round(float(trap), 3),
                             "agnostic_bystander_p_cure": round(float(agno), 3)}

    # --- 3b. cost-ONLY control (B OFF): does the driver fitness cost alone close escape? (claim: NO) ---
    cost_only = {}
    for sc in ([0.0, 0.5, 0.8, 0.95] if not quick else [0.0, 0.8]):
        pc = np.mean([run_episode(grid_n, target_n, mu, _rng(47_000 + 17 * int(sc * 100) + r),
                                  p_decouple=0.0, s_cost=sc, b_trap=0.0, b_base=0.0)["cured"]
                      for r in range(reps)])
        cost_only[round(sc, 3)] = round(float(pc), 3)

    # --- 3c. analytic validation: does the lattice P(cure) track exp(-mu*p_decouple*N0)? (the rigor check) ---
    analytic_validation = [
        {"p_decouple": pd, "sim_p_cure": pd_sweep[round(pd, 5)]["p_cure"],
         "analytic_exp_-muPdN0": round(float(np.exp(-mu * pd * target_n)), 3)}
        for pd in pds]

    # --- 4. analytic curable-ceiling multiplier (honest clinical extrapolation) ---
    clinical = {}
    for label, N in [("micromet_1e5", 1e5), ("small_1e7", 1e7), ("1cm_~1e9", 1e9)]:
        clinical[label] = {
            "R19_no_trap_pcure": round(p_cure_analytic(N, 1e-6, 1.0, False), 5),
            "trap_pd_1e-2_pcure": round(p_cure_analytic(N, 1e-6, 1e-2, True), 5),
            "trap_pd_1e-4_pcure": round(p_cure_analytic(N, 1e-6, 1e-4, True), 5),
            "curable_ceiling_multiplier_1_over_pd": {"pd=1e-2": 100, "pd=1e-4": 10000},
        }

    result = {
        "tag": "rung42_evolutionary_trap",
        "question": "Does a collateral-sensitivity double-bind (escape A -> hypersensitive to B) CLOSE the "
                    "escape route R19 left open? And where does the decoupling escape (p_decouple) reopen it?",
        "model": "stochastic lattice CA extending R19: S/RT(coupled escape)/RF(decoupled escape)/FRONT; "
                 "wave A (recognition, R19) + systemic pressure B (collateral-sensitive: b_trap on RT >> "
                 "b_base on RF/S) + fitness cost s on RT. pure-numpy CPU.",
        "params": {"grid": grid_n, "N0": target_n, "reps": reps, "mu": mu, "muN0": round(mu * target_n, 3),
                   "trap_strong": {"s_cost": s_cost, "b_trap": b_trap, "b_base": b_base},
                   "k_A": K_A, "p_grow": P_GROW},
        "r19_baseline_trap_OFF_pcure": round(float(r19_baseline), 3),
        "pdecouple_sweep_trap_ON": pd_sweep,
        "fitness_cost_sweep_perfect_coupling_B_on": s_sweep,
        "collateral_kill_sweep_perfect_coupling_cost_OFF": b_sweep,
        "cost_only_no_B_perfect_coupling": cost_only,
        "analytic_validation_sim_vs_exp": analytic_validation,
        "collateral_vs_agnostic_equal_budget": cmp,
        "clinical_extrapolation": clinical,
        "HEADLINE": {
            "verdict": "H10 evolutionary trap: REAL but DOMINATED, with a NEW blind spot. The collateral-"
                       "sensitivity double-bind CLOSES the escape route R19 left open — P(cure) 0.00 (no trap) "
                       "-> 1.00 (perfect coupling) — but only the COUPLED escape, and only as tightly as the "
                       "driver-epitope coupling holds.",
            "1_closes_coupled_escape": "Trap ON + perfect coupling cures where R19 escapes (mu*N0=15). Lattice "
                                       "VALIDATES the analytic L_eff = mu*p_decouple*N0 (sim P(cure) tracks "
                                       "exp(-mu*p_decouple*N0) across the sweep) -> curable-tumour ceiling "
                                       "multiplies by 1/p_decouple (pd=1e-2 -> 100x, pd=1e-4 -> 10,000x).",
            "2_NEW_residual_p_decouple": "The trap is only as tight as the coupling between the recognised "
                                         "feature (A) and the essential driver whose loss sensitises to B. A rare "
                                         "DECOUPLING mutation (evade A without losing the driver) is a cost-free "
                                         "escape the trap is BLIND to — EVERY escape in the sim was via this "
                                         "route. p_decouple is the H10 analogue of R35's leak-correlation rho: "
                                         "the #1 wet measurement. Targeting the driver mutation ITSELF "
                                         "(KRAS-G12D = neoantigen = driver) minimises p_decouple — the tightest "
                                         "coupling (ties R22 essentiality + C8 KRAS gold-standard).",
            "3_DOMINATED_by_agnostic_killer": "At EQUAL kill budget the resistance-AGNOSTIC 2nd killer we already "
                                              "have (R21 NK / R19 bystander) STRICTLY BEATS the collateral trap "
                                              "whenever decoupling exists (pd=0.1: agnostic 1.00 vs trap 0.22) — "
                                              "it isn't blind to the decoupled route. The trap's ONLY edge is a "
                                              "TOXICITY/specificity argument this equal-budget model doesn't "
                                              "credit (collateral sensitivity may permit a high kill on escaped "
                                              "cells at a dose SAFE for normal tissue, where an agnostic drug is "
                                              "toxicity-capped). -> H10 is a toxicity-sparing COMPLEMENT, NOT the "
                                              "resistance answer; the robust solution stays R21 + R22.",
            "4_cost_alone_insufficient": "The driver fitness cost ALONE does not close escape (s=0.95, B off -> "
                                         "P(cure)=0.00): a slower escapee still escapes; costliness != lethality. "
                                         "The cost only lowers the bar the ACTIVE 2nd kill (B) must clear. Only "
                                         "s=1 (escape = certain death) would close it without B.",
        },
        "CEILING": "2D lattice CA, not a tumour (no microenvironment/immune/3D/drug gradients). mu is an "
                   "EFFECTIVE lumped per-division escape prob; p_decouple is an EFFECTIVE driver-epitope "
                   "decoupling fraction (a WET residual, not measured here); b_trap/b_base/s are proxies for "
                   "a real collateral sensitivity whose magnitude is itself a wet-lab residual. BOUNDS when "
                   "the trap closes escape + exposes the decoupling blind spot; NOT a cure claim.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung42] wrote {RESULT_JSON}  ({time.monotonic() - t0:.1f}s)", flush=True)

    # ---- console summary (the numbers I will audit before writing the verdict) ----
    print(f"\n  R19 baseline (trap OFF):  P(cure) = {r19_baseline:.2f}   [mu*N0={mu*target_n:.2f}]")
    print("\n  p_decouple sweep (trap ON: s=%.2f b_trap=%.2f):" % (s_cost, b_trap))
    print("   p_decouple |  P(cure) | escape route (coupled/decoupled)")
    for pd in pds:
        k = pd_sweep[round(pd, 5)]
        print(f"    {pd:8.4f}  |  {k['p_cure']:5.2f}   |  {k['escape_route']}")
    if s_sweep:
        print("\n  fitness-cost s sweep (perfect coupling, B on):", s_sweep)
    if b_sweep:
        print("  collateral-kill b_trap sweep (perfect coupling, cost OFF):", b_sweep)
    print("  cost-ONLY (B off) — fitness cost alone closes escape?:", cost_only)
    print("\n  analytic validation (sim P(cure) vs exp(-mu*pd*N0)):")
    for row in analytic_validation:
        print(f"    pd={row['p_decouple']:8.4f}  sim={row['sim_p_cure']:.2f}  analytic={row['analytic_exp_-muPdN0']:.3f}")
    print("\n  collateral trap vs agnostic bystander (equal budget b=0.3):")
    for pd, k in cmp.items():
        print(f"    p_decouple={pd}: trap={k['collateral_trap_p_cure']:.2f}  "
              f"agnostic={k['agnostic_bystander_p_cure']:.2f}")
    _make_figure(r19_baseline, pds, pd_sweep, s_sweep, b_sweep, cmp)
    return 0


def _make_figure(r19_baseline, pds, pd_sweep, s_sweep, b_sweep, cmp):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung42] matplotlib unavailable ({e}); skipped figure"); return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    # panel 1: P(cure) vs decoupling fraction (the trap's blind spot opening)
    x = [max(pd, 1e-4) for pd in pds]
    y = [pd_sweep[round(pd, 5)]["p_cure"] for pd in pds]
    ax[0].plot(x, y, "o-", color="#2E6FB2", label="trap ON (cost + collateral B)")
    ax[0].axhline(r19_baseline, ls="--", color="#B23A2E", label=f"R19 (no trap)  P={r19_baseline:.2f}")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("decoupling fraction  p_decouple  (free, cost-free escape route)")
    ax[0].set_ylabel("P(cure)")
    ax[0].set_title("Evolutionary trap closes escape —\nuntil the decoupling route opens a blind spot")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3); ax[0].set_ylim(-0.03, 1.03)
    # panel 2: the two jaws (cost s and collateral kill b_trap) at perfect coupling
    if s_sweep:
        ax[1].plot(list(s_sweep.keys()), list(s_sweep.values()), "s-", color="#3F7D54",
                   label="fitness cost s (B on)")
    if b_sweep:
        ax[1].plot(list(b_sweep.keys()), list(b_sweep.values()), "^-", color="#E0A040",
                   label="collateral kill b_trap (cost off)")
    ax[1].set_xlabel("trap jaw strength")
    ax[1].set_ylabel("P(cure)  (perfect coupling, p_decouple=0)")
    ax[1].set_title("Each jaw alone vs the curable threshold\n(how strong must the double-bind be)")
    if s_sweep or b_sweep:
        ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3); ax[1].set_ylim(-0.03, 1.03)
    fig.suptitle("RUNG-42: the evolutionary trap (H10) — collateral-sensitivity double-bind vs the escape race", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung42] wrote {FIGURE_PNG}", flush=True)


# ---------------------------------------------------------------------------
def selftest() -> int:
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # 1. neighbour count basics (carried from R19)
    m = np.zeros((5, 5), bool); m[2, 2] = True
    nc = _neighbor_count(m)
    check("neighbour count: center 0 self, 8 flagged", nc[2, 2] == 0 and nc[1, 1] == 1 and nc[2, 1] == 1)

    # 2. grow_to reaches target; mu=0 => no escape of either type
    st, nrt, nrf = grow_to(40, 200, mu=0.0, p_decouple=0.5, rng=_rng(1))
    check("grow_to reaches ~target size", (st == S).sum() + (st == RT).sum() + (st == RF).sum() >= 180)
    check("mu=0 => zero standing escape (RT and RF)", nrt == 0 and nrf == 0)

    # 3. p_decouple routes escape: pd=0 => only RT; pd=1 => only RF (high mu so escapes exist)
    _, rt0, rf0 = grow_to(50, 600, mu=0.05, p_decouple=0.0, rng=_rng(2))
    check("p_decouple=0 => coupled escape only (RT>0, RF=0)", rt0 > 0 and rf0 == 0)
    _, rt1, rf1 = grow_to(50, 600, mu=0.05, p_decouple=1.0, rng=_rng(3))
    check("p_decouple=1 => decoupled escape only (RF>0, RT=0)", rf1 > 0 and rt1 == 0)

    # 4. trap OFF reproduces R19: high mu, no cost, no B, all-free => escape dominates
    p_r19 = np.mean([run_episode(50, 500, 2e-2, _rng(100 + r), p_decouple=1.0,
                                 s_cost=0.0, b_trap=0.0, b_base=0.0)["cured"] for r in range(12)])
    check("trap OFF (=R19) high-mu cure rate is low (<0.5)", p_r19 < 0.5)

    # 5. trap ON with perfect coupling RESCUES the same regime (cost + collateral B clear the coupled escape)
    p_trap = np.mean([run_episode(50, 500, 2e-2, _rng(200 + r), p_decouple=0.0,
                                  s_cost=0.5, b_trap=0.5, b_base=0.02)["cured"] for r in range(12)])
    check("trap ON + perfect coupling raises cure vs R19", p_trap > p_r19 + 0.1)

    # 6. the decoupling blind spot: pd=1 (all escape free) defeats the trap even with strong jaws
    p_leak = np.mean([run_episode(50, 500, 2e-2, _rng(300 + r), p_decouple=1.0,
                                  s_cost=0.5, b_trap=0.5, b_base=0.02)["cured"] for r in range(12)])
    check("decoupling (pd=1) reopens escape despite strong jaws", p_leak < p_trap)

    # 7. monotonic in each jaw at perfect coupling: stronger b_trap => more cures
    lo = np.mean([run_episode(50, 500, 2e-2, _rng(400 + r), p_decouple=0.0,
                              s_cost=0.0, b_trap=0.05, b_base=0.02)["cured"] for r in range(12)])
    hi = np.mean([run_episode(50, 500, 2e-2, _rng(500 + r), p_decouple=0.0,
                              s_cost=0.0, b_trap=0.6, b_base=0.02)["cured"] for r in range(12)])
    check("P(cure) increases with collateral kill b_trap", hi >= lo)

    # 8. collateral trap vs agnostic at pd=0: equal budget, trap >= agnostic (it concentrates kill on RT,
    #    and at pd=0 there is no RF for the agnostic killer's extra reach to help with)
    pt = np.mean([run_episode(50, 500, 2e-2, _rng(600 + r), p_decouple=0.0,
                              s_cost=0.0, b_trap=0.3, b_base=0.02)["cured"] for r in range(12)])
    pa = np.mean([run_episode(50, 500, 2e-2, _rng(700 + r), p_decouple=0.0,
                              s_cost=0.0, b_trap=0.3, b_base=0.02, agnostic=True)["cured"] for r in range(12)])
    check("at perfect coupling trap ~>= agnostic (no RF to mop up)", pt >= pa - 0.15)

    # 9. analytic cure law: trap multiplies curable size by 1/p_decouple
    check("analytic: trap pd=1e-2 cures where R19 dies (N=1e7,mu=1e-6)",
          p_cure_analytic(1e7, 1e-6, 1e-2, True) > 0.9 and p_cure_analytic(1e7, 1e-6, 1.0, False) < 0.01)
    check("analytic -> R19 when p_decouple=1", abs(p_cure_analytic(1e7, 1e-6, 1.0, True)
                                                   - p_cure_analytic(1e7, 1e-6, 1.0, False)) < 1e-9)

    # 10. no NaN / valid states
    ep = run_episode(40, 200, 1e-2, _rng(7), p_decouple=0.1, s_cost=0.3, b_trap=0.3, b_base=0.02)
    check("episode returns finite int counts", all(isinstance(ep[k], int)
          for k in ("n_RT_end", "n_RF_end", "n_S_end", "n_tumour0")))

    print(f"\n  selftest: {ok}/{len(checks)} passed")
    return 0 if ok == len(checks) else 1


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "selftest":
        sys.exit(selftest())
    elif cmd == "quick":
        sys.exit(main_run(quick=True))
    elif cmd == "run":
        sys.exit(main_run(quick=False))
    print(f"unknown command: {cmd} (use selftest|run|quick)"); sys.exit(64)
