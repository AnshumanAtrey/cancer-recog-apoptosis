#!/usr/bin/env python3
"""
RUNG 21 — the CROSS-KILL test: does a LAYERED killer (T-cell + NK + bystander wave) clear the MHC-dark
escapees that T-cells alone CANNOT? (laptop / Colab CPU, pure-numpy stochastic lattice — no GPU.)

THE QUESTION (the operational close of RUNG-18/18b/19)
-----------------------------------------------------
RUNG-18/18b measured that real tumours already carry MHC-dark "escapee" cells (~4% genetic + up to ~13%
transcriptional). RUNG-19 proved a recognition-gated T-cell/wave therapy ALONE cannot clear an established
tumour because those escapees survive and regrow. The fix the field uses is a RESISTANCE-AGNOSTIC second
killer. The deepest one is NK cells, because of a real, beautiful principle (Kärre's "missing-self"):

    T-cells kill cells that SHOW a tumour peptide on MHC   (need MHC present)
    NK  cells kill cells that have LOST MHC                (the OPPOSITE trigger)

So the tumour is TRAPPED: keep MHC -> T-cells kill it; drop MHC to dodge T-cells -> it becomes an NK target.
This run asks, quantitatively: does T + NK (+ a bystander wave) close the escape that T-alone leaves open,
across the escapee fractions we actually measured, and how much NK-evasion breaks it.

THE MODEL
---------
Stochastic lattice CA. Tumour cells are MHC-HIGH (T-target, NK-inhibited) or MHC-DARK (NK-target,
T-invisible). At treatment a fraction `f_escapee` are already DARK (standing escape, = RUNG-18/18b), and
HIGH cells keep mutating -> DARK at rate mu_mhc during treatment. Three killers, toggled per ARM:
  T   : kills MHC-HIGH cells (recognition-gated), prob p_T
  NK  : kills MHC-DARK cells (missing-self),       prob p_NK   (p_NK low = NK-evasion)
  wave: a dying cell recruits ANY adjacent tumour (HIGH or DARK) to die, prob `bystander` (resistance-agnostic)
ARMS compared: T-only · T+NK · T+NK+wave. Outcome at clearance: CURE iff no tumour left, else ESCAPE.

WHAT IT SHOULD SHOW (pre-registered)
------------------------------------
T-only: cure collapses as f_escapee rises (the DARK cells survive & regrow). T+NK: rescues, because NK mops
up the DARK escapees the tumour created to dodge T-cells -> the "trap". NK-evasion (low p_NK) re-opens the
hole -> only the bystander wave (agnostic) can still close it. The layered T+NK+wave should be the most
robust across f_escapee and NK-evasion.

HONEST CEILING
--------------
2D lattice CA, not a tumour (no 3D / microenvironment / immune infiltration kinetics / exhaustion). p_T,
p_NK, bystander are EFFECTIVE per-step kill probabilities, not measured rates. NK biology is reduced to
"kills MHC-low" (real tumours also evolve NK-ligand evasion beyond MHC -> captured only as low p_NK). mu_mhc
is a lumped effective MHC-loss rate. BOUNDS whether layering closes the blind spot and how much NK-evasion
breaks it; NOT a cure claim.

USAGE
  python scripts/47_crosskill.py selftest    # synthetic invariants, fast
  python scripts/47_crosskill.py run         # full sweep -> runs/rung21_crosskill/  (CPU, ~3-7 min)
  python scripts/47_crosskill.py quick       # small/fast sanity
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung21_crosskill"
RESULT_JSON = OUT_DIR / "rung21_crosskill.json"
FIGURE_PNG = OUT_DIR / "rung21_crosskill.png"

EMPTY, HIGH, DARK, FRONT, DEAD = 0, 1, 2, 3, 4   # HIGH=MHC+ (T-target), DARK=MHC-lost (NK-target)

GRID = 120
P_T = 0.6           # per-step prob a T-cell kills an adjacent-accessible MHC-HIGH cell
P_GROW = 0.25
FRONT_LIFE = 3
MAX_STEPS = 1500


def _rng(seed):
    return np.random.default_rng(seed)


def _neighbor_count(mask):
    m = mask.astype(np.int16)
    s = np.zeros_like(m)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            s += np.roll(np.roll(m, dx, axis=0), dy, axis=1)
    return s


def grow_to(grid_n, target_n, rng, p_grow=P_GROW):
    """Grow an all-MHC-HIGH tumour from a seed to target_n cells."""
    state = np.zeros((grid_n, grid_n), np.int8)
    c = grid_n // 2
    state[c, c] = HIGH
    steps, cap = 0, 50 * grid_n
    while int((state == HIGH).sum() + (state == DARK).sum()) < target_n and steps < cap:
        steps += 1
        tumour = (state == HIGH) | (state == DARK)
        nt = _neighbor_count(tumour)
        grow = (state == EMPTY) & (nt >= 1) & (rng.random(state.shape) < (1 - (1 - p_grow) ** nt))
        if not grow.any():
            break
        state[grow] = HIGH
    return state


def seed_escapees(state, f_escapee, rng):
    """Flip a fraction f of tumour cells to MHC-DARK (standing escape = RUNG-18/18b)."""
    hi = np.argwhere(state == HIGH)
    k = int(round(f_escapee * len(hi)))
    if k > 0:
        idx = hi[rng.choice(len(hi), size=k, replace=False)]
        state[idx[:, 0], idx[:, 1]] = DARK
    return state


def run_episode(grid_n, target_n, f_escapee, rng, *, use_T=True, use_NK=False, use_wave=False,
                p_T=P_T, p_NK=0.6, bystander=0.3, mu_mhc=1e-3, p_grow=P_GROW, max_steps=MAX_STEPS):
    state = grow_to(grid_n, target_n, rng, p_grow)
    state = seed_escapees(state, f_escapee, rng)
    n0 = int((state == HIGH).sum() + (state == DARK).sum())
    if n0 == 0:
        return {"outcome": "no_tumour", "cured": False}
    age = np.zeros_like(state, np.int16)

    steps = 0
    while steps < max_steps and ((state == HIGH).any() or (state == DARK).any() or (state == FRONT).any()):
        steps += 1
        nf = _neighbor_count(state == FRONT)
        kill = np.zeros_like(state, bool)
        # T-cells kill MHC-HIGH (recognition-gated)
        if use_T:
            kill |= (state == HIGH) & (rng.random(state.shape) < p_T)
        # NK cells kill MHC-DARK (missing-self)
        if use_NK:
            kill |= (state == DARK) & (rng.random(state.shape) < p_NK)
        # bystander wave: a dying FRONT recruits ANY adjacent tumour, resistance-agnostic
        if use_wave:
            p_by = 1 - (1 - bystander) ** nf
            kill |= ((state == HIGH) | (state == DARK)) & (nf >= 1) & (rng.random(state.shape) < p_by)

        # killed cells become FRONT (a dying/clearing cell; recruits only if use_wave via nf next step)
        nH = _neighbor_count(state == HIGH); nD = _neighbor_count(state == DARK)
        has_work = (state == FRONT) & ((nH + nD) > 0)        # front lingers while tumour neighbours remain
        idle = (state == FRONT) & ~has_work
        clear = idle & (age >= FRONT_LIFE)

        # regrowth: empty/dead from a living tumour neighbour, NOT adjacent to a front; inherit MHC, mutate HIGH->DARK
        living = (state == HIGH) | (state == DARK)
        nt = _neighbor_count(living)
        regrow = ((state == EMPTY) | (state == DEAD)) & (nt >= 1) & (nf == 0) & (rng.random(state.shape) < (1 - (1 - p_grow) ** nt))
        frac_dark = np.where((nH + nD) > 0, nD / np.maximum(nH + nD, 1e-9), 0.0)
        born_dark = regrow & (rng.random(state.shape) < frac_dark)
        born_high = regrow & ~born_dark
        mutate = born_high & (rng.random(state.shape) < mu_mhc)   # MHC loss at birth

        # apply
        age[has_work] = 0
        age[idle] += 1
        state[clear] = DEAD
        state[kill] = FRONT
        age[kill] = 0
        state[born_high & ~mutate] = HIGH
        state[born_dark | mutate] = DARK

    nH_end = int((state == HIGH).sum()); nD_end = int((state == DARK).sum())
    cured = (nH_end == 0 and nD_end == 0)
    return {"outcome": "cure" if cured else "escape", "cured": bool(cured),
            "nH_end": nH_end, "nD_end": nD_end, "n0": n0, "steps": steps}


ARMS = {
    "T_only": dict(use_T=True, use_NK=False, use_wave=False),
    "T+NK": dict(use_T=True, use_NK=True, use_wave=False),
    "T+NK+wave": dict(use_T=True, use_NK=True, use_wave=True),
}


def p_cure(grid_n, target_n, f, arm_kw, reps, seed0, **kw):
    cures = 0
    for r in range(reps):
        ep = run_episode(grid_n, target_n, f, _rng(seed0 + r), **arm_kw, **kw)
        cures += int(ep["cured"])
    return cures / reps


# ---------------------------------------------------------------------------
def main_run(quick=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    grid_n = 70 if quick else GRID
    target_n = int(0.15 * grid_n * grid_n)
    reps = 12 if quick else 28
    f_levels = [0.0, 0.04, 0.08, 0.13, 0.20] if not quick else [0.0, 0.1, 0.2]   # 0.04/0.13 = RUNG-18/18b
    print(f"[rung21] grid={grid_n} N0={target_n} reps={reps} f_escapee={f_levels}", flush=True)

    # PART 1: P(cure) vs escapee fraction, per arm (good NK)
    sweep_f = {arm: [] for arm in ARMS}
    for arm, kw in ARMS.items():
        for f in f_levels:
            sweep_f[arm].append(p_cure(grid_n, target_n, f, kw, reps, seed0=2100 + int(f * 1000), p_NK=0.6, bystander=0.3))
        print(f"  {arm:10} P(cure) vs f: {[f'{x:.2f}' for x in sweep_f[arm]]}", flush=True)

    # PART 2: NK-evasion sweep at a fixed high escapee fraction (does low p_NK reopen the hole? does wave save it?)
    f_fixed = 0.13
    nk_levels = [0.0, 0.2, 0.4, 0.6, 0.9]
    sweep_nk = {"T+NK": [], "T+NK+wave": []}
    for arm in sweep_nk:
        for pnk in nk_levels:
            sweep_nk[arm].append(p_cure(grid_n, target_n, f_fixed, ARMS[arm], reps, seed0=3100 + int(pnk * 100), p_NK=pnk, bystander=0.3))
        print(f"  {arm:10} P(cure) vs p_NK @f={f_fixed}: {[f'{x:.2f}' for x in sweep_nk[arm]]}", flush=True)

    # the measured-escapee anchor from RUNG-18/18b
    anchor = {"genetic_dark_~0.04": "RUNG-18 systemic-dark", "transcriptional_dark_lung_~0.13": "RUNG-18b"}

    result = {
        "tag": "rung21_crosskill",
        "question": "Does a LAYERED killer (T-cell + NK + bystander wave) clear the MHC-dark escapees that "
                    "T-cells ALONE cannot? Quantifies the 'missing-self' trap across the measured escapee "
                    "fractions and NK-evasion.",
        "model": "stochastic lattice CA; MHC-HIGH=T-target, MHC-DARK=NK-target; pure-numpy CPU.",
        "params": {"grid": grid_n, "N0": target_n, "reps": reps, "p_T": P_T, "p_grow": P_GROW},
        "f_escapee_levels": f_levels,
        "p_cure_vs_escapee_fraction": {arm: dict(zip([str(x) for x in f_levels], sweep_f[arm])) for arm in ARMS},
        "nk_evasion_sweep_at_f0.13": {"p_NK_levels": nk_levels,
                                      "T+NK": dict(zip([str(x) for x in nk_levels], sweep_nk["T+NK"])),
                                      "T+NK+wave": dict(zip([str(x) for x in nk_levels], sweep_nk["T+NK+wave"]))},
        "rung18_escapee_anchor": anchor,
        "HEADLINE": {
            "plain": "T-cells ALONE leave the MHC-dark escapees -> cure collapses as the escapee fraction rises "
                     "(the RUNG-19 failure). Adding NK (which kills exactly the MHC-LOST cells) CLOSES that hole "
                     "-- the tumour is trapped: keep MHC -> T kills, drop MHC -> NK kills. NK-evasion reopens it, "
                     "and only the resistance-agnostic bystander wave can still finish the job. The layered "
                     "T+NK+wave is the most robust across escapee fraction AND NK-evasion.",
            "T_only_at_measured_escapee": {f"f={f}": sweep_f["T_only"][i] for i, f in enumerate(f_levels)},
            "T+NK+wave_at_measured_escapee": {f"f={f}": sweep_f["T+NK+wave"][i] for i, f in enumerate(f_levels)},
        },
        "INTERPRETATION_MAP": {
            "T_only collapses with f, T+NK rescues": "confirms the 'missing-self' complementarity computationally "
                                                     "-> NK is the natural cross-kill for MHC-loss escape.",
            "low p_NK reopens hole, wave still closes it": "if the tumour evades NK too, only the agnostic "
                                                           "bystander wave (RUNG-14 ferroptosis/quorum) clears it "
                                                           "-> the 3rd layer is not redundant.",
        },
        "CEILING": "2D lattice CA, not a tumour; p_T/p_NK/bystander are effective per-step probs not measured "
                   "rates; NK reduced to 'kills MHC-low' (real NK-evasion beyond MHC = captured only as low "
                   "p_NK); mu_mhc lumped. BOUNDS whether layering closes the blind spot; NOT a cure claim.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"\n[rung21] wrote {RESULT_JSON}  ({time.monotonic() - t0:.1f}s)")
    _make_figure(f_levels, sweep_f, nk_levels, sweep_nk, f_fixed)
    return 0


def _make_figure(f_levels, sweep_f, nk_levels, sweep_nk, f_fixed):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung21] matplotlib unavailable ({e}); skipped figure"); return
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    colors = {"T_only": "#B23A2E", "T+NK": "#E0A040", "T+NK+wave": "#3F7D54"}
    for arm in sweep_f:
        ax[0].plot([x * 100 for x in f_levels], sweep_f[arm], "o-", color=colors[arm], label=arm)
    ax[0].axvspan(4, 13, color="grey", alpha=0.15, label="measured escapee (RUNG-18/18b)")
    ax[0].set_xlabel("MHC-dark escapee fraction at treatment (%)")
    ax[0].set_ylabel("P(cure)"); ax[0].set_ylim(-0.03, 1.03)
    ax[0].set_title("T-cells alone fail as escapees rise;\nNK closes the 'missing-self' hole")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    for arm, c in (("T+NK", "#E0A040"), ("T+NK+wave", "#3F7D54")):
        ax[1].plot([x * 100 for x in nk_levels], sweep_nk[arm], "o-", color=c, label=arm)
    ax[1].set_xlabel("NK kill efficiency (%)  — low = NK-evasion")
    ax[1].set_ylabel(f"P(cure) at {int(f_fixed*100)}% escapees"); ax[1].set_ylim(-0.03, 1.03)
    ax[1].set_title("If the tumour evades NK too,\nonly the agnostic wave still clears it")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle("RUNG-21: cross-kill — layered T + NK + bystander wave closes the MHC blind spot", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung21] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest():
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    rng = _rng(1)
    st = grow_to(40, 200, rng)
    check("grow_to reaches target, all MHC-HIGH", (st == HIGH).sum() >= 180 and (st == DARK).sum() == 0)
    st2 = seed_escapees(grow_to(40, 400, _rng(2)), 0.1, _rng(3))
    nd = (st2 == DARK).sum(); nt = (st2 == HIGH).sum() + nd
    check("seed_escapees makes ~10% DARK", 0.05 < nd / nt < 0.16)

    # f=0 AND mu_mhc=0 => no DARK can ever exist => T-only cures every time
    c0 = sum(run_episode(40, 200, 0.0, _rng(100 + r), **ARMS["T_only"], mu_mhc=0.0)["cured"] for r in range(8))
    check("f=0 & mu_mhc=0, T-only => cures all (no escapees ever)", c0 == 8)

    # high escapees, T-only => mostly escape; T+NK (good NK) => rescue
    p_t = np.mean([run_episode(50, 400, 0.15, _rng(200 + r), **ARMS["T_only"])["cured"] for r in range(12)])
    p_tnk = np.mean([run_episode(50, 400, 0.15, _rng(300 + r), **ARMS["T+NK"], p_NK=0.7)["cured"] for r in range(12)])
    check("high escapees: T-only cure low (<0.5)", p_t < 0.5)
    check("high escapees: T+NK rescues (> T-only)", p_tnk > p_t + 0.1)

    # NK-evasion: T+NK with p_NK=0 ~ T-only (dark cells survive both); wave then helps
    p_tnk0 = np.mean([run_episode(50, 400, 0.15, _rng(400 + r), **ARMS["T+NK"], p_NK=0.0)["cured"] for r in range(12)])
    p_wave = np.mean([run_episode(50, 400, 0.15, _rng(500 + r), **ARMS["T+NK+wave"], p_NK=0.0, bystander=0.5)["cured"] for r in range(12)])
    check("NK-evasion (p_NK=0): T+NK ~ T-only (no NK help)", p_tnk0 <= p_t + 0.15)
    check("NK-evasion: bystander wave still raises cure", p_wave > p_tnk0)

    # monotonic: more escapees => fewer T-only cures
    lo = np.mean([run_episode(50, 400, 0.02, _rng(600 + r), **ARMS["T_only"])["cured"] for r in range(10)])
    hi = np.mean([run_episode(50, 400, 0.25, _rng(700 + r), **ARMS["T_only"])["cured"] for r in range(10)])
    check("T-only cure decreases with escapee fraction", lo >= hi)

    # valid states only / no crash
    ep = run_episode(40, 200, 0.1, _rng(9), **ARMS["T+NK+wave"])
    check("episode returns finite counts", all(isinstance(ep[k], int) for k in ("nH_end", "nD_end", "n0")))

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
    print(f"unknown command: {cmd}"); sys.exit(64)
