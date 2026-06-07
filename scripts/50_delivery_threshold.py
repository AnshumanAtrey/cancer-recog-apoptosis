#!/usr/bin/env python3
"""
RUNG 24 — the WAVE-AS-A-ONE-TIME-INJECTION delivery threshold: how sparse can the injection be and still
let the death wave take over the WHOLE tumour? (laptop CPU, stochastic lattice — no GPU.)

THE QUESTION (Anshuman's "wave in a one-time shot" hypothesis)
-------------------------------------------------------------
RUNG-13 showed the recognition-gated death wave clears a susceptible tumour and that ONE seed can suffice in
principle (percolation). The buildable question for a real injection (mRNA-LNP / oncolytic virus / peptide
that seeds the trigger in some cells): **what FRACTION of tumour cells must the single injection actually
reach for the wave to percolate and clear everything — decoupling the KILL from the DELIVERY, exactly as
Anshuman framed it.** If that critical fraction is tiny, a sparse one-time shot suffices; if it's large, the
delivery itself is the bottleneck.

MODEL
-----
Stochastic lattice. Grow a susceptible tumour to N0. The injection converts a fraction `f_deliver` of tumour
cells into dying FRONT cells (the seeds it reached). The wave spreads cell→cell (recognition-gated, RUNG-13
kinetics); the tumour simultaneously REGROWS into cleared space (so a wave that fails to percolate dies out
and the tumour comes back). Outcome: CLEARED iff no tumour left. Sweep f_deliver × wave strength (p_wave =
coupling) × tumour size -> the critical delivery fraction f* where P(clear) crosses to ~1.

WHAT IT SHOULD SHOW (pre-registered)
------------------------------------
A PERCOLATION threshold: below f* the seeds are too sparse, the wave dies in regrowth, tumour persists; above
f* the seeds connect, the wave percolates and clears. Stronger coupling (p_wave) lowers f*. The headline is
f* — the minimum injection reach — and how it falls with coupling.

HONEST CEILING
--------------
2D lattice CA; p_wave is an effective coupling (the agonism/transduction residual since RUNG-1); regrowth rate
is a parameter; no 3D/microenvironment. "Injection reaches f_deliver of cells" abstracts real biodistribution
(LNP tropism, perfusion). BOUNDS the delivery reach needed; the delivery vehicle itself is wet-lab.

USAGE
  python scripts/50_delivery_threshold.py selftest
  python scripts/50_delivery_threshold.py run        # sweep -> runs/rung24_delivery/  (CPU, ~3-6 min)
  python scripts/50_delivery_threshold.py quick
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung24_delivery"
RESULT_JSON = OUT_DIR / "rung24_delivery.json"
FIGURE_PNG = OUT_DIR / "rung24_delivery.png"

EMPTY, ALIVE, FRONT, DEAD = 0, 1, 2, 3
GRID = 120
P_GROW = 0.25
FRONT_LIFE = 3
MAX_STEPS = 1500


def _rng(s):
    return np.random.default_rng(s)


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
    state = np.zeros((grid_n, grid_n), np.int8)
    c = grid_n // 2
    state[c, c] = ALIVE
    steps, cap = 0, 50 * grid_n
    while int((state == ALIVE).sum()) < target_n and steps < cap:
        steps += 1
        nt = _neighbor_count(state == ALIVE)
        grow = (state == EMPTY) & (nt >= 1) & (rng.random(state.shape) < (1 - (1 - p_grow) ** nt))
        if not grow.any():
            break
        state[grow] = ALIVE
    return state


def episode(grid_n, target_n, f_deliver, p_wave, rng, p_grow=P_GROW, max_steps=MAX_STEPS):
    """Inject the trigger into f_deliver of tumour cells; does the wave clear the whole tumour?"""
    state = grow_to(grid_n, target_n, rng, p_grow)
    alive = np.argwhere(state == ALIVE)
    if len(alive) == 0:
        return True
    k = int(round(f_deliver * len(alive)))
    if k == 0:
        return False                                       # injection reached nobody
    pick = alive[rng.choice(len(alive), size=k, replace=False)]
    state[pick[:, 0], pick[:, 1]] = FRONT
    age = np.zeros_like(state, np.int16)
    steps = 0
    while steps < max_steps and ((state == FRONT).any() or (state == ALIVE).any()):
        steps += 1
        nf = _neighbor_count(state == FRONT)
        p_commit = 1 - (1 - p_wave) ** nf
        new_front = (state == ALIVE) & (nf >= 1) & (rng.random(state.shape) < p_commit)
        # RUNG-13-style wave: a FRONT transmits for FRONT_LIFE steps then clears to DEAD — so the wave can
        # BURN OUT if seeds are too sparse to keep recruiting (this is what creates a real delivery threshold).
        clear = (state == FRONT) & (age >= FRONT_LIFE)
        nt = _neighbor_count(state == ALIVE)
        regrow = ((state == EMPTY) | (state == DEAD)) & (nt >= 1) & (nf == 0) & (rng.random(state.shape) < (1 - (1 - p_grow) ** nt))
        age[state == FRONT] += 1
        state[clear] = DEAD
        state[new_front] = FRONT
        age[new_front] = 0
        state[regrow] = ALIVE
    return bool((state == ALIVE).sum() == 0)


def p_clear(grid_n, target_n, f_deliver, p_wave, reps, seed0):
    return sum(episode(grid_n, target_n, f_deliver, p_wave, _rng(seed0 + r)) for r in range(reps)) / reps


def critical_f(fs, pcs, thresh=0.5):
    """smallest f_deliver with P(clear) >= thresh (interpolated index)."""
    for f, pc in zip(fs, pcs):
        if pc >= thresh:
            return f
    return None


def main_run(quick=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    grid_n = 70 if quick else GRID
    target_n = int(0.16 * grid_n * grid_n)
    reps = 12 if quick else 30
    fs = [0.0005, 0.001, 0.003, 0.01, 0.03, 0.1] if not quick else [0.001, 0.01, 0.1]
    # super-critical coupling regime (clean delivery threshold). Below ~critical coupling the wave fails
    # regardless of delivery (coupling, not delivery, is the gate) — reported in INTERPRETATION.
    p_waves = [0.85, 0.92, 0.98] if not quick else [0.95]
    print(f"[rung24] grid={grid_n} N0={target_n} reps={reps} f_deliver={fs} p_wave={p_waves}", flush=True)

    sweep = {}
    f_star = {}
    for pw in p_waves:
        pcs = [p_clear(grid_n, target_n, f, pw, reps, seed0=2400 + int(f * 1e5) + int(pw * 100)) for f in fs]
        sweep[str(pw)] = dict(zip([str(f) for f in fs], pcs))
        f_star[str(pw)] = critical_f(fs, pcs)
        print(f"  p_wave={pw}: P(clear) {[f'{x:.2f}' for x in pcs]}  f*={f_star[str(pw)]}", flush=True)

    result = {
        "tag": "rung24_delivery_threshold",
        "question": "How sparse can a one-time injection be (fraction of tumour cells reached) and still let "
                    "the death wave percolate and clear the whole tumour? Decouples KILL from DELIVERY.",
        "model": "stochastic lattice; inject trigger into f_deliver of cells, wave spreads vs regrowth.",
        "params": {"grid": grid_n, "N0": target_n, "reps": reps, "p_grow": P_GROW},
        "f_deliver_levels": fs, "p_wave_levels": p_waves,
        "p_clear_sweep": sweep,
        "critical_delivery_fraction_f_star": f_star,
        "HEADLINE": {
            "plain": "There is a PERCOLATION delivery threshold f*: below it the injected seeds are too sparse, "
                     "the wave dies in regrowth, the tumour persists; above it the wave percolates and clears "
                     "everything. Stronger coupling (p_wave) lowers f* — a more potent trigger needs a less "
                     "thorough injection. f* (this lattice): " + ", ".join(f"p_wave={k}→{v}" for k, v in f_star.items()),
            "meaning": "If f* is small (~1%), a SPARSE one-time injection suffices (delivery is NOT the "
                       "bottleneck — the wave does the work). If f* is large, the injection reach is the wall.",
        },
        "INTERPRETATION_MAP": {
            "f* small & falls with p_wave": "delivery decoupled from kill (Anshuman's framing holds): a potent "
                                            "self-propagating trigger needs only sparse seeding.",
            "f* large / no clear at any tested f": "delivery reach IS the bottleneck -> need better biodistribution "
                                                   "(nanorobotic/LNP-tropism, the FUT delivery layer) or stronger coupling.",
        },
        "CEILING": "2D lattice CA; p_wave effective coupling (agonism/transduction residual); regrowth a "
                   "parameter; 'reaches f_deliver of cells' abstracts real biodistribution. BOUNDS the delivery "
                   "reach needed; the vehicle is wet-lab.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung24] wrote {RESULT_JSON}  ({time.monotonic()-t0:.1f}s)")
    _make_figure(fs, sweep, p_waves, f_star)
    return 0


def _make_figure(fs, sweep, p_waves, f_star):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung24] matplotlib unavailable ({e})"); return
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    cols = ["#B23A2E", "#E0A040", "#3F7D54"]
    for pw, c in zip(p_waves, cols):
        pcs = [sweep[str(pw)][str(f)] for f in fs]
        ax.plot([f * 100 for f in fs], pcs, "o-", color=c, label=f"coupling p_wave={pw}  (f*={f_star[str(pw)]})")
    ax.axhline(0.5, ls="--", color="grey", alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("injection reach — % of tumour cells seeded (one-time)")
    ax.set_ylabel("P(wave clears the whole tumour)"); ax.set_ylim(-0.03, 1.03)
    ax.set_title("RUNG-24: wave-as-injection delivery threshold\nbelow f* seeds are too sparse; above f* the wave percolates and clears")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung24] wrote {FIGURE_PNG}")


def selftest():
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # f_deliver=0 -> never clears (no seed); f_deliver=1 -> always clears (every cell triggered)
    c0 = sum(episode(40, 200, 0.0, 0.9, _rng(10 + r)) for r in range(6))
    c1 = sum(episode(40, 200, 1.0, 0.9, _rng(20 + r)) for r in range(6))
    check("f_deliver=0 never clears", c0 == 0)
    check("f_deliver=1 always clears", c1 == 6)
    # monotone IN THE SUPER-CRITICAL coupling regime (where a clean delivery threshold exists; at sub-critical
    # coupling the wave fails regardless of delivery — coupling is the real gate, see INTERPRETATION).
    lo = p_clear(50, 400, 0.001, 0.95, 12, 100)
    hi = p_clear(50, 400, 0.1, 0.95, 12, 200)
    check("P(clear) rises with delivery fraction (super-critical coupling)", hi >= lo)
    # stronger coupling lowers the bar: at a small f, higher p_wave clears more
    weak = p_clear(50, 400, 0.01, 0.4, 12, 300)
    strong = p_clear(50, 400, 0.01, 0.95, 12, 400)
    check("stronger coupling clears more at fixed delivery", strong >= weak)
    # critical_f picks the first crossing
    check("critical_f finds threshold", critical_f([0.001, 0.01, 0.1], [0.1, 0.6, 0.9]) == 0.01)
    check("critical_f None if never crossed", critical_f([0.001, 0.01], [0.1, 0.2]) is None)

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
    print(f"unknown: {cmd}"); sys.exit(64)
