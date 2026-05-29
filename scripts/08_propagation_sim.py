#!/usr/bin/env python3
"""
Core-hypothesis test — does a SELF-PROPAGATING death wave reverse a tumour?

Tests Shriya Rai's concept ("cancer destroys itself from within via altered cell-cell
communication") in the amplified form Anshuman proposed: cancer spreads cell-to-cell, so
make the DEATH spread the same way — a chain reaction where each dying cancer cell triggers
death in neighbours, propagating through the tumour and reversing its spread.

This is an agent-based tissue model (abstract, parameterised — it tests the LOGIC of
propagation + recognition, not quantitative real-tumour kinetics). It compares death-trigger
mechanisms and asks two questions per mechanism:
  (1) does it CLEAR the tumour (propagation works)?   (2) does it SPARE healthy tissue (specificity)?

MECHANISMS (the "multiple apoptosis triggers" to compare):
  contained_apoptosis  — classic apoptosis: quiet, non-propagating (R=0). CONTROL: should NOT clear.
  diffusive_ligand     — secreted death ligand (TRAIL-like) from the seeded cells only, longer range,
                         NOT self-amplifying, recognition-gated. Tests: does a non-amplified signal suffice?
  ungated_wave         — ferroptosis-like self-sustaining wave, NOT recognition-gated. Tests: clears tumour
                         but does it LEAK into healthy tissue?
  recognition_gated    — self-sustaining wave GATED by the cancer antigen (Shriya's pillar 1). Hypothesis:
                         the only one that BOTH clears the tumour AND spares healthy tissue.

FALSIFIABLE PREDICTION: contained -> fizzles; ungated_wave -> clears but kills healthy; recognition_gated
-> clears tumour with minimal healthy death. If recognition_gated does NOT outperform on the
(cancer_killed - healthy_killed) margin, the "recognition-gated propagation" thesis is wrong.

USAGE:  python scripts/08_propagation_sim.py
REQS :  numpy (matplotlib optional for the figure). CPU, seconds.
"""

from __future__ import annotations

import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "propagation_sim"

GRID = 121                 # odd → centred tumour
SEED = 20260530
STEPS = 120
SEED_FRAC = 0.03           # fraction of cancer cells where death is initially triggered

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", stream=sys.stdout)
log = logging.getLogger("propsim")

# cell-state codes
EMPTY, HEALTHY, CANCER, DEAD = 0, 1, 2, 3


@dataclass
class Mechanism:
    name: str
    radius: int            # neighbourhood the death signal reaches
    self_sustaining: bool  # True: ALL dead cells emit (amplifying wave); False: only the seeded cells emit
    antigen_gated: bool    # True: only kills antigen+ (cancer) cells; False: kills any living cell (leaks)
    threshold: int = 1     # min dead emitting neighbours (within radius) needed to kill a killable cell


def build_tissue(rng) -> tuple[np.ndarray, np.ndarray]:
    """Central solid tumour + scattered micro-metastases, embedded in healthy tissue.
    Returns (state grid, antigen grid). Antigen=1 on cancer (the recognition marker, e.g. Trop2)."""
    state = np.full((GRID, GRID), HEALTHY, np.int8)
    yy, xx = np.mgrid[0:GRID, 0:GRID]
    c = GRID // 2
    # primary tumour: disk radius ~22
    state[(yy - c) ** 2 + (xx - c) ** 2 <= 22 ** 2] = CANCER
    # micro-metastases: small disks scattered (the "spread" we want to reverse)
    for _ in range(6):
        my, mx = rng.integers(12, GRID - 12, size=2)
        if state[my, mx] == HEALTHY:
            state[(yy - my) ** 2 + (xx - mx) ** 2 <= 5 ** 2] = CANCER
    antigen = (state == CANCER).astype(np.float32)   # recognition marker on cancer only
    return state, antigen


def neighbour_offsets(radius: int):
    offs = []
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            if dy == 0 and dx == 0:
                continue
            d = (dy * dy + dx * dx) ** 0.5
            if d <= radius + 1e-9:
                offs.append((dy, dx, d))
    return offs


def run(mech: Mechanism, rng):
    state, antigen = build_tissue(rng)
    cancer0 = int((state == CANCER).sum())
    healthy0 = int((state == HEALTHY).sum())

    # seed initial death in SEED_FRAC of cancer cells
    cancer_idx = np.argwhere(state == CANCER)
    nseed = max(1, int(SEED_FRAC * len(cancer_idx)))
    pick = cancer_idx[rng.choice(len(cancer_idx), size=nseed, replace=False)]
    state[pick[:, 0], pick[:, 1]] = DEAD
    seed_mask = np.zeros((GRID, GRID), bool)
    seed_mask[pick[:, 0], pick[:, 1]] = True
    newly_dead = state == DEAD

    offs = neighbour_offsets(mech.radius) if mech.radius > 0 else []
    history = []

    for t in range(STEPS):
        # emitters: a self-sustaining wave = EVERY dead cell re-emits; otherwise only the seeded cells.
        emitters = (state == DEAD) if mech.self_sustaining else seed_mask
        # count dead emitting neighbours within the mechanism's radius
        dead_neighbours = np.zeros((GRID, GRID), np.int16)
        if mech.radius > 0 and emitters.any():
            em = emitters.astype(np.int16)
            for dy, dx, _ in offs:
                dead_neighbours += np.roll(np.roll(em, dy, 0), dx, 1)
        # who can be killed by this mechanism
        killable = (state == CANCER)
        if not mech.antigen_gated:
            killable = killable | (state == HEALTHY)   # ungated wave can hit healthy too
        dying = killable & (dead_neighbours >= mech.threshold)
        newly_dead = dying.copy()
        state[dying] = DEAD
        cancer_dead = cancer0 - int((state == CANCER).sum())
        healthy_dead = healthy0 - int((state == HEALTHY).sum())
        history.append((cancer_dead / cancer0, healthy_dead / max(healthy0, 1)))
        if not newly_dead.any() and t > 2:
            break

    cancer_killed = history[-1][0]
    healthy_killed = history[-1][1]
    return {
        "mechanism": mech.name, "cancer_killed_frac": round(cancer_killed, 3),
        "healthy_killed_frac": round(healthy_killed, 3),
        "therapeutic_margin": round(cancer_killed - healthy_killed, 3),
        "steps": len(history), "final_state": state, "history": history,
    }


MECHANISMS = [
    # classic apoptosis: quiet/contained — only the triggered cells die (control)
    Mechanism("contained_apoptosis", radius=0, self_sustaining=False, antigen_gated=True),
    # secreted death ligand from the seeded cells only: longer reach but NO amplification → local kill
    Mechanism("diffusive_ligand",   radius=4, self_sustaining=False, antigen_gated=True),
    # ferroptosis-like self-sustaining wave, NOT recognition-gated → clears tumour but leaks to healthy
    Mechanism("ungated_wave",       radius=1, self_sustaining=True,  antigen_gated=False),
    # self-sustaining wave GATED by the cancer antigen (Shriya's pillar 1) → clears tumour, stops at healthy
    Mechanism("recognition_gated",  radius=1, self_sustaining=True,  antigen_gated=True),
]


def main() -> int:
    log.info("propagation sim — does recognition-gated self-propagating death reverse a tumour?")
    log.info("grid=%d seed_frac=%.0f%% steps=%d", GRID, SEED_FRAC * 100, STEPS)
    rng = np.random.default_rng(SEED)
    results = []
    for mech in MECHANISMS:
        r = run(mech, np.random.default_rng(SEED))   # same tissue each time (fair compare)
        results.append(r)
        log.info("[%-20s] cancer_killed=%.1f%%  healthy_killed=%.1f%%  margin=%+.2f  (steps=%d)",
                 r["mechanism"], 100 * r["cancer_killed_frac"], 100 * r["healthy_killed_frac"],
                 r["therapeutic_margin"], r["steps"])

    # ---- verdict / falsification ----
    by = {r["mechanism"]: r for r in results}
    contained = by["contained_apoptosis"]
    gated = by["recognition_gated"]
    ungated = by["ungated_wave"]
    log.info("=" * 64)
    checks = {
        "contained apoptosis fizzles (cancer_killed < 25%)": contained["cancer_killed_frac"] < 0.25,
        "recognition_gated clears tumour (cancer_killed > 80%)": gated["cancer_killed_frac"] > 0.80,
        "recognition_gated spares healthy (healthy_killed < 10%)": gated["healthy_killed_frac"] < 0.10,
        "ungated wave LEAKS to healthy (healthy_killed > recognition_gated)":
            ungated["healthy_killed_frac"] > gated["healthy_killed_frac"] + 0.05,
        "recognition_gated has the best therapeutic margin":
            gated["therapeutic_margin"] == max(r["therapeutic_margin"] for r in results),
    }
    for k, ok in checks.items():
        log.info("  [%s] %s", "✓" if ok else "✗", k)
    supported = all(checks.values())
    log.info("=" * 64)
    if supported:
        log.info("✅ HYPOTHESIS SUPPORTED (in this model): a SELF-PROPAGATING death reverses the tumour,")
        log.info("   and RECOGNITION-GATING is what makes it both lethal-to-cancer AND safe — contained")
        log.info("   apoptosis fizzles; an ungated wave clears but leaks into healthy tissue.")
        log.info("   → This is Shriya's pillar 1 (cell-cell recognition) made load-bearing, not optional.")
    else:
        log.info("⚠️ hypothesis NOT cleanly supported under these params — inspect which check failed and why.")
    log.info("CAVEAT: abstract agent-based model with assumed parameters — tests the LOGIC of")
    log.info("propagation+recognition, NOT quantitative tumour kinetics. Next: PhysiCell with real params.")

    # ---- save ----
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = [{k: v for k, v in r.items() if k not in ("final_state", "history")} for r in results]
    (OUT_DIR / "propagation_results.json").write_text(json.dumps(
        {"results": summary, "checks": checks, "hypothesis_supported": supported}, indent=2))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(["white", "#cfe8cf", "#c0392b", "#222222"])  # empty, healthy, cancer, dead
        fig, axes = plt.subplots(2, len(results), figsize=(4 * len(results), 8))
        for i, r in enumerate(results):
            axes[0, i].imshow(r["final_state"], cmap=cmap, vmin=0, vmax=3)
            axes[0, i].set_title(f"{r['mechanism']}\ncancer {100*r['cancer_killed_frac']:.0f}% / healthy {100*r['healthy_killed_frac']:.0f}%", fontsize=9)
            axes[0, i].axis("off")
            h = np.array(r["history"])
            axes[1, i].plot(h[:, 0], label="cancer killed", color="#c0392b")
            axes[1, i].plot(h[:, 1], label="healthy killed", color="#27ae60")
            axes[1, i].set_ylim(0, 1); axes[1, i].set_xlabel("step"); axes[1, i].legend(fontsize=8)
        fig.suptitle("Self-propagating death vs tumour — recognition-gating is the safe-and-lethal regime", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "propagation_sim.png", dpi=110)
        log.info("figure → %s", (OUT_DIR / "propagation_sim.png").relative_to(PROJECT_ROOT))
    except Exception as e:
        log.warning("figure skipped (%s: %s)", type(e).__name__, e)
    log.info("results → %s", (OUT_DIR / "propagation_results.json").relative_to(PROJECT_ROOT))
    return 0 if supported else 1


if __name__ == "__main__":
    sys.exit(main())
