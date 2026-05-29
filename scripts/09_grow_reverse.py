#!/usr/bin/env python3
"""
Grow-then-reverse — a tumour grows cell-to-cell, then a recognition-gated death wave
reverses it along the SAME routes. The literal form of Anshuman's idea: "the same way
cancer spread, it should kill itself" + Shriya's self-destruction-from-within.

Two phases on one tissue:
  PHASE 1 (GROW):   one cancer cell at the origin proliferates into neighbouring healthy
                    tissue (Eden-like cell-to-cell growth) into a solid tumour (+ occasional
                    distal seeding = a metastatic focus, the honest hard case).
  PHASE 2 (REVERSE): death is triggered at the ORIGIN; a RECOGNITION-GATED self-propagating
                    death wave spreads outward through the tumour (kills cancer only, stops at
                    healthy) — reversing the growth from the inside out.

Outputs (runs/grow_reverse/): a timeline PNG (key frames), the tumour-size curve (grows then
collapses), and an animated GIF of the whole grow→reverse cycle.

CAVEAT: abstract agent-based model with assumed parameters — it visualises the LOGIC of
grow-then-recognition-gated-reverse, not quantitative tumour kinetics.

USAGE:  python scripts/09_grow_reverse.py
REQS :  numpy (+ matplotlib/Pillow optional for figure/GIF). CPU, seconds.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "grow_reverse"

GRID = 161
SEED = 20260530
P_GROW = 0.18            # per-step probability a tumour-adjacent healthy cell becomes cancer
GROW_STEPS = 90
DEATH_STEPS = 90
MET_SEED_PROB = 0.012    # per-step chance of a distal metastatic focus (the honest hard case)

HEALTHY, CANCER, DEAD = 1, 2, 3
NEI8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", stream=sys.stdout)
log = logging.getLogger("growrev")


def neighbour_count(mask: np.ndarray) -> np.ndarray:
    """Count of True 8-neighbours for every cell (toroidal-safe via zero-padded rolls)."""
    acc = np.zeros(mask.shape, np.int16)
    m = mask.astype(np.int16)
    for dy, dx in NEI8:
        acc += np.roll(np.roll(m, dy, 0), dx, 1)
    return acc


def grow(rng):
    """Phase 1: Eden-like cell-to-cell tumour growth from the origin. Returns (state, frames, sizes)."""
    state = np.full((GRID, GRID), HEALTHY, np.int8)
    c = GRID // 2
    state[c, c] = CANCER
    frames, sizes = [], []
    for t in range(GROW_STEPS):
        cancer = state == CANCER
        # candidate healthy cells touching the tumour proliferate with prob P_GROW
        adj = (neighbour_count(cancer) > 0) & (state == HEALTHY)
        grow_here = adj & (rng.random((GRID, GRID)) < P_GROW)
        state[grow_here] = CANCER
        # rare distal metastatic focus
        if rng.random() < MET_SEED_PROB:
            my, mx = rng.integers(15, GRID - 15, size=2)
            if state[my, mx] == HEALTHY:
                state[my, mx] = CANCER
        sizes.append(int((state == CANCER).sum()))
        if t % 6 == 0:
            frames.append(state.copy())
    frames.append(state.copy())
    return state, frames, sizes, c


def reverse(state, origin, rng):
    """Phase 2: recognition-gated self-propagating death wave from the origin. Returns (frames, cancer_sizes, healthy_dead)."""
    c = origin
    # trigger death at the origin (the first cell to 'turn on itself')
    if state[c, c] == CANCER:
        state[c, c] = DEAD
    cancer0 = int((state == CANCER).sum()) + 1  # +1 for the just-killed origin
    healthy0 = int((state == HEALTHY).sum())
    frames, sizes, hist_healthy = [], [], []
    for t in range(DEATH_STEPS):
        dead = state == DEAD
        # recognition-gated: a dead cell kills ADJACENT CANCER cells only (stops at healthy)
        dead_adj = neighbour_count(dead) > 0
        dying = dead_adj & (state == CANCER)
        state[dying] = DEAD
        cancer_left = int((state == CANCER).sum())
        sizes.append(cancer_left)
        hist_healthy.append(healthy0 - int((state == HEALTHY).sum()))
        if t % 6 == 0:
            frames.append(state.copy())
        if not dying.any() and t > 2:
            break
    frames.append(state.copy())
    return frames, sizes, cancer0, healthy0


def render_rgb(state):
    rgb = np.zeros((*state.shape, 3), np.uint8)
    rgb[state == HEALTHY] = (207, 232, 207)
    rgb[state == CANCER] = (192, 57, 43)
    rgb[state == DEAD] = (34, 34, 34)
    return rgb


def main() -> int:
    log.info("grow-then-reverse — tumour grows cell-to-cell, recognition-gated death reverses it")
    rng = np.random.default_rng(SEED)
    state, gframes, gsizes, origin = grow(rng)
    peak = max(gsizes)
    log.info("PHASE 1 GROW: %d steps → tumour peak %d cancer cells (%.1f%% of tissue)",
             GROW_STEPS, peak, 100 * peak / (GRID * GRID))
    rframes, rsizes, cancer0, healthy0 = reverse(state, origin, rng)
    final_cancer = int((state == CANCER).sum())
    healthy_dead = healthy0 - int((state == HEALTHY).sum())
    cleared = 1 - final_cancer / max(cancer0, 1)
    log.info("PHASE 2 REVERSE: death wave from origin → cleared %.1f%% of the tumour, healthy killed %.2f%%",
             100 * cleared, 100 * healthy_dead / max(healthy0, 1))
    if final_cancer > 0:
        log.info("  %d cancer cells survive = disconnected metastatic foc/foci the contact-wave can't reach", final_cancer)

    checks = {
        "tumour grew (peak > 1500 cells)": peak > 1500,
        "reverse cleared the connected tumour (>90%)": cleared > 0.90,
        "healthy spared (<2% killed)": healthy_dead / max(healthy0, 1) < 0.02,
    }
    for k, ok in checks.items():
        log.info("  [%s] %s", "✓" if ok else "✗", k)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "grow_reverse_results.json").write_text(json.dumps({
        "grow_peak_cancer": peak, "reverse_cleared_frac": round(cleared, 3),
        "healthy_killed_frac": round(healthy_dead / max(healthy0, 1), 4),
        "final_cancer_cells": final_cancer, "checks": checks,
        "grow_sizes": gsizes, "reverse_sizes": rsizes,
    }, indent=2))

    # timeline figure + population curve
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        allframes = gframes + rframes
        pick = [gframes[0], gframes[len(gframes)//2], gframes[-1],
                rframes[len(rframes)//3], rframes[2*len(rframes)//3], rframes[-1]]
        titles = ["grow: origin", "grow: spreading", "grow: tumour (trigger)",
                  "reverse: wave", "reverse: collapsing", "reverse: cleared"]
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        for ax, fr, ti in zip(axes.ravel(), pick, titles):
            ax.imshow(render_rgb(fr)); ax.set_title(ti, fontsize=10); ax.axis("off")
        fig.suptitle("Tumour grows cell-to-cell → recognition-gated death reverses it along the same routes", fontsize=12)
        fig.tight_layout(); fig.savefig(OUT_DIR / "grow_reverse_timeline.png", dpi=110); plt.close(fig)
        # population curve
        fig2, ax = plt.subplots(figsize=(8, 4))
        ax.plot(range(len(gsizes)), gsizes, color="#c0392b", label="GROW: tumour size")
        ax.plot(range(len(gsizes), len(gsizes)+len(rsizes)), rsizes, color="#222", label="REVERSE: tumour size")
        ax.axvline(len(gsizes), ls="--", color="gray"); ax.set_xlabel("step"); ax.set_ylabel("cancer cells")
        ax.legend(); ax.set_title("Tumour size: grows, then collapses when the gated death wave fires")
        fig2.tight_layout(); fig2.savefig(OUT_DIR / "grow_reverse_curve.png", dpi=110); plt.close(fig2)
        log.info("figures → runs/grow_reverse/grow_reverse_timeline.png + grow_reverse_curve.png")
        # animated GIF
        try:
            from PIL import Image
            imgs = [Image.fromarray(render_rgb(f)).resize((322, 322), Image.NEAREST) for f in allframes]
            imgs[0].save(OUT_DIR / "grow_reverse.gif", save_all=True, append_images=imgs[1:],
                         duration=180, loop=0)
            log.info("animation → runs/grow_reverse/grow_reverse.gif (%d frames)", len(imgs))
        except Exception as e:
            log.warning("GIF skipped (%s: %s)", type(e).__name__, e)
    except Exception as e:
        log.warning("figure skipped (%s: %s)", type(e).__name__, e)

    log.info("results → runs/grow_reverse/grow_reverse_results.json")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
