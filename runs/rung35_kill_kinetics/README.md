# RUNG-35 — KILL-CIRCUIT KINETICS + LEAK-CORRELATION (stress-testing RUNG-34's load-bearing residual)

RUNG-34 concluded the 2-input AND-gated self-destruct gives a **therapeutic index ~10¹⁰** (normal false-death
4×10⁻¹¹) and flagged ONE load-bearing assumption: **leak independence**. This rung tests that assumption — and a
second hidden one — and **honestly downgrades** the headline number, then shows what wins it back. `scripts/69`.

## Two things RUNG-34 assumed (and this rung breaks)
1. **Leak is a uniform sub-threshold LEVEL.** RUNG-34 multiplied continuous leak activities (L^N) and passed them
   through the steep Hill switch → 4×10⁻¹¹. That only holds if every normal cell carries a *little* recognizer
   activity that the bistable switch filters. But off-target CRISPR cuts and toehold mis-triggers are **all-or-none
   per cell**: leak = a *fraction* of normal cells that **fully** mis-fire. In those cells the effector is FULL, so
   the Hill switch does **not** save you → false-death ≈ P(all N fire).
2. **Leaks are independent.** A shared stressor (cell-state, common delivery vehicle) makes recognizers mis-fire
   *together*.

## Result A — kinetics (reduced eARM bistable caspase switch, irreversible MOMP)
`dC/dt = u + Vmax·Cʰ/(Kʰ+Cʰ) − k_deg·C`, h=4. At zero input the switch is **bistable** (survival C≈0 **and** death
C≈0.92 both stable) → once it snaps high it stays high = **MOMP irreversibility** (RUNG-1).
- **Dynamic input threshold** u_crit (saddle-node) = **0.192** on the *drive* (RUNG-34's K=0.30 was on the
  effector *level* — different quantity).
- Mutant AND-drive (0.81) **commits in ~0.7** (1/k_deg ≈ 0.5–2 h); sustained sub-threshold graded leak (0.0025)
  **never commits**.
- **Temporal-coincidence filter:** a supra-threshold mis-fire must be **sustained ≥ T_min ≈ 0.73** (1/k_deg) to
  commit. A *brief* full-amplitude leak self-clears and does **not** kill. Kinetics add a *temporal* AND on top of
  the molecular AND.

## Result B — honest false-death (all-or-none leak + correlation ρ, one-factor Gaussian copula)
| N | L | ρ=0 (indep) | ρ=0.3 | ρ=0.6 | ρ=1 (full) | RUNG-34 graded-Hill (ref) |
|---|---|---|---|---|---|---|
| 2 | 0.05 | **2.5×10⁻³** | 7.1×10⁻³ | 1.6×10⁻² | 5.0×10⁻² | 4×10⁻¹¹ |
| 3 | 0.05 | 1.3×10⁻⁴ | 1.7×10⁻³ | 7.7×10⁻³ | 5.0×10⁻² | 1×10⁻¹⁷ |
| 4 | 0.05 | 6.3×10⁻⁶ | 5.7×10⁻⁴ | 4.7×10⁻³ | 5.0×10⁻² | 4×10⁻²⁴ |

**The N=2 AND-gate alone is NOT safe** under the conservative model: ~2.5×10⁻³ false-death even when leaks are
independent (8 orders worse than RUNG-34's 4×10⁻¹¹), rising to ~5×10⁻² as ρ→1. **Required inputs** at 5% leak to
hit a target: 1e-4 → N=4 (ρ=0) but N=7 (ρ=0.3) and *infeasible* for ρ≥0.6; 1e-9 → N=7 (ρ=0) and infeasible at any
correlation. **Correlation, not N, is the binding constraint.**

## Result C — kinetics win back what the probabilistic model lost
A coincident full mis-fire commits only if the overlap persists ≥ T_commit → reduced marginal
`L_eff = L·(1 − T_commit/T_leak)`:
| T_leak/T_commit | L_eff | N=2 (ρ=0 / 0.3) | N=3 (ρ=0 / 0.3) |
|---|---|---|---|
| ≤1 | 0 | 0 / 0 (full rescue) | 0 / 0 |
| 2 | 0.025 | 6.3×10⁻⁴ / 2.4×10⁻³ | 1.6×10⁻⁵ / 4.3×10⁻⁴ |
| 5 | 0.040 | 1.6×10⁻³ / 5.0×10⁻³ | 6.4×10⁻⁵ / 1.1×10⁻³ |
| 20 (sustained) | 0.048 | 2.3×10⁻³ / 6.6×10⁻³ | 1.1×10⁻⁴ / 1.6×10⁻³ |

If recognizer leaks are made **transient** (T_leak ≲ a few × T_commit) **and independent**, kinetics drive
effective false-death down several more orders — e.g. transient + N=3 + ρ=0 → ~10⁻⁵.

## Verdict
**Conditional feasibility — the 10¹⁰ index was optimistic.** A 2-input AND-gate alone is *not* safe under
all-or-none leak (~10⁻³, not 10⁻¹¹). Safety is recoverable but **conditional** on four now-explicit, testable
levers: (i) per-recognizer leak amplitude L, (ii) leak **duration** ≪ commitment time (kinetic filter), (iii) leak
**correlation ρ** near zero (independent delivery / cell-state triggers), (iv) **N≥3** if ρ can't be driven low.
The single most important wet measurement is the **recognizer-leak correlation in stressed normal cells** — it
decides which row of Result B you live on.

## Honest residuals
- Dimensionless time (1/k_deg ≈ 0.5–2 h, stated not fitted); reduced 1-variable caspase switch (captures the
  irreversible bistable snap of eARM/RUNG-1, not the full tBid/Bax/Apaf-1 cascade).
- The Part-C transience model (`L_eff = L·(1−T_commit/T_leak)`) is a faithful but reduced surrogate for the
  sustained-coincidence probability, not a fitted leak-kinetics measurement.
- This quantifies the **specificity logic + its dynamic robustness**, not a built circuit. Wet steps: measure L,
  the leak time-constant, and ρ across recognizers in primary normal cells.

*Result: `kill_kinetics.json`. Tests RUNG-34 (C9). New ledger row C10; CAPSTONE §14. Selftest-gated.*
