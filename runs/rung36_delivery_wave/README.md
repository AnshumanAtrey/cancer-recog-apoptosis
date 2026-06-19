# RUNG-36 — DELIVERY SEED-FRACTION × BYSTANDER-WAVE: does the wave decouple the kill from delivery?

RUNG-35's safety argument leaned on an **unquantified escape hatch**: *"cells that don't fully fire are cleared by
the bystander wave, so partial single-cell kill suffices at the tumour level."* And every rung names **delivery**
("the circuit only reaches some cells") as an open residual. This rung quantifies both at once: **if the engineered
circuit is delivered to only a seed fraction f of the tumour, and fires with per-cell probability p_kill (the R35
partial kill), does the self-propagating apoptosis wave still clear the tumour — and what is the minimum delivery
f\* for clearance, as a function of wave strength?** `scripts/70`.

## Model (stochastic lattice CA ≡ bond-percolation ignition; the R12P/B percolation framing made therapeutic)
ALIVE/DYING/DEAD lattice. **Delivery:** each cell gets the circuit w.p. f. **Ignition (t=0):** delivered cells
fire w.p. p_kill → DYING. **Wave:** a DYING cell recruits each ALIVE neighbour w.p. b (resistance-agnostic
bystander: gap-junction / Fas-FasL, R13/R21) — a recruited cell needn't be delivered-to, so a partial seed can
clear the whole tumour. Each cell DYING one step, trying each neighbour once ≡ **bond percolation, prob b** (2D
square-lattice threshold **p_c = 0.5**). Optional **resistance r** (R18/R21 escapees the wave can't recruit).
Clearance target = 99% (the residual <1% = the R19/R21 regrowth seed). Selftest-gated on the percolation limits.

## Result 1 — clearance from a SMALL seed (deliver to only 5%) vs wave strength b
| b | 0.0 | 0.35 | 0.45 | **0.50 (p_c)** | 0.60 | 0.70 | 0.80 | 0.90 |
|---|---|---|---|---|---|---|---|---|
| cleared @ f=0.05 | 0.05 | 0.36 | 0.69 | **0.83** | 0.95 | **0.99** | 0.998 | 1.00 |

Below p_c the wave dies locally (cleared ≈ f). Above p_c a 5%-seed ignites a tumour-spanning death cluster: at
b≈0.7 the wave clears **~99%** of the tumour from a **5% delivery**.

## Result 2 — minimum delivery f\* to clear ≥99% (the decoupling headline)
| b | 0.35 | 0.45 | 0.55 | 0.60 | 0.70 | 0.80 |
|---|---|---|---|---|---|---|
| f\* (p_kill=1.0) | 1.00 | 0.90 | 0.90 | 0.70 | 0.20 | **0.002** |
| f\* (p_kill=0.5) | none | none | none | none | 0.35 | **0.002** |

**f\* collapses by ~3 orders as b crosses well above p_c.** With *partial* single-cell kill (p_kill=0.5) clearance
is **impossible at any delivery** until the wave is super-critical (b≥0.7) — the precise condition under which
R35's escape hatch is true.

## Result 3 — the hard limit: wave-resistant escapees
| r (resistant frac) | 0 | 0.02 | 0.05 | 0.10 | 0.20 |
|---|---|---|---|---|---|
| cleared (b=0.8, f=0.05) | 0.998 | 0.978 | 0.948 | 0.896 | 0.786 |

Clearance caps at **≈ 1 − r** regardless of wave strength — the wave alone cannot clear cells it can't recruit.
This is exactly why R21 needed a **resistance-agnostic 2nd killer (NK)**.

## Result 4 — 3D (real tumour geometry) is *more* forgiving
Simple-cubic bond p_c ≈ 0.25 ≪ 0.5. Cleared @ f=0.05 reaches 0.93 by b=0.4, and the regime flips super-critical at
b>0.25 — so in 3D the wave is **easier** to make super-critical and the delivery requirement is even looser than
the 2D bound.

## Verdict
**The bystander wave formally decouples the kill mechanism from delivery — conditionally.** (1) Super-critical wave
(b>p_c) → minimum delivery collapses toward a few percent; **delivery stops being the bottleneck**. (2) Sub-critical
wave → near-complete delivery required; **delivery IS the bottleneck**. (3) R35's partial-kill escape hatch holds
**only** above p_c. (4) Wave-resistant escapees cap clearance at ~1−r → a resistance-agnostic 2nd killer (R21) is
mandatory for the last fraction. (5) 3D lowers p_c → more forgiving. **The engineering target shifts from "deliver
to every cell" (near-impossible) to "make the bystander coupling super-critical + add an agnostic killer for the
resistant fraction" — two separable, individually-tractable problems.**

## Honest residuals
- Lattice CA (no vasculature / 3D diffusion geometry / immune kinetics / cell motility); f, p_kill, b, r are
  **effective** probabilities, not measured efficiencies. The recognition→caspase-8 transduction potency that sets
  **b** is the same wet-lab residual R1/R13 named.
- Bond-percolation equivalence assumes a dying cell tries each neighbour once; the graded/kinetic bystander signal
  (R13) would shift the effective threshold.
- Robust claim = the **percolation collapse of f\* at b=p_c** and the **1−r resistance ceiling**, *not* an absolute
  cure %. Connects A5/A7 (death wave) + R13/R21 (wave mechanics) + R35 (partial kill).

*Result: `delivery_wave.json`. New ledger row A8; CAPSTONE §15. Selftest-gated; headline metric audited (the first
0.999 "cure" cutoff was a finite-size threshold artifact → recalibrated to a 99% clearance target).*
