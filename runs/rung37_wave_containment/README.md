# RUNG-37 — WAVE CONTAINMENT: the coupled R35×R36 safety model (a hole R36 hid)

R36 proved the bystander wave must be **super-critical** (b>p_c) to clear the tumour from a partial seed — but it
ran the wave on a **pure-tumour lattice**. The tumour is really a small island in a **sea of normal tissue**, and
the wave is **resistance-agnostic** (it kills a neighbour regardless of type). So the question R36 hid: *does that
same super-critical wave run away through normal tissue?* This rung couples R36's percolation with R35's leak on a
**two-tissue lattice** (central tumour disk = 12.6%, rest = normal). `scripts/71`.

## Model
Wave bond-prob is set by the **receiving** cell's type: `b_t` into tumour cells (free spread = super-critical),
`b_n` into normal cells (the **recognition gate** — the bystander death completes only if the receiver's own
mutation-sensor is primed, so the death signal is damped into normal cells). Tumour ignites by delivery×kill;
normal cells falsely ignite at `p_ig_n` = the R35 leak. Selftest-gated on the percolation limits.

## Result 1 — an UNGATED wave is a catastrophe by boundary spillover *alone*
| p_ig_n (R35 leak) | 0 | 1e-5 | 1e-3 | 1e-2 |
|---|---|---|---|---|
| tumour cleared | 0.999 | 0.999 | 0.999 | 0.999 |
| **normal KILLED** | **0.998** | 0.998 | 0.998 | 0.998 |

With `b_n = b_t = 0.8`, normal tissue dies **even at zero leak** — the tumour-clearing wave simply crosses the
boundary and percolates through everything. A normal false-fire only adds ignition. **The wave does not know the
tumour boundary.**

## Result 2 — a RECOGNITION-GATED wave contains it (but only well below p_c)
| b_n | 0.0 | 0.20 | 0.35 | 0.45 | **0.50 (p_c)** | 0.60 | 0.80 |
|---|---|---|---|---|---|---|---|
| tumour cleared | 0.998 | 0.998 | 0.998 | 0.999 | 0.999 | 0.999 | 0.999 |
| normal killed | 0.001 | 0.006 | 0.025 | 0.159 | **0.609** | 0.944 | 0.998 |

Tumour clearance is **flat** (b_t carries it); normal death stays small **only while b_n is sub-critical** — and
really only **well below** p_c (it jumps 0.025 → 0.16 → 0.61 across b_n = 0.35 → 0.45 → 0.50). **Safety requirement:
the recognition gate must hold b_n well under 0.5.**

## Result 3 — the wave AMPLIFIES the R35 leak (the coupling that matters)
| p_ig_n (R35) | normal killed (gate b_n=0.35) | amplification | normal killed (ungated b_n=0.8) | amplification |
|---|---|---|---|---|
| 1e-5 | 4.7e-4 | **47×** | 0.62 | 62,000× |
| 1e-4 | 2.8e-3 | 28× | 0.998 | 10,000× |
| 1e-3 | 1.6e-2 | 16× | 0.998 | 1,000× |

The wave multiplies the leak by the **mean normal-tissue death-cluster size**. A sub-critical gate **bounds** it
(~tens-fold at b_n=0.35; →1× only as b_n→0); an ungated wave **diverges** (a 1e-5 leak → total death). **So the
effective normal false-death = (R35 leak) × (normal cluster size) — R35's leak target must be *tightened* by that
factor**, and is only safe behind a gate well below p_c. The recognition gate is **as load-bearing as the leak.**

## Result 4 — the unavoidable collateral rim
Even fully contained (zero leak), the tumour wave kills a thin **rim** of normal cells at the boundary:
~**1.2%** at b_n=0.35, ~**6%** at b_n=0.45 — an intrinsic "surgical margin" that grows with b_n. Higher b_n clears
the tumour *edge* better but widens the normal rim; the gate sets where on that trade-off you sit.

## Verdict
**The wave must be bistable across tissue type: super-critical among recognised (tumour) cells, sub-critical among
unrecognised (normal) cells — i.e. a RECOGNITION-GATED bystander signal with b_n well below p_c.** This is a **new,
named safety requirement that R36's pure-tumour percolation hid**, and it sits beside R35's leak-correlation as the
**second load-bearing residual** of the kill-coupling. The same property that makes the wave clear the tumour from
a tiny seed (super-criticality) makes it lethal in normal tissue unless gated — and the gate must be tight, because
the wave *amplifies* the underlying leak by the cluster size.

## Honest residuals
- 2D lattice CA with a **sharp** tumour/normal boundary (real tumours are infiltrative/graded → rim and
  containment are softer); `b_t, b_n, p_ig_n, f` are effective probabilities, not measured.
- The recognition gate is modelled as a **reduced bond-prob into normal cells**, not a molecular mechanism.
- Robust claim = the **b_n < p_c containment requirement** + the **leak-amplification by cluster size**, not
  absolute fractions. **Wet residual:** the actual recognition-gating of the bystander death factor (does it
  require the receiver's own mutation-sensor?) — the molecular basis of b_n ≪ b_t.

*Result: `wave_containment.json`. Couples R35 (leak) + R36 (percolation). New ledger row A9; CAPSTONE §16.
Selftest-gated; narrative audited against the numbers (the "amplification ~1×" first draft was corrected to the
mean-cluster-size law).*
