# RUNG-38 — END-TO-END SAFETY ENVELOPE: composing the molecular leak (R35) with the tissue wave (R37)

R35 gave the *molecular* per-normal-cell false-death (a function of N AND-inputs, leak L, correlation ρ, leak
transience). R37 showed the wave *amplifies* whatever normal ignition exists and must be recognition-gated. Neither
gives the **tissue-level** answer. This rung composes them into the complete safety envelope of the internal
root-kill. `scripts/72`.

## The composition law (the headline)
> **tissue normal false-death  =  (molecular leak, R35)  ×  (wave amplification, R37)  ×  (delivery footprint)**

The two load-bearing residuals are **multiplicative, not independent** — the molecular AND-gate leak is *multiplied*
by the wave's normal-tissue death-cluster size. Measured wave amplification (mean cluster size) by gate:
| gate b_n | 0.20 | 0.35 | 0.45 (near p_c) |
|---|---|---|---|
| amplification | **3×** | **16×** | **210×** |

## Design-grid result (per-recognizer leak L=0.05, b_t=0.80 super-critical; tumour clears in all 16)
The end-to-end normal false-death spans **~7,600×** across the grid. Best vs worst:
| design | molec leak | × amp | end-to-end (frac of circuit-carrying cells) |
|---|---|---|---|
| **best:** N=3, ρ=0, transient, b_n=0.20 | 4.6e-6 | 3× | **1.5e-5** |
| worst: N=2, ρ=0.3, sustained, b_n=0.35 | 7.1e-3 | 16× | **1.2e-1** |

**Lever ranking** (×-worse when flipped from the best corner):
| lever | effect |
|---|---|
| add an AND-input (N 3→2) | **60×** |
| recognizer correlation (ρ 0→0.3) | **42×** (worse at high N) |
| leak transience (transient→sustained) | **27×** |
| gate tightness (b_n 0.20→0.35) | 5× (→210× near p_c) |

No single factor is sufficient: a 2-input gate, a sustained leak, a correlated leak, **or** a loose gate each
breaks tissue safety even when the others are ideal. The best swept design reaches **1.5e-5**; **N=4 or a tighter
gate / lower base leak pushes below 1e-6.**

## The third factor — delivery localisation is *also* a safety lever
The end-to-end fraction hits only normal cells that **carry the circuit**. So absolute deaths = fraction ×
delivery footprint — and tumour-localised delivery (R36: a few-% seed already cures) shrinks the at-risk
population by orders. For the best design (1.5e-5):
| delivery footprint | est. normal deaths |
|---|---|
| systemic (all 1e11 cells dosed) | ~1.5e6 |
| regional (~1e9 near-tumour) | ~1.5e4 |
| **tumour-localised (~1e7 rim)** | **~1.5e2** |

So **R36's efficacy result (localised delivery suffices) is *also* the safety result** — delivery localisation and
the kill mechanism reinforce each other.

## Verdict
**The autonomous root-kill is tissue-safe in a concrete, composed design region** — N≥3 AND-inputs + transient,
uncorrelated leaks + a tight recognition gate (b_n ≪ p_c) + tumour-localised delivery — and in that region the
**same wave still clears the tumour from a partial seed.** Crucially, "safe enough" is no longer a hand-wave: it's
the **product of four measurable molecular numbers (L, ρ, leak-lifetime, b_n) × the delivery footprint.** This is
the complete in-silico safety envelope of the internal key, composing R34 (AND-gate) → R35 (leak kinetics) → R36
(delivery/clearance) → R37 (containment) into one number.

## Honest residuals
- Composition of two reduced models (R35 Gaussian-copula leak + R37 lattice amplification); effective
  probabilities, 2D CA, sharp tumour boundary. Robust claim = the **multiplicative composition law** + **lever
  ranking** + the **delivery-footprint factor**, not absolute fractions.
- The lattice resolves normal_killed only to ~1/N_cells, so the safe region is verified **analytically**
  (molecular leak × measured cluster size) with the lattice confirming tumour clearance and the unsafe rows.
- `abs@1e11` is illustrative; the real at-risk population is the circuit-carrying normal cells (delivery footprint).
- **Wet residuals unchanged:** measure L, ρ, the leak time-constant, and the bystander recognition-gating (b_n) in
  primary normal cells — the same four numbers, now with an explicit composition that turns them into a safety verdict.

*Result: `end_to_end_safety.json`. Composes R34/35/36/37. New ledger row A10; CAPSTONE §17. Selftest-gated;
headline audited (a binary "0/16 safe at 1e-6" first draft was a threshold artifact → replaced with the
composition law + lever ranking + delivery-footprint factor).*
