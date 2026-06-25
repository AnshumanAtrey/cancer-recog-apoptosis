# RUNG-42 — THE EVOLUTIONARY TRAP (H10): does a collateral-sensitivity double-bind CLOSE the escape route R19 left open?

R19 proved the bare recognition kill (pressure A) escapes at clinical size: some cell has already lost the
recognised driver-neoantigen, and that escapee sweeps unless a resistance-AGNOSTIC 2nd killer mops it up
(R21 NK / R36 wave). H10 (Shriya's §6.3) asks a sharper variant — the **Gatenby double-bind**: set the 2nd
pressure (B = a standard therapy) so the *very mutation that escapes A makes the cell HYPERSENSITIVE to B*
(collateral sensitivity). Escape A → die to B; don't escape → die to A. A deliberate red-team: does this
close the hole, and where does it leak? `scripts/76`.

## Model — R19's growth/standing-variation lattice + SUSTAINED dual pressure
S(susceptible) / RT(coupled escape) / RF(decoupled escape) on a 130² lattice. Growth seeds Luria-Delbrück
standing variation: each S division escapes S→RT w.p. `mu(1-p_decouple)`, S→RF w.p. `mu·p_decouple`. Then
two **persistent** pressures (H10's "fitness landscape under dual selective pressures"):
- **A** = the installed recognition self-destruct circuit — kills any displayed-mutation S cell at `k_A`
  every step (persistent, NOT R19's one-shot wave; the wave was R36's question). RT/RF escaped recognition.
- **B** = systemic standard therapy — kills **RT at `b_trap` (HIGH = collateral sensitivity)**, RF/S at
  `b_base` (baseline). RT also divides slower (cost `s`) — the price of dropping the driver.
- **RT** = escaped A *and* paid the price (coupled: lost/altered the essential driver → cost + B-sensitive).
- **RF** = escaped A by a rare **DECOUPLING** mutation (evade recognition *without* losing the driver or
  sensitising to B) → cost-free, baseline B-sensitivity → the trap's BLIND SPOT.

## Result 1 — the trap CLOSES the coupled escape (and the lattice validates the analytic backbone)
At mu·N0 = 15 founders, R19 (no trap) cures 0/32. Trap ON, sweeping the decoupling fraction:
| p_decouple | 0 | 1e-3 | 1e-2 | 0.05 | 0.2 | 0.5 | 1.0 |
|---|---|---|---|---|---|---|---|
| **P(cure)** | **1.00** | 0.97 | 0.78 | 0.44 | 0.03 | 0.00 | 0.00 |
| analytic exp(−μ·p_decouple·N0) | 1.00 | 0.985 | 0.859 | 0.467 | 0.048 | 0.00 | 0.00 |

The sim **tracks the analytic L_eff = μ·p_decouple·N0** across the whole range (the rigor check, as R19
tracked Luria-Delbrück). So the trap multiplies the **curable-tumour ceiling by 1/p_decouple** (pd=1e-2 →
100×, pd=1e-4 → 10,000×). And **every single escape was via the decoupled route** (coupled column = 0 in all
32 reps at every pd): the trap shuts the coupled escape completely.

## Result 2 — the NEW load-bearing residual: the decoupling fraction p_decouple
The trap is only as tight as the **coupling between the recognised feature (A) and the essential driver
whose loss sensitises to B.** A rare decoupling mutation reopens a cost-free escape the trap can't see.
**p_decouple is the H10 analogue of R35's leak-correlation ρ** — the #1 new wet measurement. The design
lever that minimises it: **target the driver mutation ITSELF as the recognition handle** (KRAS-G12D: the
neoantigen *is* the oncogenic driver) → losing the epitope = losing the driver, the tightest possible
coupling. Ties R22 (essentiality ranking) + C8 (KRAS-G12D = MS/x-ray/T-cell gold-standard).

## Result 3 — the trap is EFFICACY-DOMINATED by the resistance-agnostic killer we already have
Collateral-sensitive trap vs R19/R21's resistance-**agnostic** 2nd killer, at **equal kill budget** (b=0.3):
| p_decouple | collateral trap P(cure) | agnostic bystander P(cure) |
|---|---|---|
| 0 | 1.00 | 1.00 |
| 0.01 | 0.78 | **1.00** |
| 0.1 | 0.22 | **1.00** |

Whenever decoupling exists, the **agnostic killer strictly wins** — it isn't blind to the decoupled route.
The trap's *only* edge is a **toxicity/specificity** argument this equal-budget model doesn't credit:
collateral sensitivity may permit a high kill on escaped cells at a dose **safe for normal tissue**, where an
agnostic drug is toxicity-capped. → **H10 is a toxicity-sparing COMPLEMENT, not the resistance answer**; the
robust solution stays **R21 (agnostic cross-kill) + R22 (multi-target essential drivers)**.

## Result 4 — fitness cost ALONE is insufficient
"Target an essential driver so escape is costly" does **not** close escape on its own: cost-only (B off),
even at s=0.95, cures 0/32. A slower escapee still escapes — **costliness ≠ lethality**. The cost only
lowers the bar the *active* 2nd kill (B) must clear (B's threshold here ≈ 0.2). Only s=1 (escape = certain
death — the feature strictly required for viability) closes it without B.

## Verdict — H10: REAL but DOMINATED, with a new blind spot (a genuine red-team result)
The evolutionary trap is a real mechanism that **closes the coupled escape route R19 left open** (1/p_decouple
ceiling multiplier, analytic-validated). But the red-team found **two honest negatives on its cleverness**:
(1) a **new residual** — the decoupling fraction p_decouple, which it's blind to; and (2) it is
**efficacy-dominated** by the resistance-agnostic 2nd killer already in the chain, edged only by a
toxicity argument. Net: it **sharpens** rather than replaces the resistance answer — drive p_decouple→0 by
recognising the driver mutation itself, keep R21's agnostic cross-kill as the robust backstop, and use
collateral sensitivity opportunistically to spare normal tissue.

## Honest residuals
- 2D lattice CA, not a tumour (no microenvironment/immune/3D/drug gradients). `mu` = effective lumped
  per-division escape prob; `b_trap/b_base/s` = proxies for a real collateral sensitivity whose magnitude is
  a wet residual; the equal-budget comparison ignores the toxicity ceiling that is the trap's actual case.
- **WET residual (new, #1 for H10):** the **decoupling fraction p_decouple** — how often a tumour can evade
  the recognised feature *without* losing the driver / sensitising to B. Decides whether the trap holds.
  Minimised by targeting the driver mutation itself; measurable as the escape-mutation spectrum under the
  recognition pressure. Joins R35's leak-correlation ρ and R37's gate on the residual list.

*Result: `rung42_evolutionary_trap.json`, figure `.png`. Extends R19 (escape race); composes with R21/R22
(resistance answer). New ledger row A12 + H10 cross-ref; CAPSTONE §21. Selftest 13/13; analytic validation +
cost-only + agnostic controls baked into the run so the JSON is reproducible from the script. Verdict written
AFTER the audit, not before — the "agnostic dominates" and "cost insufficient" negatives were the surprises.*
