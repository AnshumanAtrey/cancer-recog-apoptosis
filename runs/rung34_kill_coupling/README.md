# RUNG-34 — KILL-COUPLING: recognition → self-destruct (the step that makes the internal key a root-kill)

The recognizers (toehold R25 / allele-CRISPR R27 / DhdR-2HG R33) *detect* the cancer mutation. This rung
designs + quantifies the **coupling to death** — turning a detector into an autonomous self-destruct that kills
the cell via its **own apoptosis**, MHC-free. `scripts/68`.

## Architecture (modular self-destruct)
```
  recognizer × N  ──►  AND-gate  ──►  apoptosis effector
   toehold RNA (non-wobble) / allele-CRISPR (wobble, R27d) / DhdR D-2-HG (IDH1, R33)
   AND = split death effector (each half gated by one recognizer) | serial toeholds | dual-guide-required txn
   effector = iCasp9 (clinical CAR-T safety switch → engineered apoptosis-on-demand is PROVEN) / casp-3 / BAX
```
**Why AND, not a single recognizer:** every recognizer leaks (mis-fires in some normal cells). Wired alone to a
death switch, a 5% leak kills ~5×10⁹ of ~10¹¹ normal cells = lethal. AND-gating two **independent** recognizers
makes the joint false-positive ~L₁·L₂, and the cooperative apoptosis threshold (bistable MOMP/caspase switch,
eARM R1) filters sub-threshold leak. (RUNG-22's escape-logic, applied to KILLING.)

## Result (model: AND-gate effector → Hill apoptosis K=0.30, n=5; grounded recognizer leaks)
| design | normal false-death | tumour-kill | therapeutic index |
|---|---|---|---|
| N=1, 5% leak | **1.3×10⁻⁴** (unsafe @10¹¹ cells) | 0.996 | 7.7×10³ |
| **N=2, 5% leak** | **4×10⁻¹¹** | **0.993** | **2.5×10¹⁰** |
| N=3, 5% leak | 1.3×10⁻¹⁷ | 0.988 | 7.9×10¹⁶ |

- **A single recognizer cannot safely gate a kill; a 2-input AND-gate can** (false-death from 10⁻⁴ → 10⁻¹¹).
- The apoptosis threshold *amplifies* the mutant signal (0.81 effector → 0.99 kill), so at strong recognizers
  the kill/specificity trade-off is mild. For **weaker** recognizers (mutant-ON < ~0.8) or **N≥3**, single-cell
  kill drops — but the resistance-agnostic **bystander death wave** (R13/R21, cures across 4–13% escapees)
  clears the cells that don't fully fire, so partial single-cell kill suffices at the tumour level.

## Verdict
**Feasible root-kill.** A 2-input AND-gated self-destruct on our measured-leak recognizers kills mutant cells
via their own apoptosis with a ~10¹⁰ therapeutic index — recognizing the cancer at its *root* (the somatic
mutation / oncometabolite) and triggering the cell's own death, no immune system, no MHC. This completes the
internal key from *detector* → *therapeutic logic*.

## Honest residuals (the load-bearing one named first)
- **Independence of recognizer leaks is THE assumption.** If a stressed normal cell mis-fires several recognizers
  together (correlated failure; shared delivery / cell-state), the joint leak → L (not L²) and the AND margin
  collapses. **This is the key thing to test** (measure recognizer-leak correlation in normal cells).
- The death module must be **non-leaky** (basal iCasp9 must not trigger); the AND-reconstitution must be tight.
- **Kinetics** — does the effector cross the apoptosis threshold before the cell adapts? (next: a full eARM ODE
  of the coupled circuit, building on R1, not just this steady-state logic.)
- Multi-component **delivery** into the cell; this is a SPECIFICITY-LOGIC result, not a built/validated circuit.

*Result: `kill_coupling.json`. Connects R25/27/33 (recognizers) → R22 (AND-logic) → R1 (apoptosis) → R13/21
(bystander wave). New ledger row C9; CAPSTONE §13.*
