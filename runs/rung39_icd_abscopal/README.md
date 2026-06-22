# RUNG-39 — DANGER SIGNAL / ICD (Shriya's hypothesis H7): does the self-destruct reach METASTASES?

Every prior rung clears the **local** tumour (recognition → AND-gate kill → bystander wave → containment). But the
wave is contact/percolation-based — it cannot cross to a **distant micrometastasis that never received the
circuit**, and metastasis is what kills patients. The only mechanism that reaches untreated distant deposits is
**systemic anti-tumour immunity raised by immunogenic cell death (ICD)** — the abscopal effect. This tests
Shriya's H7 (Danger Signal / DAMPs), the most clinically decisive open hypothesis. `scripts/73`.

## The tension
Plain apoptosis is frequently **tolerogenic** (non-immunogenic, even immunosuppressive). So a "clean" internal
self-destruct could cure the primary and **leave the metastases**. H7 asks: can the kill be made immunogenic
enough to flip *local* killing into *systemic* clearance — and what does that require?

## Model (coupled ODEs)
Local clearance (R36) → a transient burst of dying-cell antigen; a fraction **`i`** (ICD immunogenicity: 0 =
tolerogenic apoptosis, 1 = fully immunogenic, e.g. pyroptotic/necroptotic) releases DAMPs. **Thresholded priming**
(DCs need enough danger signal to mature → Hill) raises a **persistent (memory)** T-cell pool that kills distant
deposits; **extinction floor** (a deposit driven below ~1 cell is eradicated). Immunosuppression `supp` scales the
T-cell kill. Selftest-gated (tolerogenic escapes, immunogenic cures, threshold exists, heavy suppression escapes
even at i=1).

## Result — abscopal cure (distant metastasis, M0=0.2) vs immunogenicity × immunosuppression
| i \ suppression | none | ×0.7 | ×0.5 | ×0.25 (heavy) |
|---|---|---|---|---|
| 0.0 (tolerogenic) | escape | escape | escape | escape |
| 0.2 | escape | escape | escape | escape |
| 0.3 | **CURE** | escape | escape | escape |
| 0.4 | CURE | CURE | CURE | escape |
| 0.6–1.0 | CURE | CURE | CURE | **escape** |

**Minimum immunogenicity i\* for abscopal cure** rises with suppression: 0.26 (none) → 0.39 (×0.5) → 0.64 (×0.35)
→ **none possible** under heavy suppression (×0.25).

Plain tolerogenic apoptosis (i≈0.05): metastasis grows back → **local-only cure**. Engineered-immunogenic kill
(i≈0.8): metastasis eradicated → **abscopal cure**.

## Verdict — three regimes
1. **Tolerogenic apoptosis → local-only cure.** Below the priming threshold no anti-tumour immunity is raised;
   distant metastases that never got the circuit **grow back**. The primary is cured but the patient is not.
2. **Immunogenic death (above threshold) + manageable suppression → systemic (abscopal) cure.** ICD raises a
   circulating, memory T-cell response that eradicates untreated deposits.
3. **Heavy immunosuppression → immunity alone fails even at max immunogenicity** → exactly where the
   resistance-agnostic **NK arm (R21, for MHC-dark deposits)** and **checkpoint-release** are required.

**DESIGN REQUIREMENT (new, named):** the internal self-destruct must be **engineered for immunogenic cell death**
(calreticulin exposure / HMGB1 / ATP release, or a pyroptotic/necroptotic flavour) — **not clean apoptosis.** This
is what turns Shriya's local *"destroy itself from within"* into a **whole-body cure**, and it adds a **third
clearance layer**: local wave (R36) + systemic ICD-immunity (R39) + NK tail (R21). H7 is therefore not "a cheap
score" — it is the pivotal **local-vs-systemic** determinant of whether this could ever be a real cure.

## Honest residuals
- Reduced immuno-oncology ODE: no explicit DC/lymph-node trafficking, T-cell exhaustion dynamics, antigen
  spreading, or MHC-loss escape (the last is precisely why R21's NK arm is needed for dark metastases). `i, kp,
  kkill, supp` are effective, not fitted. Memory persistence (slow T-cell decay) and the extinction floor are
  modelling assumptions (biologically standard, but assumptions).
- **THE KEY WET RESIDUAL:** the actual ICD immunogenicity of *this specific* kill mechanism is **unmeasured** —
  whether the engineered self-destruct exposes calreticulin / releases HMGB1 / triggers pyroptosis must be
  measured in cells.
- Robust claim = the **threshold structure** (tolerogenic fails, immunogenic-above-threshold clears) + the
  **suppression dependence**, not absolute timings.

*Result: `icd_abscopal.json`. Tests Shriya's H7; connects R36 (local wave) + R21 (NK) into systemic clearance.
Ledger H7 → DONE; CAPSTONE §18. Selftest-gated; model audited (transient-immunity regrowth bug caught → fixed
with memory persistence + extinction floor before any conclusion).*
