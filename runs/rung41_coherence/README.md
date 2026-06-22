# RUNG-41 — COHERENCE STRESS-TEST: does R39 (immunogenic death) break R37 (wave containment)?

After R34–R40 the design carries a pile of requirements. The rigorous move is not a 41st feasibility rung — it is
to ask whether the requirements **contradict**. The sharpest candidate: **R39 demands immunogenic (inflammatory,
DAMP-releasing) death to reach metastases; R37 demands a tightly-contained, recognition-gated wave to spare normal
tissue.** Immunogenic death is inherently inflammatory → raising immunogenicity *i* (for metastases) should raise
the effective normal-tissue coupling (breaking containment). Are they compatible? `scripts/75`.

## Model — compose R37's two-tissue lattice with an ICD-inflammation term
`b_n_eff(i) = b_n_base + κ·i`, where `b_n_base` = the recognition-**gated** (specific, sub-critical) coupling into
normal cells (R37), and `κ·i` = the **non-gated innate inflammatory** coupling from immunogenic death (∝ R39's
immunogenicity *i*; κ = how locally lytic/inflammatory the death mode is). The *adaptive* systemic immunity ICD
raises is antigen-specific (safe for normal cells) — it's the *innate local* inflammation (κ·i) that threatens
containment. Metastases clear iff `i ≥ i*` (R39); contained iff local normal death ≤ 5% (an R37-style margin).

## Result 1 — the compatible window (no immunosuppression, i*=0.26)
| b_n_base | κ=0 | 0.2 | 0.3 | 0.5 | 0.8 |
|---|---|---|---|---|---|
| 0.10 (tight gate) | [0.26,1.0] | [0.26,1.0] | [0.26,1.0] | [0.26,0.65] | [0.26,0.40] |
| 0.20 | [0.26,1.0] | [0.26,1.0] | [0.26,0.80] | [0.26,0.50] | [0.26,0.30] |

**The window stays OPEN across the whole inflammatory-spread range** — it *narrows* as the death gets more locally
inflammatory (upper bound falls 1.0 → 0.40) but **never closes**. No free contradiction in the realistic regime.

## Result 2 — collateral AT the metastasis-clearing dose (i = i* = 0.26)
Local normal-tissue death is **tiny** even for a very inflammatory death: 0.001 (κ=0) → 0.008 (κ=0.8) at the tight
gate; ≤0.024 at the looser gate. So at the *threshold* immunogenicity, being immunogenic costs almost no extra
local collateral. The requirements **compose comfortably** there.

## Result 3 — where it DOES break (the cornered failure mode)
R39's i* **rises with immunosuppression** (0.26 none → 0.39 ×0.5 → 0.64 ×0.35). Window open/closed (tight gate):
| i* (suppression) | κ=0.2 | κ=0.5 | κ=0.8 |
|---|---|---|---|
| 0.26 (none) | open | open | open |
| 0.39 (×0.5) | open | open | open |
| **0.64 (×0.35)** | open | open | **CLOSED** |

**The contradiction is real but cornered:** it appears only when a **maximally-lytic death (high κ)** meets
**immunosuppression (high i*)** — then the immunogenicity needed for systemic clearance exceeds the immunogenicity
that stays contained, and R37 & R39 collide.

## Verdict — the chain is coherent (with κ now on the residual list)
I tried to break the thesis; **it held in the realistic regime.** R39 and R37 **compose** for any non-suppressed
case and for moderate inflammatory spread, with negligible extra collateral at the threshold dose. The honest
caveat is the cornered failure (lytic death + immunosuppression). The **new, measurable design constraint** that
falls out: the immunogenic death must **prime systemic adaptive immunity while keeping local innate inflammatory
spread (κ) low** — favour a *calreticulin / limited-HMGB1 "eat-me + prime"* signal over full lytic pyroptotic
spillage. The effector must be **tuned**: neither pure tolerogenic apoptosis (fails R39) nor maximally-lytic
pyroptosis (risks R37 under suppression). So κ joins R35's leak-correlation and R37's gate as the third measurable
knob — and the chain R34→R41 is internally consistent.

## Honest residuals
- 2D lattice CA, sharp boundary; `b_n_eff = b_n_base + κ·i` is a **linear surrogate** for innate inflammatory
  coupling; `i*` imported from R39 (a reduced model); the 5% acceptable-margin is a stated judgement.
- Robust claim = the **window stays open in the realistic regime** + the **cornered failure (lytic + suppressed)**
  + the design constraint (systemic-priming, low local spread), not the absolute κ values.
- **WET residual:** the actual local inflammatory spread (κ) of a given engineered immunogenic-death mode — the new
  thing to measure, alongside R35's leak-correlation and R37's recognition gate.

*Result: `coherence.json`. Composes R37 (containment) + R39 (immunogenicity); a deliberate red-team that the thesis
passed in the realistic regime. New ledger row A11; CAPSTONE §20. Selftest-gated; verdict audited down from an
overstated "contradiction" to the honest "survives, with a cornered failure" once the data showed the window never
closes at no-suppression.*
