# RUNG-40 — INDEPENDENCE IS THE CURRENCY: the microenvironment sniffer (Shriya's H5) as the AND-gate's independent input

R35 showed the AND-gate's normal-cell safety is bottlenecked by leak **correlation ρ**; R38 showed the safe region
needs **N≥3 independent recognizers**. The open question: *where do independent recognizers come from?* The
mutation-sensors (allele-CRISPR R27, RNA toehold R25) all read the **same signal class** (somatic nucleic-acid
changes) → their leaks are plausibly **correlated** (shared delivery, transcriptional bursting, cell-state stress).
So you can't just stack three mutation sensors and assume ρ=0. This tests Shriya's **H5** (microenvironment
sniffer) in its real role: the AND-gate's **independence-provider**. `scripts/74`.

## Model
Block-correlated Gaussian copula → multivariate-normal orthant probability (chunked Monte-Carlo, selftest-validated
against the independent product and the R35 one-factor limit). Per-sensor leak L, within-modality correlation
ρ_within, cross-modality correlation ρ_cross. Representative leaks: mutation 5%, 2-HG 2%, microenvironment **30%**
(leaky — normal hypoxic/inflamed/ischaemic tissue shares the signal).

## Result 1 — add a 3rd input to a 2-mutation gate (base joint false-death 1.2×10⁻², ρ_within=0.5)
| 3rd input | joint false-death | × base |
|---|---|---|
| 3rd MUTATION sensor, **correlated** (ρ=0.5, L=0.05) | 4.9×10⁻³ | 0.40× |
| 3rd MUTATION sensor, hypothetically independent (L=0.05) | 6.3×10⁻⁴ | 0.05× |
| **MICROENVIRONMENT sensor, independent (L=0.30, leaky!)** | **3.6×10⁻³** | 0.30× |
| 2-HG metabolite sensor, independent (L=0.02) | 2.5×10⁻⁴ | 0.02× |

**A 30%-leaky *independent* microenvironment input (3.6×10⁻³) beats a 5%-clean *correlated* mutation input
(4.9×10⁻³).** An independent input multiplies its **full** leak factor into the joint; a correlated one barely
moves it. Independence > individual cleanliness.

## Result 2 — whole 3-input gate: orthogonal modalities vs same-modality stack
| design | joint false-death |
|---|---|
| 3× mutation, correlated ρ=0.5 | 4.9×10⁻³ |
| 3× mutation, correlated ρ=0.7 | 1.1×10⁻² |
| **mutation + 2-HG + microenvironment (orthogonal, independent)** | **3.0×10⁻⁴** |

The orthogonal-modality gate reaches the near-full product of leaks (≈ 0.05·0.02·0.30 = 3×10⁻⁴) — **16× better
than 3 correlated mutation sensors at ρ=0.5, 38× at ρ=0.7.**

## Result 3 — the advantage grows with same-modality correlation (and the honest reversal)
| ρ_within (mutation block) | 3× mutation joint | orthogonal joint | advantage |
|---|---|---|---|
| **0.0** | **1.3×10⁻⁴** | 3.0×10⁻⁴ | **0.4× (orthogonal LOSES)** |
| 0.3 | 1.7×10⁻³ | 3.0×10⁻⁴ | 5.7× |
| 0.5 | 4.9×10⁻³ | 3.0×10⁻⁴ | 16× |
| 0.7 | 1.1×10⁻² | 3.0×10⁻⁴ | 38× |
| 0.9 | 2.5×10⁻² | 3.0×10⁻⁴ | 83× |

**Honest reversal:** if mutation-sensor leaks were *truly* independent (ρ_within=0), three clean mutation sensors
(1.3×10⁻⁴) would **beat** the orthogonal design (which carries the 30%-leaky microenvironment sensor). So the
orthogonal-modality advantage is **conditional on same-modality leaks actually being correlated** — which is
exactly R35's #1 unmeasured residual. The design decision (stack mutations vs go orthogonal) is therefore **gated
on measuring ρ_within**.

## Verdict
**Independence, not individual cleanliness, is the currency of an AND-gate — and that gives Shriya's H5 a concrete
role.** The microenvironment sniffer is *not* a standalone recognizer (it leaks ~30%); it is the gate's
**independence-provider**. The winning design is an AND of **orthogonal modalities** — mutation (DNA/RNA, R27/25) +
oncometabolite (2-HG, R33) + microenvironment (hypoxia/acidosis, H5) — that fail for *unrelated physical reasons*,
so their leaks multiply. This **answers R38's open question** (the N≥3 independent recognizers come from orthogonal
physical modalities, not more of the same) and **turns R35's #1 residual into a design principle**: choose
recognizers that mis-fire independently. The catch keeps it honest: the whole advantage is conditional on the
mutation-sensor correlation R35 told us to measure.

## Honest residuals
- Leaks (0.05 / 0.02 / 0.30) and correlations are **representative, not measured**: the microenvironment standalone
  leak is literature-typical for hypoxia markers (CA9/HIF targets appear in normal hypoxic tissue) and is swept;
  the mutation-sensor **cross-correlation is the R35 wet residual** (must be measured in stressed normal cells).
- Gaussian-copula leak model is a modelling choice; robust claim = the **ranking** (independent > correlated) +
  the **orthogonal-modality design principle** + its **conditionality on ρ_within**, not absolute fractions.
- The microenvironment sensor itself (a conformational pH/hypoxia switch coupled to the kill circuit) is an
  **un-built design** (cf. the 2-HG sensor, R33).

*Result: `independence.json`. Tests Shriya's H5; resolves R38's "where do independent recognizers come from"; makes
R35's leak-correlation a design lever. Ledger H5 → DONE; CAPSTONE §19. Selftest-gated (copula orthant validated
against the independent product + correlation-increases-leak).*
