# CAPSTONE — the honest, end-to-end synthesis (RUNG 1–32)

*What we set out to test, what 32 in-silico rungs actually found, and where the irreducible wet-lab line sits.*
Every number below is a **prediction or simulation**, not biology. Honest negatives are kept as first-class.
This synthesises the public rungs in `README.md`; the source concept is private (`memory/`, gitignored).
**Part I (§0–6, RUNG 1–22)** mapped the three-stage chain + the anti-evolution layering. **Part II (§7–12, RUNG
23–32)** is the build-out: the chain became a concrete **two-KEY recognizer**, half of which (the MHC-free
internal key) is now built + genome-wide-safety-validated, while the de novo-binder half hit an honest wall and
was answered by a presentation/geometry screen.

---

## 0. The question

Shriya's concept, in one line: *a cancer cell should recognise itself (or a cancer neighbour) as abnormal and
trigger its own self-destruct (apoptosis), and that death should spread.* We asked: **does this chain hold up,
even on a computer, before any lab — and where exactly does it break?** We split it into the three stages she
named — **recognition → binding → killing/apoptosis** — and tested each to destruction.

---

## 1. The three-stage chain — what each stage gave

| Stage | Question | Result | Honest residual |
|---|---|---|---|
| **Recognition** | Is there a tag on cancer and on *nothing* healthy? | **Surface route DEAD** (every surface marker leaks into a vital organ — RUNG-5/15, 0/25 safe). **Mutation route WORKS** (neoantigens are tumour-exclusive — RUNG-11). High-mutation cancers carry ≥1 clean clonal handle in most patients (RUNG-16: MSI-CRC 99%, melanoma 81–100%, … PDAC 20–68%). | personalised (per-patient), not off-the-shelf |
| **Binding** | Will a T-cell actually grip the tag? | **Safety ↔ immunogenicity ALIGN** (RUNG-17): a tumour-exclusive handle is automatically high on agretopicity (the dominant recognition driver) — being safe *is* being recognisable. Boltz-2 independently confirms the mutants present on MHC (RUNG-20). | predicted propensity, not a proven TCR (MAGE-A3 warning) |
| **Apoptosis / spread** | If one cell dies, does death spread and clear the tumour without spilling? | **Death wave validated** (RUNG-13): snaps on, irreversible, spreads, stays contained. Arena: quorum/wave/ferroptosis lead (RUNG-14). | coupling/delivery = wet-lab |

**Stage verdict:** all three stages hold up in-silico, each with a stated residual. The chain is *coherent*.

---

## 2. The reality check that nearly broke it — and didn't

The immune route silently assumes the cancer keeps its **MHC "display window"** lit. We refused to assume and
measured it:

- **RUNG-18 (genetics, 6,319 tumours):** window intact **78%** · dimmed (one allele, still works) **18%** ·
  fully dark (route dies) only **~4%**.
- **RUNG-18b (expression, 50,719 real lung cancer cells):** transcriptionally dark **12.6%** vs 0.5% in an
  immune-cell control (metric validated) — **~2× the genetic number.** So genetics *under-counts* window-loss,
  but the window is still **lit in ~85–90%** of cancer cells, and the extra loss is the *reversible* (IFN /
  epigenetic) kind.

**Verdict:** the load-bearing assumption is broadly valid — not universal (a small permanently-dark core
exists), and genetics alone under-counts ~2×. We corrected our own confidence with data instead of assuming.

---

## 3. The hard problem — evolution — and the two-route answer

**RUNG-19 (escape race):** a single-target recognition-gated wave **cannot cure an established tumour.** Cure
collapses once expected resistant founders L = μ·N₀ cross ~1; a 1 cm tumour (~10⁹ cells) is always past that.
So the bare wave is necessary but **not sufficient.** This is the honest pivot — and it has exactly **two
escape routes**, each closed by a different layer:

```
                         ESCAPE ROUTE                          CLOSED BY            EVIDENCE
  Route 1: lose the neoantigen target (antigen loss)   →   multi-target / essential   RUNG-22
  Route 2: lose MHC entirely (window goes dark)        →   NK cross-kill + wave        RUNG-21
```

- **RUNG-21 (cross-kill):** a layered **T + NK + bystander wave** cures **100% across the measured escapee
  range (4–13%)**. The tumour is *trapped*: keep MHC → T-cells kill; drop MHC → NK kills ("missing-self").
  **All three layers are load-bearing** — remove NK and cure collapses 1.00→0.07 on the dark escapees.
- **RUNG-22 (multi-target):** escape collapses **exponentially** with the number of independent targets —
  K=3 → escape ~0 at clinical size; **one essential (un-losable) driver → escape-proof.** Rule: **≥3
  independent neoantigens OR ≥1 essential clonal driver.**

**The integrated logic (the "ideal layered system"):** a cell only fully escapes if it *independently* (a)
loses all K targets **or** loses MHC, **and** (b) evades NK, **and** (c) escapes the agnostic wave. Each is a
rare, independent failure → the probabilities **multiply** → escape driven toward zero. No single layer wins;
the *combination* leaves no route uncovered. This is exactly the multi-layer defence the field is converging
on — quantified, per layer, by our rungs.

---

## 4. The actionable design that falls out

> **Per patient (high-mutation cancer): target ≥3 clean clonal neoantigens — prefer essential drivers
> (KRAS/TP53/IDH1 class) — delivered as a personalised vaccine / TCR-T, PLUS an NK-engaging arm for the
> MHC-loss escapees, PLUS a resistance-agnostic bystander killer (ferroptosis/quorum). Each layer covers
> another's blind spot.**

- **Recognition:** clean clonal neoantigens (RUNG-11/16). Screen-priority handles: CTNNB1-S37F, EGFR-L858R,
  TP53-R248Q (RUNG-17); essential drivers IDH1-R132H / KRAS-G12D / TP53 (RUNG-12/20, presentation Boltz-confirmed).
- **Binding:** safety↔immunogenicity align → the safe handles are the recognisable ones.
- **Killing:** T-cell granzyme → caspase-3 = the cell's *own* apoptosis (Shriya's self-destruct, preserved).
- **Anti-escape:** multi-target + essential (Route 1) · NK + wave (Route 2).
- **Cancer is not one disease:** this works where mutation burden is high (melanoma, lung, MSI-CRC, bladder);
  low-TMB tumours (PDAC, most breast) lack the clean handles — for those, the autonomous / metabolic windows
  (TODO) or Shriya's MHC-independent self-destruct are the backup.

---

## 5. The irreducible wet-lab line (stated, never papered over)

Everything above is in-silico. What a computer *cannot* close, and a lab must:
1. **Does a real TCR exist** for each predicted handle? (RUNG-17/20 give propensity + presentation, not a receptor.)
2. **Can the trigger be delivered** into the body, and at what reach? (RUNG-13 says one seed suffices in
   principle; the delivery fraction is a wet question — TODO wave-as-injection.)
3. **Proteome-wide mimicry** — no MAGE-A3-style healthy look-alike. (RUNG-20 checked mut-vs-self only.)
4. **Real NK efficiency & exhaustion** in the immunosuppressive microenvironment (RUNG-21 models it as a parameter).

---

## 6. What this is, honestly

**Not** a cure, **not** biology. It is a **rigorous, falsifiable, public map** of where Shriya's recognition-
triggered self-destruct can work, where it can't, and what a complete therapy must layer to beat evolution —
with a quantitative bound at every stage and the wet-lab residual named at each. The headline finding:

> **The chain works in-silico at every stage; the bare recognition-gated wave is bounded by evolution; and a
> layered system (multi-target + essential drivers + NK cross-kill + agnostic wave) closes both escape routes
> — each layer covering another's blind spot.**

We found the door, checked it isn't locked, mapped every bolt, and specified the key. We have not walked
through it — that's the lab's step.

---

# PART II — RUNG 23–32: from "the chain works" to a buildable two-KEY recognizer

## 7. The pivot — expression leaks, so the MUTATION is the only tumour-exclusive signal (R23/25)
- **RUNG-23 (autonomous AND-gate):** every intracellular *expression/metabolic* program leaks into a vital
  organ (worst-donor) → expression-level recognition **cannot** gate safely. A v1 artifact (62% "proliferation"
  in post-mitotic cells) was caught as impossible and fixed (rule-5). **Conclusion: the somatic mutation is the
  ONLY clean tumour-exclusive signal** — which re-grounds the whole recognition layer and motivates sensing the
  mutation *itself*.
- **RUNG-25:** autonomous mutation-sensing is feasible at the RNA level for ~10/12 substitution types
  (AND-of-2 → ~1e-6 tumour-exclusive) but **fails for G·U-wobble (G>A) types** (KRAS-G12D, IDH1-R132H, TP53
  hotspots) → those need **DNA-level (CRISPR)** sensing.

## 8. The architecture that crystallized — TWO complementary KEYS
- **External key** — a *de novo binder* to the neoantigen pMHC (immune route; needs the MHC window lit).
- **Internal key** — an **autonomous intracellular mutation-sensing AND-gate self-destruct** (RNA toehold +
  allele-specific CRISPR), **MHC-FREE** → reaches the permanently-dark core the immune route can't (the ~4–13%
  residual Part I named). *This MHC-free self-destruct is the novel contribution.*

## 9. Internal key — BUILT + genome-wide safety-validated (R27/31)
- **RUNG-27d:** allele-specific CRISPR guides designed for **7/7 G>A-wobble drivers** (SpCas9-NGG rescues 5;
  SpCas9-NG rescues the 2 misses via deep SEED guides) → all DNA-addressable. Discrimination = PAM-creation or a
  PAM-proximal seed mismatch that collapses Cas9 on the WT allele.
- **RUNG-31/31b:** off-target scan of all 7 guides — coding transcriptome **and the full GRCh38 genome**
  (pigeonhole, exhaustive ≤1 mismatch, both strands). **0 cutting-competent (≤1mm) off-targets in any other
  gene, genome-wide.** Sanity-checked (each guide hits its own locus at 1mm = the SNV; KRAS→KRASP1 paralog is a
  positive control). → the internal key **does not fire in normal cells** at the resolution that matters. The
  MHC-free recognizer is built + safety-checked in-silico. (Residual: measured cutting / GUIDE-seq.)

## 10. External key — the honest negative, then the rigorous fix (R26/28)
- **de novo binders BIND but can't DISCRIMINATE subtle mutations.** IDH1-R132H/A\*01:01: strong dual-validated
  pMHC binders (PXDesign + Protenix) but **NULL** on His↔Arg discrimination across AF2-IG + ColabDesign-AF2 +
  Protenix **and** the principled **negative-design** method (ProteinMPNN two-state, R26f). Autopsy: positive
  design + a hotspot makes the binder *contact* the mutation, not *depend* on it. BRAF-V600E: mutation **buried**
  (2% exposed) + weakly presented (R26e). **Lesson: subtle/buried single substitutions are intractable for de
  novo binders** — the immune route is *bounded* for those (they're covered by the internal key instead).
- **The fix — pick targets by a SCREEN, not by prestige (RUNG-28):** MHCflurry × 16 drivers × 13 HLA-I →
  targets where a binder *can* win, via two mechanisms:
  - **PIK3CA-E545K / A\*03:01 (R29):** *presentation FLIP* — E545K creates the A\*03/A\*11 C-terminal anchor →
    WT (20 µM) isn't presented while MUT (51 nM) is (~397×) → **a strong MUT-pMHC binder is auto-specific, no
    discrimination needed.** [binder design in-flight; v1 under-powered, v2 with peptide hotspot + max batch running]
  - **KRAS-G12D / A\*11:01 (R30b):** *read-the-mutation* (Gly→Asp = presence-vs-absence) — staged from a **real
    crystal** (PDB 7OW6), but the crystal shows the G12D Asp only ~5% exposed (caveat); **gated on a free-pMHC
    fold** to confirm exposure before spending design quota.

## 11. Presentation validated against REAL mass-spec (R32)
- **Benign immunopeptidome (HLA Ligand Atlas):** PIK3CA WT `STRDPLSEITE` **ABSENT** (with A\*03/A\*11 donors
  present) → the presentation flip is corroborated by data, not just MHCflurry. KRAS WT `VVVGAGGVGK` **PRESENT**
  on normal lung/testis → no flip → the KRAS binder *must* discriminate.
- **IEDB (cancer/epitope):** both MUT peptides experimentally observed. **KRAS-G12D = gold-standard** (MS-eluted
  + x-ray + T-cell-positive 94/106 on A\*11:01/A\*03:01); PIK3CA-E545K = thinner (cellular-MHC binding only).
  → our predicted neoantigens are real; evidence levels now explicit.

## 12. Updated headline
> The recognition-gated self-destruct now has **two keys.** The **internal key** — an MHC-free, autonomous
> mutation-sensing CRISPR AND-gate — is **built and genome-wide-safety-validated** for the wobble drivers,
> covering the MHC-dark core the immune route cannot. The **external key** (de novo pMHC binder) taught an honest
> lesson — *subtle/buried single substitutions are undiscriminable* — which was answered by a presentation/
> geometry **screen** that selects targets where a binder can win (PIK3CA presentation-flip; KRAS read-the-
> mutation), with presentation now **corroborated by real mass-spec.** Part I showed the chain is *coherent*;
> Part II made *half of it a validated, safety-checked design* and put the other half on screen-chosen, data-
> backed targets — binder artifacts pending.

**Updated wet-lab line (adds to §5):** internal key — real allele-specific *cutting* efficiency + a clean
GUIDE-seq off-target profile (computational ≤1mm-clean ≠ measured); delivery of the circuit. External key —
binder *affinity/specificity* (SPR/cellular) + immunopeptidomic confirmation in the patient's own tumour.

## 13. Kill-coupling — recognition → self-destruct (the internal key becomes a root-kill, R34)
The recognizers *detect* the mutation; this closes the loop to **death**. Architecture: **recognizer × N → AND-gate
→ apoptosis effector** (iCasp9-class — the clinical CAR-T safety switch proves engineered apoptosis-on-demand).
Every recognizer leaks (~few % in normal cells); a single one wired to a death switch would kill ~5×10⁹ of 10¹¹
normal cells = lethal. **AND-gating two independent recognizers** makes the joint false-death ~L₁·L₂, and the
cooperative apoptosis threshold (eARM, §1) filters the rest:
- N=1 (5% leak): normal false-death **1.3×10⁻⁴** — unsafe. **N=2: 4×10⁻¹¹, therapeutic index ~10¹⁰, kill 0.99.**
- The threshold amplifies the mutant signal (kill stays ~0.99); cells that don't fully fire are cleared by the
  **bystander death wave** (§3) → tumour-level cure with cell-level partial kill.

→ **The internal key is now a complete root-kill logic**: recognize the somatic mutation / oncometabolite (the
only tumour-exclusive signal, §7) → AND-gate for specificity → trigger the cell's *own* apoptosis → wave clears
escapers. MHC-free, autonomous. **The load-bearing residual: independence of recognizer leaks** (correlated
mis-firing in stressed normal cells would collapse the AND margin — the key thing a lab must measure), plus a
non-leaky death module, circuit kinetics (a full coupled-ODE is next), and delivery.

## 14. Kill-circuit kinetics + leak correlation — stress-testing §13's headline number (R35)
We then attacked §13's own load-bearing residual instead of banking the 10¹⁰. Two hidden assumptions fell:
- **Leak is all-or-none, not graded.** §13's 4×10⁻¹¹ assumed leak is a uniform sub-threshold *level* the Hill
  switch filters. But off-target CRISPR cuts / toehold mis-triggers are **all-or-none per cell** → false-death ≈
  P(all N fully fire). Under that conservative model, **N=2 at 5% leak gives ~2.5×10⁻³ (ρ=0), not 4×10⁻¹¹** — and
  rises to ~5×10⁻² as leaks become **correlated** (one-factor Gaussian copula, ρ→1). *Correlation, not N, is the
  binding constraint*: at ρ≥0.6 no feasible N reaches 10⁻⁴.
- **Kinetics supply a second, orthogonal filter.** A reduced eARM bistable caspase ODE (irreversible MOMP) shows a
  supra-threshold mis-fire must be **sustained ≥ T_min (~0.7/k_deg, ~hour-scale)** to commit — a *temporal* AND on
  top of the molecular AND. If leaks are **transient and independent**, effective false-death drops several more
  orders (transient + N=3 + ρ=0 → ~10⁻⁵).

→ **Honest correction:** the kill-coupling *logic* holds, but its safety is **conditional**, not a free 10¹⁰. The
design target is now explicit and testable: per-recognizer leak amplitude L, leak **duration ≪ commitment time**,
leak **correlation ρ→0**, and **N≥3** if ρ can't be driven low. This converts §13's hand-wave residual into four
measurable engineering specs — the #1 wet measurement being recognizer-leak **correlation in stressed normal cells**.

## 15. Delivery × bystander wave — separating the kill from delivery (R36)
R35's safety rests on an escape hatch (*"partial single-cell kill is fine because the bystander wave clears the
rest"*) that also collides with the residual named on every rung: **delivery reaches only some cells.** R36 models
both as one **percolation** question — circuit delivered to a seed fraction f, fires w.p. p_kill, dying cells
recruit neighbours w.p. b (the resistance-agnostic wave) ≡ bond-percolation ignition (2D p_c=0.5):
- **Super-critical wave (b>p_c):** a **5% seed clears ~99%**; minimum delivery f\* collapses ~3 orders (→0.2%).
  *Delivery stops being the bottleneck.*
- **Sub-critical wave:** the wave dies locally → near-complete delivery required. *Delivery IS the bottleneck.*
- **Partial kill** (p_kill=0.5) clears the tumour **only** when super-critical — the precise condition that makes
  §14's partial single-cell kill acceptable at the tumour level.
- **Hard limit:** wave-resistant escapees (R18/R21) cap clearance at **~1−r** → a resistance-agnostic 2nd killer
  (NK, §3/R21) is mandatory for the last fraction. **3D geometry lowers p_c (~0.25) → even more forgiving.**

→ The therapy's hardest practical objection ("you can't deliver to every cancer cell") is **formally separable**
from the kill mechanism: the engineering target becomes *make the bystander coupling super-critical* + *add an
agnostic killer for the resistant fraction* — two individually-tractable problems instead of one impossible one.

## 16. Wave containment — the danger §15 hid (R37)
§15 says the wave must be **super-critical** to clear the tumour. But the wave is **resistance-agnostic** and the
tumour is a small island in a **sea of normal tissue** — so we asked the question §15's pure-tumour lattice hid:
does the super-critical wave run away through normal tissue? Two-tissue lattice (tumour disk in a normal sea),
wave bond-prob set by the *receiving* cell's type:
- **An ungated super-critical wave is a catastrophe by boundary spillover *alone*** — it crosses the tumour edge
  and kills ~99.8% of normal tissue **even at zero leak.** The wave doesn't know the boundary.
- **Containment requires a recognition-gated bystander signal:** super-critical into tumour (b_t), **sub-critical,
  well below p_c, into normal** (b_n) — the death completes only where the receiver's own mutation-sensor is primed.
  Then tumour clearance is unchanged while normal death stays small.
- **The wave amplifies the §14 leak** by the normal-tissue death-cluster size (~tens× even when contained; →∞ if
  ungated) ⇒ effective normal false-death = (leak)×(cluster size); §14's leak target must be *tightened* by that
  factor. An intrinsic collateral **rim** (~1–6%) dies at the boundary regardless.

→ A **second load-bearing residual** joins §14's leak-correlation: the bystander signal must be **recognition-gated**
(b_n ≪ b_t, well below p_c). The very super-criticality that lets the wave cure from a tiny seed makes it lethal in
normal tissue unless gated — and because the wave *multiplies* the leak, the gate is as load-bearing as the leak.

## 17. The end-to-end safety envelope — composing it all into one number (R38)
§14 (molecular leak) and §16 (wave containment) were treated separately; R38 composes them and finds the safety is
a **product of three factors**:
> **tissue normal false-death = (molecular leak, §14) × (wave amplification, §16) × (delivery footprint)**
- **The two residuals MULTIPLY** — the molecular AND-gate leak is scaled by the wave's normal-tissue cluster size
  (3× at a tight gate, 16× loose, 210× near p_c). End-to-end death spans **~7,600×** across the design grid.
- **Lever ranking:** adding an AND-input (N 3→2: **60×**) > recognizer correlation (ρ: **42×**) > leak transience
  (**27×**) > gate tightness (5×, →210× near p_c). No single factor suffices; you must win all of them.
- **Delivery localisation is *also* a safety lever:** the leak only hits circuit-*carrying* normal cells, so
  tumour-localised delivery (§15: a few-% seed already cures) shrinks the at-risk population by orders — §15's
  efficacy result **is** the safety result.

→ The autonomous root-kill is **tissue-safe *and* curative in a concrete, composed region**: N≥3 + transient,
uncorrelated leaks + a tight recognition-gate + localised delivery — and in that region the same wave still clears
the tumour from a partial seed. "Safe enough" is no longer a hand-wave but the **product of four measurable
molecular numbers (L, ρ, leak-lifetime, b_n) × the delivery footprint** — the complete in-silico safety envelope of
the internal key, composing §13–16 into one verdict.

## 18. Reaching the metastases — local cure vs whole-body cure (R39, Shriya's H7)
Everything above clears the **local** tumour. But the bystander wave is contact-based — it cannot reach a distant
micrometastasis that never received the circuit, and **metastasis is what kills patients.** The only mechanism
that reaches untreated distant deposits is systemic immunity raised by **immunogenic cell death (ICD)** — the
abscopal effect. R39 tests it, and finds three regimes:
- **Tolerogenic (clean) apoptosis → local-only cure.** Below the immune-priming threshold, no systemic immunity is
  raised; distant metastases **grow back.** The primary is cured but the patient is not. (Plain apoptosis is often
  *tolerogenic* — a real danger here.)
- **Immunogenic death (above threshold) → abscopal cure.** ICD (calreticulin/HMGB1/ATP, or a pyroptotic/
  necroptotic flavour) raises a persistent T-cell response that eradicates untreated deposits (i\*≈0.26, rising
  with immunosuppression).
- **Heavy immunosuppression → immunity alone fails** even at maximal immunogenicity → exactly where §3/R21's
  resistance-agnostic NK arm and checkpoint-release are required.

→ **New design requirement:** the internal self-destruct must be **engineered for immunogenic cell death, not clean
apoptosis** — that is what turns Shriya's local *"destroy itself from within"* into a **whole-body cure**, adding a
third clearance layer (local wave §15 + systemic ICD-immunity §18 + NK tail §3). H7 is not a footnote; it is the
local-vs-systemic pivot of the whole concept.

## 19. Where the independent recognizers come from — orthogonal modalities (R40, Shriya's H5)
§14/§17 demanded **N≥3 independent recognizers** with low leak-correlation, but the recognition we built is all
*one modality* (mutation sensors read the same somatic-signal class → likely **correlated** leaks). R40 resolves
this with a block-correlated leak model and a clean principle: **independence, not individual cleanliness, is the
currency of an AND-gate.** A 30%-leaky *independent* microenvironment input beats a 5%-clean *correlated* mutation
input; the winning gate ANDs **orthogonal physical modalities** — mutation (DNA/RNA) + oncometabolite (2-HG, §11's
H6) + microenvironment (hypoxia/acidosis, H5) — whose leaks multiply (joint ≈ 3×10⁻⁴), beating a correlated
mutation stack by 16–83× as correlation rises. This gives Shriya's **microenvironment sniffer (H5)** a concrete
role — *not* a standalone recognizer (it leaks), but the gate's **independence-provider** — and converts §14's #1
residual (leak correlation) into a **design principle: choose recognizers that mis-fire for unrelated reasons.**
*Honest conditionality:* the advantage holds only if same-modality leaks are genuinely correlated; if not, clean
mutation-stacking wins — so the decision is gated on the same wet measurement (ρ) §14 already flagged.

## 20. Red-team: do the requirements contradict? (R41)
§18's immunogenic-death demand and §16's containment demand pull in opposite directions — immunogenic death is
*inflammatory* (DAMP/lytic spread), and inflammation doesn't respect the recognition gate. So we deliberately tried
to **break the thesis**: model the effective normal coupling as `b_n_eff = b_n_base + κ·i` (κ = the death mode's
local inflammatory spread) and ask whether any immunogenicity both clears metastases (i ≥ i*) and stays contained.
**The thesis held in the realistic regime:** at no immunosuppression the compatible window stays open across all κ
(it narrows but never closes; local collateral ≤2.4% at the threshold dose). The contradiction is **cornered** — it
appears only when a *maximally-lytic* death (high κ) meets *immunosuppression* (which raises i* to 0.4–0.64). →
The chain §13–§19 is **internally coherent**, and a new measurable constraint falls out: the immunogenic death must
**prime systemic adaptive immunity while keeping local innate inflammatory spread κ low** (calreticulin/limited-HMGB1
"eat-me + prime", not full lytic pyroptosis). The effector is *tuned*: neither tolerogenic apoptosis (fails §18) nor
maximally-lytic pyroptosis (risks §16 under suppression). κ joins leak-correlation, the gate, and the four envelope
numbers on the list a lab must keep bounded.

---

*Rungs: see `README.md` (hypothesis catalog), `STATUS.md` (live map), `docs/HYPOTHESIS_LEDGER.md` (every
hypothesis → verdict). Next experiments: `TODO.md`. Source concept & strategy: private (`memory/`, gitignored).*
