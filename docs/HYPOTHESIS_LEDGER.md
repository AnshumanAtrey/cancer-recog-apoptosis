# HYPOTHESIS LEDGER — canonical tracker (every hypothesis → test → verdict)

**This is the single place every hypothesis lives and gets crossed off.** After each test (rung), update the
row here. Per-rung detail lives in `runs/<rung>/README.md`; the forward queue in `TODO.md`; the narrative
synthesis in `CAPSTONE.md`; the live one-page state in `STATUS.md`.

**Legend:** ✅ CONFIRMED (in-silico, with residual) · ❌ REFUTED/CLOSED · 🔄 IN-FLIGHT · ⏳ PARTIAL · ⬜ OPEN (untested)
*All ✅ are in-silico predictions with a stated wet-lab residual — never biology proven.*

---

## A. Core chain (the core concept: recognition → binding → killing) — from CAPSTONE
| # | Hypothesis | Verdict | Where |
|---|---|---|---|
| A1 | A surface marker tags cancer and nothing healthy | ❌ REFUTED — every surface marker leaks into a vital organ (0/25 safe) | R5/R15 |
| A2 | The somatic **mutation** is a tumour-exclusive tag | ✅ CONFIRMED — neoantigens tumour-exclusive; ≥1 clean clonal handle in most high-TMB patients | R11/16/17 |
| A3 | The MHC "display window" stays lit on cancer | ✅ CONFIRMED — lit in ~85–90% (genetic 78%/dim 18%/dark 4%; expression dark ~12.6%) | R18/18b |
| A4 | A binder/TCR can grip the tag; safe handles are recognizable | ✅ CONFIRMED (propensity) — safety↔immunogenicity align; mutants present on MHC | R17/R20 |
| A5 | One cell's death spreads & clears the tumour, contained | ✅ CONFIRMED (sim) — death wave snaps on, irreversible, spreads, stays contained | R13/R14 |
| A6 | Expression-LEVEL differences can gate safely | ❌ REFUTED — every expression window leaks → mutation is the ONLY tumour-exclusive signal | R23 |
| A7 | Clinical tumours need a 2nd, resistance-agnostic mechanism | ✅ CONFIRMED — cure collapses at L≈1; bystander shifts curable size ~10× | R19 |

## B. Recognition-window handles (the 12 candidate "tells") — from TODO.md
| # | Handle | Verdict | Where |
|---|---|---|---|
| B1 | Mutations (neoantigens) | ✅ main handle | R11/16/17 |
| B2 | MHC window | ✅ ~85–90% on | R18/18b |
| B3 | Surface markers (HER2/EGFR/EpCAM) | ❌ leaks | R15 |
| B4 | Glycans (Tn/sialyl-Tn) | ⬜ open (leak-test candidate) | — |
| B5 | Warburg / glucose-acid | ⬜ open | — |
| B6 | Apoptosis machinery (BCL-2 high) | ⏳ listed (BH3 mimetics) | — |
| B7 | CD47 "don't eat me" | ⬜ open (leak-test candidate) | — |
| B8 | Aneuploidy (DNA amount) | ⬜ open | — |
| B9 | Telomerase | ⬜ open | — |
| B10 | Division rate (Ki-67) | ⬜ open | — |
| B11 | PS lipid flip | ⬜ open (leak-test candidate) | — |
| B12 | Synthetic-lethal dependency (MTAP loss) | ✅ tested | R14 |

## C. The two recognition KEYS (current frontier, rungs 26–31)
| # | Hypothesis | Verdict | Where |
|---|---|---|---|
| C1 | **Internal key** — intracellular mutation-sensing AND-gate self-destruct (RNA toehold + allele-specific CRISPR) | ✅ CONFIRMED — 7/7 wobble drivers DNA-addressable | R27b/c/d |
| C2 | Internal-key guides don't fire in normal cells (no off-target) | ✅ CONFIRMED — **0 cutting-competent (≤1mm) off-targets GENOME-WIDE** (all 7 guides, full GRCh38) | R31/31b |
| C3 | **External key** — de novo pMHC binder can discriminate IDH1-R132H (His↔Arg) | ❌ REFUTED — binds but can't discriminate (NULL ×4: AF2-IG, ColabDesign, Protenix, + negative design) | R26c/d/f |
| C4 | External key works for BRAF-V600E | ❌ REFUTED — mutation buried (2% exposed) + weakly presented (3971 nM) | R26e |
| C5 | Pick binder targets by a presentation/exposure SCREEN, not prestige | ✅ CONFIRMED — screen found 3 presentation-flip + 34 up-facing read-the-mutation targets | R28 |
| C6 | **PIK3CA-E545K/A\*03:01** — presentation flip (WT unpresented) → binder auto-specific | ⏳ PARTIAL/UNCONFIRMED — v1 0/10; **v2 (hotspot B4+B6, Extended) = 1/10 Protenix-basic-only, 0/10 AF2-IG, 0 dual-passers**. rank_1 ptx_iptm 0.87 vs af2_iptm 0.14 (oracles disagree → single-model artifact). Difficulty gauge: ≤5% hardest tier on BOTH oracles. Flip is real (C8) but the 11-mer binder is at the edge of de-novo tooling at small batch. NOT a breakthrough binder. Next: larger batch / KRAS companion (C7) | R29v2 |
| C7 | **KRAS-G12D/A\*11:01** — read-the-mutation (Gly→Asp presence/absence) | 🔄 IN-FLIGHT — **GO**: free-pMHC fold shows p6 Asp **30% exposed (up-facing)**; crystal's 5% was a TCR-artifact. Gold-standard neoantigen (R32b). Design target staged (hotspot B6); PXDesign after PIK3CA queue. WT presented (R32) → MUT-vs-WT scoring = make-or-break | R30c/R32 |
| C8 | Target presentation is real (not just MHCflurry-predicted) | ✅ CONFIRMED (real MS) — benign: PIK3CA WT ABSENT (flip corroborated), KRAS WT PRESENT (no flip). MUT-in-tumour (IEDB): KRAS-G12D = gold-standard (MS-eluted+x-ray+T-cell+ on A\*11/03); PIK3CA-E545K = binding-only (thinner) | R32/32b |
| C9 | **Kill-coupling** — recognizers → AND-gate → apoptosis kills mutant, spares normal | ✅ FEASIBLE (logic), ⚠️ headline TI ~10¹⁰ **CORRECTED by R35** — that number assumed leak is a uniform sub-threshold *level* the Hill switch filters. Internal key = detector → root-kill; the *logic* holds, but the safety MARGIN is far smaller and conditional (see C10). Single recognizer (1.3×10⁻⁴) unsafe regardless | R34→R35 |
| C10 | **Kill-circuit kinetics + leak correlation** — does the AND-gate survive realistic (all-or-none, correlated) leak + irreversible kinetics? | ⏳ CONDITIONAL — under all-or-none leak, N=2 gives false-death **~2.5×10⁻³** (ρ=0) → ~5×10⁻² (ρ→1), NOT 4×10⁻¹¹. Recoverable via 4 explicit levers: leak amplitude L, leak **duration**≪commit-time (bistable-MOMP temporal filter, T_min~0.73/k_deg), **correlation ρ→0**, **N≥3**. Transient+N=3+ρ=0 → ~10⁻⁵. **Correlation ρ (not N) is the binding constraint** — the #1 wet measurement | R35 |

## D. The 15 conceptual hypotheses (`docs/hypothesis-to-consider.txt`) — cross-referenced & crossed off
| # | Conceptual hypothesis | Verdict vs our work | Where / next |
|---|---|---|---|
| H1 | Metal Overload (ferroptosis/cuproptosis) | ⏳ PARTIAL — ferroptosis tested as an alt-death mechanism in the arena | R14; metal-ligand delivery OPEN |
| H2 | Sodium Rush (NECSO/TRPM4) | ⬜ OPEN — untested (2025 pathway) | future |
| H3 | Mitochondrial Sabotage (peptide nanostructures) | ⬜ OPEN | future |
| H4 | **AND-Gate logic** (two markers) | ✅ DONE — this is a CORE spec, extensively tested | R10/R27b |
| H5 | Microenvironment Sniffer (pH/ROS/hypoxia) | ⬜ OPEN — testable now (conformational sensor) | **24h-window candidate** |
| H6 | Metabolic Waste / oncometabolite (2-HG, succinate) | ✅ FEASIBILITY CONFIRMED (R33) — 2-HG recognition gate for IDH1-R132H is feasible & EASIER than the binder route (the ~2000× concentration differential does the discrimination, not binding selectivity; needs a 2-HG-specific/non-agonist sensor). **RECOVERS IDH1** + GROUNDED: the D-2-HG-specific recognizer **DhdR** already exists (Kd 1.16µM, rejects α-KG+28 metabolites, biosensor 0.3-30mM); design = detune to ~0.1-1mM + couple to death | R33/33b |
| H7 | Danger Signal (ICD/DAMPs) | ⬜ OPEN — not yet scored; an ICD score on the kill mechanism is cheap | **24h-window candidate** |
| H8 | Gap-Junction Wave (connexin) | ✅ DONE — connexin relay/bridge tested | R12p_connexin |
| H9 | Fas-FasL Bystander | ⏳ PARTIAL — bystander cross-kill explored | R13/R21 |
| H10 | Evolutionary Trap (dual selective pressure) | ⏳ PARTIAL — escape race tested; the specific "trap" variant OPEN | R19 |
| H11 | Synthetic-Lethal Partner | ✅ DONE — MTAP synthetic-lethal tested | R14 |
| H12 | Living Drug Factory (engineered bacteria) | ⬜ OPEN (delivery/FUT) | future |
| H13 | Self-Amplifying mRNA | ⬜ OPEN (delivery/FUT) — relates to internal-key delivery | R24-adjacent |
| H14 | Nano-Origami Death Star | ⬜ OPEN (FUT) | future |
| H15 | Microwave-Activated Biorobot | ⬜ OPEN (FUT) | future |

**Reading D:** 3 of your 15 are already DONE (H4, H8, H11), 3 PARTIAL (H1, H9, H10), 9 OPEN. The OPEN-and-
testable-now set = **H5 (microenvironment), H6 (oncometabolite — strong IDH1/2-HG tie-in), H7 (ICD score)** —
these are the real 24h-window candidates (the rest are wet-lab/delivery/future). H2/3/12/14/15 need new
pathways or physical platforms (future-safe, not now).

---

## Maintenance rule
When a rung finishes: (1) set its row's verdict here (cross it off), (2) add the rung pointer, (3) if it was a
conceptual hypothesis (Section D), update its cross-ref. New hypotheses get a row here FIRST (status ⬜ OPEN),
then a TODO.md entry. `CAPSTONE.md` is the narrative (**refreshed to R32, 2026-06-19 — Part II = the two-key
build-out**); `THESIS.md` carries a SUPERSEDED banner (old DR5-bispecific thesis kept as history; current thesis
= CAPSTONE §8–12).
