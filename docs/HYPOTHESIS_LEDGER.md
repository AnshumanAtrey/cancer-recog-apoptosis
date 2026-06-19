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
| C6 | **PIK3CA-E545K/A\*03:01** — presentation flip (WT unpresented) → binder auto-specific | 🔄 IN-FLIGHT — folds landed; **design v1 0/10** (≤5% tier, under-powered + no-hotspot) → **v2: hotspot B4+B6 (up-facing) + max batch** | R29 |
| C7 | **KRAS-G12D/A\*11:01** — read-the-mutation (Gly→Asp presence/absence) | 🔄 IN-FLIGHT — real crystal target (7OW6, A\*11:01); **CAVEAT: p6 Asp only ~5% exposed in crystal (not cleanly up-facing); WT IS presented (R32) → must discriminate. GATED on free-pMHC fold to resolve Asp exposure before quota** | R30b/R32 |
| C8 | Target presentation is real (not just MHCflurry-predicted) | ✅ CONFIRMED (real MS) — benign: PIK3CA WT ABSENT (flip corroborated), KRAS WT PRESENT (no flip). MUT-in-tumour (IEDB): KRAS-G12D = gold-standard (MS-eluted+x-ray+T-cell+ on A\*11/03); PIK3CA-E545K = binding-only (thinner) | R32/32b |

## D. The 15 conceptual hypotheses (`docs/hypothesis-to-consider.txt`) — cross-referenced & crossed off
| # | Conceptual hypothesis | Verdict vs our work | Where / next |
|---|---|---|---|
| H1 | Metal Overload (ferroptosis/cuproptosis) | ⏳ PARTIAL — ferroptosis tested as an alt-death mechanism in the arena | R14; metal-ligand delivery OPEN |
| H2 | Sodium Rush (NECSO/TRPM4) | ⬜ OPEN — untested (2025 pathway) | future |
| H3 | Mitochondrial Sabotage (peptide nanostructures) | ⬜ OPEN | future |
| H4 | **AND-Gate logic** (two markers) | ✅ DONE — this is a CORE spec, extensively tested | R10/R27b |
| H5 | Microenvironment Sniffer (pH/ROS/hypoxia) | ⬜ OPEN — testable now (conformational sensor) | **24h-window candidate** |
| H6 | Metabolic Waste / oncometabolite (2-HG, succinate) | ⬜ OPEN — **note: IDH1-R132H (our target) MAKES 2-HG** → a natural tie-in | **24h-window candidate** |
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
then a TODO.md entry. `CAPSTONE.md` is the narrative (currently stale at R22 — needs R23–31 added);
`THESIS.md` is stale (old DR5-bispecific thesis, pre-dates the neoantigen dual-key pivot) — both flagged for update.
