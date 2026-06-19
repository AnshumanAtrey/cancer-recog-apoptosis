# Project STATUS — recognition-gated cancer self-apoptosis (one-page map)

**Goal:** a system that fires (kills) ONLY in cancer cells, never normal — discrimination by the cell's
own driver *mutation*. Two "keys" that recognize the mutation, then trigger death.

## The core safety logic (normal vs cancer) — TESTED
- **Expression/antigen LEVEL leaks** into normal cells (RUNG-23) → we do NOT rely on "more antigen."
- **The somatic MUTATION is the only tumour-exclusive signal** → both keys gate on the mutation, so the
  system is "tumour-exclusive by construction" (normal cells lack the mutation). (RUNG-27b)
- Per-cell gate safety quantified (RUNG-7 → grounded by RUNG-8 HLA data → RUNG-12 q_n ≤ 0.02 per handle).

## Pipeline status
| Layer | What | Status |
|---|---|---|
| **Recognition BASIS** (what distinguishes cancer) | driver mutations are addressable + tumour-exclusive | ✅ TESTED (R5/6/7/8/11/12/23/27a) |
| **Internal key** (intracellular mutation-sensing AND-gate self-destruct) | RNA toehold + allele-specific CRISPR | ✅ DONE — **7/7 wobble drivers DNA-addressable** (R27b/c/d) |
| **External key** (de novo pMHC binder = immune-route recognizer) | design a binder that's mutation-specific | 🔄 IN-FLIGHT (see below) |
| **Effector** (the killing, once recognized) | apoptosis / bystander wave / autonomous circuit / delivery | ⏳ modeled (R13/21/23/24), not the current focus |

## External key — where we are RIGHT NOW
- ❌ **IDH1-R132H** CLOSED — His↔Arg too subtle; binder binds but can't discriminate (NULL ×4 incl. negative design). R26c/d/f.
- ❌ **BRAF-V600E** CLOSED — mutation buried (2% exposed) + weakly presented. R26e.
- ✅ **Built the target-selection SCREEN we'd skipped** (R28): 16 drivers × 13 HLA → picks targets by presentation + geometry.
- 🔄 **PIK3CA-E545K / A\*03:01** (R29) — *presentation flip*. Folds landed; binder design v1 **0/10** (≤5% tier, under-powered+no-hotspot) → **v2: hotspot B4+B6 + max batch**.
- 🔄 **KRAS-G12D / A\*11:01** (R30c) — *read-the-mutation* **GO**: free-pMHC fold shows p6 Asp **30% exposed** (crystal's 5% was a TCR-artifact). Design target staged (hotspot B6); PXDesign after the PIK3CA queue.
- ✅ **Immunopeptidomics check (R32)** — real mass-spec (HLA Ligand Atlas) corroborates: PIK3CA WT **absent** from normal (flip real); KRAS WT **present** on normal lung/testis (binder must read G12D). MUT-in-tumour = next layer.
- → Two complementary shots; **currently fold-bound (hours).**

## Honest residuals (named, not hidden)
- Internal key: ON-target validated (R27c/d); off-target DONE — **0 cutting-competent (≤1mm) off-targets GENOME-WIDE** across all 7 guides (R31 coding-transcriptome + R31b full GRCh38, sanity-checked). Residual: mm2+ intergenic (low-consequence) + measured cutting (GUIDE-seq).
- External key: in-silico fold/affinity ≠ measured; immunopeptidomic presentation + binder SPR/cellular = wet-lab residual.
- Internal-key sensor MOLECULES: thermodynamic feasibility only (R27b); actual molecular design (ODesign) pending.

## Next experiments (no folds needed)
1. **CRISPR guide off-target scan** (the internal key's named residual — does any allele-specific guide cut a normal-cell locus?).
2. **ODesign** the internal-key sensor molecules (feasibility → real design).
3. When folds land: pre-crop + exposure-audit + PXDesign (PIK3CA no-hotspot, KRAS hotspot-on-p6).
