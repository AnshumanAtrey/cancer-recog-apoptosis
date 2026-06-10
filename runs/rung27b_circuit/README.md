# RUNG-27b — The autonomous mutation-sensing circuit (Shriya's MHC-free self-destruct, made buildable)

**What this is.** The counterpart to RUNG-11/16's *immune*-route map, for the route the immune system can't
reach. RUNG-23 proved every expression window leaks → the mutation is the only tumour-exclusive signal;
RUNG-25 proved single-base RNA sensing is feasible but substitution-dependent. RUNG-27b **forges the actual
key**: designs allele-specific sensors for the real recurrent clonal driver hotspots from **verified CDS**,
AND-gates two co-occurring drivers, and maps per-cancer addressability. Unlike the de novo *binder* (RUNG-26,
0 designs crossed), this route **produced designed sensors** — it's the more tractable key.

Selftest 16/16. Every hotspot codon **verified against real Ensembl MANE CDS** (`prep`). ΔΔG by ViennaRNA.

## Headline finding — the commonest drivers need DNA-level sensing

| | n | drivers |
|---|---|---|
| **RNA-toehold-sensable** (non-wobble, ΔΔG ≥ 2) | 9/16 | KRAS-G12V/G12C, NRAS-Q61R/Q61K, BRAF-V600E, IDH1-R132C, PIK3CA-H1047R, EGFR-L858R, TP53-R248W |
| **G·U-WOBBLE → DNA/CRISPR only** | 7/16 | **KRAS-G12D, KRAS-G13D, IDH1-R132H, PIK3CA-E545K, TP53-R175H/R248Q/R273H** |

**ViennaRNA independently confirms the wobble penalty** (not just our flag): median ΔΔG **0.40** kcal/mol for
wobble (false-fire 30–52% — undiscriminable; PIK3CA-E545K is even negative) vs **4.40** for non-wobble
(false-fire 1e-3…1e-5). The wobble set is exactly the **G>A transition** drivers (CpG-deamination / aging
signature) — the *most common* driver mutations. **An RNA toehold cannot discriminate them; DNA-level CRISPR
allele sensing (no wobble) can.** The autonomous gate is buildable; the *modality* is dictated by chemistry.

## Per-cancer addressability — off-the-shelf 2-driver AND-gate (MHC-free)

P(≥2 co-occurring clonal sensable drivers from different pathway groups):

| cancer | RNA-only | RNA + DNA | DNA rescue |
|---|---|---|---|
| CRC | 0.01 | **0.21** | +0.20 |
| Glioma | ~0 | **0.20** | +0.20 |
| PDAC | ~0 | **0.15** | +0.15 |
| LUAD | ~0 | 0.09 | +0.09 |
| Melanoma | 0 | 0.05 | +0.05 |
| Breast | 0 | 0.04 | +0.04 |
| Thyroid | 0 | 0 | — (single driver — no 2-input AND) |

**RNA-only is ≈0** because the top oncogene + TP53 hotspots are G>A wobble; **DNA-level sensing rescues it
entirely** (the gain *is* the coverage). This is the **off-the-shelf shared-driver FLOOR** — pre-designed
sensors, 2 pathway groups, independence-modelled, **conservative**. A personalised 2nd input (any clonal
mutation, RUNG-16-style) lifts it; the off-the-shelf advantage is that the sensors are *shared across
patients*, not bespoke.

## Why it matters (the unique value)
This route fires **inside the cell on the mutation itself, no MHC** → it reaches exactly the **~4–13%
MHC-dark escapees** (RUNG-18) that the immune route is blind to, and the AND of two somatic mutations is
**tumour-exclusive by construction** (normal cells have neither; mutual-exclusivity forces the 2 inputs from
different pathways). It's the backup layer the capstone's "ideal system" needs for the dark core.

## Ceiling (rule 3/5)
- ViennaRNA ΔG is a thermodynamic **proxy** (kinetics / RNA accessibility / genome off-targets not modelled).
- Per-cancer driver frequencies are TCGA/COSMIC **approximations** → coverage is order-of-magnitude.
- Cross-group co-occurrence modelled **independent** (TP53 co-occurs more → conservative); within a
  mutual-exclusivity group only one driver is ANDable.
- "DNA-sensable" = CRISPR allele discrimination is **established** (buildability), not a built circuit.
  Synthesis + delivery + in-cell apoptosis actuation = the wet-lab residual.

*Result: `rung27b_circuit.json` (+ designed sensor sequences per driver). Figure: `rung27b_circuit.png`.
Script + selftest: `scripts/54_mutation_circuit.py`.*
