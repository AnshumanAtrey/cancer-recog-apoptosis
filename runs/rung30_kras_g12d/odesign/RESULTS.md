# ODesign KRAS-G12D / A\*11:01 binder run — RESULTS (2026-06-27)

First successful **external-key** generation run. PXDesign was bounded at **0/80 dual-oracle passers** on this
target (AF2-IG & Protenix anti-correlate r=−0.41). ODesign — an all-atom interaction world-model (NOT AF2-based)
with explicit epitope/hotspot control — was the "different generator" lever. It ran end-to-end on a free Colab T4.

## What ran
- **Generator:** ODesign `odesign_base_prot_flex` (arXiv 2510.22304), inverse-folding via ProteinMPNN + LigandMPNN.
- **Target:** KRAS-G12D neoantigen pMHC — MHC chain A (1–108) + peptide chain B `VVVGADGVGK` (the G12D 10-mer;
  WT is `VVVGAGGVGK`), **hotspot B/4–8 centred on B/6 = the G12D Asp**.
- **Designed binder:** chain C, length **80** (the KRAS-v2 winner length).
- **Sweep:** seeds `[42,123,777,2024,31337]` × `N_sample=10` ⇒ **50 designs**. Exit 0, ~21 min.
- Artifacts: `outputs.zip` (full run: 50 CIF complexes + configs + logs), `odesign_binders.fasta` (50 binder
  sequences), `odesign_design_metrics.csv` (50 designs ranked by hotspot engagement).

## Audit (rule 5 — verified before claiming anything)
| Check | Result |
|---|---|
| Chain topology | A=108 (MHC), B=10 (peptide), C=80 (binder) — exact to spec ✓ |
| Target identity | peptide `VVVGADGVGK`, position 6 = **D** → designing against the **MUT** (not WT) ✓ |
| Binder length | all 50 = **80 residues** ✓ |
| Mode collapse | **50/50 sequences UNIQUE**; protein-like composition (A12/L11/E8/R7/V7/S6 %); no poly-X ✓ |
| **Epitope control** | median binder↔G12D-Asp distance **3.9 Å**; **44/50 contact the mutation <5 Å**; median **4** binder residues at the hotspot interface; **1/50** missed entirely ✓ |

**Top designs by hotspot engagement (score these first):** seed2024/bb7, seed777/bb9 (8 interface residues each),
seed31337/bb2, seed123/bb1 (7 each). Full ranking in `odesign_design_metrics.csv`.

## What this IS — and what it is NOT (honest framing)
- **IS:** ODesign cleared the **generation + epitope-targeting** bottleneck that bounded PXDesign — 50 diverse
  80-mer binders, 88 % of them physically docked onto the G12D Asp. The binders are *positioned to read the
  mutation*, which is the structural prerequisite for discrimination.
- **IS NOT:** a discriminating binder, and not a certified one. Two things remain **UNTESTED**:
  1. **MUT-vs-WT discrimination** — contacting position 6 (D in MUT, G in WT) is *necessary but not sufficient*;
     a binder could still grip WT equally. This is THE test (a WT-binding binder attacks normal tissue, R32).
  2. **Dual-oracle binding** — ODesign's own predicted complex is **not** an independent binding oracle. The bar
     (set in KRAS-v2) is co-certification by **Protenix AND AF2-IG**.

## Next step (the real test) — MUT-vs-WT discrimination scoring
For each top binder sequence (`odesign_binders.fasta`), fold against **MUT** pMHC
(`../staging/kras_g12d_A1101_free_mut_pmhc.pdb`) **and** **WT** pMHC (`../staging/kras_g12d_A1101_free_wt_pmhc.pdb`)
on **both** Protenix and AF2-IG. **WIN = high on MUT, low on WT, on BOTH oracles.** A dual-certified discriminating
binder = a genuine breakthrough artifact. Another 0 = a strong second-generator confirmation that the external
neoantigen binder is in-silico-bounded → recognition load rests on the validated **internal key** (R27/R33/R34–40).

## Residuals (named, not hidden)
- In-silico fold/affinity ≠ measured (binder SPR/cellular = wet-lab residual).
- ODesign internal confidence is not calibrated as a binding probability; ranking here is by *geometric* hotspot
  engagement, a proxy for "engaged the right epitope," not for affinity or specificity.
- Discrimination is structural-prediction-based; the MUT/WT ΔΔG is the load-bearing unmeasured quantity.
