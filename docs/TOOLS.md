# TOOLS — computational design/fold/MD arsenal (with HONEST verification status)

Consolidated from Anshuman's mid-2026 cross-check + tools we've actually run. **Rule 5 applies: a tool being on
this list ≠ verified-working for us.** Status column = MY confirmation level, not the source's claim. Verify the
repo installs + runs (or use the webserver) BEFORE any GPU spend (the project has burned sessions on this).

**Status:** 🟢 used/verified by us · 🔵 repo+paper confirmed by me (not yet run) · 🟡 catalog (Anshuman-provided, UNverified by me — check repo before use)

## What we ACTUALLY use (the live pipeline)
| Role | Tool | Status | Notes |
|---|---|---|---|
| Binder generation | **PXDesign** (ByteDance) | 🟢 | webserver protenix-server.com (A100s, free). Made strong pMHC binders R26c/d. Add Prediction = fold; Add Design = binder. |
| Folding / discrimination oracle | **Protenix** (ByteDance, AF3-class) | 🟢 | webserver. Folded IDH1 + BRAF + (queued) PIK3CA/KRAS pMHC. seed-match MUT/WT for specificity. |
| AF2 cross-check scoring | **ColabDesign / AF2** (mk_afdesign_model, binder protocol) | 🟢 | the `score()` in every RUNG-26 run; no-multimer config under-calls binding (within-model specificity still valid). |
| MHC presentation | **MHCflurry 2.x** | 🟢 | local venv `/tmp/mhc`; the RUNG-28 screen + BRAF/PIK3CA presentation calls. |
| Allele-specific CRISPR | own scanner (SpCas9-NGG/NG, Cas12a) | 🟢 | scripts/54–56; + off-target scanner scripts/61/62 (pigeonhole, GRCh38). |
| RFdiffusion + ProteinMPNN | dl_binder_design recipe | 🟢 | runs on a free T4 (R26b, 50 backbones). ProteinMPNN reused for negative design (R26f). |

## Confirmed by me (repo + recent paper), candidate for the fold→design stage
| Role | Tool | Status | Notes |
|---|---|---|---|
| 2nd binder generator (epitope-specified) | **ODesign** | 🔵 | github **The-Institute-for-AI-Molecular-Design/ODesign**, arXiv 2510.22304, Apache-2.0. (Anshuman's list mis-attributed it to "LG AI Research" — corrected.) All-atom, hotspot-specified; FOLD-BOUND (needs the pMHC). Use to cross-check PXDesign on PIK3CA/KRAS once folds land. |
| Phase-2 binder stability / MD | **BioEmu** (Microsoft) | 🔵 | Science 2025; pip, single-GPU, T4-friendly; monomer-native (complex off-rate approximate). Run on a binder once it clears the static bar. |

## Catalog (Anshuman-provided — UNVERIFIED by me; verify repo before spending)
| Tool (claimed origin) | Claimed role | Status |
|---|---|---|
| **SeedProteo** (ByteDance) | all-atom side-chain interface design | 🟡 plausible (same team as PXDesign) — verify repo |
| **BindCraft** (EPFL) | one-shot binder, secondary method | 🟡 real tool, but needs ≥32 GB GPU → webserver/bigger GPU only (not free T4) |
| **Genie 3** (AQLab) | fast backbone gen / screening | 🟡 verify |
| **Protenix-v2** (ByteDance) | AF3-surpassing open fold | 🟡 we use Protenix (webserver); v2 specifics unverified |
| **HelixFold3** (Baidu) / **OpenFold3** (AlQuraishi) | alt AF3-class folders (tertiary cross-check) | 🟡 verify |
| **VibeGen** (MIT/CMU) | dynamics-as-design (loop flexibility) | 🟡 future frontier; verify |
| **OVO** (MSD CZ) | industrial orchestrator (HPC/SLURM, glue) | 🟡 verify repo before relying on it |
| **ProteinDJ** | high-throughput parallel pipeline | 🟡 verify |
| **Mosaic** (Escalante Bio) | multi-objective (bind+stability+expression) | 🟡 verify |
| **SALAD** (EMBL) | ~8M-param fast backbone gen, up to 1000 res | 🟡 verify |

## Status updates (from Anshuman's cross-check — not independently confirmed)
- **Latent-X** → pivoted to an agent platform (Latent-X2 → "Latent-Y"). **BinderFlow** → no active primary repo.
  **Bolt2Design** → conceptual only, no production repo. (Treat as Anshuman-reported until verified.)

## How this maps to our plan
1. **Generate**: PXDesign (primary) + ODesign (epitope-specified cross-check) once PIK3CA/KRAS folds land.
2. **Fold/score**: Protenix + ColabDesign-AF2 (specificity = MUT-vs-WT, seed-matched).
3. **Stability**: BioEmu on any binder that clears pae≤10 / plddt≥80.
4. Orchestrators (OVO/ProteinDJ) only if we scale past hand-driving — verify-repo first; not needed yet.
