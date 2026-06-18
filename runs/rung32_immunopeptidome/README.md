# RUNG-32 — immunopeptidomics confirmation: our targets' presentation, checked against REAL mass-spec

The named residual on RUNG-28/29/30 was that presentation is **predicted** (MHCflurry), not **measured**. This
cross-checks our two active targets (+ the 3rd flip hit) against the **HLA Ligand Atlas** (rel 2020.12 — 90,428
HLA-I ligands from 227 *benign* tissue samples; `scripts/64`). It turns the prediction into a measurement.

## Result — both targets behave EXACTLY as their mechanism predicts
| peptide | benign immunopeptidome | reading |
|---|---|---|
| PIK3CA-E545K **MUT** `STRDPLSEITK` | **ABSENT** | ✓ tumour-only (not a benign self-peptide) |
| PIK3CA-E545K **WT** `STRDPLSEITE` | **ABSENT** | **presentation FLIP confirmed by mass-spec** — WT not presented on normal tissue |
| KRAS-G12D **MUT** `VVVGADGVGK` | **ABSENT** | ✓ tumour-only |
| KRAS-G12D **WT** `VVVGAGGVGK` | **FOUND** — A\*11:01, normal **Lung + Testis** | **no flip** — WT is a real self-peptide → binder MUST read G12D |
| TP53-R248W MUT `SSCMGGMNW` / WT `SSCMGGRNW` | both ABSENT | MUT tumour-only ✓ (WT on B\*57:01 — fewer B\*57 donors, less interpretable) |

## Why the absences are interpretable (not just coverage gaps)
The atlas has **4 A\*03:01 + 6 A\*11:01 donors** (of 21) and contains **43,670 A\*03/11 HLA-I 9–11mers** — so it
demonstrably *would* detect a peptide of our class/allele if it were presented. Therefore:
- **PIK3CA WT absent** despite right-allele donors = real-data support that WT-E545 isn't presented → the
  flip (MHCflurry 51 nM MUT / 20 µM WT, ~397×) is corroborated by immunopeptidomics, not just a predictor.
- **KRAS WT present** on A\*11:01 normal lung/testis = the germline peptide is genuine self → a KRAS-G12D binder
  that doesn't specifically grip the Asp would cross-react on normal tissue. **This both validates the
  read-the-mutation design (hotspot-on-p6) AND raises the specificity stakes for KRAS** (vs PIK3CA, where WT
  simply isn't there to cross-react with).

## Net
- **PIK3CA-E545K** (presentation flip) — MS-corroborated: a strong MUT-pMHC binder is auto-specific. ✅ de-risked.
- **KRAS-G12D** (read-the-mutation) — MS confirms WT IS displayed → binder discrimination is mandatory, not
  optional → hotspot-on-p6 is the right call, and specificity scoring (MUT-vs-WT) is the make-or-break.

## RUNG-32b — MUT-in-tumour confirmation (IEDB; the cancer-side complement) — `scripts/65`, `iedb_mut.json`
A benign atlas can't show a tumour mutation, so the MUT peptides were checked against **IEDB** (curated
neoantigen epitopes + MS-eluted ligands + T-cell assays). **Both are experimentally observed MHC ligands on
our target alleles — prediction → measurement:**

| MUT peptide | IEDB evidence | verdict |
|---|---|---|
| KRAS-G12D `VVVGADGVGK` | 18 MHC-ligand assays on **A\*11:01 (11) + A\*03:01 (6)**; **MS-eluted ×4**; **x-ray ×3**; **T-cell-positive 94/106** (A\*03:01/A\*11:01/A\*68:01); 9 papers | **GOLD-STANDARD validated, naturally-presented neoantigen** |
| PIK3CA-E545K `STRDPLSEITK` | 2 cellular-MHC **binding** assays (A\*11:01, A\*33:03), 1 Pos / 1 Neg; no MS/T-cell | **observed-binding, thinner** — presented-plausible, weaker evidence than KRAS |

**Implications for design:**
- **KRAS-G12D** target is *rock-solid* (real eluted + crystallized + T-cell-recognized). Two corollaries: (a) de
  novo binders are a *complementary modality* (natural TCRs already exist — e.g. the A\*11:01 work), so our
  contribution is the binder *artifact*, not first-discovery of the epitope; (b) **x-ray structures exist → use a
  real PDB as the design target** instead of (or alongside) the Protenix fold = higher-quality target. [follow-up]
- **PIK3CA-E545K** is the more *novel* shot (thin prior validation) but its presentation rests on MHCflurry + 1
  positive binding assay + the benign-absence flip — honestly the weaker-evidenced of the two on presentation,
  though the *binder* problem is easier (no discrimination needed). Both stay live; evidence levels now explicit.

## Honest residual
- Benign atlas = 21 donors / specific tissues; absence ≠ never-presented-anywhere (strong support, not proof).
- **MUT-in-tumour presentation** (is `STRDPLSEITK` / `VVVGADGVGK` actually eluted from real tumours?) needs a
  CANCER immunopeptidome (caAtlas / IEDB) — the next layer; a benign atlas can't show a tumour mutation.
- MHCflurry + atlas = presentation evidence; binder affinity/specificity (SPR/cellular) still the wet-lab line.

*Data: HLA Ligand Atlas rel 2020.12 (CC-BY 4.0), `data/immunopeptidome/` gitignored. Result: `immunopeptidome.json`.*
