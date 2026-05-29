# Step 3 — Anchor Specificity / Safety Audit: methodology (critique-hardened)

**What it proves (falsifiable):** for each candidate ANCHOR receptor, whether a real therapeutic
window exists — restricted normal expression AND low expression in the **essential parenchyma of
vital organs** — so a binder that locally clusters DR5 to trigger apoptosis would spare healthy
tissue. DR5/TNFRSF10B is the fixed *trigger*, not audited as an anchor.

This design was **adversarially stress-tested before building** (research → design → critique
workflow). The first draft was rejected; what follows is the corrected method, validated against
**real GTEx v10 + HPA data** (2026-05-29).

## Why the obvious method is wrong (and was rejected)

An absolute "TPM ≥ 10 in any vital organ → FAIL" gate **fails every surface receptor** — all are
expressed somewhere vital. On real GTEx v10 it failed 17/20 genes including FDA-approved targets
(mesothelin, claudin-18.2, Trop2, folate-receptor-α). A test that condemns everything passes the
"HER2 must fail" check trivially while carrying **zero information** — the inverse of faking
specificity, but equally useless. **Bulk dilution is also symmetric:** it can hide a vital cell type
(false-PASS) *and* make a stromal/endothelial marker look ubiquitous (false-FAIL, e.g. EPHB4).

## The two-axis metric (validated on real data)

| axis | source | meaning |
|---|---|---|
| **τ (Yanai tissue-specificity)** | GTEx v10 median TPM, clean bulk panel (drop `Cells_*` + single-nucleus sub-tissue cols → 53 tissues), on log2(TPM+1) | 1 = restricted (anchor-favourable), 0 = ubiquitous |
| **vital-parenchyma IHC** | HPA `normal_ihc_data` | max protein level (Not detected<Low<Medium<High → 0–3) in the **essential parenchymal cell types** of vital organs (cardiomyocytes, hepatocytes, pneumocytes, neurons, renal tubule/glomerulus, islet/acinar, myocytes, haematopoietic) — **not** endothelium/fibroblast/immune that are everywhere |

**Verdict (pre-registered, calibrated):**
- **FAIL** if vital-parenchyma ≥ 2 (Medium/High in an essential vital cell) **OR** τ < 0.55
- **PASS** if τ ≥ 0.70 (no-IHC-data → RNA-supported, tagged *protein-unconfirmed*; the critique's
  "absence of antibody data ≠ unsafe")
- **CAUTION** otherwise
- **FORM-DEPENDENT** flag for MUC1 (tumour glycoform) / CD44 (v6 splice) — gene-level can't resolve.

## Calibration on a labelled benchmark (not guessed thresholds)

Thresholds were chosen to **separate clinically validated targets from known-bad ones** on real data:

- **good** (validated/tolerated): CLDN18.2 τ0.92, BCMA 0.87, MSLN 0.84, FOLR1 0.83, NECTIN4 0.82, DLL3 0.79, Trop2 0.71 → all **PASS**
- **bad** (ubiquitous/known-tox): HER2 0.34, CD74 0.29, CD44 0.44, ITGB4 0.48, SDC1 0.66 → all **FAIL**

**Sensitivity 7/7, specificity 5/5.** τ ≥ 0.70 cleanly separates the two sets.

## Run-trust gate (verdicts published only if ALL pass)

1. **HER2/ERBB2 → FAIL** on cardiomyocytes (Medium IHC). *The falsification anchor.* (verified: τ0.34, heart cardiomyocyte Medium)
2. **good benchmark PASS ≥ 5/7** — proves the test discriminates, is not a uniform-FAIL machine.
3. **bad benchmark FAIL = 5/5.**
4. **DR5/TNFRSF10B non-selective** (τ < 0.55) — confirms the pipeline doesn't manufacture specificity; DR5 is correctly the trigger, not an anchor. (verified: τ0.29)

If any control fails → `RUN NOT TRUSTED`, verdicts withheld, thresholds re-examined (never tuned to a favourite).

## Correction to the critique's own assumption

The critique assumed **EPHB4** was endothelial-confined (a false-FAIL to guard). **Real HPA IHC
refutes this:** EPHB4 is Medium in cardiomyocytes and High in hepatocytes → it genuinely FAILs for an
apoptosis-triggering bispecific. We did **not** hard-code an EPHB4 pass; the data decides. (sEphB4-HSA
tolerability is mechanism-specific — a soluble decoy, not a cell-killing agonist.)

## The honest result

**None of the 10 Step-2 candidate anchors has a clean therapeutic window** — all are broadly
expressed (τ < 0.70) and/or hit vital parenchyma (Medium/High IHC): HER2 (cardiomyocytes), EPHB4
(cardiomyocytes+hepatocytes), ERBB3 (Purkinje High), DDR1 (kidney), CD44 (marrow High), etc.

This is **not** a broken test — the validated benchmark PASSes. It is a true finding: **cancer
over-expression (Step 2) does not imply a safety window.** The safe anchors come from the
tissue-**restricted** tumour-antigen class (CLDN18.2, MSLN, FOLR1, DLL3, NECTIN4, Trop2), several
relevant to our cancers (Trop2/NECTIN4 → breast; MSLN/FOLR1/DLL3 → lung; CLDN18.2 → GI).

**Recommendation for Step 4:** (a) pivot the anchor to a tissue-restricted antigen that PASSes here
and matches our tumour types; and/or (b) combinatorial **logic-gating** — require co-binding of two
individually-sub-threshold antigens (Perna & Sadelain) — so specificity is the AND of two signals.

## Data sources (authoritative, downloadable, verified live 2026-05-29)

- **GTEx v10** bulk gene median TPM: `https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`
- **HPA** normal-tissue IHC (cell-type resolved): `https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip`
- Cancer side reuses Step-2 (`targets_surface_shortlist.csv`); HPA `cancer_data.tsv.zip` available for protein cancer evidence.

## Reproduce

`python scripts/07_specificity_audit.py` (downloads GTEx+HPA to `data/specificity/`, ~14 MB) or via
`notebooks/step3_specificity_colab.ipynb`. Output: `data/specificity/specificity_audit.csv`.
