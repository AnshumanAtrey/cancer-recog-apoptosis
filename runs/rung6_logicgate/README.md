# RUNG 6 — the genetic NOT-gate (HLA-LOH) addressability test

**Question:** RUNG-5 proved no single/2-input *surface* AND gate is worst-donor-safe → a **100% per-patient
addressability gap**. The field's answer is a NOT-gate on a tumour-specific **genetic loss** (A2 Bio's Tmod
kills cells that have *lost* HLA-A\*02). **How much of the 100% gap does the genetic NOT-gate actually close?**

This is **arm (b)** of RUNG-6 — the laptop-feasible test that needs no GPU and no atlas re-fetch. It is a
**synthesis**: it joins *our* worst-donor surface gap (RUNG-5) to the *field's* published per-patient HLA-LOH
calls. Arm (a) — the GPU 3-input surface AND-NOT sweep — is the separate, compute-heavy confirmatory piece.

## Result (`rung6_hla_loh_addressability.json`, `rung6_hla_loh.png`)

| cancer | n | any-HLA-LOH ceiling | **HLA-A\*02 Tmod (deployed gate)** | gap after A\*02 gate | over-statement |
|---|---|---|---|---|---|
| NSCLC (pooled) | 587 | 27.8% [24.8–30.9] | **4.9% [3.6–6.6]** | 95.1% | 5.6× |
| breast (BRCA) | 938 | 13.5% [11.8–15.5] | **3.1% [2.3–4.1]** | 96.9% | 4.4× |
| colorectal (COREAD) | 662 | 20.4% [17.9–23.1] | **5.7% [4.4–7.4]** | 94.3% | 3.6× |

*(brackets = Jeffreys 95% CI.)*

**Finding:** the *generous* "any-HLA-LOH" ceiling (13–28%) over-states the **actually-deployed single-allele
A\*02 gate by 3.6–5.6×**. The A\*02 Tmod gate addresses only **~3–6% of patients** because it needs loss of the
*specific* sensed allele (only ~45% carry A\*02, ~39% are heterozygous, fewer lose that exact allele). **A
single-allele genetic gate barely moves the 100% gap.** Implication: closing the gap needs a **panel** of
allele-specific blockers (A\*02 + B\*07 + …) to recover the any-HLA-LOH ceiling — the next design step.

## Honest caveats (in the JSON `CEILING` field too)

1. **Different modalities / patients.** Surface gap = scRNA *mRNA* atlas; LOH = *WGS genotype* cohort. Same
   cancer *types*, not the same individuals. This is an integration argument.
2. **Patient-level, not clonal.** Subclonal HLA-LOH (TRACERx: ~76% of LUAD LOH) means the true addressable
   fraction is **lower** than reported — these are upper bounds.
3. **Lung is NSCLC-pooled**, not LUAD (RUNG-5's cancer). NSCLC here = 27.8%, which happens to match the
   TRACERx LUAD estimate (~29%).

## Provenance

- Data: `data/refs/mjimenez2023_MOESM6.xlsx`, sheet `GIE per sample` (6,319 WGS tumours).
  Martínez-Jiménez et al. 2023, *Nature Genetics*, DOI 10.1038/s41588-023-01367-1. Re-download (no auth):
  `curl -L -o data/refs/mjimenez2023_MOESM6.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01367-1/MediaObjects/41588_2023_1367_MOESM6_ESM.xlsx"`
- Code: `scripts/24_hla_loh_addressability.py` (`selftest` = 8/8; `run` = real). Log in `runs/logs/`.
