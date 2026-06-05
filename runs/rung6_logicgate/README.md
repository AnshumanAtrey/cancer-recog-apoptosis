# RUNG 6 — the genetic NOT-gate (HLA-LOH) recognition arm

**Question:** RUNG-5 proved no single/2-input *surface* AND gate is worst-donor-safe → a **100% per-patient
addressability gap**. The field's answer is a NOT-gate on a tumour-specific **genetic loss** (A2 Bio's Tmod
kills cells that have *lost* HLA-A\*02). **How much of the 100% gap does the genetic NOT-gate actually close?**

These are **laptop runs** (no GPU, seconds) — arm (b) of RUNG-6. They are a **synthesis**: they join *our*
worst-donor surface gap to the *field's* published per-patient HLA-LOH calls (Martínez-Jiménez 2023,
6,319 WGS tumours). All three were **adversarially audited** (a 5-agent review) and the fixes below applied.

> **What is measured:** NOT-arm (HLA-LOH) **availability** — a **PATIENT-LEVEL UPPER BOUND**.
> **True addressability = (HLA-LOH availability) × (broad-activator availability) × (clonal fraction).**
> Reading any single number as "addressable patients" overstates reach; the haircuts are real and quantified.

## The three runs

### `scripts/24` — single-allele gate vs ceiling
The deployed A\*02 Tmod gate addresses only **~3–6%** of patients (NSCLC 4.9%, breast 3.1%, CRC 5.7%) vs a
14–28% any-HLA-LOH ceiling — **3.6–5.6× smaller**, because it needs loss of the *specific* sensed allele.
**Clonal haircut (now computed):** with ~76% of LUAD HLA-LOH subclonal, the clonal-only NSCLC A\*02 gate is
**~1.2%, not 4.9%**. → `rung6_hla_loh_addressability.json`, `rung6_hla_loh.png`.

### `scripts/26` — blocker-PANEL curve
Reframe: Tmod gets specificity from the **negative** arm, so reach = patients with a usable lost-allele
blocker. A **~6-blocker panel** (A\*03, A\*02, A\*01, A\*24, A\*31, A\*11) takes NSCLC from 4.6% → ~19%
(**~4×** the single allele); the full ceiling needs ~22 distinct allotypes (long tail). The A\*02:01 allotype
is the **conservative floor** (the deployed blocker is A\*02-*clade* cross-reactive). Greedy tie-break is
deterministic. → `rung6_blocker_panel.json`, `rung6_blocker_panel.png`.

### `scripts/27` — pan-cancer atlas (58 types)
Highest **LOH availability**: KICH 73.8%, PANET 54.4%, CESC 40.4% — but KICH/PANET have **NO validated broad
activator** (`*` in the figure), so those ceilings are **hollow**. The **honest target list** (high LOH *and*
a validated activator): **CESC, ESCA, NSCLC, HNSC, PAAD, OV**. Lowest: CLL 0%, pilocytic 0%, liver 2% (quiet
genomes). A universal top-12 allotype panel reaches 14.1% of all patients vs 17.2% ceiling. →
`rung6_pancancer_addressability.json`, `rung6_pancancer.png`.

## What survives, honestly
- **Solid:** the single-allele gate is tiny (~3–6%, corroborated by BASECAMP-1 ~3.2% A\*02-LOH yield); a small
  panel recovers most of the per-cancer LOH ceiling; the pan-cancer biology (KICH-top/CLL-bottom) is
  mechanistically correct and reproduces the supplement's validated LILAC LOH call to rounding.
- **The contribution is integration/quantification, not a new biological finding** (the field already frames
  addressability as HLA-LOH-bounded). Honest negatives and upper-bound labels are kept first-class.

## Irreducible residuals (cannot be done on public/laptop data)
1. **Clonal LOH per cancer** — needs controlled-access TRACERx multi-region WES (the public TRACERx sheet is
   single-region). We apply the published lung multiplier illustratively only.
2. **Activator availability + density** per patient — needs proteomics/atlas cross; here it's a literature tier.
3. **Same-patient AND-NOT** — our surface gap (scRNA) and LOH calls (WGS) are different cohorts; matched by
   cancer *type*, not individual.

## Provenance
- Data: `data/refs/mjimenez2023_MOESM6.xlsx` sheet `GIE per sample`. Martínez-Jiménez 2023, *Nat Genet*,
  DOI 10.1038/s41588-023-01367-1. Re-download (no auth): `curl -L -o data/refs/mjimenez2023_MOESM6.xlsx "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01367-1/MediaObjects/41588_2023_1367_MOESM6_ESM.xlsx"`
- Selftests: `24` 8/8, `26` 12/12, `27` 6/6. Logs in `runs/logs/`.
