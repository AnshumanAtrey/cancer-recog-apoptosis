# RUNG 8 — normal-tissue HLA-I heterogeneity (grounds RUNG-7's safety parameter)

**The measurement RUNG-7 said we were missing.** RUNG-7's gate-safety result rode on one unsourced number —
the fraction of normal vital cells that are **HLA-low** (so the Tmod blocker fails and the broad activator
kills them). This run **measures** it from the CELLxGENE atlas, per vital cell type, worst-donor, and feeds it
back into RUNG-7 for a **data-grounded per-organ off-tumour-toxicity floor**.

## How to run
Open **`notebooks/rung8_hla_heterogeneity_colab.ipynb`** in Colab (a **CPU runtime is fine** — see GPU note),
run cells top to bottom. Outputs land here: `rung8_hla_heterogeneity.json` + `rung8_hla.png`. The bundle
(`rung8_run_<ts>_<sha>.zip`) is filed with `python scripts/archive_colab_run.py --commit`.

## What it computes
Per normal tissue → vital cell types (`scripts/18` VITAL_NONREGEN via `scripts/17` VITAL_AUDIT) → HLA-A/B/C
per cell, donor-resolved → per **(vital type, donor)** HLA-low fraction (UMI below threshold) + detection.
**Headline = worst-donor HLA-low fraction per vital type** (never pooled, per RUNG-5/6). HLA-A is the sensed
gene (the deployed blocker senses HLA-A\*02). The measured worst value is plugged into RUNG-7 → `data_grounded_FPR`.

## The three engineering requirements (by request)
- **Resumable across the 4-hour cap** — one Drive tile per tissue (`RUNG8_CACHE`). A disconnect loses nothing:
  re-run Cell 5 and it skips completed tissues and continues. (Only 3 genes pulled → light.)
- **Foreground-visible logging** — a background `Heartbeat` thread prints `[heartbeat] <step> | RAM` every ~20s
  (plus a flushed `[+s][rung8]` line per step), streamed live to the cell via `runlog` (`python -u`). You
  always see the current step and that it isn't stuck.
- **GPU not used, by design** — only 3 genes; the bottleneck is the Census fetch (network/disk) and the
  aggregation is a trivial numpy groupby. Stated honestly rather than bolting on idle GPU code.

## v2 fix (after the first real run, 2026-06)
The first run (`colab_runs/20260605T191641Z_ae26e7e/`) had a **methodological bug**: CELLxGENE returns `0`
for both "measured & truly zero" **and** "this dataset never measured the gene" — so datasets lacking HLA-A
were scored as 100% HLA-low (adrenal showed an impossible 0% detection; every type hit worst-donor 100% off a
single artifact donor; the coupled FPR of 29% was fake). **v2 counts only datasets that actually measured HLA**
(≥1 cell with HLA>0 anywhere), reports the **pooled** fraction among measuring datasets as the headline (robust),
and keeps per-type median/p90/worst-donor for context. **Re-running is free** (re-aggregates from the cached
Drive tiles — no refetch). The gold-standard fix (Census `feature_dataset_presence_matrix`) is noted as a
future refetch option.

## Honest ceiling
mRNA HLA ≠ surface MHC-I protein; scRNA **dropout inflates** the HLA-low fraction (headline is an **upper
bound** → conservatively over-estimates toxicity, the safe direction); HLA-I is **IFN-γ inducible** (resting
atlas may understate induced levels); scRNA resolves the **HLA-A gene**, not the **A\*02 allele**. After the v2
filter, residual zeros from *measuring* datasets are the best mRNA estimate; protein-level confirmation is the
wet-lab residual.

## Provenance
`scripts/29_hla_heterogeneity.py` (selftest 10/10, validated locally on M2). Census version pinned to
`2024-07-01` (matches RUNG-5). Reuses `scripts/17` tissue/vital conventions + the Arrow-dictionary memory-safe
fetch. Logs in `runs/logs/`; immutable per-run archive in `colab_runs/`.
