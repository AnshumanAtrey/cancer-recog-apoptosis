# RUNG 12 — pMHC structural discriminability: measure β, certify relay targets

Closes the loop opened by the RUNG-12P bridge. That bridge showed a gated relay unlocks RUNG-11's
`tcr_dependent` neoantigen handles **if** the TCR's mut-vs-WT cross-bind `β` is low — but `β` was *swept*.
This rung **measures** a per-handle `β` and re-runs the bridge with it → a **ranked, certified target list**.

## How β is measured (discriminability D = 1 − β)
Probabilistic-OR of independent discrimination mechanisms (D high if *any* fires):
- **M — MHC-level** (from RUNG-11 NetMHCpan): WT binds much worse / not presented → WT pMHC rarely on the
  surface. `clean` tier → M = 1. *(robust, fold-independent)*
- **E·P — TCR-level**: the mutated residue is solvent-**E**xposed (RSA from the ESMFold-bound conformation)
  **and** physicochemically different (**P** = charge/volume/hydropathy delta). *(P robust; E best-effort)*
- **Z — sequence**: ESM-2 embedding distance between mutant and WT peptide, normalised across handles. *(robust)*

Then `q_n = presentation_factor(wt_rank) · β` → **per-cell-safe** (q_n ≤ 0.02) / **relay-safe** (q_n ≤ 0.17,
the RUNG-12P/B 3D ceiling). Re-running coverage with the measured β gives the certified usable target set.

## How to run (Colab, **T4 GPU**, ~1–2 h, resumable)
Open `notebooks/rung12_pmhc_colab.ipynb`, **Runtime → Run all**, same Google account as before. Stages:
1. **prep** — selects the top ~32 handles by prevalence (surfaces IDH1 R132H/glioma, KRAS, BRAF V600E),
   fetches HLA α1α2 grooves from IPD-IMGT/HLA, writes `groove:peptide` ESMFold inputs. *(validated locally)*
2. **ESM-2 embeddings** → per-handle Z. *(robust core — no MSA server)*
3. **ESMFold structures** → bound pMHC PDBs → RSA exposure E. **Best-effort**: ESMFold's openfold extras
   often fail to build on Colab and its peptide docking is the soft part; if it fails, the run still produces
   a valid β from M + P + Z (E falls back to a position prior). Resumable per handle.
4. **analyze** → `rung12_pmhc.json` + `rung12_pmhc.png` (ranked targets + measured-β bridge coverage).
Bundle with `python scripts/archive_colab_run.py --commit`.

## Honest ceiling
ESMFold is a single-sequence **model** — short-peptide docking into the MHC groove is unreliable, so **E is
the softest signal**; pLDDT is reported and the β estimate leans on the fold-independent M + P + Z. **β is a
proxy, not a measured TCR Kd** — a top-ranked handle is a *prioritised hypothesis for wet-lab TCR isolation*,
not a validated target. Inherits RUNG-11 (population frequencies, mRNA→presentation) and RUNG-12P
(percolation-abstraction relay ceiling) caveats.

## Provenance
`scripts/37_pmhc_discriminability.py` (selftest 16/16; `prep` validated end-to-end against IPD-IMGT/HLA —
32/32 grooves). ESM-2 (`fair-esm`), ESMFold (best-effort), Biopython SASA. HLA grooves cached in
`data/refs/{A,B,C}_prot.fasta`. Consumes `runs/rung11_neoantigen/…` + `runs/rung12pB_relay/…`. **Next:** the
top-ranked relay-safe handles become the wet-lab shortlist; β can later be replaced by measured TCR affinities.
