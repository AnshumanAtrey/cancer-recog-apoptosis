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

## STRUCTURAL REDO — ColabFold (`notebooks/rung12_structure_colab.ipynb`)
After the ESMFold attempt failed (below), the structural arm was redone with **ColabFold / AlphaFold2-multimer**
(MSA-guided → actually docks the peptide; license-free; GPU). It folds the HLA α1α2 groove + peptide for the
**non-`clean` handles only** (clean → β=0 regardless of structure; 24 of 32 here), measures the mutated
residue's **TCR-facing exposure (RSA)**, and feeds the measured `E` into the β scoring. The RSA analysis was
**validated against a real HLA-A\*02:01 crystal (1HHK)**: it correctly reads buried anchors (P2/P9 ≈ 0.02) vs
TCR-facing residues (P4–6 ≈ 0.5). Run that notebook on **T4 GPU**, Run all; ~1.5–2.5 h, resumable.

---

## (Superseded) first attempt — `notebooks/rung12_pmhc_colab.ipynb` (ESMFold)
Open it, **Runtime → Run all**, same Google account as before. Stages:
1. **prep** — selects the top ~32 handles by prevalence (surfaces IDH1 R132H/glioma, KRAS, BRAF V600E),
   fetches HLA α1α2 grooves from IPD-IMGT/HLA, writes `groove:peptide` ESMFold inputs. *(validated locally)*
2. **ESM-2 embeddings** → per-handle Z. *(robust core — no MSA server)*
3. **ESMFold structures** → bound pMHC PDBs → RSA exposure E. **Best-effort**: ESMFold's openfold extras
   often fail to build on Colab and its peptide docking is the soft part; if it fails, the run still produces
   a valid β from M + P + Z (E falls back to a position prior). Resumable per handle.
4. **analyze** → `rung12_pmhc.json` + `rung12_pmhc.png` (ranked targets + measured-β bridge coverage).
Bundle with `python scripts/archive_colab_run.py --commit`.

## Result (real T4 run `89c7dfb`, 32 handles) — with two honest corrections
**The structural arm did NOT execute.** ESMFold's `openfold` extras failed to build on Colab (the fragility
flagged up front) → **0/32 structures**; exposure `E` fell back to a position prior throughout. And the ESM-2
`Z` signal, min-max normalized, initially over-claimed (the single biggest-change handle got `Z=1 → β=0`); the
scoring now **caps `Z`** (commit fix) so a relative-max embedding can't alone declare a handle perfectly
discriminable. So this result is a **binding (M) + physicochemical (P) + ESM-2 (Z, capped)** per-handle β —
**structural exposure is unmeasured.** (The raw run is preserved bit-for-bit in `colab_runs/`; the mirror here
is the Z-corrected re-derivation on the run's own signals.)

**With per-handle β (not swept), the picture is more sober than the bridge's optimistic swept-β estimate:**
- **9 per-cell-safe, 11 relay-safe, 2 unlocked by the relay** (vs the bridge's broad +5–30% at a uniform β=0.5).
- Most safe coverage sits in **`clean` handles** (WT not presented → genuinely safe today): PDAC 26%, glioma
  22% (IDH1 R132H), melanoma 11%.
- The relay still gives a **real melanoma unlock: 10.8% → 19.1%** (BRAF V600E handles that are relay-safe but
  not per-cell-safe). Elsewhere the marginal unlock is small.
- Notable: **KRAS-G12D/C\*08:02** (the proven clinical TCR target) sits at q_n ≈ 0.19, just *above* the relay
  ceiling → flagged borderline-risky. This is exactly where real structure (not available here) or the known
  exquisite clinical TCR would correct the generic estimate — an honest "we can't certify it from sequence alone."

**Takeaway:** the robust, fold-independent signals (binding + physicochemistry) already tier the targets and
temper the bridge's optimism. True *structural* discriminability needs heavier tooling than a robust 4h run
allows (ColabFold-multimer with MSAs, or PANDORA homology modelling) — single-sequence ESMFold can't reliably
dock a 9-mer into the groove even when it builds.

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
