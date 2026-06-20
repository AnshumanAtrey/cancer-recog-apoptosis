# RUNG-30 KRAS-G12D / A*11:01 binder design v1 (PXDesign Extended, hotspot B0:6, batch=10)

The decisive *read-the-mutation* target: KRAS WT `VVVGAGGVGK` **is** presented on normal tissue (R32), so unlike
PIK3CA the binder must discriminate the mutated residue. Designed against the FREE-fold cropped pMHC (p6 Asp 30%
exposed), hotspot forced onto the up-facing Asp (`B0:6`).

## Result (summary.csv, 10 designs)
| metric | result |
|---|---|
| AF2-IG-easy passers | **0 / 10** |
| AF2-IG passers | 0 / 10 |
| Protenix passers | 0 / 10 |
| Protenix-basic passers | **0 / 10** |
| dual-oracle passers | **0 / 10** |

No `passing-*` directory was emitted (0 passers on every oracle). Best design **rank_1**: `af2_iptm 0.37`,
`af2_ipAE 20.9` (bar ≤10), `ptx_iptm 0.77` — engages the target but clears **neither** oracle's bar. af2_iptm
across the batch: max 0.37 (cf. PIK3CA v2 max 0.46) — both weak, neither route confidently grips at batch=10.

## Difficulty gauge (`difficulty_gauge.png`, AF2 initial-guess estimate)
KRAS sits in the **≤5% hardest tier on BOTH** AF2-IG-easy and Protenix-basic — the *same* tier as PIK3CA. At a
≤5% passing rate, batch=10 has an expected yield of ~0.5 designs → **0/10 is the statistically expected outcome,
not evidence the target is intractable.** The gauge's own note names the levers: *hotspot location, binder length,
target file resolution* — and batch size.

## Verdict — honest (under-powered negative, NOT a refutation)
**KRAS-G12D binder v1 is an under-powered 0/10, not a closed door.** Two things are true at once:
1. No binder grips the MUT pMHC confidently at batch=10, so the **MUT-vs-WT discrimination test cannot run yet**
   — you can't score discrimination on a design that doesn't bind MUT. The decisive test is still *pending*, not
   *failed*.
2. This is *expected* at a ≤5% difficulty tier with only 10 designs. The target is **not** refuted (contrast the
   real refutations: IDH1 R26c/d/f had 0/12 across 4 methods *with* negative design; BRAF R26e was buried +
   weakly presented). KRAS is presence-vs-absence chemistry (Gly→Asp), the strongest discrimination case, and the
   Asp is 30% exposed — the favorable geometry is intact.

## Next (the real lever is batch)
- **KRAS v2 — large batch.** 🌐 Re-run identical settings (crop A0:1-108 + B0:1-10, hotspot B0:6, length 80–120,
  Extended) at the **largest batch the server allows (100–200)**. A ≤5% tier needs hundreds of designs to surface
  dual-passers; 10 was a scout, not a verdict.
- Only **after** a binder grips MUT: score it against the staged WT (`kras_g12d_A1101_free_wt_pmhc.pdb`) — MUT-good
  + WT-bad = the win.
- Strategic note: PIK3CA (1 single-oracle passer) and KRAS (0) are **both** ≤5%-tier → the external-binder route
  as a whole is **batch-size-bottlenecked**, while the internal key (CRISPR/2-HG, R27/R33) carries no such cost.
  The binder is the harder, compute-hungry front; the internal root-kill is the stronger contribution.

*Raw committed before this analysis (commit 6cb15a9). Ledger row C7 updated.*
