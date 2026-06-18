# RUNG-29 design v1 — PIK3CA-E545K binder, PXDesign Extended, NO hotspot → 0/10 (under-powered, fix known)

First binder batch for the screen-chosen presentation-flip target. **Honest result: 0/10 designs pass.**

## Result
| metric | value |
|---|---|
| designs | 10 (Extended mode, Protenix-Mini-Templ, length 100-120) |
| AF2-IG-easy pass | **0/10** |
| AF2-IG strict pass | 0/10 |
| Protenix-success | 0/10 |
| best af2_ipAE | 20.29 (bar ≤10; IDH1's passers were ~5.7) |
| ptx_iptm range | 0.72–0.84 (moderate, but Protenix-success=False) |

All 10 are well-folded mini-proteins (low bound↔unbound RMSD 0.2–0.5) that **do not grip the pMHC** (interface
pAE ~20–22). No strong binder.

## Why (diagnosed — NOT a dead target)
1. **Difficulty gauge: ≤5% passing tier** (hardest, both AF2-IG + Protenix — same as IDH1). At ≤5%, **10 designs →
   ~0.5 expected passers → 0 is statistically expected.** The batch was under-powered, not the target dead.
2. **No hotspot (free design).** My v1 call was "no hotspot" — correct that we must NOT pin the *buried* mutation
   (p11 Lys, 1% exposed = C-term anchor), but I wrongly omitted a hotspot entirely, so PXDesign made nice
   mini-proteins that ignore the peptide. The gauge itself says to adjust hotspot/length.

## Fix → design v2 (prescribed)
- **Hotspot = B4 + B6** — the genuinely up-facing peptide residues (SASA: p4-Asp 38%, p6-Leu 59% exposed;
  p7/p9/p10 mid; everything else incl. p11 buried). This directs the binder onto the *presented peptide* surface
  without pinning the buried mutation. (Selectivity still comes from presentation — WT isn't displayed, R32 — so
  we don't need mutation-discrimination, just a strong MUT-pMHC binder that engages the peptide.)
- **Max batch (100 designs)** — to clear the ≤5% rate (expected ~5 passers).
- Length 100-120 (as used) fine; same cropped target `staging/pik3ca_e545k_mut_cropped.pdb`.

Raw v1 archived bit-for-bit (summary.csv, task_info.json, 10 orig_designed CIFs, difficulty_gauge.png).
