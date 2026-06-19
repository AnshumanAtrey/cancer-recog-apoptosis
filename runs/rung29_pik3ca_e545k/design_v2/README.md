# RUNG-29 v2 — PIK3CA-E545K/A*03:01 binder design (PXDesign Extended, hotspot B4+B6, max batch)

Re-run of the presentation-flip target after v1 returned **0/10**. v1 was under-powered (no hotspot, ≤5% tier).
v2 added the **up-facing hotspot (B0:4 + B0:6)** from the staged MUT pMHC and ran PXDesign **Extended** with the
`Protenix-Mini-Templ` model (task_info.json).

## Result (summary.csv, 10 designs)
| metric | result |
|---|---|
| AF2-IG-easy passers | **0 / 10** |
| AF2-IG passers | 0 / 10 |
| Protenix passers | 0 / 10 |
| **Protenix-basic passers** | **1 / 10** (rank_1) |
| dual-oracle (AF2 **and** Protenix) passers | **0 / 10** |

**rank_1** (the only passer): `ptx_iptm 0.871`, `ptx_iptm_binder 0.803` (Protenix says it binds) but
`af2_iptm 0.14`, `af2_ipAE 25.5`, `AF2-IG-easy False` (AF2 says it does **not**). The two orthogonal oracles
**disagree hard** — the classic single-model-artifact signature. No design clears both.

## Difficulty gauge (`difficulty_gauge.png`, AF2 initial-guess estimate)
The submitted target sits in the **≤5% hardest tier on BOTH** AF2-IG-easy and Protenix-basic passing-rate scales
("your job is here" → ≤5% on both). Our 0/10-AF2, 1/10-Protenix outcome is exactly what a ≤5% target produces at
this batch size — the gauge and the result corroborate each other.

## Verdict — honest
**PIK3CA-E545K external binder remains UNCONFIRMED.** Even at the up-facing hotspot, this is a genuinely hard
read-the-mutation target: no dual-passing binder, and the single Protenix-basic passer fails AF2 cross-validation.
This is *not* a breakthrough binder. It does not refute the presentation-flip route (the flip itself is real and
MS-corroborated, C8/R32) — it says the *binder-design* step on this 11-mer groove is at the edge of what the
current de-novo tools can do at small batch.

**Why this is the expected outcome, not a failure of method:** the screen (C5/R28) picked PIK3CA on a
*presentation* criterion (WT unpresented → auto-specificity), which removes the MUT-vs-WT discrimination burden —
but it does **not** make the binder *easy to fold against the groove*. Foldability and discrimination are
independent axes; the screen optimized the second, the gauge measures the first.

## Next (does NOT block the internal key)
- **KRAS-G12D/A*11:01** (C7) is the staged companion target (hotspot B6, free-fold). Run it the same way; there
  the make-or-break is MUT-vs-WT scoring (KRAS WT *is* presented, R32), a different and arguably more decisive test.
- For PIK3CA specifically: a **larger batch** (the ≤5% tier needs hundreds of designs to surface dual-passers) or
  a **longer binder length** band would be the honest next lever — single-oracle passers at batch=10 don't count.

*Raw committed before this analysis (commit f2bf779). Ledger row C6 updated.*
