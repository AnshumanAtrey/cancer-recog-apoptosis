# RUNG-30 KRAS-G12D / A*11:01 binder v2 — the 80-design pre-registered test (PXDesign Extended, hotspot B0:6)

v1 was an under-powered 0/10 at a ≤5% difficulty tier. v2 ran the pre-registered **large-batch** test: 8 jobs ×
10 designs = **80 designs**, multi-seed (36461/19931/24601/13337/80085/42424/31415/27182), binder length 70–100
(concentrated at 80, the v1.5 winner), same crop `A0:1-108 B0:1-10` and hotspot `B0:6` (the up-facing G12D Asp).

## Result (80 designs pooled)
| metric | result |
|---|---|
| AF2-IG-easy passers | **0 / 80** |
| AF2-IG passers | 0 / 80 |
| Protenix passers | 0 / 80 |
| Protenix-basic passers | **1 / 80** (s3 rank1) |
| **dual-oracle passers** | **0 / 80** |

- best `af2_iptm` **0.70** (s7 rank3) — a big jump from v1's 0.37: **binders now genuinely engage** the target
  (and far from IDH1's NULL). Length 80 designs dominate the top grips.
- best `ptx_iptm` **0.86**; the one Protenix-basic passer (s3 rank1) has `af2_iptm 0.15` (AF2 says no).
- **No design passes AF2's binding bar at all** — minimum `af2_ipAE` across all 80 is **16.1** (bar ≤10).

## Why 0 dual-passers — the mechanism (not just bad luck)
The two orthogonal oracles **actively disagree** on these designs:
**Pearson r(af2_iptm, ptx_iptm) = −0.41.** Designs AF2 likes, Protenix dislikes, and vice versa.
- `af2_iptm > 0.5`: **16/80** · `ptx_iptm > 0.8`: **5/80** · **both: 0/80** (the two winner-sets do not overlap).

So a design satisfying both models is **structurally rare here**, not merely under-sampled — more seeds would keep
hitting the same anti-correlation. The single best "both-agree" candidate (s7 rank3: af2_iptm 0.70 / ptx_iptm
0.69, af2_ipAE 17.8) still fails both formal bars.

## Verdict — honest: BOUNDED at this method (not refuted like IDH1, not confirmed)
**The KRAS-G12D external binder is bounded under PXDesign (AF2-IG + Protenix) at 80 designs.** It is a *softer*
negative than IDH1/BRAF:
- IDH1 was a hard **NULL** (binders didn't engage at all, 0/12 across 4 methods *with* negative design). KRAS
  binders **do** engage (af2_iptm 0.70, ptx_iptm 0.86 individually) — the favorable geometry (Asp 30% exposed,
  Gly→Asp presence/absence) shows up as real grip.
- But **no design clears both orthogonal oracles**, and the oracles anti-correlate → our bar for a *confident*
  binder isn't met. The **MUT-vs-WT discrimination test never triggers** — you can't score discrimination on a
  design that doesn't pass the binding bar first. So discrimination remains **untested**, not failed.

## Strategic read (both neoantigen binders now tested at scale)
PIK3CA-E545K (v2: 1/10 single-oracle) and KRAS-G12D (v2: **0/80 dual**) are **both** ≤5%-tier and **both** yield
0 dual-passers at batch. **The external-binder route via PXDesign is bounded for these neoantigen pMHC targets** —
the orthogonal oracles don't co-certify any design even at scale and favorable geometry. This **reinforces the
project's strategic spine:** the internal MHC-free key (CRISPR R27 / 2-HG R33 → kill-coupling R34–37) is the
stronger, validated contribution; the external binder was always the harder, compute-hungry front.

## Next levers (if we keep pushing the binder; none block the internal key)
- **Different generator** (the honest remaining lever): RFdiffusion+ProteinMPNN or **ODesign** (the epitope-
  specified 2nd generator validated in R31b) — a different design prior might find designs *both* AF2 and Protenix
  certify. The PXDesign prior demonstrably can't here.
- **Or accept bounded** and consign KRAS/PIK3CA binders to "engage but not dual-certifiable in-silico" — the
  internal key carries the recognition load.
- Not recommended: more PXDesign seeds (the r=−0.41 anti-correlation means more of the same).

*Raw committed before this analysis (s1–s4 commit 2f8bf88, s5–s8 commit 4e8f453). Ledger row C7 updated.*
