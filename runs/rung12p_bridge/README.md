# RUNG 12P / bridge — does the gated relay UNLOCK RUNG-11's "too-risky" neoantigen handles?

Joins the two recent results. **RUNG-11**: 175 presented public-neoantigen handles, but 104 are
`tcr_dependent` — the wild-type peptide is also presented on healthy cells, so a per-cell TCR-T's safety rests
on the TCR discriminating two near-identical pMHC surfaces (MAGE-A3 failure class; per-cell bar = 0.02).
**RUNG-12P/B**: a per-hop-gated death wave tolerates a per-step false-positive up to **q_n ≈ 0.17 (3D)** before
normal-tissue kill exceeds 1%. This bridge asks: **how much addressability does switching to the relay
architecture unlock, per cancer?**

## Method (transparent proxy; the one real unknown is swept)
Per handle: `q_n = presentation_factor(wt_rank) × β × (anchor_discount if anchor)`.
- `presentation_factor(wt_rank)` — how present the WT pMHC is on normal cells (1 at strong WT binding, 0 once
  WT %rank > 2, so **clean handles → q_n = 0 → always safe**).
- `β` — the TCR's cross-bind onto a single-residue-different pMHC. **This is exactly what RUNG-12
  (AlphaFold-Multimer/ESM) would measure per handle; here it is UNKNOWN, so it is swept over [0,1].**
- anchor mutations get a 0.2 discount (WT pMHC conformationally distinct → better discrimination).

A handle is **per-cell usable** if `q_n ≤ 0.02`, **relay usable** if `q_n ≤ 0.17`.

## Result
- **Relaxation factor: ~8.6× (3D).** The relay tolerates an ~8.6× higher TCR cross-bind than a per-cell gate —
  it converts "safe only if the TCR is near-perfect (β ≤ 0.02/pf)" into "safe if the TCR is *decent*
  (β ≤ 0.17/pf)." This number is β-independent and robust.
- **Addressability unlocked at an imperfect-TCR operating point (β = 0.5):**

  | cancer | per-cell usable | relay usable | unlock |
  |---|---|---|---|
  | GLIOMA | 22.1% | **52.6%** | **+30.5%** |
  | MELANOMA | 22.8% | 31.1% | +8.3% |
  | NSCLC | 8.3% | 13.9% | +5.6% |
  | PDAC | 43.7% | 48.8% | +5.1% |
  | CRC | 19.8% | 24.4% | +4.7% |
  | BRCA | 11.4% | 15.7% | +4.3% |

- **22 handles unlocked** (relay-safe but not per-cell-safe). The standout is **IDH1 R132H** (the ~70% glioma
  driver) on A\*03:01/A\*26:01/A\*68:01/C\*07 — anchor handles whose q_n (0.07–0.10) sits between the per-cell
  bar and the relay ceiling. That single driver is why glioma's unlock is so large. Also BRAF V600E/A\*68:01,
  KRAS G12D/C\*05:01.
- **Shape (figure):** as the TCR degrades (β rises), per-cell coverage collapses fast while relay coverage
  stays high until the percolation-protected ceiling — the relay buys a wide tolerance band.

## Honest ceiling
`β` is swept **globally**, not measured per handle — pinning it down per handle is **RUNG-12's** job
(AlphaFold/ESM). `presentation_factor` is a transparent rank→presentation proxy; the relay ceiling inherits
RUNG-12P/B's percolation-abstraction caveats; frequencies inherit RUNG-11's (literature estimates, joined
datasets). This is an **integration / what-if bridge, not per-handle truth.** The robust parts are the
**relaxation factor** and the **shape**; the exact per-cancer % at a chosen β is illustrative.

## Provenance
`scripts/36_relay_neoantigen_bridge.py` (selftest 10/10). Consumes `runs/rung11_neoantigen/…json` +
`runs/rung12pB_relay/…json`; reuses `scripts/33` HLA panel + coverage(). Laptop, instant. Outputs:
`rung12p_bridge.json`, `rung12p_bridge.png`. **Next:** RUNG-12 (AlphaFold/ESM) to replace swept β with a
measured per-handle cross-bind — turning this what-if into a ranked, certified target list.
