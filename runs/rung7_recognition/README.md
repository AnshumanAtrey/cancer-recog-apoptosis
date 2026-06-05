# RUNG 7 — the AND-NOT recognition-gate DISCRIMINATION model

**Question RUNG-6 left open:** RUNG-6 *counted* how many patients have a usable genetic gate. It never asked
the literal recognition question — **does the Tmod AND-NOT gate actually separate a cancer cell from a healthy
one, and how does it fail?** This models that (bar #2 of the gate ladder) and couples it to the apoptosis
commit (RUNG-1 / EARM essence), giving the first end-to-end in-silico *recognise → commit to apoptosis* demo.

## The model
`activator = Hill(antigen density)` **AND-NOT** `blocker = Hill(HLA density)` → kill-license → bistable
apoptosis commit. Tumour cells: antigen-high, HLA lost (clonal) or retained (subclonal). Normal cells:
antigen mostly low + a leak tail, HLA high except a downregulated "HLA-low" fraction. Parameters are
literature-grounded (low-affinity CAR antigen-density discrimination; MHC-I levels; LIR-1 blocker); per-cell
distributions are **illustrative** (real joint spread needs the atlas → a Colab run).

## Result (`rung7_gate_discrimination.json`, `.png`)
Baseline: TPR (tumour killed) 48%, FPR (normal killed) **1.4%**, ROC-AUC 0.75. Two findings, both with
sensitivity sweeps:

1. **Safety is carried entirely by the blocker, not by LOH frequency.** Off-tumour toxicity floor ≈
   `P(normal HLA-low) × P(normal antigen-high)` (predicted 1.0% ≈ measured 1.4%). FPR tracks the normal
   HLA-low rate almost 1:1 (0%→0%, 5%→1.4%, 20%→5.8%) and persists across activator thresholds. **The gate
   can be no safer than HLA expression is reliable in normal tissue.**
2. **Efficacy is bounded by the clonal-LOH fraction.** TPR rises 19%→96% as clonal LOH goes 0.2→1.0 —
   subclonal HLA-retaining tumour cells keep the blocker on and escape. This is RUNG-6's clonal haircut, live.

**Two failure modes** = the actual recognition problem: **(A) false-kill** — normal cells that downregulate
the sensed HLA allele; **(B) escape** — antigen-low or subclonal HLA-retaining tumour cells.

## Honest framing
This is a **mechanistic circuit model**, not a measurement. The structure (safety = blocker reliability;
two failure routes; recognition→apoptosis coupling) is **parameter-robust**; the exact percentages are
parameter-dependent and reported with sweeps. *binding ≠ agonism* is the wet-lab residual. The model is built
so safety flows through the blocker — so the *contribution* is **quantifying the toxicity floor and showing it
is robust to the activator parameters**, plus the efficacy=clonal-LOH coupling, not the qualitative direction.

**The next measurable thing it points to:** normal-tissue HLA heterogeneity (the real safety constraint) —
measurable on the scRNA atlas (Colab). Selftest 10/10. Log in `runs/logs/`.
