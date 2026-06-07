# TODO — cancer-recog-apoptosis: next experiments

Living roadmap of what to run next. Every item is a falsifiable, in-silico test that runs on the M2 or free
Colab. Honest negatives are first-class. See `README.md` for the full hypothesis catalog and finished rungs.

**Legend:** `[ ]` open · `[x]` done · **P0** do next · **P1** soon · **P2** later/refinement · **FUT** future-safe (physics/delivery)
**Where:** 💻 M2 laptop (CPU) · ☁️ Colab CPU · ⚡ Colab GPU

---

## Done so far (recap — see README for detail)
- [x] RUNG 5–10b — surface logic gates bounded; single + combinatorial surface AND-NOT can't spare vital tissue
- [x] RUNG-11/12 — neoantigen route + AlphaFold pMHC structure (clean handles certified)
- [x] RUNG-13/14 — death wave validated + mechanism arena (quorum leads, ferroptosis/wave close)
- [x] RUNG-15/16/17 — atlas×mechanism map; clonal neoantigen burden (high-TMB broadly seedable); binding-axis immunogenicity (safety↔immunogenicity align)
- [x] RUNG-18/18b — **MHC window status**: genetically intact ~78% / dimmed ~18% / fully-dark ~4%; expression silencing adds ~2× (lung ~13% dark) → window broadly ON, genetics under-counts ~2×
- [x] RUNG-19 — evolutionary escape race: cure collapses at L=μ·N₀≈1; bystander cross-kill shifts curable size ~10×; clinical tumours NEED a resistance-agnostic 2nd mechanism
- [x] RUNG-20 (Boltz 2b) — Boltz confirms mutants present on MHC; can't discriminate mut-vs-WT by binding (saturated confidence) → discrimination is TCR-level, needs pose-RSA

---

## Recognition windows — normal cell vs cancer cell (the candidate handles)

Every row is a way to tell a cancer cell from a normal one — i.e. a potential recognition "window".
Single surface windows leak (proven). Inside windows are richer/more specific. Power = STACKING with AND-logic.

| # | Feature (plain) | Normal cell | Cancer cell | Where read | Status |
|---|---|---|---|---|---|
| 1 | **Mutations** (neoantigens = mutated protein pieces) | almost none | many, unique | inside → shown on MHC window | ✅ main handle (RUNG-11/16/17) |
| 2 | **MHC window** (display shelf) | on | mostly on, sometimes dark | surface | ✅ measured ~85–90% on (RUNG-18/18b) |
| 3 | **Surface marker proteins** (HER2/EGFR/EpCAM) | low | high — but also on organs | surface | ✅ tested → **leaks** (RUNG-15) |
| 4 | **Sugar coating** (glycans: Tn, sialyl-Tn) | normal | abnormal | surface | ❌ leak-test candidate |
| 5 | **Glucose use** (Warburg: guzzles sugar, spits acid) | calm, uses O₂ | fast, makes lactate | inside + environment | ❌ open |
| 6 | **Self-destruct machinery** (apoptosis) | intact | disabled (BCL-2 high) = weakness | inside | listed (BH3 mimetics) |
| 7 | **"Don't eat me" signal** (CD47) | low | high (hides from immune) | surface | ❌ leak-test candidate |
| 8 | **Chromosome count** (DNA amount, aneuploidy) | correct (46) | wrong, unstable | inside | ❌ open |
| 9 | **Telomerase** (immortality enzyme) | off in adult cells | switched on | inside | ❌ open |
| 10 | **Division rate** (Ki-67, replication stress) | rests | always dividing | inside | ❌ open |
| 11 | **Membrane lipid flip** (phosphatidylserine outside) | tucked inside | flipped out (stress) | surface | ❌ leak-test candidate |
| 12 | **Gene dependency** (synthetic-lethal, e.g. MTAP loss) | has both copies | lost one → addicted to partner | inside | ✅ tested (RUNG-14) |
| 13 | **Microenvironment** (pH, oxygen) | normal | acidic, hypoxic | around the cell | listed (Tier F) |
| 14 | **Physical** (stiffness, size, charge) | normal | softer/stiffer, depolarised | physical | FUT (oncotripsy/TTFields) |

---

## P0 — do next (runnable now)

- [ ] **NEW-WINDOW LEAK TEST** (the table → atlas). ☁️ Colab Census. For the top new candidate windows —
      **glycan-synthesis genes (#4)**, **CD47 / "don't eat me" (#7)**, **phosphatidylserine-flip / metabolic
      genes (#11/#5)** — measure normal-organ leakage exactly like RUNG-15 did for surface markers
      (worst-donor vital-tissue expression q_n). Output: which new windows are clean enough to STACK.
      *Reuse:* RUNG-15/34 Census loader (`scripts/40` census mode, `scripts/34` `find_leak_channels`).
      *Answers:* do any inside/surface windows beat the surface-marker leak ceiling? → the new clean handles.

- [ ] **COMBINATORIAL AND-LOGIC RECOGNITION** (stack windows for precision). 💻 M2. Take the windows that
      pass the leak test + the neoantigen window, model 2–3-input AND gates, compute per-patient
      addressability vs false-positive on vital tissue. *Answers:* does stacking windows push normal-cell
      false-positives toward zero while keeping tumour coverage? (the precision multiplier).

---

## P1 — soon

- [ ] **CAPSTONE SYNTHESIS** (the honest full story). 💻 M2 / writeup. Assemble the 3-stage chain end-to-end:
      recognition (neoantigen, surface dead) → binding (safety↔immunogenicity align) → apoptosis (wave) +
      window-status (RUNG-18/18b) + escape race (RUNG-19, needs cross-kill at clinical size). Ranked target
      shortlist + named wet-lab residuals. **Fully warranted now** — positive at every stage with stated bounds.

- [ ] **BOLTZ POSE-RSA REFINEMENT** (RUNG-20 done right). ⚡ Colab GPU. Boltz confirmed presentation but
      interface-pLDDT is the wrong ruler for TCR-level discrimination. Parse the Boltz CIF for the mutated
      residue's solvent exposure (RSA) + physicochemical change — the proper structural discrimination metric
      (matches RUNG-12's ESMFold approach, stronger model). *Reuse:* `scripts/37` `analyze_pdb` RSA logic.

- [ ] **RUNG-18b EXTEND to melanoma + bladder**. ⚡ Colab Census. They were absent in Census 2024-07-01 under
      those disease labels → try a newer Census version / specific melanoma & bladder scRNA datasets, or
      epithelial-fallback selection. Completes the route-cancer window-silencing picture (only lung is strong now).

- [ ] **PER-CANCER SURFACE AND-NOT PAIRS**. ☁️ Colab Census. RUNG-15 was pan-cancer pooled; a marker clean+high
      in ONE cancer (e.g. PSMA/prostate) could win per-cancer. Test 2-marker AND-NOT pairs per cancer type.

---

## P2 — catalog tiers not yet tested (from README hypothesis catalog)

- [ ] **C — Alternative death pathway addressability map**. ☁️ Census. Which tumours are wired for
      ferroptosis / pyroptosis / necroptosis / cuproptosis (the brake-free deaths that beat apoptosis-resistance
      & the escape race). Maps RUNG-14 `ferroptosis_wave` onto real per-cancer dependency.
- [ ] **D — BH3-mimetic dependence**. ☁️ Census. Per-tumour BCL-2/MCL-1 reliance → where lowering the
      apoptosis threshold (venetoclax-style) sensitises the recognition-gated wave.
- [ ] **E — Synthetic lethality deeper** (MTAP–PRMT5, ENO1–ENO2 collateral-deletion). ☁️ Census addressability.
- [ ] **F — Metabolic / microenvironment gates** (Warburg, pH, hypoxia as #5/#13 windows). 💻/☁️.
- [ ] **G — p53 refold** (AlphaFold ΔΔG: can a mutant p53 be folded back to active?). ⚡ GPU.
- [ ] **H — Replication-stress / mitotic-catastrophe** gating (#10 window). 💻/☁️.

---

## FUT — future-safe (physics / delivery; kept, not now)
- [ ] Oncotripsy (sound/ultrasound resonance at the cancer cell's natural frequency) — toy model in RUNG-14
- [ ] TTFields (electric-field disruption of division) — toy model in RUNG-14
- [ ] Photothermal / magnetic hyperthermia; nanorobotic nm-precision delivery; hybrids (wave + robotic injection)

---

## Standing rules (don't violate)
- Honest negatives are first-class; β / kill% / propensities are PROXIES, never verdicts; state the wet-lab residual.
- Validate (selftest) before spending compute. No subagent swarms — run the experiment.
- Commit validated results to the public repo; keep `memory/` (unpublished collaborator IP) PRIVATE/gitignored.
