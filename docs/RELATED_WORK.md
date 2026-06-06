# cancer-recog-apoptosis — Related Work

The 16 papers that justify, ground, or constrain this project. Grouped by role.

---

## Group A — Biological foundation (proves the idea is real, not speculative)

### 1. Connexin-43 in Cancer: Above and Beyond Gap Junctions (2024)
Comprehensive review of Cx43's dual role — gap junction-mediated bystander apoptosis AND non-canonical C-terminal / hemichannel signaling.
**Use:** Cite when arguing the *mechanism* by which an engineered ligand could propagate apoptosis through cancer cells.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11674308/

### 2. Cascade-Targeting Apoptosis via TRAIL Bystander + Mitochondrial Photodamage (Nano Letters, 2025)
Nanoagent expresses TRAIL on tumor cells which induce DR5-mediated apoptosis in cancer neighbors. Direct in-vivo proof-of-concept.
**Use:** Closest existing therapeutic to Shriya's hypothesis. Cite as motivation; differentiate ours via *computational design* rather than wet-lab discovery.
https://pubs.acs.org/doi/10.1021/acs.nanolett.5c00878

### 3. Radiation-induced Bystander Signalling in Cancer Therapy (Nature Reviews Cancer)
Foundational review on cell-cell apoptosis propagation via gap junctions, NO, conditioned medium.
**Use:** Required reference for any paper claiming to engineer the bystander effect.
https://www.nature.com/articles/nrc2603

---

## Group B — Structure / affinity oracle stack

### 4. AlphaFold 3 — Accurate Structure Prediction of Biomolecular Interactions (Abramson et al., Nature 2024)
All-atom complex prediction including proteins, ligands, DNA, RNA, ions.
**Use:** Ground-truth benchmark only (server-only access). Cross-validate Boltz-2 outputs against it on selected complexes.
https://pubmed.ncbi.nlm.nih.gov/38718835/

### 5. Boltz-2 — Towards Accurate and Efficient Binding Affinity Prediction (bioRxiv, June 2025) ← PRIMARY ORACLE
First open model approaching FEP-quality affinity prediction at 1000× lower cost. MIT-licensed.
**Use:** Our oracle's binding-affinity signal. Both terms of the specificity reward (cancer binding +, normal binding −).
https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1

### 6. Protenix — ByteDance's Comprehensive AlphaFold3 Reproduction (bioRxiv, Jan 2025)
ByteDance PyTorch AF3 reproduction, Apache 2.0. v1 released Feb 2026 with AF3-level performance.
**Use:** Geopolitical-hedge oracle. Cross-validate critical predictions. Defensible "sovereign stack" framing for India/iDEX context.
https://www.biorxiv.org/content/10.1101/2025.01.08.631967v1

### 7. Technical Report of HelixFold3 (Baidu, arXiv Aug 2024)
PaddleHelix AF3 replication, accuracy comparable to AF3 on conventional ligands/nucleic acids/proteins.
**Use:** Second Chinese option. Backup if Boltz-2 / Protenix unavailable.
https://arxiv.org/abs/2408.16975

### 8. ESM3 — Simulating 500M Years of Evolution with a Language Model (Science 2024)
EvolutionaryScale's 98B multimodal protein LM reasoning over sequence + structure + function.
**Use:** Candidate base for the RL policy (alternative to Llama-3.2-3B). Or as the foldability/perplexity prior in the composite reward.
https://www.science.org/doi/10.1126/science.ads0018

---

## Group C — Design tools (the policy's action space)

### 9. De Novo Design of Protein Structure and Function with RFdiffusion (Watson et al., Nature 2023)
Diffusion-based backbone generation with sub-Angstrom control over binding geometry.
**Use:** Optional pre-conditioning step — generate plausible backbones for the policy to refine, or seed initial sequence templates.
https://www.nature.com/articles/s41586-023-06415-8

### 10. One-Shot Design of Functional Protein Binders with BindCraft (Nature, 2025)
End-to-end AF2-weight-leveraging pipeline, 10–100% experimental success rates without high-throughput screening.
**Use:** Competing pipeline. Argue our value-add: (a) RL loop, (b) specificity reward, (c) end-to-end simulation tier.
https://www.nature.com/articles/s41586-025-09429-6

### 11. PepINVENT — Generative Peptide Design Beyond Natural Amino Acids (2025)
RL-guided peptide generation with non-natural amino acid support.
**Use:** Closest "RL for peptide design" prior art. Reference for our methodology section. PharmaRL-of-peptides analog.
https://pmc.ncbi.nlm.nih.gov/articles/PMC12002334/

---

## Group D — ⚠ Closest prior art (read first; differentiate explicitly)

### 12. De Novo Protein Design Enables Targeting of Intractable Oncogenic Interfaces (Baker Lab, bioRxiv Oct 2025)
RFdiffusion + ProteinMPNN binders against "undruggable" oncogenic protein-protein interfaces.
**Differentiation we own:** (a) RL loop they don't have; (b) bystander/cell-cell framing they don't have; (c) explicit cancer-vs-normal specificity reward; (d) end-to-end biological simulation tier (PySB + PhysiCell + ADMET).
https://www.biorxiv.org/content/10.1101/2025.10.22.683953v1.full.pdf

### 13. Computationally Designed High-Specificity Inhibitors of BCL2 Pro-Survival Proteins
Designed three-helix bundles with pM-nM affinity and >300× specificity to individual BCL2 family members.
**Use:** Cite for "cancer-pathway-specific designed proteins work and specificity is achievable." Apoptosis-side molecular target evidence.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5127641/

### 14. SynNotch CAR-T Cell — When Synthetic Biology and Immunology Meet Again (Frontiers Immunology 2025)
2025 review of Wendell Lim / Kole Roybal logic-gated cell therapy with AND/OR/NOT gates.
**Use:** Conceptual ancestor — "engineering cell-cell communication for cancer." Position our work as the next-step ligand-design layer.
https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1545270/full

---

## Group E — Cancer-specific target discovery (Phase 0)

### 15. CellChat v2 for Cell-Cell Communication from scRNA-seq and Spatial Transcriptomics (Nature Protocols 2024)
Standard tool for inferring L-R communication networks from single-cell data; spatial support.
**Use:** Step 2 — discover cancer-restricted L-R pairs as our training targets. R package, callable via rpy2 or CSV export.
https://www.nature.com/articles/s41596-024-01045-4

### 16. PriorCCI — Interpretable DL for Key Ligand-Receptor Interactions Between Specific Cell Types (PMC 2025)
DL framework for prioritizing cancer-vs-other-cell L-R interactions.
**Use:** Solves the *specificity* problem at the target-discovery stage. Stacks on top of CellChat to rank-filter targets.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12345837/

---

## Group F — ⚠ CLOSEST prior art for the RUNG-4/5 SELECTIVITY line (surfaced by the highest-value-move review, 2026-06-01; READ FIRST, differentiate explicitly)

> Honest framing for this whole group: **none of these is something we beat on a new mechanism.** Every one
> does combinatorial / multi-modal target selection. What none of them does — and our single defensible
> delta — is certify safety with **worst-DONOR / worst-CELL / worst-HLA-allele, fail-closed, ancestry-
> stratified** discipline instead of pooled/averaged normal-tissue thresholds, and report the **per-patient
> "no-safe-window" addressability gap** as a structural impossibility surface. If we cannot state that delta
> in one sentence for a given paper, we have no contribution over it.

### 17. ImmunoVerse — pan-cancer atlas of therapeutic T-cell targets (bioRxiv 2025.01.22.634237 / PMC12265682, 2025)
Integrates 7,188 RNA-seq + 1,771 immunopeptidomes + 208 scRNA datasets vs 17,384 normal samples / 51 tissues; puts cell-SURFACE proteins (89) AND intracellular pMHC antigens (28,446) under **explicitly "identical filtering stringency"** across 21 tumour types; reports per-cancer-type patient coverage; flags canonical CAR targets (HER2/EPCAM/MUC1) that fail their safety filter.
**This is two of our four modalities, already unified.** Our delta: they certify safety on **"max of median TPM" / "average read count" — population POOLING**, the exact hack our RUNG-5 harness forbids. We re-audit under worst-case-per-donor/cell/allele and quantify the addressability gap they cannot see.
https://www.biorxiv.org/content/10.1101/2025.01.22.634237v1

### 18. SCAN-ACT — surface-target safety-vs-coverage atlas (bioRxiv 2025.01.18.633736, 2025)
Pan-cancer scRNA-based CAR-target nomination on the surface-selectivity-vs-safety axis.
**Use:** Cite as the surface-arm precedent; differentiate on worst-case safety, same as #17.
https://www.biorxiv.org/content/10.1101/2025.01.18.633736v1

### 19. TCGA_DEPMAP / GTEX_DEPMAP — translational dependency map (He et al., Nature Cancer 2024; PMID 39009815)
Maps DepMap CRISPR dependencies onto TCGA patients, prioritises synthetic lethalities by patient OUTCOME, and builds GTEX_DEPMAP — a normal-tissue tolerability layer for therapeutic-window selection.
**This is the addressable-fraction + normal-safety atlas the DepMap arm (Option D) proposed — already in Nature Cancer.** Our delta if we ever build that arm: cross-modality (CRISPR∩RNAi∩PRISM) agreement + a CLONALITY fail-closed guard (subclonality as a reward-hacking vector) + winner's-curse shrinkage on the essentiality gap. Reproduce their numbers first; report only the SHRINKAGE delta.
https://www.nature.com/articles/s43018-024-00789-y

### 20. Vinceti / Iorio et al. — clinically-informed map of cancer dependencies + target-prioritization framework (Cancer Cell 2023)
The DepMap target-prioritization framework (common-essential filtering, biomarker association).
**Use:** The standard we benchmark any genetic-dependency arm against; vanilla re-runs are not a contribution.
https://www.cell.com/cancer-cell/fulltext/S1535-6108(23)00219-2

### 21. SiCmiR Atlas — first single-cell mature-miRNA atlas (Cai et al., Advanced Science 2026; arXiv 2508.05692)
632 datasets, 9.36M cells, 726 cell types, cell-type-resolved miRNA-target networks, via a ~977-gene mRNA→miRNA predictor trained on paired TCGA data.
**KILLS the intracellular arm's headline.** Our last-turn hypothesis ("no single-cell miRNA atlas exists, so miRNA gates can't be safety-certified") is FALSE — this fills it, using the very mRNA→miRNA trick we'd have claimed. If we ever do the intracellular audit, we USE SiCmiR for the per-cell-type leak audit; we do not claim a blind spot.
https://onlinelibrary.wiley.com/doi/10.1002/advs.202508692

### 22. RADARS — reprogrammable ADAR sensors for endogenous-transcript-gated output (Kaseniit/Qian et al., Nature Biotechnology 2022)
RNA sensors that gate a payload (incl. caspase cargo) on an endogenous transcript, to allele resolution (TP53 R248W vs R248Q); NAR Cancer 2025 adds the patient-population editability feasibility audit.
**Use:** The highest-selectivity "cell senses its own cancer state" prior art; the faithful operationalisation of Shriya's concept is ALREADY 15 years deep (see #14 SynNotch lineage too). We do not claim self-recognition as novel.
https://www.nature.com/articles/s41587-022-01493-x

### 23. Combinatorial logic-gated target search — Dannenfelser 2020; Kwon & Kang 2023; MadHitter; LogiCAR 2025
The family of methods that search AND/OR/NOT antigen combinations for tumour-vs-healthy selectivity, several at single-cell resolution (LogiCAR 2025, Kwon 2023 audit cancer-vs-healthy gates at single-cell resolution).
**We are NOT first at combinatorial gate search.** Our contribution is the anti-reward-hacking "honesty harness" (worst-donor, dropout-fail-closed, fail-closed-vital, no-multiply, held-out-donor FDR + winner's-curse shrinkage) and the addressability-gap framing — never the search itself.
Dannenfelser 2020: https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30240-1 · LogiCAR 2025: https://www.biorxiv.org/content/10.1101/2025.02.10.637420v1

### 24. Death-durability / fractional-killing prior art — Spencer 2009; Roux 2015; Legewie 2006; oncogene-withdrawal recurrence (Sci Rep 2020); CAR-T escape-reservoir ABM (Sci Rep 2024, PMC11137006)
Cell-to-cell variability in apoptosis commitment (fractional killing), bistable-threshold recovery (= anastasis), oncogene-withdrawal recurrence ODEs, and an agent-based coverage×specificity×escape×toxicity CAR-T model.
**The "does death STICK?" question (Option F) is already modelled.** Any durability appendix must cite these and claim only the synthesis (necessary-conditions checklist on our own R1/R3 engine + anastasis-signature projection) with an uncertainty envelope, not a new result.
Spencer 2009: https://www.nature.com/articles/nature08012 · CAR-T ABM 2024: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11137006/

---

## If you read only 5 (in order) — RUNG-4/5 SELECTIVITY pivot

1. **#17 ImmunoVerse 2025** — the closest competitor; the whole contribution is "they pool normal tissue, we don't."
2. **#23 LogiCAR / Kwon 2023 / Dannenfelser 2020** — we are not first at combinatorial search; concede it, own the harness.
3. **#19 TCGA_DEPMAP (He 2024)** — the genetic-NOT-gate / addressable-fraction atlas already exists; only the shrinkage delta is ours.
4. **#21 SiCmiR Atlas 2026** — kills the intracellular "blind spot" headline; proof that we must literature-check before claiming novelty.
5. **#14 SynNotch / #22 RADARS** — the logic-gated and self-recognition mechanisms are mature; ours is the rigor layer, not the mechanism.

Honest read after the 2026-06-01 review: **the white space is NOT a new mechanism.** It is the worst-case-safety,
ancestry-stratified, FDR-and-winner's-curse-controlled re-audit of known targets, and the per-patient
addressability gap reported as a first-class negative. That is narrow, methodological, and defensible — and it
is the one thing the field's pooled atlases are structurally incentivised not to do.

(The older Group A–E framing above pertains to the earlier ligand-design / bystander line; retained for provenance.)
