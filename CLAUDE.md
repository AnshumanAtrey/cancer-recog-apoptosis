# cancer-recog-apoptosis — working rules

This project INVENTS something not yet in the literature. Online research and spawned subagents can only
summarise already-published work — they cannot judge whether an un-invented thing will work. The data is
the oracle.

1. **Groundbreaking empirical testing > online research.** Default to building a runnable script that lets
   the DATA / atlas / simulation answer the question. A few targeted web lookups are fine for grounding;
   never substitute them for a real run.
2. **No subagent swarms.** Do NOT spawn Workflow/Agent fan-outs (20–30 subagents) to "verify feasibility"
   or do SOTA research — they burn tokens guessing at things that exist nowhere. Run the experiment instead.
3. **Honest negatives are first-class.** Never overclaim. Every predicted result is a HYPOTHESIS with a
   stated wet-lab residual; report threshold-sensitive numbers as ranges; let the atlas say no when it says no.
4. **Census API (recurring fix):** `cellxgene_census.get_anndata(...)` — use `obs_column_names=[...]` (and `var_column_names=[...]`), NOT the deprecated `column_names={"obs": [...]}` (FutureWarning). Same for `obs_value_filter`/`var_value_filter`.
5. **Audit a result BEFORE building on it; re-audit old runs/data too.** Every result — especially a STRONG one, and ESPECIALLY a NEGATIVE that closes a direction — is a claim about the *method* until verified. A surprising result is more often a bug than a discovery. **Cost asymmetry:** a false positive gets caught downstream (you build on it, it fails); a false negative is SILENT — it kills a real opportunity you never revisit → verify negatives *at least as hard* as positives. When a run finishes (or when revisiting an old run/dataset), before propagating any conclusion or starting the next rung:
   - **Biological sanity first (cheapest tripwire):** does the number contradict known biology? *(RUNG-23 v1's 62% proliferation "leak" into POST-MITOTIC cardiomyocytes/neurons — which cannot divide — was impossible → a method bug, not a finding.)* Can't-be-true ⇒ a flag, not a result.
   - **Were the disciplines applied?** (checklist below). Most of our bugs were *"already fixed in another rung, not carried over"* — so apply the checklist **proactively when WRITING** a new atlas script, not only when auditing after.
   - **Could a better method flip it?** If plausibly yes ⇒ HARDEN + RE-RUN before concluding, and **downgrade any claim already written** until the hardened re-run confirms (never leave an overclaim standing).
   **Method-disciplines checklist (atlas / scRNA / expression runs):**
   (a) **DEPTH-GATE** (RUNG-18b) — score only well-sequenced cells (housekeeping panel ≥ k detected) so dropout doesn't deflate the signal / fake a "dead" reading.
   (b) **ROBUST DONOR STAT** (RUNG-8) — per-donor distribution (median / p90), NEVER the single worst donor (one outlier fakes a leak); report which donor/cell-type drives any extreme.
   (c) **DATASET-MEASURING FILTER** (RUNG-8) — Census returns 0 for "not assayed" too; only count datasets that detect the gene somewhere.
   (d) **SELFTEST the math first** with an artifact-resistance check; report the tumour-vs-control CONTRAST, not dropout-sensitive absolutes.
