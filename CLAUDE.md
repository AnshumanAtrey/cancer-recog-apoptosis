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
