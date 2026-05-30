# data/dr5_agonists — RUNG 2 calibration set

Frozen reference data for RUNG 2 (the DR5 clustering-geometry consistency check). Nothing here is
fit; it is the answer key + the structural inputs, committed before any model touches them.

| file | what |
|------|------|
| `dr5_agonist_ladder.csv` | **Table A** — n=9 intrinsic-potency ladder (y = `potency_rank`). The only y for the LOO / within-valency / permutation tests. |
| `tox_inversion_tableB.csv` | **Table B** — 3-row over-clustering/hepatotox endpoint. Illustrative only (tertiary, not credited). |
| `label_provenance.md` | blinding certification + leakage notes (adversary C1 fix). |
| `geometry_features.csv` | **generated** by `scripts/13_clustering_sim.py` (the clustering sim outputs g1–g6). Not hand-edited. |
| `geometry_from_boltz.csv` | **generated** on Colab by the Boltz-2 resolution probe (seed spread on the bivalents). Optional; literature geometry is used by default because Boltz-2 is unreliable on this receptor. |

## The honest one-liner

Across the **real** DR5-agonist ladder, **valency alone explains ~90% of potency** (Spearman ≈ 0.9,
computed live by `scripts/12`). RUNG 2 asks the narrow question of whether clustering *geometry* adds
anything *beyond* valency — and is wired so the default answer is the honest "no, it's a valency lookup
table," reachable only past a blinding + negative-control + oracle-noise gate. See
`../../runs/rung2_clustering/manifest_PREREG.yaml`.

## Sources (verify before any external use)

DR5Nb1 [Huet 2014, PMC4623017]; conatumumab/AMG655 [Kaplan-Lefko 2010]; lexatumumab/HGS-ETR2
[Plummer 2007; Shivange 2021]; drozitumab/Apomab [Adams 2008]; tigatuzumab/CS-1008 crosslink-independence
**CONTESTED** [Ichikawa 2001 vs comparative reviews]; TRAIL/dulanermin [DuBois 2023]; eftozanermin/ABBV-621
[Wei 2021]; INBRX-109/ozekibart [DuBois 2023 / Patnaik]; IGM-8444/aplitabart [Wang 2021]; TAS266
[Papadopoulos 2015]. **To re-verify before publication:** IGM-8444 INN = aplitabart (not "zamerovimab");
the cross-assay nature of the rank (different cell lines) is the dominant uncertainty.
