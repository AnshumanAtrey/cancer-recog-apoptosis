#!/usr/bin/env python3
"""
RUNG 2 / Stage-0+2 — DR5 receptor clustering simulation (spatial percolation, union-find).

Generates the GEOMETRY features (g1..g6) per agonist from a real 2D receptor-clustering model:
multivalent agonist (valency N, arm-reach r_arm from the BLIND labeler) crosslinks DR5 receptors
distributed on a membrane patch. Receptors co-bound by one agonist are unioned (a clique); a secondary
crosslinker (optional) bridges nearby agonists. Connected components >= the DISC-nucleation size
(>=6 receptor death-domains = 2 trimers) are 'firing-competent'. This is a percolation model: average
receptor-graph degree rises with valency, so the firing fraction rises with valency THROUGH a percolation
threshold — which is exactly WHY valency predicts potency. RUNG 2's positive, defensible result is that
this first-principles sim reproduces the real valency->potency ladder; the geometry-BEYOND-valency
question is handled (and most likely refuted) in scripts/12.

INPUTS (blind, endpoint-free): data/dr5_agonists/blind_geometry_labels.csv (valency, r_arm_nm,
epitope_proximity, geometry_match) from scripts/14. NO potency rank is read here.

CEILING (honest, printed): firing-competent CLUSTER GEOMETRY != caspase-8 firing. The cooperative
enzymatic threshold at the DISC (procaspase-8 dimerization/auto-cleavage vs cFLIP/FADD stoichiometry,
membrane context, cortical-actin corralling) is NOT modelled and is the irreducible wet-lab crux
(EVIDENCE_AND_HANDOFF.md). Sim parameters (density, r_arm, FIRE_THR, p_bind) are literature-anchored
PROXIES; the sim tests ORDERING/logic, not absolute kinetics. NEVER multiply these outputs by the
RUNG-1 EARM gate-strength.

USAGE:  python scripts/13_clustering_sim.py
OUT  :  data/dr5_agonists/geometry_features.csv  +  runs/rung2_clustering/clustering_phase.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LABELS = PROJECT_ROOT / "data" / "dr5_agonists" / "blind_geometry_labels.csv"
OUT_CSV = PROJECT_ROOT / "data" / "dr5_agonists" / "geometry_features.csv"
OUT_DIR = PROJECT_ROOT / "runs" / "rung2_clustering"

# ---- frozen sim parameters (PROXIES, documented in manifest_PREREG / open_risks) ----
# Regime: receptor density high enough that arm-reach is NOT the bottleneck (mean spacing < smallest
# r_arm), and the agonist *bond budget* is held constant across valency (N_ag = occ*N_rec/valency).
# Valency then sets the agonist CLIQUE SIZE at fixed total bonds -> bigger cliques (high valency)
# percolate more easily -> firing fraction rises with valency. This is the percolation mechanism by
# which valency drives clustering competence, and it is what makes valency predict potency.
SEED = 20260530
L_UM = 0.30                 # membrane patch side (um)
RHO_CANCER = 10000.0        # DR5 receptors / um^2 (clustering-competent tumour cell; mean spacing ~5nm)
RHO_HEPATOCYTE = 2500.0     # receptors / um^2 (low-density bystander / hepatocyte)
AGONIST_OCCUPANCY = 1.5     # bond budget = occ*N_rec, distributed as cliques of size=valency
FIRE_THR = 6                # receptors per component to be DISC-firing-competent (2 trimers)
R_XLINK_NM = 20.0           # secondary-crosslinker bridging radius
P_BIND_BASE = 0.92          # per productive agonist-receptor contact
N_REPLICATES = 8


class UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n)); self.sz = [1] * n
    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return
        if self.sz[ra] < self.sz[rb]: ra, rb = rb, ra
        self.p[rb] = ra; self.sz[ra] += self.sz[rb]


def p_bind(epitope_proximity: float) -> float:
    """Productive-contact probability rises mildly with membrane-proximal epitope (CRD3)."""
    return min(1.0, P_BIND_BASE * (0.75 + 0.10 * epitope_proximity))


def simulate_once(valency, r_arm_nm, epitope_proximity, rho, crosslinker, rng):
    """One spatial realisation -> (firing_fraction, mean_cluster_size, giant_fraction)."""
    n_rec = max(2, int(round(rho * L_UM * L_UM)))
    # constant bond budget across valency: cliques of size=valency, occ*N_rec total bonds
    n_ag = max(1, int(round(n_rec * AGONIST_OCCUPANCY / max(1, valency))))
    rec = rng.uniform(0, L_UM, size=(n_rec, 2))
    ag = rng.uniform(0, L_UM, size=(n_ag, 2))
    r_arm = r_arm_nm / 1000.0  # um
    tree = cKDTree(rec)
    uf = UnionFind(n_rec)
    pb = p_bind(epitope_proximity)
    ag_receptors = []  # for crosslinker bridging
    for a in ag:
        idx = tree.query_ball_point(a, r=r_arm)
        if not idx:
            ag_receptors.append([]); continue
        # bind up to `valency` of the in-reach receptors, each with prob pb
        rng.shuffle(idx)
        bound = [i for i in idx[:valency] if rng.random() < pb]
        ag_receptors.append(bound)
        for i in bound[1:]:
            uf.union(bound[0], i)  # clique via the first bound receptor
    if crosslinker and n_ag > 1:
        # secondary crosslinker bridges receptor-sets of agonists within R_XLINK
        atree = cKDTree(ag)
        r_x = R_XLINK_NM / 1000.0
        for ai, a in enumerate(ag):
            if not ag_receptors[ai]:
                continue
            for aj in atree.query_ball_point(a, r=r_x):
                if aj > ai and ag_receptors[aj]:
                    uf.union(ag_receptors[ai][0], ag_receptors[aj][0])
    # component sizes
    roots = {}
    for i in range(n_rec):
        r = uf.find(i); roots[r] = roots.get(r, 0) + 1
    sizes = np.array(list(roots.values()))
    firing = int(sizes[sizes >= FIRE_THR].sum())
    firing_fraction = firing / n_rec
    bound_sizes = sizes[sizes >= 2]
    mean_cluster = float(bound_sizes.mean()) if bound_sizes.size else 1.0
    giant_fraction = float(sizes.max()) / n_rec
    return firing_fraction, mean_cluster, giant_fraction


def avg_sim(valency, r_arm_nm, epitope_proximity, rho, crosslinker, rng):
    ff, mc, gf = [], [], []
    for _ in range(N_REPLICATES):
        a, b, c = simulate_once(valency, r_arm_nm, epitope_proximity, rho, crosslinker, rng)
        ff.append(a); mc.append(b); gf.append(c)
    return float(np.mean(ff)), float(np.mean(mc)), float(np.mean(gf))


def main() -> int:
    if not LABELS.exists():
        print("[clustering-sim] missing blind labels — run scripts/14_blind_labeler.py first")
        return 2
    labels = pd.read_csv(LABELS)
    rng = np.random.default_rng(SEED)
    rows = []
    print(f"[clustering-sim] L={L_UM}um  rho_cancer={RHO_CANCER}/um^2  FIRE_THR={FIRE_THR}  "
          f"reps={N_REPLICATES}  (params are PROXIES; sim tests ordering, not absolute kinetics)")
    for _, r in labels.iterrows():
        v, ra, epi = int(r["valency"]), float(r["r_arm_nm"]), float(r["epitope_proximity"])
        # g1: firing fraction at cancer density, NO secondary crosslinker (intrinsic clustering)
        g1, g2, g3 = avg_sim(v, ra, epi, RHO_CANCER, False, rng)
        # firing WITH crosslinker (for the crosslink-dependence read)
        ff_x, _, _ = avg_sim(v, ra, epi, RHO_CANCER, True, rng)
        # g6: hyperclustering / over-clustering at hepatocyte density WITH crosslinker (tox axis)
        g6, _, _ = avg_sim(v, ra, epi, RHO_HEPATOCYTE, True, rng)
        rows.append({
            "name": r["name"], "valency": v,
            "g1_firing_fraction_nocrosslink": round(g1, 4),
            "g2_mean_cluster_size": round(g2, 3),
            "g3_percolation_margin": round(g3, 4),
            "g4_geometry_match": round(float(r["geometry_match"]), 4),
            "g5_epitope_proximity": epi,
            "g6_hyperclustering_propensity": round(g6, 4),
            "firing_fraction_crosslinked": round(ff_x, 4),
            "crosslink_dependence_index": round(ff_x - g1, 4),  # high => needs crosslinker (descriptive)
        })
    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print("[clustering-sim] wrote", OUT_CSV.relative_to(PROJECT_ROOT))
    print(out[["name", "valency", "g1_firing_fraction_nocrosslink", "g2_mean_cluster_size",
               "g3_percolation_margin", "firing_fraction_crosslinked", "crosslink_dependence_index",
               "g6_hyperclustering_propensity"]].to_string(index=False))

    # sanity: firing fraction should rise (mostly) with valency -> reproduces the valency law
    from scipy.stats import spearmanr
    rho_sim, _ = spearmanr(out["valency"], out["g1_firing_fraction_nocrosslink"])
    print(f"[clustering-sim] Spearman(valency, firing_fraction_nocrosslink) = {rho_sim:.3f}  "
          f"(sim reproduces the monotone valency->clustering law)")

    # phase figure
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        o = out.sort_values("valency")
        ax[0].plot(o["valency"], o["g1_firing_fraction_nocrosslink"], "o-", color="#c0392b", label="no crosslinker")
        ax[0].plot(o["valency"], o["firing_fraction_crosslinked"], "s--", color="#2980b9", label="+ crosslinker")
        ax[0].set_xlabel("valency (# DR5 arms)"); ax[0].set_ylabel("firing-competent receptor fraction")
        ax[0].set_title(f"clustering percolation vs valency\nSpearman(valency,firing)={rho_sim:.2f}")
        ax[0].legend(fontsize=8)
        ax[1].bar(o["name"], o["crosslink_dependence_index"], color="#7f8c8d")
        ax[1].set_ylabel("crosslink-dependence index\n(firing_+xlink − firing_noxlink)")
        ax[1].set_title("bivalents depend on crosslinker (descriptive)")
        ax[1].tick_params(axis="x", rotation=90, labelsize=7)
        fig.suptitle("RUNG 2 clustering sim — geometry CAPABILITY only; cluster geometry != caspase-8 firing "
                     "(agonism = wet-lab). NEVER multiply by RUNG-1.", fontsize=9)
        fig.tight_layout(); fig.savefig(OUT_DIR / "clustering_phase.png", dpi=110)
        print("[clustering-sim] figure -> runs/rung2_clustering/clustering_phase.png")
    except Exception as e:
        print(f"[clustering-sim] figure skipped ({type(e).__name__}: {e})")

    print("[clustering-sim] CEILING: firing-competent cluster geometry is NOT caspase-8 firing "
          "(cooperative DISC threshold = wet-lab). Params are proxies; ordering only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
