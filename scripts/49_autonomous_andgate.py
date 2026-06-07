#!/usr/bin/env python3
"""
RUNG 23 — the AUTONOMOUS intracellular AND-gate: a synthetic gene circuit that senses TWO intracellular
cancer signals (signal A AND signal B) and self-destructs — NO MHC, NO immune system. (Census; CPU.)

THE IDEA (Shriya's ORIGINAL concept, the un-crowded route)
----------------------------------------------------------
The immune route (T-cells / vaccine) dies on the MHC-DARK core (RUNG-18: ~4% permanently dark, the blind
spot). That core has NO MHC -> nothing to present -> T-cells blind. The escape: a recognition gate that lives
INSIDE the cell and needs no MHC at all — a synthetic gene circuit "signal A AND signal B -> apoptosis". AND
is the RIGHT operator here (specificity: a normal cell with only one signal must NOT fire), and it is
BUILDABLE inside a cell (synthetic TFs / toehold switches / split-effectors) — unlike a vaccine, which can
only do OR. This run asks the make-or-break question on REAL data:

    Is there a PAIR of intracellular transcriptional programs (A,B) such that (A high AND B high) fires in
    TUMOUR cells but in ~ZERO vital normal cells (worst-donor), where each single program LEAKS?

If yes -> an MHC-independent, AND-specific, autonomous self-destruct for exactly the cells nothing else
reaches. If no -> the honest negative (intracellular signals leak into regenerating vital tissue too, and AND
can't clean it) — which itself bounds the autonomous route.

PROGRAMS (intracellular signals; each a small marker panel scored by detection fraction)
---------------------------------------------------------------------------------------
PROLIF (cell division), MYC (oncogene drive), E2F (cell-cycle), GLYCOLYSIS (Warburg), HYPOXIA, WNT
(driver pathway), STEMNESS, TELOMERASE. Single proliferation/oncogene signals are KNOWN to leak into
regenerating tissue (gut crypt, bone marrow, skin) — the AND-gate's job is to clean that up.

METHOD (safety discipline carried from RUNG-5/8)
-----------------------------------------------
Per cell, a program "fires" if it detects >= FIRE_FRAC of its panel genes. Gate fire = AND/OR/single of two
programs. Per (vital cell type, donor): leak = fraction firing the gate; headline = WORST-DONOR leak (the
RUNG-5/8 worst-donor safety bar). Coverage = tumour-cell gate-fire fraction. A CLEAN autonomous gate = high
coverage AND worst-donor vital leak ~0. Only DATASETS that measured the panel are counted (Census 0 ==
not-assayed too; the RUNG-8 correction).

HONEST CEILING
--------------
mRNA != protein/activity (a transcriptional-program score is a PROXY for the real intracellular state a
circuit would sense); scRNA dropout deflates fire (-> coverage is a lower bound, leak too — read the
tumour-vs-vital CONTRAST). "Buildable" assumes a sensor exists for each program (synthetic-biology residual;
the binding/sensing step is wet-lab). Delivery of the circuit to every cell is the big residual (TODO). This
BOUNDS whether a clean intracellular AND-pair exists; it is not a built circuit.

USAGE
  python scripts/49_autonomous_andgate.py selftest    # synthetic, validates gate + leak/coverage + worst-donor
  RUNG23_CACHE=/content/drive/MyDrive/cancer-recon python scripts/49_autonomous_andgate.py run   # Colab Census
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import threading
import time
from itertools import combinations
from pathlib import Path

import numpy as np

_T0 = time.monotonic()
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung23_autonomous"
RESULT_JSON = OUT_DIR / "rung23_autonomous.json"
FIGURE_PNG = OUT_DIR / "rung23_autonomous.png"

# intracellular transcriptional programs (marker panels). Kept small + canonical.
PROGRAMS = {
    "PROLIF": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "BIRC5"],
    "MYC": ["MYC", "NPM1", "ODC1", "NCL", "SRM"],
    "E2F": ["E2F1", "MCM2", "MCM6", "CDC6", "CCNE1"],
    "GLYCOLYSIS": ["SLC2A1", "LDHA", "HK2", "PKM", "PGK1", "ENO1"],
    "HYPOXIA": ["VEGFA", "CA9", "SLC2A1", "NDRG1", "BNIP3"],
    "WNT": ["AXIN2", "LGR5", "ASCL2", "NKD1", "RNF43"],
    "STEMNESS": ["SOX2", "PROM1", "ALDH1A1", "POU5F1"],
    "TELOMERASE": ["TERT", "TERC", "DKC1"],
}
ALL_GENES = sorted({g for panel in PROGRAMS.values() for g in panel})
FIRE_FRAC = 0.5                 # a program "fires" in a cell if >= this fraction of its panel is detected
DET = 1                         # gene detected if UMI >= DET
PER_DONOR_CAP = 600
MIN_CELLS_PER_DONOR = 30
LEAK_SAFE = 0.01                # a gate is "vital-safe" if worst-donor leak <= this
COVER_MIN = 0.30                # ... and "effective" if tumour coverage >= this


def log(msg):
    print(f"[+{time.monotonic()-_T0:7.1f}s] [rung23] {msg}", flush=True)


# ---------------------------------------------------------------------------
#  CORE LOGIC (selftestable, no Census) — program fire, gate fire, leak/coverage
# ---------------------------------------------------------------------------
def program_fire(counts: np.ndarray, gene_order: list[str]) -> dict:
    """counts [N x G] in gene_order. Returns {program: bool[N]} where a program fires if >=FIRE_FRAC detected."""
    idx = {g: j for j, g in enumerate(gene_order)}
    out = {}
    for prog, panel in PROGRAMS.items():
        cols = [idx[g] for g in panel if g in idx]
        if not cols:
            out[prog] = np.zeros(counts.shape[0], bool)
            continue
        det = (counts[:, cols] >= DET).sum(axis=1)
        out[prog] = det >= np.ceil(FIRE_FRAC * len(cols))
    return out


def gate_fire(pf: dict, a: str, b: str | None, op: str) -> np.ndarray:
    if op == "single":
        return pf[a]
    if op == "AND":
        return pf[a] & pf[b]
    if op == "OR":
        return pf[a] | pf[b]
    raise ValueError(op)


def _measuring_datasets(counts, gene_order, dataset):
    """A dataset measured the panel iff it detects >=1 of ALL_GENES somewhere (Census 0 == not-assayed too)."""
    any_det = (counts >= DET).any(axis=1)
    return {d for d in set(dataset.tolist()) if any_det[dataset == d].any()}


def evaluate_all(counts, gene_order, label_vital, donor, dataset, is_tumour):
    """Evaluate every single program + every AND/OR pair. Returns ranked clean gates."""
    pf = program_fire(counts, gene_order)
    measuring = _measuring_datasets(counts, gene_order, dataset)
    meas = np.array([d in measuring for d in dataset], bool)
    is_vital = (label_vital != None) & meas                # vital normal cells from measuring datasets
    tum = is_tumour & meas
    progs = list(PROGRAMS)
    results = []

    def metrics(fire):
        cov = float(fire[tum].mean()) if tum.any() else 0.0
        worst = 0.0
        for d in np.unique(donor[is_vital]):
            m = is_vital & (donor == d)
            if m.sum() >= MIN_CELLS_PER_DONOR:
                worst = max(worst, float(fire[m].mean()))
        return round(cov, 4), round(worst, 4)

    for p in progs:
        cov, leak = metrics(pf[p])
        results.append({"gate": p, "op": "single", "coverage": cov, "worst_leak": leak,
                        "safe_effective": leak <= LEAK_SAFE and cov >= COVER_MIN})
    for a, b in combinations(progs, 2):
        for op in ("AND", "OR"):
            cov, leak = metrics(gate_fire(pf, a, b, op))
            results.append({"gate": f"{a} {op} {b}", "op": op, "a": a, "b": b,
                            "coverage": cov, "worst_leak": leak,
                            "safe_effective": leak <= LEAK_SAFE and cov >= COVER_MIN})
    results.sort(key=lambda r: (-int(r["safe_effective"]), r["worst_leak"], -r["coverage"]))
    return results, {"n_datasets_measuring": len(measuring), "n_vital_cells": int(is_vital.sum()),
                     "n_tumour_cells": int(tum.sum())}


# ---------------------------------------------------------------------------
#  Census fetch (reuses the RUNG-8/29 vital pull + RUNG-5 tumour pattern)
# ---------------------------------------------------------------------------
class Heartbeat:
    def __init__(self, interval=20):
        self.interval, self.label, self._stop = interval, "starting", False

    def set(self, label):
        self.label = label; log(label)

    def _run(self):
        while not self._stop:
            for _ in range(self.interval * 2):
                if self._stop:
                    return
                time.sleep(0.5)
            if not self._stop:
                print(f"[+{time.monotonic()-_T0:7.1f}s] [heartbeat] {self.label}", flush=True)

    def start(self):
        threading.Thread(target=self._run, daemon=True).start(); return self

    def stop(self):
        self._stop = True


def _load(name, mod):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / mod)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def _codes(col):
    import pyarrow as pa
    col = col.combine_chunks()
    if not pa.types.is_dictionary(col.type):
        col = col.dictionary_encode()
    return col.indices.to_numpy(zero_copy_only=False), [str(x) for x in col.dictionary.to_pylist()]


def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    HB = Heartbeat().start()
    import cellxgene_census
    d4 = _load("d4", "17_logicgate_data.py")
    HB.set(f"opening Census {d4.CENSUS_VERSION} ...")
    census = cellxgene_census.open_soma(census_version=d4.CENSUS_VERSION)
    exp = census["census_data"]["homo_sapiens"]

    def fetch(vf, cap_label):
        tbl = exp.obs.read(value_filter=vf, column_names=["soma_joinid", "cell_type", "donor_id", "dataset_id"]).concat()
        if tbl.num_rows == 0:
            return None
        jid = tbl.column("soma_joinid").to_numpy()
        ct_codes, ct_vals = _codes(tbl.column("cell_type"))
        ds_codes, ds_vals = _codes(tbl.column("dataset_id"))
        dn_codes, dn_vals = _codes(tbl.column("donor_id"))
        gid = ds_codes.astype(np.int64) * (int(dn_codes.max()) + 1) + dn_codes
        keep = []
        for g in np.unique(gid):
            keep.extend(np.where(gid == g)[0][:PER_DONOR_CAP].tolist())
        keep = np.sort(np.array(keep, np.int64))
        sel = jid[keep]
        HB.set(f"{cap_label}: materialising {len(ALL_GENES)} genes for {len(sel):,} cells ...")
        ad = cellxgene_census.get_anndata(census, organism="Homo sapiens", obs_coords=sel.tolist(),
                                          var_value_filter=f"feature_name in {ALL_GENES}",
                                          column_names={"obs": ["cell_type", "donor_id", "dataset_id"]})
        vn = list(ad.var["feature_name"]) if "feature_name" in ad.var else list(ad.var_names)
        X = np.asarray(ad.X.todense() if hasattr(ad.X, "todense") else ad.X)
        counts = np.zeros((X.shape[0], len(ALL_GENES)), np.int32)
        for j, g in enumerate(ALL_GENES):
            if g in vn:
                counts[:, j] = X[:, vn.index(g)].astype(np.int32)
        ct = ad.obs["cell_type"].astype(str).to_numpy()
        donor = np.array([f"{a}::{b}" for a, b in zip(ad.obs["dataset_id"].astype(str), ad.obs["donor_id"].astype(str))], object)
        dataset = np.array([str(a) for a in ad.obs["dataset_id"].astype(str)], object)
        return counts, ct, donor, dataset

    # vital normal cells across tissues
    all_c, all_lab, all_don, all_ds, all_tum = [], [], [], [], []
    for tissue in d4.NORMAL_TISSUES:
        r = fetch(f"is_primary_data == True and disease == 'normal' and tissue_general == '{tissue}'", f"normal {tissue}")
        if r is None:
            continue
        counts, ct, donor, dataset = r
        lab = np.array([next((v for k, v in d4.VITAL_AUDIT.items() if k in c.lower()), None) for c in ct], object)
        all_c.append(counts); all_lab.append(lab); all_don.append(donor); all_ds.append(dataset)
        all_tum.append(np.zeros(len(ct), bool))
    # tumour cells (malignant)
    r = fetch(f"is_primary_data == True and cell_type in ['malignant cell','neoplastic cell']", "malignant")
    if r is not None:
        counts, ct, donor, dataset = r
        all_c.append(counts); all_lab.append(np.array([None] * len(ct), object))
        all_don.append(donor); all_ds.append(dataset); all_tum.append(np.ones(len(ct), bool))
    HB.stop()

    counts = np.vstack(all_c); label = np.concatenate(all_lab); donor = np.concatenate(all_don)
    dataset = np.concatenate(all_ds); is_tum = np.concatenate(all_tum)
    log(f"total {counts.shape[0]:,} cells; vital {int((label!=None).sum()):,}; tumour {int(is_tum.sum()):,}")

    results, meta = evaluate_all(counts, ALL_GENES, label, donor, dataset, is_tum)
    clean = [r for r in results if r["safe_effective"]]
    best_and = next((r for r in results if r["op"] == "AND" and r["safe_effective"]), None)

    out = {
        "tag": "rung23_autonomous_andgate",
        "question": "Is there an intracellular AND-pair (A high AND B high) firing in tumour cells but ~zero "
                    "vital normal cells (worst-donor) -> an MHC-independent autonomous self-destruct gate?",
        "programs": PROGRAMS, "fire_frac": FIRE_FRAC, "leak_safe": LEAK_SAFE, "cover_min": COVER_MIN,
        "meta": meta,
        "all_gates_ranked": results[:40],
        "safe_effective_gates": clean,
        "best_AND_gate": best_and,
        "HEADLINE": (f"{len(clean)} safe&effective gates (worst-donor vital leak ≤{LEAK_SAFE}, tumour coverage "
                     f"≥{COVER_MIN}). Best AND: {best_and['gate'] if best_and else 'NONE'} "
                     f"(cov {best_and['coverage'] if best_and else '-'}, leak {best_and['worst_leak'] if best_and else '-'}). "
                     f"If an AND beats its single programs on leak -> autonomous MHC-free gate is buildable for the dark core."),
        "CEILING": "mRNA != protein/activity (program score is a PROXY for the state a circuit senses); scRNA "
                   "dropout deflates fire (read the tumour-vs-vital CONTRAST, not absolutes); 'buildable' "
                   "assumes a sensor per program (synthetic-biology + delivery = wet-lab residual). BOUNDS "
                   "whether a clean intracellular AND-pair exists; not a built circuit.",
    }
    RESULT_JSON.write_text(json.dumps(out, indent=2))
    log(f"wrote {RESULT_JSON}")
    log(f"safe&effective: {[r['gate'] for r in clean][:8]}")
    _make_figure(results)
    return 0


def _make_figure(results):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log(f"matplotlib unavailable ({e})"); return
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for r in results:
        c = {"single": "#888", "AND": "#3F7D54", "OR": "#E0A040"}[r["op"]]
        ax.scatter(r["worst_leak"] * 100, r["coverage"] * 100, color=c, s=28, alpha=0.8)
    ax.axvline(LEAK_SAFE * 100, ls="--", color="#B23A2E", label=f"safe leak ≤{LEAK_SAFE*100:.0f}%")
    ax.axhline(COVER_MIN * 100, ls=":", color="grey", label=f"effective cov ≥{COVER_MIN*100:.0f}%")
    for op, c in [("single", "#888"), ("AND", "#3F7D54"), ("OR", "#E0A040")]:
        ax.scatter([], [], color=c, label=op)
    ax.set_xlabel("worst-donor vital-tissue leak (%)  ← safer"); ax.set_ylabel("tumour coverage (%)")
    ax.set_title("RUNG-23: autonomous intracellular gates — specificity × coverage\n(top-left green = clean buildable AND-gate)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_xscale("symlog", linthresh=1)
    fig.tight_layout(); fig.savefig(FIGURE_PNG, dpi=130)
    log(f"wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    go = ALL_GENES

    def cell(fire_progs):
        """build a cell expressing exactly the panels of the named programs."""
        c = np.zeros(len(go), np.int32)
        idx = {g: j for j, g in enumerate(go)}
        for p in fire_progs:
            for g in PROGRAMS[p]:
                c[idx[g]] = 5
        return c

    # program_fire detects exactly the intended programs
    cm = np.array([cell(["PROLIF"]), cell(["PROLIF", "GLYCOLYSIS"]), cell([])])
    pf = program_fire(cm, go)
    check("PROLIF fires where its panel present", pf["PROLIF"][0] and pf["PROLIF"][1] and not pf["PROLIF"][2])
    check("GLYCOLYSIS fires only in cell 1", (not pf["GLYCOLYSIS"][0]) and pf["GLYCOLYSIS"][1])

    # gate logic
    check("AND requires both", gate_fire(pf, "PROLIF", "GLYCOLYSIS", "AND").tolist() == [False, True, False])
    check("OR requires either", gate_fire(pf, "PROLIF", "GLYCOLYSIS", "OR").tolist() == [True, True, False])
    check("single = program", gate_fire(pf, "PROLIF", None, "single").tolist() == [True, True, False])

    # leak/coverage + worst-donor + AND-cleans-leak scenario:
    # vital donor V1: PROLIF only (proliferating crypt) -> single PROLIF leaks, but AND(PROLIF,GLYCOLYSIS) does NOT
    # tumour: PROLIF AND GLYCOLYSIS -> AND covers tumour, spares vital
    rows_c, rows_lab, rows_don, rows_ds, rows_tum = [], [], [], [], []
    for _ in range(40):
        rows_c.append(cell(["PROLIF"])); rows_lab.append("cardiomyocyte"); rows_don.append("dsV::V1"); rows_ds.append("dsV"); rows_tum.append(False)
    for _ in range(40):
        rows_c.append(cell(["PROLIF", "GLYCOLYSIS"])); rows_lab.append(None); rows_don.append("dsT::T1"); rows_ds.append("dsT"); rows_tum.append(True)
    C = np.array(rows_c); L = np.array(rows_lab, object); D = np.array(rows_don, object)
    DS = np.array(rows_ds, object); TUM = np.array(rows_tum, bool)
    res, meta = evaluate_all(C, go, L, D, DS, TUM)
    g = {r["gate"]: r for r in res}
    check("single PROLIF LEAKS into vital (worst_leak==1.0)", abs(g["PROLIF"]["worst_leak"] - 1.0) < 1e-9)
    check("AND(PROLIF,GLYCOLYSIS) spares vital (leak==0)", abs(g["PROLIF AND GLYCOLYSIS"]["worst_leak"] - 0.0) < 1e-9)
    check("AND covers tumour (coverage==1.0)", abs(g["PROLIF AND GLYCOLYSIS"]["coverage"] - 1.0) < 1e-9)
    check("AND is safe&effective, single PROLIF is NOT",
          g["PROLIF AND GLYCOLYSIS"]["safe_effective"] and not g["PROLIF"]["safe_effective"])

    # dataset-measuring exclusion: a dataset detecting nothing is dropped
    Cz = np.vstack([C, np.zeros((20, len(go)), np.int32)])
    Lz = np.concatenate([L, np.array([None] * 20, object)])
    Dz = np.concatenate([D, np.array(["dsZ::Z1"] * 20, object)])
    DSz = np.concatenate([DS, np.array(["dsZ"] * 20, object)])
    TUMz = np.concatenate([TUM, np.ones(20, bool)])
    _, metaz = evaluate_all(Cz, go, Lz, Dz, DSz, TUMz)
    check("non-measuring dataset excluded (tumour count unchanged)", metaz["n_tumour_cells"] == meta["n_tumour_cells"])

    print(f"\n  selftest: {ok}/{len(checks)} passed")
    return 0 if ok == len(checks) else 1


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "selftest":
        sys.exit(selftest())
    elif cmd == "run":
        sys.exit(main_run())
    print(f"unknown: {cmd}"); sys.exit(64)
