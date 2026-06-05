#!/usr/bin/env python3
"""
RUNG 8 — normal-tissue HLA-I expression heterogeneity (the measurement RUNG-7 said we were missing).

WHY
---
RUNG-7's entire safety result rides on ONE unsourced parameter: the fraction of NORMAL cells that are
"HLA-low" (so the Tmod LIR-1 blocker fails to fire and the broad activator kills them). This run MEASURES that
fraction, per vital cell type, from the CELLxGENE atlas, and feeds it straight back into RUNG-7's model to get
a DATA-GROUNDED per-organ off-tumour-toxicity floor for the HLA-A*02 gate.

WHAT IT COMPUTES
----------------
For each of the 9 normal tissues, for the non-regenerating VITAL cell types (scripts/18 VITAL_NONREGEN,
mapped via scripts/17 VITAL_AUDIT), pull HLA-A / HLA-B / HLA-C per cell (donor-resolved), then per
(vital cell type, donor): the HLA-low fraction (UMI below a threshold) + detection fraction. Headline =
the WORST-DONOR HLA-low fraction per vital type (consistent with RUNG-5/6 worst-donor safety). The blocker
in the deployed gate senses HLA-A, so HLA-A is primary; B/C are reported for context.

ENGINEERING (the three hard requirements)
-----------------------------------------
* RESUMABLE: one Drive tile per tissue (RUNG8_CACHE). A 4-hour-cap disconnect loses nothing — re-run and it
  skips completed tissues and continues. (Only 3 genes are pulled, so this is far lighter than RUNG-5.)
* FOREGROUND-VISIBLE LOGGING: a background HEARTBEAT thread prints "[heartbeat] <current step> | elapsed | RAM"
  every ~20s, so even during a long opaque Census query you SEE it is alive (not stuck). Every step also logs
  a flushed line. Run via runlog (python -u) so the subprocess stream reaches the Colab cell live.
* GPU: NOT used and NOT needed. We pull only 3 genes; the bottleneck is the Census fetch (network/disk), and
  the per-donor aggregation is a trivial numpy groupby. Stated honestly rather than bolting on idle GPU code.

HONEST CEILING
--------------
mRNA HLA != surface MHC-I protein (the blocker senses protein); scRNA DROPOUT inflates the apparent HLA-low
fraction, so the headline is an UPPER BOUND on true HLA-low (conservative -> over-estimates toxicity, the safe
direction). HLA-I is IFN-gamma INDUCIBLE, so resting atlas expression can understate induced levels. scRNA
resolves the HLA-A GENE, not the A*02 ALLELE. These are stated, not hidden.

USAGE
  python scripts/29_hla_heterogeneity.py selftest         # LOCAL: synthetic, no Census, validates the math
  RUNG8_CACHE=/content/drive/MyDrive/cancer-recon python scripts/29_hla_heterogeneity.py run   # Colab real run
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import threading
import time
from pathlib import Path

import numpy as np

_T0 = time.monotonic()
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung8_hla"
RESULT_JSON = OUT_DIR / "rung8_hla_heterogeneity.json"
FIGURE_PNG = OUT_DIR / "rung8_hla.png"

HLA_GENES = ["HLA-A", "HLA-B", "HLA-C"]
LOW_THRESHOLDS = (1, 2)              # a cell is 'HLA-low' at k if UMI < k  (k=1 => undetected)
PER_DONOR_CAP = 600                  # cap vital cells per (type, donor): bounds RAM, keeps donors powered
MIN_CELLS_PER_DONOR = 30             # a donor must have >= this many cells to enter the worst-donor max
CACHE = Path(os.environ["RUNG8_CACHE"]) if os.environ.get("RUNG8_CACHE") else None


def log(msg):
    print(f"[+{time.monotonic() - _T0:7.1f}s] [rung8] {msg}", flush=True)


def _ram_gb():
    import resource
    m = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return m / 1e9 if m > 1e7 else m / 1e6      # bytes on macOS, KB on Linux/Colab


# ---------------------------------------------------------------------------
#  HEARTBEAT — foreground visibility so a long/opaque step never looks "stuck"
# ---------------------------------------------------------------------------
class Heartbeat:
    def __init__(self, interval=20):
        self.interval = interval
        self.label = "starting"
        self._stop = False

    def set(self, label):
        self.label = label
        log(label)

    def _run(self):
        while not self._stop:
            for _ in range(self.interval * 2):       # responsive stop, ~0.5s granularity
                if self._stop:
                    return
                time.sleep(0.5)
            if not self._stop:
                print(f"[+{time.monotonic() - _T0:7.1f}s] [heartbeat] {self.label} | RAM {_ram_gb():.1f}GB",
                      flush=True)

    def start(self):
        threading.Thread(target=self._run, daemon=True).start()
        return self

    def stop(self):
        self._stop = True


HB = Heartbeat()


# ---------------------------------------------------------------------------
#  reuse the tested data conventions (scripts/17 = tissues / vital map / census)
# ---------------------------------------------------------------------------
def _load(name, mod):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / mod)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _vital_label(cell_type_str: str, vital_audit: dict):
    """Map a Census cell_type string to a canonical VITAL label, or None if not a protected vital type."""
    cl = cell_type_str.lower()
    return next((v for key, v in vital_audit.items() if key in cl), None)


def _codes(col):
    """Arrow column -> (int codes, small unique list) WITHOUT one Python str per cell (the brain-OOM fix)."""
    import pyarrow as pa
    col = col.combine_chunks()
    if not pa.types.is_dictionary(col.type):
        col = col.dictionary_encode()
    return col.indices.to_numpy(zero_copy_only=False), [str(x) for x in col.dictionary.to_pylist()]


# ---------------------------------------------------------------------------
#  per-tissue fetch (vital cells only, donor-capped, memory-safe via Arrow codes)
# ---------------------------------------------------------------------------
def pull_vital_hla(census, d4, tissue, ti, n_tissue):
    import cellxgene_census
    exp = census["census_data"]["homo_sapiens"]
    vf = (f"is_primary_data == True and disease == 'normal' and tissue_general == '{tissue}'")
    HB.set(f"[{ti + 1}/{n_tissue}] {tissue}: reading obs metadata (cell_type, donor) ...")
    tbl = exp.obs.read(value_filter=vf,
                       column_names=["soma_joinid", "cell_type", "donor_id", "dataset_id"]).concat()
    total = tbl.num_rows
    if total == 0:
        log(f"[{ti + 1}/{n_tissue}] {tissue}: 0 normal cells matched"); return None
    jid = tbl.column("soma_joinid").to_numpy()
    ct_codes, ct_vals = _codes(tbl.column("cell_type"))
    # map the SMALL unique cell-type set -> vital label (or None); per-cell label via index (cheap)
    mapped_unique = np.array([_vital_label(c, d4.VITAL_AUDIT) for c in ct_vals], dtype=object)
    is_vital = np.array([m is not None for m in mapped_unique])[ct_codes]
    n_vital = int(is_vital.sum())
    HB.set(f"[{ti + 1}/{n_tissue}] {tissue}: {total:,} normal cells, {n_vital:,} vital -> capping per donor")
    if n_vital == 0:
        return {"tissue": tissue, "counts": np.zeros((0, 3), np.int32),
                "label": np.array([], object), "donor": np.array([], object)}
    labels = mapped_unique[ct_codes]                                   # per-cell vital label or None
    ds_codes, ds_vals = _codes(tbl.column("dataset_id"))
    dn_codes, dn_vals = _codes(tbl.column("donor_id"))
    donor_gid = ds_codes.astype(np.int64) * (int(dn_codes.max()) + 1) + dn_codes

    # per (vital label, donor) keep first PER_DONOR_CAP -> donor-aware, bounds RAM
    keep = []
    vidx = np.where(is_vital)[0]
    lab_v = labels[vidx]
    for lab in np.unique(lab_v):
        lab_mask = vidx[lab_v == lab]
        for g in np.unique(donor_gid[lab_mask]):
            keep.extend(lab_mask[donor_gid[lab_mask] == g][:PER_DONOR_CAP].tolist())
    keep = np.sort(np.array(keep, dtype=np.int64))
    sel_jid = jid[keep]
    HB.set(f"[{ti + 1}/{n_tissue}] {tissue}: materialising HLA-A/B/C for {len(sel_jid):,} vital cells ...")
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", obs_coords=sel_jid.tolist(),
        var_value_filter=f"feature_name in {HLA_GENES}",
        column_names={"obs": ["cell_type", "donor_id", "dataset_id"]})
    # order columns as HLA_GENES
    vnames = list(ad.var["feature_name"]) if "feature_name" in ad.var else list(ad.var_names)
    order = [vnames.index(g) for g in HLA_GENES if g in vnames]
    X = ad.X[:, order]
    counts = np.asarray(X.todense() if hasattr(X, "todense") else X).astype(np.int32)
    if counts.shape[1] != 3:                                            # a gene missing in this slice -> pad
        full = np.zeros((counts.shape[0], 3), np.int32)
        present = [g for g in HLA_GENES if g in vnames]
        for j, g in enumerate(present):
            full[:, HLA_GENES.index(g)] = counts[:, j]
        counts = full
    donor = np.array([f"{ds}::{dn}" for ds, dn in
                      zip(ad.obs["dataset_id"].astype(str), ad.obs["donor_id"].astype(str))], dtype=object)
    lab = np.array([_vital_label(c, d4.VITAL_AUDIT) for c in ad.obs["cell_type"].astype(str)], dtype=object)
    log(f"[{ti + 1}/{n_tissue}] {tissue}: kept {counts.shape[0]:,} vital cells "
        f"(types: {sorted(set(x for x in lab if x))}) RAM {_ram_gb():.1f}GB")
    return {"tissue": tissue, "counts": counts, "label": lab, "donor": donor}


# ---------------------------------------------------------------------------
#  aggregation: per vital type HLA-low among DATASETS THAT ACTUALLY MEASURED HLA
# ---------------------------------------------------------------------------
#  CRITICAL (v2 fix): CELLxGENE returns 0 for BOTH "measured & truly zero" AND "this dataset never measured
#  the gene". Counting the second as "HLA-low" inflated every number (adrenal showed 0% detection = its
#  datasets simply lack HLA-A; every type hit worst-donor 100% off a single artifact donor). We now treat a
#  dataset as HLA-MEASURING iff it has >=1 cell with HLA>0 ANYWHERE, drop cells from non-measuring datasets,
#  and report a robust per-donor DISTRIBUTION (median/p90/worst) instead of the brittle single worst donor.
#  Residual ceiling: a dataset that measured HLA but where a cell type is genuinely ~0 is correctly kept
#  (real biology, e.g. immune-privileged neurons); mRNA dropout still inflates low (upper bound). The
#  gold-standard fix (Census feature_dataset_presence_matrix) needs a refetch and is noted, not silently done.
def aggregate(counts, label, donor, gene_idx=0):
    """Returns (per_type, coverage_meta). gene_idx 0 = HLA-A (the sensed gene)."""
    a = counts[:, gene_idx].astype(np.int64)
    dataset = np.array([str(d).split("::")[0] for d in donor], dtype=object)
    all_ds = set(dataset.tolist())
    measuring = {d for d in all_ds if a[dataset == d].max(initial=0) > 0}    # dataset detects HLA somewhere
    is_meas = np.array([d in measuring for d in dataset], dtype=bool)
    meta = {"n_datasets_total": len(all_ds),
            "n_datasets_excluded_unmeasured": len(all_ds) - len(measuring),
            "n_cells_excluded_unmeasured": int((~is_meas).sum())}

    out = {}
    for t in sorted(set(x for x in label if x)):
        m = (label == t) & is_meas
        if m.sum() == 0:
            out[t] = {"n_cells_measured": 0, "note": "no cells from HLA-measuring datasets (excluded)"}
            continue
        c = a[m]; dn = donor[m]
        donor_low = []
        for d in set(dn.tolist()):
            cc = c[dn == d]
            if len(cc) >= MIN_CELLS_PER_DONOR:
                donor_low.append(float((cc < 1).mean()))
        dl = np.array(sorted(donor_low)) if donor_low else np.array([])
        out[t] = {
            "n_cells_measured": int(m.sum()),
            "n_datasets_measuring": int(len(set(dataset[m].tolist()))),
            "n_informative_donors": int(len(dl)),
            "pooled_low_k1": round(float((c < 1).mean()), 4),
            "pooled_detect": round(float((c >= 1).mean()), 4),
            "median_donor_low_k1": round(float(np.median(dl)), 4) if len(dl) else None,
            "p90_donor_low_k1": round(float(np.percentile(dl, 90)), 4) if len(dl) else None,
            "worst_donor_low_k1": round(float(dl.max()), 4) if len(dl) else None,
        }
    return out, meta


def couple_to_rung7(worst_hla_low_overall: float):
    """Plug the MEASURED worst-vital HLA-low fraction into RUNG-7's gate model -> data-grounded FPR."""
    try:
        r7 = _load("r7", "28_andnot_gate_discrimination.py")
        p = dict(r7.PARAMS, normal_hla_low_frac=float(worst_hla_low_overall))
        res = r7.evaluate(p, np.random.default_rng(r7.SEED))
        return {"normal_hla_low_frac_measured": round(float(worst_hla_low_overall), 4),
                "data_grounded_FPR": res["fpr"], "data_grounded_TPR": res["tpr"],
                "note": "RUNG-7 model re-run with the MEASURED worst-vital HLA-low fraction in place of the "
                        "assumed 5%. FPR = the gate's off-tumour toxicity at this measured blocker-failure rate."}
    except Exception as e:
        return {"error": f"could not couple to RUNG-7: {type(e).__name__}: {e}"}


# ---------------------------------------------------------------------------
def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    HB.start()
    d4 = _load("d4", "17_logicgate_data.py")
    census = None                                   # LAZY: opened only if a tile is missing (network-free re-aggregate)
    tissues = d4.NORMAL_TISSUES
    tile_dir = CACHE if CACHE else (OUT_DIR / "tiles")
    tile_dir.mkdir(parents=True, exist_ok=True)
    log(f"cache/tiles -> {tile_dir}  (resumable: completed tissues are skipped on re-run)")

    all_counts, all_label, all_donor = [], [], []
    for ti, tissue in enumerate(tissues):
        tile = tile_dir / f"rung8_{tissue.replace(' ', '_')}.npz"
        if tile.exists():
            d = np.load(tile, allow_pickle=True)
            log(f"[{ti + 1}/{len(tissues)}] {tissue}: RESUMED from tile ({d['counts'].shape[0]:,} cells)")
            all_counts.append(d["counts"]); all_label.append(d["label"]); all_donor.append(d["donor"])
            continue
        if census is None:                          # first missing tile -> open Census now (only if needed)
            import cellxgene_census
            HB.set(f"opening CELLxGENE Census {d4.CENSUS_VERSION} ...")
            census = cellxgene_census.open_soma(census_version=d4.CENSUS_VERSION)
        res = pull_vital_hla(census, d4, tissue, ti, len(tissues))
        if res is None:
            continue
        np.savez_compressed(tile, counts=res["counts"], label=res["label"], donor=res["donor"])
        log(f"[{ti + 1}/{len(tissues)}] {tissue}: tile checkpointed -> {tile.name} (safe to disconnect now)")
        all_counts.append(res["counts"]); all_label.append(res["label"]); all_donor.append(res["donor"])

    HB.set("all tissues done — aggregating per (vital type, donor) ...")
    counts = np.vstack(all_counts) if all_counts else np.zeros((0, 3), np.int32)
    label = np.concatenate(all_label) if all_label else np.array([], object)
    donor = np.concatenate(all_donor) if all_donor else np.array([], object)
    per_type, coverage = aggregate(counts, label, donor, gene_idx=0)  # HLA-A (the sensed gene)
    log(f"dataset coverage: {coverage['n_datasets_excluded_unmeasured']}/{coverage['n_datasets_total']} "
        f"datasets excluded as non-HLA-measuring ({coverage['n_cells_excluded_unmeasured']:,} cells dropped)")

    # headline = the vital type with the highest POOLED HLA-low fraction (among HLA-measuring datasets);
    # pooled (not single-worst-donor) is the robust, RUNG-7-relevant 'fraction of normal cells HLA-low'.
    valid = {t: v for t, v in per_type.items() if v.get("n_cells_measured", 0) > 0 and v.get("pooled_low_k1") is not None}
    worst_type = max(valid, key=lambda t: valid[t]["pooled_low_k1"]) if valid else None
    worst_pooled = valid[worst_type]["pooled_low_k1"] if worst_type else 0.0
    coupling = couple_to_rung7(worst_pooled)

    result = {
        "tag": "rung8_hla_heterogeneity",
        "question": "What fraction of NORMAL vital cells are HLA-low (so the Tmod blocker fails)? — grounds "
                    "RUNG-7's single load-bearing parameter with atlas data.",
        "census_version": d4.CENSUS_VERSION, "hla_genes": HLA_GENES, "sensed_gene": "HLA-A",
        "n_cells_total": int(counts.shape[0]),
        "dataset_coverage": coverage,
        "per_vital_type": per_type,
        "worst_vital_type": worst_type,
        "worst_vital_pooled_HLA_low_k1": round(float(worst_pooled), 4),
        "rung7_coupling": coupling,
        "CEILING": "v2 fix: only DATASETS that actually measured HLA are counted (Census returns 0 for "
                   "'not measured' too -> v1 conflated unmeasured with HLA-low; see dataset_coverage). Residual: "
                   "mRNA HLA != surface MHC-I protein; scRNA DROPOUT still inflates HLA-low (UPPER BOUND -> "
                   "conservatively over-estimates toxicity); HLA-I is IFN-gamma inducible (resting atlas may "
                   "understate induced); HLA-A GENE not the A*02 ALLELE. Headline = POOLED among measuring "
                   "datasets (robust); per-type worst/p90/median donor also reported. Gold-standard fix = "
                   "Census feature_dataset_presence_matrix (needs refetch).",
        "INTERPRETATION": "This replaces RUNG-7's assumed 5% normal HLA-low fraction with a measured, "
                          "per-vital-tissue number. The worst vital type sets the gate's off-tumour-toxicity "
                          "floor (data_grounded_FPR). If the measured HLA-low fraction is well below 5%, the "
                          "blocker is more reliable than assumed (good); if higher, RUNG-7's safety was "
                          "optimistic and the gate is more toxic than modelled.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    log(f"wrote {RESULT_JSON}")
    log(f"WORST vital type = {worst_type}  pooled HLA-low(k1) = {worst_pooled:.1%}  "
        f"-> data-grounded gate FPR = {coupling.get('data_grounded_FPR', 'n/a')}")
    log(f"{'vital type':22} {'n_meas':>8} {'donors':>6} {'pooled_low':>10} {'p90_donor':>9} {'detect':>7}")
    for t, v in sorted(valid.items(), key=lambda kv: -kv[1]["pooled_low_k1"]):
        log(f"  {t:22} {v['n_cells_measured']:>8,} {v['n_informative_donors']:>6} "
            f"{v['pooled_low_k1']:>10.1%} {(v['p90_donor_low_k1'] or 0):>9.1%} {v['pooled_detect']:>7.1%}")
    HB.stop()
    _make_figure(valid)
    return 0


def _make_figure(per_type):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        log(f"matplotlib unavailable ({e}); skipped figure"); return
    valid = {t: v for t, v in per_type.items() if v.get("pooled_low_k1") is not None}
    if not valid:
        log("no measured vital types -> no figure"); return
    items = sorted(valid.items(), key=lambda kv: kv[1]["pooled_low_k1"], reverse=True)
    names = [f"{t} (n={v['n_cells_measured']:,})" for t, v in items]
    pooled = [v["pooled_low_k1"] * 100 for _, v in items]
    p90 = [(v["p90_donor_low_k1"] or 0) * 100 for _, v in items]
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(10, max(3, 0.5 * len(names) + 1.5)))
    ax.barh(y, p90, color="#C1432B", label="p90-donor HLA-low (high-end)")
    ax.barh(y, pooled, height=0.5, color="#2B6CB0", label="pooled HLA-low (headline)")
    ax.axvline(5, ls="--", color="grey", label="RUNG-7 assumed 5%")
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("% of vital cells HLA-A-low (UMI=0) among HLA-MEASURING datasets — UPPER BOUND (dropout)")
    ax.set_title("RUNG-8 (v2): normal-tissue HLA-low fraction per vital cell type\n"
                 "(measuring datasets only; grounds RUNG-7's safety parameter)", fontsize=11)
    ax.legend(fontsize=8); ax.grid(axis="x", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIGURE_PNG, dpi=130)
    log(f"wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    """Synthetic, no Census — validates the aggregation math, worst-donor logic, and RUNG-7 coupling."""
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # two HLA-MEASURING donors per type (each has some HLA>0); cardiomyocyte D2 is all-zero (real silence)
    rng = np.random.default_rng(8)
    rows = []
    def block(label, donor, n, hla_low):
        a = np.zeros(n, int) if hla_low else rng.integers(3, 30, n)
        for i in range(n):
            rows.append((label, donor, a[i]))
    block("cardiomyocyte", "dsA::D1", 100, False)    # measuring dataset, HLA-high
    block("cardiomyocyte", "dsA::D2", 100, True)     # SAME dataset dsA (measuring) but this donor all-zero
    block("neuron", "dsB::D3", 100, False)
    block("neuron", "dsB::D4", 100, False)
    label = np.array([r[0] for r in rows], object)
    donor = np.array([r[1] for r in rows], object)
    counts = np.array([[r[2], r[2], r[2]] for r in rows], np.int32)

    agg, cov = aggregate(counts, label, donor, gene_idx=0)
    check("aggregate returns (per_type, coverage) tuple", isinstance(cov, dict) and "n_datasets_total" in cov)
    check("both vital types aggregated", set(agg) == {"cardiomyocyte", "neuron"})
    check("cardiomyocyte pooled HLA-low(k1) == 50% (D2 zero, in a MEASURING dataset -> kept)",
          abs(agg["cardiomyocyte"]["pooled_low_k1"] - 0.5) < 1e-9)
    check("cardiomyocyte worst-donor HLA-low(k1) == 100% (the all-zero donor D2)",
          abs(agg["cardiomyocyte"]["worst_donor_low_k1"] - 1.0) < 1e-9)
    check("neuron pooled HLA-low ~ 0 (all detected)", agg["neuron"]["pooled_low_k1"] < 0.05)
    check("worst >= pooled (per type)",
          all(agg[t]["worst_donor_low_k1"] >= agg[t]["pooled_low_k1"] - 1e-9 for t in agg))

    # THE v2 FIX: a dataset with NO HLA anywhere (unmeasured) is excluded, NOT counted as HLA-low.
    rows3 = ([("kidney_tubule", "dsMeas::D1", v) for v in rng.integers(3, 30, 80)]   # measuring
             + [("kidney_tubule", "dsUnmeas::D9", 0) for _ in range(80)])            # never measured -> exclude
    l3 = np.array([r[0] for r in rows3], object); d3 = np.array([r[1] for r in rows3], object)
    c3 = np.array([[r[2]] * 3 for r in rows3], np.int32)
    agg3, cov3 = aggregate(c3, l3, d3, gene_idx=0)
    check("unmeasured dataset excluded (1 of 2 datasets dropped)", cov3["n_datasets_excluded_unmeasured"] == 1)
    check("kidney_tubule pooled HLA-low ~ 0 (the unmeasured all-zero dataset was NOT counted as low)",
          agg3["kidney_tubule"]["pooled_low_k1"] < 0.05)

    # MIN_CELLS gate: a tiny donor must not enter the informative-donor set
    rows2 = [("neuron", "dsX::tiny", 7)] * 5 + [("neuron", "dsX::big", v) for v in rng.integers(5, 30, 80)]
    lab2 = np.array([r[0] for r in rows2], object); dn2 = np.array([r[1] for r in rows2], object)
    cnt2 = np.array([[r[2]] * 3 for r in rows2], np.int32)
    agg2, _ = aggregate(cnt2, lab2, dn2, gene_idx=0)
    check("tiny donor (<MIN_CELLS) excluded from informative donors", agg2["neuron"]["n_informative_donors"] == 1)

    # RUNG-7 coupling runs and FPR rises with the measured HLA-low fraction
    c_lo = couple_to_rung7(0.01).get("data_grounded_FPR", None)
    c_hi = couple_to_rung7(0.20).get("data_grounded_FPR", None)
    check("RUNG-7 coupling returns an FPR", isinstance(c_lo, float))
    check("data-grounded FPR rises with measured HLA-low frac", (c_hi or 0) > (c_lo or 0))

    # heartbeat starts/stops cleanly
    hb = Heartbeat(interval=1).start(); hb.set("selftest heartbeat"); time.sleep(0.1); hb.stop()
    check("heartbeat starts and stops without error", True)

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", nargs="?", default="run", choices=["run", "selftest"])
    args = ap.parse_args()
    sys.exit(selftest() if args.mode == "selftest" else main_run())
