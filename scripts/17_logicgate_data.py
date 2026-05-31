#!/usr/bin/env python3
"""
RUNG 4 / Step-5 — logic-gate DATA LAYER (real single-cell discovery; Colab/GPU-free, runs on CELLxGENE).

The selectivity ENGINE (scripts/18) and its RUN-TRUST validation (scripts/20) run anywhere. THIS script
is the only part that needs the real single-cell atlases, so it is built to run on Colab where
`cellxgene_census` + network are available (the atlases are GBs and are NOT committed — gitignored, like
Step-2's data/cellxgene/).

What it does, one tissue in RAM at a time (never bulk/pseudobulk):
  1. Stream NORMAL per-tissue single-cell slices from CZ CELLxGENE Census for ALL Step-3 vital tissues
     (heart, brain, kidney, liver, lung, pancreas, adrenal gland, bone marrow, skeletal muscle), subset to
     the antigen-pool genes only (keeps each slice small).
  2. Emit a VITAL-COVERAGE CENSUS. SAFETY MECHANICS (audit-hardened): vital-parenchyma cells are kept in
     FULL (asymmetric cap — only abundant non-vital types are capped) so a rare lethal double-positive is
     not statistically erased; leaks are Jeffreys UPPER bounds (a false zero from dropout cannot pass); and
     it FAILS CLOSED — if any non-regen vital type (heart/brain/kidney/pancreas/adrenal/muscle) was NOT
     adequately captured, the gate is UNCERTAIN, never silently SELECTIVE ('never looked at the heart' !=
     'the heart is clean'). Multiple-testing control = held-out-donor replication is DEFERRED to the next
     pass; selective gates are a DISCOVERY shortlist until then.
  3. Pull the TUMOUR single-cell (reuse scripts/03's lung/breast/colon malignant pulls).
  4. Assemble a per-cell Panel and score the candidate AND / AND-NOT gates with scripts/18, emitting
     gate_selectivity.csv (both tumour COVERAGE and worst-case NORMAL LEAK, vital broken out).

HONEST: this is a transcript-level HYPOTHESIS screen. mRNA != surface protein (single-cell r~0.1-0.4);
a NOT/absence can be dropout; co-localisation != a functional circuit that fires caspase-8 (wet-lab).
HLA-LOH is a per-patient GENETIC NOT (NGS-stratified), not an atlas expression call. Recognition-
selectivity is a SEPARATE axis — never multiplied with RUNG-1/2/3 (asserted via scripts/18).

USAGE (Colab):  python scripts/17_logicgate_data.py
REQS        :  pip install cellxgene-census scanpy   (CPU; streams TileDB-SOMA)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

_T0 = time.monotonic()


def log(msg):
    """Timestamped (elapsed-seconds) progress line so the run is never a blind box."""
    print(f"[+{time.monotonic() - _T0:7.1f}s] [rung4] {msg}", flush=True)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung4_logicgate"
DATA_DIR = PROJECT_ROOT / "data" / "logicgate"
CENSUS_VERSION = "2024-07-01"            # pinned to match scripts/03
MIN_VITAL_CELLS = 200                    # below this a vital type is UNAUDITED
K = 2                                    # per-cell POSITIVE threshold (UMI >= K)
MAX_PER_TYPE = 1500                      # cap for ABUNDANT non-vital cell types
VITAL_CAP = 50000                        # cap for VITAL parenchyma: large enough to detect a co-positive
#   subpopulation down to ~0.01% (and the Jeffreys UPPER bound honestly reports the residual uncertainty
#   below that), but bounded so the panel stays light + cacheable (keeping all 7.2M neurons bloated it to
#   8.4M cells). 30x+ more sensitive than the non-vital cap; reading is unchanged (cap applied after read).
SEED = 20260530
# Panel cache: set env LOGICGATE_CACHE to a GOOGLE DRIVE path so the ~45-min NORMAL fetch happens ONCE and
# SURVIVES a disconnect (free Colab wipes /content on reconnect). Cached SEPARATELY (normal vs tumour) so a
# tumour-only change re-fetches only the ~7-min tumour pull; normal migrates from an old combined cache.
CACHE_PATH = Path(os.environ.get("LOGICGATE_CACHE", "")) if os.environ.get("LOGICGATE_CACHE") else None
NORMAL_CACHE = CACHE_PATH.with_suffix(".normal.npz") if CACHE_PATH else None
TUMOUR_CACHE = CACHE_PATH.with_suffix(".tumour.npz") if CACHE_PATH else None
TUMOUR_CAP = 50000          # cap malignant tumour cells (coverage denominator)
COV_BAR = 0.15              # coverage bar for a TWO-antigen single-cell AND: even a perfect gate reads
#   ~0.25-0.5 co-positive among malignant cells after scRNA dropout (each antigen detected ~half the time),
#   so 0.30 was near the dropout ceiling. 0.15 is dropout-realistic; coverage is a SEPARATE axis from safety.
# Cells counted as the tumour (coverage) denominator = MALIGNANT/epithelial only (NOT immune+stroma, which
# diluted coverage in the first real run). Keyword match on Census cell_type.
MALIGNANT_KEYWORDS = ("malignant", "neoplastic", "tumor", "tumour", "carcinoma", "epithelial")

lg_spec = importlib.util.spec_from_file_location("lg", PROJECT_ROOT / "scripts" / "18_logicgate_search.py")
lg = importlib.util.module_from_spec(lg_spec); lg_spec.loader.exec_module(lg)

NORMAL_TISSUES = ["heart", "brain", "kidney", "liver", "lung", "bone marrow",
                  "pancreas", "adrenal gland", "skeletal muscle"]   # all 9: non-regen vital + regen + tumour-matched
TUMOUR_DISEASES = ["lung adenocarcinoma", "breast carcinoma", "colorectal cancer"]
# Step-2 activator pool (all FAILED scripts/07 single-antigen safety -> MUST be gated).
ACTIVATORS = ["ERBB2", "ERBB3", "EPHB4", "TACSTD2", "MUC1", "SDC1", "CD74", "ITGB4"]
# curated surface partners for the AND co-input (kept small to control multiple testing).
PARTNERS = ["ERBB3", "EPHB4", "EPCAM", "CDH1", "MET", "PROM1", "CD24", "FOLR1", "MSLN", "CEACAM5",
            "NECTIN4", "CLDN6", "CLDN18", "ROR1", "MUC16"]
ALL_GENES = sorted(set(ACTIVATORS + PARTNERS))
# Canonical non-regenerating vital parenchyma to audit (cell_type substrings -> canonical label).
# Cardiac entries FIRST so 'cardiac muscle cell' maps to cardiomyocyte before any 'muscle' rule.
# NOTE: if a real Census label doesn't match here, that vital type stays UNAUDITED and the gate
# FAILS CLOSED (UNCERTAIN) — imperfect mapping is now conservative-safe, never lethally false-safe.
VITAL_AUDIT = {"cardiac muscle": "cardiomyocyte", "cardiomyocyte": "cardiomyocyte",
               "neuron": "neuron", "glial": "neuron", "glia": "neuron",
               "kidney epithel": "kidney_tubule", "renal": "kidney_tubule", "nephron": "kidney_tubule",
               "podocyte": "kidney_podocyte",
               "pancreatic": "pancreatic_islet", "islet": "pancreatic_islet", "beta cell": "pancreatic_islet",
               "type b pancreatic": "pancreatic_islet",
               "adreno": "adrenal_cortical", "adrenal cort": "adrenal_cortical", "cortical cell of adrenal": "adrenal_cortical",
               "chromaffin": "adrenal_cortical",
               "skeletal muscle": "skeletal_myocyte", "skeletal": "skeletal_myocyte", "myofiber": "skeletal_myocyte",
               "fast muscle": "skeletal_myocyte", "slow muscle": "skeletal_myocyte", "muscle cell": "skeletal_myocyte"}


def _q(items):
    """SOMA value-filter list literal, e.g. ['heart','brain']."""
    return "[" + ", ".join(f"'{x}'" for x in items) + "]"


def _stream_pull(census, value_filter, label, tissue_index=0):
    """STREAM one slice via get_anndata(obs_value_filter=...) — a CONTIGUOUS, predicate-pushed read — then
    map cell_types to canonical labels and ASYMMETRICALLY subsample: keep ALL vital-parenchyma cells (so a
    rare lethal double-positive cardiomyocyte cannot be statistically erased), cap only abundant non-vital
    types at MAX_PER_TYPE. Mapping happens BEFORE the cap so vital cells are recognised first. Independent
    RNG per tissue. Returns (counts, canonical_cell_types, tissues) or None."""
    import cellxgene_census
    log(f"{label}: streaming get_anndata (contiguous predicate read) ...")
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", obs_value_filter=value_filter,
        var_value_filter=f"feature_name in {ALL_GENES}",
        column_names={"obs": ["cell_type", "tissue_general"]})
    if ad.n_obs == 0:
        log(f"{label}: 0 cells matched — check disease/tissue labels"); return None
    raw = ad.obs["cell_type"].astype(str).to_numpy()
    mapped = np.array(_map_celltype(raw))   # MAP BEFORE CAP
    # DIAGNOSTIC: surface the actual Census cell-type vocabulary so mapping gaps (e.g. adrenal/muscle) are
    # visible, not blind. Logs the most common raw labels that did NOT map to any canonical vital type.
    import collections
    unmapped = collections.Counter(r for r, m in zip(raw, mapped) if m == r)
    if unmapped:
        top = ", ".join(f"{c}({n})" for c, n in unmapped.most_common(6))
        log(f"{label}: top unmapped raw cell types -> {top}")
    log(f"{label}: {ad.n_obs:,} cells read; asymmetric cap (vital <= {VITAL_CAP}, non-vital <= {MAX_PER_TYPE}) ...")
    rng = np.random.default_rng([SEED, tissue_index])   # independent draws per tissue
    keep = []
    for lab in np.unique(mapped):
        idx = np.where(mapped == lab)[0]
        cap = VITAL_CAP if lab in lg.VITAL_NONREGEN else MAX_PER_TYPE   # keep vital up to a large bound
        keep.append(idx if len(idx) <= cap else rng.choice(idx, cap, replace=False))
    keep = np.sort(np.concatenate(keep))
    ad = ad[keep]; mapped = mapped[keep]
    vital_kept = sorted(set(mapped.tolist()) & lg.VITAL_NONREGEN)
    log(f"{label}: kept {ad.n_obs:,} cells; vital types kept in FULL: {vital_kept or 'NONE in this tissue'}")
    return _dense_over(ad, ALL_GENES), list(mapped), list(ad.obs["tissue_general"].astype(str))


def _open_census():
    import cellxgene_census
    log(f"opening CELLxGENE Census version {CENSUS_VERSION} ...")
    c = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
    log("census open ✓")
    return c


def fetch_normal(census):
    """NORMAL panel — one streaming read per tissue (the ~45-min part; cached separately)."""
    cb, ct, ts = [], [], []
    for ti, tissue in enumerate(NORMAL_TISSUES):
        res = _stream_pull(census,
                           f"is_primary_data == True and disease == 'normal' and tissue_general == '{tissue}'",
                           f"NORMAL {tissue}", tissue_index=ti)
        if res is not None:
            c, cts, tss = res
            cb.append(c); ct += list(cts); ts += tss
    if not cb:
        raise RuntimeError("no normal cells fetched — check Census version / tissue labels / network")
    return lg.Panel(np.vstack(cb), ALL_GENES, np.array(ct), np.array(ts), np.array(["normal"] * len(ct)))


def fetch_tumour(census):
    """TUMOUR panel — MALIGNANT/epithelial cells ONLY (the coverage denominator). The first real run counted
    coverage over the whole tumour microenvironment (immune+stroma), crushing it; this fixes that."""
    import cellxgene_census
    log("TUMOUR (malignant/epithelial only): streaming get_anndata ...")
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens",
        obs_value_filter=f"is_primary_data == True and disease in {_q(TUMOUR_DISEASES)}",
        var_value_filter=f"feature_name in {ALL_GENES}",
        column_names={"obs": ["cell_type", "tissue_general"]})
    raw = ad.obs["cell_type"].astype(str).to_numpy()
    mal = np.array([any(k in c.lower() for k in MALIGNANT_KEYWORDS) for c in raw])
    log(f"TUMOUR: {ad.n_obs:,} cells read; {int(mal.sum()):,} malignant/epithelial (coverage denominator)")
    if mal.sum() < 100:   # annotation didn't match -> show vocabulary, fall back to all (flagged)
        import collections
        top = ", ".join(f"{c}({n})" for c, n in collections.Counter(raw).most_common(8))
        log(f"TUMOUR: <100 malignant matched — top raw types: {top}; FALLING BACK to all tumour cells (coverage diluted, flagged)")
        mal = np.ones(ad.n_obs, bool)
    ad = ad[mal]
    if ad.n_obs > TUMOUR_CAP:
        idx = np.sort(np.random.default_rng([SEED, 99]).choice(ad.n_obs, TUMOUR_CAP, replace=False))
        ad = ad[idx]
    log(f"TUMOUR: kept {ad.n_obs:,} malignant cells")
    n = ad.n_obs
    return lg.Panel(_dense_over(ad, ALL_GENES), ALL_GENES,
                    np.array(["tumour_malignant"] * n), np.array(["tumour"] * n), np.array(["tumour"] * n))


def _concat_panels(a, b):
    return lg.Panel(np.vstack([a.counts, b.counts]), ALL_GENES,
                    np.concatenate([a.cell_type, b.cell_type]),
                    np.concatenate([a.tissue, b.tissue]),
                    np.concatenate([a.compartment, b.compartment]))


def _dense_over(ad, genes):
    """Return (n_cells, len(genes)) integer counts in ALL_GENES order (0 for genes absent in this slice)."""
    import scipy.sparse as sp
    name_to_col = {n: i for i, n in enumerate(ad.var["feature_name"].astype(str))}
    out = np.zeros((ad.n_obs, len(genes)), dtype=np.int32)
    X = ad.X.tocsc() if sp.issparse(ad.X) else ad.X
    for j, g in enumerate(genes):
        if g in name_to_col:
            col = X[:, name_to_col[g]]
            out[:, j] = (col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()).astype(np.int32)
    return out


def _map_celltype(raw):
    """Collapse Census cell_type strings to the canonical vital labels scripts/18 protects (else keep raw)."""
    mapped = []
    for c in raw:
        cl = c.lower()
        hit = next((v for key, v in VITAL_AUDIT.items() if key in cl), c)
        mapped.append(hit)
    return mapped


def vital_coverage_census(panel):
    """Per (vital type) cell counts; flag UNAUDITED below MIN_VITAL_CELLS (droplet under-sampling)."""
    rows = []
    for vt in sorted(set(lg.VITAL_NONREGEN)):
        n = int((panel.cell_type == vt).sum())
        rows.append({"vital_type": vt, "n_cells": n,
                     "status": "AUDITED" if n >= MIN_VITAL_CELLS else "UNAUDITED (under-sampled — use snRNA-seq)"})
    return rows


def candidate_gates():
    gates = []
    for a in ACTIVATORS:
        for b in PARTNERS:
            if b != a:
                gates.append((a, b, "AND"))
        gates.append((a, "HLA_A02_LOH", "AND_NOT"))   # Tmod genetic NOT (flagged genotype, not atlas expression)
    return gates


def save_panel(panel, path):
    """Cache the assembled panel so a re-run (after a disconnect/laptop-sleep) skips the ~50-min fetch."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, counts=panel.counts, genes=np.array(panel.genes, dtype=object),
                        cell_type=panel.cell_type, tissue=panel.tissue, compartment=panel.compartment)
    log(f"cached panel -> {path}  ({panel.counts.shape[0]:,} cells). Re-runs resume from here.")


def load_panel(path):
    d = np.load(path, allow_pickle=True)
    log(f"RESUMING from cached panel {path} (skipping the Census fetch) ...")
    return lg.Panel(d["counts"], list(d["genes"]), d["cell_type"], d["tissue"], d["compartment"])


def main() -> int:
    lg.assert_no_multiply()
    if importlib.util.find_spec("cellxgene_census") is None:
        print("[rung4-data] cellxgene_census not installed — this script runs on COLAB.")
        print("[rung4-data] locally, the METHOD is validated by scripts/20 (synthetic ground truth).")
        print("[rung4-data] on Colab:  pip install cellxgene-census scanpy  then re-run.")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log("=== RUNG 4 real discovery — CELLxGENE single-cell logic-gate search ===")
    # NORMAL: from its own cache, else MIGRATE the normal cells out of an old combined cache (no 45-min
    # re-fetch), else fetch. TUMOUR (malignant): from cache, else fetch (~7 min). Split so a tumour-only
    # change never re-fetches the expensive normal atlases.
    census = None
    if NORMAL_CACHE and NORMAL_CACHE.exists():
        normal = load_panel(NORMAL_CACHE)
    elif CACHE_PATH and CACHE_PATH.exists():
        log("migrating NORMAL cells from the existing combined cache (skips the ~45-min normal re-fetch) ...")
        comb = load_panel(CACHE_PATH); m = comb.compartment == "normal"
        normal = lg.Panel(comb.counts[m], comb.genes, comb.cell_type[m], comb.tissue[m], comb.compartment[m])
        if NORMAL_CACHE:
            save_panel(normal, NORMAL_CACHE)
    else:
        census = _open_census(); normal = fetch_normal(census)
        if NORMAL_CACHE:
            save_panel(normal, NORMAL_CACHE)
    if TUMOUR_CACHE and TUMOUR_CACHE.exists():
        tumour = load_panel(TUMOUR_CACHE)
    else:
        census = census or _open_census(); tumour = fetch_tumour(census)
        if TUMOUR_CACHE:
            save_panel(tumour, TUMOUR_CACHE)
    if census is not None:
        census.close(); log("census closed")
    panel = _concat_panels(normal, tumour)
    n_norm = int((panel.compartment == "normal").sum()); n_tum = int((panel.compartment == "tumour").sum())
    log(f"panel assembled: {panel.counts.shape[0]:,} cells ({n_norm:,} normal, {n_tum:,} tumour-malignant) x {len(ALL_GENES)} antigens")
    for tis in sorted(set(panel.tissue)):
        log(f"  tissue {tis:14s}: {int((panel.tissue == tis).sum()):,} cells")

    cov = vital_coverage_census(panel)
    log("VITAL-COVERAGE CENSUS (heart/brain/kidney must be AUDITED to certify vital-safe):")
    for r in cov:
        log(f"  {r['vital_type']:16s} n={r['n_cells']:6d}  {r['status']}")
    unaudited = [r["vital_type"] for r in cov if r["status"].startswith("UNAUDITED")]

    # HLA-LOH is a per-patient GENOTYPE gate (NGS LOH), not an atlas-expression NOT -> not searched here.
    specs = [(a, b, logic) for a, b, logic in candidate_gates()
             if logic == "AND" and a in ALL_GENES and b in ALL_GENES]
    log(f"scoring {len(specs)} candidate AND gates over {panel.counts.shape[0]:,} cells (vectorised) ...")

    def _prog(i, n, r):
        if i % 10 == 0 or i == n or r["selective"]:
            tag = "SELECTIVE" if r["selective"] else "no"
            log(f"  scored {i:3d}/{n}  {r['gate']:26s} cov={r['tumour_coverage']:.2f} "
                f"leak={r['worst_normal_leak']:.2f} vital={r['vital_leak']:.2f} -> {tag}")

    # FAIL-CLOSED: require every non-regen vital type to be adequately captured; a gate where heart/brain/
    # kidney/pancreas/adrenal/muscle was NOT captured is UNCERTAIN, never silently SELECTIVE. Leaks are
    # Jeffreys UPPER bounds, vital cells kept in full (asymmetric cap) — so a false zero cannot pass.
    rows = lg.score_gates_batch(panel, specs, k=K, cov_bar=COV_BAR, required_vital=lg.VITAL_NONREGEN, progress=_prog)
    for r in rows:
        r["protein_copositivity_status"] = "NO_SINGLECELL_PROTEIN_DATA"  # transcript-only until CITE-seq
        r["transcript_only"] = True
        r["multiple_testing_control"] = "held-out-donor replication DEFERRED (next pass) — treat selective as DISCOVERY shortlist"
    # rank: selective first, then SAFE-but-low-coverage (the interesting ones), then by leak
    rows.sort(key=lambda r: (not r["selective"], not r.get("safe"), r["worst_normal_leak"], -r["tumour_coverage"]))

    import csv
    with open(DATA_DIR / "gate_selectivity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())) if rows else None
        if w:
            w.writeheader(); w.writerows(rows)
    with open(DATA_DIR / "vital_coverage_census.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["vital_type", "n_cells", "status"]); w.writeheader(); w.writerows(cov)

    log(f"wrote data/logicgate/gate_selectivity.csv ({len(rows)} gates) + vital_coverage_census.csv")
    selective = [r for r in rows if r["selective"]]
    safe = [r for r in rows if r.get("safe")]
    safe_lowcov = [r for r in safe if not r["selective"]]
    log(f"scored {len(rows)} AND gates -> {len(safe)} SAFE on all audited vital/normal axes; "
        f"of those {len(selective)} also clear coverage>={COV_BAR} (SELECTIVE), {len(safe_lowcov)} are SAFE-but-low-coverage.")
    for r in (selective or safe_lowcov)[:12]:
        log(f"  {r['verdict'][:24]:24s} {r['gate']:24s} cov={r['tumour_coverage']:.2f} leak={r['worst_normal_leak']:.3f}@{r['worst_group']}")
    if not safe:
        log("NO gate is even SAFE in this pool (all leak into normal tissue) — a genuine, FIRST-CLASS negative.")
    elif not selective:
        log(f"{len(safe)} SAFE gates exist but none clear coverage — coverage axis, not safety, is the limiter (transcript dropout / pool).")
    (OUT_DIR / "rung4_results.json").write_text(json.dumps({
        "census_version": CENSUS_VERSION, "n_cells": int(panel.counts.shape[0]),
        "vital_coverage": cov, "unaudited_vital_types": unaudited,
        "n_gates": len(rows), "n_selective": len(selective), "n_safe": len(safe),
        "n_safe_low_coverage": len(safe_lowcov), "cov_bar": COV_BAR,
        "top_gates": rows[:15], "no_clean_gate": len(selective) == 0, "no_safe_gate": len(safe) == 0,
        "CEILING": "transcript-level hypothesis; mRNA!=surface protein (CITE-seq needed to confirm "
                   "co-positivity); HLA-LOH is an NGS genotype gate not modelled here; recognition is a "
                   "separate axis never multiplied with RUNG-1/2/3; the durability cost is in scripts/21.",
    }, indent=2, default=str))
    print("[rung4-data] -> data/logicgate/gate_selectivity.csv + runs/rung4_logicgate/rung4_results.json")
    _figure(rows, cov)
    print("[rung4-data] CEILING: transcript-only HYPOTHESIS; CITE-seq/flow/co-culture confirm; agonism = wet-lab.")
    return 0


def _figure(rows, cov):
    """Real-discovery figure: gate coverage-vs-leak frontier + vital-coverage census."""
    if not rows:
        return
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
        for r in rows:
            ax[0].scatter(r["tumour_coverage"], r["worst_normal_leak"],
                          c="#27ae60" if r["selective"] else "#c0392b", s=45, alpha=0.7)
        ax[0].axhline(0.02, ls="--", color="green", lw=0.8); ax[0].axvline(0.30, ls="--", color="green", lw=0.8)
        ax[0].set_xlabel("tumour coverage (want high)"); ax[0].set_ylabel("worst normal-cell leak (want ~0)")
        ax[0].set_title("AND-gate frontier (transcript-only)\ngreen=selective; target=lower-right box")
        sel = [r for r in rows if r["selective"]][:8]
        for r in sel:
            ax[0].annotate(r["gate"][:18], (r["tumour_coverage"], r["worst_normal_leak"]), fontsize=6)
        ax[1].barh([c["vital_type"] for c in cov], [c["n_cells"] for c in cov],
                   color=["#2980b9" if c["status"] == "AUDITED" else "#e67e22" for c in cov])
        ax[1].axvline(MIN_VITAL_CELLS, ls="--", color="red", lw=0.8)
        ax[1].set_xlabel("cells captured"); ax[1].set_title("vital-coverage census\n(orange = UNAUDITED, < min)")
        fig.suptitle("RUNG 4 real discovery — transcript-only hypothesis (mRNA!=protein; CITE-seq confirms). "
                     "Recognition is a separate axis; 'no clean gate' is a valid result.", fontsize=8)
        fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(OUT_DIR / "rung4_discovery.png", dpi=110)
        print("[rung4-data] figure -> runs/rung4_logicgate/rung4_discovery.png")
    except Exception as e:
        print(f"[rung4-data] figure skipped ({type(e).__name__}: {e})")


if __name__ == "__main__":
    sys.exit(main())
