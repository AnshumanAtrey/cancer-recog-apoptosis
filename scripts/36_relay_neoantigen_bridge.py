#!/usr/bin/env python3
"""
RUNG 12P / bridge — does the gated relay UNLOCK RUNG-11's "too risky" neoantigen handles? (laptop, instant)

THE BRIDGE
----------
RUNG-11 found 175 presented public-neoantigen handles, but 104 are `tcr_dependent`: the WILD-TYPE peptide is
also presented on healthy cells, so a per-cell TCR-T's safety rests entirely on the TCR discriminating two
near-identical pMHC surfaces (the MAGE-A3 failure class). RUNG-12P/B found a per-hop-gated death wave tolerates
a much higher per-step false-positive than a per-cell gate: safe ceiling q_n ~ 0.17 (3D, conservative) / 0.30
(2D) vs the per-cell worst-donor bar 0.02. This run JOINS them:

  For each RUNG-11 handle, estimate its per-cell false-positive (cross-reactivity), then ask:
    - PER-CELL usable?  cross-reactivity <= 0.02 (R5 bar)  -> deployable as a TCR-T today
    - RELAY usable?     cross-reactivity <= 0.17 (3D relay ceiling) -> deployable as a per-hop relay gate
  How much ADDRESSABILITY does switching to the relay architecture UNLOCK, per cancer?

THE PROXY (transparent; the one real unknown is SWEPT, not guessed)
-------------------------------------------------------------------
A handle's per-cell false-positive q_n = P(normal cell presents WT pMHC) x P(anti-mutant receptor binds WT).
  - P(WT presented) = presentation_factor(wt_rank): 1 at wt_rank<=0.5 (strong), linearly to 0 at wt_rank>=2
    (so CLEAN handles, WT not presented, have q_n=0 -> always safe).
  - P(receptor binds WT | presented) = beta x anchor_discount, where beta is the TCR's cross-bind onto a
    single-residue-different pMHC. beta is EXACTLY what RUNG-12 (AlphaFold-Multimer/ESM) would measure per
    handle; here it is UNKNOWN so we SWEEP it over [0,1]. ANCHOR mutations (mutation at an MHC anchor ->
    WT pMHC conformationally distinct) get a discount (better discrimination), per the structural literature.

So q_n(handle, beta) = presentation_factor(wt_rank) * beta * (ANCHOR_DISCOUNT if anchor else 1). The ROBUST,
beta-independent result is the RELAXATION FACTOR: the relay raises the tolerable beta by relay_ceiling/0.02
(~8.6x in 3D). The addressability numbers are reported AS A FUNCTION of beta (honest about the unknown).

HONEST CEILING
--------------
beta is swept globally, not measured per handle (RUNG-12's job). presentation_factor is a transparent
rank->presentation proxy. The relay ceiling inherits RUNG-12P/B's percolation-abstraction caveats. Population
frequencies inherit RUNG-11's (literature point estimates, joined datasets). This is an integration/what-if
bridge, not per-handle truth -- but the relaxation factor and the SHAPE (relay sustains coverage as beta rises)
are robust.

USAGE
  python scripts/36_relay_neoantigen_bridge.py            # -> JSON + figure
  python scripts/36_relay_neoantigen_bridge.py selftest   # proxy + classification + coverage logic checks
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUNG11_JSON = PROJECT_ROOT / "runs" / "rung11_neoantigen" / "rung11_neoantigen_addressability.json"
RUNG12B_JSON = PROJECT_ROOT / "runs" / "rung12pB_relay" / "rung12pB_relay.json"
OUT_DIR = PROJECT_ROOT / "runs" / "rung12p_bridge"
RESULT_JSON = OUT_DIR / "rung12p_bridge.json"
FIGURE_PNG = OUT_DIR / "rung12p_bridge.png"

BINDER_RANK = 2.0
STRONG_RANK = 0.5
PER_CELL_BAR = 0.02           # R5/R7 worst-donor per-cell false-positive bar
ANCHOR_DISCOUNT = 0.2         # anchor mutations -> WT pMHC conformationally distinct -> ~5x better discrimination


def _load(name, mod):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / mod)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


# ---------------------------------------------------------------------------
#  The proxy (pure -> selftest exercises it directly).
# ---------------------------------------------------------------------------
def presentation_factor(wt_rank: float) -> float:
    """How present is the WT pMHC on normal cells? 1 if strong binder (<=0.5), 0 if non-binder (>=2)."""
    return float(np.clip((BINDER_RANK - wt_rank) / (BINDER_RANK - STRONG_RANK), 0.0, 1.0))


def qn_of(wt_rank: float, tier: str, anchor: bool, beta: float) -> float:
    """Per-cell false-positive (cross-reactivity) of a handle at TCR cross-bind `beta`."""
    if tier == "clean":
        return 0.0
    disc = ANCHOR_DISCOUNT if (anchor or tier == "anchor") else 1.0
    return presentation_factor(wt_rank) * beta * disc


# ---------------------------------------------------------------------------
def load_handles(d33):
    """Reconstruct per-handle records from the RUNG-11 result JSON (presented handles only)."""
    d11 = json.loads(RUNG11_JSON.read_text())
    posmap = {(g, f"{wt}{pos}{mut}"): pos for (g, _acc, pos, wt, mut, _prev) in d33.DRIVERS}
    handles = []
    for _drv, info in d11["per_driver"].items():
        gene, mut_label, prev = info["gene"], info["mutation"], info["cancer_prev"]
        pos = posmap.get((gene, mut_label))
        for allele, b in info["by_allele"].items():
            handles.append({"gene": gene, "pos": pos, "mut_label": mut_label, "allele": allele,
                            "tier": b["tier"], "wt_rank": b["wt_rank"], "anchor": b["anchor"],
                            "cancer_prev": prev})
    return handles


def classify(handles, beta, relay_ceiling):
    """Tag each handle usable_percell / usable_relay at cross-bind beta."""
    for h in handles:
        q = qn_of(h["wt_rank"], h["tier"], h["anchor"], beta)
        h["qn"] = q
        h["usable_percell"] = q <= PER_CELL_BAR
        h["usable_relay"] = q <= relay_ceiling
    return handles


# ---------------------------------------------------------------------------
def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d33 = _load("d33", "33_neoantigen_addressability.py")
    if not RUNG11_JSON.exists() or not RUNG12B_JSON.exists():
        print(f"[bridge] need {RUNG11_JSON.name} and {RUNG12B_JSON.name}; run RUNG-11 and RUNG-12P/B first.")
        return 2
    d12 = json.loads(RUNG12B_JSON.read_text())
    relay_2d = d12["relay_safe_q_n_ceiling_at_1pct"]
    relay_3d = d12["sensitivity_3D"]["safe_q_n_ceiling_at_1pct"]
    relay_ceiling = relay_3d                                   # conservative primary
    relax_factor = round(relay_ceiling / PER_CELL_BAR, 1)

    handles = load_handles(d33)
    print(f"[bridge] {len(handles)} RUNG-11 presented handles; relay ceiling q_n={relay_ceiling} (3D), "
          f"per-cell bar {PER_CELL_BAR} -> relaxation {relax_factor}x")

    # sweep the unknown TCR cross-bind beta; record per-cancer coverage for per-cell vs relay
    betas = [round(b, 3) for b in np.linspace(0.0, 1.0, 21)]
    sweep = []
    for beta in betas:
        classify(handles, beta, relay_ceiling)
        cov_pc = d33.coverage(handles, d33.HLA_PANEL, lambda h: h["usable_percell"])
        cov_rl = d33.coverage(handles, d33.HLA_PANEL, lambda h: h["usable_relay"])
        sweep.append({"beta": beta,
                      "percell": {c: cov_pc[c]["central"] for c in d33.CANCERS},
                      "relay": {c: cov_rl[c]["central"] for c in d33.CANCERS}})

    # representative "imperfect TCR" operating point for the headline
    beta_star = 0.5
    classify(handles, beta_star, relay_ceiling)
    cov_pc = d33.coverage(handles, d33.HLA_PANEL, lambda h: h["usable_percell"])
    cov_rl = d33.coverage(handles, d33.HLA_PANEL, lambda h: h["usable_relay"])
    unlocked = sorted([h for h in handles if h["usable_relay"] and not h["usable_percell"]],
                      key=lambda h: -max(h["cancer_prev"].values()))
    unlocked_list = [{"handle": f"{h['gene']} {h['mut_label']} / {h['allele']}", "tier": h["tier"],
                      "wt_rank": h["wt_rank"], "qn": round(h["qn"], 3),
                      "cancers": h["cancer_prev"]} for h in unlocked[:20]]
    n_clean = sum(1 for h in handles if h["tier"] == "clean")

    result = {
        "tag": "rung12p_relay_neoantigen_bridge",
        "question": "How much of RUNG-11's neoantigen addressability does the RUNG-12P gated relay UNLOCK vs a "
                    "per-cell TCR-T, as a function of the (unknown, RUNG-12) TCR mut-vs-WT cross-bind beta?",
        "inputs": {"rung11_handles": len(handles), "n_clean_always_safe": n_clean,
                   "per_cell_bar": PER_CELL_BAR, "relay_ceiling_3d": relay_3d, "relay_ceiling_2d": relay_2d,
                   "anchor_discount": ANCHOR_DISCOUNT},
        "RELAXATION_FACTOR": f"{relax_factor}x (relay tolerates beta up to relay_ceiling/0.02 higher than per-cell)",
        "proxy": "q_n = presentation_factor(wt_rank) * beta * (anchor_discount if anchor). beta swept (RUNG-12 "
                 "would measure it per handle). CLEAN handles (WT not presented) -> q_n=0 -> always safe.",
        "beta_sweep_coverage": sweep,
        "headline_at_beta": {
            "beta": beta_star,
            "per_cancer_percell_vs_relay": {c: {"percell": round(cov_pc[c]["central"], 4),
                                                "relay": round(cov_rl[c]["central"], 4),
                                                "unlock": round(cov_rl[c]["central"] - cov_pc[c]["central"], 4)}
                                            for c in d33.CANCERS},
            "n_handles_unlocked": len(unlocked),
            "unlocked_handles_top": unlocked_list,
        },
        "CEILING": "beta is SWEPT globally not measured per handle (RUNG-12 AlphaFold/ESM would pin it down); "
                   "presentation_factor is a transparent rank->presentation proxy; relay ceiling inherits "
                   "RUNG-12P/B percolation-abstraction caveats; frequencies inherit RUNG-11's. Robust parts: the "
                   "relaxation factor (relay_ceiling/per_cell_bar) and the SHAPE (relay sustains coverage as beta "
                   "rises). Not per-handle truth -- an integration/what-if bridge.",
        "INTERPRETATION": "CLEAN handles are safe in both architectures. The relay's value is on tcr_dependent "
                          "handles: it converts 'safe only if the TCR is near-perfect (beta<=0.02/pf)' into "
                          "'safe if the TCR is decent (beta<=0.17/pf)'. Whether a given handle clears even the "
                          "relaxed bar is the RUNG-12 structural question -- but the relay widens the target "
                          "phase space ~8.6x.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[bridge] wrote {RESULT_JSON}")

    print(f"\n  RELAXATION: relay tolerates ~{relax_factor}x higher TCR cross-bind than per-cell gate")
    print(f"\n  per-cancer coverage at beta={beta_star} (imperfect TCR):   per-cell -> relay  (unlock)")
    for c in d33.CANCERS:
        pc, rl = cov_pc[c]["central"], cov_rl[c]["central"]
        print(f"    {c:9} {pc:6.1%} -> {rl:6.1%}   (+{rl-pc:5.1%})")
    print(f"\n  {len(unlocked)} handles unlocked by the relay at beta={beta_star} (relay-safe but NOT per-cell-safe). top:")
    for u in unlocked_list[:8]:
        print(f"    {u['handle']:30} {u['tier']:13} wt_rank={u['wt_rank']:.2f} qn={u['qn']:.2f}")

    _make_figure(result, d33.CANCERS, beta_star)
    return 0


def _make_figure(result, cancers, beta_star):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[bridge] matplotlib unavailable ({e})"); return
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    # panel 1: per-cancer per-cell vs relay coverage at beta*
    hb = result["headline_at_beta"]["per_cancer_percell_vs_relay"]
    x = np.arange(len(cancers)); w = 0.38
    pc = [hb[c]["percell"] * 100 for c in cancers]
    rl = [hb[c]["relay"] * 100 for c in cancers]
    ax[0].bar(x - w / 2, pc, w, label="per-cell TCR-T usable", color="#C1432B")
    ax[0].bar(x + w / 2, rl, w, label="RELAY-gate usable", color="#1B5E20")
    ax[0].set_xticks(x); ax[0].set_xticklabels(cancers, rotation=30, ha="right")
    ax[0].set_ylabel("% patients addressable (safely)")
    ax[0].set_title(f"Usable neoantigen addressability at imperfect TCR (beta={beta_star})\n"
                    f"relay unlocks the gap")
    ax[0].legend(fontsize=8); ax[0].grid(axis="y", alpha=0.3)

    # panel 2: coverage vs beta (the cross-bind sweep) for representative cancers
    sweep = result["beta_sweep_coverage"]
    betas = [s["beta"] for s in sweep]
    for c, col in [("PDAC", "#3B7DD8"), ("CRC", "#8E44AD"), ("MELANOMA", "#E67E22")]:
        ax[1].plot(betas, [s["relay"][c] * 100 for s in sweep], "-", color=col, label=f"{c} relay")
        ax[1].plot(betas, [s["percell"][c] * 100 for s in sweep], "--", color=col, alpha=0.6, label=f"{c} per-cell")
    rc = result["inputs"]["relay_ceiling_3d"]
    ax[1].set_xlabel("TCR mut-vs-WT cross-bind  beta  (RUNG-12 would measure this)")
    ax[1].set_ylabel("% patients usable")
    ax[1].set_title("Relay sustains usable coverage as TCR discrimination degrades\n(solid=relay, dashed=per-cell)")
    ax[1].legend(fontsize=7, ncol=3); ax[1].grid(alpha=0.3)

    fig.suptitle("RUNG-12P bridge: the gated relay unlocks RUNG-11's 'too-risky' tcr_dependent neoantigen handles",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[bridge] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # presentation_factor monotone + endpoints
    check("pf=0 for non-binder WT (rank>=2)", presentation_factor(2.5) == 0.0)
    check("pf=1 for strong WT (rank<=0.5)", presentation_factor(0.3) == 1.0)
    check("pf monotone decreasing in wt_rank", presentation_factor(0.8) > presentation_factor(1.5))

    # qn: clean -> 0 regardless of beta
    check("clean handle q_n=0 (WT not presented)", qn_of(3.0, "clean", False, 1.0) == 0.0)
    # tcr_dependent, strong WT, beta=0.5 -> q_n = 1*0.5 = 0.5
    check("tcr_dependent strong-WT q_n = pf*beta", abs(qn_of(0.5, "tcr_dependent", False, 0.5) - 0.5) < 1e-9)
    # anchor discount lowers q_n
    check("anchor discount lowers q_n", qn_of(0.5, "anchor", True, 0.5) < qn_of(0.5, "tcr_dependent", False, 0.5))

    # classify: relay set always superset of per-cell set (relay ceiling > per-cell bar)
    H = [{"gene": "K", "pos": 12, "mut_label": "G12D", "allele": "A", "tier": "tcr_dependent",
          "wt_rank": 0.6, "anchor": False, "cancer_prev": {"PDAC": 0.4}},
         {"gene": "T", "pos": 175, "mut_label": "R175H", "allele": "A", "tier": "clean",
          "wt_rank": 2.5, "anchor": False, "cancer_prev": {"PDAC": 0.1}}]
    classify(H, beta=0.1, relay_ceiling=0.173)
    pcset = {id(h) for h in H if h["usable_percell"]}
    rlset = {id(h) for h in H if h["usable_relay"]}
    check("relay-usable set superset of per-cell-usable", pcset <= rlset)
    # at beta=0.1: tcr_dependent q_n = pf(0.6)*0.1. pf(0.6)=(2-0.6)/1.5=0.933 -> 0.093 -> relay-safe(<=0.173) not per-cell(<=0.02)
    h0 = H[0]
    check("leaky handle unlocked by relay at low beta", h0["usable_relay"] and not h0["usable_percell"])
    # clean handle always safe
    check("clean handle usable in both", H[1]["usable_percell"] and H[1]["usable_relay"])
    # at high beta the leaky handle falls out of relay too
    classify(H, beta=0.9, relay_ceiling=0.173)
    check("at high beta leaky handle not even relay-safe", not H[0]["usable_relay"])

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="RUNG-12P bridge: relay unlock of RUNG-11 neoantigen handles")
    ap.add_argument("mode", nargs="?", default="run", choices=["run", "selftest"])
    args = ap.parse_args()
    sys.exit(selftest() if args.mode == "selftest" else main_run())
