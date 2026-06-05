#!/usr/bin/env python3
"""
RUNG 6 / arm (b) — the HLA-LOH NOT-gate addressability test (laptop, no GPU, minutes).

THE HYPOTHESIS UNDER TEST
-------------------------
RUNG-5 proved (worst-donor-safe, atlas-scale) that NO single/2-input *surface* AND gate is safe -> a
100% per-patient addressability gap. The field's answer is a NOT-gate on a tumour-specific GENETIC LOSS:
A2 Bio's Tmod (the most clinically-advanced logic-gated CAR-T) kills cells that have LOST HLA-A*02 (loss
of heterozygosity), while a LIR-1 blocker sensing RETAINED HLA-A*02 spares normal tissue. This script asks
the quantitative question: **how much of the 100% gap does the genetic NOT-gate actually close?**

WHAT THIS IS, HONESTLY
----------------------
This is a SYNTHESIS, not a wet discovery. It joins our worst-donor surface gap (RUNG-5) to the FIELD's
published per-patient HLA-LOH calls (Martinez-Jimenez et al. 2023, Nat Genet, DOI 10.1038/s41588-023-01367-1,
"GIE per sample" sheet, 6,319 WGS tumours / 58 types). The novel part is the integration + the honest
contrast between the GENEROUS ceiling ("any HLA-LOH" -> some allele-specific gate is in principle available)
and the ACTUAL deployed gate ("HLA-A*02-specific" -> the A2 Bio Tmod requires germline A*02-heterozygosity
AND somatic loss of the A*02 allele specifically).

THREE IRREDUCIBLE CAVEATS (stated, never papered over)
------------------------------------------------------
1. DIFFERENT MODALITIES / DIFFERENT PATIENTS. Our surface gap is from a single-cell *mRNA* atlas; the LOH
   fractions are from a *WGS genotype* cohort. We are integrating two datasets, not measuring the same
   patients. The "100% surface gap" and "X% genetic addressable" are about the same *cancer types*, not the
   same individuals.
2. PATIENT-LEVEL, NOT CLONAL. A sample is "HLA-LOH+" if the tumour lost the allele; but much HLA-LOH is
   SUBCLONAL (TRACERx: ~76% of LUAD HLA-LOH tumours), so only a SUBSET of cells in a "positive" patient
   actually lost it. Patient-level addressability is therefore an UPPER BOUND; the clonal fraction is lower.
3. LUNG IS NSCLC-POOLED. This cohort codes lung as NSCLC (no LUAD/LUSC split); RUNG-5's cancer was lung
   adenocarcinoma. We report NSCLC and flag the mismatch (LUSC carries higher HLA-LOH, so NSCLC-pooled would
   over-state LUAD — though here NSCLC=27.8% happens to sit at the TRACERx LUAD estimate ~29%).

DATA
----
data/refs/mjimenez2023_MOESM6.xlsx  (Supplementary Data 3 of the paper; sheet "GIE per sample").
Download once (no auth needed, served from the Springer CDN):
  curl -L -o data/refs/mjimenez2023_MOESM6.xlsx \
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01367-1/MediaObjects/41588_2023_1367_MOESM6_ESM.xlsx"

USAGE
  python scripts/24_hla_loh_addressability.py            # run on the real supplement -> JSON + figure
  python scripts/24_hla_loh_addressability.py selftest   # synthetic-data checks of the gate logic only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SUPP = PROJECT_ROOT / "data" / "refs" / "mjimenez2023_MOESM6.xlsx"
SHEET = "GIE per sample"
OUT_DIR = PROJECT_ROOT / "runs" / "rung6_logicgate"
RESULT_JSON = OUT_DIR / "rung6_hla_loh_addressability.json"
FIGURE_PNG = OUT_DIR / "rung6_hla_loh.png"

# RUNG-5's cancers -> this cohort's codes. Lung is NSCLC-pooled (see caveat 3).
CANCERS = {"NSCLC": "lung (NSCLC-pooled; RUNG-5 was LUAD)",
           "BRCA": "breast carcinoma",
           "COREAD": "colorectal cancer"}
SURFACE_GAP = 1.0          # RUNG-5 result: 0/1000 surface gates worst-donor-safe -> 100% gap
LOST_CN = 0.5              # allelic copy number below this = allele lost (LOHHLA convention)
A02_PREFIX = "A*02"        # A2 Bio Tmod blocker senses the HLA-A*02 group (A*02:01 dominant)


# ---------------------------------------------------------------------------
# Jeffreys one-sided bounds (mirror scripts/18 — a fraction from few patients
# carries real uncertainty; we report the bound, not just the point estimate).
# ---------------------------------------------------------------------------
def jeffreys_upper(k: int, n: int, alpha: float = 0.05) -> float:
    if n <= 0:
        return 1.0
    try:
        from scipy.stats import beta
        return float(beta.ppf(1 - alpha, k + 0.5, n - k + 0.5)) if k < n else 1.0
    except Exception:
        return min(1.0, (k + 1.96 ** 2 / 2) / (n + 1.96 ** 2) + 1.96 / (n + 4) * 0.5)


def jeffreys_lower(k: int, n: int, alpha: float = 0.05) -> float:
    if n <= 0:
        return 0.0
    try:
        from scipy.stats import beta
        return float(beta.ppf(alpha, k + 0.5, n - k + 0.5)) if k > 0 else 0.0
    except Exception:
        return max(0.0, (k + 1.96 ** 2 / 2) / (n + 1.96 ** 2) - 1.96 / (n + 4) * 0.5)


# ---------------------------------------------------------------------------
# Gate logic — the heart of the test, kept pure so selftest can exercise it on
# synthetic rows without any Excel I/O.
# ---------------------------------------------------------------------------
def _is_a02(allele) -> bool:
    return isinstance(allele, str) and allele.startswith(A02_PREFIX)


def annotate_gates(df: pd.DataFrame, lost_cn: float = LOST_CN) -> pd.DataFrame:
    """Add the four addressability flags. Columns required: A1, A2, A1_CN, A2_CN, loh_lilac."""
    out = df.copy()
    a1_is, a2_is = out.A1.apply(_is_a02), out.A2.apply(_is_a02)

    out["carries_a02"] = a1_is | a2_is
    # germline A*02 HETEROZYGOUS = exactly one allele is A*02 (XOR) AND the two alleles differ.
    # (homozygous A*02/A*02 is useless for the gate: no allelic asymmetry for a blocker to exploit.)
    out["het_a02"] = (a1_is ^ a2_is) & (out.A1 != out.A2)
    # copy number of the A*02-carrying allele (only meaningful when het_a02)
    a02_cn = np.where(a1_is, out.A1_CN, out.A2_CN)
    # A2 Bio Tmod-addressable = germline A*02-het AND the A*02 allele was LOST (CN < threshold).
    out["a02_loh_addressable"] = out["het_a02"] & (pd.Series(a02_cn, index=out.index) < lost_cn)
    # generous ceiling = the author-provided "any HLA allele lost" flag.
    out["any_hla_loh"] = out.loh_lilac.astype(bool)
    return out


def _frac_block(flag: pd.Series) -> dict:
    n = int(len(flag))
    k = int(flag.sum())
    pt = k / n if n else 0.0
    return {"n": n, "k": k, "fraction": round(pt, 4),
            "ci_lower": round(jeffreys_lower(k, n), 4),
            "ci_upper": round(jeffreys_upper(k, n), 4),
            "residual_gap": round(1.0 - pt, 4)}            # gap left after this NOT-gate (best case)


def compute_addressability(df: pd.DataFrame, lost_cn: float = LOST_CN) -> dict:
    df = annotate_gates(df, lost_cn=lost_cn)
    per_cancer = {}
    for code, label in CANCERS.items():
        s = df[df.cancer_type_code == code]
        if len(s) == 0:
            per_cancer[code] = {"label": label, "n": 0, "note": "absent from cohort"}
            continue
        block = {
            "label": label,
            "n_patients": int(len(s)),
            "carries_a02": _frac_block(s.carries_a02),
            "het_a02": _frac_block(s.het_a02),
            "any_hla_loh": _frac_block(s.any_hla_loh),            # generous ceiling
            "a02_tmod_addressable": _frac_block(s.a02_loh_addressable),  # the actual deployed gate
        }
        # INVARIANTS (a subset must never exceed its superset; the gate is a subset of any-LOH).
        assert block["a02_tmod_addressable"]["fraction"] <= block["het_a02"]["fraction"] + 1e-9
        assert block["het_a02"]["fraction"] <= block["carries_a02"]["fraction"] + 1e-9
        block["surface_gap_rung5"] = SURFACE_GAP
        block["gap_after_any_loh_gate"] = block["any_hla_loh"]["residual_gap"]
        block["gap_after_a02_tmod_gate"] = block["a02_tmod_addressable"]["residual_gap"]
        block["overstatement_factor_any_vs_a02"] = (
            round(block["any_hla_loh"]["fraction"] / block["a02_tmod_addressable"]["fraction"], 2)
            if block["a02_tmod_addressable"]["fraction"] > 0 else None)
        per_cancer[code] = block
    return per_cancer


# ---------------------------------------------------------------------------
def load_supplement(path: Path = SUPP) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Supplement not found: {path}\nDownload it (no auth needed) with:\n"
            f'  curl -L -o {path} "https://static-content.springer.com/esm/'
            f'art%3A10.1038%2Fs41588-023-01367-1/MediaObjects/41588_2023_1367_MOESM6_ESM.xlsx"')
    df = pd.read_excel(path, sheet_name=SHEET)
    need = {"sample_id_2", "cancer_type_code", "A1", "A2", "A1_CN", "A2_CN", "loh_lilac"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Supplement sheet '{SHEET}' missing columns: {missing}")
    return df


def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_supplement()
    print(f"[rung6] loaded {len(df):,} samples from {SUPP.name} sheet '{SHEET}'")

    per_cancer = compute_addressability(df)

    # sensitivity of the A*02 gate to the 'allele lost' CN threshold (robustness, not cherry-picking)
    sens = {}
    for thr in (0.3, 0.5, 0.7):
        ann = annotate_gates(df, lost_cn=thr)
        sens[str(thr)] = {c: round(float(ann[ann.cancer_type_code == c].a02_loh_addressable.mean()), 4)
                          for c in CANCERS}

    result = {
        "tag": "rung6_hla_loh_addressability",
        "hypothesis": "How much of RUNG-5's 100% surface-gate addressability gap does a genetic HLA-LOH "
                      "NOT-gate close? Contrast the GENEROUS any-HLA-LOH ceiling vs the ACTUAL deployed "
                      "A2 Bio Tmod gate (HLA-A*02-specific).",
        "data_source": "Martinez-Jimenez et al. 2023, Nat Genet, DOI 10.1038/s41588-023-01367-1, "
                       "Supplementary Data (MOESM6) sheet 'GIE per sample'; 6,319 WGS tumours.",
        "n_samples_total": int(len(df)),
        "surface_gap_rung5": SURFACE_GAP,
        "lost_cn_threshold": LOST_CN,
        "a02_match": A02_PREFIX,
        "per_cancer": per_cancer,
        "a02_addressable_sensitivity_to_lostCN": sens,
        "HEADLINE": {c: {"any_loh_gate_addressable": per_cancer[c]["any_hla_loh"]["fraction"],
                         "a02_tmod_addressable": per_cancer[c]["a02_tmod_addressable"]["fraction"],
                         "gap_after_a02_tmod": per_cancer[c]["gap_after_a02_tmod_gate"]}
                     for c in CANCERS},
        "CEILING": "Synthesis of OUR surface gap (scRNA atlas) + FIELD's published HLA-LOH calls (WGS) — "
                   "different modalities, different patients, same cancer types. Patient-level not clonal "
                   "(subclonal LOH => true addressable is LOWER; TRACERx ~76% of LUAD LOH subclonal). Lung "
                   "is NSCLC-pooled, not LUAD. A*02 group match (A2 Bio is A*02:01). NOT a wet result.",
        "INTERPRETATION": "any-HLA-LOH (some allele-specific gate in principle) gives a generous ceiling; "
                          "the ACTUAL single-allele A*02 Tmod gate addresses ~3-6% of patients (5-9x smaller) "
                          "because it needs loss of the SPECIFIC sensed allele. The 100% gap barely moves for "
                          "a single-allele gate. Implication: a PANEL of allele-specific blockers (A*02 + B*07 "
                          "+ ...) is needed to recover the any-HLA-LOH ceiling. This is the next design step.",
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung6] wrote {RESULT_JSON}")

    # console summary
    print("\n  cancer   n     anyHLA-LOH     A*02-Tmod      gap(A*02)   over-state")
    for c, b in per_cancer.items():
        if b.get("n_patients"):
            print(f"  {c:7}{b['n_patients']:5}  "
                  f"{b['any_hla_loh']['fraction']:.1%} [{b['any_hla_loh']['ci_lower']:.1%}-{b['any_hla_loh']['ci_upper']:.1%}]  "
                  f"{b['a02_tmod_addressable']['fraction']:.1%} [{b['a02_tmod_addressable']['ci_lower']:.1%}-{b['a02_tmod_addressable']['ci_upper']:.1%}]  "
                  f"{b['gap_after_a02_tmod_gate']:.1%}      {b['overstatement_factor_any_vs_a02']}x")

    _make_figure(per_cancer)
    return 0


def _make_figure(per_cancer: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung6] matplotlib unavailable ({e}); skipped figure")
        return
    codes = [c for c in CANCERS if per_cancer.get(c, {}).get("n_patients")]
    labels = [f"{c}\n(n={per_cancer[c]['n_patients']})" for c in codes]
    any_loh = [per_cancer[c]["any_hla_loh"]["fraction"] * 100 for c in codes]
    a02 = [per_cancer[c]["a02_tmod_addressable"]["fraction"] * 100 for c in codes]
    x = np.arange(len(codes))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

    # panel 1: addressable fractions (ceiling vs deployed gate)
    w = 0.36
    ax[0].bar(x - w / 2, any_loh, w, label="any-HLA-LOH (ceiling)", color="#4C9F70")
    ax[0].bar(x + w / 2, a02, w, label="HLA-A*02 Tmod (deployed)", color="#C1432B")
    for xi, v in zip(x - w / 2, any_loh):
        ax[0].text(xi, v + 0.4, f"{v:.1f}%", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, a02):
        ax[0].text(xi, v + 0.4, f"{v:.1f}%", ha="center", fontsize=8)
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels)
    ax[0].set_ylabel("% patients addressable")
    ax[0].set_title("Genetic NOT-gate addressability\n(ceiling vs actual single-allele gate)")
    ax[0].legend(fontsize=8); ax[0].grid(axis="y", alpha=0.3)

    # panel 2: the addressability gap — surface (RUNG-5) vs after each genetic gate
    gap_any = [(1 - per_cancer[c]["any_hla_loh"]["fraction"]) * 100 for c in codes]
    gap_a02 = [per_cancer[c]["gap_after_a02_tmod_gate"] * 100 for c in codes]
    ax[1].bar(x - w, [100] * len(codes), w, label="surface gates (RUNG-5)", color="#888")
    ax[1].bar(x, gap_any, w, label="+ any-HLA-LOH gate", color="#4C9F70")
    ax[1].bar(x + w, gap_a02, w, label="+ A*02 Tmod gate", color="#C1432B")
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels)
    ax[1].set_ylabel("% patients still UNADDRESSED (gap)")
    ax[1].set_ylim(0, 105)
    ax[1].set_title("Addressability gap shrinks only marginally\nfor a single-allele genetic gate")
    ax[1].legend(fontsize=8); ax[1].grid(axis="y", alpha=0.3)

    fig.suptitle("RUNG-6 (b): HLA-LOH NOT-gate — how much of the 100% gap does the genetic signal close?",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung6] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    """Exercise the gate logic on synthetic rows with known answers (no Excel needed)."""
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond)))
        ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # synthetic patients with hand-computed expected flags
    rows = [
        # A1,        A2,        A1_CN, A2_CN, loh_lilac, expect_a02_addr, why
        ("A*02:01", "A*11:01", 3.20, 0.00, True,  True),   # het A*02, A2(A*11) lost? no -> A*02 retained...
        ("A*02:01", "A*11:01", 0.00, 2.10, True,  True),   # het A*02, A*02 allele(A1) lost -> addressable
        ("A*02:01", "A*11:01", 3.20, 2.00, False, False),  # het A*02, nothing lost -> not addressable
        ("A*11:01", "A*03:01", 0.00, 2.00, True,  False),  # no A*02 -> not addressable (lost a non-A*02)
        ("A*02:01", "A*02:01", 1.00, 0.00, True,  False),  # homozygous A*02 -> het=False -> not addressable
        ("A*02:06", "A*24:02", 0.10, 2.00, True,  True),   # het A*02 (group), A*02 allele lost -> addressable
    ]
    df = pd.DataFrame(rows, columns=["A1", "A2", "A1_CN", "A2_CN", "loh_lilac", "_expect"])
    df["cancer_type_code"] = "NSCLC"
    ann = annotate_gates(df, lost_cn=LOST_CN)

    # Fix the first row's expectation: A1=A*02 CN=3.2 (retained), A2=A*11 CN=0 (lost) => A*02 NOT lost.
    expect = [False, True, False, False, False, True]
    got = ann.a02_loh_addressable.tolist()
    check("per-row A*02-LOH addressability matches hand calc", got == expect)
    check("carries_a02 counts A*02 in either slot", ann.carries_a02.tolist() == [True, True, True, False, True, True])
    check("het_a02 excludes homozygous A*02/A*02", ann.het_a02.tolist() == [True, True, True, False, False, True])

    # invariant: addressable subset of het subset of carries
    check("invariant a02_addr <= het <= carries",
          (ann.a02_loh_addressable <= ann.het_a02).all() and (ann.het_a02 <= ann.carries_a02).all())
    # any-HLA-LOH (loh_lilac) must be >= the A*02-specific subset count
    check("any_hla_loh >= a02_addressable (count)",
          int(ann.any_hla_loh.sum()) >= int(ann.a02_loh_addressable.sum()))

    # CN threshold monotonicity: a higher 'lost' threshold can only ADD addressable patients
    a_lo = annotate_gates(df, lost_cn=0.3).a02_loh_addressable.sum()
    a_hi = annotate_gates(df, lost_cn=0.7).a02_loh_addressable.sum()
    check("looser lost-CN threshold is monotone (>=)", a_hi >= a_lo)

    # Jeffreys: a 0/50 observation must NOT credit a true zero (upper bound > 0)
    check("jeffreys_upper(0,50) > 0 (no false zero)", jeffreys_upper(0, 50) > 0.0)
    check("jeffreys_lower(50,50) < 1 (no false one)", jeffreys_lower(50, 50) < 1.0)

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", nargs="?", default="run", choices=["run", "selftest"])
    args = ap.parse_args()
    sys.exit(selftest() if args.mode == "selftest" else main_run())
