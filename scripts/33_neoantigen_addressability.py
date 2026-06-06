#!/usr/bin/env python3
"""
RUNG 11 — Public-neoantigen Addressability x Discriminability map (laptop, no GPU, minutes).

WHY THIS RUNG EXISTS (the pivot)
--------------------------------
RUNG 5->10b exhausted ONE recognition axis: "which genes are ON/OFF" (the EXPRESSION axis, on the
scRNA atlas). We proved it is bounded -- surface logic gates can't be worst-donor-safe, no surface
NOT-blocker even exists, and the genetic HLA-LOH NOT-gate only reaches ~15-28% of patients. That map
is real, but it is ONE axis.

This rung opens the SECOND axis -- the one AlphaFold/ESM live on: "is this protein MUTATED?" A somatic
mutation is, BY CONSTRUCTION, absent from every healthy cell in the body. It is the tumour-EXCLUSIVE
signal the expression axis could never produce. Concretely: PUBLIC neoantigens -- peptides from
RECURRENT driver hotspot mutations (KRAS G12D/V/C, TP53 R175H, PIK3CA H1047R, BRAF V600E, ...) presented
on common HLA-I alleles. These are a live clinical frontier (KRAS-G12D/C*08:02 drove metastatic CRC
regression; KRAS-G12V/A*11:01; TP53-R175H/A*02:01; CTNNB1; PIK3CA).

THE TWO QUESTIONS THIS SCRIPT ANSWERS EMPIRICALLY
-------------------------------------------------
1. ADDRESSABILITY -- what fraction of patients (per cancer type) carry a driver hotspot whose mutant
   peptide is PRESENTED on an HLA allele they actually carry? = sum over (mutation x allele) of
   P(has mutation in this cancer) x P(carries the restricting allele). Compared head-to-head with the
   HLA-LOH ceiling (15-28%) and the deployed single-allele A*02 gate (3-6%) from RUNG-6.
2. DISCRIMINABILITY (the safety axis) -- for each presented mutant peptide, does its WILD-TYPE
   counterpart ALSO get presented? If WT is NOT presented, the WT pMHC is simply absent from healthy
   cells and cross-reactivity is structurally impossible at the pMHC level ("clean"). If both are
   presented, safety rests on the TCR telling two near-identical surfaces apart -- the MAGE-A3 failure
   mode ("tcr_dependent", routes to the RUNG-12 AlphaFold/ESM structural check). Anchor-residue
   mutations (P2 / C-terminus) flip MHC binding itself, giving MHC-level discrimination ("anchor").

WHAT THIS IS, HONESTLY
----------------------
A PREDICTION-driven SYNTHESIS, not a wet discovery. It joins three public ingredients: (a) canonical
protein sequences (UniProt), (b) MHC-I presentation prediction (MHCflurry 2.0 -- affinity + processing),
(c) literature driver-mutation frequencies + HLA carrier frequencies. The novelty is the integration and
the explicit addressability-x-discriminability contrast against the expression-axis ceiling. Every hit
is a HYPOTHESIS.

FIVE IRREDUCIBLE CAVEATS (stated, never papered over)
-----------------------------------------------------
1. PREDICTED, NOT MEASURED. MHCflurry presentation_score != actual surface presentation. Ground truth =
   immunopeptidomics (MS). A presented-by-prediction peptide may not be on the real cell surface.
2. PROCESSING IS MODELLED COARSELY. Proteasomal cleavage / TAP transport are captured by MHCflurry's
   processing predictor, but real processing is context- and proteasome-type-dependent.
3. TCR-EXISTENCE RESIDUAL. A presented mutant pMHC is NECESSARY but NOT SUFFICIENT -- a TCR (or TCR-mimic
   antibody) with the needed mut-vs-WT discrimination must EXIST or be engineerable. That is the wet-lab /
   structure residual; RUNG-12 (AlphaFold-Multimer / ESM) probes it. The "tcr_dependent" tier flags
   exactly where this residual dominates.
4. FREQUENCIES ARE LITERATURE POINT ESTIMATES, POPULATION-DEPENDENT. Mutation prevalences (COSMIC/TCGA/
   cBioPortal) and HLA carrier frequencies (AFND/NMDP, here European-weighted) vary by study and ancestry.
   East/South-Asian populations raise A*11:01 / A*24:02 substantially -> KRAS coverage is HIGHER there.
   Coverage is an estimate joining independent datasets, NOT a measured same-patient cohort intersection.
5. CLASS I / CD8 ONLY. Class-II / CD4 neoantigens (also therapeutically real) are out of scope here.

USAGE
  python scripts/33_neoantigen_addressability.py            # real run (needs mhcflurry) -> JSON + figure
  python scripts/33_neoantigen_addressability.py selftest   # offline logic checks (mock predictor; no deps)

REAL RUN SETUP (one time, laptop, ~2 min)
  pip install mhcflurry pandas numpy matplotlib scipy
  mhcflurry-downloads fetch models_class1_presentation
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REFS = PROJECT_ROOT / "data" / "refs"
OUT_DIR = PROJECT_ROOT / "runs" / "rung11_neoantigen"
RESULT_JSON = OUT_DIR / "rung11_neoantigen_addressability.json"
FIGURE_PNG = OUT_DIR / "rung11_neoantigen.png"

# binder/presented thresholds (MHCflurry affinity percentile rank; standard field cutoffs)
BINDER_RANK = 2.0          # affinity_percentile <= 2  -> "presented/binder"
STRONG_RANK = 0.5          # <= 0.5 -> "strong binder"
CLEAR_NONBINDER = 4.0      # WT %rank > this = CLEARLY not presented (robust "clean" call; >2 alone is fragile)
PEPTIDE_LENS = (8, 9, 10, 11)

# HLA-LOH benchmarks from RUNG-6 (the expression-axis ceiling we are trying to beat)
HLA_LOH_CEILING = (0.135, 0.278)   # any-HLA-LOH availability, breast..NSCLC (patient-level UPPER bound)
A02_TMOD_DEPLOYED = (0.031, 0.057) # the actual single-allele A*02 Tmod gate, breast..colorectal

# ---------------------------------------------------------------------------
# DRIVER HOTSPOTS -- recurrent, "public" mutations. position is 1-based residue in the UniProt canonical
# sequence; wt/mut are single-letter. cancer_prev = fraction of ALL patients in that cancer type carrying
# THIS exact mutation (literature point estimates; COSMIC / TCGA pan-cancer / cBioPortal + KRAS reviews).
# Mutations sharing (gene,position) are MUTUALLY EXCLUSIVE at that codon (handled in coverage()).
# ---------------------------------------------------------------------------
DRIVERS = [
    # gene,   acc,       pos,  wt,  mut, cancer_prev
    ("KRAS",  "P01116",   12, "G", "D", {"PDAC": 0.40, "CRC": 0.13, "NSCLC": 0.04, "OV": 0.02}),
    ("KRAS",  "P01116",   12, "G", "V", {"PDAC": 0.22, "CRC": 0.09, "NSCLC": 0.03}),
    ("KRAS",  "P01116",   12, "G", "C", {"NSCLC": 0.13, "CRC": 0.03, "PDAC": 0.01}),
    ("KRAS",  "P01116",   12, "G", "R", {"PDAC": 0.16}),
    ("TP53",  "P04637",  175, "R", "H", {"OV": 0.03, "PDAC": 0.02, "CRC": 0.02, "BRCA": 0.02, "NSCLC": 0.02, "HNSC": 0.03}),
    ("TP53",  "P04637",  248, "R", "Q", {"OV": 0.02, "PDAC": 0.02, "CRC": 0.02, "BRCA": 0.015}),
    ("TP53",  "P04637",  273, "R", "H", {"OV": 0.02, "CRC": 0.02, "NSCLC": 0.02}),
    ("TP53",  "P04637",  220, "Y", "C", {"OV": 0.02, "BRCA": 0.015, "NSCLC": 0.015}),
    ("PIK3CA","P42336", 1047, "H", "R", {"BRCA": 0.13, "CRC": 0.04, "HNSC": 0.05}),
    ("PIK3CA","P42336",  545, "E", "K", {"BRCA": 0.08, "CRC": 0.05, "HNSC": 0.04}),
    ("BRAF",  "P15056",  600, "V", "E", {"MELANOMA": 0.45, "CRC": 0.10, "NSCLC": 0.02}),
    ("EGFR",  "P00533",  858, "L", "R", {"NSCLC": 0.11}),
    ("IDH1",  "O75874",  132, "R", "H", {"GLIOMA": 0.70, "HCC": 0.01}),
    ("CTNNB1","P35222",   37, "S", "F", {"HCC": 0.06}),
    ("CTNNB1","P35222",   45, "S", "P", {"HCC": 0.04}),
    ("NRAS",  "P01111",   61, "Q", "R", {"MELANOMA": 0.13}),
    ("NRAS",  "P01111",   61, "Q", "K", {"MELANOMA": 0.07}),
]
CANCERS = ["PDAC", "CRC", "NSCLC", "BRCA", "MELANOMA", "GLIOMA", "HCC", "OV"]

# common HLA-I alleles with CARRIER (phenotype) frequency = P(individual has >=1 copy). European-weighted
# (AFND / NMDP); flagged population-dependent. East/South-Asian shift A*11:01/A*24:02 up (helps KRAS).
HLA_PANEL = {
    "HLA-A*02:01": 0.43, "HLA-A*01:01": 0.24, "HLA-A*03:01": 0.22, "HLA-A*11:01": 0.20,
    "HLA-A*24:02": 0.22, "HLA-A*26:01": 0.08, "HLA-A*32:01": 0.08, "HLA-A*68:01": 0.10,
    "HLA-B*07:02": 0.18, "HLA-B*08:01": 0.14, "HLA-B*35:01": 0.12, "HLA-B*44:02": 0.12,
    "HLA-B*44:03": 0.10, "HLA-B*15:01": 0.10, "HLA-B*40:01": 0.10,
    "HLA-C*07:01": 0.26, "HLA-C*07:02": 0.24, "HLA-C*08:02": 0.06, "HLA-C*06:02": 0.16, "HLA-C*05:01": 0.12,
}

# the oracle: clinically-validated public-neoantigen handles the pipeline MUST recover as "presented".
KNOWN_POSITIVES = [("KRAS", "G12D", "HLA-C*08:02"),
                   ("KRAS", "G12V", "HLA-A*11:01"),
                   ("TP53", "R175H", "HLA-A*02:01")]
# known-DIFFICULT (hydrophobic, weak class-I presenter) -> a sanity negative, not expected to score high.
KNOWN_HARD = [("BRAF", "V600E", "HLA-A*02:01")]


# ---------------------------------------------------------------------------
# Jeffreys bounds (mirror scripts/24 -- a fraction from finite estimates carries uncertainty).
# ---------------------------------------------------------------------------
def jeffreys_upper(k: int, n: int, alpha: float = 0.05) -> float:
    if n <= 0:
        return 1.0
    try:
        from scipy.stats import beta
        return float(beta.ppf(1 - alpha, k + 0.5, n - k + 0.5)) if k < n else 1.0
    except Exception:
        return min(1.0, (k + 1.92) / (n + 3.84) + 0.98 / (n + 4))


# ---------------------------------------------------------------------------
# Sequence + peptide logic (pure -- selftest exercises it with a mock sequence, no network).
# ---------------------------------------------------------------------------
def fetch_uniprot(acc: str) -> str:
    """Canonical UniProt sequence, cached to data/refs/uniprot_<acc>.fasta. (network once per accession)"""
    cache = REFS / f"uniprot_{acc}.fasta"
    if cache.exists():
        txt = cache.read_text()
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        try:
            with urllib.request.urlopen(url, timeout=30) as r:    # noqa: S310 (trusted host)
                txt = r.read().decode()
        except Exception:
            # macOS Python often lacks a CA bundle -> fall back to curl (system cert store)
            import subprocess
            txt = subprocess.check_output(["curl", "-fsSL", url], text=True, timeout=60)
        REFS.mkdir(parents=True, exist_ok=True)
        cache.write_text(txt)
    seq = "".join(l.strip() for l in txt.splitlines() if l and not l.startswith(">"))
    if not seq:
        raise ValueError(f"empty sequence for {acc}")
    return seq


def gen_registers(seq: str, pos: int, wt: str, mut: str, lens=PEPTIDE_LENS) -> list[dict]:
    """All k-mers (k in lens) that CONTAIN the mutated residue. Returns paired WT/MUT peptides + the
    mutated residue's index/anchor status within each register. pos is 1-based."""
    i = pos - 1
    if not (0 <= i < len(seq)):
        raise IndexError(f"position {pos} out of range for length {len(seq)}")
    if seq[i] != wt:
        raise AssertionError(f"UniProt residue {pos} is '{seq[i]}', expected WT '{wt}'")
    mut_seq = seq[:i] + mut + seq[i + 1:]
    out = []
    for L in lens:
        for start in range(max(0, i - L + 1), min(i, len(seq) - L) + 1):
            j = i - start                                     # 0-based index of mutated residue in peptide
            if not (0 <= j < L):
                continue
            pep_wt = seq[start:start + L]
            pep_mut = mut_seq[start:start + L]
            if len(pep_mut) != L or pep_wt == pep_mut:
                continue
            anchor = (j == 1) or (j == L - 1)                 # P2 or C-terminus = primary HLA-I anchors
            out.append({"len": L, "j": j, "p_in_pep": j + 1, "anchor": anchor,
                        "pep_wt": pep_wt, "pep_mut": pep_mut})
    return out


def classify_tier(mut_presented: bool, wt_presented: bool, anchor: bool) -> str:
    """The discriminability (safety) tier of a presented handle."""
    if not mut_presented:
        return "not_presented"
    if not wt_presented:
        return "clean"            # WT pMHC absent from healthy surface -> pMHC-level cross-reactivity impossible
    if anchor:
        return "anchor"           # mutation at MHC anchor -> MHC-level discrimination assists the TCR
    return "tcr_dependent"        # both pMHC on surfaces -> safety rests entirely on the TCR (MAGE-A3 risk class)


TIER_ORDER = {"clean": 0, "anchor": 1, "tcr_dependent": 2, "not_presented": 3}
SAFE_TIERS = ("clean",)
BROAD_TIERS = ("clean", "anchor", "tcr_dependent")


def acc_tiers(*tiers):
    """Predicate: handle's tier is in `tiers`."""
    return lambda h: h["tier"] in tiers


def acc_clean_robust(h) -> bool:
    """Robust SAFE: clean AND the WT is CLEARLY a non-binder (not just barely over rank 2)."""
    return (h["tier"] == "clean" and h.get("wt_rank", 0.0) > CLEAR_NONBINDER
            and h.get("mut_rank", 99.0) <= BINDER_RANK)


def evaluate_driver(label: str, registers: list[dict], alleles: list[str], scores: dict) -> dict:
    """For each allele, pick the best-presented MUT register and classify its discriminability.
    `scores` maps (peptide, allele) -> {"rank": affinity_percentile, "pres": presentation_score}."""
    per_allele = {}
    for a in alleles:
        best = None
        for reg in registers:
            sm = scores.get((reg["pep_mut"], a))
            if sm is None:
                continue
            if best is None or sm["rank"] < best["mut_rank"]:
                sw = scores.get((reg["pep_wt"], a), {"rank": 100.0, "pres": 0.0})
                mut_presented = sm["rank"] <= BINDER_RANK
                wt_presented = sw["rank"] <= BINDER_RANK
                best = {"allele": a, "len": reg["len"], "p_in_pep": reg["p_in_pep"], "anchor": reg["anchor"],
                        "pep_mut": reg["pep_mut"], "pep_wt": reg["pep_wt"],
                        "mut_rank": round(sm["rank"], 3), "wt_rank": round(sw["rank"], 3),
                        "mut_pres": round(sm.get("pres", 0.0), 3), "wt_pres": round(sw.get("pres", 0.0), 3),
                        "mut_presented": mut_presented, "wt_presented": wt_presented,
                        "strong": sm["rank"] <= STRONG_RANK}
        if best is None:
            continue
        best["tier"] = classify_tier(best["mut_presented"], best["wt_presented"], best["anchor"])
        per_allele[a] = best
    return per_allele


# ---------------------------------------------------------------------------
# Addressability / coverage. Per cancer: a patient benefits from mutation m iff they carry m AND >=1 allele
# that presents m (in the included tiers). Mutations at the same (gene,codon) are MUTUALLY EXCLUSIVE -> sum;
# different codons/genes treated INDEPENDENT -> 1 - product. Reports a [central, naive-sum-upper] range.
# ---------------------------------------------------------------------------
def coverage(handles: list[dict], hla_freq: dict, accept) -> dict:
    """handles: list of {gene, pos, mut_label, allele, tier, cancer_prev, ...}. `accept` is a predicate
    selecting which handles count. Returns per-cancer coverage."""
    if isinstance(accept, (tuple, list)):     # convenience: a tier tuple
        accept = acc_tiers(*accept)
    out = {}
    for cancer in CANCERS:
        # per-mutation P(presented by >=1 carried allele) over accepted handles
        by_mut = {}   # (gene,pos,mut_label) -> {"p_mut": x, "alleles": set}
        for h in handles:
            if not accept(h):
                continue
            p_mut = h["cancer_prev"].get(cancer, 0.0)
            if p_mut <= 0:
                continue
            key = (h["gene"], h["pos"], h["mut_label"])
            rec = by_mut.setdefault(key, {"p_mut": p_mut, "alleles": set()})
            rec["alleles"].add(h["allele"])
        # group mutually-exclusive codon variants; combine groups independently
        groups = {}   # (gene,pos) -> benefit (summed over its mut variants)
        naive_sum = 0.0
        for (gene, posn, _mut), rec in by_mut.items():
            p_present = 1.0 - np.prod([1.0 - hla_freq[a] for a in rec["alleles"]])
            benefit = rec["p_mut"] * float(p_present)
            groups[(gene, posn)] = groups.get((gene, posn), 0.0) + benefit
            naive_sum += benefit
        central = 1.0 - float(np.prod([1.0 - b for b in groups.values()])) if groups else 0.0
        out[cancer] = {"central": round(central, 4),
                       "upper_naive_sum": round(min(1.0, naive_sum), 4),
                       "n_mutations_targetable": len(by_mut)}
    return out


# ---------------------------------------------------------------------------
# Real MHCflurry predictor (per (peptide, allele) presentation). Selftest never calls this.
# ---------------------------------------------------------------------------
def mhcflurry_scores(peptides: list[str], alleles: list[str]) -> dict:
    try:
        from mhcflurry import Class1PresentationPredictor
    except Exception as e:
        raise RuntimeError(
            "mhcflurry not installed. Run:\n"
            "  pip install mhcflurry\n"
            "  mhcflurry-downloads fetch models_class1_presentation\n"
            f"(import error: {e})")
    predictor = Class1PresentationPredictor.load()
    peptides = sorted({p for p in peptides if 8 <= len(p) <= 15})
    scores, supported = {}, []
    for a in alleles:
        try:
            df = predictor.predict(peptides=peptides, alleles={a: [a]},
                                   include_affinity_percentile=True, verbose=0)
        except Exception as e:
            print(f"[rung11]   allele {a} unsupported by MHCflurry, skipped ({type(e).__name__})")
            continue
        supported.append(a)
        for _, row in df.iterrows():
            scores[(row["peptide"], a)] = {"rank": float(row["affinity_percentile"]),
                                           "pres": float(row.get("presentation_score", 0.0))}
    return scores, supported


# ---------------------------------------------------------------------------
def main_run() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[rung11] {len(DRIVERS)} driver hotspots x {len(HLA_PANEL)} HLA alleles; fetching UniProt seqs...")

    # 1) sequences + registers
    seq_cache, reg_by_driver, all_peptides = {}, {}, set()
    for gene, acc, pos, wt, mut, _prev in DRIVERS:
        if acc not in seq_cache:
            seq_cache[acc] = fetch_uniprot(acc)
        regs = gen_registers(seq_cache[acc], pos, wt, mut)
        reg_by_driver[(gene, pos, mut)] = regs
        for r in regs:
            all_peptides.add(r["pep_mut"]); all_peptides.add(r["pep_wt"])
    print(f"[rung11] generated {len(all_peptides):,} unique peptides "
          f"(mut+WT, {min(PEPTIDE_LENS)}-{max(PEPTIDE_LENS)}mers spanning each hotspot)")

    # 2) presentation prediction (MHCflurry)
    alleles = list(HLA_PANEL)
    try:
        scores, supported = mhcflurry_scores(sorted(all_peptides), alleles)
    except RuntimeError as e:
        print(f"[rung11] {e}")
        return 3
    print(f"[rung11] scored {len(scores):,} (peptide,allele) pairs over {len(supported)}/{len(alleles)} alleles")

    # 3) per-driver evaluation -> flat handle list
    handles, per_driver = [], {}
    for gene, acc, pos, wt, mut, prev in DRIVERS:
        mut_label = f"{wt}{pos}{mut}"
        per_allele = evaluate_driver(mut_label, reg_by_driver[(gene, pos, mut)], supported, scores)
        per_driver[f"{gene} {mut_label}"] = {
            "gene": gene, "mutation": mut_label, "cancer_prev": prev,
            "by_allele": {a: {k: b[k] for k in ("tier", "len", "p_in_pep", "anchor",
                                                "mut_rank", "wt_rank", "mut_presented", "wt_presented",
                                                "pep_mut", "pep_wt", "strong")}
                          for a, b in per_allele.items() if b["tier"] != "not_presented"},
        }
        for a, b in per_allele.items():
            if b["tier"] == "not_presented":
                continue
            handles.append({"gene": gene, "pos": pos, "mut_label": mut_label, "allele": a,
                            "tier": b["tier"], "cancer_prev": prev, "anchor": b["anchor"],
                            "mut_rank": b["mut_rank"], "wt_rank": b["wt_rank"], "strong": b["strong"]})

    # 4) addressability. SAFE reported as a RANGE because the clean call is threshold-sensitive:
    #    LOWER = robust clean (WT clearly non-binding, rank > CLEAR_NONBINDER);
    #    UPPER = lenient clean (WT just over the binder cutoff, rank > 2). BROAD = all presented handles.
    cov_safe_strict = coverage(handles, HLA_PANEL, acc_clean_robust)
    cov_safe = coverage(handles, HLA_PANEL, SAFE_TIERS)
    cov_broad = coverage(handles, HLA_PANEL, BROAD_TIERS)

    # GOLD handles: the most defensible targets -- mutant is a STRONG binder AND WT is CLEARLY off the MHC.
    gold = sorted(
        [{"driver": f"{h['gene']} {h['mut_label']}", "allele": h["allele"], "anchor": h["anchor"],
          "mut_rank": h["mut_rank"], "wt_rank": h["wt_rank"],
          "cancers": {c: p for c, p in h["cancer_prev"].items()}}
         for h in handles if h["strong"] and h["tier"] == "clean" and h["wt_rank"] > CLEAR_NONBINDER],
        key=lambda g: (g["mut_rank"], -g["wt_rank"]))

    # 5) oracle validation -- did we recover the known clinical handles?
    def _recovered(gene, mut_label, allele):
        d = per_driver.get(f"{gene} {mut_label}", {}).get("by_allele", {}).get(allele)
        return None if d is None else {"presented": d["mut_presented"], "tier": d["tier"],
                                       "mut_rank": d["mut_rank"], "peptide": d["pep_mut"]}
    validation = {f"{g} {m} / {a}": _recovered(g, m, a) for g, m, a in KNOWN_POSITIVES}
    hard_check = {f"{g} {m} / {a}": _recovered(g, m, a) for g, m, a in KNOWN_HARD}
    n_recovered = sum(1 for v in validation.values() if v and v["presented"])

    # tier census + fragility audit (clean handles whose WT sits in the marginal [2, CLEAR_NONBINDER) band)
    tier_counts = {}
    for h in handles:
        tier_counts[h["tier"]] = tier_counts.get(h["tier"], 0) + 1
    n_clean = sum(1 for h in handles if h["tier"] == "clean")
    n_clean_fragile = sum(1 for h in handles if h["tier"] == "clean" and h["wt_rank"] <= CLEAR_NONBINDER)
    n_gold = len(gold)

    result = {
        "tag": "rung11_neoantigen_addressability",
        "axis": "SEQUENCE/neoantigen (is this protein MUTATED?) -- the tumour-EXCLUSIVE recognition axis, "
                "distinct from the EXPRESSION axis (RUNG 5-10b) which proved bounded.",
        "question": "Per cancer type, what fraction of patients carry a driver hotspot whose MUTANT peptide "
                    "is presented on an HLA allele they carry (ADDRESSABILITY), and how cleanly can the "
                    "mutant be discriminated from wild-type at the pMHC level (DISCRIMINABILITY/safety)?",
        "method": {
            "presentation_predictor": "MHCflurry 2.0 Class1PresentationPredictor (affinity %rank + processing)",
            "binder_rank": BINDER_RANK, "strong_rank": STRONG_RANK, "peptide_lengths": list(PEPTIDE_LENS),
            "sequences": "UniProt canonical, validated wt residue == sequence[pos-1]",
            "anchor_positions": "P2 and C-terminus (primary HLA-I anchors)",
        },
        "n_drivers": len(DRIVERS), "n_alleles_panel": len(HLA_PANEL), "n_alleles_supported": len(supported),
        "n_handles_presented": len(handles), "tier_counts": tier_counts,
        "clean_fragility": {"n_clean": n_clean, "n_clean_fragile_wt_2to4": n_clean_fragile,
                            "n_gold_strong_mut_clear_wt": n_gold,
                            "note": "a 'clean' call (WT not presented) flips to tcr_dependent if the WT %rank "
                                    "cutoff tightens; fragile = WT in [2,4). GOLD = strong mutant binder "
                                    "(rank<=0.5) AND WT clearly off MHC (rank>4) -> the defensible targets."},
        "benchmarks_expression_axis": {
            "hla_loh_any_ceiling": HLA_LOH_CEILING, "a02_tmod_deployed": A02_TMOD_DEPLOYED,
            "note": "RUNG-6 expression-axis numbers we are trying to beat (patient-level upper bounds).",
        },
        "addressability_SAFE_range": {"strict_clean_robust": cov_safe_strict, "lenient_clean": cov_safe,
                                      "meaning": "SAFE = WT pMHC absent from healthy cells (cross-reactivity "
                                                 "structurally impossible at pMHC level). Reported strict..lenient."},
        "addressability_BROAD_all_presented": cov_broad,
        "gold_handles_strong_mut_clear_wt": gold,
        "oracle_validation_known_positives": validation,
        "oracle_validation_n_recovered": f"{n_recovered}/{len(KNOWN_POSITIVES)}",
        "sanity_known_hard": hard_check,
        "per_driver": per_driver,
        "CAVEATS": "PREDICTED not measured (MHCflurry != immunopeptidomics); processing modelled coarsely; "
                   "TCR-existence is the wet-lab residual (tcr_dependent tier = where it dominates, -> RUNG-12 "
                   "AlphaFold/ESM); mutation & HLA frequencies are population-dependent literature estimates "
                   "joining independent datasets (not a same-patient cohort); class-I/CD8 only.",
        "INTERPRETATION_MAP": {
            "A_positive": "neoantigen SAFE addressability > HLA-LOH ceiling -> a tumour-exclusive handle that "
                          "generalises; the positive we hunted.",
            "B_safer": "comparable addressability but cleaner discrimination (clean/anchor tiers dominate) -> "
                       "same reach, fewer cross-reactivity fatalities.",
            "C_tcr_bound": "most public neoantigens fall in tcr_dependent -> decisive map of WHERE TCR "
                           "engineering (RUNG-12) must focus; the WT pMHC is on healthy cells.",
        },
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung11] wrote {RESULT_JSON}")

    # console summary
    print(f"\n  oracle (known clinical handles recovered): {n_recovered}/{len(KNOWN_POSITIVES)}")
    for k, v in validation.items():
        print(f"    {k:28} -> {'PRESENTED' if v and v['presented'] else 'MISSED  '} "
              f"({'rank %.2f, %s' % (v['mut_rank'], v['tier']) if v else 'n/a'})")
    print(f"\n  tier census (presented handles): {tier_counts}")
    print(f"  clean handles: {n_clean} ({n_clean_fragile} fragile WT in [2,4)); GOLD (strong mut + clear WT): {n_gold}")
    print(f"\n  cancer     SAFE(clean) strict..lenient   BROAD(all)   [HLA-LOH ceiling {HLA_LOH_CEILING[0]:.0%}-{HLA_LOH_CEILING[1]:.0%}]")
    for c in CANCERS:
        ss, sl, b = cov_safe_strict[c], cov_safe[c], cov_broad[c]
        print(f"  {c:9} {ss['central']:6.1%} .. {sl['central']:6.1%}      "
              f"{b['central']:6.1%}    n_mut={b['n_mutations_targetable']}")

    _make_figure(cov_safe_strict, cov_safe, cov_broad, tier_counts, validation)
    return 0


def _make_figure(cov_safe_strict, cov_safe, cov_broad, tier_counts, validation) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung11] matplotlib unavailable ({e}); skipped figure")
        return
    x = np.arange(len(CANCERS))
    safe_lo = [cov_safe_strict[c]["central"] * 100 for c in CANCERS]
    safe_hi = [cov_safe[c]["central"] * 100 for c in CANCERS]
    broad = [cov_broad[c]["central"] * 100 for c in CANCERS]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))

    w = 0.38
    ax[0].bar(x - w / 2, broad, w, label="BROAD (all presented; needs TCR eng.)", color="#3B7DD8")
    # SAFE shown as a range: solid bar to strict (robust clean), error-cap up to lenient clean
    ax[0].bar(x + w / 2, safe_lo, w, label="SAFE strict (WT clearly off MHC)", color="#1B5E20")
    ax[0].errorbar(x + w / 2, safe_lo, yerr=[[0] * len(CANCERS), np.array(safe_hi) - np.array(safe_lo)],
                   fmt="none", ecolor="#66BB6A", elinewidth=6, alpha=0.6, capsize=3,
                   label="SAFE lenient (WT just over cutoff)")
    ax[0].axhspan(HLA_LOH_CEILING[0] * 100, HLA_LOH_CEILING[1] * 100, color="#C1432B", alpha=0.15)
    ax[0].axhline(HLA_LOH_CEILING[1] * 100, color="#C1432B", ls="--", lw=1,
                  label=f"HLA-LOH ceiling {HLA_LOH_CEILING[0]:.0%}-{HLA_LOH_CEILING[1]:.0%} (RUNG-6)")
    for xi, v in zip(x - w / 2, broad):
        ax[0].text(xi, v + 0.6, f"{v:.0f}", ha="center", fontsize=7)
    ax[0].set_xticks(x); ax[0].set_xticklabels(CANCERS, rotation=30, ha="right")
    ax[0].set_ylabel("% patients addressable")
    ax[0].set_title("Neoantigen (sequence axis) addressability\nvs the expression-axis HLA-LOH ceiling")
    ax[0].legend(fontsize=6.5, loc="upper right"); ax[0].grid(axis="y", alpha=0.3)

    # discriminability tier census + oracle
    tiers = ["clean", "anchor", "tcr_dependent"]
    colors = {"clean": "#1B5E20", "anchor": "#F9A825", "tcr_dependent": "#C1432B"}
    counts = [tier_counts.get(t, 0) for t in tiers]
    ax[1].bar(tiers, counts, color=[colors[t] for t in tiers])
    for i, v in enumerate(counts):
        ax[1].text(i, v + 0.3, str(v), ha="center", fontsize=9)
    n_rec = sum(1 for v in validation.values() if v and v["presented"])
    ax[1].set_ylabel("# presented (mutation x allele) handles")
    ax[1].set_title(f"Discriminability tiers (safety)\noracle: {n_rec}/{len(validation)} known clinical handles recovered")
    ax[1].grid(axis="y", alpha=0.3)

    fig.suptitle("RUNG-11: public-neoantigen addressability x discriminability "
                 "(tumour-EXCLUSIVE recognition axis)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung11] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    """Exercise sequence/register/tier/coverage logic on a mock sequence + mock predictor (no deps, offline)."""
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # --- register generation on a known mini-protein. pos=5 (1-based), wt='G'->mut='D'
    seq = "AAAAGAAAA"          # residue 5 (index 4) is G
    regs = gen_registers(seq, 5, "G", "D", lens=(8, 9))
    # all registers must (a) contain the mutation, (b) differ from WT only at one position
    diffs_ok = all(sum(c1 != c2 for c1, c2 in zip(r["pep_wt"], r["pep_mut"])) == 1 for r in regs)
    check("registers differ from WT at exactly one position", diffs_ok)
    contains = all("D" in r["pep_mut"] and r["pep_wt"][r["j"]] == "G" for r in regs)
    check("mutated residue present at index j in every register", contains)
    # length-9 register starting at 0 places G at index 4 -> P5 (TCR-facing, not anchor)
    r9 = [r for r in regs if r["len"] == 9][0]
    check("9-mer full-span: mutation at P5 is non-anchor", r9["p_in_pep"] == 5 and not r9["anchor"])
    # wrong WT residue must raise
    try:
        gen_registers(seq, 5, "Q", "D"); raised = False
    except AssertionError:
        raised = True
    check("gen_registers asserts on wrong WT residue", raised)

    # --- anchor detection: a long-enough sequence yields a register with the mutation at P2 (anchor).
    # (a short mock can't: the register start is capped by len-L, so the mutation never reaches P2.)
    seqa = "A" * 7 + "G" + "A" * 7    # length 15, index 7 ('G'), pos 8 -> start=6 gives j=1 (P2)
    ra = gen_registers(seqa, 8, "G", "D", lens=(9,))
    has_anchor = any(r["anchor"] and r["p_in_pep"] == 2 for r in ra)
    check("anchor flagged when mutation lands at P2", has_anchor)

    # --- tier classification truth table
    check("tier clean  (mut yes, WT no)",  classify_tier(True, False, False) == "clean")
    check("tier anchor (both, anchor)",    classify_tier(True, True, True) == "anchor")
    check("tier tcr    (both, non-anchor)",classify_tier(True, True, False) == "tcr_dependent")
    check("tier none   (mut no)",          classify_tier(False, False, False) == "not_presented")

    # --- evaluate_driver with a MOCK scorer: build registers, inject ranks, check best-pick + tier
    regs2 = gen_registers("AAAAGAAAA", 5, "G", "D", lens=(9,))
    rg = regs2[0]
    mock = {(rg["pep_mut"], "HLA-A*02:01"): {"rank": 0.3, "pres": 0.9},   # mut strong binder
            (rg["pep_wt"], "HLA-A*02:01"): {"rank": 50.0, "pres": 0.1}}   # WT non-binder -> clean
    ev = evaluate_driver("G5D", regs2, ["HLA-A*02:01"], mock)
    check("evaluate_driver picks presented mut + tier clean when WT not presented",
          ev["HLA-A*02:01"]["tier"] == "clean" and ev["HLA-A*02:01"]["mut_presented"])

    # --- coverage math. one mutation p_mut=0.4 presented by one allele carrier-freq 0.5 -> 0.20
    handles = [{"gene": "X", "pos": 1, "mut_label": "A1B", "allele": "HLA-A*02:01",
                "tier": "clean", "cancer_prev": {"PDAC": 0.4}, "anchor": False, "mut_rank": 0.3}]
    cov = coverage(handles, {"HLA-A*02:01": 0.5}, SAFE_TIERS)
    check("coverage single handle = p_mut * carrier_freq", abs(cov["PDAC"]["central"] - 0.20) < 1e-9)
    # adding a second allele for the SAME mutation only raises coverage (more carriers presented)
    handles2 = handles + [{"gene": "X", "pos": 1, "mut_label": "A1B", "allele": "HLA-A*11:01",
                           "tier": "clean", "cancer_prev": {"PDAC": 0.4}, "anchor": False, "mut_rank": 0.3}]
    cov2 = coverage(handles2, {"HLA-A*02:01": 0.5, "HLA-A*11:01": 0.5}, SAFE_TIERS)
    check("second presenting allele raises coverage", cov2["PDAC"]["central"] > cov["PDAC"]["central"])
    # mutually-exclusive codon variants are SUMMED (not double-discounted): two muts at same codon,
    # p 0.3 & 0.2, each presented w.p.1 (carrier 1.0) -> central = 0.3+0.2 = 0.5 (product-of-one-group)
    hx = [{"gene": "K", "pos": 12, "mut_label": "G12D", "allele": "A", "tier": "clean",
           "cancer_prev": {"PDAC": 0.3}, "anchor": False, "mut_rank": 0.1},
          {"gene": "K", "pos": 12, "mut_label": "G12V", "allele": "A", "tier": "clean",
           "cancer_prev": {"PDAC": 0.2}, "anchor": False, "mut_rank": 0.1}]
    covx = coverage(hx, {"A": 1.0}, SAFE_TIERS)
    check("same-codon variants summed (mutually exclusive)", abs(covx["PDAC"]["central"] - 0.5) < 1e-9)
    # SAFE (clean only) <= BROAD (all presented), always
    hb = hx + [{"gene": "K", "pos": 12, "mut_label": "G12C", "allele": "A", "tier": "tcr_dependent",
                "cancer_prev": {"PDAC": 0.1}, "anchor": False, "mut_rank": 0.4}]
    safe_c = coverage(hb, {"A": 1.0}, SAFE_TIERS)["PDAC"]["central"]
    broad_c = coverage(hb, {"A": 1.0}, BROAD_TIERS)["PDAC"]["central"]
    check("SAFE(clean) <= BROAD(all presented)", safe_c <= broad_c + 1e-12)

    # --- Jeffreys sanity
    check("jeffreys_upper(0,100) > 0 (no false zero)", jeffreys_upper(0, 100) > 0.0)

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="RUNG-11 public-neoantigen addressability x discriminability")
    ap.add_argument("mode", nargs="?", default="run", choices=["run", "selftest"])
    args = ap.parse_args()
    sys.exit(selftest() if args.mode == "selftest" else main_run())
