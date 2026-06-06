#!/usr/bin/env python3
"""
RUNG 12 — per-handle pMHC structural discriminability: replace the SWEPT beta with a MEASURED one.

WHAT THIS CLOSES
----------------
The RUNG-12P bridge (scripts/36) showed a gated relay unlocks RUNG-11's `tcr_dependent` neoantigen handles
IF the TCR's mut-vs-WT cross-bind `beta` is low enough -- but `beta` was SWEPT (unknown). This rung MEASURES a
per-handle `beta` proxy from structure + sequence, turning the what-if into a RANKED, CERTIFIED target list.

THE SIGNALS (per handle: mutant pMHC vs wild-type pMHC)
------------------------------------------------------
Discriminability D in [0,1] = how distinguishable the mutant pMHC is from WT to a TCR. beta (cross-bind) = 1-D.
D combines independent discrimination MECHANISMS (probabilistic OR -> D high if ANY discriminates):
  - M  MHC-level: WT binds much worse / not presented (from RUNG-11 NetMHCpan mut_rank vs wt_rank). clean->1.
  - E*P  TCR-level: the mutated residue is solvent-EXPOSED (E = RSA from the ESMFold-bound conformation) AND
         physicochemically DIFFERENT (P = charge/volume/hydropathy delta). A buried or conservative change is
         invisible to the TCR.
  - Z  sequence-context: ESM-2 embedding distance between mut and WT peptide.
  D = 1 - (1-M)(1-E*P)(1-Z).  Structural change (peptide RMSD mut vs WT) + pLDDT are reported for transparency.

Then per handle: q_n = presentation_factor(wt_rank) * beta  ->  per-cell-safe (<=0.02) / relay-safe (<=0.17) ->
re-run the bridge coverage with MEASURED beta (no longer swept).

STRUCTURE ENGINE: ESMFold (single-sequence, NO MSA server -> robust on Colab, fast on a T4). pMHC = HLA alpha1-
alpha2 groove (mature 1-182, fetched from IPD-IMGT/HLA) + ':' + peptide (ESMFold multimer chainbreak).

HONEST CEILING
--------------
ESMFold is a single-sequence MODEL; peptide docking into the groove can be imperfect -> E (exposure) is the
softest signal (we report pLDDT and lean on M+P+Z which don't depend on fold quality). beta is a PROXY, not a
measured TCR Kd; a surviving handle is a prioritised HYPOTHESIS for wet-lab TCR isolation, not a validated
target. Inherits RUNG-11 (frequencies) + RUNG-12P (relay ceiling) caveats. mRNA->presentation already flagged.

USAGE
  python scripts/37_pmhc_discriminability.py selftest
  python scripts/37_pmhc_discriminability.py prep      # select handles, fetch HLA grooves, write FASTAs+manifest
  python scripts/37_pmhc_discriminability.py analyze    # consume PDBs + ESM deltas -> measured beta + bridge
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REFS = PROJECT_ROOT / "data" / "refs"
RUNG11_JSON = PROJECT_ROOT / "runs" / "rung11_neoantigen" / "rung11_neoantigen_addressability.json"
RUNG12B_JSON = PROJECT_ROOT / "runs" / "rung12pB_relay" / "rung12pB_relay.json"
OUT_DIR = PROJECT_ROOT / "runs" / "rung12_pmhc"
WORK = Path(os.environ.get("RUNG12_WORK", str(OUT_DIR / "work")))   # FASTAs, PDBs, embeddings (Drive on Colab)
MANIFEST = OUT_DIR / "rung12_manifest.json"
RESULT_JSON = OUT_DIR / "rung12_pmhc.json"
FIGURE_PNG = OUT_DIR / "rung12_pmhc.png"

N_TOP = int(os.environ.get("RUNG12_N", "32"))     # handles to fold (x2 for mut+WT)
GROOVE_LEN = 182                                  # alpha1-alpha2 (peptide-binding groove) length
PER_CELL_BAR = 0.02
RELAY_CEILING_DEFAULT = 0.173                     # 3D from RUNG-12P/B (overridden from its JSON if present)
BINDER_RANK, STRONG_RANK = 2.0, 0.5
MATURE_MOTIF = re.compile("SHS[MLKQ][RK]YF")      # conserved HLA-I alpha1 N-terminus -> mature start

# physicochemical scales
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5}
VOLUME = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4, "G": 60.1,
          "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0,
          "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2,
      "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
# Tien et al. 2013 theoretical max ASA (for RSA normalisation)
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223, "G": 104, "H": 224,
          "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
          "Y": 263, "V": 174}


def _load(name, mod):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / mod)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


# ---------------------------------------------------------------------------
#  Pure scoring (selftest exercises these directly).
# ---------------------------------------------------------------------------
def presentation_factor(wt_rank: float) -> float:
    return float(np.clip((BINDER_RANK - wt_rank) / (BINDER_RANK - STRONG_RANK), 0.0, 1.0))


def physchem_delta(wt_aa: str, mut_aa: str) -> dict:
    dq = abs(CHARGE.get(wt_aa, 0) - CHARGE.get(mut_aa, 0)) / 2.0
    dv = abs(VOLUME[wt_aa] - VOLUME[mut_aa]) / (227.8 - 60.1)
    dh = abs(KD[wt_aa] - KD[mut_aa]) / 9.0
    P = float(np.clip(max(dq, dv, dh), 0.0, 1.0))      # max -> any large change makes it discriminable
    return {"d_charge": round(dq, 3), "d_volume": round(dv, 3), "d_hydropathy": round(dh, 3), "P": round(P, 3)}


def mhc_level_discrimination(tier: str, mut_rank: float, wt_rank: float) -> float:
    """WT presents much worse than mutant -> the WT pMHC is rarely on the surface (MHC-level discrimination)."""
    if tier == "clean":
        return 1.0
    return float(np.clip((wt_rank - mut_rank) / 2.0, 0.0, 1.0))


Z_WEIGHT = 0.4   # ESM-2 embedding distance is a SOFT, relative sequence hint (correlated with P) — it must
                 # NOT by itself declare a handle perfectly discriminable. Cap its OR-contribution.


def combine_discriminability(M: float, E: float, P: float, Z: float) -> dict:
    """Probabilistic-OR of independent discrimination mechanisms. Returns D and beta=1-D.
    M (MHC-binding differential) and E*P (exposure x physicochemistry) are the load-bearing signals; Z (ESM-2)
    is capped (Z_WEIGHT) so a relative-max embedding distance can't single-handedly force beta=0."""
    tcr = E * P
    z_eff = Z_WEIGHT * float(np.clip(Z, 0.0, 1.0))
    D = 1.0 - (1.0 - M) * (1.0 - tcr) * (1.0 - z_eff)
    D = float(np.clip(D, 0.0, 1.0))
    return {"M": round(M, 3), "E": round(E, 3), "P": round(P, 3), "tcr_EP": round(tcr, 3),
            "Z": round(Z, 3), "Z_eff": round(z_eff, 3), "D": round(D, 3), "beta": round(1.0 - D, 3)}


def qn_measured(wt_rank: float, beta: float) -> float:
    return presentation_factor(wt_rank) * beta


# ---------------------------------------------------------------------------
#  HLA groove fetch (alpha1-alpha2, mature 1-182).
# ---------------------------------------------------------------------------
def _read_fasta(text):
    recs, name, seq = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                recs[name] = "".join(seq)
            parts = line.split()
            name = parts[1] if len(parts) > 1 else line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if name:
        recs[name] = "".join(seq)
    return recs


def fetch_locus(locus: str) -> dict:
    cache = REFS / f"{locus}_prot.fasta"
    if not cache.exists():
        url = f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/{locus}_prot.fasta"
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()       # noqa: S310
        except Exception:
            import subprocess
            txt = subprocess.check_output(["curl", "-fsSL", url], text=True, timeout=120)
        REFS.mkdir(parents=True, exist_ok=True); cache.write_text(txt)
    return _read_fasta(cache.read_text())


def hla_groove(allele: str) -> str:
    """allele like 'HLA-A*02:01' -> mature alpha1-alpha2 (182 aa)."""
    a = allele.replace("HLA-", "")
    locus = a[0]
    recs = fetch_locus(locus)
    keys = sorted([k for k in recs if k.startswith(a)])
    if not keys:
        raise ValueError(f"no IMGT sequence for {allele}")
    s = recs[keys[0]].replace("*", "").replace("X", "")
    m = MATURE_MOTIF.search(s)
    if not m:
        raise ValueError(f"mature start motif not found for {allele}")
    groove = s[m.start():m.start() + GROOVE_LEN]
    if len(groove) < GROOVE_LEN - 5:
        raise ValueError(f"groove too short for {allele}: {len(groove)}")
    return groove


# ---------------------------------------------------------------------------
#  Structure analysis (Biopython; degrade gracefully if unavailable).
# ---------------------------------------------------------------------------
def analyze_pdb(pdb_path: Path, pep_index_1based: int):
    """ESMFold output: chain A = MHC groove, chain B = peptide. Return RSA of the mutated peptide residue,
    mean peptide pLDDT (B-factor). None if Biopython missing or parse fails."""
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except Exception:
        return None
    try:
        st = PDBParser(QUIET=True).get_structure("p", str(pdb_path))
        model = next(iter(st))
        chains = list(model)
        pep_chain = min(chains, key=lambda c: sum(1 for _ in c.get_residues()))   # peptide = shortest chain
        ShrakeRupley().compute(model, level="R")
        residues = [r for r in pep_chain.get_residues() if r.id[0] == " "]
        if not (1 <= pep_index_1based <= len(residues)):
            return None
        res = residues[pep_index_1based - 1]
        aa3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
                  "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
                  "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
        aa = aa3to1.get(res.resname, None)
        rsa = float(res.sasa) / MAXASA.get(aa, 200.0) if aa else None
        plddt = float(np.mean([a.bfactor for a in pep_chain.get_atoms()]))
        return {"rsa": round(min(rsa, 1.0), 3) if rsa is not None else None,
                "pep_pLDDT": round(plddt, 1), "pep_len": len(residues)}
    except Exception:
        return None


def peptide_rmsd(pdb_mut: Path, pdb_wt: Path):
    """Superpose on MHC CA, return peptide-CA RMSD (mut vs WT). None on failure."""
    try:
        from Bio.PDB import PDBParser, Superimposer
    except Exception:
        return None
    try:
        p = PDBParser(QUIET=True)
        def chains(path):
            m = next(iter(p.get_structure("x", str(path))))
            cs = sorted(m, key=lambda c: -sum(1 for _ in c.get_residues()))
            mhc = [r["CA"] for r in cs[0].get_residues() if "CA" in r]
            pep = [r["CA"] for r in cs[-1].get_residues() if "CA" in r]
            return mhc, pep
        mhc_m, pep_m = chains(pdb_mut); mhc_w, pep_w = chains(pdb_wt)
        n = min(len(mhc_m), len(mhc_w))
        si = Superimposer(); si.set_atoms(mhc_m[:n], mhc_w[:n])
        rot, tran = si.rotran
        npep = min(len(pep_m), len(pep_w))
        cm = np.array([a.coord for a in pep_m[:npep]])
        cw = np.array([a.coord for a in pep_w[:npep]]) @ rot + tran
        return round(float(np.sqrt(((cm - cw) ** 2).sum(1).mean())), 3) if npep else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
def select_handles(d11, n_top):
    """Top presented handles by max cancer prevalence (naturally surfaces IDH1-R132H glioma, KRAS PDAC, etc.)."""
    posmap = None
    handles = []
    for _drv, info in d11["per_driver"].items():
        gene, mut_label, prev = info["gene"], info["mutation"], info["cancer_prev"]
        for allele, b in info["by_allele"].items():
            j = b["p_in_pep"] - 1
            handles.append({"id": f"{gene}_{mut_label}_{allele.replace('HLA-','').replace('*','').replace(':','')}",
                            "gene": gene, "mut_label": mut_label, "allele": allele,
                            "pep_mut": b["pep_mut"], "pep_wt": b["pep_wt"], "p_in_pep": b["p_in_pep"],
                            "wt_aa": b["pep_wt"][j], "mut_aa": b["pep_mut"][j],
                            "wt_rank": b["wt_rank"], "mut_rank": b["mut_rank"], "tier": b["tier"],
                            "anchor": b["anchor"], "cancer_prev": prev, "max_prev": max(prev.values())})
    handles.sort(key=lambda h: -h["max_prev"])
    return handles[:n_top]


def main_prep() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True); WORK.mkdir(parents=True, exist_ok=True)
    if not RUNG11_JSON.exists():
        print("[rung12] need RUNG-11 result first"); return 2
    d11 = json.loads(RUNG11_JSON.read_text())
    handles = select_handles(d11, N_TOP)
    grooves, ok = {}, 0
    for h in handles:
        try:
            if h["allele"] not in grooves:
                grooves[h["allele"]] = hla_groove(h["allele"])
            h["groove"] = grooves[h["allele"]]
            # write ESMFold inputs: groove:peptide (mut and WT)
            (WORK / f"{h['id']}_mut.fasta").write_text(f">{h['id']}_mut\n{h['groove']}:{h['pep_mut']}\n")
            (WORK / f"{h['id']}_wt.fasta").write_text(f">{h['id']}_wt\n{h['groove']}:{h['pep_wt']}\n")
            ok += 1
        except Exception as e:
            h["groove_error"] = str(e)
            print(f"[rung12] {h['id']}: groove fetch FAILED ({e})")
    MANIFEST.write_text(json.dumps({"handles": handles, "groove_len": GROOVE_LEN, "n": len(handles),
                                    "n_grooves_ok": ok}, indent=2))
    print(f"[rung12] prepped {ok}/{len(handles)} handles -> {WORK} (mut+WT FASTAs) | manifest {MANIFEST}")
    print(f"[rung12] next (Colab GPU): ESMFold each *.fasta -> <id>_{{mut,wt}}.pdb; ESM-2 embed peptides; then `analyze`")
    return 0


def main_analyze() -> int:
    d33 = _load("d33", "33_neoantigen_addressability.py")
    if not MANIFEST.exists():
        print("[rung12] run `prep` (and the Colab fold/embed) first"); return 2
    man = json.loads(MANIFEST.read_text())
    relay_ceiling = RELAY_CEILING_DEFAULT
    if RUNG12B_JSON.exists():
        relay_ceiling = json.loads(RUNG12B_JSON.read_text())["sensitivity_3D"]["safe_q_n_ceiling_at_1pct"]
    # optional ESM-2 embedding deltas written by the notebook: {id: cosine_distance}
    emb = {}
    ef = WORK / "esm_deltas.json"
    if ef.exists():
        emb = json.loads(ef.read_text())

    rows = []
    for h in man["handles"]:
        if "groove" not in h:
            continue
        j = h["p_in_pep"]
        struct_m = analyze_pdb(WORK / f"{h['id']}_mut.pdb", j)
        rmsd = peptide_rmsd(WORK / f"{h['id']}_mut.pdb", WORK / f"{h['id']}_wt.pdb")
        # exposure E: prefer measured RSA; else fall back to position-based prior (anchor buried, tcr-facing up)
        if struct_m and struct_m.get("rsa") is not None:
            E = struct_m["rsa"]; e_src = "ESMFold RSA"
        else:
            E = 0.2 if h["anchor"] else 0.7; e_src = "position prior (no structure)"
        P = physchem_delta(h["wt_aa"], h["mut_aa"])["P"]
        M = mhc_level_discrimination(h["tier"], h["mut_rank"], h["wt_rank"])
        Z = float(np.clip(emb.get(h["id"], 0.0), 0.0, 1.0))
        disc = combine_discriminability(M, E, P, Z)
        beta = disc["beta"]
        qn = qn_measured(h["wt_rank"], beta)
        rows.append({**{k: h[k] for k in ("id", "gene", "mut_label", "allele", "tier", "anchor",
                                          "wt_rank", "mut_rank", "cancer_prev")},
                     "wt_aa": h["wt_aa"], "mut_aa": h["mut_aa"], "physchem_P": P, "E": E, "E_source": e_src,
                     "pep_pLDDT": (struct_m or {}).get("pep_pLDDT"), "pep_rmsd_mut_wt": rmsd,
                     "esm_Z": Z, **disc, "qn_measured": round(qn, 4),
                     "per_cell_safe": qn <= PER_CELL_BAR, "relay_safe": qn <= relay_ceiling})

    # bridge with MEASURED per-handle beta: build coverage-handles and classify by measured safety
    cov_handles = []
    posmap = {(g, f"{wt}{pos}{mut}"): pos for (g, _a, pos, wt, mut, _p) in d33.DRIVERS}
    for r in rows:
        cov_handles.append({"gene": r["gene"], "pos": posmap.get((r["gene"], r["mut_label"])),
                            "mut_label": r["mut_label"], "allele": r["allele"], "tier": r["tier"],
                            "cancer_prev": r["cancer_prev"], "per_cell_safe": r["per_cell_safe"],
                            "relay_safe": r["relay_safe"]})
    cov_pc = d33.coverage(cov_handles, d33.HLA_PANEL, lambda h: h["per_cell_safe"])
    cov_rl = d33.coverage(cov_handles, d33.HLA_PANEL, lambda h: h["relay_safe"])

    ranked = sorted(rows, key=lambda r: (-int(r["relay_safe"]), -r["D"], -max(r["cancer_prev"].values())))
    n_struct = sum(1 for r in rows if r["E_source"].startswith("ESMFold"))
    result = {
        "tag": "rung12_pmhc_discriminability",
        "question": "Per top neoantigen handle, how structurally discriminable is mutant pMHC from wild-type "
                    "(measured beta), and how does the relay-usable target set look with MEASURED (not swept) beta?",
        "n_handles": len(rows), "n_with_structure": n_struct, "relay_ceiling_used": relay_ceiling,
        "per_cell_bar": PER_CELL_BAR,
        "n_per_cell_safe": sum(1 for r in rows if r["per_cell_safe"]),
        "n_relay_safe": sum(1 for r in rows if r["relay_safe"]),
        "n_unlocked_by_relay": sum(1 for r in rows if r["relay_safe"] and not r["per_cell_safe"]),
        "coverage_measured": {c: {"per_cell": round(cov_pc[c]["central"], 4),
                                  "relay": round(cov_rl[c]["central"], 4)} for c in d33.CANCERS},
        "ranked_targets": ranked,
        "CEILING": "ESMFold single-sequence MODEL (peptide docking can be imperfect -> E/RSA is the softest "
                   "signal; pLDDT reported, weight on M+P+Z which are fold-independent). beta is a PROXY not a "
                   "TCR Kd; a top handle is a prioritised HYPOTHESIS for wet-lab TCR isolation. Inherits RUNG-11 "
                   "frequency + RUNG-12P relay-ceiling caveats.",
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung12] wrote {RESULT_JSON}  ({n_struct}/{len(rows)} with ESMFold structure)")
    print(f"  per-cell-safe {result['n_per_cell_safe']} | relay-safe {result['n_relay_safe']} | "
          f"unlocked-by-relay {result['n_unlocked_by_relay']}")
    print("\n  top targets (relay-safe first, by discriminability D):")
    for r in ranked[:12]:
        tag = "PERCELL" if r["per_cell_safe"] else ("RELAY" if r["relay_safe"] else "risky")
        print(f"    {r['id']:34} {r['tier']:13} D={r['D']:.2f} beta={r['beta']:.2f} qn={r['qn_measured']:.3f} "
              f"E={r['E']:.2f} P={r['physchem_P']:.2f} -> {tag}")
    _make_figure(result, d33.CANCERS)
    return 0


def _make_figure(result, cancers):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung12] matplotlib unavailable ({e})"); return
    rows = result["ranked_targets"][:20][::-1]
    names = [r["id"].replace("_", " ") for r in rows]
    betas = [r["beta"] for r in rows]
    colors = ["#1B5E20" if r["per_cell_safe"] else ("#F9A825" if r["relay_safe"] else "#C1432B") for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(14, max(4, 0.42 * len(rows) + 1.5)))
    y = np.arange(len(rows))
    ax[0].barh(y, [1 - b for b in betas], color=colors)
    ax[0].set_yticks(y); ax[0].set_yticklabels(names, fontsize=6.5)
    ax[0].set_xlabel("discriminability D = 1 - beta (green=per-cell-safe, amber=relay-only, red=risky)")
    ax[0].set_title("Measured per-handle pMHC discriminability\n(mutant vs wild-type)")
    ax[0].axvline(1 - result["relay_ceiling_used"], ls=":", color="grey")
    ax[0].grid(axis="x", alpha=0.3)

    cm = result["coverage_measured"]
    x = np.arange(len(cancers)); w = 0.38
    ax[1].bar(x - w / 2, [cm[c]["per_cell"] * 100 for c in cancers], w, label="per-cell usable", color="#C1432B")
    ax[1].bar(x + w / 2, [cm[c]["relay"] * 100 for c in cancers], w, label="relay usable", color="#1B5E20")
    ax[1].set_xticks(x); ax[1].set_xticklabels(cancers, rotation=30, ha="right")
    ax[1].set_ylabel("% patients usable (MEASURED beta)")
    ax[1].set_title("Usable addressability with MEASURED beta\n(top-handle subset)")
    ax[1].legend(fontsize=8); ax[1].grid(axis="y", alpha=0.3)
    fig.suptitle("RUNG-12: pMHC structural discriminability -> ranked, certified relay targets", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung12] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest() -> int:
    checks, ok = [], 0

    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # physicochemical delta: charge-flipping G->D big; conservative I->L tiny
    check("G->D (adds charge+volume) P large", physchem_delta("G", "D")["P"] > 0.4)
    check("I->L (conservative) P small", physchem_delta("I", "L")["P"] < 0.2)
    check("charge flip detected E->K", physchem_delta("E", "K")["d_charge"] == 1.0)

    # MHC-level discrimination
    check("clean tier -> M=1", mhc_level_discrimination("clean", 0.5, 3.0) == 1.0)
    check("WT much weaker -> M high", mhc_level_discrimination("tcr_dependent", 0.1, 1.9) > 0.5)
    check("WT ~ mut -> M~0", mhc_level_discrimination("tcr_dependent", 0.6, 0.65) < 0.1)

    # discriminability combination: any mechanism firing -> D high; none -> D low (beta high)
    dlo = combine_discriminability(0.0, 0.7, 0.05, 0.0)   # exposed but tiny chem change, WT presents -> low D
    dhi = combine_discriminability(0.0, 0.9, 0.9, 0.0)    # exposed + big chem change -> high D
    check("exposed+big-change -> high D (low beta)", dhi["beta"] < 0.2)
    check("exposed+tiny-change -> low D (high beta)", dlo["beta"] > 0.6)
    check("MHC discrimination alone rescues D", combine_discriminability(1.0, 0.0, 0.0, 0.0)["D"] == 1.0)
    check("ESM-2 Z ALONE cannot force beta=0 (capped)", combine_discriminability(0.0, 0.0, 0.0, 1.0)["beta"] > 0.4)

    # presentation_factor + qn
    check("clean (wt_rank>2) -> pf=0 -> qn=0", qn_measured(2.5, 0.9) == 0.0)
    check("qn = pf*beta", abs(qn_measured(0.5, 0.4) - 0.4) < 1e-9)

    # HLA groove extraction on a synthetic IMGT-style record
    syn = ">HLA:HLA9 A*99:99:01:01 365 bp\n" + "M" * 24 + "SHSMRYF" + "G" * 200 + "\n"
    (REFS).mkdir(parents=True, exist_ok=True)
    (REFS / "Z_prot.fasta").write_text(syn)
    rec = _read_fasta(syn)
    check("fasta parse keys allele name", "A*99:99:01:01" in rec)
    s = rec["A*99:99:01:01"]; m = MATURE_MOTIF.search(s)
    check("mature motif at idx 24", m.start() == 24)
    groove = s[m.start():m.start() + GROOVE_LEN]
    check("groove length 182", len(groove) == GROOVE_LEN and groove.startswith("SHSMRYF"))
    (REFS / "Z_prot.fasta").unlink()

    # handle selection from a synthetic RUNG-11-like JSON (IDH1 R132H glioma should top by prevalence)
    d11 = {"per_driver": {
        "IDH1 R132H": {"gene": "IDH1", "mutation": "R132H", "cancer_prev": {"GLIOMA": 0.70},
                       "by_allele": {"HLA-A*03:01": {"tier": "anchor", "p_in_pep": 4, "anchor": True,
                                                     "pep_mut": "AAAHAAAAA", "pep_wt": "AAARAAAAA",
                                                     "wt_rank": 1.0, "mut_rank": 0.5, "len": 9, "strong": True,
                                                     "mut_presented": True, "wt_presented": True}}},
        "KRAS G12D": {"gene": "KRAS", "mutation": "G12D", "cancer_prev": {"PDAC": 0.40},
                      "by_allele": {"HLA-C*08:02": {"tier": "tcr_dependent", "p_in_pep": 3, "anchor": False,
                                                    "pep_mut": "GADGVGKSAL", "pep_wt": "GAGGVGKSAL",
                                                    "wt_rank": 1.05, "mut_rank": 0.05, "len": 10, "strong": True,
                                                    "mut_presented": True, "wt_presented": True}}}}}
    hs = select_handles(d11, 5)
    check("handle selection orders by max prevalence (IDH1 glioma first)", hs[0]["gene"] == "IDH1")
    check("handle carries mutated residue indices", hs[0]["wt_aa"] == "R" and hs[0]["mut_aa"] == "H")

    total = len(checks)
    print(f"\nselftest: {ok}/{total} checks passed")
    return 0 if ok == total else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="RUNG-12 pMHC structural discriminability")
    ap.add_argument("mode", nargs="?", default="analyze", choices=["prep", "analyze", "selftest"])
    args = ap.parse_args()
    sys.exit({"prep": main_prep, "analyze": main_analyze, "selftest": selftest}[args.mode]())
