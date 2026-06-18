#!/usr/bin/env python3
"""
Generalized pMHC staging: fold CIFs -> PXDesign upload target + MUT/WT scoring PDBs + meta.json.
Reusable across targets (BRAF used scripts/58; this is the parameterized version for PIK3CA/KRAS/...).

Usage:
  python3 scripts/63_stage_target.py --run runs/rung29_pik3ca_e545k --name pik3ca_e545k \
      --pep_mut STRDPLSEITK --pep_wt STRDPLSEITE --p 11 --hotspot none --allele HLA-A*03:01
  python3 scripts/63_stage_target.py --run runs/rung30_kras_g12d --name kras_g12d \
      --pep_mut VVVGADGVGK --pep_wt VVVGAGGVGK --p 6 --hotspot p6 --allele HLA-A*03:01

Reads <run>/folds/{MUT,WT}/<name>_*_sample_0.cif, writes <run>/staging/.
hotspot 'none' = design a strong binder to the up-facing surface (selectivity by presentation, e.g. PIK3CA).
hotspot 'pN'   = pin the design hotspot on peptide position N (read-the-mutation, e.g. KRAS p6).
"""
import os, sys, json, glob, argparse, warnings
import numpy as np
warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BACKBONE = {"N", "CA", "C", "O"}
GROOVE_PREFIX = "SHSMRYF"


def load(cif, pep_len):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    m = next(iter(MMCIFParser(QUIET=True).get_structure("x", cif)))
    groove = peptide = None
    for ch in m:
        res = [r for r in ch if r.id[0] == " "]
        seq = "".join(seq1(r.resname) if len(r.resname) == 3 else "X" for r in res)
        if seq.startswith(GROOVE_PREFIX):
            groove = (ch.id, res, seq)
        elif len(res) == pep_len:
            peptide = (ch.id, res, seq)
    assert groove and peptide, f"chain id failed in {cif}: groove={bool(groove)} pep={bool(peptide)}"
    return groove, peptide


def write_pdb(path, chain_specs):
    serial, lines = 1, []
    for cid, residues in chain_specs:
        for i, r in enumerate(residues):
            for a in r:
                x, y, z = a.get_coord(); nm = a.get_name()
                name = f" {nm:<3s}" if len(nm) < 4 else nm
                el = (a.element or nm[0]).strip()[:2]
                lines.append(f"ATOM  {serial:5d} {name:<4s}{'':1s}{r.resname:>3s} {cid:1s}{i+1:4d}{'':1s}   "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
                serial += 1
        lines.append(f"TER   {serial:5d}      {residues[-1].resname:>3s} {cid:1s}{len(residues):4d}"); serial += 1
    lines.append("END"); open(path, "w").write("\n".join(lines) + "\n")


def crop_groove(groove_res, peptide_res, radius=10.0):
    pep = np.array([a.get_coord() for r in peptide_res for a in r])
    kept = []
    for r in groove_res:
        ga = np.array([a.get_coord() for a in r])
        if np.linalg.norm(ga[:, None, :] - pep[None, :, :], axis=-1).min() <= radius:
            kept.append(r)
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--pep_mut", required=True); ap.add_argument("--pep_wt", required=True)
    ap.add_argument("--p", type=int, required=True); ap.add_argument("--hotspot", default="none")
    ap.add_argument("--allele", default="HLA-A*03:01")
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run); OUT = os.path.join(run, "staging"); os.makedirs(OUT, exist_ok=True)
    mut_cif = sorted(glob.glob(f"{run}/folds/MUT/*_sample_0.cif"))[0]
    wt_cif = sorted(glob.glob(f"{run}/folds/WT/*_sample_0.cif"))[0]
    gM, pM = load(mut_cif, len(a.pep_mut)); gW, pW = load(wt_cif, len(a.pep_wt))
    assert pM[2] == a.pep_mut, f"MUT peptide {pM[2]} != {a.pep_mut}"
    assert pW[2] == a.pep_wt, f"WT peptide {pW[2]} != {a.pep_wt}"
    assert pM[2][a.p - 1] != pW[2][a.p - 1], f"p{a.p} identical in MUT/WT"

    write_pdb(f"{OUT}/{a.name}_mut_pmhc.pdb", [("A", gM[1]), ("B", pM[1])])
    write_pdb(f"{OUT}/{a.name}_wt_pmhc.pdb", [("A", gW[1]), ("B", pW[1])])
    crop = crop_groove(gM[1], pM[1], 10.0)
    write_pdb(f"{OUT}/{a.name}_mut_cropped.pdb", [("A", crop), ("B", pM[1])])

    REPO_STAGE = f"/content/cancer-recog-apoptosis/{a.run}/staging"
    meta = {"target": a.name, "allele": a.allele, "pep_mut": a.pep_mut, "pep_wt": a.pep_wt,
            "p_in_pep": a.p, "mut_aa": pM[2][a.p - 1], "wt_aa": pW[2][a.p - 1],
            "design_hotspot": a.hotspot,
            "design_note": ("NO hotspot — strong binder to up-facing surface; selectivity by presentation"
                            if a.hotspot == "none" else
                            f"hotspot on peptide {a.hotspot} (read-the-mutation)"),
            "mut_pdb": f"{REPO_STAGE}/{a.name}_mut_pmhc.pdb",
            "wt_pdb": f"{REPO_STAGE}/{a.name}_wt_pmhc.pdb",
            "cropped_design_target": f"{REPO_STAGE}/{a.name}_mut_cropped.pdb",
            "source": f"{a.run}/folds/{{MUT,WT}}/*_sample_0.cif"}
    json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=2)
    print(f"[staged] {OUT}")
    print(f"  groove {len(gM[1])} aa | peptide {pM[2]} (p{a.p}: {pM[2][a.p-1]} vs WT {pW[2][a.p-1]}) | "
          f"cropped groove kept {len(crop)} aa | hotspot={a.hotspot}")


if __name__ == "__main__":
    main()
