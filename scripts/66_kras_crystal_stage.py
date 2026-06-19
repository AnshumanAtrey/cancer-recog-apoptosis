#!/usr/bin/env python3
"""
RUNG-30b — stage the KRAS-G12D / HLA-A*11:01 design target from a REAL CRYSTAL STRUCTURE (PDB 7OW6, 2.64 A;
Sim et al PMID 36088370) instead of a fold. Higher quality + the better-validated allele (A*11:01: IEDB
MS-eluted + T-cell-positive, RUNG-32b). 7OW6 is a TCR-pMHC complex; we strip the TCR + b2m, keep the
A*11:01 groove (chain A) + the KRAS-G12D peptide VVVGADGVGK (chain C), crop to peptide + 10 A groove, and
SASA-confirm the hotspot (KRAS needs read-the-mutation -> hotspot ON the up-facing G12D Asp).
"""
import os, json, warnings
import numpy as np
warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CIF = os.path.join(ROOT, "data/pdb_kras/7OW6.cif")
OUT = os.path.join(ROOT, "runs/rung30_kras_g12d/staging_crystal")
PEP = "VVVGADGVGK"; P_MUT = 6   # G12D at peptide p6 (Asp)
BACKBONE = {"N", "CA", "C", "O"}


def load():
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    m = next(iter(MMCIFParser(QUIET=True).get_structure("x", CIF)))
    groove = peptide = None
    for ch in m:
        res = [r for r in ch if r.id[0] == " "]
        seq = "".join(seq1(r.resname) if len(r.resname) == 3 else "" for r in res)
        if "SHSMRYF" in seq[:10] and groove is None:
            groove = res
        elif seq == PEP and peptide is None:
            peptide = res
    assert groove and peptide, f"groove={bool(groove)} peptide={bool(peptide)}"
    return m, groove, peptide


def write_pdb(path, chain_specs):
    serial, lines = 1, []
    for cid, residues in chain_specs:
        for i, r in enumerate(residues):
            for a in r:
                if a.get_name() not in ("N", "CA", "C", "O", "CB") and len(a.get_name()) > 0:
                    pass
                x, y, z = a.get_coord(); nm = a.get_name()
                if a.element == "H":
                    continue
                name = f" {nm:<3s}" if len(nm) < 4 else nm
                el = (a.element or nm[0]).strip()[:2]
                lines.append(f"ATOM  {serial:5d} {name:<4s}{'':1s}{r.resname:>3s} {cid:1s}{i+1:4d}{'':1s}   "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
                serial += 1
        lines.append(f"TER   {serial:5d}      {residues[-1].resname:>3s} {cid:1s}{len(residues):4d}"); serial += 1
    lines.append("END"); open(path, "w").write("\n".join(lines) + "\n")


def crop_groove(groove, peptide, radius=10.0):
    pep = np.array([a.get_coord() for r in peptide for a in r])
    return [r for r in groove
            if np.linalg.norm(np.array([a.get_coord() for a in r])[:, None, :] - pep[None, :, :], axis=-1).min() <= radius]


def peptide_sasa(pmhc_pdb):
    """SASA on the ISOLATED pMHC (groove+peptide only, no TCR/b2m) -> true presented exposure.
    (Computing on the full crystal would have the TCR occlude the peptide = artifactually buried.)"""
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.SeqUtils import seq1
    model = next(iter(PDBParser(QUIET=True).get_structure("x", pmhc_pdb)))
    ShrakeRupley().compute(model, level="A")
    peptide = [r for r in model["B"] if r.id[0] == " "]
    maxs = {"VAL":160,"ALA":129,"GLY":104,"ASP":193,"LYS":236}
    out = []
    for i, r in enumerate(peptide):
        sc = [a for a in r if a.get_name() not in BACKBONE]
        sasa = sum(a.sasa for a in sc) if sc else 0.0
        rel = sasa / maxs.get(r.resname, 180) * 100
        out.append((i + 1, r.resname, seq1(r.resname), round(sasa, 1), round(rel),
                    "UP" if rel >= 25 else ("mid" if rel >= 12 else "buried")))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    m, groove, peptide = load()
    from Bio.SeqUtils import seq1
    assert "".join(seq1(r.resname) for r in peptide) == PEP
    assert peptide[P_MUT - 1].resname == "ASP", f"p{P_MUT} not ASP"
    pmhc_pdb = os.path.join(OUT, "kras_g12d_A1101_mut_pmhc.pdb")
    write_pdb(pmhc_pdb, [("A", groove), ("B", peptide)])
    crop = crop_groove(groove, peptide, 10.0)
    write_pdb(os.path.join(OUT, "kras_g12d_A1101_mut_cropped.pdb"), [("A", crop), ("B", peptide)])
    sasa = peptide_sasa(pmhc_pdb)   # isolated pMHC (no TCR) = true presented exposure
    print("KRAS-G12D peptide VVVGADGVGK on A*11:01 (crystal 7OW6, TCR-stripped) — sidechain SASA:")
    for pos, rn, aa, s, rel, cls in sasa:
        tag = " <-- p6 G12D (the read-the-mutation hotspot)" if pos == P_MUT else ""
        print(f"  p{pos:>2} {rn} {aa}: {s:5.1f} A^2 ({rel:>3}%) {cls}{tag}")
    p6 = next(x for x in sasa if x[0] == P_MUT)
    REPO_STAGE = "/content/cancer-recog-apoptosis/runs/rung30_kras_g12d/staging_crystal"
    meta = {"target": "kras_g12d_A1101", "allele": "HLA-A*11:01", "pep_mut": PEP, "pep_wt": "VVVGAGGVGK",
            "p_in_pep": P_MUT, "structure": "PDB 7OW6 (2.64A, crystal, TCR stripped)",
            "design_hotspot": f"B{P_MUT}", "hotspot_residue": f"p{P_MUT} Asp (G12D), only {p6[4]}% exposed in this crystal conformation",
            "design_note": "read-the-mutation: hotspot ON the G12D Asp (p6). KRAS WT IS presented (R32) so the "
                           "binder MUST discriminate G12D from G12 (Gly has no sidechain).",
            "mut_pdb": f"{REPO_STAGE}/kras_g12d_A1101_mut_pmhc.pdb",
            "cropped_design_target": f"{REPO_STAGE}/kras_g12d_A1101_mut_cropped.pdb",
            "EXPOSURE_CAVEAT": f"p6 G12D Asp is only ~{p6[4]}% solvent-exposed in 7OW6 (TCR-stripped) = NOT cleanly "
                      "up-facing, contradicting the screen's position-heuristic. The crystals are either TCR-BOUND "
                      "(7OW6/7PB2 = induced-fit conformation) or free with p5-p6 DISORDERED (7OW4) -> the free-pMHC "
                      "Asp exposure is genuinely uncertain. TCRs DO recognize G12D (gold-standard, R32b) so it IS "
                      "readable (likely via induced fit/charge), but a de novo binder targeting a partially-buried "
                      "charged residue is harder (echoes IDH1/BRAF). RECOMMEND: also fold the FREE pMHC (Protenix) "
                      "to compare conformations before/with the design.",
            "p6_sasa": {"resname": p6[1], "sasa_A2": p6[3], "rel_pct": p6[4], "class": p6[5]}}
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=2)
    print(f"\n[staged] {OUT}  (groove {len(groove)} aa, cropped to {len(crop)} aa; hotspot B{P_MUT} = G12D Asp)")


if __name__ == "__main__":
    main()
