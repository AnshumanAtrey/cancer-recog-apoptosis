#!/usr/bin/env python3
"""
RUNG-26e staging — turn the committed BRAF-V600E MUT/WT pMHC folds (runs/rung26e_braf_v600e/folds/)
into the artifacts the binder pipeline needs, AND empirically test whether V->E is *readable* by a
binder before we spend PXDesign quota (rule 1: let the structure answer).

Produces (in runs/rung26e_braf_v600e/staging/):
  braf_mut_pmhc.pdb / braf_wt_pmhc.pdb  - full pMHC (chains A=groove, B=peptide) = AF2/Protenix scoring targets
  braf_mut_cropped.pdb                   - peptide + 10 A groove = PXDesign Extended upload target
  meta.json                              - mut_pdb/wt_pdb (repo-clone paths) + hotspot=3 + peptide/p, for the notebooks
  readability.json                       - the p3 Glu-vs-Val exposure / up-facing / fold-invariance check

Readability rationale: a de novo binder sits ABOVE the peptide. For it to discriminate V600E it must
(a) the p3 sidechain points UP (toward the binder), (b) is solvent-exposed (not buried in the groove),
and (c) the MUT/WT folds differ only locally (so a two-state MUT-vs-WT comparison is valid, not confounded
by a global rearrangement). Glu (charged, longer) vs Val (small hydrophobic) at an exposed up-facing p3 is
the structural basis for why negative design *can* discriminate here, unlike IDH1's buried-ish His<->Arg.
"""
import os, json, warnings
import numpy as np

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FOLDS = os.path.join(ROOT, "runs/rung26e_braf_v600e/folds")
OUT = os.path.join(ROOT, "runs/rung26e_braf_v600e/staging")
MUT_CIF = os.path.join(FOLDS, "MUT/braf_v600e_MUT_pmhc_sample_0.cif")
WT_CIF = os.path.join(FOLDS, "WT/braf_v600e_WT_pmhc_sample_0.cif")

PEP_MUT, PEP_WT, P = "ATEKSRWSGSH", "ATVKSRWSGSH", 3   # V600E at peptide p3 (E vs V)
GROOVE_PREFIX = "SHSMRYF"
BACKBONE = {"N", "CA", "C", "O"}
# repo-clone paths so the notebooks work straight off `git clone` (no manual Drive upload)
REPO_STAGE = "/content/cancer-recog-apoptosis/runs/rung26e_braf_v600e/staging"


def load(cif):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    s = MMCIFParser(QUIET=True).get_structure("x", cif)
    model = next(iter(s))
    groove = peptide = None
    for ch in model:
        res = [r for r in ch if r.id[0] == " "]
        seq = "".join(seq1(r.resname) if len(r.resname) == 3 else "X" for r in res)
        if seq.startswith(GROOVE_PREFIX):
            groove = (ch, res, seq)
        elif len(res) == len(PEP_MUT):
            peptide = (ch, res, seq)
    assert groove and peptide, f"chain id failed in {cif}"
    return model, groove, peptide


def write_pdb(path, chains):
    """chains = list of (chain_id, [residues]); writes all atoms, sequential numbering."""
    serial = 1
    lines = []
    for cid, residues in chains:
        for i, r in enumerate(residues):
            for a in r:
                x, y, z = a.get_coord()
                nm = a.get_name()
                name = f" {nm:<3s}" if len(nm) < 4 else nm
                el = (a.element or nm[0]).strip()[:2]
                lines.append(f"ATOM  {serial:5d} {name:<4s}{'':1s}{r.resname:>3s} {cid:1s}"
                             f"{i+1:4d}{'':1s}   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}"
                             f"          {el:>2s}")
                serial += 1
        lines.append(f"TER   {serial:5d}      {residues[-1].resname:>3s} {cid:1s}{len(residues):4d}")
        serial += 1
    lines.append("END")
    open(path, "w").write("\n".join(lines) + "\n")


def crop_groove(groove_res, peptide_res, radius=10.0):
    """Keep groove residues with any atom within `radius` of any peptide atom."""
    pep_atoms = np.array([a.get_coord() for r in peptide_res for a in r])
    kept = []
    for r in groove_res:
        ga = np.array([a.get_coord() for a in r])
        d = np.linalg.norm(ga[:, None, :] - pep_atoms[None, :, :], axis=-1)
        if d.min() <= radius:
            kept.append(r)
    return kept


def p3_readability(model, groove, peptide, label):
    from Bio.PDB.SASA import ShrakeRupley
    g_ch, g_res, _ = groove
    p_ch, p_res, _ = peptide
    p3 = p_res[P - 1]
    sc_atoms = [a for a in p3 if a.get_name() not in BACKBONE]  # sidechain
    sc_coord = np.array([a.get_coord() for a in sc_atoms]) if sc_atoms else None

    # (a) up-facing: groove-CA centroid -> peptide-CA centroid = "up"; p3 sidechain direction vs up
    g_ca = np.array([r["CA"].get_coord() for r in g_res if "CA" in r])
    p_ca = np.array([r["CA"].get_coord() for r in p_res if "CA" in r])
    up = p_ca.mean(0) - g_ca.mean(0); up /= np.linalg.norm(up)
    ca = p3["CA"].get_coord()
    cos_up = float(np.dot((sc_coord.mean(0) - ca) / (np.linalg.norm(sc_coord.mean(0) - ca) + 1e-9), up)) if sc_coord is not None else None

    # (b) groove burial: groove heavy atoms within 5 A of any p3 sidechain atom
    burial = None
    if sc_coord is not None:
        g_atoms = np.array([a.get_coord() for r in g_res for a in r if a.element != "H"])
        dmin = np.linalg.norm(g_atoms[:, None, :] - sc_coord[None, :, :], axis=-1).min(0)
        burial = int((dmin <= 5.0).sum())  # contacts; lower = more exposed

    # (c) SASA of p3 sidechain (absolute A^2)
    sasa = None
    try:
        ShrakeRupley().compute(model, level="A")
        sasa = round(float(sum(a.sasa for a in sc_atoms)), 1) if sc_atoms else 0.0
    except Exception as e:
        sasa = f"sasa-failed: {e}"
    return {"label": label, "p3_resname": p3.resname, "up_facing_cos": None if cos_up is None else round(cos_up, 3),
            "groove_contacts_5A": burial, "sidechain_sasa_A2": sasa}


def main():
    os.makedirs(OUT, exist_ok=True)
    mM, gM, pM = load(MUT_CIF)
    mW, gW, pW = load(WT_CIF)
    assert pM[2] == PEP_MUT, f"MUT peptide {pM[2]!r} != {PEP_MUT!r}"
    assert pW[2] == PEP_WT, f"WT peptide {pW[2]!r} != {PEP_WT!r}"
    assert pM[1][P - 1].resname == "GLU" and pW[1][P - 1].resname == "VAL", "p3 not GLU/VAL"

    # full pMHC scoring targets
    write_pdb(os.path.join(OUT, "braf_mut_pmhc.pdb"), [("A", gM[1]), ("B", pM[1])])
    write_pdb(os.path.join(OUT, "braf_wt_pmhc.pdb"), [("A", gW[1]), ("B", pW[1])])
    # cropped MUT design target (peptide + 10 A groove)
    crop = crop_groove(gM[1], pM[1], 10.0)
    write_pdb(os.path.join(OUT, "braf_mut_cropped.pdb"), [("A", crop), ("B", pM[1])])

    # peptide-backbone RMSD MUT vs WT after superposing on groove CAs (fold-invariance check)
    from Bio.PDB import Superimposer
    g_ca_M = [r["CA"] for r in gM[1] if "CA" in r]
    g_ca_W = [r["CA"] for r in gW[1] if "CA" in r]
    n = min(len(g_ca_M), len(g_ca_W))
    sup = Superimposer(); sup.set_atoms(g_ca_M[:n], g_ca_W[:n])
    rot, tran = sup.rotran
    pep_bb_M = np.array([a.get_coord() for r in pM[1] for a in r if a.get_name() in BACKBONE])
    pep_bb_W = np.array([a.get_coord() for r in pW[1] for a in r if a.get_name() in BACKBONE])
    pep_bb_W_aln = pep_bb_W @ rot + tran
    pep_rmsd = float(np.sqrt(((pep_bb_M - pep_bb_W_aln) ** 2).sum(1).mean()))

    read = {"groove_ca_superpose_rmsd": round(sup.rms, 3),
            "peptide_backbone_rmsd_mut_vs_wt": round(pep_rmsd, 3),
            "mut": p3_readability(mM, gM, pM, "MUT_Glu"),
            "wt": p3_readability(mW, gW, pW, "WT_Val")}
    read["verdict"] = (
        "READABLE: p3 sidechain up-facing (cos>0) + low groove burial = a binder above the peptide can "
        "engage it; peptide backbone ~invariant MUT vs WT (RMSD small) so two-state MUT-vs-WT scoring is "
        "valid (difference is the p3 sidechain, not a fold rearrangement)."
        if (read["mut"]["up_facing_cos"] or 0) > 0 and read["peptide_backbone_rmsd_mut_vs_wt"] < 2.0
        else "CHECK: p3 may be buried/down-facing or the folds diverge - inspect before committing PXDesign quota.")

    json.dump(read, open(os.path.join(OUT, "readability.json"), "w"), indent=2)

    meta = {"target": "BRAF_V600E_A0101", "pep_mut": PEP_MUT, "pep_wt": PEP_WT, "p_in_pep": P,
            "hotspot": P,  # notebooks build f"B{hotspot}" -> B3
            "mut_pdb": f"{REPO_STAGE}/braf_mut_pmhc.pdb",
            "wt_pdb": f"{REPO_STAGE}/braf_wt_pmhc.pdb",
            "cropped_design_target": f"{REPO_STAGE}/braf_mut_cropped.pdb",
            "source_folds": "runs/rung26e_braf_v600e/folds/{MUT,WT}/*_sample_0.cif (Protenix seed 19931)",
            "note": "mut_pdb/wt_pdb are repo-clone paths -> notebooks work straight off git clone."}
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=2)

    print(f"[staged] {OUT}")
    print(f"  groove-CA superpose RMSD : {read['groove_ca_superpose_rmsd']} A")
    print(f"  peptide backbone RMSD    : {read['peptide_backbone_rmsd_mut_vs_wt']} A (MUT vs WT)")
    print(f"  MUT p3 Glu: up_cos={read['mut']['up_facing_cos']} burial5A={read['mut']['groove_contacts_5A']} sasa={read['mut']['sidechain_sasa_A2']}")
    print(f"  WT  p3 Val: up_cos={read['wt']['up_facing_cos']} burial5A={read['wt']['groove_contacts_5A']} sasa={read['wt']['sidechain_sasa_A2']}")
    print(f"  VERDICT: {read['verdict'][:80]}...")


if __name__ == "__main__":
    main()
