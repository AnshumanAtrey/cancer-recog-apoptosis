#!/usr/bin/env python3
"""
RUNG-26f - ProteinMPNN two-state negative-design harness (the principled fix RUNG-26d's
autopsy proved we need: positive design + a hotspot makes the binder CONTACT the mutation
but cannot make binding DEPEND on it; only explicit negative design - score/optimize the
binder against MUT *and against* WT - can discriminate).

This file is the M2/CPU half: it (a) SCORES a binder's preference for the MUT vs WT pMHC
context using ProteinMPNN conditional likelihoods, and (b) GENERATES new binder sequences
on a fixed backbone, two-state-selected by that preference. The Colab notebook reuses the
exact same functions on the full Drive backbone set and AF2-confirms the survivors.

Mechanic (no coordinates change between states - that is the point):
  Same backbone. The peptide's mutated position p4 is HIS in the MUT target and ARG in the
  WT target (IDH1 R132H: WT=Arg -> mut=His). ProteinMPNN uses backbone atoms for geometry
  and the *fixed* (visible) target residue identities as decoding context, so the binder
  chain's conditional NLL changes with that one identity. We design/score chain C (binder),
  fixing chains A (HLA-A*01:01 groove) + B (peptide).
    score = mean NLL over binder residues  (LOWER = ProteinMPNN finds the binder more
            likely given that target context = better "fit")
    dscore = NLL_WT - NLL_MUT  (POSITIVE => binder prefers the MUT context => discriminating)
  Averaged over K random decoding orders to damp permutation noise (report mean +/- sd).

VALIDATION GATE (rule 10 - prove the tool before spending GPU): run on the existing
PXDesign rank_1/2/3 binders, which AF2 + Protenix already showed are NON-specific
(bind MUT ~= WT). The harness must (i) run clean and (ii) return |dscore| ~ 0 for them -
i.e. it must NOT hallucinate discrimination where three other models found none. A harness
that "discovers" specificity in a known-null binder is broken, not groundbreaking.

Targets: peptide MUT=IIGHHAYGDQY / WT=IIGRHAYGDQY, mutation at p4 (His<-Arg); HLA-A*01:01.
"""
import os, sys, json, argparse, tempfile, warnings
import numpy as np

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MPNN_DIR = os.path.join(ROOT, ".tools", "ProteinMPNN")
sys.path.insert(0, MPNN_DIR)

# --- target definition (IDH1 R132H / HLA-A*01:01) -------------------------------------
PEP_MUT = "IIGHHAYGDQY"   # p4 = H (the R132H His)
PEP_WT  = "IIGRHAYGDQY"   # p4 = R (germline Arg)
P_MUT   = 4               # 1-indexed mutated peptide position
MUT_RESNAME, WT_RESNAME = "HIS", "ARG"
GROOVE_PREFIX = "SHSMRYF"  # HLA-A heavy-chain N-terminus, to identify the groove chain
BACKBONE_ATOMS = ["N", "CA", "C", "O"]

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
THREE2ONE = {  # standard 20 + a couple of common alts -> 1-letter
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","SEC":"C","HSD":"H","HSE":"H",
}


def read_complex(cif_path):
    """Return (groove, peptide, binder), each a list of (resname, {atom: (x,y,z)})
    with only backbone atoms, chains identified by content not by label."""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure("x", cif_path)
    model = next(iter(s))
    chains = []
    for ch in model:
        residues = []
        for r in ch:
            if r.id[0] != " ":
                continue
            atoms = {}
            for a in r:
                if a.get_name() in BACKBONE_ATOMS:
                    atoms[a.get_name()] = tuple(float(v) for v in a.get_coord())
            if all(bb in atoms for bb in BACKBONE_ATOMS):
                residues.append((r.resname.strip().upper(), atoms))
        if residues:
            seq = "".join(THREE2ONE.get(rn, "X") for rn, _ in residues)
            chains.append({"id": ch.id, "seq": seq, "res": residues})
    groove = peptide = binder = None
    for c in chains:
        if c["seq"].startswith(GROOVE_PREFIX):
            groove = c
        elif len(c["res"]) == len(PEP_MUT):
            peptide = c
    rest = [c for c in chains if c is not groove and c is not peptide]
    if len(rest) == 1:
        binder = rest[0]
    if not (groove and peptide and binder):
        raise RuntimeError(f"chain id failed in {cif_path}: "
                           f"{[(c['id'], len(c['res'])) for c in chains]}")
    # validate the mutation register
    assert peptide["seq"] == PEP_MUT, f"peptide {peptide['seq']!r} != MUT {PEP_MUT!r}"
    assert peptide["res"][P_MUT - 1][0] == MUT_RESNAME, \
        f"p{P_MUT} is {peptide['res'][P_MUT-1][0]}, expected {MUT_RESNAME}"
    return groove, peptide, binder


def _atom_line(serial, atom, resname, chain, resseq, xyz):
    name = f" {atom:<3s}" if len(atom) < 4 else atom
    el = atom[0]
    return (f"ATOM  {serial:5d} {name}{'':1s}{resname:>3s} {chain:1s}{resseq:4d}{'':1s}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")


def write_backbone_pdb(path, groove, peptide, binder, pep_p4_resname):
    """Clean backbone-only PDB: chain A=groove, B=peptide, C=binder. The peptide's p4
    residue is renamed to pep_p4_resname (HIS for MUT, ARG for WT); coords unchanged."""
    serial = 1
    lines = []
    chain_specs = [("A", groove["res"], None),
                   ("B", peptide["res"], pep_p4_resname),
                   ("C", binder["res"], None)]
    for chain, residues, p4name in chain_specs:
        for i, (resname, atoms) in enumerate(residues):
            rn = resname
            if p4name is not None and (i + 1) == P_MUT:
                rn = p4name
            for atom in BACKBONE_ATOMS:
                lines.append(_atom_line(serial, atom, rn, chain, i + 1, atoms[atom]))
                serial += 1
        lines.append(f"TER   {serial:5d}      {residues[-1][0]:>3s} {chain:1s}{len(residues):4d}")
        serial += 1
    lines.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


_MODEL = {"m": None}


def load_model(ckpt_name="v_48_020.pt", device=None):
    if _MODEL["m"] is not None:
        return _MODEL["m"]
    import torch
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    from protein_mpnn_utils import ProteinMPNN
    ckpt_path = os.path.join(MPNN_DIR, "vanilla_model_weights", ckpt_name)
    ckpt = torch.load(ckpt_path, map_location=device)
    model = ProteinMPNN(ca_only=False, num_letters=21, node_features=128, edge_features=128,
                        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                        augment_eps=0.0, k_neighbors=ckpt["num_edges"])
    model.to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    _MODEL["m"] = (model, device)
    return _MODEL["m"]


def score_binder(pdb_path, n_orders=16, seed=0):
    """Mean ProteinMPNN NLL over the binder chain (C), fixing target chains A+B,
    averaged over n_orders random decoding permutations. Returns (mean, sd)."""
    import torch
    from protein_mpnn_utils import tied_featurize, parse_PDB, _scores
    model, device = load_model()
    pdb = parse_PDB(pdb_path)
    name = pdb[0]["name"]
    chain_id_dict = {name: (["C"], ["A", "B"])}  # design C, fix A+B
    (X, S, mask, lengths, chain_M, chain_encoding_all, _cl, _vl, _ml, _mcl,
     chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, _tp, pssm_coef,
     pssm_bias, pssm_log_odds_all, bias_by_res_all, _tb) = tied_featurize(
        pdb, device, chain_id_dict, None, None, None, None, None, ca_only=False)
    mask_for_loss = mask * chain_M * chain_M_pos
    torch.manual_seed(seed)
    vals = []
    with torch.no_grad():
        for _ in range(n_orders):
            randn = torch.randn(chain_M.shape, device=device)
            log_probs = model(X, S, mask, chain_M * chain_M_pos, residue_idx,
                              chain_encoding_all, randn)
            sc = _scores(S, log_probs, mask_for_loss).cpu().numpy()[0]
            vals.append(float(sc))
    return float(np.mean(vals)), float(np.std(vals)), int(mask_for_loss.sum().item())


def two_state_score(cif_path, n_orders=16):
    """Score the binder in cif_path against MUT and WT pMHC contexts."""
    groove, peptide, binder = read_complex(cif_path)
    with tempfile.TemporaryDirectory() as td:
        mut_pdb = os.path.join(td, "mut.pdb")
        wt_pdb = os.path.join(td, "wt.pdb")
        write_backbone_pdb(mut_pdb, groove, peptide, binder, MUT_RESNAME)
        write_backbone_pdb(wt_pdb, groove, peptide, binder, WT_RESNAME)
        mut_m, mut_s, nb = score_binder(mut_pdb, n_orders=n_orders, seed=0)
        wt_m, wt_s, _ = score_binder(wt_pdb, n_orders=n_orders, seed=0)
    return {
        "binder_len": len(binder["res"]),
        "binder_scored_residues": nb,
        "nll_mut": round(mut_m, 4), "nll_mut_sd": round(mut_s, 4),
        "nll_wt": round(wt_m, 4), "nll_wt_sd": round(wt_s, 4),
        "dscore_wt_minus_mut": round(wt_m - mut_m, 4),
        # noise band: pooled sd of the two independent estimates
        "dscore_noise_sd": round((mut_s**2 + wt_s**2) ** 0.5, 4),
    }


def _featurize(pdb_path):
    from protein_mpnn_utils import tied_featurize, parse_PDB
    _model, device = load_model()
    pdb = parse_PDB(pdb_path)
    name = pdb[0]["name"]
    chain_id_dict = {name: (["C"], ["A", "B"])}  # design C, fix A+B
    feats = tied_featurize(pdb, device, chain_id_dict, None, None, None, None, None,
                           ca_only=False)
    return feats, device


def generate_two_state(cif_path, n_candidates=64, temperature=0.2, n_orders=4,
                       seed=0, top_k=8):
    """Sample n_candidates binder sequences on this fixed backbone conditioned on the MUT
    context, then two-state-score each (paired decoding orders) and rank by
    dscore = NLL_wt - NLL_mut (positive => prefers MUT => discriminating). This is the
    greedy negative-design version: oversample under MUT, keep the WT-disfavouring tail.
    The Colab notebook runs this on every Drive backbone and AF2-confirms the top_k."""
    import torch
    from protein_mpnn_utils import _scores, _S_to_seq
    model, device = load_model()
    groove, peptide, binder = read_complex(cif_path)
    with tempfile.TemporaryDirectory() as td:
        mut_pdb = os.path.join(td, "mut.pdb"); wt_pdb = os.path.join(td, "wt.pdb")
        write_backbone_pdb(mut_pdb, groove, peptide, binder, MUT_RESNAME)
        write_backbone_pdb(wt_pdb, groove, peptide, binder, WT_RESNAME)
        fM, _ = _featurize(mut_pdb)
        fW, _ = _featurize(wt_pdb)
    X, S_mut, mask, _len, chain_M, chain_enc = fM[0], fM[1], fM[2], fM[3], fM[4], fM[5]
    chain_M_pos, omit_AA_mask, residue_idx = fM[10], fM[11], fM[12]
    pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res = fM[15], fM[16], fM[17], fM[18]
    S_wt = fW[1]
    N = n_candidates

    def rep(t):
        return t.repeat(*([N] + [1] * (t.dim() - 1)))

    X_, S_mut_, mask_, chain_M_, chain_enc_, chain_M_pos_, residue_idx_ = map(
        rep, (X, S_mut, mask, chain_M, chain_enc, chain_M_pos, residue_idx))
    S_wt_ = rep(S_wt)
    omit_AA_mask_, pssm_coef_, pssm_bias_, pssm_log_odds_all_, bias_by_res_ = map(
        rep, (omit_AA_mask, pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res))
    pssm_log_odds_mask_ = (pssm_log_odds_all_ > 0.0).float()

    omit_AAs_np = np.array([aa in "X" for aa in ALPHABET], dtype=np.float32)
    bias_AAs_np = np.zeros(len(ALPHABET), dtype=np.float32)

    torch.manual_seed(seed)
    with torch.no_grad():
        randn = torch.randn(chain_M_.shape, device=device)
        sd = model.sample(X_, randn, S_mut_, chain_M_, chain_enc_, residue_idx_,
                          mask=mask_, temperature=temperature, omit_AAs_np=omit_AAs_np,
                          bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos_,
                          omit_AA_mask=omit_AA_mask_, pssm_coef=pssm_coef_,
                          pssm_bias=pssm_bias_, pssm_multi=0.0, pssm_log_odds_flag=0,
                          pssm_log_odds_mask=pssm_log_odds_mask_, pssm_bias_flag=0,
                          bias_by_res=bias_by_res_)
        S_sample = sd["S"]
        cmb = chain_M_.bool()
        S_eval_mut = torch.where(cmb, S_sample, S_mut_)
        S_eval_wt = torch.where(cmb, S_sample, S_wt_)
        mask_for_loss = mask_ * chain_M_ * chain_M_pos_
        nll_mut = torch.zeros(N, device=device)
        nll_wt = torch.zeros(N, device=device)
        for _ in range(n_orders):
            r = torch.randn(chain_M_.shape, device=device)  # paired order for both states
            lpm = model(X_, S_eval_mut, mask_, chain_M_ * chain_M_pos_, residue_idx_,
                        chain_enc_, r)
            lpw = model(X_, S_eval_wt, mask_, chain_M_ * chain_M_pos_, residue_idx_,
                        chain_enc_, r)
            nll_mut += _scores(S_eval_mut, lpm, mask_for_loss)
            nll_wt += _scores(S_eval_wt, lpw, mask_for_loss)
        nll_mut /= n_orders
        nll_wt /= n_orders
        dscore = (nll_wt - nll_mut).cpu().numpy()
        nll_mut_np = nll_mut.cpu().numpy()
        nll_wt_np = nll_wt.cpu().numpy()
        seqs = [_S_to_seq(S_sample[b], chain_M_[b]) for b in range(N)]

    order = np.argsort(-dscore)  # best (most MUT-preferring) first
    cands = [{
        "seq": seqs[b],
        "nll_mut": round(float(nll_mut_np[b]), 4),
        "nll_wt": round(float(nll_wt_np[b]), 4),
        "dscore_wt_minus_mut": round(float(dscore[b]), 4),
    } for b in order[:top_k]]
    return {
        "n_candidates": N, "temperature": temperature, "n_orders": n_orders,
        "dscore_max": round(float(dscore.max()), 4),
        "dscore_mean": round(float(dscore.mean()), 4),
        "dscore_std": round(float(dscore.std()), 4),
        "top": cands,
    }


def selftest():
    """Sanity-check writer/reader/identity on rank_1 without invoking torch."""
    cif = os.path.join(ROOT, "runs/rung26d_pxdesign_hsB4_idh1/top_designs/rank_1.cif")
    groove, peptide, binder = read_complex(cif)
    assert peptide["seq"] == PEP_MUT
    assert peptide["res"][P_MUT - 1][0] == "HIS"
    assert groove["seq"].startswith(GROOVE_PREFIX)
    assert 60 <= len(binder["res"]) <= 130
    with tempfile.TemporaryDirectory() as td:
        wt = os.path.join(td, "wt.pdb")
        write_backbone_pdb(wt, groove, peptide, binder, WT_RESNAME)
        # re-read the WT pdb via ProteinMPNN's own parser to confirm p4 flipped to R
        from protein_mpnn_utils import parse_PDB
        d = parse_PDB(wt)[0]
        pep_seq = d["seq_chain_B"]
        assert pep_seq == PEP_WT, f"written WT peptide {pep_seq!r} != {PEP_WT!r}"
        assert d["seq_chain_C"] == binder["seq"], "binder seq corrupted in write/read"
    print("[selftest] OK: chain ID, mutation register, His->Arg swap, MPNN re-parse all pass")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--generate", action="store_true",
                    help="negative-design generation demo on rank_1 backbone")
    ap.add_argument("--n", type=int, default=64, help="candidates to sample (generate mode)")
    ap.add_argument("--temp", type=float, default=0.2)
    ap.add_argument("--orders", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(ROOT, "runs/rung26f_negdesign/validate.json"))
    args = ap.parse_args()

    if args.selftest:
        selftest()
        sys.exit(0)

    if args.generate:
        cif = os.path.join(ROOT, "runs/rung26d_pxdesign_hsB4_idh1/top_designs/rank_1.cif")
        gen = generate_two_state(cif, n_candidates=args.n, temperature=args.temp,
                                 n_orders=max(4, args.orders // 4))
        out = os.path.join(ROOT, "runs/rung26f_negdesign/generate_demo.json")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        payload = {
            "tag": "rung26f_negdesign_generate_demo",
            "backbone": "rung26d rank_1 (PXDesign dual-passer)",
            "target": "IDH1_R132H_A0101  pep MUT=IIGHHAYGDQY / WT=IIGRHAYGDQY  p4",
            "noise_floor_nll": "~0.03-0.05 (from validate.json) - dscore must clear this to mean anything",
            **gen,
        }
        with open(out, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"generate demo: dscore max={gen['dscore_max']} mean={gen['dscore_mean']} "
              f"std={gen['dscore_std']}  (noise floor ~0.03-0.05)")
        for c in gen["top"][:3]:
            print(f"  dscore={c['dscore_wt_minus_mut']:+.4f}  nll_mut={c['nll_mut']} "
                  f"nll_wt={c['nll_wt']}  {c['seq']}")
        print(f"[written] {out}")
        sys.exit(0)

    designs = ["rank_1", "rank_2", "rank_3"]
    base = os.path.join(ROOT, "runs/rung26d_pxdesign_hsB4_idh1/top_designs")
    results = []
    for d in designs:
        cif = os.path.join(base, f"{d}.cif")
        if not os.path.exists(cif):
            print(f"[skip] {cif} missing"); continue
        r = two_state_score(cif, n_orders=args.orders)
        r["design"] = d
        results.append(r)
        print(f"{d}: NLL_mut={r['nll_mut']}+/-{r['nll_mut_sd']}  "
              f"NLL_wt={r['nll_wt']}+/-{r['nll_wt_sd']}  "
              f"dscore(WT-MUT)={r['dscore_wt_minus_mut']} (noise sd {r['dscore_noise_sd']})")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    payload = {
        "tag": "rung26f_negdesign_validate",
        "method": "ProteinMPNN v_48_020 two-state conditional NLL; design C, fix A+B; "
                  "MUT p4=HIS vs WT p4=ARG; mean over decoding orders",
        "n_decoding_orders": args.orders,
        "target": "IDH1_R132H_A0101  pep MUT=IIGHHAYGDQY / WT=IIGRHAYGDQY  p4",
        "interpretation": "dscore>0 => binder prefers MUT context (discriminating). "
                          "These rank_1/2/3 are AF2+Protenix-confirmed NON-specific binders, "
                          "so the GATE is |dscore| ~ 0 (harness must not invent discrimination).",
        "results": results,
    }
    with open(args.out, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"[written] {args.out}")
