#!/usr/bin/env python3
"""
RUNG-28 — neoantigen TARGET-SELECTION screen (the step we skipped before hand-picking IDH1/BRAF).

Both hand-picked external-binder targets FAILED for target-geometry reasons, not method failure:
  IDH1-R132H  — His<->Arg subtle + binder wraps conserved core (RUNG-26c/d NULL x3 models)
  BRAF-V600E  — mutation is a BURIED p3 anchor (2% exposed) AND the peptide barely presents on
                A*01:01 (MHCflurry MUT 3971 nM / pres 0.025) (RUNG-26e audit)

A de-novo binder can discriminate a neoantigen only when the target satisfies a GEOMETRY/PRESENTATION
criterion we never checked:
  (A) the MUT peptide is STRONGLY PRESENTED (else there's nothing on the cell surface to bind), AND
  (B) either the mutated residue is SOLVENT-EXPOSED/up-facing (binder reads it directly) OR the
      mutation FLIPS presentation (MUT presents, WT doesn't) -> selectivity by presentation, binder
      just needs to bind MUT well.

This screen does the CHEAP half — (A) + the presentation-flip part of (B) — with MHCflurry across the
common HLA-I alleles, over the real driver panel (CDS-derived mutant peptides). The exposure half of (B)
needs a fold (Protenix) + SASA, run ONLY on the presentation winners (scripts/58 machinery).

Run with the mhcflurry venv:  /tmp/mhc/bin/python scripts/59_neoantigen_presentation_screen.py
"""
import os, json, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CDS_DIR = os.path.join(ROOT, "data/refs/cds")
OUT = os.path.join(ROOT, "runs/rung28_neoag_screen")

CODON = {  # RNA codon -> aa (from scripts/54; stops omitted -> translation halts)
 "AAA":"K","AAC":"N","AAG":"K","AAU":"N","ACA":"T","ACC":"T","ACG":"T","ACU":"T","AGA":"R","AGC":"S",
 "AGG":"R","AGU":"S","AUA":"I","AUC":"I","AUG":"M","AUU":"I","CAA":"Q","CAC":"H","CAG":"Q","CAU":"H",
 "CCA":"P","CCC":"P","CCG":"P","CCU":"P","CGA":"R","CGC":"R","CGG":"R","CGU":"R","CUA":"L","CUC":"L",
 "CUG":"L","CUU":"L","GAA":"E","GAC":"D","GAG":"E","GAU":"D","GCA":"A","GCC":"A","GCG":"A","GCU":"A",
 "GGA":"G","GGC":"G","GGG":"G","GGU":"G","GUA":"V","GUC":"V","GUG":"V","GUU":"V","UAC":"Y","UAU":"Y",
 "UCA":"S","UCC":"S","UCG":"S","UCU":"S","UGC":"C","UGG":"W","UGU":"C","UUA":"L","UUC":"F","UUG":"L","UUU":"F"}

# (gene, label, aa_pos, wt_codon, mut_codon) — verified driver hotspots (scripts/54)
DRIVERS = [
    ("KRAS","G12D",12,"GGU","GAU"), ("KRAS","G12V",12,"GGU","GUU"), ("KRAS","G12C",12,"GGU","UGU"),
    ("KRAS","G13D",13,"GGC","GAC"), ("NRAS","Q61R",61,"CAA","CGA"), ("NRAS","Q61K",61,"CAA","AAA"),
    ("BRAF","V600E",600,"GUG","GAG"), ("IDH1","R132H",132,"CGU","CAU"), ("IDH1","R132C",132,"CGU","UGU"),
    ("PIK3CA","E545K",545,"GAG","AAG"), ("PIK3CA","H1047R",1047,"CAU","CGU"), ("EGFR","L858R",858,"CUG","CGG"),
    ("TP53","R175H",175,"CGC","CAC"), ("TP53","R248Q",248,"CGG","CAG"), ("TP53","R273H",273,"CGU","CAU"),
    ("TP53","R248W",248,"CGG","UGG"),
]

# common HLA-I (broad global coverage; A & B, the high-expression class-I)
ALLELES = ["HLA-A*01:01","HLA-A*02:01","HLA-A*03:01","HLA-A*11:01","HLA-A*24:02","HLA-A*26:01",
           "HLA-B*07:02","HLA-B*08:01","HLA-B*15:01","HLA-B*35:01","HLA-B*40:01","HLA-B*44:02","HLA-B*57:01"]

KMERS = (8, 9, 10, 11)
AFF_PRESENTED = 500.0   # nM; MUT must be a real binder
PRES_PRESENTED = 0.5    # presentation score; likely surface-presented


def translate(cds):
    aa = []
    for i in range(0, len(cds) - 2, 3):
        c = cds[i:i+3]
        if c not in CODON:
            break  # stop codon
        aa.append(CODON[c])
    return "".join(aa)


def neoag_peptides():
    """For every driver: WT and MUT protein, all 8-11mers spanning the mutated residue."""
    rows = []
    for gene, label, pos, wt_cod, mut_cod in DRIVERS:
        cds = open(os.path.join(CDS_DIR, f"{gene}.txt")).read().strip()
        s = (pos - 1) * 3
        assert cds[s:s+3] == wt_cod, f"{gene} codon{pos} {cds[s:s+3]} != {wt_cod}"
        wt_aa, mut_aa = CODON[wt_cod], CODON[mut_cod]
        prot = translate(cds)
        assert prot[pos-1] == wt_aa, f"{gene} prot[{pos}]={prot[pos-1]} != {wt_aa}"
        mut_prot = prot[:pos-1] + mut_aa + prot[pos:]
        pi = pos - 1
        for k in KMERS:
            for start in range(max(0, pi - k + 1), min(pi, len(prot) - k) + 1):
                wt_pep, mut_pep = prot[start:start+k], mut_prot[start:start+k]
                if len(mut_pep) == k and mut_pep != wt_pep:
                    rows.append(dict(gene=gene, label=label, k=k, offset=pi - start,
                                     wt_pep=wt_pep, mut_pep=mut_pep))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = neoag_peptides()
    mut_peps = sorted({r["mut_pep"] for r in rows})
    wt_peps = sorted({r["wt_pep"] for r in rows})
    print(f"{len(rows)} (driver,window) pairs -> {len(mut_peps)} MUT + {len(wt_peps)} WT unique peptides "
          f"x {len(ALLELES)} alleles")

    from mhcflurry import Class1PresentationPredictor
    P = Class1PresentationPredictor.load()

    def score(peps):  # {(peptide, allele): (affinity, presentation)}
        out = {}
        for al in ALLELES:
            df = P.predict(peptides=peps, alleles=[al], verbose=0)
            for _, r in df.iterrows():
                out[(r["peptide"], al)] = (float(r["affinity"]), float(r["presentation_score"]))
        return out
    mut_sc, wt_sc = score(mut_peps), score(wt_peps)

    # best MUT window per (driver, allele), with its matched WT
    best = {}
    for r in rows:
        for al in ALLELES:
            ma, mp = mut_sc[(r["mut_pep"], al)]
            wa, wp = wt_sc[(r["wt_pep"], al)]
            key = (r["label"], al)
            cand = dict(driver=r["label"], gene=r["gene"], allele=al, k=r["k"], offset=r["offset"],
                        mut_pep=r["mut_pep"], wt_pep=r["wt_pep"],
                        mut_affinity_nM=round(ma, 1), mut_presentation=round(mp, 3),
                        wt_affinity_nM=round(wa, 1), wt_presentation=round(wp, 3),
                        affinity_ratio_wt_over_mut=round(wa / max(ma, 1e-6), 2),
                        presentation_drop=round(mp - wp, 3))
            if key not in best or mp > best[key]["mut_presentation"]:
                best[key] = cand
    hits = sorted(best.values(), key=lambda c: -c["mut_presentation"])

    presented = [h for h in hits if h["mut_affinity_nM"] <= AFF_PRESENTED and h["mut_presentation"] >= PRES_PRESENTED]
    flip = [h for h in presented if h["affinity_ratio_wt_over_mut"] >= 5.0]  # MUT presents, WT much weaker

    print(f"\n=== {len(presented)} (driver,allele) pairs where MUT is strongly presented "
          f"(<= {AFF_PRESENTED} nM AND pres >= {PRES_PRESENTED}) ===")
    for h in presented[:25]:
        flipt = "  <FLIP wt%.0fx" % h["affinity_ratio_wt_over_mut"] + ">" if h in flip else ""
        print(f"  {h['driver']:7s} {h['allele']:12s} {h['mut_pep']:11s} "
              f"MUT {h['mut_affinity_nM']:7.1f}nM/{h['mut_presentation']:.2f}  "
              f"WT {h['wt_affinity_nM']:8.1f}nM/{h['wt_presentation']:.2f}{flipt}")
    print(f"\n=== {len(flip)} of those also show a PRESENTATION FLIP (WT >= 5x weaker) "
          f"= differential-presentation targets (binder needs only to bind MUT) ===")
    for h in flip[:15]:
        print(f"  {h['driver']:7s} {h['allele']:12s} {h['mut_pep']:11s} "
              f"MUT {h['mut_affinity_nM']:.0f}nM  WT {h['wt_affinity_nM']:.0f}nM  ({h['affinity_ratio_wt_over_mut']:.0f}x)")

    json.dump({"tag": "rung28_neoag_presentation_screen", "n_drivers": len(DRIVERS), "alleles": ALLELES,
               "thresholds": {"mut_affinity_nM": AFF_PRESENTED, "mut_presentation": PRES_PRESENTED, "flip_ratio": 5.0},
               "presented": presented, "presentation_flip": flip, "all_best_per_driver_allele": hits},
              open(os.path.join(OUT, "screen.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/screen.json")


if __name__ == "__main__":
    main()
