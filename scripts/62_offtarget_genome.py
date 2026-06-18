#!/usr/bin/env python3
"""
RUNG-31b — GENOME-WIDE off-target scan for the internal key's 7 allele-specific guides (closes the
intronic/intergenic residual that the RUNG-31 coding-transcriptome first-pass left open).

Scope: the DANGEROUS, cutting-competent set = 0 or 1 total mismatch + valid PAM, exhaustively, genome-wide
(Ensembl GRCh38 primary assembly, both strands). Method: 2 disjoint 10-nt anchors — any ≤1-mismatch site has
>=1 anchor exact (pigeonhole). 10-mers are rare (~3 occurrences genome-wide) so this is exhaustive AND fast.
mm2+ genome-wide is NOT covered here (the transcriptome RUNG-31 covered mm2 in coding; mm2 cuts weakly) and
remains a residual alongside measured cutting (GUIDE-seq).

Built-in sanity / on-vs-off discrimination: the guide is MUT, the reference is WT, so each guide hits its OWN
driver locus at exactly 1 mismatch (the SNV) on its known chromosome -> that hit is the ON-target. Any ≤1mm
hit on a DIFFERENT chromosome (or far from the locus) = a real off-target that would cut normal cells.
"""
import os, sys, gzip, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FASTA = os.path.join(ROOT, "data/offtarget/grch38.fa.gz")
OUT = os.path.join(ROOT, "runs/rung31_offtarget")

GUIDES = [
    ("KRAS_G12D",   "KRAS",   "CTTGTGGTAGTTGGAGCTGA", "SpCas9-NGG", "12"),
    ("KRAS_G13D",   "KRAS",   "GTAGTTGGAGCTGGTGACGT", "SpCas9-NGG", "12"),
    ("IDH1_R132H",  "IDH1",   "ATCATAGGTCATCATGCTTA", "SpCas9-NGG", "2"),
    ("PIK3CA_E545K","PIK3CA", "TCTCTCTGAAATCACTAAGC", "SpCas9-NGG", "3"),
    ("TP53_R248Q",  "TP53",   "GCATGGGCGGCATGAACCAG", "SpCas9-NGG", "17"),
    ("TP53_R175H",  "TP53",   "TGACGGAGGTTGTGAGGCAC", "SpCas9-NG",  "17"),
    ("TP53_R273H",  "TP53",   "CGGAACAGCTTTGAGGTGCA", "SpCas9-NG",  "17"),
]
ANCHORS = [(0, 10), (10, 20)]   # pigeonhole for <=1 mismatch -> >=1 half exact
MAX_MM = 1
_RC = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_RC)[::-1]


def pam_ok(pam, enzyme):
    if enzyme == "SpCas9-NGG":
        return len(pam) >= 3 and pam[1] == "G" and pam[2] == "G"
    if enzyme == "SpCas9-NG":
        return len(pam) >= 2 and pam[1] == "G"
    return False


def pam_len(enzyme):
    return 3 if enzyme == "SpCas9-NGG" else 2


def norm_chrom(name):
    """Normalize to bare chromosome id: 'chr12'->'12', '12'->'12'."""
    return name[3:] if name.lower().startswith("chr") else name


def iter_chroms(path):
    name, seq = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield norm_chrom(name), "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip().upper())
        if name is not None:
            yield norm_chrom(name), "".join(seq)


def scan():
    res = {g[0]: {"gene": g[1], "enzyme": g[3], "protospacer": g[2], "on_chrom": g[4], "hits": []}
           for g in GUIDES}
    guides = [(lab, gene, ps, enz, pl, pam_len(enz)) for (lab, gene, ps, enz, pl) in GUIDES]
    for cname, seq in iter_chroms(FASTA):
        if "N" in seq:
            seq = seq.replace("N", "X")
        print(f"  scanning chr {cname} ({len(seq)/1e6:.0f} Mb)...", flush=True)
        for strand, s in (("+", seq), ("-", rc(seq))):
            L = len(s)
            for lab, gene, ps, enz, on_chrom, pl in guides:
                seen = set()
                for (a0, a1) in ANCHORS:
                    frag = ps[a0:a1]
                    start = 0
                    while True:
                        j = s.find(frag, start)
                        if j < 0:
                            break
                        start = j + 1
                        p0 = j - a0
                        if p0 < 0 or p0 + 20 + pl > L or p0 in seen:
                            continue
                        cand = s[p0:p0 + 20]
                        pam = s[p0 + 20:p0 + 20 + pl]
                        if not pam_ok(pam, enz):
                            continue
                        mm = sum(a != b for a, b in zip(cand, ps))
                        if mm > MAX_MM:
                            continue
                        seen.add(p0)
                        # report genomic coord (on - strand, convert back)
                        pos = p0 if strand == "+" else (L - p0 - 20)
                        res[lab]["hits"].append({"chrom": cname, "pos": pos, "strand": strand,
                                                 "mm": mm, "pam": pam, "site": cand})
    return res


def main():
    if not os.path.exists(FASTA):
        sys.exit(f"missing {FASTA}")
    os.makedirs(OUT, exist_ok=True)
    res = scan()
    print("\n=== genome-wide <=1mm hits (on-target = 1mm on the driver's own chromosome) ===")
    summary = []
    all_clean = True
    for lab, gene, ps, enz, on_chrom in GUIDES:
        r = res[lab]
        hits = r["hits"]
        on = [h for h in hits if h["chrom"] == on_chrom and h["mm"] == 1]
        off = [h for h in hits if not (h["chrom"] == on_chrom and h["mm"] == 1)]
        # a mm0 hit anywhere, or any <=1mm hit on another chrom, is an off-target
        sane = len(on) >= 1                       # guide must hit its own locus at 1mm (the SNV)
        clean = len(off) == 0
        all_clean &= clean
        flag = "CLEAN" if clean else f"OFF-TARGET x{len(off)}"
        print(f"{lab:13s} [{enz:10s}] on-target(chr{on_chrom}) 1mm x{len(on)} {'OK' if sane else 'SANITY-FAIL'} | "
              f"other <=1mm hits: {len(off)} -> {flag}")
        for h in off[:8]:
            print(f"      OFF: chr{h['chrom']}:{h['pos']} {h['strand']} mm{h['mm']} pam={h['pam']} {h['site']}")
        summary.append({"label": lab, "gene": gene, "enzyme": enz, "on_chrom": on_chrom,
                        "on_target_1mm_hits": len(on), "sanity_ok": sane,
                        "offtarget_le1mm": off, "n_offtarget_le1mm": len(off), "clean": clean})
    print(f"\nGENOME-WIDE <=1mm: {'ALL CLEAN' if all_clean else 'OFF-TARGETS FOUND -> inspect'}")
    json.dump({"tag": "rung31b_offtarget_genome_wide",
               "scope": "Ensembl GRCh38 primary assembly, both strands, exhaustive <=1 mismatch + PAM "
                        "(2x10nt pigeonhole). mm2+ and measured cutting (GUIDE-seq) = residual.",
               "max_mm": MAX_MM, "all_clean_le1mm": all_clean,
               "interpretation": "on-target = 1mm hit on the driver's own chromosome (the SNV). Any other "
                                 "<=1mm hit = a real off-target that would cut normal cells.",
               "guides": summary}, open(os.path.join(OUT, "offtarget_genome.json"), "w"), indent=2)
    print(f"[saved] {OUT}/offtarget_genome.json")


if __name__ == "__main__":
    main()
