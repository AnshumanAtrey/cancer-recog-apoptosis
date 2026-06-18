#!/usr/bin/env python3
"""
RUNG-31 — off-target first-pass for the internal key's allele-specific guides (RUNG-27c/d).

The internal CRISPR key is "done" (7/7 wobble drivers ON-target allele-specific) but we flagged a NAMED
residual: a guide that ALSO cuts some OTHER locus would fire the self-destruct in NORMAL cells = breaks the
whole normal-vs-cancer premise. This scans the 7 designed guides against the human CODING transcriptome
(Ensembl GRCh38 cDNA) for near-matches with a valid PAM in OTHER genes.

Method (PIGEONHOLE, exhaustive for <=2 mismatches): split the 20 nt protospacer into 3 parts; any site with
<=2 mismatches has >=1 part exact (2 mismatches can't touch all 3 parts). So we exact-search each part (fast
substring), reconstruct the full 20 nt at the implied offset, count TOTAL mismatches anywhere, and require a
valid PAM (NGG for SpCas9, NG for SpCas9-NG) immediately 3'. Both strands, de-duplicated. This catches
seed-mismatch off-targets too (unlike a seed-exact filter -- which, fatally, also misses the on-target because
the allele-specific SNV SITS IN the seed, so the MUT seed never exactly matches the WT reference).

Built-in sanity: the guide hits its OWN gene in the WT reference at ~1 mismatch (the designed SNV) -> proves
the machinery works + locates the on-target locus/PAM. Own-gene hits are reported as sanity, excluded from the
off-target tally. DANGEROUS off-target = a 0- or 1-mismatch hit in ANOTHER gene (would cut normal cells).

HONEST SCOPE: coding transcriptome + <=2 mismatches = the most consequential off-targets (cutting an expressed
gene). mm3+, and intronic/intergenic/regulatory sites, need a genome-wide tool (Cas-OFFinder on GRCh38) = the
residual.
"""
import os, sys, gzip, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FASTA = os.path.join(ROOT, "data/offtarget/cdna.fa.gz")
OUT = os.path.join(ROOT, "runs/rung31_offtarget")

# 7 designed allele-specific guides (MUT-targeting protospacers, 5'->3'), RUNG-27c (NGG) + 27d (NG)
GUIDES = [
    ("KRAS_G12D",   "KRAS",   "CTTGTGGTAGTTGGAGCTGA", "SpCas9-NGG"),
    ("KRAS_G13D",   "KRAS",   "GTAGTTGGAGCTGGTGACGT", "SpCas9-NGG"),
    ("IDH1_R132H",  "IDH1",   "ATCATAGGTCATCATGCTTA", "SpCas9-NGG"),
    ("PIK3CA_E545K","PIK3CA", "TCTCTCTGAAATCACTAAGC", "SpCas9-NGG"),
    ("TP53_R248Q",  "TP53",   "GCATGGGCGGCATGAACCAG", "SpCas9-NGG"),
    ("TP53_R175H",  "TP53",   "TGACGGAGGTTGTGAGGCAC", "SpCas9-NG"),
    ("TP53_R273H",  "TP53",   "CGGAACAGCTTTGAGGTGCA", "SpCas9-NG"),
]
ANCHORS = [(0, 7), (7, 14), (14, 20)]   # pigeonhole: <=2 mismatches -> >=1 anchor exact
MAX_MM = 2             # exhaustive up to 2 mismatches (the dangerous, cutting-competent set)
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


def iter_transcripts(path):
    """Yield (gene_symbol, sequence) from an Ensembl cDNA fasta.gz."""
    gene, seq = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if gene is not None:
                    yield gene, "".join(seq)
                gene = "?"
                for tok in line.split():
                    if tok.startswith("gene_symbol:"):
                        gene = tok.split(":", 1)[1]
                seq = []
            else:
                seq.append(line.strip().upper())
        if gene is not None:
            yield gene, "".join(seq)


def scan():
    results = {g[0]: {"gene": g[1], "enzyme": g[3], "protospacer": g[2],
                      "on_target_min_mm": None, "offtarget_by_mm": {0: 0, 1: 0, 2: 0},
                      "offtarget_genes": {}} for g in GUIDES}
    guides = [(label, gene, ps, enz, pam_len(enz)) for (label, gene, ps, enz) in GUIDES]

    n = 0
    for gsym, seq in iter_transcripts(FASTA):
        if "N" in seq:
            seq = seq.replace("N", "X")  # never matches ACGT
        for strand, s in (("+", seq), ("-", rc(seq))):
            L = len(s)
            for label, gene, ps, enz, pl in guides:
                r = results[label]
                seen = set()                               # dedupe offsets found via multiple anchors
                for (a0, a1) in ANCHORS:
                    frag = ps[a0:a1]
                    start = 0
                    while True:
                        j = s.find(frag, start)
                        if j < 0:
                            break
                        start = j + 1
                        ps_start = j - a0                  # implied protospacer 5' start
                        if ps_start < 0 or ps_start + 20 + pl > L or ps_start in seen:
                            continue
                        cand = s[ps_start:ps_start + 20]
                        pam = s[ps_start + 20:ps_start + 20 + pl]
                        if not pam_ok(pam, enz):
                            continue
                        mm = sum(a != b for a, b in zip(cand, ps))
                        if mm > MAX_MM:
                            continue
                        seen.add(ps_start)
                        if gsym == gene:
                            if r["on_target_min_mm"] is None or mm < r["on_target_min_mm"]:
                                r["on_target_min_mm"] = mm
                        else:
                            r["offtarget_by_mm"][mm] += 1
                            if gsym not in r["offtarget_genes"] or mm < r["offtarget_genes"][gsym]:
                                r["offtarget_genes"][gsym] = mm
        n += 1
    return results, n


def main():
    if not os.path.exists(FASTA):
        sys.exit(f"missing {FASTA} — download Ensembl GRCh38 cdna.all.fa.gz first")
    os.makedirs(OUT, exist_ok=True)
    results, n = scan()
    print(f"scanned {n} transcripts\n")
    summary = []
    sanity_ok = True
    for label, _, ps, enz in GUIDES:
        r = results[label]
        on = r["on_target_min_mm"]
        ot = r["offtarget_by_mm"]
        closest = next((mm for mm in (0, 1, 2) if ot[mm] > 0), None)
        worst_genes = sorted(r["offtarget_genes"].items(), key=lambda kv: kv[1])[:6]
        sane = on is not None and on <= 1          # guide must hit its own gene at <=1 mm (the SNV)
        sanity_ok &= sane
        print(f"{label:13s} [{enz:10s}] on-target({r['gene']}) min_mm={on} {'OK' if sane else 'FAIL'} | "
              f"off mm0/1/2 = {ot[0]}/{ot[1]}/{ot[2]} | closest off mm={closest} | top: {worst_genes}")
        summary.append({"label": label, "gene": r["gene"], "enzyme": enz,
                        "on_target_min_mm": on, "sanity_ok": sane, "offtarget_by_mm": ot,
                        "closest_offtarget_mm": closest, "n_offtarget_genes": len(r["offtarget_genes"]),
                        "closest_offtarget_genes": dict(worst_genes)})
    print(f"\nSANITY (all guides hit own gene at <=1mm): {'PASS' if sanity_ok else 'FAIL -> do not trust off-target counts'}")
    payload = {"tag": "rung31_offtarget_firstpass",
               "scope": "Ensembl GRCh38 coding transcriptome; pigeonhole exhaustive <=2 mismatches + PAM, both strands; "
                        "mm3+ and intronic/intergenic = residual (Cas-OFFinder on GRCh38).",
               "n_transcripts": n, "max_mm": MAX_MM, "sanity_all_pass": sanity_ok,
               "interpretation": "on_target_min_mm ~1 (guide vs WT reference = designed SNV) = sanity OK. "
                                 "DANGEROUS = a mm0/mm1 off-target in ANOTHER gene (would cut normal cells).",
               "guides": summary}
    json.dump(payload, open(os.path.join(OUT, "offtarget.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/offtarget.json")


if __name__ == "__main__":
    main()
