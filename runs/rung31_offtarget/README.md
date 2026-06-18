# RUNG-31 — internal-key off-target first-pass: the 7 allele-specific guides are clean of dangerous (≤1mm) off-targets in expressed genes

The internal CRISPR key was "done" on ON-target allele-specificity (RUNG-27c/d, 7/7 wobble drivers) but carried
a NAMED safety residual: a guide that *also* cuts another locus would fire the self-destruct in NORMAL cells →
breaks the normal-vs-cancer premise. This scans the 7 designed guides against the human coding transcriptome.

## Result — dangerous category EMPTY for all 7
| guide | enzyme | on-target (sanity) | off mm0 / mm1 / mm2 |
|---|---|---|---|
| KRAS-G12D | SpCas9-NGG | 1 mm ✓ | 0 / 0 / **1** (KRASP1, pseudogene) |
| KRAS-G13D | SpCas9-NGG | 1 mm ✓ | 0 / 0 / 0 |
| IDH1-R132H | SpCas9-NGG | 1 mm ✓ | 0 / 0 / 0 |
| PIK3CA-E545K | SpCas9-NGG | 1 mm ✓ | 0 / 0 / **1** (unnamed) |
| TP53-R248Q | SpCas9-NGG | 1 mm ✓ | 0 / 0 / 0 |
| TP53-R175H | **SpCas9-NG** | 1 mm ✓ | 0 / 0 / 0 |
| TP53-R273H | **SpCas9-NG** | 1 mm ✓ | 0 / 0 / 0 |

**No 0- or 1-mismatch off-target in any other gene** — the cutting-competent danger set is empty. The only
≤2mm hits are KRAS-G12D→**KRASP1** (a known KRAS pseudogene; mm2) and one unnamed mm2 for PIK3CA. The
**SpCas9-NG guides** (TP53 R175H/R273H), which we flagged for off-target liability from the relaxed NG PAM,
show **zero** ≤2mm off-targets in expressed genes — the concern doesn't materialize at this resolution.

## Why this is trustworthy (rule 5)
- **Sanity PASSED:** every guide hits its OWN gene at exactly 1 mismatch (= the designed SNV; the reference is
  WT) → the scan machinery works. (A first seed-exact version FAILED this — the SNV sits in the seed, so a
  seed-exact filter misses both on- and off-target; switched to a pigeonhole search exhaustive for ≤2 mm.)
- **Positive control:** KRAS-G12D correctly surfaces **KRASP1** (the KRAS pseudogene paralog) at mm2 → the scan
  finds real paralogous near-matches, so the empty mm0/mm1 set is a real absence, not a blind scan.

## Method
Pigeonhole: split each 20 nt protospacer into 3 parts; any ≤2-mismatch site has ≥1 part exact → exact-search
each part, reconstruct the 20mer, count total mismatches, require a valid PAM (NGG / NG) immediately 3'. Both
strands, de-duplicated. `scripts/61_crispr_offtarget_scan.py`; data = Ensembl GRCh38 cdna.all (gitignored).

## Honest residuals
- **mm3+** not exhaustively scanned (low cutting risk); **intronic/intergenic/regulatory** off-targets need a
  genome-wide tool (**Cas-OFFinder on GRCh38**) — the wet-lab-grade residual.
- Coding-transcriptome ≈ exonic genomic sequence; introns/UTRs-of-genomic-DNA not represented.
- The mm2 KRASP1 hit (pseudogene, weakly cutting at 2 mm) is low-risk but noted.
- Computational near-match ≠ measured cutting; GUIDE-seq / amplicon sequencing is the experimental residual.

**Net:** the internal key's allele-specific guides clear the dangerous (≤1mm, expressed-gene) off-target bar —
the normal-vs-cancer safety of the internal key holds at first-pass; genome-wide confirmation is the next tier.
