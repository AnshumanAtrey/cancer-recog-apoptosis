# RUNG-27d — Alt-enzyme CRISPR rescue (SpCas9-NG, Cas12a) closes RUNG-27c's 2 misses → 7/7 wobble drivers DNA-addressable

RUNG-27c designed allele-specific SpCas9-**NGG** guides for 5/7 G>A-wobble drivers and flagged 2 misses —
**TP53-R175H** (no NGG PAM) and **TP53-R273H** (SNV in NGG seed position 16 = too distal) — as "need SpCas9-NG /
Cas12a." RUNG-27d scans those engineered enzymes on the same CDS-local windows.

## Result — both misses rescued; combined 7/7

| driver | SpCas9-NGG (R27c) | **SpCas9-NG** | Cas12a (TTTV) | rescued by |
|---|---|---|---|---|
| TP53-R175H | NONE (no PAM) | **SEED, seed pos 2**, `TGACGGAGGTTGTGAGGCAC` (+) | none | **SpCas9-NG** |
| TP53-R273H | DISTAL (seed 16) | **SEED, seed pos 1** (deepest), `CGGAACAGCTTTGAGGTGCA` (+) | MID (seed 7) | **SpCas9-NG** |

**Combining SpCas9-NGG + NG + Cas12a → 7/7 wobble drivers** (KRAS-G12D, KRAS-G13D, IDH1-R132H, PIK3CA-E545K,
TP53-R175H/R248Q/R273H) have a designed allele-specific guide. The wobble drivers RNA toehold can't sense
(RUNG-27b) are now **fully DNA-addressable** for the autonomous MHC-free self-destruct circuit.

## Why it works (and why it's expected, not a fluke)
Both rescues are **SEED**, not PAM-creating — correct biology: a G>A transition *destroys* a G and *makes* an
A, so it can't CREATE an NG (needs G) or a TTTV (needs T) PAM. SpCas9-NG's relaxed **NG** PAM has ~every-other-
base density, so a PAM almost always sits where the SNV falls in the PAM-proximal seed (a single seed mismatch
collapses Cas9 on the WT allele). That is exactly the engineered purpose of SpCas9-NG — so finding seed guides
here is the expected outcome, validated by the deep seed positions (1 and 2 = strongest single-mismatch
discrimination). Self-caught consistency: each guide's mutant base (A) sits at the reported seed position
(R273H: 3'-most/PAM-adjacent; R175H: 2nd-from-PAM).

## Honest framing & ceiling (rule 3/5)
- This is a **sequence-design completion**, NOT new biology — modest, expected closure of the R27c gap.
- **SpCas9-NG's relaxed PAM raises genome-wide OFF-target liability** vs NGG → "PAM + seed available" is
  necessary, NOT sufficient; real allele-specific cutting + a clean off-target profile = the wet-lab residual.
  NG/Cas12a on-target efficiency is context-dependent and generally below NGG → NGG guides (the 5 from R27c)
  remain preferable where they exist; NG is the fallback that *unlocks* R175H/R273H, at an off-target cost.
- CDS-local (Ensembl MANE); SEED near an exon edge needs intron-aware genomic confirmation.
- enAsCas12a (TTYN) / SpRY (near-PAMless) would relax further but worsen the off-target trade → not scanned.

*Result: `rung27d_alt_crispr.json` (per-driver, per-enzyme guides). Script + selftest 6/6: `scripts/56_crispr_alt_enzymes.py` (reuses scripts/54 drivers/CDS + scripts/55 window/NGG scan).*
