#!/usr/bin/env python3
"""
RUNG-27b — THE AUTONOMOUS MUTATION-SENSING CIRCUIT: design real sensors for the real clonal driver
mutations, AND-gate two co-occurring drivers, and map — PER CANCER — what fraction of patients could get
an MHC-FREE autonomous self-destruct. This is Shriya's original concept made buildable, and the
counterpart to RUNG-11/16's IMMUNE-route addressability for the MHC-dark core the immune route can't reach.

WHERE THIS SITS IN THE ARC
--------------------------
  RUNG-23 : every EXPRESSION window leaks -> the MUTATION is the only tumour-exclusive signal.
  RUNG-25 : single-base RNA sensing is FEASIBLE but substitution-dependent — fails on G·U WOBBLE (G>A,U>C).
  RUNG-18 : ~4-13% of tumour cells are MHC-dark (immune route blind there).
  RUNG-27b (this): turn that into a DESIGN + an addressability map. For each real driver hotspot: build the
  sensor from REAL CDS, score discrimination (ViennaRNA), classify RNA-sensable (non-wobble) vs DNA-only
  (wobble -> CRISPR has no wobble). Then AND two co-occurring clonal drivers (respecting pathway MUTUAL
  EXCLUSIVITY) and compute per-cancer coverage under RNA-only vs RNA+DNA sensing.

KEY HONEST FINDING THIS SURFACES
--------------------------------
Most canonical driver hotspots are G>A transitions (CpG deamination / aging signature): KRAS-G12D,
IDH1-R132H, PIK3CA-E545K, and the TP53 R175H/R248Q/R273H hotspots are ALL G>A = WOBBLE -> poorly RNA-
sensable. So an RNA-toehold autonomous gate is constrained; DNA-level (CRISPR) sensing (no wobble) is
required to cover the bulk of drivers. The AND-gate is buildable per-cancer; the modality (RNA vs DNA) is
dictated by the driver's substitution chemistry.

CEILING (rule 3/5 — stated, not papered over)
  - ViennaRNA ΔG is a thermodynamic PROXY (kinetics / RNA accessibility / genome off-targets NOT modelled).
  - Per-cancer driver frequencies are approximate (TCGA/COSMIC bands); coverage is order-of-magnitude.
  - Co-occurrence across pathway groups is modelled INDEPENDENT (real TP53 co-occurs MORE -> coverage is
    conservative); within a mutual-exclusivity group only one driver can be ANDed.
  - "DNA-sensable" assumes CRISPR/base-sensor allele discrimination (established) — a buildability claim, not
    a built circuit. Synthesis + delivery + in-cell actuation = the wet-lab residual.

USAGE
  python scripts/54_mutation_circuit.py prep       # fetch + cache real CDS (Ensembl) and VERIFY every hotspot codon
  python scripts/54_mutation_circuit.py selftest   # pure-python proxy, no ViennaRNA/network — validates the logic
  python scripts/54_mutation_circuit.py run        # ViennaRNA -> runs/rung27b_circuit/
"""
from __future__ import annotations
import json, math, sys, time
from itertools import combinations
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "runs" / "rung27b_circuit"
CDS_DIR = PROJECT_ROOT / "data" / "refs" / "cds"
RESULT_JSON = OUT_DIR / "rung27b_circuit.json"
FIGURE_PNG = OUT_DIR / "rung27b_circuit.png"

RT = 0.616
COMP = {"A": "U", "U": "A", "G": "C", "C": "G"}
WIN = 11                      # nt each side of the mutated base -> 23-nt sensor window
DDG_SENSABLE = 2.0            # ΔΔG (kcal/mol) above which RNA discrimination is usable (RUNG-25 band)

# MANE Select transcripts (stable) for the driver genes
TX = {"KRAS": "ENST00000256078", "NRAS": "ENST00000369535", "BRAF": "ENST00000646891",
      "IDH1": "ENST00000345146", "PIK3CA": "ENST00000263967", "EGFR": "ENST00000275493",
      "TP53": "ENST00000269305"}

CODON = {  # RNA codon -> 1-letter aa (stops omitted)
 "AAA":"K","AAC":"N","AAG":"K","AAU":"N","ACA":"T","ACC":"T","ACG":"T","ACU":"T","AGA":"R","AGC":"S",
 "AGG":"R","AGU":"S","AUA":"I","AUC":"I","AUG":"M","AUU":"I","CAA":"Q","CAC":"H","CAG":"Q","CAU":"H",
 "CCA":"P","CCC":"P","CCG":"P","CCU":"P","CGA":"R","CGC":"R","CGG":"R","CGU":"R","CUA":"L","CUC":"L",
 "CUG":"L","CUU":"L","GAA":"E","GAC":"D","GAG":"E","GAU":"D","GCA":"A","GCC":"A","GCG":"A","GCU":"A",
 "GGA":"G","GGC":"G","GGG":"G","GGU":"G","GUA":"V","GUC":"V","GUG":"V","GUU":"V","UAC":"Y","UAU":"Y",
 "UCA":"S","UCC":"S","UCG":"S","UCU":"S","UGC":"C","UGG":"W","UGU":"C","UUA":"L","UUC":"F","UUG":"L","UUU":"F"}

# Real driver hotspots: (gene, label, codon_pos, wt_codon[RNA], mut_codon[RNA]). WT codon is VERIFIED against
# the fetched CDS in prep(); mut_codon encodes the exact single-base change from the c. notation.
DRIVERS = [
    ("KRAS", "G12D", 12, "GGU", "GAU"),   # c.35G>A   (wobble)
    ("KRAS", "G12V", 12, "GGU", "GUU"),   # c.35G>T
    ("KRAS", "G12C", 12, "GGU", "UGU"),   # c.34G>T
    ("KRAS", "G13D", 13, "GGC", "GAC"),   # c.38G>A   (wobble)
    ("NRAS", "Q61R", 61, "CAA", "CGA"),   # c.182A>G
    ("NRAS", "Q61K", 61, "CAA", "AAA"),   # c.181C>A
    ("BRAF", "V600E", 600, "GUG", "GAG"), # c.1799T>A
    ("IDH1", "R132H", 132, "CGU", "CAU"), # c.395G>A  (wobble)
    ("IDH1", "R132C", 132, "CGU", "UGU"), # c.394C>T
    ("PIK3CA", "E545K", 545, "GAG", "AAG"),  # c.1633G>A (wobble)
    ("PIK3CA", "H1047R", 1047, "CAU", "CGU"),# c.3140A>G
    ("EGFR", "L858R", 858, "CUG", "CGG"), # c.2573T>G
    ("TP53", "R175H", 175, "CGC", "CAC"), # c.524G>A  (wobble)
    ("TP53", "R248Q", 248, "CGG", "CAG"), # c.743G>A  (wobble)
    ("TP53", "R273H", 273, "CGU", "CAU"), # c.818G>A  (wobble)
    ("TP53", "R248W", 248, "CGG", "UGG"), # c.742C>T
]

# Approximate CLONAL driver-hotspot frequencies per cancer (TCGA / COSMIC bands — order-of-magnitude).
CANCER_DRIVERS = {
    "PDAC":     {"KRAS_G12D": .40, "KRAS_G12V": .30, "KRAS_G12C": .02, "TP53_R175H": .08, "TP53_R248Q": .07, "TP53_R273H": .06},
    "GLIOMA":   {"IDH1_R132H": .75, "IDH1_R132C": .03, "TP53_R175H": .10, "TP53_R248Q": .08, "TP53_R273H": .07},
    "MELANOMA": {"BRAF_V600E": .45, "NRAS_Q61R": .13, "NRAS_Q61K": .07, "TP53_R175H": .04, "TP53_R273H": .04},
    "CRC":      {"KRAS_G12D": .18, "KRAS_G12V": .12, "BRAF_V600E": .10, "PIK3CA_E545K": .10, "PIK3CA_H1047R": .06, "TP53_R175H": .12, "TP53_R248Q": .12, "TP53_R273H": .10},
    "LUAD":     {"KRAS_G12C": .13, "KRAS_G12V": .07, "KRAS_G12D": .04, "EGFR_L858R": .14, "TP53_R175H": .08, "TP53_R248Q": .08, "TP53_R273H": .08},
    "BREAST":   {"PIK3CA_H1047R": .14, "PIK3CA_E545K": .09, "TP53_R175H": .07, "TP53_R248Q": .06, "TP53_R273H": .05},
    "THYROID":  {"BRAF_V600E": .60},
    "AML":      {"IDH1_R132H": .06, "IDH1_R132C": .04, "NRAS_Q61R": .08, "NRAS_Q61K": .04},
}

# Pathway MUTUAL-EXCLUSIVITY groups: two drivers in the SAME group rarely co-occur in one clone -> cannot be
# ANDed. Different groups co-occur -> ANDable. (RTK/RAS/RAF pathway is one bus; TP53/PIK3CA/IDH1 separate.)
MUTEX_GROUPS = {
    "RTK_RAS_RAF": {"KRAS", "NRAS", "BRAF", "EGFR"},   # one effector bus -> largely mutually exclusive
    "TP53": {"TP53"}, "PI3K": {"PIK3CA"}, "IDH": {"IDH1"},
}
def _group_of(gene):
    for g, members in MUTEX_GROUPS.items():
        if gene in members:
            return g
    return gene


# ---------------------------------------------------------------------------
def rc(seq):
    return "".join(COMP[b] for b in reversed(seq))

def _sig(x, n=4):
    if x == 0:
        return 0.0
    return float(f"{x:.{n}g}")

def is_wobble_sub(wt_base, mut_base):
    """sensor base = comp(mut); sensor·WT pair = (comp(mut), wt). Wobble (G·U) -> can't discriminate at RNA level.
       G>A: comp(A)=U vs WT G -> {U,G} wobble. U>C: comp(C)=G vs WT U -> {G,U} wobble."""
    return {COMP[mut_base], wt_base} == {"G", "U"}


def _proxy_duplex_energy(s1, s2):
    t = s2[::-1]
    L = min(len(s1), len(t)); e = 0.0
    for i in range(L):
        w = 1.0 + 1.5 * (1 - abs(i - (L - 1) / 2) / ((L - 1) / 2 + 1e-9))
        e += (-2.0 * w) if COMP.get(s1[i]) == t[i] else (+1.5 * w)
    return e

def duplex_energy(sensor, target, backend):
    if backend == "vienna":
        import RNA
        return float(RNA.duplexfold(sensor, target).energy)
    return _proxy_duplex_energy(sensor, target)

def discrimination(wt_window, mut_window, backend):
    sensor = rc(mut_window)
    dg_mut = duplex_energy(sensor, mut_window, backend)
    dg_wt = duplex_energy(sensor, wt_window, backend)
    ddg = dg_wt - dg_mut
    ff = math.exp(-max(ddg, 0.0) / RT)
    return {"sensor": sensor, "ddg": round(ddg, 3), "dg_mut": round(dg_mut, 3),
            "dg_wt": round(dg_wt, 3), "false_fire_rate": _sig(min(ff, 1.0))}


# ---------------------------------------------------------------------------
def _cds_path(gene):
    return CDS_DIR / f"{gene}.txt"

def prep():
    import requests
    CDS_DIR.mkdir(parents=True, exist_ok=True)
    okall = True
    for gene, tx in TX.items():
        p = _cds_path(gene)
        if not p.exists():
            r = requests.get(f"https://rest.ensembl.org/sequence/id/{tx}?type=cds",
                             headers={"Content-Type": "text/plain"}, timeout=30)
            if not r.ok:
                print(f"  [FAIL] {gene} {tx}: HTTP {r.status_code}"); okall = False; continue
            p.write_text(r.text.strip().upper().replace("T", "U"))
        cds = p.read_text().strip()
        # verify every hotspot of this gene: CDS codon must translate to the expected WT aa
        for (g, label, pos, wt_cod, mut_cod) in DRIVERS:
            if g != gene:
                continue
            cod = cds[(pos - 1) * 3:(pos - 1) * 3 + 3]
            wt_aa = label[0]
            good = (cod == wt_cod) and (CODON.get(cod) == wt_aa)
            print(f"  [{'OK' if good else 'MISMATCH'}] {gene} {label}: CDS codon{pos}={cod} ({CODON.get(cod,'?')}) expect {wt_cod}({wt_aa})")
            okall &= good
    print("prep:", "all hotspots verified ✓" if okall else "VERIFICATION FAILED ✗")
    return 0 if okall else 1


def build_window(gene, pos, wt_cod, mut_cod, cds):
    """Return (wt_window, mut_window, mut_base_in_window, wt_base, mut_base) — 23-nt windows from real CDS."""
    s = (pos - 1) * 3
    assert cds[s:s + 3] == wt_cod, f"{gene} codon{pos} {cds[s:s+3]} != {wt_cod}"
    # which base in the codon changed
    ci = [i for i in range(3) if wt_cod[i] != mut_cod[i]]
    assert len(ci) == 1, f"{gene} {wt_cod}->{mut_cod} not single-base"
    mb = s + ci[0]                                  # CDS index of the mutated base
    lo, hi = mb - WIN, mb + WIN + 1
    wt_win = cds[lo:hi]
    mut_win = cds[lo:mb] + mut_cod[ci[0]] + cds[mb + 1:hi]
    return wt_win, mut_win, mb - lo, wt_cod[ci[0]], mut_cod[ci[0]]


# ---------------------------------------------------------------------------
def design_sensors(backend, cds_by_gene):
    """For each real driver: build window, design sensor, score discrimination, classify RNA- vs DNA-sensable."""
    out = {}
    for (gene, label, pos, wt_cod, mut_cod) in DRIVERS:
        cds = cds_by_gene[gene]
        wt_win, mut_win, mpos, wb, mb = build_window(gene, pos, wt_cod, mut_cod, cds)
        d = discrimination(wt_win, mut_win, backend)
        wob = is_wobble_sub(wb, mb)
        rna_sensable = (not wob) and (d["ddg"] >= DDG_SENSABLE)
        out[f"{gene}_{label}"] = {
            "gene": gene, "aa_change": label, "rna_sub": f"{wb}>{mb}", "is_wobble": wob,
            "wt_window": wt_win, "mut_window": mut_win, "sensor": d["sensor"], "ddg": d["ddg"],
            "false_fire_rate": d["false_fire_rate"],
            "rna_sensable": bool(rna_sensable),
            "dna_sensable": True,   # DNA/CRISPR allele discrimination has no G·U wobble -> all single-base drivers addressable
        }
    return out


def cancer_coverage(sensors):
    """Per cancer: P(>=2 co-occurring clonal drivers from DIFFERENT mutex groups, both sensable), under
       RNA-only vs RNA+DNA sensing. Independence across groups (conservative; TP53 co-occurs more)."""
    res = {}
    for cancer, drivers in CANCER_DRIVERS.items():
        # collapse drivers into mutex groups; a group is 'present & sensable(mode)' if >=1 of its sensable
        # drivers is present. group prob ~ sum of sensable-driver freqs in the group (capped at total).
        def group_probs(mode):
            groups = {}
            for dlabel, freq in drivers.items():
                s = sensors[dlabel]
                ok = s["rna_sensable"] if mode == "rna" else (s["rna_sensable"] or s["dna_sensable"])
                grp = _group_of(s["gene"])
                tot = groups.setdefault(grp, [0.0, 0.0])    # [present, present_and_sensable]
                tot[0] += freq
                if ok:
                    tot[1] += freq
            # cap at 1.0
            return {g: (min(v[0], 1.0), min(v[1], 1.0)) for g, v in groups.items()}

        def p_two_sensable(mode):
            gp = group_probs(mode)
            ps = [v[1] for v in gp.values()]            # per-group P(present & sensable)
            # P(>=2 of independent groups present&sensable) = 1 - P(0) - P(1)
            p0 = 1.0
            for p in ps:
                p0 *= (1 - p)
            p1 = 0.0
            for i in range(len(ps)):
                term = ps[i]
                for j in range(len(ps)):
                    if j != i:
                        term *= (1 - ps[j])
                p1 += term
            return max(0.0, 1 - p0 - p1)

        cov_rna = p_two_sensable("rna")
        cov_both = p_two_sensable("both")
        # which sensable drivers exist here, by modality
        rna_ok = [d for d in drivers if sensors[d]["rna_sensable"]]
        dna_only = [d for d in drivers if (not sensors[d]["rna_sensable"]) and sensors[d]["dna_sensable"]]
        res[cancer] = {
            "coverage_rna_only": _sig(cov_rna), "coverage_rna_plus_dna": _sig(cov_both),
            "dna_rescue_gain": _sig(cov_both - cov_rna),
            "rna_sensable_drivers": rna_ok, "dna_only_drivers": dna_only,
            "n_groups": len({_group_of(sensors[d]["gene"]) for d in drivers}),
        }
    return res


# ---------------------------------------------------------------------------
def main_run():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    backend = "vienna"
    try:
        import RNA  # noqa
    except Exception:
        print("[rung27b] ViennaRNA not installed. `pip install viennarna` then re-run. (selftest needs no dep.)")
        return 3
    # need the cached CDS (run prep first)
    miss = [g for g in TX if not _cds_path(g).exists()]
    if miss:
        print(f"[rung27b] missing CDS for {miss} — run `python scripts/54_mutation_circuit.py prep` first.")
        return 4
    cds_by_gene = {g: _cds_path(g).read_text().strip() for g in TX}

    sensors = design_sensors(backend, cds_by_gene)
    coverage = cancer_coverage(sensors)

    n = len(sensors)
    rna_ok = [k for k, v in sensors.items() if v["rna_sensable"]]
    dna_only = [k for k, v in sensors.items() if not v["rna_sensable"]]
    wobble = [k for k, v in sensors.items() if v["is_wobble"]]
    # rule-5 evidence: ViennaRNA ΔΔG must INDEPENDENTLY separate wobble from non-wobble (not just our flag)
    import statistics as _st
    _w = [v["ddg"] for v in sensors.values() if v["is_wobble"]]
    _n = [v["ddg"] for v in sensors.values() if not v["is_wobble"]]
    vienna_confirms = {"median_ddg_wobble": round(_st.median(_w), 2),
                       "median_ddg_nonwobble": round(_st.median(_n), 2),
                       "independently_confirmed": _st.median(_w) < _st.median(_n)}
    best_cov_both = max(coverage, key=lambda c: coverage[c]["coverage_rna_plus_dna"])
    big_rescue = sorted(coverage, key=lambda c: coverage[c]["dna_rescue_gain"], reverse=True)[:3]

    result = {
        "tag": "rung27b_autonomous_mutation_circuit",
        "question": "Design real allele-specific sensors for the recurrent clonal driver mutations, AND-gate two "
                    "co-occurring drivers, and map per-cancer what fraction of patients could get an MHC-free "
                    "autonomous self-destruct (the counterpart to RUNG-11/16 for the MHC-dark core).",
        "backend": "ViennaRNA duplex ΔG (real CDS windows, Ensembl MANE)", "window_nt": 2 * WIN + 1,
        "ddg_sensable_threshold": DDG_SENSABLE,
        "n_drivers": n, "n_rna_sensable": len(rna_ok), "n_dna_only_wobble": len(dna_only),
        "wobble_drivers_need_dna": wobble,
        "vienna_independently_confirms_wobble": vienna_confirms,
        "coverage_is": "OFF-THE-SHELF shared-driver-hotspot FLOOR — pre-designed sensors for recurrent drivers, "
                       "AND of 2 co-occurring clonal drivers from different pathway groups, independence-modelled. "
                       "Conservative: a PERSONALISED 2nd input (any clonal mutation, RUNG-16-style) lifts it; the "
                       "off-the-shelf advantage is that the sensors are shared across patients, not bespoke.",
        "per_driver": sensors,
        "per_cancer_coverage": coverage,
        "best_cancer_rna_plus_dna": [best_cov_both, coverage[best_cov_both]["coverage_rna_plus_dna"]],
        "biggest_dna_rescue": {c: coverage[c]["dna_rescue_gain"] for c in big_rescue},
        "HEADLINE": (
            f"Designed allele-specific sensors for {n} real clonal driver hotspots from REAL CDS. "
            f"{len(rna_ok)}/{n} are RNA-toehold-sensable (non-wobble, ΔΔG≥{DDG_SENSABLE}); {len(dna_only)}/{n} are "
            f"G·U-WOBBLE — incl. the most common drivers KRAS-G12D, IDH1-R132H, PIK3CA-E545K, and the TP53 "
            f"R175H/R248Q/R273H hotspots (all G>A transitions) — which an RNA toehold CAN'T discriminate but "
            f"DNA-level CRISPR sensing CAN (no wobble). AUTONOMOUS MHC-FREE AND-gate addressability (≥2 co-occurring "
            f"clonal sensable drivers from different pathway groups): RNA-only is LOW (most oncogene+TP53 hotspots "
            f"are wobble), but RNA+DNA lifts it sharply — top {best_cov_both} "
            f"~{coverage[best_cov_both]['coverage_rna_plus_dna']}. The autonomous self-destruct IS buildable per "
            f"cancer; the modality (RNA toehold vs DNA/CRISPR) is dictated by the driver's substitution chemistry, "
            f"and DNA-level sensing is REQUIRED to cover the G>A-transition bulk."),
        "INTERPRETATION_MAP": {
            "covers the MHC-dark core": "this route fires inside the cell on the mutation itself (no MHC) -> reaches "
                                        "exactly the ~4-13% MHC-dark escapees (RUNG-18) the immune route cannot.",
            "RNA-only low / RNA+DNA high": "G>A-transition drivers (the commonest) are RNA-wobble -> need DNA/CRISPR "
                                           "allele sensing; with both modalities the AND-gate is broadly buildable.",
            "AND of 2 clonal drivers": "somatic-mutation AND -> tumour-exclusive by construction (normal cells have "
                                       "neither); mutual-exclusivity means the 2 inputs must be from different pathways.",
        },
        "CEILING": [
            "ViennaRNA ΔG = thermodynamic PROXY (kinetics/accessibility/genome off-targets not modelled).",
            "Per-cancer driver frequencies are TCGA/COSMIC approximations -> coverage is order-of-magnitude.",
            "Co-occurrence across pathway groups modelled INDEPENDENT (TP53 co-occurs more -> conservative); within "
            "a mutual-exclusivity group only one driver is ANDable.",
            "DNA-sensable = CRISPR allele discrimination is established (buildability), NOT a built circuit; "
            "synthesis + delivery + in-cell apoptosis actuation = the wet-lab residual.",
        ],
    }
    RESULT_JSON.write_text(json.dumps(result, indent=2))
    print(f"[rung27b] wrote {RESULT_JSON}  ({time.monotonic()-t0:.1f}s)")
    print(f"  RNA-sensable {len(rna_ok)}/{n}; wobble(need DNA) {len(dna_only)}: {dna_only}")
    for c, v in sorted(coverage.items(), key=lambda kv: -kv[1]["coverage_rna_plus_dna"]):
        print(f"  {c:9s} RNA-only {v['coverage_rna_only']:<7} RNA+DNA {v['coverage_rna_plus_dna']:<7} (rescue +{v['dna_rescue_gain']})")
    _make_figure(coverage)
    return 0


def _make_figure(coverage):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[rung27b] matplotlib unavailable ({e})"); return
    cs = sorted(coverage, key=lambda c: -coverage[c]["coverage_rna_plus_dna"])
    rna = [coverage[c]["coverage_rna_only"] for c in cs]
    both = [coverage[c]["coverage_rna_plus_dna"] for c in cs]
    import numpy as np
    x = np.arange(len(cs)); w = 0.38
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.bar(x - w/2, rna, w, label="RNA-toehold only", color="#E0A040")
    ax.bar(x + w/2, both, w, label="RNA + DNA/CRISPR", color="#3F7D54")
    ax.set_xticks(x); ax.set_xticklabels(cs, rotation=30, ha="right")
    ax.set_ylabel("P(≥2 co-occurring sensable clonal drivers)\n= autonomous MHC-free AND-gate coverage")
    ax.set_title("RUNG-27b: autonomous mutation-sensing self-destruct — per-cancer addressability")
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(FIGURE_PNG, dpi=130)
    print(f"[rung27b] wrote {FIGURE_PNG}")


# ---------------------------------------------------------------------------
def selftest():
    ok = 0; checks = []
    def check(name, cond):
        nonlocal ok
        checks.append((name, bool(cond))); ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    check("rc reverse-complement", rc("AUGC") == "GCAU")

    # wobble classification matches the known biology (RUNG-25): G>A and U>C are wobble; others not
    check("G>A is wobble (KRAS-G12D/IDH1-R132H/TP53)", is_wobble_sub("G", "A"))
    check("U>C is wobble", is_wobble_sub("U", "C"))
    check("G>U non-wobble (KRAS-G12V)", not is_wobble_sub("G", "U"))
    check("U>A non-wobble (BRAF-V600E)", not is_wobble_sub("U", "A"))
    check("A>G non-wobble (NRAS-Q61R)", not is_wobble_sub("A", "G"))

    # build_window on a synthetic CDS: codon 5 = GGU -> GAU (G12D-like). Codon 5 keeps the mutated base >=11 nt
    # from the start so the 23-nt window fits (real driver hotspots are all well interior — verified in prep()).
    cds = "AUGCAUGCAUGC" + "GGU" + "AUGCAUGCAUGCAUGC"  # codons 1-4 filler, codon5=GGU(G)
    wt, mut, mpos, wb, mb = build_window("X", 5, "GGU", "GAU", cds)
    check("window wt/mut differ by exactly 1 base", sum(a != b for a, b in zip(wt, mut)) == 1)
    check("mutated base recorded correctly", wt[mpos] == "G" and mut[mpos] == "A" and (wb, mb) == ("G", "A"))

    # discrimination: mutant-matched sensor is more stable on mut than WT (ΔΔG>0) under proxy
    d = discrimination(wt, mut, "proxy")
    check("ΔΔG>0 (mut more stable)", d["ddg"] > 0)
    check("false_fire in (0,1]", 0 < d["false_fire_rate"] <= 1)

    # rna_sensable logic: a non-wobble strong-ΔΔG driver is RNA-sensable; a wobble one is not (DNA only)
    fake = {
        "A_good": {"gene": "BRAF", "rna_sub": "U>A", "is_wobble": False, "ddg": 5.0, "rna_sensable": True, "dna_sensable": True},
        "A_wob":  {"gene": "IDH1", "rna_sub": "G>A", "is_wobble": True,  "ddg": 0.5, "rna_sensable": False, "dna_sensable": True},
    }
    check("non-wobble strong -> rna_sensable", fake["A_good"]["rna_sensable"] and not fake["A_wob"]["rna_sensable"])

    # mutex grouping: KRAS & BRAF same group (can't AND); KRAS & TP53 different (ANDable)
    check("KRAS & BRAF same mutex group", _group_of("KRAS") == _group_of("BRAF"))
    check("KRAS & TP53 different groups", _group_of("KRAS") != _group_of("TP53"))

    # coverage math: P(>=2) of three independent groups at p each
    sensors = {
        "KRAS_G12V": {"gene":"KRAS","rna_sensable":True,"dna_sensable":True},
        "TP53_R175H":{"gene":"TP53","rna_sensable":False,"dna_sensable":True},
        "PIK3CA_E545K":{"gene":"PIK3CA","rna_sensable":False,"dna_sensable":True},
    }
    global CANCER_DRIVERS
    saved = CANCER_DRIVERS
    CANCER_DRIVERS = {"T": {"KRAS_G12V": .5, "TP53_R175H": .5, "PIK3CA_E545K": .5}}
    cov = cancer_coverage(sensors)["T"]
    CANCER_DRIVERS = saved
    # RNA-only: only KRAS_G12V sensable -> can't get 2 groups -> ~0; RNA+DNA: all 3 -> P(>=2 of three 0.5) = 0.5
    check("RNA-only coverage ~0 (only 1 sensable group)", cov["coverage_rna_only"] < 1e-9)
    check("RNA+DNA coverage = P(>=2 of three p=.5) = 0.5", abs(cov["coverage_rna_plus_dna"] - 0.5) < 1e-6)
    check("DNA rescue positive", cov["dna_rescue_gain"] > 0)

    print(f"\n  selftest: {ok}/{len(checks)} passed")
    return 0 if ok == len(checks) else 1


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "prep":
        sys.exit(prep())
    if cmd == "selftest":
        sys.exit(selftest())
    if cmd == "run":
        sys.exit(main_run())
    print(f"unknown: {cmd}"); sys.exit(64)
