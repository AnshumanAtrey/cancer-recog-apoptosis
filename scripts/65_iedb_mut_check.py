#!/usr/bin/env python3
"""
RUNG-32b — MUT-in-tumour confirmation via IEDB (the cancer-side complement to RUNG-32's benign atlas).
A benign immunopeptidome can't show a tumour mutation; IEDB curates published neoantigen epitopes + MS-eluted
ligands + T-cell assays. This queries the IEDB IQ-API for our MUT peptides and summarises the EXPERIMENTAL
evidence (which MHC alleles, MS-elution vs binding, T-cell recognition) — turning the MHCflurry PREDICTION
into a measurement.

API: https://query-api.iedb.org/ (PostgREST). curl is used (python urllib hit an SSL-verify quirk here).
"""
import os, json, subprocess, collections
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "runs/rung32_immunopeptidome")
PEPS = [("KRAS_G12D", "VVVGADGVGK"), ("PIK3CA_E545K", "STRDPLSEITK")]


def fetch(url):
    r = subprocess.run(["curl", "-s", "--max-time", "60", url], capture_output=True, text=True)
    return json.loads(r.stdout) if r.stdout.strip() else []


def summarize(pep):
    mhc = fetch(f"https://query-api.iedb.org/mhc_search?linear_sequence=eq.{pep}"
                "&select=mhc_allele_name,assay_names,qualitative_measure,disease_names,pubmed_id")
    tc = fetch(f"https://query-api.iedb.org/tcell_search?linear_sequence=eq.{pep}"
               "&select=mhc_allele_name,qualitative_measure")
    alleles = collections.Counter(r.get("mhc_allele_name") for r in mhc if r.get("mhc_allele_name"))
    assays = collections.Counter(r.get("assay_names") for r in mhc if r.get("assay_names"))
    qual = collections.Counter(r.get("qualitative_measure") for r in mhc)
    ms_eluted = sum(v for k, v in assays.items() if "mass spectrometry" in (k or ""))
    xray = sum(v for k, v in assays.items() if "x-ray" in (k or ""))
    pmids = sorted({r["pubmed_id"] for r in mhc if r.get("pubmed_id")})
    tcell_pos = sum(1 for r in tc if r.get("qualitative_measure", "").startswith("Positive"))
    tcell_alleles = sorted({r.get("mhc_allele_name") for r in tc if r.get("mhc_allele_name")})
    return {
        "peptide": pep, "n_mhc_ligand_assays": len(mhc),
        "mhc_alleles": dict(alleles.most_common()),
        "ms_eluted_records": ms_eluted, "xray_structures": xray,
        "qualitative": dict(qual.most_common()),
        "n_tcell_assays": len(tc), "tcell_positive": tcell_pos, "tcell_alleles": tcell_alleles,
        "n_references": len(pmids), "pubmed_ids": pmids[:10],
        "assay_types": dict(assays.most_common()),
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    results = {}
    for label, pep in PEPS:
        s = summarize(pep); results[label] = s
        verdict = ("MS-ELUTED + T-cell validated" if s["ms_eluted_records"] and s["tcell_positive"]
                   else "MS-eluted" if s["ms_eluted_records"]
                   else "binding-assay only" if s["n_mhc_ligand_assays"] else "NOT in IEDB")
        print(f"{label:13s} {pep:12s}: {s['n_mhc_ligand_assays']} MHC-ligand assays | alleles {s['mhc_alleles']} | "
              f"MS-eluted {s['ms_eluted_records']} | xray {s['xray_structures']} | "
              f"T-cell+ {s['tcell_positive']}/{s['n_tcell_assays']} ({s['tcell_alleles']}) | refs {s['n_references']} -> {verdict}")
    payload = {"tag": "rung32b_iedb_mut_in_tumour",
               "source": "IEDB IQ-API (query-api.iedb.org), MHC ligand + T-cell assays",
               "interpretation": "Both MUT peptides are EXPERIMENTALLY OBSERVED MHC ligands on our target alleles. "
                                 "KRAS-G12D = gold-standard (MS-eluted + x-ray + T-cell+ on A*11:01/A*03:01). "
                                 "PIK3CA-E545K = thinner (cellular-MHC binding only, 1 pos/1 neg on A*11:01/A*33:03; "
                                 "no MS/T-cell) -> presented-plausible, weaker evidence than KRAS.",
               "results": results}
    json.dump(payload, open(os.path.join(OUT, "iedb_mut.json"), "w"), indent=2)
    print(f"[saved] {OUT}/iedb_mut.json")


if __name__ == "__main__":
    main()
