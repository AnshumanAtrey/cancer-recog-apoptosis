#!/usr/bin/env python3
"""
Step 1 — Boltz-2 oracle smoke test (CLOUD-GRADE).

Tests whether Boltz-2 can SEPARATE a real DR5 binder from a non-binder, using
the predicted INTERFACE confidence (ipTM) of two protein-protein complexes:

  (A) DR5 ECD + TRAIL ectodomain (residues 114-281)  ← positive control, real binder
  (B) DR5 ECD + scrambled TRAIL ectodomain            ← negative control, non-binder

WHY ipTM, NOT AFFINITY:
  Boltz-2's affinity head only accepts a SMALL-MOLECULE ligand as the binder
  (SMILES / CCD). A protein or peptide chain is categorically rejected
  ("Chain B is not a ligand! Affinity is currently only supported for ligands.").
  Our binder is a protein, so the correct binding proxy is the structural
  interface confidence ipTM — exactly how AlphaFold-Multimer is used to score
  protein-protein interactions.
  Refs: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
        Passaro et al., Boltz-2 (bioRxiv 2025.06.14.659707)

METRIC:
  ipTM = "Predicted TM score when aggregating at the interfaces", range [0,1],
  higher = more confident interface. From confidence_[name]_model_0.json.
  (AF-Multimer rule of thumb: ipTM > 0.6 confident, 0.4-0.6 uncertain, < 0.4 weak.)

DECISION RULE:
  PASS if ipTM(positive) - ipTM(negative) >= 0.15  → oracle discriminates → Step 2.
  Otherwise FAIL → see ASSESSMENT.md Day-1 kill criteria; pivot oracle stack.

BIOLOGY CAVEAT (logged, not fatal):
  Native TRAIL binds DR5 in the groove between TWO trimer protomers; we model the
  pairwise interface with a single TRAIL chain, which AF-class models usually
  recover from training data (PDB 1D4V/1DU3/1D0G are pre-cutoff). If ipTM(positive)
  is unexpectedly low (< 0.40) the pivot is to model the TRAIL homotrimer (3 chains)
  so the true groove is present.

RESUMABILITY:
  Each complex caches in runs/step1_boltz/<sub>/. If a confidence JSON already
  exists for a complex it is SKIPPED. To force a full rerun, delete
  runs/step1_boltz/state.json and the per-complex subdirs.

REQUIREMENTS:
  GPU (T4 16GB is enough for two ~300-residue complexes at diffusion_samples=1).
  pip install boltz ; ~10GB weights auto-download on first run.

USAGE:
  python scripts/01_boltz_smoketest.py
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SEQ_DIR = PROJECT_ROOT / "data" / "sequences"
RUN_DIR = PROJECT_ROOT / "runs" / "step1_boltz"
STATE_PATH = RUN_DIR / "state.json"

POSITIVE = RUN_DIR / "positive"
NEGATIVE = RUN_DIR / "negative"

# Decision threshold — see ASSESSMENT.md. ipTM margin between binder and non-binder.
MARGIN_THRESHOLD = 0.15
# Below this absolute ipTM, even a "winning" positive is a weak interface → trimer pivot.
WEAK_INTERFACE_IPTM = 0.40

# ---------- logging ----------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,                    # Colab swallows stderr; send logs to stdout
)
log = logging.getLogger("step1")


# ---------- helpers ----------
def read_fasta_sequence(path: Path) -> str:
    lines = path.read_text().strip().splitlines()
    seq_lines = [ln.strip() for ln in lines if not ln.startswith(">")]
    return "".join(seq_lines).upper().replace(" ", "")


def write_boltz_yaml(yaml_path: Path, receptor_seq: str, binder_seq: str, name: str) -> None:
    """Write a Boltz-2 input YAML for a protein-protein complex (NO affinity property).

    Two protein chains A (receptor) and B (binder). With --use_msa_server the MSAs
    are fetched automatically. We do NOT request affinity because the binder is a
    protein, not a small molecule.
    """
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    yaml_path.write_text(
        f"# Boltz-2 input for {name} (protein-protein; ipTM is the binding proxy)\n"
        f"sequences:\n"
        f"  - protein:\n      id: A\n      sequence: {receptor_seq}\n"
        f"  - protein:\n      id: B\n      sequence: {binder_seq}\n"
    )


def have_boltz() -> bool:
    return shutil.which("boltz") is not None


def find_confidence_jsons(out_dir: Path) -> list[Path]:
    """All confidence_*.json files Boltz wrote under this complex's out_dir."""
    return sorted(out_dir.rglob("confidence_*.json"))


def run_boltz(input_yaml: Path, out_dir: Path, diffusion_samples: int = 1) -> int:
    """Run boltz CLI, mirroring output to stdout AND a per-complex boltz.log."""
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "boltz.log"
    cmd = [
        "boltz", "predict", str(input_yaml),
        "--use_msa_server",
        "--diffusion_samples", str(diffusion_samples),
        "--out_dir", str(out_dir),
    ]
    log.info("invoking: %s", " ".join(cmd))
    log.info("boltz log mirrored to %s", log_path.relative_to(PROJECT_ROOT))
    t0 = time.time()
    with open(log_path, "w") as logf:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            logf.write(line)
        rc = proc.wait()
    log.info("boltz exited rc=%d in %.1fs", rc, time.time() - t0)
    return rc


def tail_log(log_path: Path, n: int = 30) -> None:
    if not log_path.exists():
        log.warning("no boltz.log at %s", log_path)
        return
    lines = log_path.read_text().splitlines()
    log.error("--- last %d lines of %s ---", min(n, len(lines)), log_path.name)
    for ln in lines[-n:]:
        print(f"    {ln}", flush=True)
    log.error("--- end of %s ---", log_path.name)


def list_out_dir(out_dir: Path) -> None:
    log.error("out_dir contents:")
    for p in sorted(out_dir.rglob("*")):
        print(f"    {p.relative_to(out_dir)}", flush=True)


@dataclass
class InterfaceResult:
    name: str
    iptm: Optional[float]
    ptm: Optional[float]
    complex_plddt: Optional[float]
    confidence_score: Optional[float]
    raw_path: Optional[str]
    status: str   # "DONE" | "CACHED" | "MISSING" | "FAILED"


def parse_best_confidence(out_dir: Path, name: str, cached: bool) -> InterfaceResult:
    """Parse confidence JSON(s); if several samples, keep the one with the highest ipTM."""
    jsons = find_confidence_jsons(out_dir)
    if not jsons:
        return InterfaceResult(name, None, None, None, None, None, "MISSING")
    best = None
    best_path = None
    for p in jsons:
        d = json.loads(p.read_text())
        iptm = d.get("iptm")
        if iptm is None:
            continue
        if best is None or float(iptm) > float(best.get("iptm", -1)):
            best = d
            best_path = p
    if best is None:                       # JSONs existed but had no iptm field
        d = json.loads(jsons[0].read_text())
        best, best_path = d, jsons[0]
    return InterfaceResult(
        name=name,
        iptm=best.get("iptm"),
        ptm=best.get("ptm"),
        complex_plddt=best.get("complex_plddt"),
        confidence_score=best.get("confidence_score"),
        raw_path=str(best_path),
        status="CACHED" if cached else "DONE",
    )


def process_complex(name: str, receptor: str, binder: str, out_dir: Path) -> InterfaceResult:
    """Idempotent: returns cached result if a confidence JSON already exists, else runs boltz."""
    log.info("=" * 60)
    log.info("[%s] start (out_dir=%s)", name, out_dir.relative_to(PROJECT_ROOT))
    if find_confidence_jsons(out_dir):
        log.info("[%s] SKIP — confidence JSON already present", name)
        return parse_best_confidence(out_dir, name, cached=True)

    yaml_path = out_dir / "input.yaml"
    write_boltz_yaml(yaml_path, receptor, binder, name)
    log.info("[%s] wrote YAML → %s", name, yaml_path.relative_to(PROJECT_ROOT))

    rc = run_boltz(yaml_path, out_dir)
    if rc != 0:
        log.error("[%s] boltz failed with rc=%d", name, rc)
        tail_log(out_dir / "boltz.log")
        list_out_dir(out_dir)
        return InterfaceResult(name, None, None, None, None, None, "FAILED")

    res = parse_best_confidence(out_dir, name, cached=False)
    if res.status == "MISSING":
        log.error("[%s] boltz rc=0 but no confidence JSON found — inspect log below", name)
        tail_log(out_dir / "boltz.log")
        list_out_dir(out_dir)
    return res


def save_state(positive: InterfaceResult, negative: InterfaceResult, margin: Optional[float]) -> None:
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps({
        "metric": "iptm",
        "positive": asdict(positive),
        "negative": asdict(negative),
        "iptm_margin": margin,
        "margin_threshold": MARGIN_THRESHOLD,
        "weak_interface_iptm": WEAK_INTERFACE_IPTM,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }, indent=2))
    log.info("state saved → %s", STATE_PATH.relative_to(PROJECT_ROOT))


# ---------- main ----------
def main() -> int:
    log.info("cancer-recon-apoptosis — Step 1 — Boltz-2 oracle smoke test (CLOUD, ipTM-based)")

    if not have_boltz():
        log.error("`boltz` CLI not found in PATH. Install: pip install boltz")
        log.error("For local plumbing test instead: python scripts/01_local_smoketest.py")
        return 2

    try:
        dr5_ecd   = read_fasta_sequence(SEQ_DIR / "dr5_ecd_human.fasta")
        trail_ecd = read_fasta_sequence(SEQ_DIR / "trail_ecd_human.fasta")
        scrambled = read_fasta_sequence(SEQ_DIR / "scrambled_ecd_control.fasta")
    except FileNotFoundError as e:
        log.error("FASTA missing: %s", e)
        return 4
    log.info("DR5_ECD len=%d  TRAIL_ECD len=%d  Scrambled len=%d",
             len(dr5_ecd), len(trail_ecd), len(scrambled))

    pos = process_complex("POSITIVE", dr5_ecd, trail_ecd, POSITIVE)
    neg = process_complex("NEGATIVE", dr5_ecd, scrambled,  NEGATIVE)

    log.info("=" * 60)
    log.info("RESULTS (metric = ipTM, interface confidence, range [0,1])")
    for r in (pos, neg):
        log.info("[%s] status=%s ipTM=%s pTM=%s complex_plddt=%s confidence=%s",
                 r.name, r.status, r.iptm, r.ptm, r.complex_plddt, r.confidence_score)

    if pos.iptm is None or neg.iptm is None:
        log.error("One or both complexes have no ipTM — cannot decide. See logs above.")
        save_state(pos, neg, None)
        return 3

    margin = float(pos.iptm) - float(neg.iptm)
    log.info("ipTM margin = ipTM(positive) - ipTM(negative) = %+.3f", margin)
    log.info("PASS threshold: margin >= %.2f", MARGIN_THRESHOLD)
    save_state(pos, neg, margin)

    if margin >= MARGIN_THRESHOLD:
        log.info("✅ PASS — oracle separates DR5 binder from non-binder. Proceed to Step 2.")
        if float(pos.iptm) < WEAK_INTERFACE_IPTM:
            log.warning("note: positive ipTM=%.3f is < %.2f (weak interface). It still beats the "
                        "negative, but for stronger signal model the TRAIL homotrimer (3 chains) "
                        "so the true DR5 groove is present.", float(pos.iptm), WEAK_INTERFACE_IPTM)
        return 0

    log.error("❌ FAIL — ipTM margin %.3f < %.2f. Oracle does not clearly discriminate.", margin, MARGIN_THRESHOLD)
    log.error("Pivot options (ASSESSMENT.md Day-1 kill criteria):")
    log.error("  (i)   model TRAIL homotrimer (3 chains) — presents the true DR5-binding groove")
    log.error("  (ii)  cross-check with AlphaFold 3 Server")
    log.error("  (iii) raise diffusion_samples (5) and take best-of-N ipTM")
    return 1


if __name__ == "__main__":
    sys.exit(main())
