#!/usr/bin/env python3
"""
Step 1a — Boltz-2 structure PREDICTION for the oracle smoke test (CLOUD-GRADE).

Predicts protein-protein complexes of DR5 ECD (chain A) with several binders
(chain B) and saves Boltz's raw confidence outputs. It does NOT decide pass/fail
— that is done by scripts/02_interface_metrics.py, which computes the proper
interface metrics (interface pLDDT, pDockQ, interface PAE) from these outputs.

WHY THE SPLIT:
  An earlier version decided on raw ipTM and "failed" because ipTM is fooled by
  confident *nonspecific* docking (a scrambled binder scored HIGHER ipTM than the
  real one, while having much lower pLDDT). The field uses pLDDT-based interface
  metrics to tell real binders from non-binders, so prediction and decision are
  now separate steps.

WHY ipTM/affinity are NOT used here:
  - Boltz-2's affinity head only accepts a small-molecule ligand binder; our binder
    is a protein → categorically rejected.
  - ipTM alone is a poor binder/non-binder classifier (see above). See script 02.

COMPLEXES (binder = chain B against DR5 ECD = chain A):
  positive         TRAIL ectodomain (real DR5 binder, PDB-validated)
  negative         scrambled TRAIL ectodomain (composition-matched, unfoldable)
  negative_folded  hen lysozyme (well-folded NON-binder → isolates interface
                   specificity from mere foldability; the rigorous control)

RESUMABILITY:
  Each complex caches in runs/step1_boltz/<name>/. If a confidence JSON already
  exists for a complex it is SKIPPED. Delete runs/step1_boltz/<name>/ to redo one.

REQUIREMENTS:
  GPU (T4 16GB suffices). pip install boltz ; ~10GB weights auto-download once.

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

RECEPTOR_FASTA = "dr5_ecd_human.fasta"

# (subdir, binder_fasta, description). Order: positive first.
COMPLEXES = [
    ("positive",        "trail_ecd_human.fasta",       "TRAIL ectodomain — real DR5 binder"),
    ("negative",        "scrambled_ecd_control.fasta",  "scrambled TRAIL — composition-matched non-binder"),
    ("negative_folded", "lysozyme_control.fasta",       "hen lysozyme — well-folded non-binder (specificity control)"),
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,                    # Colab swallows stderr; logs go to stdout
)
log = logging.getLogger("step1a")


def read_fasta_sequence(path: Path) -> str:
    lines = path.read_text().strip().splitlines()
    return "".join(ln.strip() for ln in lines if not ln.startswith(">")).upper().replace(" ", "")


def write_boltz_yaml(yaml_path: Path, receptor_seq: str, binder_seq: str, name: str) -> None:
    """Two protein chains A (receptor) and B (binder). No affinity property."""
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    yaml_path.write_text(
        f"# Boltz-2 input for {name} (protein-protein; analysed via interface metrics)\n"
        f"sequences:\n"
        f"  - protein:\n      id: A\n      sequence: {receptor_seq}\n"
        f"  - protein:\n      id: B\n      sequence: {binder_seq}\n"
    )


def have_boltz() -> bool:
    return shutil.which("boltz") is not None


def find_confidence_jsons(out_dir: Path) -> list[Path]:
    return sorted(out_dir.rglob("confidence_*.json"))


def run_boltz(input_yaml: Path, out_dir: Path, diffusion_samples: int = 1) -> int:
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
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1)
        for line in proc.stdout:
            sys.stdout.write(line); sys.stdout.flush(); logf.write(line)
        rc = proc.wait()
    log.info("boltz exited rc=%d in %.1fs", rc, time.time() - t0)
    return rc


def tail_log(log_path: Path, n: int = 30) -> None:
    if not log_path.exists():
        log.warning("no boltz.log at %s", log_path); return
    lines = log_path.read_text().splitlines()
    log.error("--- last %d lines of %s ---", min(n, len(lines)), log_path.name)
    for ln in lines[-n:]:
        print(f"    {ln}", flush=True)
    log.error("--- end of %s ---", log_path.name)


@dataclass
class PredResult:
    name: str
    description: str
    iptm: Optional[float]
    ptm: Optional[float]
    complex_plddt: Optional[float]
    complex_iplddt: Optional[float]
    complex_ipde: Optional[float]
    confidence_score: Optional[float]
    raw_path: Optional[str]
    status: str   # DONE | CACHED | MISSING | FAILED


def parse_best_confidence(out_dir: Path, name: str, desc: str, cached: bool) -> PredResult:
    jsons = find_confidence_jsons(out_dir)
    if not jsons:
        return PredResult(name, desc, None, None, None, None, None, None, None, "MISSING")
    best, best_path = None, None
    for p in jsons:
        d = json.loads(p.read_text())
        if d.get("iptm") is None:
            continue
        if best is None or float(d["iptm"]) > float(best["iptm"]):
            best, best_path = d, p
    if best is None:
        best, best_path = json.loads(jsons[0].read_text()), jsons[0]
    return PredResult(
        name=name, description=desc,
        iptm=best.get("iptm"), ptm=best.get("ptm"),
        complex_plddt=best.get("complex_plddt"),
        complex_iplddt=best.get("complex_iplddt"),
        complex_ipde=best.get("complex_ipde"),
        confidence_score=best.get("confidence_score"),
        raw_path=str(best_path),
        status="CACHED" if cached else "DONE",
    )


def process_complex(name: str, desc: str, receptor: str, binder: str, out_dir: Path) -> PredResult:
    log.info("=" * 60)
    log.info("[%s] %s (out_dir=%s)", name, desc, out_dir.relative_to(PROJECT_ROOT))
    if find_confidence_jsons(out_dir):
        log.info("[%s] SKIP — confidence JSON already present", name)
        return parse_best_confidence(out_dir, name, desc, cached=True)
    yaml_path = out_dir / "input.yaml"
    write_boltz_yaml(yaml_path, receptor, binder, name)
    log.info("[%s] wrote YAML → %s", name, yaml_path.relative_to(PROJECT_ROOT))
    rc = run_boltz(yaml_path, out_dir)
    if rc != 0:
        log.error("[%s] boltz failed rc=%d", name, rc); tail_log(out_dir / "boltz.log")
        return PredResult(name, desc, None, None, None, None, None, None, None, "FAILED")
    res = parse_best_confidence(out_dir, name, desc, cached=False)
    if res.status == "MISSING":
        log.error("[%s] boltz rc=0 but no confidence JSON — inspect log", name)
        tail_log(out_dir / "boltz.log")
    return res


def main() -> int:
    log.info("cancer-recog-apoptosis — Step 1a — Boltz-2 PREDICTION (decision in script 02)")
    if not have_boltz():
        log.error("`boltz` not in PATH. Install: pip install boltz")
        return 2
    try:
        receptor = read_fasta_sequence(SEQ_DIR / RECEPTOR_FASTA)
    except FileNotFoundError as e:
        log.error("receptor FASTA missing: %s", e); return 4
    log.info("receptor DR5_ECD len=%d", len(receptor))

    results = []
    for sub, fasta, desc in COMPLEXES:
        path = SEQ_DIR / fasta
        if not path.exists():
            log.warning("[%s] binder FASTA %s not found — skipping this complex", sub, fasta)
            continue
        binder = read_fasta_sequence(path)
        log.info("[%s] binder %s len=%d", sub, fasta, len(binder))
        results.append(process_complex(sub, desc, receptor, binder, RUN_DIR / sub))

    log.info("=" * 60)
    log.info("PREDICTION SUMMARY (raw confidence; decision = script 02)")
    for r in results:
        log.info("[%s] status=%s iptm=%s plddt=%s iplddt=%s ipde=%s",
                 r.name, r.status, r.iptm, r.complex_plddt, r.complex_iplddt, r.complex_ipde)

    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(
        {"complexes": [asdict(r) for r in results],
         "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2))
    log.info("raw state saved → %s", STATE_PATH.relative_to(PROJECT_ROOT))

    missing = [r.name for r in results if r.status in ("MISSING", "FAILED")]
    if missing:
        log.error("complexes without usable output: %s — re-run this cell", missing)
        return 3
    if not any(r.name == "positive" for r in results):
        log.error("no positive complex predicted — check data/sequences/"); return 4
    log.info("✓ all complexes predicted. Next: python scripts/02_interface_metrics.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
