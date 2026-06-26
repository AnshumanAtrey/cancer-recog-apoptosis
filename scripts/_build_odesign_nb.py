#!/usr/bin/env python3
"""
Builder for notebooks/binder_odesign_kras_colab.ipynb — the KRAS-G12D / A*11:01 binder retry via ODesign
(the "different generator" lever; PXDesign was bounded at 0/80 because AF2-IG & Protenix anti-correlate).

Why a builder (rule 7): a GPU/Colab notebook CANNOT be run locally, so a single bad cell = a wasted user
round-trip. This script AST-CHECKS every code cell (magics `!`/`%`/`#@` neutralised to `pass`) before the
.ipynb is written, so the notebook is guaranteed to PARSE. Rule-7 hygiene baked in: GPU-guard at the top of
the install AND the run cell; a heartbeat thread (output-dir file count + GPU mem + elapsed) since inference
is one long subprocess; NO condacolab (pip only).

WORKFLOW (the project Drive convention — matches the other cancer-recon notebooks):
- Heavy assets (the ODesign checkpoints + the CCD file) are CACHED in YOUR Drive at
  `/content/drive/MyDrive/cancer-recon/odesign_assets/` -> they download ONCE; every later run reads them
  from Drive (no re-download, no re-pasting Google-Drive file IDs).
- The small run INPUTS (target pMHC PDB + input JSON) wget straight from the public repo (no upload).
- So after the first run, the only thing the notebook needs is a GPU.

Run:  python scripts/_build_odesign_nb.py     # writes the .ipynb + AST-checks every cell
"""
from __future__ import annotations
import ast
import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "binder_odesign_kras_colab.ipynb"

# ---- markdown cells -------------------------------------------------------
MD_INTRO = """# KRAS-G12D / A\\*11:01 binder retry — **ODesign** (the different-generator lever)

**Why this run.** PXDesign bounded the KRAS-G12D external binder at **0/80 dual-oracle passers** because AF2-IG
and Protenix **anti-correlate** (r=−0.41) — no design satisfies both. ODesign is an **all-atom interaction
world-model, NOT AF2-based** → a genuinely different generative prior, with **explicit epitope/hotspot control**
(force the binder onto the G12D Asp). The one principled shot left at the external key.
*(Repo: The-Institute-for-AI-Molecular-Design/ODesign, Apache-2.0, arXiv 2510.22304.)*

**Fully automatic — no uploads, no Google-Drive file IDs:**
- Target pMHC PDB + input JSON: **`wget` from the public repo**.
- **CCD file** (528MB): **`wget` from this repo's GitHub Release** (`odesign-ccd-v20240608`).
- ODesign **checkpoints**: from HuggingFace.
- All heavy assets are **cached in YOUR Drive** (`MyDrive/cancer-recon/odesign_assets/`) on the first run, so
  every later run reads them from Drive — no re-download.
- ⇒ Set a GPU and **Run-all**. That's it.

> **Honest banner (rule 7).** ODesign's stack is frozen for **Python 3.10**; Colab is 3.12, so the install
> builds an isolated **Python-3.10 venv** (`uv venv`) and runs inference under it. Verified the cp310 wheels
> exist (rdkit 2023.3.1 / torch 2.3.1+cu121 / pyg). I still couldn't execute the CUDA install on the M2, so if
> a wheel fights on first run, the maintainers' **container path** (`RUN_GUIDE.md` Path A) is the rock-solid fallback.

**Runtime:** set **Runtime → Change runtime type → T4 GPU** BEFORE running — the first cell refuses to proceed
without a GPU (a CPU fallback runs ~50× slower and silently wastes hours)."""

MD_RESTART = """### The install builds an isolated Python-3.10 venv (`/content/odv`)
ODesign's stack is frozen for Python 3.10, but Colab runs 3.12 — so the install puts everything in a **separate
3.10 venv** (`uv venv`) and leaves Colab's own Python alone. ⇒ **no "restart runtime" prompt, nothing to re-run.**
Inference runs under `/content/odv/bin/python`. The install is ~6–9 min. (Rule 7: no `condacolab` — the venv is
cleaner.) If a pyg/torch wheel ever fights, the container (RUN_GUIDE Path A) is the fallback."""

MD_SCORE = """## ✅ Next: the real test — MUT-vs-WT discrimination scoring (separate step)
ODesign wrote designed binders (`outputs/.../*.cif`) against the **MUT** pMHC. Binding isn't the win — it must
**discriminate**. For each top binder sequence:
1. Fold it against the **MUT** pMHC (`kras_g12d_A1101_free_mut_pmhc.pdb`) **and** the **WT** pMHC
   (`kras_g12d_A1101_free_wt_pmhc.pdb`) — on **both** Protenix and AF2-IG (the bound was dual-certification).
   Both PDBs are in the repo under `runs/rung30_kras_g12d/staging/` (the notebook can wget them too).
2. **Win = high on MUT, low on WT, on BOTH oracles.** A binder that grips MUT and WT equally fails (it would
   attack the WT peptide on normal tissue, R32).
3. Use the existing scoring path (`notebooks/binder_score_colab.ipynb` / `binder_specificity_*`).

**Either outcome is informative:** a dual-certified discriminating binder = a breakthrough artifact; another 0
= a strong second-generator confirmation the external-binder route is bounded in-silico → recognition rests on
the validated **internal key** (R27/R33/R34–40)."""

# ---- code cells (each is pure-python OR magics-on-their-own-lines) ---------
C_GPU_PRE = r'''# --- GPU GUARD (before install; torch not present yet, so probe nvidia-smi) ---
import subprocess
r = subprocess.run(["nvidia-smi"], capture_output=True, text=True)
if r.returncode != 0:
    raise SystemExit("NO GPU detected. Runtime -> Change runtime type -> T4 GPU, then re-run this cell.")
print(r.stdout[:600])
print("GPU present — proceeding to install.")'''

C_INSTALL = r'''%cd /content
!rm -rf /content/ODesign
!git clone --depth 1 https://github.com/The-Institute-for-AI-Molecular-Design/ODesign.git
%cd /content/ODesign
# --- correct ODesign's under-specified requirements.txt IN PLACE, then install ONCE, so the single resolve
#     enforces EVERY frozen pin simultaneously. VALIDATED locally on a real py3.10 + linux-cp310 resolve:
#     numpy 1.26.3 / scipy 1.15.2 / biopython 1.83 all stay bit-identical (-> torch 2.3.1 + rdkit 2023.3.1
#     ABI unbroken), and prody backtracks to the 2.3.1 compatible with all of them. Corrections found by a
#     local import-scan of their source (rule 7), NOT guesswork:
#       - biotite: their pin ==1.0.1 LACKS biotite.interface.rdkit.from_mol (added in >=1.2.0)  -> 1.2.0
#       - prody + addict: imported by ProteinMPNN/invfold but ABSENT from requirements.txt       -> append
#     flash_attn is deliberately NOT added: every use_flash/use_deepspeed_evo_attention/use_lma config
#     defaults false and is never set true (verified in configs), so its 20-30min nvcc compile buys nothing.
import re, pathlib
_req = pathlib.Path("/content/ODesign/requirements.txt")
_lines = [l for l in _req.read_text().splitlines()
          if not re.match(r"^(biotite|prody|addict)([=<>!~ ]|$)", l.strip())]
_req.write_text("\n".join(_lines + ["biotite==1.2.0", "prody", "addict"]) + "\n")
print("requirements.txt corrected: biotite==1.2.0, +prody, +addict (flash_attn skipped by design)")
# ODesign's stack is frozen for Py3.10 (rdkit==2023.3.1 / torch==2.3.1 have NO 3.12 wheels; Colab is 3.12).
# Build a Python-3.10 venv with uv and install the CORRECTED stack THERE; Colab's own Python is left
# untouched (no restart prompt). NO condacolab (rule 7). [~6-9 min]
!pip install -q uv
!uv venv --python 3.10 --seed /content/odv
!/content/odv/bin/pip install -q -r requirements.txt -f https://data.pyg.org/whl/torch-2.3.1+cu121.html'''

C_GPU_POST = r'''# --- GPU GUARD (check the 3.10 VENV's torch sees CUDA; fail LOUD before spending time) ---
import subprocess
chk = subprocess.run(["/content/odv/bin/python", "-c",
                      "import torch; print('venv torch', torch.__version__, '| CUDA', torch.version.cuda, "
                      "'| GPU', torch.cuda.get_device_name(0)); assert torch.cuda.is_available()"],
                     capture_output=True, text=True)
print(chk.stdout.strip()); print(chk.stderr[-800:] if chk.returncode else "")
assert chk.returncode == 0, "venv torch/CUDA check FAILED — see stderr (re-run install, or use container Path A)"'''

C_DRIVE = r'''# --- Mount YOUR Drive for a persistent asset cache (checkpoints + CCD download ONCE) ---
import os
try:
    from google.colab import drive
    drive.mount("/content/drive")
    ASSETS = "/content/drive/MyDrive/cancer-recon/odesign_assets"
    print("Drive mounted — asset cache:", ASSETS)
except Exception as e:
    ASSETS = "/content/odesign_assets"
    print("no Drive (", type(e).__name__, ") — assets EPHEMERAL (re-download each session):", ASSETS)
for d in (ASSETS + "/ckpt", ASSETS + "/data", "/content/ODesign/ckpt", "/content/ODesign/data"):
    os.makedirs(d, exist_ok=True)'''

C_CKPT_CCD = r'''# --- Checkpoints + CCD: read from Drive cache if present; else fetch ONCE and cache to Drive ---
import os, glob, shutil, subprocess
RELEASE = "https://github.com/AnshumanAtrey/cancer-recog-apoptosis/releases/download/odesign-ccd-v20240608"

# 1) checkpoints — ONLY the 3 that protein-design inference loads (verified vs
#    src/utils/inference/infer_runner.py: {infer_model_name}.pt @L53, oinvfold_{modality}.ckpt @L105,
#    v_48_020.pt @L88). The other 6 in get_odesign_ckpt.sh are ligand/NA-only -> skip (less disk, less
#    time, lower free-tier eviction risk). Per-file Drive cache.
HF = "https://huggingface.co/The-Institute-for-AI-Molecular-Design"
CKPTS = [("odesign_base_prot_flex.pt", HF + "/ODesign/resolve/main/ckpt/odesign_base_prot_flex.pt"),
         ("oinvfold_protein.ckpt", HF + "/OInvFold/resolve/main/oinvfold_protein.ckpt"),
         ("v_48_020.pt", "https://github.com/dauparas/ProteinMPNN/raw/main/vanilla_model_weights/v_48_020.pt")]
for name, url in CKPTS:
    dst, cache = "/content/ODesign/ckpt/" + name, ASSETS + "/ckpt/" + name
    if os.path.exists(cache):
        shutil.copy(cache, dst); print("ckpt (Drive cache):", name)
    else:
        print("ckpt download:", name, "..."); subprocess.run(["wget", "-q", "-O", dst, url])
        if os.path.exists(dst) and os.path.getsize(dst) > 100_000:
            shutil.copy(dst, cache)
    assert os.path.exists(dst) and os.path.getsize(dst) > 100_000, "checkpoint missing/too small: " + name
print("checkpoints ready (3 protein-design files) ✓")

# 2) CCD files — Drive cache, else wget from THIS repo's GitHub Release (automatic, no file IDs)
CIF = "/content/ODesign/data/components.v20240608.cif"
PKL = "/content/ODesign/data/components.v20240608.cif.rdkit_mol.pkl"
cif_c, pkl_c = ASSETS + "/data/" + os.path.basename(CIF), ASSETS + "/data/" + os.path.basename(PKL)
if os.path.exists(cif_c) and os.path.exists(pkl_c):
    shutil.copy(cif_c, CIF); shutil.copy(pkl_c, PKL)
    print("CCD: loaded from Drive cache ✓")
else:
    print("CCD: downloading from GitHub Release (one-time, ~528MB) ...")
    subprocess.run(["wget", "-q", "-O", CIF, RELEASE + "/components.v20240608.cif"])
    subprocess.run(["wget", "-q", "-O", PKL, RELEASE + "/components.v20240608.cif.rdkit_mol.pkl"])
    assert os.path.getsize(CIF) > 1_000_000 and os.path.getsize(PKL) > 1_000_000, "CCD download failed"
    shutil.copy(CIF, cif_c); shutil.copy(PKL, pkl_c)   # cache to Drive for next time
    print("CCD: downloaded from release + cached to Drive ✓")
print("ckpt dir:", os.listdir("/content/ODesign/ckpt"))'''

C_FETCH_DIRS = r'''# Pull the validated INPUT directly from the public repo (no manual upload — the files are tracked).
import os
os.makedirs("/content/ODesign/examples/protein_design/prot_binding_prot", exist_ok=True)
RAW = "https://raw.githubusercontent.com/AnshumanAtrey/cancer-recog-apoptosis/main"
print("fetching inputs from", RAW)'''

C_FETCH_DL = r'''!wget -q -O /content/ODesign/data/kras_g12d_A1101_free_mut_cropped.pdb {RAW}/runs/rung30_kras_g12d/staging/kras_g12d_A1101_free_mut_cropped.pdb
!wget -q -O /content/ODesign/examples/protein_design/prot_binding_prot/kras.json {RAW}/runs/rung30_kras_g12d/odesign/kras_odesign_input.json'''

C_FETCH_VERIFY = r'''import os, json
spec = json.load(open("/content/ODesign/examples/protein_design/prot_binding_prot/kras.json"))
print("input JSON OK:", json.dumps(spec, indent=2))
p = "/content/ODesign/data/kras_g12d_A1101_free_mut_cropped.pdb"
assert os.path.exists(p) and os.path.getsize(p) > 1000, "target PDB fetch failed (check the RAW url / branch)"
print("target PDB OK:", os.path.getsize(p), "bytes")'''

C_RUNLIB = r'''# --- inference runner: reset-guard + GPU-guard + heartbeat + CAPTURED output (rule 7) ---
# stderr is MERGED into stdout and both streamed live AND kept in a ring buffer, so a crash traceback is
# ALWAYS printed in the copyable cell output (the old `subprocess.run` left tracebacks in Colab's red box,
# which got lost on copy). num_workers=0 in the smoke gives a clean main-process traceback (no dataloader wrapper).
import os, time, threading, subprocess, collections

def run_inference(exp_name, seeds, n_sample, num_workers):
    if not os.path.isdir("/content/ODesign") or not os.path.exists("/content/odv/bin/python"):
        raise SystemExit("Runtime was reset (free-tier eviction wiped /content). Runtime -> Run all from the "
                         "top — the Drive cache makes the checkpoint/CCD steps fast this time.")
    import torch
    assert torch.cuda.is_available(), "NO GPU — abort (CPU would be ~50x slower and silent)."
    out = "/content/ODesign/outputs"
    os.makedirs(out, exist_ok=True)
    n0 = sum(len(fs) for _, _, fs in os.walk(out))
    cmd = ["/content/odv/bin/python", "./scripts/inference.py",
           "exp=train_odesign_base_prot_flex",
           "data_root_dir=./data", "ckpt_root_dir=./ckpt",
           "exp.infer_model_name=odesign_base_prot_flex",
           "exp.design_modality=protein",
           "exp.input_json_path=./examples/protein_design/prot_binding_prot/kras.json",
           "exp.exp_name=" + exp_name,
           "exp.seeds=" + seeds,
           "exp.model.sample_diffusion.N_sample=" + str(n_sample),
           "exp.use_msa=false", "exp.num_workers=" + str(num_workers),
           "exp.model.inference_noise_schedulers.coordinate.partial_diffusion.enable=false",
           "exp.model.inference_noise_schedulers.coordinate.partial_diffusion.snr=0.1"]
    print("RUN:", " ".join(cmd), flush=True)
    stop = threading.Event()
    def hb():
        t0 = time.time()
        while not stop.is_set():
            n = sum(len(fs) for _, _, fs in os.walk(out))
            try:
                free, total = torch.cuda.mem_get_info(); used = (total - free) / 1e9
            except Exception:
                used = -1.0
            print(f"[hb] {time.time()-t0:6.0f}s  outputs={n} files  GPU_used={used:.1f}GB", flush=True)
            stop.wait(30)
    th = threading.Thread(target=hb, daemon=True); th.start()
    ring = collections.deque(maxlen=120)
    proc = subprocess.Popen(cmd, cwd="/content/ODesign", stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        ring.append(line)
        print(line, end="", flush=True)
    proc.wait()
    stop.set(); time.sleep(0.1)
    n1 = sum(len(fs) for _, _, fs in os.walk(out))
    print(f"\ninference exit code: {proc.returncode} | new output files: {n1 - n0}", flush=True)
    if proc.returncode != 0:
        print("\n===== LAST 120 LINES (traceback) =====\n" + "".join(ring), flush=True)
    return proc.returncode == 0 and n1 > n0

print("run_inference() ready.")'''

C_SMOKE = r'''# --- SMOKE: 1 seed x 2 samples, num_workers=0 — PROVE the full pipeline (~5 min) before the 50-design sweep.
# Validate-before-spend (rule 5/10): we have NEVER yet seen ODesign emit a single design, so confirm the deps
# resolve + the generator actually writes outputs before committing GPU-hours that a free-tier eviction could kill.
ok = run_inference("kras_smoke", "[42]", 2, 0)
assert ok, "SMOKE FAILED — read the '===== LAST 120 LINES (traceback) =====' block printed above."
print("\n✅ SMOKE PASSED — ODesign wrote outputs. Safe to run the full sweep below.")'''

C_RUN = r'''# --- FULL sweep: 5 seeds x 10 samples = ~50 designs (run ONLY after the smoke passed) ---
ok = run_inference("kras_g12d_odesign_v1", "[42,123,777,2024,31337]", 10, 4)
print("full sweep produced outputs:", ok)'''

C_PERSIST = r'''# Persist outputs to your Drive BEFORE the runtime dies (rule 7)
import shutil, os
dst = ASSETS.rsplit("/", 1)[0] + "/odesign_kras_g12d_v1"   # MyDrive/cancer-recon/odesign_kras_g12d_v1
os.makedirs(dst, exist_ok=True)
shutil.make_archive(dst + "/outputs", "zip", "/content/ODesign/outputs")
print("saved:", dst + "/outputs.zip")'''


def code(src):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}


def md(src):
    return {"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)}


CELLS = [
    md(MD_INTRO),
    code(C_GPU_PRE),
    md(MD_RESTART),
    code(C_INSTALL),
    code(C_GPU_POST),
    code(C_DRIVE),
    code(C_CKPT_CCD),
    code(C_FETCH_DIRS),
    code(C_FETCH_DL),
    code(C_FETCH_VERIFY),
    code(C_RUNLIB),
    code(C_SMOKE),
    code(C_RUN),
    code(C_PERSIST),
    md(MD_SCORE),
]


def ast_check(cells):
    """Rule 7: AST-parse every code cell with magics (!,%,#@) neutralised to `pass` (preserving indent)."""
    errs = 0
    for i, c in enumerate(cells):
        if c["cell_type"] != "code":
            continue
        lines = []
        for ln in "".join(c["source"]).splitlines():
            st = ln.lstrip()
            if st.startswith(("!", "%", "#@")):
                lines.append((len(ln) - len(st)) * " " + "pass")
            else:
                lines.append(ln)
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as e:
            errs += 1
            print(f"  [FAIL] code cell #{i}: {e}")
    if errs == 0:
        print(f"  [PASS] all {sum(c['cell_type']=='code' for c in cells)} code cells AST-parse cleanly")
    return errs


def main():
    if ast_check(CELLS):
        raise SystemExit("cell(s) failed AST check — fix before writing the notebook")
    nb = {"cells": CELLS,
          "metadata": {"accelerator": "GPU",
                       "colab": {"provenance": [], "gpuType": "T4"},
                       "kernelspec": {"display_name": "Python 3", "name": "python3"},
                       "language_info": {"name": "python"}},
          "nbformat": 4, "nbformat_minor": 0}
    NB_PATH.parent.mkdir(parents=True, exist_ok=True)
    NB_PATH.write_text(json.dumps(nb, indent=1))
    print(f"  wrote {NB_PATH}  ({len(CELLS)} cells)")


if __name__ == "__main__":
    main()
