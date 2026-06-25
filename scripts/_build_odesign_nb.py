#!/usr/bin/env python3
"""
Builder for notebooks/binder_odesign_kras_colab.ipynb — the KRAS-G12D / A*11:01 binder retry via ODesign
(the "different generator" lever; PXDesign was bounded at 0/80 because AF2-IG & Protenix anti-correlate).

Why a builder (rule 7): a GPU/Colab notebook CANNOT be run locally, so a single bad cell = a wasted user
round-trip. This script AST-CHECKS every code cell (magics `!`/`%`/`#@` neutralised to `pass`) before the
.ipynb is written, so the notebook is guaranteed to PARSE. It also encodes the rule-7 hygiene that has cost
us sessions before: GPU-guard at the top of the install AND the run cell; a heartbeat thread (output-dir file
count + GPU mem + elapsed) since inference is one long subprocess; NO condacolab (pip only); persist outputs
to Drive before the runtime dies; and the CCD files fetched by file-ID (NOT a folder download — that folder
also holds an 850 GB training tarball).

Run:  python scripts/_build_odesign_nb.py     # writes the .ipynb + AST-checks every cell
"""
from __future__ import annotations
import ast
import json
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent.parent / "notebooks" / "binder_odesign_kras_colab.ipynb"

MD = []  # not used as a list of all; we interleave below

# ---- markdown cells -------------------------------------------------------
MD_INTRO = """# KRAS-G12D / A\\*11:01 binder retry — **ODesign** (the different-generator lever)

**Why this run.** PXDesign bounded the KRAS-G12D external binder at **0/80 dual-oracle passers** because AF2-IG
and Protenix **anti-correlate** (r=−0.41) — no design satisfies both. ODesign is an **all-atom interaction
world-model, NOT AF2-based** → a genuinely different generative prior, with **explicit epitope/hotspot control**
(force the binder onto the G12D Asp). It is the one principled shot left at the external key.
*(Repo: The-Institute-for-AI-Molecular-Design/ODesign, Apache-2.0, arXiv 2510.22304.)*

> **Honest banner (rule 7).** This is a heavy CUDA-12.1 stack (torch 2.3.1 + pyg + deepspeed + a HF checkpoint
> + a Google-Drive CCD file). I could **not** dry-run the CUDA wheels on the M2 (platform-specific), so treat
> the first run as possibly needing install iteration. If pip fights you, the maintainers' **container path**
> (`Dockerfile`/`odesign.def` in the repo, see `RUN_GUIDE.md` Path A) is the lower-risk option on a CUDA box.

**Runtime:** set **Runtime → Change runtime type → T4 GPU** (or better) BEFORE running. The first cell will
refuse to proceed without a GPU (a CPU fallback would run ~50× slower and silently waste hours).

**What's the actual test.** ODesign only *generates* binders against the **MUT** pMHC. Binding ≠ winning — the
binder must **discriminate**: high on MUT, low on **WT**, on **both** oracles. The MUT-vs-WT scoring is the last
cell's pointer (a separate scoring step), not this notebook. 0 discriminating binders here too = a strong
2nd-generator confirmation the external route is bounded → recognition rests on the internal key (already the
validated contribution)."""

MD_RESTART = """### If Colab prompts to RESTART after the install
The full `requirements.txt` pins `numpy==1.26.3` / `protobuf==3.20.2` etc., which can downgrade Colab's base and
trigger a restart prompt. **That's fine** — click restart, then **re-run from the next cell down** (skip the
install cell). Do NOT re-run the install. (Rule 7: never `condacolab` here — it would wipe the pip installs.)"""

MD_CCD = """### CCD files — ONE manual step (and a trap to avoid)
ODesign needs the all-atom Chemical Component Dictionary: **`components.v20240608.cif`** and
**`components.v20240608.cif.rdkit_mol.pkl`**, from their Google Drive folder:
`https://drive.google.com/drive/folders/1wPmwIrC3G52q1JFY0RXY95tjKDl7YEln`

⚠️ **Do NOT `gdown --folder`** that link — the same folder holds an **850 GB** `odesign_full_data.tar.gz`
training tarball. Instead, open the folder, right-click **each** of the two CCD files → *Share → Copy link*,
take the ID from `.../d/<FILE_ID>/view`, and paste the two IDs in the cell below."""

MD_SCORE = """## ✅ Next: the real test — MUT-vs-WT discrimination scoring (separate step)
ODesign wrote designed binders (`outputs/.../*.cif`) against the **MUT** pMHC. A binder that binds isn't the
win — it must **discriminate**. For each top ODesign binder sequence:
1. Fold it against the **MUT** pMHC (`kras_g12d_A1101_free_mut_pmhc.pdb`) **and** the **WT** pMHC
   (`kras_g12d_A1101_free_wt_pmhc.pdb`) — on **both** Protenix and AF2-IG (the bound was dual-certification).
2. **Win = high on MUT, low on WT, on BOTH oracles.** A binder that grips MUT and WT equally fails (it would
   attack the WT peptide on normal tissue, R32).
3. Use the existing scoring path (`notebooks/binder_score_colab.ipynb` / `binder_specificity_*`).

**Either outcome is informative:** a dual-certified discriminating binder = a breakthrough artifact; another 0
= a strong second-generator confirmation the external-binder route is bounded in-silico → the recognition load
sits entirely on the **internal key** (R27/R33/R34–40), which is already the validated contribution.

Persist `outputs/` to Drive (cell above) BEFORE the runtime dies, then bring the winners into the scoring run."""

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
!pip install -q -r requirements.txt -f https://data.pyg.org/whl/torch-2.3.1+cu121.html'''

C_GPU_POST = r'''# --- GPU GUARD (after install; fail LOUD before spending time) ---
import torch
assert torch.cuda.is_available(), "NO GPU after install — set runtime to T4 and re-run from here."
print("torch", torch.__version__, "| CUDA", torch.version.cuda, "| GPU", torch.cuda.get_device_name(0))'''

C_CKPT = r'''%cd /content/ODesign
!bash ./ckpt/get_odesign_ckpt.sh ./ckpt
!ls -lh ./ckpt'''

C_CCD_IDS = r'''# Paste the two Google Drive file IDs (see the markdown above). Do NOT download the 850 GB tar.
CIF_FILE_ID = ""   # <-- components.v20240608.cif
PKL_FILE_ID = ""   # <-- components.v20240608.cif.rdkit_mol.pkl
import os
assert CIF_FILE_ID and PKL_FILE_ID, "paste BOTH Google Drive file IDs above before running the next cell"
os.makedirs("/content/ODesign/data", exist_ok=True)
print("ok — IDs set")'''

C_CCD_DL = r'''!pip install -q -U gdown
!gdown {CIF_FILE_ID} -O /content/ODesign/data/components.v20240608.cif
!gdown {PKL_FILE_ID} -O /content/ODesign/data/components.v20240608.cif.rdkit_mol.pkl
!ls -lh /content/ODesign/data'''

C_FETCH_DIRS = r'''# Pull the validated input DIRECTLY from the public repo (no manual upload — the files are tracked).
import os
os.makedirs("/content/ODesign/data", exist_ok=True)
os.makedirs("/content/ODesign/examples/protein_design/prot_binding_prot", exist_ok=True)
RAW = "https://raw.githubusercontent.com/AnshumanAtrey/cancer-recog-apoptosis/main"
print("fetching from", RAW)'''

C_FETCH_DL = r'''!wget -q -O /content/ODesign/data/kras_g12d_A1101_free_mut_cropped.pdb {RAW}/runs/rung30_kras_g12d/staging/kras_g12d_A1101_free_mut_cropped.pdb
!wget -q -O /content/ODesign/examples/protein_design/prot_binding_prot/kras.json {RAW}/runs/rung30_kras_g12d/odesign/kras_odesign_input.json'''

C_FETCH_VERIFY = r'''import os, json
spec = json.load(open("/content/ODesign/examples/protein_design/prot_binding_prot/kras.json"))
print("input JSON OK:", json.dumps(spec, indent=2))
p = "/content/ODesign/data/kras_g12d_A1101_free_mut_cropped.pdb"
assert os.path.exists(p) and os.path.getsize(p) > 1000, "target PDB fetch failed (check the RAW url / branch)"
print("target PDB OK:", os.path.getsize(p), "bytes")'''

C_RUN = r'''# --- INFERENCE with GPU guard + heartbeat (inference is one long subprocess; rule 7) ---
import torch, os, time, threading, subprocess
assert torch.cuda.is_available(), "NO GPU — abort (CPU would be ~50x slower and silent)."
os.chdir("/content/ODesign")
OUT = "/content/ODesign/outputs"
os.makedirs(OUT, exist_ok=True)

def heartbeat(stop):
    t0 = time.time()
    while not stop.is_set():
        n = sum(len(fs) for _, _, fs in os.walk(OUT))
        try:
            free, total = torch.cuda.mem_get_info(); used = (total - free) / 1e9
        except Exception:
            used = -1.0
        print(f"[hb] {time.time()-t0:6.0f}s  outputs={n} files  GPU_used={used:.1f}GB", flush=True)
        stop.wait(30)

cmd = ["python", "./scripts/inference.py",
       "exp=train_odesign_base_prot_flex",
       "data_root_dir=./data", "ckpt_root_dir=./ckpt",
       "exp.infer_model_name=odesign_base_prot_flex",
       "exp.design_modality=protein",
       "exp.input_json_path=./examples/protein_design/prot_binding_prot/kras.json",
       "exp.exp_name=kras_g12d_odesign_v1",
       "exp.seeds=[42,123,777,2024,31337]",
       "exp.model.sample_diffusion.N_sample=10",
       "exp.use_msa=false", "exp.num_workers=4",
       "exp.model.inference_noise_schedulers.coordinate.partial_diffusion.enable=false",
       "exp.model.inference_noise_schedulers.coordinate.partial_diffusion.snr=0.1"]
print("RUN:", " ".join(cmd), flush=True)
stop = threading.Event()
th = threading.Thread(target=heartbeat, args=(stop,), daemon=True); th.start()
proc = subprocess.run(cmd)
stop.set(); time.sleep(0.1)
print("inference.py exit code:", proc.returncode,
      "| outputs:", sum(len(fs) for _, _, fs in os.walk(OUT)), "files")
assert proc.returncode == 0, "inference failed — read the traceback above (often a CCD/ckpt path or a pyg/torch mismatch)."'''

C_PERSIST = r'''# Persist outputs to YOUR Drive BEFORE the runtime dies (rule 7)
from google.colab import drive
drive.mount("/content/drive")
import shutil, os
dst = "/content/drive/MyDrive/odesign_kras_g12d_v1"
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
    code(C_CKPT),
    md(MD_CCD),
    code(C_CCD_IDS),
    code(C_CCD_DL),
    code(C_FETCH_DIRS),
    code(C_FETCH_DL),
    code(C_FETCH_VERIFY),
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
        src = "\n".join(lines)
        try:
            ast.parse(src)
        except SyntaxError as e:
            errs += 1
            print(f"  [FAIL] code cell #{i}: {e}")
    if errs == 0:
        print(f"  [PASS] all {sum(c['cell_type']=='code' for c in cells)} code cells AST-parse cleanly")
    return errs


def main():
    errs = ast_check(CELLS)
    if errs:
        raise SystemExit(f"{errs} cell(s) failed AST check — fix before writing the notebook")
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
