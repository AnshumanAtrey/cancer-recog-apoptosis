# KRAS-G12D binder retry via ODesign — run guide (the honest "different generator" lever)

**Why ODesign (not more PXDesign, not BindCraft):** KRAS v2 was bounded because AF2-IG and Protenix **anti-correlate**
(r=−0.41) — no design satisfies both. A generator that optimizes for AF2 (BindCraft) would inherit that problem.
**ODesign is an all-atom interaction *world-model*, not AF2-based** — a genuinely different generative prior, with
explicit **epitope/hotspot specification** (exactly our need: force the binder onto the G12D Asp). It is the one
principled shot left at the external binder. *(Repo: The-Institute-for-AI-Molecular-Design/ODesign, Apache-2.0,
arXiv 2510.22304; verified R31b.)*

> **Honest caveat (rule 7):** this is a heavy CUDA-12.1 research stack (torch 2.3.1 + pyg + deepspeed + a HF
> checkpoint + a Google-Drive CCD file). I could **not** truly dry-run the CUDA wheels on the M2 (platform-specific),
> so the steps below are ODesign's **documented** recipe + our **validated** input — not a "tested" notebook.
> Expect possible install iteration on first run. The **container path is the maintainers' recommended, lower-risk
> option** if pip fights you.

## What's already prepped (validated on M2)
- **Input JSON:** `kras_odesign_input.json` (this dir) — target = pMHC (chain A/1-108 MHC + chain B/1-10 peptide),
  designed binder length **80** (the v2 winner), **hotspot `B/4,B/5,B/6,B/7,B/8`** centred on **B/6 = the G12D Asp**.
  Validated against the staged PDB (peptide `VVVGADGVGK`, res 6 = D). ✓
- **Target structure:** `runs/rung30_kras_g12d/staging/kras_g12d_A1101_free_mut_cropped.pdb` (gitignored — upload it
  to the run box). Rename/copy to `./data/kras_g12d_A1101_free_mut_cropped.pdb` (matches `ref_file` in the JSON).

## Path A — container (recommended, lower-risk; needs a CUDA GPU box)
```bash
git clone https://github.com/The-Institute-for-AI-Molecular-Design/ODesign.git && cd ODesign
# get the protein-binder checkpoint (flexible receptor):
bash ./ckpt/get_odesign_ckpt.sh ./ckpt          # or grab odesign_base_prot_flex.pt from the HF repo
# download components.v20240608.cif (+ .rdkit_mol.pkl) from their Google Drive -> ./data/   (one-time, all-atom CCD)
cp /path/to/kras_g12d_A1101_free_mut_cropped.pdb ./data/
cp /path/to/kras_odesign_input.json ./examples/protein_design/prot_binding_prot/kras.json
# build + run the container (Dockerfile / odesign.def both shipped):
docker run --gpus all -it --rm --shm-size=8g -v $(pwd)/ckpt:/app/ODesign/ckpt \
  -v $(pwd)/data:/app/ODesign/data -v $(pwd)/outputs:/app/ODesign/outputs <odesign-image> \
  bash inference_demo.sh
```
Edit `inference_demo.sh` first:
- `infer_model_name="odesign_base_prot_flex"`  (protein binder, flexible receptor)
- `input_json_path="./examples/protein_design/prot_binding_prot/kras.json"`
- `seeds="[42,123,777,2024,31337]"`  · `N_sample=10`   → ~50 designs (multi-seed, like the KRAS v2 sweep)
- `exp_name="kras_g12d_odesign_v1"`

## Path B — Colab (built + AST-checked): `notebooks/binder_odesign_kras_colab.ipynb`
Open that notebook in Colab (rebuild with `python scripts/_build_odesign_nb.py`; every code cell is AST-checked).
Set **Runtime → T4 GPU** first. It runs, in order: GPU-guard (pre-install, probes `nvidia-smi`) → clone +
**build a Python-3.10 venv** (`uv venv --python 3.10 /content/odv`) and `pip install -r requirements.txt -f
https://data.pyg.org/whl/torch-2.3.1+cu121.html` **into it** (ODesign's stack is frozen for Py3.10; Colab is
3.12 → rdkit==2023.3.1 / torch==2.3.1 have no 3.12 wheels; NO condacolab, rule 7). **3 corrections to ODesign's
under-specified requirements.txt are applied IN PLACE before that single install** (so one resolve enforces every
frozen pin at once — validated locally on a real py3.10/linux-cp310 resolve: numpy 1.26.3 / scipy 1.15.2 /
biopython 1.83 all stay bit-identical → torch/rdkit ABI unbroken; prody backtracks to the compatible 2.3.1):
their `biotite==1.0.1` is bumped to **1.2.0** (1.0.1 lacks `biotite.interface.rdkit.from_mol`, added ≥1.2.0), and
**prody + addict** are appended (imported by ProteinMPNN/invfold but absent from reqs). `flash_attn` is
deliberately skipped — every `use_flash`/`use_deepspeed_evo_attention`/`use_lma` config defaults false, so its
20-30min nvcc compile buys nothing. All found by a **local import-scan** of their source (rule 7). → GPU-guard
checks the venv's torch →
`ckpt/get_odesign_ckpt.sh` (Drive-cached) → **CCD wget from the GitHub Release** (Drive-cached) → inputs wget from
the repo → **SMOKE run** (1 seed × 2 samples, `num_workers=0`, ~5 min) that must pass before → the **FULL sweep**
(`/content/odv/bin/python scripts/inference.py`, seeds `[42,123,777,2024,31337]`, N_sample=10). Both go through
one `run_inference()` helper: **stderr merged into stdout, streamed live AND ring-buffered**, so any crash
traceback is always printed in the copyable cell output (never lost in Colab's red box), with a **heartbeat
thread** (output-dir file count + GPU mem every 30 s) → **persist `outputs/` to Drive**.
- **Fully automatic — no uploads, no Google-Drive file IDs.** The **CCD file** (528MB) `wget`s from this
  repo's **GitHub Release** `odesign-ccd-v20240608` (the two CCD files re-hosted as release assets, since git
  rejects >100MB; the 238GB training tar is NOT included). Checkpoints from HuggingFace; inputs from the repo.
  All heavy assets are **cached in YOUR Drive** (`MyDrive/cancer-recon/odesign_assets/`) on the first run, so
  every later run reads them from Drive. ⇒ After setting a GPU, just **Run-all**.
- **Honest (rule 7):** the CUDA-12.1 wheels could not be M2-dry-run (platform-specific), so first-run install
  iteration is possible; numpy/protobuf pins may force a Colab restart (then re-run from the cell *after*
  install — do NOT re-run install). If pip fights you, fall back to **Path A** (container) on a CUDA box.

## The make-or-break: MUT-vs-WT scoring (this is the real test, ODesign only generates)
ODesign outputs binder sequences/structures against the **MUT** pMHC. A binder that *binds* isn't the win — it must
**discriminate**. For each top ODesign binder:
1. Fold it against the **MUT** pMHC and the **WT** pMHC (`runs/rung30_kras_g12d/staging/kras_g12d_A1101_free_wt_pmhc.pdb`)
   — Protenix *and* AF2-IG (both oracles, since the bound was dual-certification).
2. **Win = high on MUT, low on WT, on BOTH oracles.** A binder that grips MUT and WT equally fails (it would attack
   the WT peptide presented on normal tissue, R32).
3. If ODesign **also** yields 0 dual-certified discriminating binders → the neoantigen-pMHC external binder is
   genuinely bounded in-silico across two independent generators → recognition load sits entirely on the **internal
   key** (R27/R33/R34–40), which is already the validated contribution.

## Honest expectation
The prior is a long shot (KRAS is a ≤5%-difficulty, single-residue-discrimination target; PXDesign got 0/80). ODesign's
*different prior* + *all-atom epitope control* is a real reason it *could* succeed where AF2-based design didn't — but
treat it as one clean shot, not a sure thing. Either outcome is informative: a dual-certified discriminating binder is
a breakthrough artifact; another 0 is a strong second-generator confirmation that the route is bounded.
