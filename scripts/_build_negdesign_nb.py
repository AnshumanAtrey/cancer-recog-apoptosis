#!/usr/bin/env python3
"""Builder for notebooks/negdesign_two_state_colab.ipynb (RUNG-26f).
Emits a valid .ipynb and AST-checks every code cell. Run: python3 scripts/_build_negdesign_nb.py
Keeping the builder (not hand-edited JSON) so the notebook stays regenerable + lint-clean."""
import ast, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "notebooks", "negdesign_two_state_colab.ipynb")

CELLS = []
def md(text):  CELLS.append(("markdown", text))
def code(text): CELLS.append(("code", text))

md(r"""# RUNG-26f — ProteinMPNN two-state **negative design** + AF2 confirmation

The principled fix RUNG-26d's autopsy proved we need. PXDesign/RFdiffusion are **positive** design:
a hotspot makes the binder *contact* the mutation but cannot make binding *depend* on it, so binders
to IDH1-R132H bind MUT≈WT (NULL across AF2-IG, ColabDesign-AF2, Protenix). **Negative design** scores
the binder against the MUT pMHC **and against the WT pMHC** and keeps only sequences the model finds
much *less* likely on WT — discrimination by construction.

**Mechanic** (same backbone for both states — that's the point): the peptide's mutated position is
`HIS` in the MUT target / `ARG` in the WT target (IDH1 R132H). ProteinMPNN scores the binder chain's
conditional NLL given the fixed target identities, so:
- `dscore = NLL_wt − NLL_mut` &nbsp; (**positive ⇒ binder prefers MUT ⇒ discriminating**).

This notebook **reuses the M2-validated harness** (`scripts/57_negdesign_proteinmpnn.py`) — already gated
on the known-non-specific rank_1/2/3 (|dscore|≈0, within noise → the scorer doesn't hallucinate
discrimination). Here it (1) **generates** new binders two-state-selected by dscore on every available
backbone (GPU-fast sampling), then (2) **AF2-confirms** the top candidates by folding each vs MUT and WT
pMHC (the same `score()` used in all prior RUNG-26 runs).

**Bars:** MPNN dscore is a *cheap screen* (noise floor ~0.03 unpaired). The real gate is AF2:
`mutant_specific = binder on MUT (pae_interaction≤10, binder_plddt≥80) AND discrimination (pae_wt − pae_mut) ≥ 3`.

**Repoint to BRAF** (V600E, V→E — the chemically tractable case) by editing the `CFG` block in Cell 2
once the BRAF MUT/WT pMHC folds land. T4-OK.
""")

code(r"""#@title Cell 1 — clone repo + ProteinMPNN + GPU-guard + ColabDesign/AF2 params  [~6 min]
import os, glob, subprocess
# torch (ProteinMPNN, Cell 2) + jax (AF2, Cell 3) must SHARE one T4 -> stop jax pre-grabbing
# all VRAM, and reduce torch fragmentation. MUST be set before jax/torch import CUDA.
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')
os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')   # jax allocates on demand
os.environ.setdefault('XLA_PYTHON_CLIENT_MEM_FRACTION', '.8')
from pathlib import Path
REPO = Path('/content/cancer-recog-apoptosis')
if REPO.exists():
    get_ipython().system('cd {REPO} && git fetch origin && git reset --hard origin/main')
else:
    get_ipython().system('git clone https://github.com/AnshumanAtrey/cancer-recog-apoptosis.git {REPO}')
os.chdir(REPO)
# ProteinMPNN must live where scripts/57 expects it: <repo>/.tools/ProteinMPNN
MPNN = REPO / '.tools' / 'ProteinMPNN'
if not MPNN.exists():
    get_ipython().system('git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git {MPNN}')
get_ipython().system('nvidia-smi -L')
from google.colab import drive; drive.mount('/content/drive')
def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode: print('  !', (r.stderr or '')[-300:])
    return r.returncode
sh('pip install -q gemmi git+https://github.com/sokrypton/ColabDesign.git@v1.1.1')
if not glob.glob('params/params_model_*.npz'):
    sh('apt-get -qq install -y aria2 >/dev/null; mkdir -p params && cd params && '
       '(aria2c -q -x16 https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar || '
       'wget -q https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar) && '
       'tar -xf alphafold_params_2022-12-06.tar && rm -f alphafold_params_2022-12-06.tar')
import jax, torch
assert any(d.platform == 'gpu' for d in jax.devices()), 'NO GPU — Runtime>Change runtime type>T4, rerun'
print('jax GPU:', jax.devices(), '| torch cuda:', torch.cuda.is_available())
from colabdesign import mk_afdesign_model  # noqa: F401
print('[CELL 1] done')
""")

code(r"""#@title Cell 2 — ProteinMPNN two-state GENERATION across backbones -> rank by dscore
import os, glob, json, importlib.util
import numpy as np

# ====================== CONFIG — edit here to repoint the target ======================
CFG = dict(
    tag='rung26f_negdesign_idh1',
    pep_mut='IIGHHAYGDQY', pep_wt='IIGRHAYGDQY', p=4,           # IDH1 R132H (His<-Arg) at p4
    mut_resname='HIS', wt_resname='ARG', groove_prefix='SHSMRYF',
    meta='/content/drive/MyDrive/cancer-recon/rung26b_rfdiff/meta.json',  # mut_pdb/wt_pdb/hotspot
    drive_backbones='/content/drive/MyDrive/cancer-recon/rung26b_rfdiff',  # RFdiffusion PDBs (weak; opt-in)
    n=64, temp=0.25, orders=8, top_k=8, chunk=8,                # chunk bounds GPU mem (T4-safe)
    max_drive=0,    # 0 = only the 3 good repo PXDesign backbones. The 155 RFdiffusion d_*.pdb never
                    # scored as binders (pae 25-27) -> AF2-confirming their candidates wastes GPU.
                    # Raise (e.g. 12) only on a bigger GPU / more compute units.
)
# --- BRAF V600E (V->E at p3) — uncomment when the BRAF MUT/WT folds land ---
# CFG.update(tag='rung26f_negdesign_braf', pep_mut='ATEKSRWSGSH', pep_wt='ATVKSRWSGSH', p=3,
#            mut_resname='GLU', wt_resname='VAL',
#            meta='/content/drive/MyDrive/cancer-recon/rung26e_braf/meta.json',
#            drive_backbones='/content/drive/MyDrive/cancer-recon/rung26e_braf')
# ======================================================================================

REPO = '/content/cancer-recog-apoptosis'
spec = importlib.util.spec_from_file_location('nd57', f'{REPO}/scripts/57_negdesign_proteinmpnn.py')
nd = importlib.util.module_from_spec(spec); spec.loader.exec_module(nd)
# override the target so the validated harness scores OUR neoantigen
nd.PEP_MUT, nd.PEP_WT, nd.P_MUT = CFG['pep_mut'], CFG['pep_wt'], CFG['p']
nd.MUT_RESNAME, nd.WT_RESNAME, nd.GROOVE_PREFIX = CFG['mut_resname'], CFG['wt_resname'], CFG['groove_prefix']
_m, _dev = nd.load_model(); print('ProteinMPNN on', _dev)

# discover backbones: repo PXDesign complexes (always present) + Drive RFdiffusion PDBs (if any)
import torch
cands_cif = sorted(glob.glob(f'{REPO}/runs/rung26d_pxdesign_hsB4_idh1/top_designs/*.cif'))
drive_all = sorted(glob.glob(f"{CFG['drive_backbones']}/**/*.pdb", recursive=True)) if os.path.isdir(CFG['drive_backbones']) else []
drive_pdb = drive_all[:CFG['max_drive']]
if len(drive_all) > len(drive_pdb):
    print(f"NOTE: {len(drive_all)} drive backbones found, using {len(drive_pdb)} (max_drive={CFG['max_drive']}); "
          f"{len(drive_all) - len(drive_pdb)} weak RFdiffusion backbones DROPPED (raise CFG['max_drive'] to include).")
backbones, skipped = [], []
for bb in cands_cif + drive_pdb:
    try:
        nd.read_complex(bb); backbones.append(bb)
    except Exception as e:
        skipped.append((os.path.basename(bb), str(e)[:60]))
print(f'usable backbones: {len(backbones)}  (repo cif {len(cands_cif)}, drive used {len(drive_pdb)}, skipped {len(skipped)})')
for s in skipped[:5]: print('  skip', s)

OUTDIR = f'{REPO}/runs/rung26f_negdesign'; os.makedirs(OUTDIR, exist_ok=True)
all_cands = []
for i, bb in enumerate(backbones):
    g = nd.generate_two_state(bb, n_candidates=CFG['n'], temperature=CFG['temp'],
                              n_orders=CFG['orders'], top_k=CFG['top_k'], chunk=CFG['chunk'])
    for c in g['top']:
        c['backbone'] = os.path.basename(bb)
    all_cands += g['top']
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    print(f"[{i+1}/{len(backbones)}] {os.path.basename(bb):28s} dscore max={g['dscore_max']:+.4f} mean={g['dscore_mean']:+.4f}")
    # incremental save (resumable / crash-safe)
    all_cands.sort(key=lambda c: -c['dscore_wt_minus_mut'])
    json.dump({'tag': CFG['tag'], 'config': CFG, 'n_backbones': len(backbones),
               'candidates_ranked_by_dscore': all_cands},
              open(f'{OUTDIR}/candidates.json', 'w'), indent=2)
print(f"\\nTOP negative-design candidates (by MPNN dscore):")
for c in all_cands[:CFG['top_k']]:
    print(f"  dscore={c['dscore_wt_minus_mut']:+.4f}  {c['backbone']}  {c['seq']}")
print(f"[saved] {OUTDIR}/candidates.json")

# release ProteinMPNN (torch) GPU memory so Cell 3's jax/AF2 has room on the shared T4
import gc, torch
nd.load_model.__globals__['_MODEL']['m'] = None
gc.collect(); torch.cuda.empty_cache()
print('freed ProteinMPNN GPU memory; ready for Cell 3 (AF2)')
""")

code(r"""#@title Cell 3 — AF2 two-state CONFIRMATION of the top candidates (the real gate)
import os, json, gc, sys
# defensively free any ProteinMPNN (torch) GPU memory before jax/AF2 grabs the T4
try:
    import torch
    _nd = sys.modules.get('nd57')
    if _nd is not None:
        _nd.load_model.__globals__['_MODEL']['m'] = None
    gc.collect(); torch.cuda.empty_cache()
except Exception as _e:
    print('torch teardown skipped:', _e)
from colabdesign import mk_afdesign_model, clear_mem
import jax
assert any(d.platform == 'gpu' for d in jax.devices()), 'NO GPU — switch to T4 and rerun'
REPO = '/content/cancer-recog-apoptosis'; OUTDIR = f'{REPO}/runs/rung26f_negdesign'
C = json.load(open(f'{OUTDIR}/candidates.json')); CFG = C['config']
M = json.load(open(CFG['meta'])); MUT_PDB, WT_PDB, HOT = M['mut_pdb'], M['wt_pdb'], f"B{M['hotspot']}"
print('MUT_PDB', os.path.exists(MUT_PDB), '| WT_PDB', os.path.exists(WT_PDB), '| hotspot', HOT)

def score(pdb, seq):
    clear_mem()
    m = mk_afdesign_model(protocol='binder', use_multimer=False, num_recycles=3, data_dir='.')
    m.prep_inputs(pdb_filename=pdb, chain='A,B', binder_len=len(seq), hotspot=HOT,
                  rm_target_seq=False, ignore_missing=True)
    m.predict(seq=seq, verbose=False); log = m.aux['log']
    return {'pae_interaction': round(float(log['i_pae']) * 31.0, 3),
            'binder_plddt': round(float(log['plddt']) * 100.0, 1)}

TOPN = 12  # AF2-confirm the best dscore candidates
cand = C['candidates_ranked_by_dscore'][:TOPN]
print(f'AF2 two-state on {len(cand)} top-dscore negative-design candidates (bar: pae_mut<=10, plddt>=80, disc>=3)')
results = []
for k, c in enumerate(cand):
    seq = c['seq']
    mut = score(MUT_PDB, seq); wt = score(WT_PDB, seq)
    disc = round(wt['pae_interaction'] - mut['pae_interaction'], 3)
    is_binder = mut['pae_interaction'] <= 10.0 and mut['binder_plddt'] >= 80.0
    spec = bool(is_binder and disc >= 3.0)
    r = {'rank': k + 1, 'backbone': c['backbone'], 'mpnn_dscore': c['dscore_wt_minus_mut'],
         'seq': seq, 'mut': mut, 'wt': wt, 'af2_discrimination_wt_minus_mut': disc,
         'is_binder': is_binder, 'mutant_specific': spec}
    results.append(r)
    print(f"  #{k+1} dscore={c['dscore_wt_minus_mut']:+.4f} | MUT pae={mut['pae_interaction']} "
          f"plddt={mut['binder_plddt']} | WT pae={wt['pae_interaction']} | disc={disc} "
          f"| {'SPECIFIC' if spec else ('binds' if is_binder else 'no-bind')}")
    json.dump({'tag': CFG['tag'] + '_af2_confirm', 'meta': {'mut_pdb': MUT_PDB, 'wt_pdb': WT_PDB, 'hotspot': HOT},
               'bars': {'pae_bind': 10.0, 'plddt_min': 80.0, 'discrimination_min': 3.0},
               'n_specific': sum(x['mutant_specific'] for x in results), 'results': results},
              open(f'{OUTDIR}/af2_confirm.json', 'w'), indent=2)
n_spec = sum(x['mutant_specific'] for x in results)
print(f"\\n=== {n_spec}/{len(results)} AF2-confirmed mutant-specific ===")
print(f"[saved] {OUTDIR}/af2_confirm.json  — git add/commit from the repo to persist")
""")

md(r"""## Reading the result (honest bars)

- **MPNN dscore** (Cell 2) is a cheap screen. On IDH1 the M2 demo showed the two-state tail grows weakly
  with oversampling (best-of-8 ≈0 → best-of-24 ≈ +0.024) but stays near the ~0.03 noise floor — His↔Arg
  is a *hard* discrimination. A positive dscore is **necessary, not sufficient**.
- **AF2 confirmation** (Cell 3) is the gate that matters: `mutant_specific` requires the binder to actually
  bind MUT (pae≤10, plddt≥80) **and** lose ≥3 pae units on WT. This is the same metric that returned NULL
  for every positive-design binder.
- **Outcomes:**
  - `n_specific ≥ 1` on IDH1 → negative design cracked even the hard case (big result — verify the winner
    with the Protenix webserver MUT/WT, seed-matched, before any claim).
  - `n_specific = 0` on IDH1 → expected; the method is built + GPU-validated and **repoints to BRAF V600E**
    (V→E, where an unsatisfied buried charge on WT-Val gives negative design a real handle). Edit `CFG` in Cell 2.
- **Wet-lab residual (unchanged):** in-silico likelihood/structure ≠ measured affinity or specificity;
  proteome off-target, expression, and SPR/cellular validation remain.
""")

nb = {
    "cells": [
        ({"cell_type": "markdown", "metadata": {}, "source": (t if t.endswith("\n") else t + "\n").splitlines(keepends=True)}
         if k == "markdown" else
         {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
          "source": (t if t.endswith("\n") else t + "\n").splitlines(keepends=True)})
        for (k, t) in CELLS
    ],
    "metadata": {"accelerator": "GPU", "colab": {"provenance": []},
                 "kernelspec": {"name": "python3", "display_name": "Python 3"},
                 "language_info": {"name": "python"}},
    "nbformat": 4, "nbformat_minor": 0,
}

# AST-gate every code cell (strip IPython magics / shell lines first)
def astable(src):
    out = []
    for ln in src.splitlines():
        s = ln.strip()
        if s.startswith('#@title') or s.startswith('!') or s.startswith('%'):
            out.append('pass')
        else:
            out.append(ln)
    return "\n".join(out)

for i, (k, t) in enumerate(CELLS):
    if k == "code":
        ast.parse(astable(t))
print("[AST] all code cells parse clean")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(nb, open(OUT, "w"), indent=1)
print(f"[written] {OUT}  ({len(CELLS)} cells)")
