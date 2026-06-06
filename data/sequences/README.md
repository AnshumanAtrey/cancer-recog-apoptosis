# Reference Sequences for Smoke Tests

These FASTA files hold the canonical sequences used by `scripts/01_*_smoketest.py`.

| File | What it is | Role |
|---|---|---|
| `dr5_ecd_human.fasta` | DR5 (TNFRSF10B) extracellular domain | Receptor (chain A) — the cancer-cell target |
| `trail_ecd_human.fasta` | Full TRAIL (TNFSF10) ectodomain, residues 114–281 | **Positive binder (chain B)** — real DR5 binder, PDB-validated interface |
| `scrambled_ecd_control.fasta` | Composition-matched shuffle of the TRAIL ectodomain (seed 20260529) | **Negative binder (chain B)** — non-binder; should score LOW |
| `dr4_ecd_human.fasta` | DR4 (TNFRSF10A) extracellular domain | "Normal homolog" stand-in for Step-3 specificity checks |
| `trail_dr5_binding_loop.fasta` | 20-mer TRAIL fragment (legacy) | Retained for reference; NOT used by the current smoke test (see note) |
| `scrambled_control.fasta` | Shuffle of the 20-mer (legacy) | Retained for reference; NOT used by the current smoke test |

## Smoke test logic (ipTM, not affinity)

The smoke test (`scripts/01_boltz_smoketest.py`) predicts two protein–protein complexes
with Boltz-2 and compares their **interface confidence ipTM** (range [0,1], higher =
more confident interface — the same metric AlphaFold-Multimer uses):

- `DR5_ECD + TRAIL_ECD`   (positive, real binder)
- `DR5_ECD + Scrambled`   (negative, non-binder)

**PASS if `ipTM(positive) − ipTM(negative) ≥ 0.15`** → oracle discriminates → Step 2.
Otherwise → see [ASSESSMENT.md](../../docs/ASSESSMENT.md) Day-1 kill criteria.

### Why not affinity / not the 20-mer loop?

- **Affinity head is small-molecule only.** Boltz-2 rejects a protein binder for affinity
  ("Chain B is not a ligand! Affinity is currently only supported for ligands."). Our binder
  is a protein, so ipTM (structural interface confidence) is the correct proxy.
- **Full ectodomain, not a 20-mer.** A bare linear peptide lacks its native scaffold and
  usually won't dock confidently, giving an uninformative test. The full TRAIL ectodomain is
  a known DR5-binding fold present in Boltz's training data (PDB 1D4V/1DU3/1D0G).

## Sequence sources

- **DR5 ECD:** UniProt O14763 (TNFRSF10B), residues 56–183 (cysteine-rich ligand-binding region).
- **DR4 ECD:** UniProt O00220 (TNFRSF10A), analogous region.
- **TRAIL ECD:** UniProt P50591 (TNFSF10), residues 114–281 — the soluble ectodomain used in
  the DR5-complex crystal structures. Verified by exact slice of the canonical 281-aa sequence.
- **Scrambled ECD:** deterministic shuffle (`random.Random(20260529)`) of the TRAIL ectodomain,
  identical amino-acid composition.

## Biology caveat (documented in the script)

Native TRAIL binds DR5 in the groove between two trimer protomers. The smoke test models the
pairwise interface with a single TRAIL chain, which AF-class models typically recover from
training data. If the positive ipTM comes back weak (< 0.40), the pivot is to model the TRAIL
**homotrimer** (3 chains) so the true groove is present — the script prints this guidance.
