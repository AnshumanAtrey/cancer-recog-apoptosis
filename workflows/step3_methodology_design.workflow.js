export const meta = {
  name: 'step3-specificity-methodology',
  description: 'Research + design + adversarially verify the Step-3 anchor specificity-audit methodology (GTEx/HPA tumour-vs-normal therapeutic window) before building',
  phases: [
    { title: 'Research', detail: 'parallel: GTEx access, HPA access, target-selection methodology, candidate clinical ground-truth' },
    { title: 'Synthesize', detail: 'merge into a falsifiable Step-3 methodology spec' },
    { title: 'Critique', detail: 'adversarial: is it circular / fakeable / non-falsifiable?' },
  ],
}

const CTX = `PROJECT: cancer-recon-apoptosis (github.com/AnshumanAtrey/cancer-recon-apoptosis).
We are designing a HYBRID bispecific "recognition ligand": one end binds a cancer-OVER-EXPRESSED
cell-surface receptor (the ANCHOR = specificity), the other end clusters the death receptor DR5
(TNFRSF10B) to trigger apoptosis. Operationalises Shriya Rai's "recognise cancer -> self-destruct".

STEP 2 (done) gave cancer-enriched surface receptor candidates from scRNA-seq (lung/breast/colon
tumour vs normal). STEP 3 = ANCHOR SPECIFICITY/SAFETY AUDIT: for each candidate, prove (or disprove)
a real therapeutic window = HIGH on cancer AND LOW on vital healthy organs, so a binder spares
healthy tissue. The cautionary example: HER2/ERBB2 topped Step 2 but is on heart (cardiomyocytes) ->
trastuzumab cardiotoxicity. A correct method MUST independently flag this; if it says HER2 is heart-safe
it is wrong. This is real cancer research — it must be falsifiable and must not fake specificity.

CANDIDATE ANCHORS: ERBB2(HER2), ERBB3(HER3), EPHB4, CD44, DDR1, MUC1, SDC1, CD74, ITGB4, ADGRG1.
TRIGGER (fixed): TNFRSF10B (DR5).
ENVIRONMENT for the build: Python on Colab CPU (no GPU); pandas/requests/scanpy available; we value
authoritative, programmatically-downloadable data + reproducible parsing. Keep it light enough for Colab.`

const RESEARCH_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    source: { type: 'string' },
    summary: { type: 'string' },
    access_methods: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: { name: {type:'string'}, url: {type:'string'}, format: {type:'string'},
        code_snippet: {type:'string'}, notes: {type:'string'} },
      required: ['name','url','format','code_snippet','notes'] } },
    key_facts: { type: 'array', items: { type: 'string' } },
    gotchas: { type: 'array', items: { type: 'string' } },
    citations: { type: 'array', items: { type: 'string' } },
  },
  required: ['source','summary','access_methods','key_facts','gotchas','citations'],
}

phase('Research')
const [gtex, hpa, method, biology] = await parallel([
  () => agent(CTX + `

YOUR TASK (GTEx access): Find the EXACT programmatic way to get GTEx normal-tissue median gene
expression (TPM/nTPM) per tissue for a list of gene symbols. Cover BOTH: (a) the GTEx Portal API
(v2) — exact base URL, endpoint path(s) for median gene expression by tissue, query params (gene id
type — does it need Ensembl/GENCODE id?), response JSON shape; (b) the bulk "gene median TPM" GCT
file — exact current download URL (v8 or v10), gzip/gct format, how to parse with pandas (skiprows,
columns). List the ~54 GTEx tissues and which are VITAL organs (heart, brain regions, lung, liver,
kidney, etc.). Give a minimal working code snippet for at least one method. Verify URLs resolve.`,
    { label: 'gtex-access', phase: 'Research', schema: RESEARCH_SCHEMA, agentType: 'general-purpose' }),

  () => agent(CTX + `

YOUR TASK (Human Protein Atlas access): Find EXACT current download URLs + file formats from
proteinatlas.org/about/download for: (1) normal_tissue.tsv.zip — protein IHC expression in normal
tissues (columns: Gene, Gene name, Tissue, Cell type, Level [Not detected/Low/Medium/High],
Reliability); (2) rna_tissue_consensus.tsv.zip — consensus RNA nTPM per normal tissue; (3)
pathology.tsv.zip — protein IHC in CANCERS (so normal vs cancer protein from the SAME IHC method =
directly comparable). Give exact URLs, the precise column names, how the expression Level / nTPM is
encoded, and a pandas parsing snippet. Note gene-symbol vs Ensembl id usage and any size caveats.
Verify the download URLs are current and resolve.`,
    { label: 'hpa-access', phase: 'Research', schema: RESEARCH_SCHEMA, agentType: 'general-purpose' }),

  () => agent(CTX + `

YOUR TASK (target-selection METHODOLOGY): How do rigorous published CAR-T / ADC / T-cell-engager /
bispecific TARGET-SELECTION studies quantify a tumour-vs-normal therapeutic window and on-target/
off-tumour safety? Cover: (1) the tissue-specificity index TAU (give the exact formula and its 0-1
interpretation); (2) tumour:normal expression ratio definitions and what fold/threshold is treated
as "selective"; (3) the concept of a VITAL/essential normal-tissue exclusion list and why a target
high in any vital organ is disqualified regardless of tumour level; (4) RNA-vs-protein discordance
caveats (why protein/IHC matters for a surface binder); (5) the specific PITFALLS that make such an
analysis FAKE or misleading (e.g., cherry-picking tissues, ignoring protein, comparing incomparable
units, no threshold pre-registration). Cite concrete papers/frameworks (e.g. CAR-T target safety
analyses, ADC target criteria, the 'paradigm of target selection' literature). The summary field
should state what a DEFENSIBLE, non-circular safety verdict looks like.`,
    { label: 'methodology', phase: 'Research', schema: RESEARCH_SCHEMA, agentType: 'general-purpose' }),

  () => agent(CTX + `

YOUR TASK (candidate clinical GROUND-TRUTH — this is our falsification check): For EACH candidate
anchor — ERBB2/HER2, ERBB3/HER3, EPHB4, CD44, DDR1, MUC1, SDC1, CD74, ITGB4, ADGRG1 — summarise from
the literature: (a) its KNOWN normal-tissue expression (which healthy organs, especially vital ones);
(b) any KNOWN clinical on-target/off-tumour TOXICITY when targeted (e.g. HER2 -> cardiotoxicity;
EGFR -> skin/GI; CD44 -> hematopoietic/widespread; MUC1 -> epithelial). This is the ground-truth our
GTEx/HPA analysis must REPRODUCE — if our computed result contradicts well-known toxicity (e.g. calls
HER2 heart-safe), the analysis is wrong. In key_facts, give one line per candidate: "GENE: normal
tissues = ...; known tox = ...". Cite. In summary, state which candidates are widely considered to
have the BEST vs WORST tumour-selectivity per the literature.`,
    { label: 'candidate-ground-truth', phase: 'Research', schema: RESEARCH_SCHEMA, agentType: 'general-purpose' }),
])

phase('Synthesize')
const METHODOLOGY_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    what_we_prove: { type: 'string' },
    what_would_falsify: { type: 'string' },
    data_sources: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: { name:{type:'string'}, role:{type:'string'}, access_url:{type:'string'}, parse_notes:{type:'string'} },
      required: ['name','role','access_url','parse_notes'] } },
    metrics: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: { name:{type:'string'}, definition:{type:'string'}, threshold:{type:'string'} },
      required: ['name','definition','threshold'] } },
    vital_tissues: { type: 'array', items: { type: 'string' } },
    decision_rule: { type: 'string' },
    candidates: { type: 'array', items: { type: 'string' } },
    sanity_checks: { type: 'array', items: { type: 'string' } },
    pitfalls_to_avoid: { type: 'array', items: { type: 'string' } },
    script_outline: { type: 'array', items: { type: 'string' } },
  },
  required: ['what_we_prove','what_would_falsify','data_sources','metrics','vital_tissues','decision_rule','candidates','sanity_checks','pitfalls_to_avoid','script_outline'],
}
const methodology = await agent(CTX + `

YOUR TASK: synthesise the four research outputs below into ONE rigorous, FALSIFIABLE Step-3
methodology spec for scripts/07_specificity_audit.py (Python/Colab). It must: pick concrete data
sources + exact access URLs; define metrics with PRE-REGISTERED numeric thresholds (tau, tumour:normal
ratio, vital-organ rule); give a clear PASS/CAUTION/FAIL decision rule; list explicit SANITY CHECKS
that validate the method against clinical ground-truth (e.g. "HER2 must be flagged FAIL/CAUTION due to
heart; EGFR-class skin tox") — if a sanity check fails, the run is not trusted; enumerate pitfalls to
avoid; and give a concrete script_outline (ordered steps). Keep it implementable on Colab CPU.

GTEX RESEARCH: ${JSON.stringify(gtex)}
HPA RESEARCH: ${JSON.stringify(hpa)}
METHODOLOGY RESEARCH: ${JSON.stringify(method)}
CANDIDATE GROUND-TRUTH: ${JSON.stringify(biology)}`,
  { label: 'synthesis', phase: 'Synthesize', schema: METHODOLOGY_SCHEMA })

phase('Critique')
const CRITIQUE_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    is_falsifiable: { type: 'boolean' },
    is_circular: { type: 'boolean' },
    circularity_explanation: { type: 'string' },
    rigor_issues: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: { issue:{type:'string'}, severity:{type:'string'}, fix:{type:'string'} },
      required: ['issue','severity','fix'] } },
    required_changes: { type: 'array', items: { type: 'string' } },
    verdict: { type: 'string' },
  },
  required: ['is_falsifiable','is_circular','circularity_explanation','rigor_issues','required_changes','verdict'],
}
const critique = await agent(CTX + `

YOUR TASK: adversarially stress-test the Step-3 methodology below. We must PROVE a real cancer-vs-
healthy therapeutic window — not manufacture one. Attack it hard: (1) Is the verdict FALSIFIABLE —
could the data actually make a candidate FAIL, or is every candidate engineered to pass? (2) Is it
CIRCULAR — does it assume what it's trying to prove (e.g. using the same scRNA that selected the
candidate to also "validate" it)? (3) Where could it be FAKED or mislead — cherry-picked tissues,
ignoring protein, incomparable units, post-hoc thresholds, missing vital organs? (4) Does it include
a real ground-truth check (HER2->heart etc.) that would expose a broken pipeline? (5) Are the thresholds
defensible and pre-registered, not tuned to get a nice answer? List concrete required_changes the
builder MUST make. Default to skepticism; only call it sound if it genuinely is.

METHODOLOGY: ${JSON.stringify(methodology)}
CANDIDATE GROUND-TRUTH: ${JSON.stringify(biology)}`,
  { label: 'adversarial-critique', phase: 'Critique', schema: CRITIQUE_SCHEMA })

return { methodology, critique }
