# Lecture 18 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 18 (Cancer Genomics: Integrated Capstone).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-18.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-18/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Tumour / normal palette: tumour red, normal cobalt; subclones progressive tints.
- Substitution colours: C>A cobalt, C>G amber, C>T red, T>A grey, T>C teal, T>G violet (COSMIC convention).
- Tier colours: Tier I green, Tier II amber-light, Tier III grey, Tier IV muted.
- Biomarker colours: TMB cobalt, MSI amber, HRD teal.
- Typography: Inter for UI chrome; JetBrains Mono for HGVS variants, gene names, signature IDs (SBS1...), VAF values, CCF values, percentages.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters.
- **Re-simulate / Run / Reset** — for stochastic / triggered computations.
- **Click-to-select** — for variant choice in interpretation artifacts.
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Sequencing Design Chooser → recommended design + tradeoffs explanation.
- Somatic Variant Caller Simulator → TP/FP/FN per subclone at the current purity × depth; detection-limit curve.
- Mutational Signature Decomposer → recovered exposure weights vs ground truth; fit quality (cosine similarity); suggested K.
- Gene Fusion Visualiser → fusion annotation (parent genes, breakpoints, domains, clinical action).
- Subclonal Reconstruction Explorer → recovered subclones (CCF + variant assignment) vs ground truth; phylogenetic tree.
- AMP/ASCO/CAP Oncology Tier Classifier → Tier I–IV classification + matched therapy action.
- Targeted Therapy Matcher → ranked therapy plan for the variant list.

### 1.5 Feasibility gate on user input (required where input is free-form)

- All seven artifacts use structured input (preset variants, sliders, dropdowns) — no free-form text parsing needed.
- Numeric inputs clamp to valid ranges.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted; otherwise none).
- No external data files > 100 KB; all knowledge-base entries (variants, drugs, COSMIC signatures, gene lists) hardcoded in-browser.
- All artifacts explicitly state they are **educational only** and not clinical-grade tools.

---

## 2. Artifact #1 — Sequencing Design Chooser

**File**: `artifacts/lecture-18/01-sequencing-design.html`
**Lecture anchor**: §2.5 The 2025 practical workflow
**EE framing reinforced**: design follows clinical question + sample + budget constraints.

### Teaching purpose

Student describes a clinical use case; artifact recommends the sequencing design (tumour-only panel / tumour-normal panel / WES / WGS / ctDNA / combo) with a reasoned explanation.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Clinical scenario:                                          │
│   Use case: [ newly-diagnosed advanced solid tumour ▾ /     │
│               early-stage resected (MRD monitoring) /       │
│               rare-tumour research / leukaemia /            │
│               unknown primary / hereditary cancer risk ]    │
│   Sample availability: [ fresh-frozen tumour ▾ / FFPE only /│
│                          blood only / FFPE + blood ]         │
│   Budget: [ constrained ▾ / moderate / research-grade ]      │
│   Turnaround time needed: [ urgent (<7d) ▾ / standard /     │
│                              research ]                      │
├─────────────────────────────────────────────────────────────┤
│ Recommendation:                                             │
│   Use: **Tumour-normal paired panel (FoundationOne-style,   │
│         324 genes)**                                         │
│   Rationale:                                                 │
│   - Advanced solid tumour → need actionable-gene coverage    │
│   - FFPE + blood available → paired design feasible          │
│   - Constrained budget → WES/WGS too expensive                │
│   - Standard TAT (7–14 d) → panel achievable                  │
│ Alternative: tumour-only panel if blood unavailable          │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Panel covers the actionable gene list at high depth.   │   │
│ │ Switch to "MRD monitoring" → recommendation shifts to  │   │
│ │   ctDNA assay.                                          │  │
│ │ Switch to "research grade" → recommendation shifts to  │   │
│ │   WGS + ctDNA combination.                              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Four dropdowns: use case, sample availability, budget, turnaround time.

### What they see

- Recommendation card (design + rationale).
- Alternative suggestion if primary has a caveat.
- Outcome banner.

### Target aha

Design is deterministic given the inputs. Clinical question + sample + budget + TAT jointly drive the choice; the decision is auditable.

### Technical notes

- Pure JS.
- Hardcoded decision table: (use case × sample × budget × TAT) → (primary design, alternative, rationale text).
- ~24 scenarios covered with reasonable approximations.

### Acceptance criteria

- [ ] Default (advanced solid, FFPE + blood, moderate budget, standard TAT) → tumour-normal panel.
- [ ] MRD monitoring preset → ctDNA.
- [ ] Research grade + fresh-frozen → WGS + ctDNA combo.
- [ ] Rationale text updates dynamically.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Somatic Variant Caller Simulator

**File**: `artifacts/lecture-18/02-somatic-caller.html`
**Lecture anchor**: §3.3 Low-VAF subclonal detection
**EE framing reinforced**: detection-theoretic limits under purity × depth × error floor.

### Teaching purpose

Simulate a tumour-normal pair at user-controlled tumour purity, sequencing depth, subclonal structure, and sequencer error rate. Run a simplified Mutect2-like somatic caller. Show per-subclone detection rate (TP/FP/FN) as a function of the inputs.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulation parameters:                                      │
│   Tumour purity: [──●── 0.40 ]                              │
│   Tumour coverage: [──●── 150× ]  (log)                      │
│   Normal coverage: [──●── 30× ]                              │
│   Subclones: [──●── 3 ]                                     │
│   Sequencer error rate: [──●── 0.001 ]                      │
│   FFPE orientation bias: [──●── 0.0 (off) ]                  │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Per-subclone detection:                                     │
│   Clone A (CCF 100%, clonal): TP 250 / FP 5 / FN 0           │
│   Clone B (CCF 50%): TP 120 / FP 2 / FN 10                   │
│   Clone C (CCF 15%): TP 15 / FP 1 / FN 40 ← many missed       │
├─────────────────────────────────────────────────────────────┤
│ Detection-limit curve:                                      │
│   sensitivity vs CCF threshold                              │
│   current purity/depth marked                               │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Clonal + mid-CCF well-detected; low-CCF mostly missed. │   │
│ │ Double tumour coverage to 300× → Clone C TP goes from  │   │
│ │   15 to 32 (still not all).                             │   │
│ │ Add FFPE orientation bias → Clone C FP explodes.        │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Tumour purity: 0.1 – 1.0, default 0.40.
- Tumour coverage: 30 – 1000× (log), default 150×.
- Normal coverage: 10 – 200× (log), default 30×.
- Subclones: 1 – 5, default 3.
- Sequencer error rate: 0.0001 – 0.01 (log), default 0.001.
- FFPE orientation bias: 0 – 0.05, default 0.

### What they see

- Per-subclone detection counts (TP/FP/FN).
- Detection-limit curve (sensitivity vs CCF) with current-setting marker.
- Outcome banner.

### Target aha

Low-CCF subclones require high purity AND high depth; neither alone is enough. FFPE orientation bias degrades subclone detection sharply because the error floor rises. Signal-processing-style detection limits are directly visible.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate variants per subclone; simulate reads with sequencer error; apply orientation bias if enabled.
- Simplified Mutect2-like caller: per-position log-likelihood; threshold on LR; filter orientation bias.
- Detection-limit curve: analytical / semi-empirical.

### Acceptance criteria

- [ ] Default (purity 0.4, 150×) detects clonal + mid-CCF well; low-CCF partially.
- [ ] Increasing depth to 500× recovers more low-CCF.
- [ ] Adding orientation bias spikes FPs.
- [ ] Reducing purity to 0.1 makes even mid-CCF hard.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — Mutational Signature Decomposer

**File**: `artifacts/lecture-18/03-signature-decomposer.html`
**Lecture anchor**: §3.4 Mutational signatures
**EE framing reinforced**: NMF decomposition of a 96-dim spectrum into interpretable signatures.

### Teaching purpose

Generate a synthetic 96-bin mutation spectrum as a mixture of user-selected COSMIC signatures at user-specified exposures. Run a simplified NMF decomposition with a user-chosen K. Show recovered exposures vs ground truth.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Ground-truth mixture:                                       │
│   Signature 1: [ SBS1 ▾ ] at exposure [──●── 0.40 ]          │
│   Signature 2: [ SBS7a ▾ ] at exposure [──●── 0.35 ]         │
│   Signature 3: [ SBS3 ▾ ] at exposure [──●── 0.25 ]          │
│   Total mutations simulated: [──●── 10,000 ]                 │
│ [Simulate spectrum]                                          │
├─────────────────────────────────────────────────────────────┤
│ Fitted NMF:                                                 │
│   K (number of signatures): [──●── 3 ]                      │
│ [Run NMF]                                                   │
├─────────────────────────────────────────────────────────────┤
│ Comparison:                                                 │
│   Truth 96-bin spectrum vs fitted spectrum                  │
│   Recovered exposures: 0.38 / 0.36 / 0.26 (cosine sim 0.97) │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ K = 3 recovers the planted signatures well.             │   │
│ │ Try K = 4 → extra signature picks up noise; cosine     │   │
│ │   similarity barely improves (overfitting).             │   │
│ │ Try K = 2 → merges two signatures; cosine drops.        │   │
│ │ Simulate fewer mutations (1000) → noise dominates.      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- 3 signature dropdowns + exposure sliders.
- Mutation-count slider.
- K slider (# signatures to fit).
- Simulate / Run buttons.

### What they see

- Ground-truth vs fitted 96-bin bar charts.
- Recovered exposure values.
- Cosine similarity fit quality.
- Outcome banner.

### Target aha

NMF decomposition of the mutation spectrum is straightforward; the hard part is picking K. Too few signatures → merge; too many → overfit noise. Canonical model-selection problem from L16.

### Technical notes

- Pure JS.
- COSMIC signatures SBS1, SBS2, SBS3, SBS4, SBS5, SBS6, SBS7a, SBS13, SBS40 hardcoded (96-dim vectors).
- Mixture = weighted sum; simulate Poisson per bin at count `total × weight × component_bin_probability`.
- NMF: multiplicative-update algorithm, ~50 iterations.
- Cosine similarity between truth and recovered signatures.

### Acceptance criteria

- [ ] Default (3-signature mixture, 10k muts, K=3) → cosine sim ≥ 0.95.
- [ ] K=2 on a 3-signature truth → cosine drops noticeably.
- [ ] K=5 on a 3-signature truth → extra signatures fit noise; fit improves marginally.
- [ ] 1000 mutations → noisy; cosine drops.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Gene Fusion Visualiser

**File**: `artifacts/lecture-18/04-fusion-visualiser.html`
**Lecture anchor**: §4.3 Gene fusions from RNA-seq
**EE framing reinforced**: fusion = modular recombination; kinase domain + oligomerisation domain → constitutive activation.

### Teaching purpose

Pick a preset canonical fusion; visualise parent genes, the breakpoint, the resulting fusion transcript, the fusion protein domain structure, and the clinical relevance.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Fusion preset: [ BCR-ABL1 ▾ / EML4-ALK / NTRK3-ETV6 /       │
│                  TMPRSS2-ERG / FGFR3-TACC3 ]                 │
├─────────────────────────────────────────────────────────────┤
│ Parent genes:                                               │
│   EML4 (chr 2p21): 23 exons, TAPE + coiled-coil domains     │
│   ALK (chr 2p23): 29 exons, kinase + extracellular domain    │
│                                                             │
│ Rearrangement event:                                        │
│   chr 2p inversion joining EML4 exon 13 → ALK exon 20       │
│                                                             │
│ Fusion transcript:                                          │
│   [EML4 exons 1–13] [ALK exons 20–29]                        │
│                                                             │
│ Fusion protein:                                             │
│   EML4 coiled-coil → dimerisation → ALK kinase constitutive │
│                                                             │
│ Clinical action:                                            │
│   EML4-ALK in ~5% of lung adenocarcinomas                    │
│   Tier I: crizotinib, alectinib, lorlatinib (FDA-approved)  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ EML4-ALK activation mechanism: coiled-coil forces      │  │
│ │   dimerisation → ALK kinase phosphorylates itself.     │   │
│ │ Switch to BCR-ABL1 → same story: BCR's coiled-coil +   │  │
│ │   ABL1 kinase. Different genes, same mechanism pattern. │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Fusion preset dropdown (5 presets).

### What they see

- Parent-gene summary.
- Rearrangement event description.
- Fusion transcript + protein structure diagram.
- Clinical-action panel.
- Outcome banner.

### Target aha

Fusions are modular. The same "oligomerisation domain + kinase domain" pattern appears in BCR-ABL1, EML4-ALK, TMPRSS2-ERG. The fusion activates a kinase by forcing dimerisation; targeted therapy inhibits the activated kinase. A canonical precision-oncology story repeats across oncogenes.

### Technical notes

- Pure JS.
- Hardcoded metadata per fusion: parent genes, domains, breakpoint, clinical significance, approved drugs.
- SVG visualisation: parent-gene bars → junction arrow → fusion bar.

### Acceptance criteria

- [ ] Each preset renders a meaningful visualisation.
- [ ] Clinical-action text per preset is correct.
- [ ] Mechanism pattern (dimerisation + kinase) highlighted when applicable.
- [ ] Opens pre-rendered with EML4-ALK default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — Subclonal Reconstruction Explorer

**File**: `artifacts/lecture-18/05-subclonal-reconstruction.html`
**Lecture anchor**: §5.1–5.2 Subclonal architecture and phylogenetic reconstruction
**EE framing reinforced**: blind source separation / deconvolution from bulk sequencing.

### Teaching purpose

Simulate a tumour with user-specified subclonal structure (purity, number of subclones, branching tree). Simulate bulk sequencing; compute per-variant VAF. Run a PyClone-style clustering to recover subclones. Compare recovered vs ground-truth.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Tumour structure:                                           │
│   Tumour purity: [──●── 0.50 ]                              │
│   Number of subclones: [──●── 3 ]                           │
│   Branching topology: [ linear ▾ / branching ]               │
│   Variants per subclone: [──●── 100 ]                        │
│   Tumour coverage: [──●── 100× ]                             │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ VAF distribution (simulated):                               │
│   histogram showing modes per subclone                      │
├─────────────────────────────────────────────────────────────┤
│ PyClone-style recovery (CCF space):                         │
│   recovered clusters with estimated CCF                     │
│   recovered vs true assignment accuracy                     │
├─────────────────────────────────────────────────────────────┤
│ Phylogenetic tree:                                          │
│   trunk + branches + leaves with driver annotations         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 3 subclones recovered at CCF 100% / 48% / 22%           │  │
│ │   (true: 100% / 50% / 20%; well-matched).               │  │
│ │ Lower purity to 0.20 → CCFs harder to estimate; clusters│  │
│ │   may merge.                                            │   │
│ │ Increase subclones to 5 → recovery starts breaking.     │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Tumour purity: 0.1 – 1.0, default 0.50.
- Subclones: 1 – 5, default 3.
- Topology: linear / branching, default branching.
- Variants per subclone: 50 – 500, default 100.
- Coverage: 30 – 500×, default 100×.

### What they see

- VAF histogram.
- CCF-space clustering visualisation.
- Phylogenetic tree.
- Outcome banner.

### Target aha

Subclonal reconstruction is source separation: bulk mixture → per-source CCF + variant assignment. Purity × coverage × subclone-count together determine recoverability. At some level of complexity, reconstruction breaks down.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate variants per subclone, compute VAF from (purity × CCF × 0.5) + noise.
- PyClone-like clustering via Dirichlet Process mixture or simple EM.
- Phylogenetic tree inferred from cluster CCFs + assignment (nested-subset logic).

### Acceptance criteria

- [ ] Default (purity 0.5, 3 subclones, 100×) recovers 3 clusters within 5% CCF.
- [ ] Purity 0.2 → recovery degrades visibly.
- [ ] 5 subclones → clusters start merging.
- [ ] Tree topology matches for simple cases (≤ 3 subclones).
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — AMP/ASCO/CAP Oncology Tier Classifier

**File**: `artifacts/lecture-18/06-oncology-tier.html`
**Lecture anchor**: §6.1 AMP/ASCO/CAP tier system
**EE framing reinforced**: action-centric classification; different question from germline ACMG/AMP.

### Teaching purpose

Input: variant (gene + HGVS + cancer type) + evidence level. Artifact applies AMP/ASCO/CAP rules → returns Tier I–IV + recommended therapy action + citations.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Variant characterisation:                                   │
│   Gene: [ EGFR ▾ ]                                           │
│   Variant: [ L858R ▾ / T790M / exon19del / L861Q / ... ]     │
│   Cancer type: [ non-small cell lung cancer ▾ / colorectal / │
│                  melanoma / breast / ... ]                   │
│   Evidence level: [ A (FDA-approved) ▾ / B (guidelines) /    │
│                     C (different-cancer-type) / D (trial) ]  │
├─────────────────────────────────────────────────────────────┤
│ Classification:                                             │
│   Tier I-A — FDA-approved targeted therapy in this cancer    │
│             type for this variant.                           │
│   Drug: erlotinib, gefitinib, osimertinib                   │
│   Citation: OncoKB Level 1; FDA label, 2013 onwards          │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ EGFR L858R in NSCLC is the canonical Tier I-A example. │   │
│ │ Change cancer type to glioblastoma → Tier II-C          │  │
│ │   (drug FDA-approved for different cancer).             │   │
│ │ Change variant to a rare missense → Tier III unless     │   │
│ │   functional evidence present.                          │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Gene dropdown.
- Variant dropdown (populated per gene with common variants).
- Cancer-type dropdown.
- Evidence-level dropdown.

### What they see

- Tier classification + level.
- Recommended therapy + citation.
- Outcome banner.

### Target aha

Oncology variant tiering is a deterministic lookup: variant × cancer type × evidence-level → tier. Contrast with ACMG/AMP germline (evidence-aggregation for pathogenicity).

### Technical notes

- Pure JS.
- Hardcoded knowledge base: (gene, variant, cancer type) → (tier, drug(s), evidence level) for ~30 canonical driver-variant-cancer combinations.
- Rule application: exact match → Tier I-A; same variant different cancer → Tier II-C; unknown match → Tier III.

### Acceptance criteria

- [ ] EGFR L858R in NSCLC → Tier I-A; drugs listed.
- [ ] EGFR L858R in glioblastoma → Tier II-C (cross-tumour).
- [ ] Rare PIK3CA missense → Tier III.
- [ ] Common gnomAD SNP → Tier IV.
- [ ] Opens pre-rendered with EGFR L858R / NSCLC default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Targeted Therapy Matcher

**File**: `artifacts/lecture-18/07-therapy-matcher.html`
**Lecture anchor**: §7.1 Targeted therapy decisions
**EE framing reinforced**: therapy selection as a multi-variant integration problem.

### Teaching purpose

Input: a list of somatic variants + cancer type (from a preset "panel run"). Artifact maps each variant to potential therapies; integrates biomarkers (TMB, MSI, HRD) if present; returns a ranked treatment plan.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Case summary (preset):                                      │
│   Cancer: [ lung adenocarcinoma ▾ / colorectal MSI-H /      │
│             BRCA+ ovarian / rare tumour unknown primary ]    │
│   Variants detected (from preset panel run):                │
│     - EGFR L858R (Tier I-A)                                 │
│     - TP53 R175H (Tier III, prognostic)                      │
│     - PTEN loss (Tier II-D)                                  │
│   Biomarkers:                                               │
│     - TMB: 4.8 mut/Mb (low)                                 │
│     - MSI: stable                                           │
│     - HRD: negative                                          │
├─────────────────────────────────────────────────────────────┤
│ Ranked treatment plan:                                      │
│   1. **First line**: osimertinib (FDA-approved EGFR TKI for │
│      NSCLC, Tier I-A). Standard of care.                     │
│   2. Alternative: erlotinib / gefitinib (older-generation    │
│      TKIs).                                                  │
│   3. If progression on first line: check for resistance      │
│      mutation (T790M); if present → osimertinib continues;   │
│      if not → other lines (chemo).                           │
│   4. PTEN loss → consider PI3K pathway trial.                │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Tier I-A dominant → targeted therapy first line.        │  │
│ │ Biomarkers (TMB, MSI, HRD) not informative here.        │   │
│ │ Switch preset to "MSI-H colorectal" → immunotherapy     │  │
│ │   first-line (pembrolizumab).                           │  │
│ │ Switch to "BRCA+ ovarian" → PARP inhibitor first-line.  │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Preset dropdown (4 canonical cases).

### What they see

- Case summary (variants + biomarkers).
- Ranked treatment plan.
- Outcome banner.

### Target aha

Therapy matching integrates multiple signals (tier variants + biomarkers + cancer type). The mapping is deterministic given the inputs; the human judgment is in ordering therapies + considering patient-specific factors.

### Technical notes

- Pure JS.
- Hardcoded 4 presets with realistic variant + biomarker profiles.
- Therapy-matching logic: variant tier → approved drugs; biomarkers → tumour-agnostic indications; integration per priority rules.

### Acceptance criteria

- [ ] EGFR+ lung preset → osimertinib first line.
- [ ] MSI-H colorectal preset → pembrolizumab first line.
- [ ] BRCA+ ovarian preset → olaparib / niraparib first line.
- [ ] Rare tumour preset → clinical-trial focus.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Substitution-type colours identical across Artifact #3 and lecture Figure 5.
- Tier colours identical across Artifacts #6, #7 and lecture Figure 9.
- Biomarker colours identical across Artifact #7 and lecture Figure 10.
- Every artifact emits an **outcome banner** per convention §1.4.
- Educational-disclaimer text appears in every artifact.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
- [ ] Educational disclaimer visible somewhere in the artifact.
