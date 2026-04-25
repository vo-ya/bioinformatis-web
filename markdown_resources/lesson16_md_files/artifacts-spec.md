# Lecture 16 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 16 (ML in Genomics: Architectures, Pitfalls, Frontiers).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-16.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-16/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Architecture family colours: CNN amber, transformer cobalt, VAE violet, GNN teal, equivariant amber-dark.
- Data-set / train-test distinction: train cobalt, test amber.
- Leakage indicators: red tint on leaky elements; green tint on clean.
- Typography: Inter for UI chrome; JetBrains Mono for tensor shapes, accession-like IDs, numerical values, architecture-hyperparameter labels.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters.
- **Run / Re-simulate / Reset** — for trigger-based computations.
- **Select / Match** — for matching / checklist artifacts.
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Architecture-to-Problem Matcher → correctness of match + explanation of inductive-bias fit.
- Data Leakage Demonstrator → per-split-strategy validation F1; gap between naive and corrected splits.
- DNA-LM Tokeniser Explorer → vocabulary size, token count, effective receptive field in bp per tokenisation choice.
- Cell Foundation Model Zero-Shot Explorer → classification accuracy under leave-one-cell vs leave-one-study protocol.
- Inductive Bias Explorer → validation F1 per architecture at each training-set size.
- Self-Supervised Pretraining Simulator → downstream F1 with vs without pretraining as a function of labelled budget.
- Genomics ML Pitfall Checklist → severity-ranked list of likely issues in the student's proposed project.

### 1.5 Feasibility gate on user input (required where input is free-form)

- The Architecture-to-Problem Matcher and Pitfall Checklist accept structured choices from a dropdown; no free-form text needed.
- Numeric sliders clamp to valid ranges with inline clamp-feedback.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted; otherwise none).
- No external data files > 100 KB; all simulations use synthetic data generated in-browser.

---

## 2. Artifact #1 — Architecture-to-Problem Matcher

**File**: `artifacts/lecture-16/01-architecture-matcher.html`
**Lecture anchor**: §1.7 Common threads
**EE framing reinforced**: architecture choice = inductive-bias matching.

### Teaching purpose

Present a sequence of genomic problems; student chooses the best-matched architecture from a fixed list; feedback explains why each choice is or isn't correct.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Problem (1 of 7):                                           │
│   "Predict gene expression from 100 kb of surrounding DNA." │
├─────────────────────────────────────────────────────────────┤
│ Pick architecture:                                          │
│   [ ] CNN (pileup-style)                                    │
│   [ ] Dilated conv + transformer                            │
│   [ ] VAE with NB likelihood                                │
│   [ ] GNN (message-passing)                                 │
│   [ ] Equivariant transformer                               │
│   [ ] Transformer (base)                                    │
│   [ ] MLP                                                    │
├─────────────────────────────────────────────────────────────┤
│ Explanation panel:                                          │
│   Correct choice: "dilated conv + transformer"              │
│   Why: multi-scale structure (motif → enhancer → gene);     │
│     sparse long-range interactions; receptive-field mgmt.    │
│   Example model: Enformer, Borzoi.                          │
├─────────────────────────────────────────────────────────────┤
│ [Next problem]  [Score: 2/3 correct]                        │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Radio-button selection per problem.
- Next / Previous buttons.
- Score tracker.

### What they see

- Problem description.
- Architecture-choice list.
- Post-answer explanation panel with correct choice + reasoning.
- Cumulative score.

### Target aha

Each problem's architecture choice is determined by data structure, not architectural novelty. Student gets quick, actionable intuition for which prior fits which data.

### Technical notes

- Pure JS.
- Hardcoded list of ~7 problems: variant calling / long-range regulation / single-cell counts / molecular property / protein structure / regulatory motif discovery / rare-variant burden.
- Each problem has a correct architecture + an explanation paragraph; student sees the explanation on any choice (whether correct or not).
- Score counter updates live.

### Acceptance criteria

- [ ] Seven problems present, navigable forward/backward.
- [ ] Correct/incorrect feedback is specific.
- [ ] Explanation panel always covers: what problem structure is, why the chosen architecture fits or doesn't.
- [ ] Score counter correct.
- [ ] Opens pre-rendered on problem 1.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Data Leakage Demonstrator

**File**: `artifacts/lecture-16/02-data-leakage.html`
**Lecture anchor**: §4.1–4.6
**EE framing reinforced**: data leakage as violated independence in correlated samples; the need for blocked cross-validation.

### Teaching purpose

Simulate a genomic regression task. Show how validation F1 changes dramatically under random / chromosome / family-aware split protocols, using a synthetic correlated dataset.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulation parameters:                                      │
│   # samples:        [──●── 2000 ]                           │
│   Correlation structure: [ spatial ▾ / family / batch ]     │
│   Correlation strength:  [──●── 0.7 ]                       │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Split strategy comparison:                                  │
│ [ random 80/20 ▾ / chromosome-split / family-aware /         │
│   leave-one-study-out ]                                     │
├─────────────────────────────────────────────────────────────┤
│ Validation F1 per split:                                    │
│   Random:        0.94 (inflated)                            │
│   Chromosome:    0.71 (corrected)                           │
│   Family-aware:  0.65 (corrected)                           │
│   Leave-one-study: 0.58 (hardest, honest)                    │
├─────────────────────────────────────────────────────────────┤
│ Visualisation:                                              │
│   chromosome bar with train/test positions shown per split  │
│   correlation matrix between train and test samples         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Random split F1 = 0.94; honest split F1 = 0.65.         │  │
│ │ Gap of 29 percentage points = pure leakage.              │   │
│ │ Strengthening correlation → widens the gap.               │  │
│ │ Lowering correlation → both converge.                     │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # samples: 500 – 10,000 (log), default 2000.
- Correlation structure: spatial (chromosome) / family / batch.
- Correlation strength: 0 – 1.0, default 0.7.
- Re-simulate.

### What they see

- Per-split-strategy F1 scores.
- Chromosome bar visualisation for chosen split.
- Sample-sample correlation matrix colour-coded.
- Outcome banner.

### Target aha

Gap of 20–30 percentage points is normal on correlated data. Random splits over-fit; blocked splits are honest. The gap grows with correlation strength and shrinks toward zero as correlation weakens.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate genotypes with block-wise correlation (simulating LD blocks or family relatedness).
- Train a simple logistic regression on a binary outcome.
- Evaluate under each split strategy; report F1.

### Acceptance criteria

- [ ] Default (spatial correlation, strength=0.7) → random F1 > chromosome F1 by > 0.15.
- [ ] Switching correlation structure to "family" changes the ordering relative to chromosome.
- [ ] Correlation strength = 0 → all splits give similar F1.
- [ ] Visualisation updates with each change.
- [ ] Opens pre-rendered with default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — DNA-LM Tokeniser Explorer

**File**: `artifacts/lecture-16/03-dna-lm-tokeniser.html`
**Lecture anchor**: §5.4 Why DNA is NOT like text
**EE framing reinforced**: tokenisation choice = representation-scale trade-off.

### Teaching purpose

Compare tokenisation strategies for DNA — per-nucleotide, 3-mer, 6-mer, BPE-style — on vocabulary size, typical token count per sequence, and effective receptive field per 512-token window.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ DNA input:                                                  │
│   Preset: [ short motif ▾ / promoter region / long seq ]     │
│   Or paste your own (≤ 5 kb).                                │
├─────────────────────────────────────────────────────────────┤
│ Tokenisation choice:                                        │
│ [ per-nucleotide ▾ / 3-mer / 6-mer / BPE ]                   │
├─────────────────────────────────────────────────────────────┤
│ Visualisation:                                              │
│   Sequence as tokens (one row per token; each token a       │
│   coloured box).                                            │
├─────────────────────────────────────────────────────────────┤
│ Statistics:                                                 │
│   Vocabulary size: 4 / 64 / 4096 / 512                       │
│   Total tokens for this sequence: 5000 / 1666 / 833 / 1200   │
│   Effective receptive field at 512 tokens: 512 bp / 1.5 kb / │
│     3 kb / ≈ 1.8 kb                                           │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Per-nucleotide: finest resolution, slowest model.      │   │
│ │ 6-mer: fixed ~3 kb receptive field with 4096-word vocab │   │
│ │ BPE: data-adaptive; hits a rough middle ground.         │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- DNA preset dropdown or paste box.
- Tokenisation choice dropdown.

### What they see

- Visual tokenisation of the sequence.
- Per-tokeniser statistics table.
- Outcome banner.

### Target aha

Tokenisation is a design choice with measurable trade-offs — vocabulary, receptive field, efficiency. There's no universally best choice; per-task experimentation is needed.

### Technical notes

- Pure JS.
- Implement four tokenisers: per-nucleotide, 3-mer sliding, 6-mer sliding, simplified BPE (pretrained subword vocab from a hardcoded small table).
- Receptive field at 512 tokens = 512 × avg(tokens-per-bp)^-1.

### Acceptance criteria

- [ ] Each tokenisation choice produces correct vocabulary + token count.
- [ ] Visualisation updates on tokenisation change.
- [ ] Preset sequences render without errors.
- [ ] Opens pre-rendered on promoter preset with 6-mer.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Cell Foundation Model Zero-Shot Explorer

**File**: `artifacts/lecture-16/04-cell-foundation.html`
**Lecture anchor**: §6.5 How to read a foundation-model paper honestly
**EE framing reinforced**: zero-shot accuracy is a function of evaluation protocol; leave-one-study-out is the honest test.

### Teaching purpose

Simulate a scGPT-like foundation model's classifier. Evaluate under leave-one-cell-out vs leave-one-study-out protocols. Show how the "zero-shot accuracy" headline number shifts with the split.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Preset single-cell cohort:                                  │
│   # studies: [──●── 5 ]                                     │
│   Cells per study: [──●── 500 ]                              │
│   Cell-type count: [──●── 8 ]                                │
│   Study-specific batch effect: [──●── 0.4 ]                 │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Evaluation protocol:                                        │
│ [ leave-one-cell-out ▾ / leave-one-study-out ]               │
├─────────────────────────────────────────────────────────────┤
│ Results:                                                    │
│   LOO-cell accuracy:   92%                                  │
│   LOO-study accuracy:  71%                                  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ LOO-cell: evaluation leaks via study-specific batch     │  │
│ │   → inflated.                                           │   │
│ │ LOO-study: each test study is truly held out → honest.  │  │
│ │ The honest accuracy is usually 15–25 pp lower.          │   │
│ │ Reducing batch effect → gap shrinks; increasing → widens.│  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # studies: 2 – 20, default 5.
- Cells per study: 100 – 2000, default 500.
- Cell-type count: 2 – 20, default 8.
- Batch effect strength: 0 – 1.0, default 0.4.
- Protocol dropdown.

### What they see

- Accuracy per protocol side-by-side.
- Gap annotation.
- Outcome banner.

### Target aha

Foundation-model zero-shot numbers are a function of the evaluation protocol. Without leave-one-study-out, accuracy numbers are inflated by batch-effect memorisation.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate cells: each cell has a cell-type signature + study-specific batch offset.
- "Classifier" is a simplified nearest-neighbour in an embedding; embeddings = ideal signature + batch effect.
- LOO-cell: one cell held out, classify it. LOO-study: one study held out, classify its cells without any cells of that study in training.

### Acceptance criteria

- [ ] Default (5 studies, batch=0.4) → LOO-cell ≥ LOO-study by ≥ 15 pp.
- [ ] Batch effect → 0 collapses the gap.
- [ ] Batch effect → 1.0 widens gap to 30+ pp.
- [ ] Protocol dropdown updates results instantly.
- [ ] Opens pre-rendered with default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — Inductive Bias Explorer

**File**: `artifacts/lecture-16/05-inductive-bias.html`
**Lecture anchor**: §2.2 CNN vs MLP vs Transformer
**EE framing reinforced**: architecture choice is data-regime-dependent.

### Teaching purpose

Train three architectures (MLP, CNN, Transformer) on a synthetic pileup-like task at user-controlled training-set size. Show validation accuracy per architecture. Student experiences data-efficiency ordering directly.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Training-set size: [──●── 1000 ]  (log)                     │
│ Pileup height × width: [──●── 20 × 30 ]                     │
│ Task difficulty: [──●── medium ▾ / easy / hard ]             │
│ [Train all three]                                           │
├─────────────────────────────────────────────────────────────┤
│ Training progress bar per architecture (animated)           │
├─────────────────────────────────────────────────────────────┤
│ Validation F1:                                              │
│   MLP:                0.68                                  │
│   CNN:                0.94                                  │
│   Transformer:         0.77                                  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ At 1000 samples: CNN dominates.                         │  │
│ │ Raise to 100000 → transformer approaches CNN.           │  │
│ │ Architecture choice is data-regime-dependent; "choose   │  │
│ │ transformer always" is a large-data mirage.             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Training set size: 100 – 100,000 (log), default 1000.
- Pileup dimensions.
- Task difficulty.
- Train button.

### What they see

- Per-architecture F1 after training.
- Training progress animation.
- Outcome banner.

### Target aha

At small data, inductive bias wins. At large data, architectures converge. "Use a transformer because it's new" is poor engineering when your data is small and structured.

### Technical notes

- Pure JS, seeded PRNG.
- Training simulated via a simplified proxy: compute theoretical F1 from (architecture × data size × task difficulty) using a calibrated curve. Animate "progress" as gradient-descent-like log-loss decrease.
- No real backpropagation needed; use hardcoded learning curves.

### Acceptance criteria

- [ ] Default (1000 samples) → CNN F1 > MLP F1; CNN F1 > Transformer F1.
- [ ] Raising sample count to 100k → Transformer F1 rises toward CNN.
- [ ] Progress animation completes in ≤ 400 ms per architecture.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — Self-Supervised Pretraining Simulator

**File**: `artifacts/lecture-16/06-ssl-pretraining.html`
**Lecture anchor**: §3.3 Self-supervised pretraining
**EE framing reinforced**: pretraining as representation learning from unlabelled data.

### Teaching purpose

Show the downstream-label efficiency gain of pretraining. Train a small transformer from scratch vs pretrain + fine-tune. Plot validation F1 vs labelled-example count.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Pretraining corpus size: [──●── 1M tokens ]  (toggleable)    │
│ Downstream labelled budget: [──●── 1000 ]                   │
│ Pretraining objective: [ MLM ▾ / next-token ]               │
│ [Train from scratch]  [Pretrain + fine-tune]                │
├─────────────────────────────────────────────────────────────┤
│ Validation F1 curves:                                       │
│   (two lines, one per condition, as a function of labelled   │
│    examples)                                                 │
│   Scratch at labelled=1000:         0.62                     │
│   Pretrained + fine-tune at 1000:   0.89                     │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Pretrained + 1k labels matches scratch-trained at 100k  │  │
│ │   labels: ~100× sample efficiency.                      │   │
│ │ Remove pretraining → downstream collapses at small N.   │   │
│ │ This is why foundation models matter most in genomics: │   │
│ │   labels are scarce.                                    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Pretraining corpus size: 0 – 10M tokens (log), default 1M.
- Downstream labelled budget: 100 – 100,000 (log), default 1000.
- Pretraining objective: MLM or next-token.
- Train buttons.

### What they see

- Per-condition F1 curve as a function of labelled examples.
- Outcome banner quantifying the sample-efficiency gain.

### Target aha

Pretraining is ~10–100× sample-efficient compared to training from scratch, and the efficiency gap widens with smaller labelled budgets. In genomics, where labels are scarce, pretraining is not optional — it's the default.

### Technical notes

- Pure JS.
- Simulate training via hardcoded F1-vs-budget curves for (scratch / pretrained) × (MLM / next-token). Curves calibrated to realistic numbers.
- Visual: two line plots on a shared log-log axis.

### Acceptance criteria

- [ ] Default → pretrained F1 > scratch F1 by > 0.2 at 1000 labels.
- [ ] Increasing labelled budget to 100,000 → both converge.
- [ ] Zero pretraining corpus → pretrained = scratch.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Genomics ML Pitfall Checklist

**File**: `artifacts/lecture-16/07-pitfall-checklist.html`
**Lecture anchor**: §7.4 Convergence with mainstream AI (or wrap-up)
**EE framing reinforced**: genomics ML design bugs are recognisable; the mistakes are predictable.

### Teaching purpose

Student describes a proposed genomics ML project via structured dropdowns (architecture, data source, split strategy, evaluation metric, baseline, likelihood choice, pretraining strategy). The checklist flags likely pitfalls with severity tags and explanatory hints.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Define your proposed project:                               │
│   Architecture:         [ CNN ▾ / Transformer / VAE / ... ] │
│   Data source:          [ SRA / dbGaP / GEO / GTEx / ... ]   │
│   Split strategy:       [ random / chromosome / family /     │
│                           leave-one-study ]                 │
│   Labelled budget:      [──●── 10,000 ]                     │
│   Loss function:        [ MSE ▾ / log-loss / NB / … ]        │
│   Baseline:             [ specialised SOTA ▾ / naive ]       │
│   Pretraining strategy: [ none ▾ / DNA-LM / protein-LM ]    │
│ [Audit]                                                     │
├─────────────────────────────────────────────────────────────┤
│ Flags:                                                      │
│   🔴 Random split on homologous sequences → leakage risk    │
│   🟡 MSE loss for count data → wrong likelihood              │
│   🟡 Naive baseline only → insufficient comparison          │
│   🟢 Leave-one-study-out split chosen → good               │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 2 red flags, 2 yellow, 1 green.                          │  │
│ │ Most critical: fix the split strategy before any training│  │
│ │ Second priority: pick NB loss for single-cell counts.    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Seven dropdowns / sliders describing the proposed project.
- Audit button.

### What they see

- Severity-ranked list of flags (red / yellow / green).
- Outcome banner prioritising which issues matter most.

### Target aha

Common genomics ML bugs are recognisable from a checklist. Most of them (leakage, wrong likelihood, weak baselines) recur across projects. Habit: audit your project setup before training.

### Technical notes

- Pure JS.
- Hardcoded decision table: per-combination-of-choices → list of flags with severity + explanation.
- Examples of pitfalls covered: random split on homologous/correlated data; MSE for counts; weak baselines; mismatched architecture to data type; pretraining on test-contaminated corpora; leave-one-sample in related-samples cohort.

### Acceptance criteria

- [ ] Default project configuration flags ≥ 2 issues.
- [ ] Changing "split strategy" from random to chromosome downgrades one red to green.
- [ ] Changing "loss" from MSE to NB on single-cell counts downgrades another.
- [ ] Fully-green configuration is achievable with correct choices.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Architecture-family colours identical across artifacts.
- "Train vs test" colour convention (cobalt / amber) identical across leakage-demo and foundation-model artifacts.
- Every artifact emits an **outcome banner** per convention §1.4.
- KaTeX (where used) loaded consistently.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls beyond KaTeX (where used).
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
