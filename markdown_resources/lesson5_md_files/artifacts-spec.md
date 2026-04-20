# Lecture 5 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 5 (Bulk RNA-seq).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-05.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-05/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified in the per-artifact section.
- The file must render standalone when opened directly in a browser.
- Artifact is embedded in the lecture page via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` near the end of `<body>`.** Missing this silently forces the iframe to a default CSS height with an internal scrollbar — caught in Lecture 4 and now checked in C6.

### 1.2 Visual design

- Use the design tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- DNA / RNA bases, when rendered, use the `--base-*` palette exactly (with `--base-u` when the context is mRNA/RNA).
- Typography: **Inter** for UI chrome; **JetBrains Mono** for sequences, k-mers, transcript IDs, count values, abundance values.
- Default state is instructive: the artifact opens showing a meaningful example, no user input required.
- Controls grouped in a panel above or to the left of the visualization.
- No animations longer than ~400 ms. Motion only when it carries information.

### 1.3 Interaction model

- **Input / sequences / thresholds** — editable text fields, sliders, or dropdowns, validated against the artifact's input alphabet.
- **Step / Play / Pause / Reset** — where the artifact shows an algorithm with discrete steps (STAR simulator, EM visualizer).
- **Speed** — optional slider, 0.25×–4×. Default 1×.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact in this lecture answers its own question at the end. The student should never be left to infer the result from the final animation state. Concretely:

- If the artifact walks a pipeline, it shows the **output** at each stage (BAM row, count row, abundance estimate) — not just the process.
- If the artifact simulates alignment, it shows the **final alignment** explicitly with a verdict: "✓ aligned with 2 splices" or "✗ no valid alignment".
- If the artifact computes a compatibility class, it shows the **class set** explicitly as a list of transcript IDs.
- If the artifact runs EM, it shows the **converged abundance estimate** and the total log-likelihood margin.
- If the artifact computes a normalisation, it shows the **before-vs-after values** side by side.

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts that accept user input (read sequences, count matrices) must pre-flight the input and report *why* it is or isn't runnable before the calculation runs. Rejected inputs get an inline explanation.

### 1.6 Pedagogical constraint

Every artifact must produce a **specific realization** — the "target aha moment" named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed and should be revised.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry or analytics.
- No external data files larger than ~50 KB.

---

## 2. Artifact #1 — RNA-seq Pipeline Walkthrough

**File**: `artifacts/lecture-05/01-pipeline-walkthrough.html`
**Lecture anchor**: §1.4 The RNA-seq pipeline at a glance
**EE framing reinforced**: pipeline-as-stages; each stage's output is the next stage's input.

### Teaching purpose

Turn the abstract "RNA-seq pipeline" from Figure 2 into concrete data at each stage. Student should leave able to name what goes into and out of every stage of both the align-then-count path and the pseudoalign path.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Pipeline path: [Align-then-count ▾ / Pseudoalign]           │
│ Sample: [Default · Bulk RNA-seq, 4 genes, 20 reads]         │
├─────────────────────────────────────────────────────────────┤
│ Stages (click to inspect):                                  │
│  [FASTQ]  →  [QC]  →  [Align]  →  [BAM]  →  [Counts]  →  [DE│
│                                                             │
│  (selected stage's data shown below)                        │
├─────────────────────────────────────────────────────────────┤
│ Stage inspector:                                            │
│   Input to this stage: [FASTQ reads rendered with Q scores] │
│   Operation: [description]                                  │
│   Output: [next-stage data]                                 │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Pipeline produced 4-gene count matrix               │   │
│ │   Gene1: 12 · Gene2: 8 · Gene3: 3 · Gene4: 5         │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Pipeline path selector**: Align-then-count / Pseudoalign.
- **Sample presets**: 3 presets (default, with-ambiguous-reads, low-coverage).
- **Stage selector**: click any stage box to see what the data looks like there.
- **Reset**: returns to the initial preset.

### What they see

- Top: a horizontal pipeline flow, matching Figure 2's structure.
- Middle: stage inspector showing the data at the clicked stage — FASTQ text at stage 1, BAM-like entries at stage 3/4, a small count table at stage 5.
- Bottom: final count matrix rendered as a small table; outcome banner with gene-level totals.

### Target aha moment

Switch from "Align-then-count" to "Pseudoalign" on the same input. Watch the middle stages change — BAM disappears, transcript-level abundances appear instead, and `tximport` is the bridge to gene-level counts. The final count matrix is very close but not identical. Student sees that the two paths converge at the endpoint but traverse very different data structures in between.

### Technical notes

- Pure JS, no libraries. All data is hardcoded per preset.
- Stage transitions are instant (no animation); click-to-inspect is the interaction.

### Acceptance criteria

- [ ] Default pipeline walkthrough shows reasonable data at every stage.
- [ ] Switching paths (align vs pseudoalign) updates the middle stages correctly.
- [ ] Each preset produces a different, coherent count-matrix output.
- [ ] Opens standalone with the default stage pre-selected.
- [ ] HTML parses; inline JS passes `node --check`.

---

## 3. Artifact #2 — STAR Seed-and-Extend Simulator

**File**: `artifacts/lecture-05/02-star-simulator.html`
**Lecture anchor**: §2.2 STAR and the spliced-read HMM
**EE framing reinforced**: HMM with long-skip transitions at canonical motifs.

### Teaching purpose

Make STAR's seed-cluster-extend process visible. Student should see how MMPs are found, clustered, and extended across introns.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Read: [editable 100-bp sequence]                            │
│ Reference: [two exons + intron, preset]                     │
│ Options: [✓] canonical-GT-AG-bonus                          │
│          [✓] allow splicing                                 │
│          [──●── min MMP length = 20 bp]                     │
│ [Step] [Play] [Reset]                                       │
├─────────────────────────────────────────────────────────────┤
│ Algorithm progress:                                         │
│   Phase 1 · Seed search — 3 MMPs found                      │
│   Phase 2 · Clustering — 2 clusters (one discarded)         │
│   Phase 3 · Extending — splicing at GT-AG motif             │
├─────────────────────────────────────────────────────────────┤
│ Visualisation (SVG): read above, reference below with MMPs  │
│                       highlighted and cluster arcs drawn   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Aligned with 1 splice at pos 1050 (GT-AG motif)    │   │
│ │   Match score: 95 · Splice penalty: 0                │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Read input**: editable 100–150 bp sequence (presets available).
- **Canonical motif bonus toggle**: when off, splicing is penalised everywhere; when on, GT-AG motifs allow free splicing.
- **Allow splicing toggle**: when off, STAR runs in unspliced mode; shows how splicing breaks.
- **Min MMP length slider**: 10–30 bp.
- **Step / Play / Reset** for the three-phase algorithm.

### What they see

- SVG visualisation of the read and reference, MMPs highlighted as the algorithm finds them, clustering arcs drawn as the cluster step runs, splice jumps rendered in the extension phase.
- A small trace log showing each phase's output.
- Outcome banner reporting the final alignment score and number of splices.

### Target aha moment

Toggle off "canonical-GT-AG-bonus" on a preset read that spans a junction. The alignment either fails or produces a much lower score because the splice is now penalised. Student sees that the motif recognition is what makes spliced alignment tractable — without it, the aligner would have to try every possible intron length at every position.

### Technical notes

- Pure JS. Presets embed read + reference + intron position.
- MMP search implemented naively for clarity (linear scan over reference for each seed).
- Smith-Waterman extension with splice bonus at GT-AG motifs.

### Acceptance criteria

- [ ] Default preset produces correct seed-cluster-extend visualisation with one splice.
- [ ] Custom reads are accepted and processed.
- [ ] Toggling canonical-motif-bonus visibly changes alignment quality.
- [ ] Step/Play/Reset controls work.
- [ ] Opens showing the default preset pre-aligned.
- [ ] HTML parses; inline JS passes `node --check`.

---

## 4. Artifact #3 — Compatibility Class Viewer

**File**: `artifacts/lecture-05/03-compatibility-classes.html`
**Lecture anchor**: §3.1 The compatibility-class insight
**EE framing reinforced**: sketch-based set membership; pseudoalignment's core abstraction.

### Teaching purpose

Make compatibility classes concrete. Student types a read (or picks a preset), sees it decomposed into k-mers, sees each k-mer's transcript-set from the hash, and sees the intersected compatibility class.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Transcriptome: [Preset · 3 transcripts, ~80 bp each]        │
│   T1: ACGTACGTACGTACGTACGTACGTACGTACGT... (editable)       │
│   T2: CGTACGTACGTACGTGGGCCGTACGTACGTACG... (editable)       │
│   T3: TAGTACGTACGTACGTAAATACGTACGTACGTA... (editable)       │
│ k: [──●── 11 bp]                                            │
├─────────────────────────────────────────────────────────────┤
│ Read: [editable · default: ACGTACGTACGTACGT]                │
│ [Pseudoalign]                                               │
├─────────────────────────────────────────────────────────────┤
│ k-mer decomposition and transcript-set membership:          │
│  k-mer         | belongs to                                 │
│  ACGTACGTACG    | {T1, T2}                                  │
│  CGTACGTACGT    | {T1, T2}                                  │
│  GTACGTACGTA    | {T1, T2, T3}                              │
│  ...                                                        │
│                                                             │
│  Intersection = compatibility class = {T1, T2}              │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Compatibility class: {T1, T2}                       │   │
│ │   Read compatible with 2 of 3 transcripts             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Transcriptome editor**: three transcripts editable; up to 200 bp each.
- **k slider**: 5–21 bp (default 11).
- **Read input**: editable 15–30 bp.
- **Presets**: "unique match {T1}", "multi-match {T1,T2,T3}", "no match {}", "novel k-mer".
- **Pseudoalign button**: triggers the computation.

### What they see

- The transcriptome with per-transcript colour coding.
- The read's k-mers listed in a table with their transcript-set tags.
- The intersection result shown live as set operations proceed.
- Outcome banner stating the compatibility class.

### Target aha moment

Load the "novel k-mer" preset. One of the read's k-mers doesn't appear in any transcript. The intersection of sets including {} is {} — the read drops. Student sees that a single non-matching k-mer is enough to disqualify the whole read. Contrast with a read that has ambiguous k-mers in each position but still has a non-empty intersection.

### Technical notes

- Pure JS. k-mer hash is a simple Map<string, Set<transcript>>.
- Intersection computed as set intersection; rendered incrementally k-mer by k-mer.
- Feasibility gate: invalid characters in transcripts or read produce inline errors.

### Acceptance criteria

- [ ] Default preset produces a non-trivial compatibility class.
- [ ] Custom transcriptome and reads work.
- [ ] Empty compatibility class is rendered correctly (with an explanation).
- [ ] k slider changes the decomposition in real time.
- [ ] Opens showing a default preset with the class already computed.
- [ ] HTML parses; inline JS passes `node --check`.

---

## 5. Artifact #4 — EM Iteration Visualizer

**File**: `artifacts/lecture-05/04-em-visualizer.html`
**Lecture anchor**: §3.4 Expectation-Maximization for read-to-transcript assignment
**EE framing reinforced**: iterative soft-decision decoding.

### Teaching purpose

Let the student step through EM iterations on a small problem, watch the fractional votes redistribute, and watch the transcript abundance estimates converge.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Scenario: [Preset ▾ · 3 transcripts, 10 reads, default]     │
│   T1, ℓ=500 | T2, ℓ=1500 | T3, ℓ=2000                       │
│   Reads: 10 with compatibility classes listed               │
│ [Step EM] [Play] [Pause] [Reset] [Run to convergence]       │
│                                                             │
│ Iteration: [5]                                              │
├─────────────────────────────────────────────────────────────┤
│ Fractional vote matrix (reads × transcripts):               │
│   r1  T1: 0.48  T2: 0.52                                    │
│   r2  T2: 1.00                                              │
│   ...                                                       │
│                                                             │
│ Current abundance estimate:                                 │
│   θ = (0.22, 0.58, 0.20)                                    │
│   Bar chart showing the three bars.                         │
│                                                             │
│ Log-likelihood: -42.3 (climbing)                            │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Converged at iter 12 · log-L margin = 8.7          │   │
│ │   Final θ = (0.21, 0.60, 0.19)                       │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Preset selector**: 3 scenarios (default balanced, heavy-skew, fully-ambiguous).
- **Custom mode**: textarea for entering read compatibility classes manually.
- **Step EM button**: one iteration per click.
- **Play/Pause**: ~500ms per iteration.
- **Run to convergence**: skip animation, go straight to final state.
- **Reset**: restart from uniform initial abundance.

### What they see

- A matrix of fractional votes updated each iteration, with cells shaded by vote magnitude.
- A bar chart of the abundance estimate θ below the matrix, updated each iteration.
- A running log-likelihood curve showing convergence.
- Outcome banner with converged θ and log-likelihood margin vs. alternatives.

### Target aha moment

Load the "fully-ambiguous" scenario — where every read is compatible with all transcripts. Watch: the iterations *don't converge* to a unique solution — any θ that satisfies the constraints is a local maximum. The student sees that compatibility classes need to span multiple subsets for EM to discriminate between transcripts.

### Technical notes

- Pure JS. EM implemented straightforwardly (E-step computes votes, M-step normalises).
- Convergence detected by Δθ < 10⁻⁶.
- Custom input: parse comma-separated compatibility classes like `T1,T2|T2|T1,T2,T3|T3`.
- Feasibility gate: reject empty or malformed input with inline error.

### Acceptance criteria

- [ ] Default preset converges in ~10–20 iterations to reasonable θ.
- [ ] "Fully-ambiguous" preset shows non-convergence / degenerate solution.
- [ ] Step / Play / Reset controls all work.
- [ ] Log-likelihood climbs monotonically each iteration.
- [ ] Opens showing the default preset at iteration 0 (uniform θ).
- [ ] HTML parses; inline JS passes `node --check`.

---

## 6. Artifact #5 — Normalisation Calculator

**File**: `artifacts/lecture-05/05-normalization.html`
**Lecture anchor**: §4.2 CPM, TPM, FPKM — what each normalises for
**EE framing reinforced**: unit normalisation as a removal of systematic biases.

### Teaching purpose

Given a small count matrix with gene lengths and per-sample depths, compute CPM, FPKM, TPM, and DESeq2-style size-factor-adjusted values, all at once. Student should see where the units differ and when each is the right choice.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Count matrix (editable):                                    │
│   Gene   Len(bp)   S1      S2      S3                       │
│   G1     500       120     150     100                      │
│   G2     2000      400     480     350                      │
│   G3     5000      200     250     180                      │
│   G4     1200      80      95      70                       │
│   G5     800       60      72      50                       │
│   [+ add gene] [+ add sample] [reset to preset]             │
├─────────────────────────────────────────────────────────────┤
│ Normalised values (side-by-side):                           │
│   Gene | CPM                | TPM                |FPKM      │
│   G1   | S1=1200 S2=1500... | S1=456 S2=495...  | S1=2400...│
│   G2   | ...                 | ...               | ...      │
│                                                             │
│ Size factors (DESeq2):                                      │
│   S1 = 0.91 · S2 = 1.08 · S3 = 1.00                         │
│                                                             │
│ Size-factor-adjusted counts:                                │
│   G1: S1=131 S2=139 S3=100                                  │
│   ...                                                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 5 genes × 3 samples normalised                        │   │
│ │ G3 has highest CPM in S2, but 2nd-highest TPM         │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Count matrix editor**: spreadsheet-like, click a cell to edit.
- **Add/remove genes and samples**.
- **Reset to preset**: restore the default 5-gene, 3-sample example.
- **Unit toggle**: highlight CPM, FPKM, TPM, or size-factor-adjusted counts by clicking a unit label (sorts the gene list by that unit).

### What they see

- The raw count matrix.
- A side-by-side rendering of all four normalisations, computed live on input change.
- Size factors computed via median-of-ratios, shown below.
- Outcome banner highlighting one or two genes where the rank differs across units.

### Target aha moment

Start with the default preset. Highest-CPM gene is G2 in S2. Now toggle to TPM: the highest TPM is a different gene, because TPM corrects for transcript length. G3 (5000 bp) dominated CPM because it's long and collects lots of reads; in TPM, G1 (short and abundant per-base) pulls ahead. Student sees that the three units *rank genes differently* — which is the whole point.

### Technical notes

- Pure JS. Median-of-ratios computed live.
- CPM, FPKM, TPM formulas as in §4.2 of the lecture.
- Feasibility gate: non-numeric or negative counts, non-positive lengths → inline error.

### Acceptance criteria

- [ ] Default preset produces matching CPM, TPM, FPKM, and size-factor values.
- [ ] Editing a cell updates all four normalisation outputs in real time.
- [ ] Adding a gene or sample produces sensible output.
- [ ] Switching unit highlights shows different gene rankings.
- [ ] Opens showing the default preset with all four normalisations visible.
- [ ] HTML parses; inline JS passes `node --check`.

---

## 7. Artifact #6 — Poisson vs Negative Binomial Explorer

**File**: `artifacts/lecture-05/06-nb-vs-poisson.html`
**Lecture anchor**: §4.4 The count distribution and why Poisson isn't enough
**EE framing reinforced**: overdispersion as multiplicative noise on top of Poisson.

### Teaching purpose

Let the student see how NB differs from Poisson interactively, and compare both against a real RNA-seq mean-variance scatter. Reader should leave with an intuitive grasp of when Poisson breaks and why NB replaces it.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Parameters:                                                 │
│   mean μ: [──●── 50]                                        │
│   dispersion α: [──●── 0.2]                                 │
│   Observations to simulate: [──●── 200]                     │
│ [Re-sample]                                                 │
├─────────────────────────────────────────────────────────────┤
│ Discrete distributions at μ = 50:                           │
│   ┌──────────────┐                                          │
│   │  Poisson     │   ▁▂▃▅███▅▃▂▁  mean 50, var 50           │
│   │              │                                          │
│   │  Neg Binom   │   ▁▂▃▄▅▆██▇▆▅▄▃▂  mean 50, var 550       │
│   └──────────────┘                                          │
├─────────────────────────────────────────────────────────────┤
│ Empirical mean-variance plot (real RNA-seq):                │
│   [log-log scatter of ~2000 gene points]                    │
│   [Poisson diagonal overlaid]                               │
│   Points visibly sit ABOVE the Poisson line.                │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Poisson: var = mean (50)                              │   │
│ │ NB: var = mean + α·mean² = 50 + 0.2·2500 = 550        │   │
│ │ Real RNA-seq dispersion ≈ 0.02–0.3 depending on gene │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Mean slider**: 1–500 (log scale).
- **Dispersion α slider**: 0 to 1.0 (α = 0 collapses NB to Poisson).
- **Sample size slider**: 50–1000.
- **Re-sample button**: draws new samples and re-renders.

### What they see

- Two histograms side by side, both at the same mean μ. Poisson is always tight; NB widens as α increases.
- A real mean-variance scatter from an embedded RNA-seq dataset (~2000 points), with the Poisson reference line drawn.
- Outcome banner showing the analytical variance for both distributions and the typical empirical dispersion range.

### Target aha moment

Crank α up from 0 to 0.3. The NB histogram widens dramatically while the Poisson stays put. Simultaneously, the empirical mean-variance scatter fits the NB curve well above Poisson. Student sees that (1) Poisson's variance-equals-mean assumption is violated by real data, and (2) α is not a small correction — it's the dominant variance contributor at moderate-to-high expression.

### Technical notes

- Pure JS. Sample from Poisson using Knuth's method; sample from NB using gamma-Poisson composition.
- Embedded RNA-seq mean-variance data: ~2000 (mean, variance) pairs extracted from a GTEx-like dataset. ~50 KB.
- Feasibility gate: α < 0 → inline error.

### Acceptance criteria

- [ ] Default (μ = 50, α = 0.2) shows clear NB widening over Poisson.
- [ ] α = 0 collapses NB to Poisson visually.
- [ ] Real RNA-seq scatter sits above the Poisson line.
- [ ] Re-sampling produces visibly different but consistent histograms.
- [ ] Opens showing the default parameters with both distributions rendered.
- [ ] HTML parses; inline JS passes `node --check`.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Transcript IDs are rendered in JetBrains Mono consistently: `T1`, `T2`, etc.
- Count values and abundance values use JetBrains Mono.
- The Compatibility Class Viewer (#3), T-DBG concept in lecture, and EM Visualizer (#4) use the same transcript-set colour palette.
- The Normalisation Calculator (#5) and NB Explorer (#6) use consistent histogram styling.
- Every artifact emits an **outcome banner** per convention §1.4.

## 9. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>` (C6 check).
- [ ] Outcome banner or equivalent verdict line visible at the end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
