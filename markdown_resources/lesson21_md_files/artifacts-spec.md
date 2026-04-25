# Lecture 21 (proposed L7) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 21 (HMMs, Profile HMMs, Gene Finding).
> **Companion files**: `lecture-21.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-21/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- HMM-state colours: emission heat in cobalt-orange; transitions grey.
- Gene-finding state classes: intergenic grey, exon green, intron amber, UTR muted-cobalt, splice-site red.
- Chromatin-state colours: active-promoter red, enhancer amber, transcribed green, polycomb purple, heterochromatin black-grey.
- Typography: Inter for chrome; JetBrains Mono for state labels, codons, scores.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget

Seven interactive tools.

| # | Title | Anchor |
|---|---|---|
| 1 | HMM playground | §1 |
| 2 | Viterbi stepper | §2.3 |
| 3 | Baum-Welch convergence | §2.4 |
| 4 | Profile HMM builder | §3 |
| 5 | GENSCAN-style gene finder simulator | §4 |
| 6 | End-to-end annotation walkthrough | §7 |
| 7 | ChromHMM segmentation | §5.1 |

---

## Artifact #1 — HMM Playground

**File**: `artifacts/lecture-21/01-hmm-playground.html`
**Anchor**: §1

### Teaching purpose

Build a 2-state CpG island HMM. Tune transition + emission parameters; see Viterbi-decoded segmentation update live.

### UI layout

- Two-state HMM diagram with editable transition probabilities.
- Per-state emission distribution: 4 nucleotide bars (+ optional 16-bin dinucleotide).
- 1 kb input sequence (preset toggleable; user-editable for advanced).
- Output: Viterbi-decoded state-per-position sequence visualisation.
- Statistics: number of CpG islands called, total fraction of genome in island state.
- Outcome banner: "with these parameters: 12 islands found, average length 280 bp."

### Target aha

Simple 2-state HMM segments a sequence; the segmentation depends on transitions (how often switch states) and emissions (how strongly each state prefers different observations).

### Acceptance criteria

- Default opens with realistic CpG island parameters.
- Decoding updates within 200 ms.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — Viterbi Stepper

**File**: `artifacts/lecture-21/02-viterbi-stepper.html`
**Anchor**: §2.3

### Teaching purpose

Step through Viterbi decoding cell-by-cell on the trellis. Each step shows max-plus computation and pointer storage.

### UI layout

- 3-state HMM (preset).
- 8-character input sequence.
- Trellis lattice: 3 rows × 8 columns.
- Step-forward button advances one cell.
- Current cell highlighted; pre-existing cells coloured by $\log \delta$ value.
- Backtrace button at the end shows the optimal path.
- Outcome banner: "step 6 / 24: at cell ($t = 3, S = $ State 2): $\log \delta_3(2) = -4.2$ from predecessor State 1."

### Target aha

Viterbi is a max over all paths, computed bottom-up via the trellis. Pointers store argmax for backtrace.

### Acceptance criteria

- Step controls work; backtrace paints the path correctly.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Baum-Welch Convergence

**File**: `artifacts/lecture-21/03-baum-welch.html`
**Anchor**: §2.4

### Teaching purpose

Watch Baum-Welch (EM) iterate from random initialisation; see log-likelihood monotonically increase to a local optimum.

### UI layout

- Generative HMM defined (3-state with known parameters).
- Synthetic observations generated.
- 5 random initial parameter sets, each tracked separately.
- Log-likelihood trajectories plotted vs iteration.
- Run / pause / step buttons.
- Outcome banner: "after 50 iterations: 3 of 5 inits converged to global optimum; 2 trapped at local optima with -8% lower likelihood."

### Target aha

EM is local. Multiple restarts are essential; pick the highest-likelihood result.

### Acceptance criteria

- Visible monotonicity in each curve.
- 30+ iteration default.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Profile HMM Builder

**File**: `artifacts/lecture-21/04-profile-hmm.html`
**Anchor**: §3

### Teaching purpose

Construct a profile HMM from a small MSA; visualise the 3-state-per-column architecture; align a query sequence.

### UI layout

- Pre-loaded 5-row MSA.
- Generated profile HMM displayed with match emissions, insert states, delete states.
- Query sequence input (preset / custom).
- Viterbi alignment of query to model.
- Output: alignment showing match / insert / delete states traversed; HMMER-style bit score and E-value.
- Outcome banner: "query aligns with bit score 47, E = 1e-12 → strong domain match."

### Target aha

Profile HMM = MSA → 3-state-per-column model. Viterbi alignment scores any new sequence against the family.

### Acceptance criteria

- 6-row MSA building correct profile.
- Query alignment visualisation animates.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — Gene Finder Simulator

**File**: `artifacts/lecture-21/05-gene-finder.html`
**Anchor**: §4

### Teaching purpose

Run a simplified GENSCAN-style HMM on a 10 kb genomic input. Visualise the gene segmentation.

### UI layout

- 10 kb input genomic region (toggle: "with hints" vs "without hints").
- Multi-state HMM: intergenic, promoter, exon-initial / internal / final, intron-phase 0/1/2, UTR-5/3.
- Viterbi-decoded segmentation as a horizontal track with colour-coded states.
- Slider: intron-length distribution mean (controls gene-finding sensitivity to intron length).
- Statistics: gene count, total CDS length, average exon size.
- Outcome banner: "found 3 genes; total CDS 1.4 kb; mean exon 320 bp."

### Target aha

Eukaryotic gene finding is multi-state HMM segmentation. Length distributions and emission patterns drive the result.

### Acceptance criteria

- Default input contains 3 simulated genes; tool finds all 3.
- Toggling "with hints" improves boundary precision.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — End-to-End Annotation Walkthrough

**File**: `artifacts/lecture-21/06-annotation-walkthrough.html`
**Anchor**: §7

### Teaching purpose

Walk a 1 kb mock contig through assembly → ORF prediction → Pfam scan → BLAST function transfer → KEGG pathway.

### UI layout

- Step-by-step pipeline (5 stages, each clickable).
- At each step: input visualised, tool's action explained, output shown.
- Final output: gene annotation card with name, function, Pfam domain, pathway.
- Outcome banner: "annotated contig contains 1 gene → 'putative protein kinase' → KEGG pathway 'MAPK signaling'."

### Target aha

Annotation is a synthesis: each tool contributes one piece; the final annotation is the union.

### Acceptance criteria

- 5 steps, each with non-trivial output.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — ChromHMM Segmentation

**File**: `artifacts/lecture-21/07-chromhmm.html`
**Anchor**: §5.1

### Teaching purpose

Run a simplified 15-state ChromHMM on a 1 Mb region with synthetic ChIP-seq tracks.

### UI layout

- 6 ChIP-seq tracks (H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K27me3, H3K9me3).
- ChromHMM emission matrix displayed.
- 15-state segmentation track over the 1 Mb region.
- Toggle: which states to highlight (active vs poised vs repressed).
- Outcome banner: "active promoters 2.4% of genome; enhancers 4.1%; transcribed 18%."

### Target aha

ChromHMM integrates multi-track ChIP-seq into interpretable chromatin states; the per-state emission matrix tells you what each state represents.

### Acceptance criteria

- Realistic emission matrix.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- HMM-state visualisations consistent across #1, #2, #3.
- Profile HMM topology in #4 matches Figure 3 in lecture.
- ChromHMM colours in #7 match figure 9 in lecture.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
