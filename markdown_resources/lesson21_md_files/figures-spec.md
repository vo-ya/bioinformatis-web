# Lecture 21 (proposed L7) — Figures Specification

> **Scope**: Static diagrams for Lecture 21 (HMMs, Profile HMMs, Gene Finding).
> **Companion files**: `lecture-21.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- HMM-state colours: emission heat in cobalt-orange gradient; transitions grey.
- Gene-finding state classes: intergenic grey, exon green, intron amber, UTR muted-cobalt, splice-site red.
- Chromatin-state colours: active-promoter red, enhancer amber, transcribed green, polycomb purple, heterochromatin black-grey.
- Typography: Inter for chrome; JetBrains Mono for state labels, codons, scores.

## Figure budget

Twelve figures.

| # | Title | Part |
|---|---|---|
| 1 | HMM graph and joint-probability factorisation | Part 1 |
| 2 | Viterbi trellis with backtrace | Part 2 |
| 3 | Profile HMM as MSA → 3-state-per-column model | Part 3 |
| 4 | GENSCAN HMM architecture | Part 4 |
| 5 | HMM ↔ deep learning landscape | Part 6 |
| 6 | Forward-backward as belief propagation | Part 2 |
| 7 | Baum-Welch convergence | Part 2 |
| 8 | Splice-site PWM logo | Part 4 |
| 9 | ChromHMM emission matrix | Part 5 |
| 10 | Pair-HMM in variant calling | Part 5 |
| 11 | Pfam domain hit on a protein | Part 3 |
| 12 | End-to-end annotation pipeline | Part 7 |

---

## Figure 1 — HMM graph and joint-probability factorisation

**File**: `diagrams/lecture-21/01-hmm-graph.svg`
**ViewBox**: `0 0 1200 540`

Top half: a 4-state HMM as a directed graph. States as nodes; transitions as labelled edges with probability values. Each state has a small bar chart showing emission distribution.

Bottom half: the unrolled-in-time graphical model for $T = 6$ observations. States $z_1, ..., z_6$ in a horizontal chain; observations $x_t$ hanging below each state. Edges showing conditional dependencies.

Annotation: "$P(X, Z) = \pi_{z_1} \prod a_{z_{t-1}, z_t} \prod b_{z_t}(x_t)$ — Markov + emission factorisation."

---

## Figure 2 — Viterbi trellis with backtrace

**File**: `diagrams/lecture-21/02-viterbi-trellis.svg`
**ViewBox**: `0 0 1200 600`

A 4-state HMM unrolled over $T = 8$ time steps.

- Trellis lattice: 4 rows × 8 columns.
- Each cell shows $\log \delta_t(i)$ as colour intensity.
- Edges between cells show transitions.
- Viterbi backtrace path drawn in bold red from end to start.
- Side panel: pseudocode of Viterbi recurrence.

Annotation: "Viterbi = max-plus DP on the trellis; same machinery as pairwise alignment."

---

## Figure 3 — Profile HMM as MSA → 3-state-per-column model

**File**: `diagrams/lecture-21/03-profile-hmm.svg`
**ViewBox**: `0 0 1200 600`

Top: a 5-row × 8-column MSA (toy globin alignment).

Bottom: the corresponding profile HMM with 8 match-state columns:

- Each match state shows its emission distribution as a small 20-bar amino-acid histogram.
- Insert states drawn as self-loops above each match.
- Delete states drawn as bypass arrows below.
- B (begin) and E (end) states at the ends.

Side annotation: "Match emissions = column composition; transitions = gap statistics."

---

## Figure 4 — GENSCAN HMM architecture

**File**: `diagrams/lecture-21/04-genscan.svg`
**ViewBox**: `0 0 1200 720`

A high-level state diagram of GENSCAN:

- Intergenic state (grey).
- Promoter / 5'-UTR (muted cobalt).
- Initial exon, internal exons, final exon (green).
- Intron-phase 0 / 1 / 2 (amber).
- 3'-UTR + polyA signal (muted cobalt).
- Single-exon gene (green).

Arrows show canonical transitions. Length distributions noted next to states (geometric for intron, gamma for exon).

Annotation: "Viterbi on this HMM segments a 100 Mb chromosome into a list of genes in minutes."

---

## Figure 5 — HMM vs deep learning landscape

**File**: `diagrams/lecture-21/05-hmm-vs-dl.svg`
**ViewBox**: `0 0 1200 540`

2D scatter:

- X-axis: training-set size (log scale, $10^2$–$10^9$).
- Y-axis: prediction accuracy.

Two curves:

- HMM (cobalt): rises early, plateaus around $10^4$ examples.
- Deep learning (red): rises late, exceeds HMM around $10^5$.

Crossover point marked with annotation: "HMM regime: limited data + interpretability. DL regime: abundant data + flexibility. Most genomic-classification tasks sit near the crossover."

---

## Figure 6 — Forward-backward as belief propagation

**File**: `diagrams/lecture-21/06-forward-backward.svg`
**ViewBox**: `0 0 1200 540`

Chain of hidden states $z_1, ..., z_T$ with observations.

- Forward messages $\alpha_t(i)$ propagate left-to-right (cobalt arrows).
- Backward messages $\beta_t(i)$ propagate right-to-left (amber arrows).
- Per-position posterior $\gamma_t(i) = \alpha_t(i) \beta_t(i) / P(X)$ shown as combined messages.

Annotation: "Forward-backward IS belief propagation on a chain. Felsenstein's pruning is BP on a tree."

---

## Figure 7 — Baum-Welch convergence

**File**: `diagrams/lecture-21/07-baum-welch.svg`
**ViewBox**: `0 0 1200 540`

Plot of $\log P(X \mid \theta^{(k)})$ vs iteration $k$ for 5 different random initialisations:

- 3 converge to the global optimum.
- 2 to a local optimum at lower likelihood.
- Y-axis: log-likelihood.
- Annotation: "EM monotonically increases likelihood. Multiple restarts mitigate local optima."

---

## Figure 8 — Splice-site PWM logo

**File**: `diagrams/lecture-21/08-splice-pwm.svg`
**ViewBox**: `0 0 1200 480`

Two side-by-side splice-site logos:

- Donor (5' splice site): position-weight matrix from -3 to +6 around the GT consensus. Y-axis: information content (bits). Tall G and T at positions 1 and 2 of the intron.
- Acceptor (3' splice site): PWM from -14 to +3, including polypyrimidine tract upstream of the AG consensus.

Annotation: "PWMs = simplest profile HMM (no inserts, no deletes). Modern SpliceAI generalizes these to deep CNN."

---

## Figure 9 — ChromHMM emission matrix

**File**: `diagrams/lecture-21/09-chromhmm-emission.svg`
**ViewBox**: `0 0 1200 720`

Heat map of 15-state ChromHMM model:

- Rows: 15 states (active promoter, weak promoter, strong enhancer, weak enhancer, ...).
- Columns: 9 histone marks (H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K27me3, H3K9me3, H3K9ac, H4K20me1, CTCF).

Heat map values: $b_i(\text{mark}) = P(\text{mark observed} | \text{state})$.

Side panel: the genome-wide proportion of each state in human ESCs.

Annotation: "States are distinguished by their mark-co-occurrence patterns. ChromHMM finds these by Baum-Welch on whole-genome ChIP-seq."

---

## Figure 10 — Pair-HMM in variant calling

**File**: `diagrams/lecture-21/10-pair-hmm.svg`
**ViewBox**: `0 0 1200 540`

Three states: Match (M), Insert (I), Delete (D). Edges:

- M → M: $1 - 2\delta - \tau$ (continue match).
- M → I: $\delta$ (insertion).
- M → D: $\delta$ (deletion).
- I → I: $\epsilon$ (extend insertion).
- I → M: $1 - \epsilon$.
- D → D: $\epsilon$ (extend deletion).
- D → M: $1 - \epsilon$.
- M → END: $\tau$.

Annotation: "Pair-HMM is the engine behind GATK HaplotypeCaller and DeepVariant. Computes per-read alignment likelihood with explicit indel handling."

---

## Figure 11 — Pfam domain hit on a protein

**File**: `diagrams/lecture-21/11-pfam-hit.svg`
**ViewBox**: `0 0 1200 540`

A protein sequence of ~600 residues drawn as a horizontal bar.

- Identified Pfam domains highlighted as coloured boxes:
  - Kinase domain (250 aa, position 200-450), e.g., PF00069.
  - SH2 domain (100 aa, position 50-150), e.g., PF00017.
  - PH domain (100 aa, position 480-580).
- HMMER E-value annotation per domain.

Side panel: the profile-HMM logos for the three domains, showing position-specific conservation.

Annotation: "domain assignment is the first-pass functional annotation for any new protein."

---

## Figure 12 — End-to-end annotation pipeline

**File**: `diagrams/lecture-21/12-pipeline.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart for bacterial genome annotation:

1. Raw reads (FASTQ).
2. Assembly (SPAdes / Flye).
3. Repeat masking.
4. ORF prediction (Prodigal).
5. Protein translation.
6. Pfam scan (HMMER).
7. UniProt search (DIAMOND).
8. KEGG pathway mapping.
9. Final annotation table.

Each step annotated with typical runtime + tool name. Side panel: equivalent eukaryotic pipeline (AUGUSTUS / BRAKER + RNA-seq hints).
