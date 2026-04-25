# Lecture 20 (proposed L4) — Figures Specification

> **Scope**: Static diagrams for Lecture 20 (MSA, Phylogenetics, Comparative Genomics).
> **Companion files**: `lecture-20.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Sequence colours: amino-acid hydrophobicity gradient (matches L15 protein conventions).
- Tree colours: edges grey; bootstrap labels green if ≥ 70, amber 50–70, red < 50.
- Selection regime colours: positive selection red, purifying cobalt, neutral grey.
- Synteny colours: conserved-block green, rearranged amber, inverted red.
- Typography: Inter for chrome; JetBrains Mono for sequences, codons, residue scores.

## Figure budget

Thirteen figures.

| # | Title | Part |
|---|---|---|
| 1 | Pairwise → MSA generalization | Part 1 |
| 2 | Progressive alignment step-by-step | Part 2 |
| 3 | Neighbor-Joining algorithm step-by-step | Part 3 |
| 4 | Molecular clock dating | Part 4 |
| 5 | dN/dS distribution under selection regimes | Part 5 |
| 6 | Synteny dot plot (human–mouse) | Part 6 |
| 7 | Standard phylogenetics workflow | Part 7 |
| 8 | MSA quality benchmark comparison | Part 2 |
| 9 | Felsenstein's pruning algorithm | Part 3 |
| 10 | Ortholog vs paralog distinction | Part 6 |
| 11 | Long-branch attraction failure mode | Part 3 |
| 12 | Tree of life (vertebrates excerpt) | Part 7 |
| 13 | RNA secondary structure and SCFG model | Part 6.5 |

---

## Figure 1 — Pairwise → MSA generalization

**File**: `diagrams/lecture-20/01-pairwise-to-msa.svg`
**ViewBox**: `0 0 1200 540`

Three panels:

- Left: pairwise DP matrix (10×10 cells), traceback diagonal highlighted.
- Middle: 5-sequence DP tensor in 5D — visualised as nested cubes; total state count $10^{10}$ labelled; "exact algorithm impossible" annotation.
- Right: progressive alignment showing decomposition into 4 pairwise/profile-vs-sequence steps; total cost $\sim N^2 L^2$.

---

## Figure 2 — Progressive alignment step-by-step

**File**: `diagrams/lecture-20/02-progressive-msa.svg`
**ViewBox**: `0 0 1200 600`

5-sequence MSA build:

- Step 1: align two closest sequences (sea cucumber vs lamprey) → first 2-row alignment.
- Step 2: profile from step 1 + chicken → 3-row alignment.
- Step 3: profile from step 2 + mouse → 4-row.
- Step 4: profile from step 3 + human → final 5-row alignment.

Conserved positions highlighted in colour (functional residues like the heme-binding histidine).

---

## Figure 3 — Neighbor-Joining step-by-step

**File**: `diagrams/lecture-20/03-nj-tree.svg`
**ViewBox**: `0 0 1200 720`

Six-taxon example:

- Top: input pairwise distance matrix (6×6 colour-coded heat map).
- Middle: Q-matrix at iteration 1, with smallest entry circled.
- Iterations 2-3: progressively-collapsed matrices.
- Bottom: final NJ tree with branch lengths.
- Annotation: "Q-matrix correction removes molecular-clock assumption."

---

## Figure 4 — Molecular clock dating

**File**: `diagrams/lecture-20/04-clock-dating.svg`
**ViewBox**: `0 0 1200 540`

Phylogenetic tree of 5 mammalian species (human, chimpanzee, mouse, rat, dog).

- Branch lengths in substitutions-per-site.
- Calibration point: chimpanzee-human split fixed at 6 Mya.
- Time axis on bottom showing absolute dates.
- Posterior 95% credible intervals on internal nodes (relaxed clock).

---

## Figure 5 — dN/dS distribution under selection regimes

**File**: `diagrams/lecture-20/05-dnds.svg`
**ViewBox**: `0 0 1200 540`

Three panels of $\omega$ distributions:

- Panel 1: ω distribution under purifying selection — sharp mode at 0.1; long tail truncated at 1.
- Panel 2: ω distribution under positive selection — bimodal; main mass < 1, secondary mass > 1.
- Panel 3: real-world examples — HA (broad with positive tail), histones (narrow at 0), housekeeping (centered around 0.1).

Bottom annotation: "ω = dN/dS gives quantitative selection detection at gene and site level."

---

## Figure 6 — Synteny dot plot

**File**: `diagrams/lecture-20/06-synteny.svg`
**ViewBox**: `0 0 1200 720`

Two-genome dot plot:

- X-axis: human chromosome 17 (0–80 Mb).
- Y-axis: mouse chromosome 11 (0–120 Mb).
- Dots: homolog gene pairs.
- Visible diagonal stripes = conserved synteny blocks.
- Off-diagonal points = rearrangements (small).

Side annotations call out specific gene clusters (HOXB cluster, p53 region, CDKN1B). Top inset: dot-plot overview at full chromosome × chromosome scale.

---

## Figure 7 — Standard phylogenetics workflow

**File**: `diagrams/lecture-20/07-workflow.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart:

1. Collect homologs (BLAST / OrthoFinder).
2. MAFFT L-INS-i alignment.
3. Trim with trimAl / Gblocks.
4. IQ-TREE model selection + ML tree + 1000 bootstraps.
5. Visualise on iTOL / FigTree.
6. PAML for selection / BEAST for dating.

Side panel: alternative tools at each step.

---

## Figure 8 — MSA quality benchmark comparison

**File**: `diagrams/lecture-20/08-msa-benchmark.svg`
**ViewBox**: `0 0 1200 540`

Bar chart of MSA quality (sum-of-pairs accuracy on BAliBASE) for six tools:

- ClustalW (1994).
- ClustalO (2011).
- MUSCLE (2004).
- MAFFT FFT-NS-2 (2002).
- MAFFT L-INS-i (2007).
- PROBCONS (2005).

Y-axis: average accuracy across BAliBASE reference sets.
Annotation at the bottom: "MAFFT L-INS-i and PROBCONS are top; speed-accuracy trade-off is real."

---

## Figure 9 — Felsenstein's pruning algorithm

**File**: `diagrams/lecture-20/09-felsenstein.svg`
**ViewBox**: `0 0 1200 600`

A 5-leaf binary tree.

- Leaves: observed residues (e.g., 4 trees with leaves A, A, G, T, A).
- At each internal node: dynamic-programming computation: $L(\text{node}) = \prod_{\text{children}} \sum_{\text{state}} P(\text{child} \mid \text{parent state}, t) L(\text{child})$.
- Bottom-up flow shown as arrows.
- Final root likelihood as the marginal $\sum_{\text{state}} P(\text{root state}) L(\text{root}, \text{state})$.

Annotation: "Belief propagation on the tree."

---

## Figure 10 — Ortholog vs paralog distinction

**File**: `diagrams/lecture-20/10-ortho-paralog.svg`
**ViewBox**: `0 0 1200 540`

Toy gene-family tree spanning 3 species (human, mouse, dog) and 2 sub-families (alpha, beta) created by an ancient duplication:

- Internal duplication node clearly marked.
- All speciation nodes after duplication.
- Pairs labelled: Human-α / Mouse-α = orthologs (speciation only). Human-α / Human-β = paralogs (duplication). Human-α / Mouse-β = "ohnologs" or distant homologs.

Side panel: classification table summarising relationships.

---

## Figure 11 — Long-branch attraction failure mode

**File**: `diagrams/lecture-20/11-lba.svg`
**ViewBox**: `0 0 1200 540`

Two trees side by side:

- True tree: A-B (close pair); C-D (close pair); long branch from A and from C to deep root.
- Parsimony-inferred tree (incorrect): A-C close; B-D close; topology error.

Annotation: "long-branch attraction occurs when fast-evolving lineages share many substitutions by chance — parsimony incorrectly groups them."

---

## Figure 12 — Tree of life (vertebrate excerpt)

**File**: `diagrams/lecture-20/12-tree-of-life.svg`
**ViewBox**: `0 0 1200 800`

Simplified vertebrate phylogeny:

- Lampreys → fish (cartilaginous + bony) → tetrapods (amphibia + reptilia + birds + mammals).
- Mammals expanded: monotremes → marsupials → placentals (rodents, primates, carnivores, etc.).
- Branch lengths approximate to time since divergence.
- Calibration points (fossils) marked.
- Annotation: "modern tree-of-life work uses 100s of marker genes + relaxed clock; Hug et al. 2016 placed archaea as sister to eukaryotes."

---

## Figure 13 — RNA secondary structure and SCFG model

**File**: `diagrams/lecture-20/13-rna-structure.svg`
**ViewBox**: `0 0 1200 720`

Three panels:

- Panel 1 (top-left): tRNA cloverleaf diagram with anticodon arm (red), acceptor stem (cobalt), D-arm (green), TΨC arm (amber); secondary structure clearly drawn.
- Panel 2 (top-right): same tRNA in dot-bracket notation: `(((((..(((....))).(((....))).....(((....)))))))).` aligned with the sequence beneath.
- Panel 3 (bottom): profile-CM (covariance model) for the tRNA family — match-pair states (paired columns), match-singlet states (loop columns), insert and delete states. SCFG production rules schematically shown.

Annotation: "SCFGs handle RNA's nested base-pairing structure exactly via the inside-outside algorithm. Pseudoknots make folding NP-hard; classical RNA folding ignores them."

Side panel: example Infernal `cmsearch` output line: family ID, E-value, bit score, strand, structure conservation.
