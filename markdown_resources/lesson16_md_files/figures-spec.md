# Lecture 16 — Figures Specification

> **Scope**: Static diagrams for Lecture 16 (ML in Genomics: Architectures, Pitfalls, Frontiers).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-16.md`.

---

## 0. Conventions for This Lecture

- Figures are custom SVG; content is diagram / table / schematic heavy rather than plot-heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Architecture diagrams: rounded rectangles for tensor / module blocks with tensor-shape annotations in JetBrains Mono.
- Colour conventions per architecture family: CNN amber, transformer cobalt, VAE violet, GNN teal, equivariant structure amber-dark. Consistent across all figures.
- Typography: Inter for UI labels; JetBrains Mono for tensor shapes, dataset counts, code-like details (k-mer lengths, loss names).
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | DeepVariant pileup-CNN | Part 1 | Custom SVG |
| 2 | Enformer architecture | Part 1 | Custom SVG |
| 3 | Architecture-to-problem map | Part 1 | Custom SVG |
| 4 | Inductive-bias table | Part 2 | Custom SVG |
| 5 | Labelled-data pyramid in genomics | Part 3 | Custom SVG |
| 6 | Self-supervised pretraining flow | Part 3 | Custom SVG |
| 7 | Data leakage — random vs chromosome split | Part 4 | Custom SVG |
| 8 | DNA-LM landscape | Part 5 | Custom SVG |
| 9 | Cell foundation model architecture | Part 6 | Custom SVG |
| 10 | CNN vs Transformer vs MLP comparison | Part 2 | Custom SVG |
| 11 | Leakage case studies — before / after splits | Part 4 | Custom SVG |
| 12 | Multimodal foundation-model landscape | Part 7 | Custom SVG |

---

## Figure 1 — DeepVariant pileup-CNN

**File**: `diagrams/lecture-16/01-deepvariant-cnn.svg`
**Lecture anchor**: §1.2 Pileups → Convolutional Neural Network
**ViewBox**: `0 0 1200 520`

### Purpose

Show the DeepVariant architecture: pileup image input → CNN → genotype classification.

### Content

**Left panel — pileup image**. A 2D grid of ~15 rows × 20 columns. Rows = reads. Columns = genomic positions around a candidate variant (highlighted with a vertical amber band at column 10). Each cell colour-coded by base (A cobalt, C amber, G green, T red) with a small darkening for lower-quality reads. Strand shown via thin arrow on the left of each row.

**Middle panel — CNN stack**. Rounded rectangles labelled:

- "Conv 3×3 × 32" with output shape "15 × 20 × 32".
- "Conv 3×3 × 64" → "7 × 10 × 64" (after max-pool).
- "Conv 3×3 × 128" → "3 × 5 × 128".
- "Dense 512".
- "Dense 4 (softmax)".

Arrows between blocks.

**Right panel — output**. Four horizontal bars labelled "hom-ref / het / hom-alt / no-variant" with probabilities 0.02 / 0.95 / 0.02 / 0.01. The "het" bar highlighted.

**Bottom caption**. "Translation-invariant local feature extraction matches pileup structure: mismatches are local 2D patterns; column-in-window shouldn't matter."

### Style notes

- Pileup image: small colour-coded cells.
- CNN blocks: amber outlines.
- Output bars: cobalt.

---

## Figure 2 — Enformer architecture

**File**: `diagrams/lecture-16/02-enformer.svg`
**Lecture anchor**: §1.3 Long-range sequence regulation
**ViewBox**: `0 0 1200 520`

### Purpose

Show Enformer's dilated-conv + transformer hybrid architecture for 100 kb regulatory prediction.

### Content

**Top band — Input**. A horizontal strip labelled "DNA sequence, 100 kb × 4 one-hot". Colour-coded bases.

**Middle band — Dilated conv stack**. A horizontal row of blocks labelled "dilation 1", "dilation 2", "dilation 4", "dilation 8", "..." up to "dilation 512". Each block annotated with its receptive-field contribution. Arrows show output feeding forward. Total receptive field: ~2 kb.

**Middle-lower band — Transformer blocks**. 11 transformer blocks stacked, each labelled "multi-head self-attention + FFN". Tensor shape: down-sampled to "1536 positions × 768 channels". Annotation: "attention routes from distal enhancers to gene body".

**Bottom band — Output heads**. Three output rectangles:

- "Chromatin accessibility (per cell type)" — 1536 × 684.
- "CAGE gene expression (per cell type)" — 1536 × 638.
- "ChIP-seq TF binding" — 1536 × ~2000.

**Right-side annotation box**. "Dilated conv = efficient local-to-mid scale; attention = sparse long-range routing. Full 2D attention at 100 kb would be O(10¹⁰) — infeasible."

### Style notes

- Dilated conv blocks: cobalt gradient (dark → light with increasing dilation).
- Transformer blocks: cobalt-solid.
- Output heads: teal.

---

## Figure 3 — Architecture-to-problem map

**File**: `diagrams/lecture-16/03-architecture-map.svg`
**Lecture anchor**: §1.6 Protein 3D (wrap-up of Part 1)
**ViewBox**: `0 0 1200 600`

### Purpose

A 5-row summary table connecting genomic problems to architectures and inductive biases.

### Content

**5-column table**. Columns:

1. **Problem**
2. **Input representation**
3. **Architecture**
4. **Key inductive bias**
5. **Example model**

**Rows**:

- Variant calling / pileup image / 2D CNN / translation invariance + locality / **DeepVariant** (amber).
- Long-range regulation / one-hot DNA, 100 kb / dilated conv + transformer / multi-scale + sparse long-range / **Enformer** / **Borzoi** (cobalt).
- Single-cell expression / cells × genes counts / VAE + NB likelihood / NB likelihood + latent prior / **scVI** / **totalVI** (violet).
- Molecular property / atom-bond graph / GNN (message-passing) / permutation invariance + local aggregation / **GIN** / **SchNet** (teal).
- Protein structure / MSA + pair representation / equivariant transformer / axial attention + SE(3) equivariance / **AlphaFold2** (amber-dark).

Each row colour-tagged with the architecture family's signature colour.

### Style notes

- Header row: bold Inter with dark background.
- Cells in alternating tint for readability.
- Model names in JetBrains Mono.

---

## Figure 4 — Inductive-bias table

**File**: `diagrams/lecture-16/04-inductive-biases.svg`
**Lecture anchor**: §2.6 Summary table
**ViewBox**: `0 0 1200 600`

### Purpose

Cross-reference between architectures and the inductive biases they encode, with data-regime guidance.

### Content

**3-column table with ~8 rows**. Columns:

1. **Architecture**
2. **Key inductive bias**
3. **Data regime where it wins**

**Rows** (one per architecture):

- MLP / none / very large, no structure.
- CNN / translation invariance + locality / small-to-medium, position-independent local features.
- Dilated conv + transformer / multi-scale + long-range / structured, moderate-scale sequences.
- VAE + NB likelihood / NB count distribution + latent prior / count data with batch effects.
- GNN / permutation invariance + locality / molecular graphs.
- Equivariant transformer / SO(3)-equivariance / 3D geometric data.
- Transformer (base) / weak locality, attention over all positions / very large data.
- Plain fully-connected / none / never in high-dim data.

### Style notes

- Each row colour-coded.
- "Where it wins" column in italic.
- Annotation at bottom: "Match inductive bias to data structure. Mismatch is the #1 cause of genomics ML failure."

---

## Figure 5 — Labelled-data pyramid in genomics

**File**: `diagrams/lecture-16/05-data-pyramid.svg`
**Lecture anchor**: §3.1 Where labels come from in genomics
**ViewBox**: `0 0 1080 540`

### Purpose

Show the inverted-pyramid / funnel shape of genomics label availability.

### Content

**Inverted pyramid** (widest at top, narrowest at bottom). Label each tier with count and example:

1. Top (widest): "Unlabelled sequences (UniRef, SRA, NCBI) — ~10¹⁰ entries". Label quality: "implicit only".
2. "Bulk RNA-seq (SRA, GEO) — ~10⁶ samples". Label quality: "tissue / condition tags, often noisy".
3. "GWAS associations — ~10⁵ SNP × trait pairs". Label quality: "p-values, effect sizes".
4. "Single-cell transcriptomes (HCA, cellxgene) — ~10⁸ cells". Label quality: "cluster labels, imperfect".
5. "ClinVar variants — ~10⁶". Label quality: "expert-curated pathogenicity".
6. Bottom (narrowest): "Deep mutational scans / MPRA — ~10⁴ variants per gene, single-digit thousands of genes". Label quality: "experimental ground truth, gold".

**Right side** — annotation box. "Genomics ML is a self-supervised pretraining + fine-tuning field by necessity: the top is where scale lives; the bottom is where labels live."

### Style notes

- Pyramid in cobalt gradient (lightest at top, darkest at bottom).
- Counts in JetBrains Mono.
- Label-quality tags in italic.

---

## Figure 6 — Self-supervised pretraining flow

**File**: `diagrams/lecture-16/06-ssl-flow.svg`
**Lecture anchor**: §3.3 Self-supervised pretraining
**ViewBox**: `0 0 1200 480`

### Purpose

Left-to-right workflow of self-supervised pretraining → fine-tuning for a downstream task.

### Content

**Panel 1 (left)**: A large cloud labelled "Unlabelled sequence corpus (~650M proteins / ~3B DNA bp)". Inside, a few example sequences shown as strips.

**Panel 2**: A transformer block with a masked-token pretraining objective: one sequence with ~15% of tokens replaced by `[MASK]`; the model predicts the original. Arrow labelled "MLM pretraining, millions of steps".

**Panel 3**: "Pretrained model ready". Icon of a saved model checkpoint.

**Panel 4**: A second, smaller cloud labelled "Labelled downstream task (~50k samples)". Examples: "variant pathogenic / benign", "promoter / non-promoter".

**Panel 5**: Fine-tuning arrow → final deployment model.

**Bottom caption**. "Pretraining uses 100,000× more data than fine-tuning. Representations transfer. Downstream labelled-sample needs drop by ~10–100×."

### Style notes

- Clouds for data stores.
- Transformer block with rectangular shape.
- Arrows connecting panels.

---

## Figure 7 — Data leakage — random vs chromosome split

**File**: `diagrams/lecture-16/07-leakage-splits.svg`
**Lecture anchor**: §4.1 Why genomics leaks worse
**ViewBox**: `0 0 1200 480`

### Purpose

Contrast random-position splitting (leaky) vs chromosome splitting (clean).

### Content

**Top panel — Random split**. A chromosome bar with tiny alternating amber (train) and cobalt (test) cells distributed across it, like a chess pattern. Annotation arrow: "test positions nearby train positions → LD bleeds".

**Bottom panel — Chromosome split**. A chromosome bar with chr1–chr21 entirely amber (train) and chr22 entirely cobalt (test). Annotation arrow: "test is spatially separated → no LD bleed".

**Right-side summary table**. Three rows:

- **Random split**: "test F1 = 0.95 (inflated)".
- **Chromosome split**: "test F1 = 0.68 (honest)".
- **Family-aware split**: "test F1 = 0.65 (even stricter)".

### Style notes

- Amber for train; cobalt for test.
- Chromosome bar styled as a horizontal rectangle.
- Inflation arrows in red for the random-split case.

---

## Figure 8 — DNA-LM landscape

**File**: `diagrams/lecture-16/08-dna-lm-timeline.svg`
**Lecture anchor**: §5.1 What DNA LMs are
**ViewBox**: `0 0 1200 520`

### Purpose

Show the chronological evolution of DNA language models, by year and by context length.

### Content

**Axes**: X = year (2021–2024). Y = max context length in bp (log scale, 10² – 10⁷).

**Points** (each a labelled dot):

- 2021: DNABERT — 512 tokens × 6-mer = ~3 kb effective; 100M params. Bottom-left.
- 2023: Nucleotide Transformer — 1000 tokens × 6-mer = ~6 kb; 500M–2.5B params. Mid-left.
- 2023: HyenaDNA — 1 Mb context; 100M params. Top-right.
- 2024: Evo — 131k context; 7B params. Upper-right.

**Shaded region "transformer-feasible"** along the x-axis below 10 kb; **shaded region "sub-quadratic models"** above 100 kb.

**Annotation**. "Context length is the main axis of competition; sub-quadratic architectures (Hyena) extend beyond transformer-feasible."

### Style notes

- Dots sized by parameter count.
- Region shading: `--bg-subtle` vs slightly-tinted cobalt.
- Labels in JetBrains Mono.

---

## Figure 9 — Cell foundation model architecture

**File**: `diagrams/lecture-16/09-cell-foundation-model.svg`
**Lecture anchor**: §6.1 The pitch
**ViewBox**: `0 0 1200 520`

### Purpose

A schematic of a cell foundation model (scGPT / Geneformer / scFoundation style).

### Content

**Left — Input**. A cell's expression vector over ~20,000 genes. Above, small illustration of three tokenisation choices: rank-value, binned-expression pairs, real-valued embedding. Only one chosen at a time per model.

**Middle — Transformer stack**. A vertical stack of ~12 transformer blocks labelled "self-attention + FFN". Per-block shape: "n_genes × hidden_dim".

**Right — Output heads**. Three possible downstream heads:

- Cell-type classifier.
- Perturbation response prediction.
- Gene-gene interaction attribution.

**Bottom band — Pretraining corpus details**. "Pretraining: 30M+ cells from CellxGene / HCA; cross-tissue / cross-species / cross-technology".

**Right-side annotation box (coloured red/orange)**. "Caveats: specialised methods (scVI, Harmony, CellTypist) still win on specialised downstream tasks. Foundation-model utility over well-tuned specialised methods is unsettled."

### Style notes

- Transformer blocks: cobalt.
- Output heads: teal (one per column).
- Caveats box: amber border.

---

## Figure 10 — CNN vs Transformer vs MLP comparison

**File**: `diagrams/lecture-16/10-architecture-comparison.svg`
**Lecture anchor**: §2.2 CNN vs MLP vs Transformer on pileup-like data
**ViewBox**: `0 0 1200 520`

### Purpose

Quantify the data-efficiency ordering of MLP < CNN < Transformer (untuned) < Transformer (local-attention-biased) on a pileup-like task across training-set sizes.

### Content

**Line plot**. X-axis: training set size, log scale from 10³ to 10⁷. Y-axis: validation F1 on a synthetic pileup task, range 0.5–1.0.

**Four curves**:

- MLP: starts ~0.70 at 10³, plateaus ~0.90 at 10⁷. (Persistent underfitting.) Grey.
- CNN: jumps to 0.98 by 10⁴, flat at 0.995 thereafter. Amber.
- Transformer (vanilla, absolute positional embeddings): starts 0.70 at 10³, climbs through 0.90 at 10⁵, 0.99 at 10⁷. Cobalt.
- Transformer (2D positional + local-attention inductive bias): starts 0.85 at 10³, matches CNN at 10⁵, pulls ahead to 0.999 at 10⁷. Cobalt-dark.

**Crossover annotations**: "CNN ahead below 10⁵; transformer with good priors ahead above 10⁶".

**Bottom caption**. "Inductive bias is data-regime-dependent. The right prior buys 10×–100× sample efficiency."

### Style notes

- Curves with distinct markers per architecture.
- Log-scale x-axis with gridlines at each decade.

---

## Figure 11 — Leakage case studies — before / after splits

**File**: `diagrams/lecture-16/11-leakage-cases.svg`
**Lecture anchor**: §4.6 Real-world case studies
**ViewBox**: `0 0 1080 480`

### Purpose

Real-world examples where tightening train/test splits dramatically changed the reported headline accuracy.

### Content

**Grouped bar chart**. Three groups on the x-axis (case studies). Per group, two bars: "As reported" (tall, cobalt) and "After tightened split" (short, amber).

- Case 1: "Protein function prediction (2015 era)". As reported: 95%. Tightened: 45%.
- Case 2: "Cancer driver prediction (2016 era)". As reported: 99%. Tightened: 70%.
- Case 3: "Enformer chromosome-split (2021)". As reported: 0.82 Pearson. Tightened: 0.74 Pearson.

**Y-axis**: performance metric (0–1.0). Each group uses its own scale noted inline.

**Bottom annotation band**. "Reproducing published genomics-ML papers should start by auditing the split. The corrected number is usually the honest one."

### Style notes

- Tall bars: cobalt.
- Short bars: amber.
- Gap between bars with a difference label (e.g. "−50pp").

---

## Figure 12 — Multimodal foundation-model landscape

**File**: `diagrams/lecture-16/12-multimodal.svg`
**Lecture anchor**: §7.1 Multimodal foundation models
**ViewBox**: `0 0 1200 560`

### Purpose

Map existing foundation models in a 2D space of (modalities covered) × (parameter count); highlight the frontier.

### Content

**Axes**. X-axis: modalities covered (1 → 2 → 3+). Y-axis: parameter count (log, 10⁸ → 10¹¹).

**Dots** (each labelled):

- DNABERT (1 modality: DNA, ~10⁸). Bottom-left.
- Nucleotide Transformer (1 modality: DNA, ~10⁹). Mid-left.
- Geneformer (1 modality: single-cell, ~10⁷). Bottom-left.
- scGPT (1 modality: single-cell, ~10⁸). Bottom-left.
- Evo (1 modality: DNA, ~10¹⁰). Upper-left.
- ESM-2 (1 modality: protein, ~10¹⁰). Upper-left.
- AlphaFold3 (3 modalities: protein + nucleic + ligand, ~10⁹). Upper-right.
- ESM-3 (3 modalities: seq + struct + function, ~10¹¹). Top-right.

**Shaded "frontier" region** in the top-right (high-modality, high-param).

**Arrow** from bottom-left to top-right labelled "The direction of progress".

**Bottom annotation band**. "Data — matched multi-modal measurements — is the bottleneck, not model scale."

### Style notes

- Dots sized by parameter count.
- Model-family colour-coding (DNA cobalt; protein amber; cell teal; multi-modal violet).
- Shaded frontier region: `--bg-subtle` with a soft outline.
