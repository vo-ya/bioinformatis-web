# Lecture 7 — Single-Cell RNA-seq Fundamentals

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-07.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Explain what bulk RNA-seq (Lecture 5) averages away, and name two biological questions that can only be answered at single-cell resolution.
2. Sketch a 10x Genomics droplet-based single-cell RNA-seq workflow end-to-end, and explain what cell barcodes and UMIs are doing.
3. Describe why UMIs solve the PCR-duplicate problem that bulk RNA-seq tolerates, and relate UMIs to error-correcting tags.
4. Read a knee plot and set an EmptyDrops-style threshold separating real cells from ambient barcodes.
5. Name three QC filters that every single-cell analysis applies (doublets, mitochondrial fraction, ambient RNA) and explain what each removes.
6. Explain why bulk normalisation methods fail on sparse single-cell data, and describe the alternative (SCTransform / log-normalisation + HVG selection).
7. Describe the dimensionality-reduction pipeline: PCA → UMAP/t-SNE → graph-based clustering (Leiden/Louvain) → marker gene identification.
8. Annotate a cluster of cells as a cell type using either marker-based or reference-based methods, and state the tradeoff between them.

---

## Part 1 — Why Single-Cell and the Droplet Chemistry (≈45 min)

### 1.1 What bulk RNA-seq hides (≈8 min)

Lecture 5 covered bulk RNA-seq: take millions of cells, extract RNA, sequence. The output is an expression profile averaged across the whole tube. For a tissue like liver or blood, that average is a weighted mean over many distinct cell types — hepatocytes and macrophages and endothelial cells all contributing to the same number.

Three biological questions the bulk average can't answer:

- **Which cell type expresses the gene I care about?** A single bulk value for "TNF in whole blood" doesn't say whether TNF is being made by monocytes, T cells, or neutrophils — and that matters for therapy.
- **Do rare cells exist that express something nobody else does?** A tumor-infiltrating regulatory T cell at 0.5% of the cells contributes 0.5% to any bulk average. Its distinct expression program is invisible.
- **Are cells continuously transitioning, or discrete?** Development, immune activation, differentiation — these look like smooth gradients when you sample one cell at a time, but averaging turns the gradient into a single static value.

Single-cell RNA-seq (scRNA-seq) answers these by measuring transcript abundance in each cell separately. The output is a **cell × gene count matrix**: thousands of cells as rows, ~20,000 genes as columns, entries are molecule counts.

> **Intuition box**: Bulk is a mean pool over thousands of cells. Single-cell is every individual sample. The cost: the per-cell signal is orders of magnitude noisier, because each cell yields ~10⁴ mRNAs at best and you sample only a fraction. The payoff: cellular heterogeneity becomes visible.

### 1.2 Cell-type heterogeneity and the "average cell" fallacy (≈7 min)

The empirical claim that drove the field: **almost no tissue is homogeneous.** Even a supposedly pure cell population — say, CD4+ T cells sorted to 99% purity from peripheral blood — turns out to contain five or ten functionally distinct states when you look cell-by-cell. Naïve T cells, effectors, central memory, regulatory, exhausted, and a few transitional states whose existence wasn't even hypothesised before scRNA-seq.

This heterogeneity means the bulk average is not a good representative of any single cell. It is, in fact, a value that no single cell has. Treating the average as "what the cell is doing" is the **ergodic fallacy** — assuming the ensemble mean tells you the per-sample behaviour.

**FIGURE — Figure #1: Single-cell vs bulk — what the average hides** → `diagrams/lecture-07/01-bulk-vs-single-cell.svg`
*Top: a tissue with 4 visibly distinct cell types, shown schematically. Middle: the bulk expression value for each gene, rendered as a single averaged trace. Bottom: the per-cell expression for three genes, showing that distinct cell types express distinct gene sets.*

> **EE framing**: Bulk RNA-seq is a low-resolution integrator — it sums over all cells and reports the mean. Single-cell is a high-resolution sampler that resolves individual signal sources. The same measurement trade-off exists in radio astronomy (beamforming averages sky-emissions vs interferometry that resolves sources), microscopy (bulk fluorescence vs single-molecule imaging), and signal processing (power-spectrum estimation vs time-domain sample inspection). Going single-cell trades per-sample SNR for resolution — and most biological questions turn out to need the resolution.

### 1.3 The 10x Genomics droplet chemistry (≈15 min)

10x Genomics' Chromium platform dominates single-cell RNA-seq in 2024. It's a droplet-based method: cells are encapsulated in oil-water emulsion droplets alongside a barcoded gel bead, each cell's mRNA is reverse-transcribed into labelled cDNA inside its own droplet, then everything is pooled back together and sequenced in bulk.

The chemistry in four steps:

1. **Load.** A suspension of ~1000–50,000 cells is loaded alongside a large excess of gel beads, each bead carrying millions of copies of a 16–18 bp cell barcode (the same barcode across all copies on one bead; different barcodes across beads). A microfluidic chip co-encapsulates one cell + one bead + oil into each droplet (called a GEM, for gel-bead-in-emulsion).
2. **Lyse + prime.** Inside each droplet, the cell is lysed. mRNAs bind to oligo-dT primers on the bead. Each primer carries: the bead's cell barcode + a 10–12 bp unique molecular identifier (UMI) + a poly-T tail.
3. **Reverse transcribe.** Reverse transcriptase extends each primer into a cDNA copy of the mRNA, with the cell barcode and UMI locked to the front.
4. **Pool + sequence.** Droplets are broken, all cDNA pooled, PCR-amplified, and sequenced on Illumina. Every read carries its cell barcode (which cell it came from), its UMI (which mRNA molecule it came from), and its transcript sequence.

Typical Chromium v3 run: 10,000 cells, 50,000 reads per cell, ~3000 genes detected per cell. Total cost in 2024: ~$3,500 per sample plus sequencing.

**FIGURE — Figure #2: 10x Genomics droplet chemistry** → `diagrams/lecture-07/02-droplet-chemistry.svg`
*A schematic of the Chromium microfluidic chip: three input channels (cells, gel beads, oil) converging into a droplet former; a single GEM shown in cross-section with a cell, a bead, and the oil boundary; the RT step inside the droplet producing labelled cDNA.*

> **Warning box**: Not all single-cell protocols are 10x. Smart-seq2 and Smart-seq3 are plate-based methods that sequence each cell individually in a well — higher gene coverage per cell, lower throughput, no UMIs. Drop-seq is the open-source droplet ancestor of 10x. CEL-seq / MARS-seq use barcoded primers with a different chemistry. Each platform has distinct sensitivity, throughput, and cost profiles, and analysis pipelines are not interchangeable. When reading a paper, check the platform before comparing numbers.

### 1.4 Cell barcodes and UMIs (≈15 min)

The cell barcode says **which cell** the read came from. The UMI says **which mRNA molecule** within that cell.

A cell barcode is 16 bp (10x Chromium v3), drawn from a fixed list of ~6 million valid barcodes. After sequencing, reads are grouped by barcode — all reads sharing the same barcode came from the same cell. Some of those barcodes will be error-corrected to the nearest valid barcode in the list (Hamming-1 distance).

A UMI is 10–12 bp. It's a random sequence attached to each primer on the bead before the cell is even captured. When reverse transcription runs, each individual mRNA in the cell gets tagged with a different UMI. After PCR amplification, you end up with many copies of the same cDNA — each copy carrying the original UMI from the one mRNA molecule it came from.

**The critical move**: during quantification, collapse reads that share both cell barcode *and* UMI. They came from the same original mRNA molecule; counting them separately would inflate the count. The final count per gene per cell is the number of *distinct UMIs*, not the number of reads.

**FIGURE — Figure #3: Cell barcode + UMI structure and deduplication** → `diagrams/lecture-07/03-barcode-umi.svg`
*Top: a bead with primers laid out — each primer has (cell barcode, UMI, poly-T) structure. Middle: inside a droplet, 3 mRNA molecules bind 3 different primers (same cell barcode, different UMIs); PCR amplifies each to 50 copies. Bottom: after sequencing, 150 reads collapse back to 3 distinct (barcode, UMI) pairs — the true count.*

**EMBED — Artifact #1: UMI Deduplication Demo** → `artifacts/lecture-07/01-umi-dedup.html`
*A small set of reads with cell barcodes and UMIs. Watch PCR amplification duplicate them; then watch the count collapse back to the true molecule count by grouping on (barcode, UMI).*

> **EE framing**: A UMI is an error-correcting tag placed on each molecule before the amplification channel. The PCR step is a noisy amplifier (positive feedback with multiplicative noise — callback to Lecture 1's PCR framing). Without a pre-amplification tag, you cannot distinguish "I sequenced one molecule 50 times" from "I sequenced 50 distinct molecules once." The UMI is the serial number stamped on each input before it enters the amplifier; counting distinct serial numbers at the output recovers the original molecule count regardless of amplification bias. This is the same idea as a MAC address or a TCP sequence number — unique identifiers that survive a noisy channel.

> **Historical pointer**: The first UMI design was Kivioja et al. 2012, applied to bulk RNA-seq. The killer application turned out to be single-cell, where PCR amplification is enormous (hundreds of copies per molecule) and without UMIs the counts would be fiction. 10x Genomics' founding insight in 2015 was combining droplet microfluidics (from Drop-seq, Macosko et al. 2015) with commercial-grade manufacturing of UMI-bearing gel beads. The resulting instrument shipped in 2016 and within three years had overtaken every competing platform.

---

## Part 2 — Data Pipeline and QC (≈60 min)

### 2.1 From BCL to a count matrix (≈10 min)

The sequencing instrument (Illumina NovaSeq, typically) emits **BCL** (Binary Base-Call) files directly from the flow cell. A complete pipeline from BCL to a scRNA-seq count matrix:

```
BCL files (Illumina flow cell)
    ↓ bcl2fastq
FASTQ files (reads split by sample, barcodes on Read 1, transcript on Read 2)
    ↓ cellranger count  (or STARsolo / alevin-fry / kb-python)
BAM (aligned reads)  +  filtered_feature_bc_matrix.h5 (the count matrix)
    ↓ QC + normalisation (§2.3-§3)
AnnData / Seurat object (cells × genes, with metadata)
    ↓ clustering + annotation (§4-§5)
Cell-type-labeled expression atlas
```

Each stage has well-defined outputs. A clean pipeline records intermediate artefacts (BCL → FASTQ, FASTQ → BAM, BAM → matrix) so you can re-run any downstream step without re-sequencing or re-aligning.

### 2.2 Quantification tools: Cell Ranger, STARsolo, alevin-fry (≈15 min)

The core job of the quantifier: align reads to the transcriptome (or pseudoalign, callback Lecture 5 §3), group them by (cell barcode, UMI), count distinct UMIs per gene per cell, and emit a sparse count matrix.

Three dominant tools in 2024:

- **Cell Ranger** (10x Genomics' own). Historically the "reference" pipeline. Under the hood, uses STAR (Lecture 5 §2.2) with 10x-specific barcode handling. Slower than modern alternatives, outputs are the de facto standard others are compared against.
- **STARsolo** (part of STAR). Adds 10x / Drop-seq / Smart-seq barcode handling to vanilla STAR. Much faster than Cell Ranger at similar accuracy. The go-to for anyone who isn't locked into Cell Ranger output format.
- **alevin-fry** (from the Salmon team). Pseudoalignment-based (Lecture 5 §3); fastest of the three. Recent enough that not everyone trusts it yet, but benchmarks match STARsolo within ~2% on standard datasets.

Output formats:

- The count matrix is saved as an **HDF5** file (`.h5`) or in 10x's Matrix Market triplet format (`.mtx` + `barcodes.tsv` + `features.tsv`). Both are *sparse* formats — a 10,000-cell × 20,000-gene dataset has ~5% non-zero entries; storing the full dense matrix wastes an order of magnitude of disk.
- The object that downstream analysis tools (Seurat, Scanpy) consume is the **AnnData** (Python) or **SeuratObject** (R) container — essentially "the sparse count matrix plus per-cell and per-gene metadata."

**FIGURE — Figure #4: scRNA-seq quantification pipeline** → `diagrams/lecture-07/04-quant-pipeline.svg`
*A horizontal flow: BCL → FASTQ → (Cell Ranger / STARsolo / alevin-fry branches) → BAM + count matrix (h5) → AnnData/Seurat object. Each stage annotated with its file format.*

### 2.3 Empty droplets and the knee plot (≈10 min)

Every 10x run generates more barcodes than it has cells. Two reasons: (1) some droplets contain only ambient RNA (lysed debris floating in the buffer) with no real cell — they still produce reads, tagged with their bead's barcode; (2) some droplets contain just a bead with no cell at all.

The empty-droplet barcode counts are not zero — they contain small amounts of ambient RNA, typically 10–1000 UMIs per barcode. Real cells usually have 1000–100,000 UMIs. A **knee plot** — barcodes sorted by total UMI count, plotted on log-log — shows a characteristic sharp transition (the "knee") separating the two regimes.

**FIGURE — Figure #5: The knee plot and EmptyDrops threshold** → `diagrams/lecture-07/05-knee-plot.svg`
*Log-log plot: x-axis barcode rank (sorted), y-axis UMI count per barcode. The canonical knee shape: a plateau of high-count "real cell" barcodes on the left, a steep drop in the middle, and a long tail of low-count "empty" barcodes on the right. EmptyDrops' ambient-profile test marks the boundary.*

**EMBED — Artifact #2: Knee Plot Explorer** → `artifacts/lecture-07/02-knee-plot.html`
*A simulated knee plot with adjustable cell count and ambient RNA level. Slide the threshold; see cells kept vs discarded. Run EmptyDrops' ambient-profile logic to place the threshold automatically.*

> **Intuition box**: The knee is a discriminability boundary between two populations — cells and empty droplets — drawn naturally by log-count distribution. On real data the boundary is rarely as clean as a textbook example; 10–30% of "cells" on the right side of the knee are partially-lysed or low-quality. EmptyDrops (Lun et al. 2019) refines the simple knee cut by modelling the ambient profile as a multinomial and testing each barcode against it — cells whose gene distribution looks significantly different from ambient are kept, even if their total UMI count is moderate.

### 2.4 Doublet detection (≈10 min)

When two cells accidentally end up in the same droplet (two cells captured together), their mRNA gets the same cell barcode. The resulting "cell" appears to express the gene programs of both underlying cells — a **doublet**.

Doublet rate scales with cell-loading concentration. 10x's recommended loading produces about 4–8% doublets at 10,000 targeted cells. Load more cells, you get a higher doublet rate (this is a per-chip hardware limitation, not a choice you can skip).

Doublet-detection algorithms work by recognising the "mixture" signature:

- **Scrublet** (Wolock et al. 2019) simulates synthetic doublets by averaging random pairs of real cells; then trains a classifier to distinguish real cells from the synthetic ones; applies the classifier to your data.
- **DoubletFinder** works similarly — simulate, classify.
- **scDblFinder** uses a slightly different simulation + iterative refinement.

All three produce a per-cell doublet score. Thresholds (typically 0.3–0.4) mark high-scoring cells as likely doublets and filter them out.

**FIGURE — Figure #6: Doublet detection strategy** → `diagrams/lecture-07/06-doublet-detection.svg`
*Top: two real cells with distinct expression profiles are computationally combined to make a simulated doublet. Middle: the combined profile overlaid on a UMAP of real cells — it lands between real clusters. Bottom: distribution of doublet scores, with a threshold marked.*

> **EE framing**: Doublet detection is outlier detection in a sparse high-dimensional feature space. Real cells lie on well-separated manifolds (the cell-type clusters); doublets project to the *midpoints between manifolds*, a region that real cells of any type rarely occupy. A classifier trained on the simulated midpoint samples learns to recognise this geometry. The same structural move appears in anomaly-detection pipelines (normal vs. mixed-class intruder detection in network traffic, outlier images in vision datasets).

> **Warning box**: Doublet rate scales linearly with the number of cells loaded per chip. "Load more cells to save money per cell" is a false economy past the 6–8% doublet rate — you save on per-cell cost but lose power from cells you have to discard.

### 2.5 Mitochondrial fraction and ambient RNA (≈15 min)

Two more QC filters that every pipeline runs:

**Mitochondrial fraction.** The percentage of reads mapping to mitochondrial genes (MT-*). Healthy cells have ~5–15% mitochondrial expression. A cell whose mitochondrial fraction is >25% is usually stressed, apoptotic, or had its plasma membrane compromised during dissociation — the cytoplasmic mRNAs leaked out, leaving mostly mitochondrial RNA (protected inside the mitochondrial membrane). Filtering at a per-dataset threshold (commonly 10–20%) removes these low-quality cells.

**Ambient RNA.** In a droplet, some of the mRNA you measure didn't come from your cell — it came from the solution around the cell, released by other lysed cells during dissociation. This ambient contamination is present in every cell to some degree, and shows up as "everyone expresses everything at low levels." Tools that correct for it:

- **SoupX** (Young & Behjati 2020) estimates the ambient-RNA profile from empty droplets and subtracts a fraction of it from each cell's count.
- **CellBender** (Fleming et al. 2023) uses a deep learning model (VAE-flavored) that jointly decomposes the observed counts into (cell signal + ambient background + technical noise).

**FIGURE — Figure #7: QC metrics dashboard** → `diagrams/lecture-07/07-qc-metrics.svg`
*Four panels showing per-cell QC metrics: (a) total UMIs per cell — distribution with knee, (b) genes-per-cell — histogram with low-gene cells in red, (c) MT% — histogram with threshold line, (d) ambient-RNA contamination fraction per cell (post-SoupX).*

**EMBED — Artifact #3: QC Filter Playground** → `artifacts/lecture-07/03-qc-filter.html`
*A simulated scRNA-seq dataset with known-quality cells. Slide thresholds on total UMIs, gene count, and MT% and watch which cells pass. Target aha: different thresholds produce radically different cell populations.*

> **Warning box**: Over-filtering is as damaging as under-filtering. If you set MT% threshold at 5%, you'll exclude cells from tissues with naturally high mitochondrial content (cardiac myocytes, fast-twitch muscle, certain neurons). Always inspect the MT distribution before setting a threshold — the "right" value is dataset-specific.

---

## Part 3 — Normalisation and Feature Selection (≈25 min)

### 3.1 Why bulk normalisation fails on sparse data (≈8 min)

Lecture 5 §4.3 introduced DESeq2's median-of-ratios size factor. It works well for bulk data because (a) each sample has thousands of well-measured genes to compute a robust median over, and (b) library depths are high enough that per-gene counts are typically positive.

Single-cell data breaks both assumptions. A typical scRNA-seq cell has ~3,000 non-zero genes out of ~20,000 total — most genes are zero. The median-of-ratios estimator on mostly-zero vectors becomes unstable: the "median" is dominated by whichever genes happen to be non-zero in both samples, which is a biased subset.

**FIGURE — Figure #8: Sparse count matrix — zeros everywhere** → `diagrams/lecture-07/08-sparse-matrix.svg`
*A 50-cell × 200-gene slice of a real scRNA-seq count matrix rendered as a heatmap. ~95% of cells are dark (zero); the visible structure is patches of expression specific to cell groups.*

> **EE framing**: A scRNA-seq count matrix is a heavily undersampled signal per cell. Each cell measures on the order of 10⁴ molecules out of ~10⁶ in the actual cell; the per-gene signal is shot-noise-limited. The zeros in the matrix are a mixture of "truly not expressed" and "expressed but not captured" — a mixture the normalisation step has to handle without pretending it can distinguish the two.

### 3.2 Log-normalisation and SCTransform (≈10 min)

Two dominant approaches for single-cell normalisation:

**Log-normalisation** (Seurat's `NormalizeData`). Two steps:

1. Divide each cell's counts by its total UMIs (per-cell size factor).
2. Multiply by a scaling factor (typically 10,000), add 1, take natural log.

Result: log-CPM-equivalent values, with the pseudo-count of 1 preventing log(0). Fast, simple, works well for most downstream operations. The price: it ignores the mean-variance relationship of count data (variance scales with mean even after log).

**SCTransform** (Hafemeister & Satija 2019, Choudhary & Satija 2022). A more principled alternative. For each gene, fit a regularised negative-binomial generalised linear model predicting counts from total UMIs per cell. The **Pearson residuals** from this model are used as the normalised expression values.

This is statistically more rigorous than log-norm (handles the mean-variance relationship explicitly, stabilises variance across expression levels). The cost: slower, and the residuals are on a different scale than the raw counts — some downstream steps need adjustment.

In 2024, log-norm is still the default in Seurat and Scanpy workflows because it's fast and "good enough" for most analyses. SCTransform is used when variance-stabilisation matters (e.g., integrating datasets from different batches, §Lecture 8).

**EMBED — Artifact #4: Normalisation Comparison** → `artifacts/lecture-07/04-normalisation.html`
*A simulated scRNA-seq dataset with three cell types at different total UMI counts. Compare raw counts, log-normalisation, and SCTransform-style Pearson residuals side-by-side. See how the distributions change across the three representations.*

> **Warning box**: A zero count in scRNA-seq does not mean "not expressed." A cell with low total UMIs might express a gene that a deeper cell of the same type would show — it just wasn't sampled. This is why scRNA-seq clustering is on *relative* expression patterns, not absolute counts, and why inter-cell expression comparisons are always noisy for low-count genes.

### 3.3 Highly variable gene selection (≈7 min)

A typical scRNA-seq dataset has ~20,000 gene features per cell. Most of them are boring — expressed uniformly across all cell types, contributing no discrimination signal. Keeping all 20,000 wastes compute and, counterintuitively, hurts clustering (high-dimensional random walks drown out the real signal).

The fix: **feature selection** — keep only the top 2,000–5,000 **highly variable genes (HVGs)** for downstream analysis. Two common methods:

- **Variance / mean-based** (Seurat's `FindVariableFeatures` vst method): fit a mean-variance trend across all genes (similar to Lecture 6's dispersion trend); pick the genes with the highest *residual* variance (variance above what the trend predicts).
- **Pearson-residual-based** (SCTransform): the genes with the largest variance in Pearson residuals across cells — essentially the same operation as above, but on the residual scale.

The top 2,000 HVGs capture most of the discriminating signal between cell types. The remaining 18,000 genes contribute noise to clustering and UMAP. Feature selection isn't optional — it's what makes the downstream steps work.

---

## Part 4 — Dimensionality Reduction and Clustering (≈50 min)

### 4.1 Why reduce dimensions (≈5 min)

After HVG selection, each cell is a 2,000-dimensional point. Distance between cells in this space is meaningful (cells with similar expression are close) but we can't visualise 2,000D, and most clustering algorithms that work in 2D or 10D collapse in 2,000D — the **curse of dimensionality**.

> **Intuition box**: In high dimensions, every pair of random points is approximately the same distance apart. "Nearest neighbour" stops being meaningful past ~20–50 dimensions. Any distance-based clustering at 2,000D will be dominated by noise in the many-dimensional marginals. Reducing dimension first is not optional — it's a precondition for any geometric algorithm to work.

Single-cell analysis tackles this with two successive reductions:

1. **PCA** (principal component analysis) — linear, produces ~30–50 components capturing 50–80% of total variance.
2. **UMAP or t-SNE** — nonlinear, produces a 2D visualisation from the PCA components.

Clustering happens on the PCA components, not on the 2D UMAP. UMAP is purely for visualisation.

### 4.2 PCA on the normalised matrix (≈15 min)

PCA on a cells × genes matrix X (after normalisation, HVG selection, scaling to zero mean and unit variance per gene):

1. Compute the covariance matrix C = X^T X / (n-1) (or equivalently, compute SVD of X directly).
2. Find the eigenvectors (principal components) of C — these are directions in gene space along which cells vary maximally.
3. Project cells onto the top k components (k = 30–50 typically): cell i in the new representation is just X[i, :] · V_{:, 1:k}.

The result is a cells × k matrix where k is 30–50 instead of 2,000. Each new dimension is a linear combination of the original genes, ordered by how much variance it captures.

**FIGURE — Figure #9: PCA on scRNA-seq — scree plot and projection** → `diagrams/lecture-07/09-pca-projection.svg`
*Left: scree plot showing eigenvalues of the first 50 PCs — steep drop over the first few, plateau after ~30. Right: cells projected into PC1 × PC2 space, with three cell types visible as distinct clouds.*

> **EE framing**: PCA is SVD of the normalised expression matrix. It produces a low-rank approximation of X where the first k left-singular vectors are the principal components. Same operation as any eigendecomposition-based compression (JPEG basis, MP3 frames, face recognition eigenfaces). The first N PCs are the optimal linear rank-N approximation by Frobenius norm — the same sense in which the first N DFT coefficients are the optimal Fourier-basis approximation.

Practical choice: how many PCs to keep? The scree plot shows diminishing returns past some point. A common heuristic is "keep enough PCs to capture 50–80% of variance" — usually 30–50. More formally, **parallel analysis** compares your scree to one from permuted data and keeps only PCs above the permutation baseline.

### 4.3 UMAP and t-SNE for visualisation (≈15 min)

Once you have cells in 30–50D PC space, you need a 2D visualisation so a biologist can see the clusters. Two tools dominate:

**t-SNE** (van der Maaten & Hinton 2008). Older. Optimises a 2D embedding such that pairwise distances in the low-dimensional space match pairwise distances in the high-dimensional space, where "distance" is measured via Gaussian kernels in high-D and t-distributions in low-D. Preserves local structure (near neighbours stay near) at the cost of distorting global structure (distances between distant clusters are meaningless in t-SNE).

**UMAP** (McInnes et al. 2018). Newer. Based on algebraic topology — treats the data as a fuzzy topological structure and finds a low-dimensional embedding that preserves that structure. In practice: faster than t-SNE, slightly better at preserving global structure (clusters sit at roughly-meaningful distances from each other), similar local-neighbourhood preservation. UMAP is the 2024 default.

**FIGURE — Figure #10: UMAP embedding of a multi-cell-type dataset** → `diagrams/lecture-07/10-umap-embedding.svg`
*A 2D UMAP projection showing ~8,000 cells organised into ~10 visible clusters of different sizes, colour-coded by cluster assignment. Each cluster annotated with its dominant cell-type marker.*

**EMBED — Artifact #5: PCA + UMAP Embedder** → `artifacts/lecture-07/05-pca-umap-embedder.html`
*A simulated 3-cell-type dataset (1,000 cells × 100 genes). Apply PCA; see the scree plot. Then run UMAP on top-10 PCs; see the 2D embedding. Compare to a PCA-only 2D projection. Varying the UMAP perplexity-equivalent parameter reshapes the embedding.*

> **EE framing**: UMAP is a nonlinear manifold-learning algorithm that preserves local topology. Conceptually: treat each cell as a vertex, build a k-nearest-neighbours graph in PC space, then lay out that graph in 2D such that local connectivity is preserved. The analog in signal processing is manifold-based dimensionality reduction (LLE, Isomap) used for speech / gesture / image analysis. UMAP's formal foundation in algebraic topology gives it better global structure preservation than t-SNE's purely-geometric objective.

### 4.4 Graph-based clustering: Leiden and Louvain (≈10 min)

Clustering 10,000 cells in 30D PC space is done graph-theoretically, not with k-means or hierarchical clustering (which would not respect manifold structure):

1. Build a **k-nearest-neighbours graph** (kNN, k=10–30). Each cell is a node; edges connect each cell to its k nearest neighbours in PC space.
2. Optional: build a **shared-nearest-neighbours graph** (SNN) — edge weight is the Jaccard overlap of k-NN neighbour sets. Denoises the raw kNN graph.
3. Run **Leiden** or **Louvain** community detection — both maximise graph modularity (informally: find partitions where edges within clusters are dense and edges between clusters are sparse).

Leiden (Traag, Waltman, van Eck 2019) is strictly better than Louvain — it guarantees well-connected clusters, whereas Louvain can produce disconnected communities. Leiden is the 2024 default.

A critical hyperparameter: the **resolution**. Higher resolution → more clusters (finer subclustering). Lower resolution → fewer, larger clusters. This choice is not automatic; it depends on what biological question you're asking. Common practice: scan resolutions 0.1 to 2.0, pick the one where clusters align with known cell-type markers.

**EMBED — Artifact #6: Leiden Clustering Visualizer** → `artifacts/lecture-07/06-leiden-clustering.html`
*A simulated UMAP embedding of ~500 cells. Build the kNN graph (adjustable k); visualise it overlaid on the UMAP. Run Leiden at adjustable resolution; watch clusters form and dissolve. Target aha: resolution is a knob; more isn't always better.*

> **Discussion prompt**: You run Leiden at resolution 1.0 and get 8 clusters. You raise to 2.0 and get 15 clusters (several of the original 8 split into two or three). Which is "right"? (Neither is automatically right. Raise resolution only if the splits reveal distinct marker genes and map to biologically meaningful cell states. If the sub-clusters have no distinguishing markers, the split is noise. Resolution should follow the biology, not the other way around.)

### 4.5 Marker gene identification (≈5 min)

Once cells are clustered, the final analysis step is **finding markers** — genes that are differentially expressed in a given cluster compared to all other clusters. Markers are what let you assign cell-type identities (§5).

For each cluster *c* and each gene *g*, run a DE test comparing cells in *c* vs cells not in *c*. The output is a log-fold change + adjusted p-value per (cluster, gene). Sort by log-FC; the top genes are the cluster's markers.

Two DE test choices for this:

- **Wilcoxon rank-sum test** (Seurat default, Scanpy default). Nonparametric. Fast. Doesn't assume any count distribution. Robust to scRNA-seq's weird mean-variance relationship.
- **MAST** (Finak et al. 2015). A hurdle model — jointly models the fraction of cells expressing the gene and the expression level in expressing cells. Slower, sometimes more sensitive for low-count markers.

Wilcoxon is the default for speed. MAST is used when the markers really matter.

**FIGURE — Figure #11: Marker gene heatmap** → `diagrams/lecture-07/11-marker-heatmap.svg`
*Rows = cells (grouped by cluster), columns = top 5 marker genes per cluster. Darkness indicates expression. Visible block-diagonal structure: each cluster expresses its own markers distinctly from other clusters.*

> **Intuition box**: A cluster's markers are the genes that say "I'm different from everyone else." If cluster A has as its top marker MS4A1 (= CD20, a B-cell marker), cluster A is probably B cells. If no gene reliably distinguishes cluster A from cluster B, those two clusters are probably over-split at the current resolution.

---

## Part 5 — Cell-Type Annotation (≈20 min)

### 5.1 Reference-based vs marker-based annotation (≈12 min)

Once markers are identified, the final step is **cell-type annotation** — assigning a biological label to each cluster (T cell, hepatocyte, endothelial cell, etc.). Two paradigms:

**Marker-based annotation** (manual). For each cluster, look at its top markers and consult a reference of known cell-type markers:

- B cell markers: MS4A1 (CD20), CD79A, CD19.
- Cytotoxic T cell: CD8A, CD8B, GZMB, PRF1.
- Hepatocyte: ALB, AFP, APOA1, TF.
- Macrophage: CD68, CSF1R, MARCO.

Assign the cluster the label whose markers best match the cluster's top genes. Works well when cell types are well-characterised (immune cells, major organ cell types). Fails when analysing rare or novel populations without documented markers.

**Reference-based annotation** (automated). Given a pre-annotated reference atlas (e.g., Tabula Sapiens — human tissues, 500k cells, cell types annotated by experts), compute the most-similar reference cell for each query cell, transfer the reference's cell-type label:

- **SingleR** (Aran et al. 2019). Per-cell correlation to reference-cell-type average profiles. Fast, stable. The R ecosystem default.
- **Azimuth** (Stuart et al. 2022). Web application running a pre-trained reference for common tissues (PBMC, bone marrow, heart). Upload your data, get annotations in 20 minutes.
- **scArches / scPoli** (Lotfollahi et al. 2022). Deep-learning transfer of cell-type labels via latent-space mapping.

Reference-based methods are faster and less labour-intensive; they work well for well-studied tissues. They fail silently on novel cell types (no reference match means the query is labeled as the closest wrong type). Always validate reference-based annotations against markers.

> **Discussion prompt**: You annotate a dataset with Azimuth against a PBMC reference. One cluster of 300 cells gets labeled "activated NK cell" with high confidence. But the cluster's top marker is CEACAM8 (granulocyte marker), which should not appear in any lymphocyte. What's going on? (The Azimuth reference doesn't contain granulocytes. The 300 cells are in fact granulocytes — a common mistake when running PBMC-reference annotators on samples that were supposed to be PBMC-only but weren't. This is why reference-based annotation should always be paired with marker checking.)

### 5.2 Cell-type ontologies and atlases (≈8 min)

For consistent reporting, cell-type annotations should reference a shared vocabulary:

- **Cell Ontology (CL)** (Bard et al. 2005, ongoing). The formal ontology of animal cell types. Every cell type has a unique CL identifier (e.g., CL:0000814 for "mature NK T cell"). Hierarchical — parent-child relationships between cell types and lineages.
- **The Human Cell Atlas (HCA)**. An effort to map every cell type in the human body (HCA Consortium, 2017–). As of 2024, ~30 tissue-scale datasets published. HCA data is the reference for most modern annotation pipelines.
- **Tabula Sapiens / Tabula Muris** (Tabula Sapiens Consortium 2022). Single-institution pan-tissue atlases with consistent annotation. Smaller than HCA but more consistently annotated.

Use CL identifiers in your metadata wherever possible. Vocabulary consistency is what makes cross-study comparison work.

> **Historical pointer**: The Human Cell Atlas launched in 2017 with the goal of defining "every cell type in every human tissue." By 2024 the field has ~50M single-cell profiles across hundreds of tissue datasets, plus emerging spatial and multi-omic modalities. The biological reference for what a cell *is* has shifted from textbook morphology-based taxonomy to single-cell-molecular-signature-based taxonomy. For working bioinformatics engineers in 2024+, the reference to build against is HCA, not histology.

---

## Wrap-up (≈10 min)

### What you should take away

- **Single-cell RNA-seq measures per-cell transcript abundance.** Bulk averages hide cell-type composition; single-cell resolves it. The cost is per-cell signal-to-noise; the payoff is heterogeneity becomes visible.
- **UMIs solve PCR duplicate counting via error-correcting tags.** Each mRNA gets a unique serial number before the amplification channel; final counts are distinct UMIs, not reads. Essential at the amplification levels scRNA-seq uses.
- **QC is non-optional: empty droplets, doublets, MT%, ambient RNA.** Every pipeline applies all four. Filtering cuts 20–40% of observed "cells" before downstream analysis.
- **Normalisation must respect sparsity.** Bulk methods (median-of-ratios) fail on mostly-zero matrices. Log-normalisation or Pearson residuals, followed by HVG selection, is the 2024 default.
- **Dimensionality reduction is two-stage: PCA → UMAP.** PCA gives a linear 30–50D representation; UMAP embeds to 2D for visualisation. Clustering runs on PCs, not on UMAP coordinates.
- **Leiden graph clustering is the 2024 default for cell grouping.** Resolution is a free parameter; tune it to biology, not to a number.
- **Annotation: marker-based or reference-based.** Reference-based is fast and automated; marker-based is manual but catches reference gaps. Do both.

### Next lecture

Advanced single-cell: pseudotime and trajectory inference, RNA velocity (a state-space estimator on spliced/unspliced ratios), batch integration (Harmony, scVI as variational EM), multi-modal single-cell (CITE-seq, multi-ome), spatial transcriptomics.

### Homework

1. Download the 10x PBMC 3k dataset (`pbmc_10k_v3` from 10x public datasets). Run a standard Scanpy or Seurat pipeline end-to-end: QC, normalise, HVG, PCA, UMAP, Leiden. Report: cells after QC, number of clusters at resolution 1.0, approximate cell-type assignment for each cluster.
2. Make a knee plot of the PBMC dataset's raw `raw_feature_bc_matrix.h5`. Mark the EmptyDrops threshold. How many barcodes are above the knee? How does that compare to the number of cells 10x reports as valid?
3. Run Scrublet on the same dataset. What's the estimated doublet rate? How does it compare to the 6-8% expected for the loading concentration?
4. Compute the dispersion trend (as in Lecture 6 §3.3) but now on the per-gene variance across single cells. How does it compare to the bulk mean-dispersion trend? Is the empirical-Bayes shrinkage idea still useful here?
5. Annotate the PBMC clusters two ways: (a) SingleR against a reference, (b) manually, using canonical markers. How many clusters do the two methods agree on? Where do they disagree, and which is right?

### Recommended reading

- Macosko, E. Z., Basu, A., Satija, R., et al. (2015). Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. *Cell* 161, 1202–1214. (The Drop-seq paper, the open-source precursor to 10x.)
- Zheng, G. X. Y., Terry, J. M., Belgrader, P., et al. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications* 8, 14049. (The Chromium 10x paper.)
- Lun, A. T. L., Riesenfeld, S., Andrews, T., et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology* 20, 63.
- Hafemeister, C., &amp; Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology* 20, 296. (The SCTransform paper.)
- Traag, V. A., Waltman, L., &amp; van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9, 5233.
- McInnes, L., Healy, J., &amp; Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv:1802.03426*.
- Aran, D., Looney, A. P., Liu, L., et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology* 20, 163–172. (The SingleR paper.)
- The Tabula Sapiens Consortium (2022). The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science* 376, eabl4896.
- Scanpy tutorial: <https://scanpy.readthedocs.io/en/stable/tutorials.html>
- Seurat tutorial: <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Why Single-Cell + Droplet Chemistry          | 45&nbsp;min | 0:45 |
| Part 2 — Data Pipeline and QC                          | 60&nbsp;min | 1:45 |
| Part 3 — Normalisation and Feature Selection           | 25&nbsp;min | 2:10 |
| Part 4 — Dimensionality Reduction and Clustering       | 50&nbsp;min | 3:00 |
| Part 5 — Cell-Type Annotation                          | 20&nbsp;min | 3:20 |
| Wrap-up                                                 | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
