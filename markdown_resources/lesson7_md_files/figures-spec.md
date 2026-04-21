# Lecture 7 — Figures Specification

> **Scope**: Static diagrams for Lecture 7 (Single-Cell RNA-seq Fundamentals).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-07.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. scRNA-seq content is pipeline-/algorithm-heavy; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Line-art first; fills only where they earn the ink.
- Base colours (`--base-a/t/g/c/u`) only where actual nucleotide characters appear (Figure 3's barcode/UMI structure).
- Cell-type and cluster labels in scatter plots use the accent palette variants — three or four distinct colour families per scatter, consistent between Figures 9–11.
- Monospace (JetBrains Mono) for: barcodes, UMIs, gene symbols, UMI counts, numeric values. Inter for everything else.
- Arrows follow the shared `<marker id="arrow-accent">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content.

## Figure Budget

Eleven figures for a ~3h 30min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Single-cell vs bulk — what the average hides | Part 1 | Custom SVG |
| 2 | 10x Genomics droplet chemistry | Part 1 | Custom SVG |
| 3 | Cell barcode + UMI structure and deduplication | Part 1 | Custom SVG |
| 4 | scRNA-seq quantification pipeline | Part 2 | Custom SVG |
| 5 | The knee plot and EmptyDrops threshold | Part 2 | Custom SVG |
| 6 | Doublet detection strategy | Part 2 | Custom SVG |
| 7 | QC metrics dashboard | Part 2 | Custom SVG |
| 8 | Sparse count matrix — zeros everywhere | Part 3 | Custom SVG |
| 9 | PCA on scRNA-seq — scree plot and projection | Part 4 | Custom SVG |
| 10 | UMAP embedding of a multi-cell-type dataset | Part 4 | Custom SVG |
| 11 | Marker gene heatmap | Part 4 | Custom SVG |

---

## Figure 1 — Single-cell vs bulk, what the average hides

**File**: `diagrams/lecture-07/01-bulk-vs-single-cell.svg`
**Lecture anchor**: §1.2 Cell-type heterogeneity and the "average cell" fallacy
**ViewBox**: `0 0 900 440`

### Purpose

Visual core of the lecture's motivation. Show a tissue's cell-type composition, the bulk averaging operation, and the per-cell expression pattern side-by-side. Readers should see why the bulk average hides biologically meaningful structure.

### Content

Three stacked bands:

**Band 1 — A tissue schematic.** A simple cartoon of a tissue slice (rectangular box with ~60 cells drawn as circles, colour-coded into 4 distinct cell types — say, hepatocyte-like, immune-like, endothelial-like, rare progenitor-like). Numbers next to each type: "45%", "30%", "22%", "3%".

**Band 2 — Bulk RNA-seq.** An arrow from the tissue to a single test-tube icon labeled "bulk mRNA extraction". Below it, a flat bar chart of 5 example gene expression levels, each shown as one averaged value.

**Band 3 — Single-cell RNA-seq.** An arrow from the tissue to a "droplet encapsulation" step, producing a cell × gene matrix visualisation. Same 5 genes shown as a heatmap: cells grouped by type, each type showing distinct expression patterns. Annotations point out the rare cell type — visible in the heatmap, invisible in the bulk chart.

Caption: "Bulk RNA-seq loses cell-type composition. Single-cell RNA-seq recovers it."

### Style notes

- Cell-type colours: distinct accent-family hues. Suggest `--accent`, `--success`, `--warning`, `--base-u` (violet) for the 4 types.
- Tissue outline: `--fg-muted` 1.5px stroke.
- Arrows between bands: `--fg-muted` 1.5px with standard marker.
- Heatmap cells: `--accent-bg` to `--accent` gradient for expression intensity.

---

## Figure 2 — 10x Genomics droplet chemistry

**File**: `diagrams/lecture-07/02-droplet-chemistry.svg`
**Lecture anchor**: §1.3 The 10x Genomics droplet chemistry
**ViewBox**: `0 0 900 460`

### Purpose

Show the microfluidic workflow — cells + beads + oil → droplets containing labelled cDNA. Readers should be able to describe the four steps from the figure alone.

### Content

**Left — Chromium chip schematic.** Three input channels converging at a junction: (a) cell suspension (top), (b) gel-bead suspension (middle), (c) oil (bottom). Labelled "Chromium chip — droplet generator". Output channel on the right produces a stream of droplets.

**Right — Single droplet cross-section.** A circular droplet with: a cell inside (labelled, with mRNAs floating in the cytoplasm), a gel bead on the side (with primers extending from its surface), and the oil boundary around everything. Inside the droplet, the reverse-transcription reaction is shown as arrows from mRNAs to primers producing labelled cDNA.

**Bottom band — Four step summary.** 4 panels laid out horizontally showing the four chemistry steps: (1) Load, (2) Lyse + prime, (3) Reverse transcribe, (4) Pool + sequence. Small schematic per step.

### Style notes

- Cell: a filled circle in `--base-a` (red family).
- Gel bead: a filled circle in `--bg-muted` with visible primer extensions in `--accent`.
- Droplet outline: `--fg` 2px.
- Oil: outer region `--bg-inset`.
- Primers: thin `--accent` lines extending from bead, each labeled with (CB + UMI + poly-T) structure tag.

---

## Figure 3 — Cell barcode + UMI structure and deduplication

**File**: `diagrams/lecture-07/03-barcode-umi.svg`
**Lecture anchor**: §1.4 Cell barcodes and UMIs
**ViewBox**: `0 0 960 460`

### Purpose

The central figure tying together UMI biology and the amplification-channel EE framing. Readers should leave able to explain why UMIs give true molecule counts despite PCR bias.

### Content

**Top band — Primer structure on a bead.** A single primer with its three regions labelled: Cell barcode (16 bp), UMI (12 bp), poly-T (30+ T's). Each region shown as a coloured bar with its bp length annotated. A few more primers shown in greyed-out form to indicate "millions of these per bead, all sharing the cell barcode, with distinct UMIs."

**Middle band — Inside a droplet.** 3 mRNA molecules (distinct colours) bind to 3 different primers on the bead. Each primer has the same cell barcode but different UMIs. Arrows show RT extension producing 3 distinct cDNAs: (CB + UMI₁ + mRNA₁), (CB + UMI₂ + mRNA₂), (CB + UMI₃ + mRNA₃).

**Lower-middle band — PCR amplification.** The 3 cDNAs are amplified to ~50 copies each. Visualise as 150 cDNA copies arranged in three coloured groups.

**Bottom band — Deduplication.** After sequencing, 150 reads collapse on (cell barcode, UMI) to recover 3 distinct molecule counts. A small table showing the dedup operation: input 150 reads → output 3 UMI groups → count = 3.

### Style notes

- Cell barcode: `--accent` coloured bar.
- UMI: `--warning` coloured bar (distinct from CB).
- Poly-T: `--fg-muted` coloured bar.
- mRNA molecules: three distinct base colours (`--base-a`, `--base-t`, `--base-g`).
- PCR duplicates: same colour as parent mRNA, reduced opacity.

---

## Figure 4 — scRNA-seq quantification pipeline

**File**: `diagrams/lecture-07/04-quant-pipeline.svg`
**Lecture anchor**: §2.2 Quantification tools
**ViewBox**: `0 0 1100 340`

### Purpose

Horizontal flow showing the pipeline from BCL files to an AnnData/SeuratObject, with the three tool-choice branches (Cell Ranger, STARsolo, alevin-fry) visible.

### Content

Left-to-right stages:

1. **BCL files** (Illumina flow cell icon).
2. **bcl2fastq** → FASTQ files.
3. **Quantification** — a three-branch box: Cell Ranger / STARsolo / alevin-fry (all three listed as parallel options with tool-name tags).
4. **Output** — BAM (for alignment-based) + count matrix (h5 or mtx).
5. **AnnData / Seurat** — cell × gene container, ready for QC.

Each stage labeled with its input format and output format. Each tool branch annotated with runtime characteristics ("fastest", "reference-standard", "pseudoalignment"). Annotation callout: "format choice locks you into specific downstream ecosystems."

### Style notes

- Boxes: `--bg-muted` fill with `--border-strong` outline; rounded corners.
- The quantifier box is the widest — show as a grouped set of three smaller boxes nested inside a parent "Quantification" frame.
- Arrows between stages: `--accent` 1.5px with standard arrow marker.
- Format labels on each arrow: JetBrains Mono 10, `--fg-muted`.

---

## Figure 5 — The knee plot and EmptyDrops threshold

**File**: `diagrams/lecture-07/05-knee-plot.svg`
**Lecture anchor**: §2.3 Empty droplets and the knee plot
**ViewBox**: `0 0 900 440`

### Purpose

The canonical scRNA-seq diagnostic plot. Readers should understand both the shape and the threshold-placement logic.

### Content

Log-log plot with:

- X-axis: barcode rank (sorted descending by UMI count). Range 1 to ~1,000,000.
- Y-axis: UMI count per barcode. Range 10⁰ to 10⁵.
- The knee curve: sharp transition from a "real cells" plateau on the left to a steep drop at ~5,000 barcodes, and a long tail of low-count empty barcodes on the right.
- Marked threshold line at the knee (cell/background boundary).
- Annotations on the plateau: "real cells (~10,000)".
- Annotations in the transition: "ambiguous zone (EmptyDrops test)".
- Annotations in the tail: "empty droplets (ambient only)".

Inset: EmptyDrops' test — a small panel showing the ambient gene profile (distribution across all genes) vs a candidate barcode's gene profile. Caption: "significantly different = real cell".

### Style notes

- Curve: `--accent` 2px solid.
- Threshold line: `--warning` 1.5px dashed.
- Axis ticks + labels in JetBrains Mono 10 `--fg-muted`.
- Inset: `--bg-inset` filled background, `--border-strong` outline.

---

## Figure 6 — Doublet detection strategy

**File**: `diagrams/lecture-07/06-doublet-detection.svg`
**Lecture anchor**: §2.4 Doublet detection
**ViewBox**: `0 0 900 400`

### Purpose

Show the classifier-training idea behind Scrublet / DoubletFinder: simulate doublets as averages of pairs, use them as labeled positives.

### Content

Three panels side-by-side:

**Panel 1 — Two real cells.** Two distinct cell types shown as UMAP-like points, each with a mini expression profile vector next to it.

**Panel 2 — Synthetic doublet.** The two real-cell vectors are combined (averaged). The resulting "doublet" vector is shown. An arrow indicates "synthetic doublet simulation".

**Panel 3 — Classifier decision boundary.** A 2D UMAP embedding showing real cells organized into 4-5 cell-type clusters, with synthetic doublets overlaid — they land in the "between-cluster" gap regions. A decision boundary (curved contour) is drawn separating high-doublet-score regions from low.

Bottom panel: a histogram of doublet scores with a threshold (e.g., 0.3) marked. Cells above threshold are flagged as predicted doublets.

### Style notes

- Real cells: filled circles in cell-type-specific colours (consistent with Figures 10–11).
- Synthetic doublets: open circles with striped pattern, in `--warning` colour.
- Decision boundary: `--accent` 1.5px dashed.
- Histogram bars: `--bg-muted` fill with `--accent` stroke; threshold line in `--warning`.

---

## Figure 7 — QC metrics dashboard

**File**: `diagrams/lecture-07/07-qc-metrics.svg`
**Lecture anchor**: §2.5 Mitochondrial fraction and ambient RNA
**ViewBox**: `0 0 960 460`

### Purpose

Four-panel QC dashboard showing the standard metrics and their filtering. Readers should recognise each plot type and know what threshold decisions it supports.

### Content

2×2 grid of small panels:

**(a) Total UMIs per cell.** Histogram with x-axis log-scale 10²–10⁵, showing the bimodal distribution seen on the knee plot from a different angle. Threshold line marked.

**(b) Genes-per-cell.** Histogram with x-axis 200–8000. Low-gene cells highlighted in `--warning`.

**(c) MT fraction %.** Histogram with x-axis 0–50%. Most cells 5–15%; long tail of high-MT cells in `--error`. Threshold line at 15%.

**(d) Ambient RNA contamination fraction (per cell, post-SoupX).** Histogram showing the distribution of per-cell ambient fractions 0–40%. Some cells high-contamination.

Each panel labelled with metric name and typical threshold.

### Style notes

- All panels share the same axis-label style and bar colours.
- Good-cell regions in `--accent-bg`, bad-cell regions in `--warning-bg` or `--error-bg`.
- Each panel has a small legend: "keep" / "remove" colour key.

---

## Figure 8 — Sparse count matrix, zeros everywhere

**File**: `diagrams/lecture-07/08-sparse-matrix.svg`
**Lecture anchor**: §3.1 Why bulk normalisation fails on sparse data
**ViewBox**: `0 0 900 440`

### Purpose

Show visually that a scRNA-seq count matrix is dominated by zeros. Readers should leave with a strong visual memory that sparsity is a feature of the data, not a choice.

### Content

A heatmap tile: 50 cells (rows) × 200 genes (columns). Cells are organized by cell-type group (shown as row-side colour bars). Genes are organized by category (shown as column headers).

Cell intensity coloured by expression level:

- Zero: pure white (or `--bg-muted` for subtle contrast).
- Low expression: pale `--accent-bg`.
- Medium: `--accent` half-intensity.
- High: full `--accent`.

Visible structure: 2–3 visible blocks where one cell type has distinct non-zero expression in a specific gene group. The rest of the matrix is mostly blank (zero).

Annotations: percentages. "95% of entries are zero" (prominent callout). "Real signal visible as sparse blocks."

### Style notes

- Cells: small squares (8×8 px each, so the full matrix is 400×1600 — but scale to viewBox).
- Row colour bar: 4 cell-type colours, consistent with Figures 1, 6, 10.
- Column headers: in 8px Inter, rotated 90° if needed.
- Non-zero blocks: clear contrast with zero background.

---

## Figure 9 — PCA on scRNA-seq, scree plot and projection

**File**: `diagrams/lecture-07/09-pca-projection.svg`
**Lecture anchor**: §4.2 PCA on the normalised matrix
**ViewBox**: `0 0 900 420`

### Purpose

Show the two canonical PCA outputs — scree plot and 2D projection — that every scRNA-seq workflow produces.

### Content

**Left panel — scree plot.** X-axis: PC rank (1 to 50). Y-axis: fraction of variance explained per PC (log scale).
- Steep drop from PC1 (~8%) to PC5 (~2%).
- Plateau from PC10 onward (~0.5% each).
- A "kept" threshold at PC30 — everything left is kept for downstream.

**Right panel — PC1 × PC2 projection.** Scatter of ~1000 cells. Three visibly distinct clusters visible in PC1/PC2 space (representing 3 cell types). Cells colour-coded by true cell type.

Annotations: scree kneepoint highlighted. Cumulative variance (e.g., 65% through PC30) labelled.

### Style notes

- Scree bars or dots: `--accent`.
- Kept threshold: `--warning` vertical line.
- Scatter clusters: distinct colours matching Figure 10 convention.

---

## Figure 10 — UMAP embedding of a multi-cell-type dataset

**File**: `diagrams/lecture-07/10-umap-embedding.svg`
**Lecture anchor**: §4.3 UMAP and t-SNE for visualisation
**ViewBox**: `0 0 900 480`

### Purpose

The canonical scRNA-seq result visualisation. Readers should recognise this view from any paper.

### Content

A single 2D UMAP scatter with:

- ~8000 cells, each as a 2px circle.
- Colour-coded into 10 Leiden clusters (distinct colour per cluster).
- Cluster labels overlaid in the centre of each cluster.
- Cluster names: immune-system appropriate (e.g., "CD4 T", "CD8 T", "B cells", "NK", "monocytes", "DC", "neutrophils", "platelets", "erythrocytes", "rare X").

Annotations: scale bar, "UMAP1" / "UMAP2" axis labels, legend.

Side inset: a small callout showing a zoomed-in boundary between two related clusters (e.g., CD4 vs CD8 T cells) — showing smooth transition vs sharp separation.

### Style notes

- Cluster colours: 10 distinct accent-family hues. Maintain consistency across any scRNA-seq figure (7, 10, 11) that shows clusters.
- Labels in Inter 11, white-haloed for readability against coloured backgrounds.
- Axes: minimal. Just "UMAP1" / "UMAP2" labels; no axis numbers (they're arbitrary).

---

## Figure 11 — Marker gene heatmap

**File**: `diagrams/lecture-07/11-marker-heatmap.svg`
**Lecture anchor**: §4.5 Marker gene identification
**ViewBox**: `0 0 960 460`

### Purpose

Show the canonical way cluster markers are visualised — a heatmap with cells grouped by cluster and top-5 markers per cluster on the columns.

### Content

- Rows: cells, grouped by cluster (cluster colour bar on left).
- Columns: top 5 marker genes per cluster (50 genes total for 10 clusters).
- Cell colour: expression intensity (same gradient as Figure 8).
- Block-diagonal structure: each cluster's cells show strong expression in their own marker columns, near-zero in other clusters' markers.

Label each cluster's marker block (e.g., "B cell markers: MS4A1, CD79A, CD19, ...").

### Style notes

- Clear cluster boundaries as dark horizontal lines between groups.
- Column group boundaries as thin vertical lines between cluster-specific marker blocks.
- Gene symbols in JetBrains Mono 8, rotated 90°.

---

## Cross-Figure Consistency Notes

- **Cell-type colours**: maintain the same colour for the same cell type across Figures 1, 6, 10, 11.
- **UMAP coordinates**: Figures 6 and 10 should use coordinate-consistent layouts so a reader can trace the same cells.
- **Count-matrix heatmaps**: Figures 8 and 11 should use the same expression-intensity gradient (`--bg-muted` to `--accent`).
- **Scree / histogram axes**: log-scale axes labeled in JetBrains Mono, with consistent tick density.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All eleven figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient (except the count-matrix intensity gradient in Figures 8, 11), drop shadow, glow, or 3D effect.
- [ ] All sequences, barcodes, UMIs, gene symbols, numeric values in JetBrains Mono; all other labels in Inter.
- [ ] Base colours (`--base-*`) only where nucleotides are shown (Figure 3 primer + mRNA structure).
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
- [ ] All `&`, `<`, `>` in text content are XML-escaped.
