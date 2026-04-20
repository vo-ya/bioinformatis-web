# Lecture 5 — Figures Specification

> **Scope**: Static diagrams for Lecture 5 (Bulk RNA-seq).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-05.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. RNA-seq content is algorithmic / schematic; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Line-art first; fill only where the fill earns the ink.
- Base colors (`--base-a/t/g/c/u`) only where actual nucleotide characters are shown. RNA figures use `--base-u` (violet) where the context is mRNA.
- Genomic sequences, k-mers, and transcript IDs use **JetBrains Mono**; all other labels use **Inter**.
- Arrows follow the shared `<marker id="arrow-accent">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content.

## Figure Budget

Eleven figures for a ~3h 30min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Transcription and splicing — gene to isoforms | Part 1 | Custom SVG |
| 2 | RNA-seq pipeline at a glance | Part 1 | Custom SVG |
| 3 | Reads spanning an exon–intron–exon junction | Part 2 | Custom SVG |
| 4 | STAR seed-cluster-extend with spliced extension | Part 2 | Custom SVG |
| 5 | HISAT2's hierarchical graph FM-index | Part 2 | Custom SVG |
| 6 | Compatibility classes — transcripts, reads, k-mers | Part 3 | Custom SVG |
| 7 | Target de Bruijn graph and pseudoalignment walk | Part 3 | Custom SVG |
| 8 | EM iterations for read-to-transcript assignment | Part 3 | Custom SVG |
| 9 | TPM vs FPKM vs CPM for a small example | Part 4 | Custom SVG |
| 10 | DESeq2 size factors — median of ratios | Part 4 | Custom SVG |
| 11 | Poisson vs Negative Binomial for RNA-seq counts | Part 4 | Custom SVG |

---

## Figure 1 — Transcription and splicing, gene to isoforms

**File**: `diagrams/lecture-05/01-splicing-isoforms.svg`
**Lecture anchor**: §1.2 Transcription, splicing, and isoforms
**ViewBox**: `0 0 900 420`

### Purpose

Establish, in one image, the three biological facts that drive every algorithm in the rest of the lecture: transcription, splicing, and alternative isoforms. Readers should come away knowing why a mature mRNA is shorter than its genomic locus and why the same gene can produce multiple transcripts.

### Content

Three horizontal bands stacked vertically.

**Band 1 — Genomic locus.** A horizontal bar representing a gene's genomic region (~40 kb), with four labeled segments: five exons (filled `--accent-bg`) and four introns (unfilled). Tick labels in JetBrains Mono showing approximate coordinates (`chr7:1..40k`).

**Band 2 — Transcription arrow.** A thin arrow labeled "transcription" connects the genomic locus to a pre-mRNA bar of the same structure (exons + introns preserved). This is the nascent transcript.

**Band 3 — Splicing and isoforms.** A "splicing" arrow from the pre-mRNA produces two alternative mature mRNAs:
- Isoform A: contains all 5 exons, introns removed.
- Isoform B: exon 3 is skipped (alternative splicing); contains exons 1, 2, 4, 5.

Show splice-site motifs GT and AG at two intron boundaries. Caption under each isoform: "Isoform A · 2.1 kb" and "Isoform B · 1.6 kb".

### Style notes

- Exons drawn as filled rectangles `--accent-bg` fill with `--accent` outline.
- Introns drawn as thin lines `--fg-muted` with no fill.
- Transcription + splicing arrows: `--fg-muted` 1.5-weight with arrow markers.
- Splice motifs GT / AG in JetBrains Mono 10, `--accent`.

---

## Figure 2 — RNA-seq pipeline at a glance

**File**: `diagrams/lecture-05/02-rnaseq-pipeline.svg`
**Lecture anchor**: §1.4 The RNA-seq pipeline at a glance
**ViewBox**: `0 0 1100 420`

### Purpose

Road map for the whole lecture. Every subsequent section corresponds to one stage. Readers come back here to locate which stage any given topic is about.

### Content

Seven boxes laid out in two parallel paths, merging at the end:

```
  [FASTQ]
    ↓
  [QC + trim]
    ↓
  ┌─────────────────┬─────────────────┐
  │ Align-then-count│ Pseudoalign     │
  │ STAR, HISAT2    │ Kallisto, Salmon│
  │ (§2)            │ (§3)            │
  │      ↓          │      ↓          │
  │    [BAM]        │ [abundances]    │
  │      ↓          │      ↓          │
  │ [feature counts]│ [tximport]      │
  │      ↓          │      ↓          │
  └──────┬──────────┴──────┬──────────┘
         └── [gene counts] ──┘
                  ↓
          [norm + DE]  (Lecture 6)
```

Each box: kicker line above, tool names below in JetBrains Mono 10 `--fg-muted`. Section cross-references (§2, §3, §4) visible as small tags.

### Style notes

- Boxes: `--bg-muted` fill, `--border-strong` 1.5px stroke, 6px rounded corners.
- Align-then-count path gets `--accent-bg` tint.
- Pseudoalign path gets `--accent-bg` tint of a slightly different shade (still within palette).
- Arrows `--accent` 1.5-weight with the shared marker pattern.

---

## Figure 3 — Reads spanning an exon–intron–exon junction

**File**: `diagrams/lecture-05/03-junction-reads.svg`
**Lecture anchor**: §2.1 The exon-intron structure problem
**ViewBox**: `0 0 900 380`

### Purpose

Make the "jumping gap" problem viscerally concrete. Readers should see exactly why a 150 bp read can't be aligned contiguously to a 40 kb genomic locus.

### Content

Top: a genomic region with three exons (labeled E1, E2, E3) separated by two introns. Exons drawn as filled blocks; introns as thin lines with loop motifs suggesting "5 kb intron here". Scale: genome span shown as ~12 kb.

Middle: six short reads aligned to the genome, each colour-coded by alignment difficulty:

- Reads r1, r2 (green) — fully inside E1. Easy: continuous alignment.
- Reads r3, r4 (amber) — span E1–E2 junction. Their left half lands in E1, right half in E2, with the intron gap bridged in the middle. Show the reads as *split across the intron* with a dashed arc indicating the splice jump.
- Read r5 (amber) — spans E2–E3 junction similarly.
- Read r6 (red) — spans both E1–E2 and E2–E3 junctions. Very challenging; multi-junction read.

Annotations: number of reads falling in each category; "This is what the aligner has to sort out" label.

### Style notes

- Exons: `--accent-bg` filled with `--accent` outline.
- Intron loops: `--fg-muted` 1px dashed.
- Easy reads (r1, r2): `--success` thin stroke.
- Junction reads (r3, r4, r5): `--warning` thin stroke.
- Multi-junction read (r6): `--error` thin stroke.
- Splice jumps (dashed arcs): `--accent` 1px dashed.

---

## Figure 4 — STAR seed-cluster-extend with spliced extension

**File**: `diagrams/lecture-05/04-star-mechanics.svg`
**Lecture anchor**: §2.2 STAR and the spliced-read HMM
**ViewBox**: `0 0 900 460`

### Purpose

Visualise STAR's three-phase algorithm so the reader can narrate it from the figure alone.

### Content

Four stacked bands showing a single read's journey through STAR:

**Band 1 — Input.** A 150-bp read drawn in JetBrains Mono. Label: "Read · 150 bp".

**Band 2 — Seed search (MMPs).** The read below, with three colour-coded MMP segments highlighted: MMP1 (30 bp, genomic coord ~10500), MMP2 (25 bp, genomic coord ~15800), MMP3 (40 bp, genomic coord ~16200). Each MMP drawn as a highlighted sub-range of the read, with an arrow pointing to its genomic location on a reference strip below.

**Band 3 — Seed clustering.** The three MMPs on a genomic coordinate axis. MMP2 and MMP3 cluster (within ~500 bp of each other); MMP1 is far away (~5 kb). The cluster {MMP2, MMP3} is circled; MMP1 is discarded or tried separately.

**Band 4 — Spliced extension.** The clustered MMPs with the read extending across the gap at canonical GT-AG motifs. The extension jumps the intron (shown as a grey zig-zag region) at zero cost because the motif is canonical. Final alignment visualised as a continuous match across two exons with a splice break in the middle.

Side annotation: "Splice jump is FREE at GT-AG motifs; penalised elsewhere."

### Style notes

- MMPs coloured distinctly: MMP1 `--warning`, MMP2 `--accent`, MMP3 `--accent-bright`.
- Clustering circle: `--accent` dashed 1.2px.
- Splice motifs GT / AG in JetBrains Mono 11 weight 600, `--accent`.
- Intron zig-zag: `--fg-subtle` 1px.

---

## Figure 5 — HISAT2's hierarchical graph FM-index

**File**: `diagrams/lecture-05/05-hisat2-hgfm.svg`
**Lecture anchor**: §2.3 HISAT2 and the hierarchical graph FM-index
**ViewBox**: `0 0 920 420`

### Purpose

Show HISAT2's two-layer index design and the graph-edge extension. Readers should see why HISAT2's memory is 4× smaller than STAR's.

### Content

Three stacked bands:

**Band 1 — Linear reference genome.** A long horizontal bar labeled "reference · 3.1 Gb" (for human). Small tick marks along it representing genomic coordinates.

**Band 2 — Global + local FM-indexes.** Below the reference, a wide "global index" box labeled "~8 GB FM-index over the whole reference". Below that, a series of smaller "local index" boxes, each labeled "~57 kb region, ~300 KB index". A visual callout shows how a read search starts in the global index and hands off to a local index for the last few bases.

**Band 3 — Graph overlay.** On top of the reference, draw a set of alternate-path edges for known splice junctions: arcs connecting the end of one exon to the start of a non-adjacent exon. Each arc is labeled "annotated splice". The FM-index is shown as able to traverse these edges during search — reads that match across a splice junction follow the arc.

Side comparison panel: "STAR: ~30 GB suffix array. HISAT2: ~8 GB HGFM. Cost: ~4× speedup on index load."

### Style notes

- Global index: filled `--bg-muted` rectangle with `--border-strong` outline.
- Local index: small filled `--bg-inset` rectangles below.
- Splice arcs: `--accent` 1.5px dashed with small arrowhead at the far end.
- Comparison panel: enclosed in a thin `--accent` bordered box.

---

## Figure 6 — Compatibility classes: transcripts, reads, k-mers

**File**: `diagrams/lecture-05/06-compatibility-classes.svg`
**Lecture anchor**: §3.1 The compatibility-class insight
**ViewBox**: `0 0 900 440`

### Purpose

Establish the central abstraction of pseudoalignment. Reader should see exactly how a read's k-mers define the set of transcripts it's compatible with, and should understand why this information is sufficient for quantification.

### Content

**Top band — three transcripts.** T1, T2, T3 as labeled horizontal bars (~300 bp each). Overlap regions visible: T1 and T2 share a middle segment; T2 and T3 share their right tail.

**Middle band — k-mer tables.** For each transcript, a small strip of k-mer tags (k = 11 here for legibility). Colour-code each k-mer by its transcript-set:
- Solid green if only in T1.
- Solid blue if only in T2.
- Blue-green diagonal if in both T1 and T2.
- Red-blue diagonal if in both T2 and T3.
- Black stripes if in all three.

**Bottom band — four example reads and their compatibility classes.** Each read drawn with its k-mer colour pattern:
- Read A has all-green k-mers → compatibility class {T1}
- Read B has all-(blue-red) k-mers → compatibility class {T2, T3}
- Read C has all-black k-mers → {T1, T2, T3}
- Read D has a k-mer not in any transcript → {} (drop)

### Style notes

- Transcripts drawn as filled `--accent-bg` bars with `--accent` outline.
- k-mer tags in JetBrains Mono 9 with coloured backgrounds per transcript-set.
- Reads shown with per-k-mer coloured cells.
- Compatibility class labels in Inter 11 weight 600, colour-coded to match.

---

## Figure 7 — Target de Bruijn graph and pseudoalignment walk

**File**: `diagrams/lecture-05/07-tdbg-kallisto.svg`
**Lecture anchor**: §3.2 Kallisto and the Target de Bruijn Graph
**ViewBox**: `0 0 920 440`

### Purpose

Connect the T-DBG concept to the de Bruijn graph from Lecture 3, and show how a read's walk through the graph produces the compatibility class.

### Content

**Top band — three small transcripts.** T1, T2, T3 as short sequences (~30 bp each) in JetBrains Mono. Two transcripts share a middle segment.

**Middle band — Target de Bruijn graph.** A de Bruijn graph constructed from the three transcripts at k = 5. Nodes are 4-mers; edges are 5-mers. Each node is tagged with the *set of transcripts* whose sequence contains it (the tag is a small coloured square: solid = single transcript, striped = multiple transcripts).

**Bottom band — a read's walk.** A single read (say 15 bp) decomposed into its k-mers. The k-mers are traced through the graph as a path; at each step, the transcript-set tag is read; the intersection of all visited tags is the read's compatibility class.

Annotation box: "Pseudoalignment = walk + intersect. No base-level alignment needed."

### Style notes

- Transcripts in JetBrains Mono 12 with `--base-*` per-base colouring (only if displaying raw RNA; otherwise neutral).
- DBG nodes as small circles `--bg-muted` fill, `--fg` outline.
- Edges `--fg-muted` 1.2px.
- Transcript-set tags as small coloured squares beside each node.
- Read walk highlighted with a thick `--accent` edge colouring.

---

## Figure 8 — EM iterations for read-to-transcript assignment

**File**: `diagrams/lecture-05/08-em-iterations.svg`
**Lecture anchor**: §3.4 Expectation-Maximization for read-to-transcript assignment
**ViewBox**: `0 0 1020 480`

### Purpose

Visualise EM converging on a toy example. Readers should be able to track both the fractional-vote matrix (reads × transcripts) and the abundance estimate (θ) across iterations, and should see them stabilise.

### Content

A wide figure showing EM across 5 iterations, arranged as columns.

**Header row:** T1, T2, T3 column headers.

**Left column — Inputs.** Four reads r1..r4 with their compatibility classes listed, and the transcript lengths ℓ = (500, 1500, 2000).

**Each of 5 iteration columns:**
- Small matrix of fractional votes zᵣ,ₜ (reads as rows, transcripts as columns, cells shaded `--accent-bg` with intensity proportional to the vote).
- Below the matrix, the summed vote per transcript, divided by length, renormalised → θₜ⁽ᵏ⁾.
- A small bar chart of θ for that iteration.

**Right column — Converged.** Final θ estimate; the likelihood value; and a note like "converged at iter 5 (Δθ < 10⁻⁶)".

### Style notes

- Fractional-vote cells shaded `--accent-bg` with opacity proportional to vote value.
- θ bars in `--accent` outline.
- Iteration arrows `--fg-muted` connecting successive columns.
- Title kicker: "EM · read-to-transcript assignment · toy example".

---

## Figure 9 — TPM vs FPKM vs CPM for a small example

**File**: `diagrams/lecture-05/09-normalization-units.svg`
**Lecture anchor**: §4.2 CPM, TPM, FPKM — what each normalises for
**ViewBox**: `0 0 900 440`

### Purpose

Make the difference between three normalisation units concrete. The reader should see that the three units rank genes differently depending on transcript length and sample composition.

### Content

A table-style figure with six genes in rows and four sample columns: raw counts, CPM, FPKM, TPM. Each row coloured by the gene's length (short, medium, long). Highlight cells where the rank order between samples differs between units.

Small footnote panel: "CPM corrects for depth only. FPKM corrects for depth and length, but doesn't sum to 10⁶. TPM corrects for length first, then normalises to a per-sample sum of 10⁶."

### Style notes

- Gene length shown as a small coloured tag beside the gene name: green (short), amber (medium), red (long).
- Values in JetBrains Mono 10.
- Cells where the rank differs between units highlighted with a thin `--accent` border.
- Footnote in `--bg-inset` tinted panel.

---

## Figure 10 — DESeq2 size factors, median of ratios

**File**: `diagrams/lecture-05/10-size-factors.svg`
**Lecture anchor**: §4.3 Size factors and library composition
**ViewBox**: `0 0 900 440`

### Purpose

Walk through the median-of-ratios estimator step by step so the reader can implement it by hand after reading the figure.

### Content

Four stacked panels:

**Panel 1 — Raw counts.** Four samples × eight genes count matrix. Show totals per sample (varying depths).

**Panel 2 — Geometric mean per gene.** For each gene, compute the geometric mean of its counts across samples. Show a small column of geometric means.

**Panel 3 — Per-sample ratios.** For each sample, divide each gene's count by the gene's geometric mean. Show a new matrix of ratios. Highlight one sample in detail.

**Panel 4 — Size factors = median of ratios.** For each sample, take the median of its ratio column. Show the four size factors as a small bar chart.

**Panel 5 — Normalised counts.** Divide original counts by size factors. Show the result. A small annotation: "Gene rankings now reflect abundance, not depth or composition."

### Style notes

- Matrix cells: JetBrains Mono 10, right-aligned.
- Ratios cells: colour-coded (pale blue for ratio < 1, pale red for > 1).
- Median-per-sample cell: `--accent-bg` fill, `--accent` border.
- Arrows connecting panels: `--fg-muted` 1.2px.

---

## Figure 11 — Poisson vs Negative Binomial for RNA-seq counts

**File**: `diagrams/lecture-05/11-poisson-vs-nb.svg`
**Lecture anchor**: §4.4 The count distribution and why Poisson isn't enough
**ViewBox**: `0 0 900 440`

### Purpose

Show both the *shape* difference between Poisson and NB distributions at the same mean, and the *empirical* evidence that real RNA-seq data needs NB (the mean-variance plot).

### Content

**Top panel — Two discrete distributions side by side**, both at mean μ = 50.
- Left subplot: Poisson(50). Tight; variance = 50; visible as a bell-shaped concentration of probability.
- Right subplot: NegBin(50, α = 0.3). Visibly wider; variance = 50 + 0.3·50² = 800. Right tail much heavier than Poisson.

**Bottom panel — Empirical mean-variance plot.** X-axis: mean count per gene (log scale, 1 to 10000). Y-axis: variance (log scale, 1 to 10⁷). Scatter of ~2000 gene points. Overlay the Poisson line var = mean (slope 1 on log-log). Points systematically sit *above* this line, visibly curving upward — the hallmark of overdispersion. Annotated arrow: "at μ ≈ 100, typical variance is ~500 — 5× Poisson."

### Style notes

- Poisson bars: `--accent-bg` fill, `--accent` outline.
- NB bars: same fill/outline but slightly different — e.g. `--accent-bright` outline to distinguish.
- Poisson reference line: `--fg-muted` 1.5px dashed.
- Scatter: `--fg` 1px circles, r=2.
- Annotation arrow: `--warning` 1.2px.

---

## Cross-Figure Consistency Notes

- **Transcripts** are drawn as filled `--accent-bg` bars across Figures 1, 6, 7, 9. Use the same aspect ratio and outline style.
- **k-mers and compatibility classes** in Figures 6 and 7 should use consistent colouring: same green-blue-red palette for transcript-set tags.
- **Genomic loci** in Figures 1, 3, 5 should share the exon-filled / intron-line convention.
- **Bars and count data** in Figures 9, 10 should use JetBrains Mono for values and Inter for headers.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All eleven figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient, drop shadow, glow, or 3D effect.
- [ ] All sequences, k-mers, counts, and abundance values are in JetBrains Mono; all other labels in Inter.
- [ ] Base colors appear only where actual nucleotides are shown (primarily Figure 1 for splicing motifs GT/AG).
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
- [ ] All `&`, `<`, `>` in text content are XML-escaped.
