# Lecture 3 — Figures Specification

> **Scope**: Static diagrams for Lecture 3 (DNA Sequence Assembly).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-03.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. No photographs in Lecture 3 — the content is algorithmic / conceptual.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Every figure is **line-art first**: stroke-only by default, fill only where the fill earns the ink.
- Base colors (`--base-a/t/g/c`) are used only when actual nucleotide characters are shown. For generic alphabet diagrams (graphs of k-mer structure, coverage plots, metric plots), use the neutral palette.
- Any DNA sequence, k-mer, contig, or graph-node label uses **JetBrains Mono**; all other labels use **Inter**.
- Arrows follow the shared `<marker id="arrow-accent">` pattern from `diagram-style-guide.md` §7.

## Figure Budget

Eleven figures for a ~3h 25min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | De novo vs resequencing | Part 1 | Custom SVG |
| 2 | The assembly pipeline | Part 1 | Custom SVG |
| 3 | Read types for assembly | Part 1 | Custom SVG |
| 4 | Coverage distribution | Part 2 | Custom SVG |
| 5 | k-mer spectrum | Part 2 | Custom SVG |
| 6 | De Bruijn graph | Part 3 | Custom SVG |
| 7 | Overlap graph | Part 3 | Custom SVG |
| 8 | Graph topology signatures | Part 4 | Custom SVG |
| 9 | Repeat collapse | Part 4 | Custom SVG |
| 10 | Scaffolding with mate pairs | Part 4 | Custom SVG |
| 11 | N50 metric | Part 5 | Custom SVG |

---

## Figure 1 — De Novo vs Resequencing

**File**: `diagrams/lecture-03/01-denovo-vs-resequencing.svg`
**Lecture anchor**: §1.2 De novo vs resequencing
**ViewBox**: `0 0 860 360`

### Purpose

Establish the semantic difference between the two problems. The lecture repeatedly returns to the distinction, and readers should be able to recall which stack applies from a single visual.

### Content

Two side-by-side panels separated by a thin vertical rule.

**Left panel — Resequencing.** A long horizontal reference bar, pre-drawn in `--fg-muted`, labeled `REFERENCE`. Short reads drawn as thin accent-coloured bars above the reference, each with a dashed leader line pointing to its aligned position. Kicker above: `RESEQUENCING · READS ← REFERENCE`.

**Right panel — De novo assembly.** The same short reads, but no reference bar. Instead, the reads are partially overlapping each other, with overlap-regions highlighted. A faint dashed horizontal contig line beneath them (labeled `reconstructed contig`) showing the assembly emerging from the overlaps. Kicker above: `DE NOVO · CONTIG ← READS`.

### Style notes

- Reference bar in resequencing panel: `--bg-muted` fill, `--fg-muted` 2.5-weight stroke, labeled with 5' and 3' ends.
- Reads: `--accent` 1.5-weight outline, no fill.
- Overlap highlighting in de novo panel: soft `--accent-bg` fill on the overlapping segments, small `·∼·` tick between the reads.
- Caption at bottom centred: "Same reads, different problem. Resequencing aligns to a known template; de novo assembly reconstructs the template from the reads."

---

## Figure 2 — The Assembly Pipeline

**File**: `diagrams/lecture-03/02-assembly-pipeline.svg`
**Lecture anchor**: §1.3 The typical assembly pipeline
**ViewBox**: `0 0 1100 260`

### Purpose

A canonical six-to-seven-box flow diagram matching the structure of Lecture 2's pipeline figure. Readers come back to this figure to locate which stage any given section is talking about.

### Content

Seven boxes, left-to-right:

1. **Reads** — input, kicker `FASTQ`
2. **QC / Trimming** — `fastp`, `trim_galore`
3. **Error Correction** — `BayesHammer`
4. **Graph Construction** — de Bruijn / overlap
5. **Graph Cleanup** — tips, bubbles, tangles
6. **Contig Construction** — non-branching paths
7. **Scaffolding** — paired / mate-pair / Hi-C
8. *(Optional, drawn as a branch)* **Polishing** — `Pilon`, `Racon`

Tapered cobalt arrows between boxes. Below each box, a thin `--fg-muted` caption with the tools typically used (Inter 10).

### Style notes

- Same box treatment as Lecture 2's pipeline figure: `--bg-muted` fill, `--border-strong` 1.5px stroke, 6px rounded corners.
- The optional "Polishing" box is drawn slightly below the main line with a dashed feeder arrow from "Scaffolding" to indicate it's optional.
- Kicker above the whole flow: `ASSEMBLY PIPELINE · ~30 MIN FOR A BACTERIAL GENOME`.

---

## Figure 3 — Read Types for Assembly

**File**: `diagrams/lecture-03/03-read-types.svg`
**Lecture anchor**: §1.4 Reads and read types for assembly
**ViewBox**: `0 0 860 360`

### Purpose

Compare the four relevant read categories on a common scale — length, accuracy, pair structure. Readers should see at a glance why long reads matter for repeats.

### Content

Four horizontal tracks stacked vertically, all drawn on the same base-pair scale bar along the bottom (0 → 50 kb, with tick marks at 1, 10, 50 kb).

1. **Illumina short-read**: 150 bp, single. A single small tight bar on the left side of the scale.
2. **Illumina paired-end**: 2×150 bp with a 500 bp insert. Two small bars linked by a dashed arc labeled `insert ≈ 500 bp`.
3. **PacBio HiFi**: 15 kb, single, high accuracy. A long bar with a crisp outline.
4. **ONT raw**: 50 kb, single, low accuracy. A very long bar with a dashed/fuzzy outline to imply noise.

Each track has a left-side label in Inter small-caps: `ILLUMINA SHORT`, `ILLUMINA PE`, `PACBIO HIFI`, `ONT RAW`. Right-side: a small annotation with typical error rate in `--fg-muted` (e.g. `≤0.5%`, `0.1%`, `5–15%`).

### Style notes

- Base-pair scale bar: 1px `--border-strong` with tick marks and JetBrains Mono labels.
- Read bars: 2.5-weight `--accent` stroke.
- ONT dashed outline: `stroke-dasharray="4 2"`.
- Insert arc for paired-end: thin `--fg-muted` dashed curve.

---

## Figure 4 — Coverage Distribution

**File**: `diagrams/lecture-03/04-coverage-distribution.svg`
**Lecture anchor**: §2.1 Coverage and average coverage
**ViewBox**: `0 0 860 420`

### Purpose

Visualize coverage both spatially (reads pile up on a genome) and statistically (histogram vs Poisson). Readers should see that 30× mean does not mean uniform 30×.

### Content

**Top half** — spatial. A horizontal reference bar with ~15 reads drawn as horizontal bars at various positions, offset vertically so overlap is visible. Below the reference, a filled-area "coverage depth" plot showing per-position coverage (height oscillates between ~5 and ~15). Annotations point out a low-coverage dip and a high-coverage peak.

**Bottom half** — statistical. A histogram of per-position coverage (bars), overlaid with the theoretical Poisson probability mass at the observed mean. X-axis: depth. Y-axis: fraction of positions. The Poisson curve matches well in the bulk but the histogram has a small spike at depth = 0 (the "dropouts" that the Poisson doesn't predict).

### Style notes

- Reads: thin `--accent` bars.
- Coverage depth plot: filled area, `--bg-muted` fill, `--fg` outline.
- Histogram bars: filled `--accent-bg`, `--accent` outline.
- Poisson curve: dashed `--fg-muted`.
- Annotations in Inter 10, `--fg-muted`.

---

## Figure 5 — k-mer Spectrum

**File**: `diagrams/lecture-03/05-kmer-spectrum.svg`
**Lecture anchor**: §2.3 Error correction
**ViewBox**: `0 0 860 360`

### Purpose

Show that true k-mers and erroneous k-mers form two cleanly separable populations in the frequency histogram. Make the error-correction threshold viscerally obvious.

### Content

A histogram of k-mer frequency. X-axis: frequency (how many times each k-mer appears, from 1 to ~60). Y-axis: count of distinct k-mers with that frequency, **log scale**.

Two populations visible:

- A very tall bar at frequency 1 (labeled `erroneous k-mers` with an annotation pointing at it).
- A bump centered around frequency 30 (labeled `true k-mers · ~30× coverage`).

A dashed vertical line at frequency ≈ 3 labeled `error threshold · discard below`.

An inset (upper right) shows a magnified view of the transition region with Inter annotation explaining "k-mers below threshold are assumed to be single-read errors."

### Style notes

- Log-y axis with tick marks at 10⁰, 10¹, 10², 10³, 10⁴, 10⁵ in JetBrains Mono.
- Histogram bars: `--accent` fill.
- Error-threshold line: `--error` stroke 1.5, dashed (`stroke-dasharray="4 3"`).
- Inset: `--bg-muted` background with thin `--border` stroke.

---

## Figure 6 — De Bruijn Graph

**File**: `diagrams/lecture-03/06-debruijn-graph.svg`
**Lecture anchor**: §3.2 De Bruijn graphs
**ViewBox**: `0 0 860 440`

### Purpose

The canonical construction figure. Given three input reads and k = 4, show the k-mer decomposition and the resulting graph. This is the figure readers will scroll back to during homework.

### Content

**Top band — input reads.** Three reads in JetBrains Mono: `ACGTCC`, `CGTCCA`, `GTCCAT`. Kicker: `INPUT READS`.

**Middle band — k-mers.** Every k-mer (k=4) extracted from each read, laid out horizontally under its parent read. Arrows from each k-mer splitting it into its (k−1)-mer prefix and suffix:

```
ACGT → ACG | CGT
CGTC → CGT | GTC
GTCC → GTC | TCC
TCCA → TCC | CCA
CCAT → CCA | CAT
```

**Bottom band — de Bruijn graph.** Six nodes — `ACG`, `CGT`, `GTC`, `TCC`, `CCA`, `CAT` — drawn as circles, connected by directed edges in a left-to-right chain. Each edge is labeled with the k-mer that produced it.

Annotation bubble: "One node per distinct (k−1)-mer; one edge per distinct k-mer."

### Style notes

- Reads in JetBrains Mono 14, `--fg`.
- k-mer arrows (→) in `--fg-muted`, thin.
- Graph nodes: circles, `--bg-muted` fill, `--fg` stroke 1.5.
- Graph edges: cobalt arrows with labels in JetBrains Mono 10.
- Node labels in JetBrains Mono 12, weight 500.

---

## Figure 7 — Overlap Graph

**File**: `diagrams/lecture-03/07-overlap-graph.svg`
**Lecture anchor**: §3.4 Overlap graphs
**ViewBox**: `0 0 860 420`

### Purpose

Contrast with Figure 6. Show that in an overlap graph each read is one node and edges encode suffix-prefix overlaps. Highlights that this graph is smaller in node count but harder to search.

### Content

Five reads listed on the left side in JetBrains Mono:
```
R1: ACGTCC
R2: CGTCCA
R3: GTCCAT
R4: CCATAG
R5: CATAGC
```

Middle panel: a pairwise overlap table showing which read's suffix matches which read's prefix, with overlap length highlighted.

Right panel: the overlap graph — 5 nodes labeled R1..R5, directed edges labeled with the overlap length (e.g. "5 bp", "4 bp"). Graph is almost linear — R1 → R2 → R3 → R4 → R5 — which is the ideal (no branching) case.

### Style notes

- Read labels JetBrains Mono 12.
- Overlap table: `--bg-muted` background, cells with numbers in `--fg`; matching cells highlighted `--accent-bg`.
- Graph nodes: circles, `--bg-muted` fill. Edges: thin cobalt with overlap-length labels in `--fg-muted` Inter 10.

---

## Figure 8 — Graph Topology Signatures

**File**: `diagrams/lecture-03/08-topology-signatures.svg`
**Lecture anchor**: §4.2 Graph topology: tips, bubbles, and tangles
**ViewBox**: `0 0 860 320`

### Purpose

Three canonical patterns side by side. Each reader should come away able to sketch the three from memory.

### Content

Three equal panels:

**Panel 1 — Tip.** A linear path with a short dead-end branch attached at one node (2–3 edges long, then terminates). Kicker `TIP`, caption `single-read error · remove if shorter than 2k`.

**Panel 2 — Bubble.** Two branch nodes connected by two parallel paths of roughly equal length. Kicker `BUBBLE`, caption `SNP or small variant · keep higher-coverage path, collapse`.

**Panel 3 — Tangle.** A dense subgraph with three or four in-edges and three or four out-edges, multiple internal cycles. Kicker `TANGLE`, caption `repeat · unresolvable without long reads`.

### Style notes

- All three panels use the same graph vocabulary: nodes as small circles, edges as thin cobalt arrows.
- Each kicker above its panel in small-caps Inter 11 weight 600 `--accent`.
- Captions below in Inter 10 `--fg-muted`.

---

## Figure 9 — Repeat Collapse

**File**: `diagrams/lecture-03/09-repeat-collapse.svg`
**Lecture anchor**: §4.3 Repeats and ambiguity
**ViewBox**: `0 0 860 360`

### Purpose

Explain visually why a repeat collapses in a de Bruijn graph. The reader should see the two copies merging into one node and the resulting ambiguity on the upstream/downstream connections.

### Content

**Top panel — genome layout.** A horizontal genome bar with two identical repeats `R` drawn as shaded blocks (soft `--accent-bg` fill) at positions X and Y, with distinct flanking regions labeled `A` (left of first copy), `B` (between the two copies), `C` (right of second copy).

**Bottom panel — de Bruijn graph.** The two R copies collapse into a single graph node/path labeled `R`. Coming into R: two in-edges from regions `A` and `B`. Coming out of R: two out-edges to regions `B` and `C`.

Annotation: "Two possible orderings consistent with the graph: A-R-B-R-C and A-R-C-B-R."

### Style notes

- Genome bar as in Figure 1.
- Repeat blocks: `--accent-bg` fill, `--accent` stroke.
- Graph nodes: circles `--bg-muted` fill.
- Ambiguous paths labeled with a `?` marker in `--warning`.

---

## Figure 10 — Scaffolding with Mate Pairs

**File**: `diagrams/lecture-03/10-scaffolding-mate-pairs.svg`
**Lecture anchor**: §4.4 Gaps and scaffolding
**ViewBox**: `0 0 980 320`

### Purpose

Show how paired-end / mate-pair reads bridge unresolved gaps between contigs to produce a scaffold.

### Content

**Top band — contigs.** Three contigs drawn as horizontal bars labeled `CONTIG 1`, `CONTIG 2`, `CONTIG 3`. Gaps between them shown as zig-zag breaks.

**Middle band — mate-pair reads.** Several arcs connecting one point on `CONTIG 1` to a point on `CONTIG 2` (labeled `2 kb insert`), similar arcs connecting `CONTIG 2` to `CONTIG 3` (labeled `10 kb insert`).

**Bottom band — resulting scaffold.** A single scaffold bar `SCAFFOLD` showing contigs in order with `N×??` placeholder regions between them (where the gap length is estimated from the linking-mate insert size).

### Style notes

- Contig bars: `--bg-muted` fill, `--fg` 2.5-weight stroke.
- Gaps between contigs drawn as zig-zag `∿∿∿` in `--fg-subtle`.
- Mate-pair arcs: thin cobalt dashed curves.
- Scaffold N-placeholder regions: shown as `N×` labels inside a slightly-differently-shaded bar section.

---

## Figure 11 — N50 Metric

**File**: `diagrams/lecture-03/11-n50-metric.svg`
**Lecture anchor**: §5.2 Assembly metrics
**ViewBox**: `0 0 860 420`

### Purpose

Make N50 concrete with a small worked example. The reader should be able to compute N50 by hand after seeing the figure.

### Content

**Top half — sorted contig bars.** Horizontal bars representing contigs, sorted by length longest-first, with length labels. Example list: [4500, 3200, 2800, 2100, 1500, 900, 700, 600, 400, 300]. Total 17,000. Half = 8,500.

**Bottom half — cumulative-length curve.** X-axis: contig index (1..10). Y-axis: cumulative total length. The curve rises sharply first (large contigs), then levels off. A horizontal reference line at 50% (= 8,500) intersects the curve at the third contig. A vertical drop line from that intersection points to the corresponding contig in the top half. The length of that contig (2,800) is labeled `N50 = 2,800`.

A second horizontal reference line at 50% of an assumed genome size (say 20,000, → 10,000) shows how NG50 differs — it intersects the curve at a later contig, giving a smaller `NG50 = 2,100`.

### Style notes

- Contig bars: `--accent` 2.5-weight stroke, no fill, with length labels in JetBrains Mono 10.
- Cumulative curve: thick cobalt line.
- N50 reference line: dashed `--accent` at y = 50%.
- NG50 reference line: dashed `--fg-muted` at y = 50% of genome.
- Intersection labeled with a cobalt dot and a JetBrains Mono value.

---

## Cross-Figure Consistency Notes

- **Reference bars and reads** look the same across Figures 1 and 4: `--bg-muted` filled reference, `--accent` outlined reads, consistent with the Lecture 2 convention.
- **Graph nodes and edges** look the same across Figures 6, 7, 8, 9: circles with `--bg-muted` fill and `--fg` stroke, cobalt directional edges.
- **Kicker labels** (uppercase, small-caps Inter 11 weight 600) are used for section headers inside every figure.
- **Base colors** are used only where actual nucleotide characters are shown (Figures 6, 7). Graph-structure figures (8, 9, 10) use only the neutral palette.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All eleven figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient, drop shadow, glow, or 3D effect.
- [ ] All sequences and graph labels are in JetBrains Mono; all other labels in Inter.
- [ ] Base colors appear only where actual nucleotides are shown.
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
