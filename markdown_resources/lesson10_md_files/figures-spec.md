# Lecture 10 — Figures Specification

> **Scope**: Static diagrams for Lecture 10 (Methylation, Hi-C, and 3D Genome Organisation).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-10.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG; content is chemistry-, matrix-, and linear-algebra-heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Base-palette colours (`--base-a/t/g/c`) only where literal nucleotide letters appear (Figures 1, 3, 4).
- Matrices rendered as pixel grids with a consistent accent/amber heatmap; never rainbow.
- Monospace (JetBrains Mono) for: DNA letters, eigenvector coordinates, numerical labels. Inter for UI labels. Source Serif for math / equations where used.
- Arrows follow the shared `<marker id="arrow-accent">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content.

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | DNA methylation — 5mC on a CpG | Part 1 | Custom SVG |
| 2 | CpG islands across a promoter | Part 1 | Custom SVG |
| 3 | Bisulfite conversion chemistry | Part 2 | Custom SVG |
| 4 | Bismark three-letter alignment | Part 2 | Custom SVG |
| 5 | Long-read direct methylation calling | Part 2 | Custom SVG |
| 6 | Beta distributions for methylation proportions | Part 3 | Custom SVG |
| 7 | 3C → Hi-C family tree | Part 4 | Custom SVG |
| 8 | Hi-C library protocol | Part 4 | Custom SVG |
| 9 | Raw vs normalised contact matrix | Part 4 | Custom SVG |
| 10 | TADs and A/B compartments | Part 5 | Custom SVG |
| 11 | A/B compartments as first eigenvector | Part 6 | Custom SVG |
| 12 | TAD calling as change-point detection | Part 6 | Custom SVG |

---

## Figure 1 — DNA methylation, 5mC on a CpG

**File**: `diagrams/lecture-10/01-methylation-chemistry.svg`
**Lecture anchor**: §1.1 What methylation is
**ViewBox**: `0 0 960 420`

### Purpose

Introduce 5-methylcytosine as the modification and CpG as the context. Reader leaves with a concrete picture of the covalent change.

### Content

Three stacked bands:

**Top — DNA helix detail.** A short stretch (10–15 bp) of DNA double helix rendered schematically. Highlight one CpG dinucleotide position with a callout box. Base letters rendered in base-palette colours (A red, C blue, G amber, T green). The cytosine in the highlighted CpG has a small methyl group drawn protruding from its 5-carbon (a small `CH₃` label with a bond line).

**Middle — before and after.** Two adjacent panels:
- Left: `5' ...ACG**C**GTT... 3'` with a filled circle above the C (methylated).
- Right: `5' ...ACG**C**GTT... 3'` with an open circle above the C (unmethylated).
- Annotate "same sequence, two epigenetic states".

**Bottom — structural formulas.** Skeletal formulas of cytosine and 5-methylcytosine side by side. Atom labels (N, C) in standard chemistry convention. Methyl group in accent colour to draw the eye. Arrow between them labelled "DNMT adds" and "TET oxidises" for forward and reverse reactions.

### Style notes

- Helix strands: `--fg` 1.2px.
- Base letters: base-palette per nucleotide.
- Methyl group: `--accent` with thicker stroke (2px) for emphasis.
- Callout boxes: `--border-strong` 0.8px dashed.

---

## Figure 2 — CpG islands across a promoter

**File**: `diagrams/lecture-10/02-cpg-islands.svg`
**Lecture anchor**: §1.2 CpG islands and the methylation landscape
**ViewBox**: `0 0 1080 360`

### Purpose

Show the genome-scale methylation pattern: sparse methylated CpGs across most of the genome, with a dense unmethylated CpG island at the promoter. Reader leaves recognising the pattern and its regulatory implication.

### Content

**Top band — gene model.** A ~10 kb genomic region. Gene body with exons (coloured rectangles) and introns (lines), TSS arrow at the left end. Annotate gene name (e.g. `GENE_X`).

**Middle band — CpG density.** Per-base CpG density plot. High peak at the promoter (inside the CpG island box); sparse low-level CpGs everywhere else. Shaded "CpG island" box overlaying the promoter region.

**Bottom band — methylation status per CpG.** One coloured circle per CpG site, positioned at its genomic coordinate. Filled = methylated; open = unmethylated. Pattern:
- Inside the CpG island: ~30 CpGs, nearly all open (unmethylated).
- Outside: ~15 sparse CpGs, nearly all filled (methylated).
- One intermediate region showing partial methylation (some open, some filled).

Annotations on the right: "island: unmethylated → active promoter", "bulk genome: methylated → structural/silenced".

### Style notes

- Gene exons: `--accent-bg` filled, `--accent` outline.
- CpG island box: `--fg-subtle` dashed outline, semi-transparent `--accent-bg` fill.
- Filled circles (methylated): `--fg` solid.
- Open circles (unmethylated): `--fg` outline only.
- CpG density curve: `--accent` fill under, 1.5px stroke.

---

## Figure 3 — Bisulfite conversion chemistry

**File**: `diagrams/lecture-10/03-bisulfite-chemistry.svg`
**Lecture anchor**: §2.1 Bisulfite conversion chemistry
**ViewBox**: `0 0 1080 420`

### Purpose

Core chemistry reference. Reader leaves able to state: before, middle, after — and which bases change when.

### Content

Four stacked horizontal strands, same DNA sequence (~30 bp) shown each time.

**(a) Reference.** `5' AATCGTCGAAGCGTCGAAGCGTCGAA 3'` with CpGs highlighted. Each C in a CpG gets a colour annotation showing its methylation state: half methylated (filled) and half not (open) for a mix.

**(b) After bisulfite treatment.** Same sequence but:
- Unmethylated Cs → U (highlighted in warning colour).
- Methylated Cs → C (unchanged).
- Non-CpG Cs → U (unmethylated non-CpG Cs converted too, though they are fewer in mammalian genomes).

**(c) After PCR amplification.** U → T. So:
- Unmethylated Cs → T (warning highlighting).
- Methylated Cs → C.

**(d) Summary arrow box.** A mapping diagram: `unmethylated C → T` (warning), `methylated C → C` (accent). Caption: "every surviving C in the read = a methylated C in the original".

### Style notes

- Strand background: `--bg-muted` thin line.
- Letters: JetBrains Mono 12.
- Changed bases: `--warning` background highlight.
- Unchanged methylated C: `--accent` background.
- Arrow markers between stages: `--accent` 1.5px.

---

## Figure 4 — Bismark three-letter alignment

**File**: `diagrams/lecture-10/04-bismark-alignment.svg`
**Lecture anchor**: §2.3 Three-letter alignment
**ViewBox**: `0 0 1080 460`

### Purpose

Show the key algorithmic trick — convert both reference and read into three-letter space, align, then decode back. Reader leaves with a procedural mental model of the algorithm.

### Content

Three rows (algorithm stages), each showing a reference/read pair.

**Row 1 — original.** Reference (top): `ATCGTCG...` with Cs visible. Bisulfite read (bottom): `ATTGTTG...` with Ts where unmethylated Cs were. Label: "naive aligner sees 3 mismatches — fails".

**Row 2 — C→T converted both.** Reference: `ATTGTTG...` (every C turned to T, in silico). Read: `ATTGTTG...` (already has Ts). Label: "three-letter alignment succeeds — exact match".

**Row 3 — decode back.** The read is placed at the matched position in original coordinates. At each original-C position, check what the read has:
- Read has C → methylated (mark with filled circle).
- Read has T → unmethylated (mark with open circle).
Output: per-CpG methylation call.

### Style notes

- Row boxes: `--bg-muted` fill with `--border` outline.
- Original Cs highlighted in `--accent` background; Ts in `--warning`.
- After conversion, all bases uniform `--fg-muted`.
- Methylation-call circles: as per Figure 2 convention.

---

## Figure 5 — Long-read direct methylation calling

**File**: `diagrams/lecture-10/05-longread-methylation.svg`
**Lecture anchor**: §2.4 Long-read direct methylation calling
**ViewBox**: `0 0 1080 440`

### Purpose

Show the two long-read platforms calling methylation from raw signal. Reader leaves able to name the signal modality each uses.

### Content

Two side-by-side panels, plus a shared output pipeline below.

**Left — Oxford Nanopore.** A protein pore (schematic) with DNA threading through. Current-vs-time trace below the pore showing steady state broken by one deeper deflection at a 5mC base. Label: "5mC produces distinctive current signature".

**Right — PacBio HiFi.** A DNA polymerase icon with a nascent strand. Interpulse duration (IPD) vs time trace showing a longer pause at a 5mC position. Label: "5mC → polymerase slowdown → longer IPD".

**Bottom pipeline.** Both signals feed into a neural-network basecaller → outputs per-position (base, modification probability). Annotations: "ONT Dorado / Guppy", "PacBio SMRTLink". Output: a bit vector of methylation calls.

### Style notes

- Pore: `--bg-muted` with `--border-strong` 1.5px outline.
- Current trace: 1px `--fg`, deflection in `--accent`.
- Polymerase: `--accent` filled blob.
- IPD trace: 1px `--fg`, pause highlighted in `--accent`.
- Basecaller box: rounded rectangle with `--accent-bg` fill.

---

## Figure 6 — Beta distributions for methylation proportions

**File**: `diagrams/lecture-10/06-beta-distribution.svg`
**Lecture anchor**: §3.2 Beta distribution for proportion data
**ViewBox**: `0 0 960 420`

### Purpose

Introduce the Beta density visually. Reader leaves able to read the curves and associate shape parameters with biological regimes.

### Content

A single plot area with x-axis `p` from 0 to 1 and y-axis density.

Four Beta curves overlaid:

1. `Beta(1, 1)` — flat uniform (dashed `--fg-muted`).
2. `Beta(2, 8)` — peaked near 0.2 (cobalt; labelled "likely unmethylated").
3. `Beta(8, 2)` — peaked near 0.8 (red; labelled "likely methylated").
4. `Beta(20, 20)` — peaked near 0.5 (amber; labelled "intermediate / mixed").

Each curve labelled with its $(\alpha, \beta)$ and a biological interpretation.

Axis labels and scale; legend with coloured swatches.

Small callout box to the right: "posterior update: Beta(α+m, β+n−m) given m methylated of n reads".

### Style notes

- Curves: filled area under each curve in 0.3 alpha for easy reading.
- Line strokes: 2px each.
- Grid: `--border` 0.3px horizontal lines at density 0.5, 1.0, etc.

---

## Figure 7 — 3C → Hi-C family tree

**File**: `diagrams/lecture-10/07-3c-family-tree.svg`
**Lecture anchor**: §4.2 The 3C family
**ViewBox**: `0 0 1080 360`

### Purpose

Show the genealogy of proximity-capture assays. Reader leaves able to place each assay on the dimensionality scale.

### Content

Five boxes connected by arrows, left-to-right timeline:

1. **3C** (2002). Label: "one vs one". Small cartoon: two loci; one pair read out.
2. **4C** (2006). Label: "one vs all". Cartoon: one viewpoint; all contacts.
3. **5C** (2006). Label: "many vs many (preselected)". Cartoon: grid of selected loci.
4. **Hi-C** (2009). Label: "all vs all". Cartoon: full contact matrix.
5. **Micro-C** (2015). Label: "all vs all, nucleosome-resolution". Cartoon: denser matrix with finer grid.

Above each box: resolution achievable. Below each: paper citation (year + first author).

A dimension axis at the bottom showing "pairs measured per experiment" growing left-to-right by orders of magnitude.

### Style notes

- Boxes: rounded rectangles, `--bg-muted` fill, `--border-strong` outline.
- Arrows between: `--accent` 2px with standard marker.
- Dimension axis: `--fg-muted` line with log-scale tick marks.
- Icons: simple 12px cartoons in `--fg`.

---

## Figure 8 — Hi-C library protocol

**File**: `diagrams/lecture-10/08-hic-protocol.svg`
**Lecture anchor**: §4.3 The Hi-C library protocol
**ViewBox**: `0 0 1080 440`

### Purpose

Six-step schematic making the chemistry concrete. Reader leaves able to name each step.

### Content

Six panels in a 2×3 grid (or single row):

1. **Crosslink.** Nucleus with chromatin; formaldehyde icon; crosslinks forming between nearby fragments.
2. **Digest.** Restriction enzyme icon (scissors) cutting DNA; produce sticky ends.
3. **Fill-in + biotin-label.** Klenow polymerase filling in ends with biotinylated dCTP; biotin shown as orange triangles on the ends.
4. **Ligate proximity.** T4 ligase joining adjacent fragment ends (preferentially ones in 3D proximity); a chimeric fragment drawn spanning two originally-distant loci with a biotin junction in the middle.
5. **Shear + pull down.** Sonication icon; streptavidin beads pulling down biotin-tagged fragments.
6. **Paired-end sequence.** Illumina instrument; two reads mapping to two distinct genomic positions; one entry `C[i, j] += 1` added to the output matrix.

Each panel labelled with step number + name + key reagent.

### Style notes

- Panel boxes: `--bg-muted` fill.
- Reagent icons: simple shapes in `--accent`.
- Biotin: `--warning` triangle.
- DNA: `--fg` curves.
- Arrows between panels: `--accent` 1.5 with marker.

---

## Figure 9 — Raw vs normalised contact matrix

**File**: `diagrams/lecture-10/09-contact-matrix.svg`
**Lecture anchor**: §4.4 The contact matrix as raw output
**ViewBox**: `0 0 1080 460`

### Purpose

Side-by-side comparison of before / after ICE normalisation. Reader leaves recognising the diagonal, distance decay, and the effect of normalisation.

### Content

Two heatmap panels, shared colourbar:

**Left — Raw.** A 50×50 or 80×80 grid contact matrix. Features:
- Bright diagonal.
- Smooth distance-decay away from the diagonal.
- Visible banding artefacts: 2–3 rows/columns darker or brighter than their neighbours (bin-specific biases).
- Log colour scaling annotated.

**Right — ICE-normalised.** Same genomic window. Features:
- Diagonal still visible but less dominant (because we often detrend distance after ICE).
- Banding gone.
- Triangular TAD blocks visible along the diagonal.
- Faint off-diagonal checkerboard (compartment structure).

Annotations on the left: "banding = bin bias" (arrow). On the right: "TADs" (arrow pointing at triangular blocks), "compartments" (arrow pointing at off-diagonal checker).

Shared colour scale at the bottom: pale → dark via `--accent` gradient; log-scaled counts.

### Style notes

- Heatmap cells: `--accent` gradient from `--bg-muted` (0) to deep `--accent` (max).
- Diagonal visibility: built into the data, not a special highlight.
- Bordered cells in the colourbar for reference.

---

## Figure 10 — TADs and A/B compartments

**File**: `diagrams/lecture-10/10-tads-compartments.svg`
**Lecture anchor**: §5.2 TADs
**ViewBox**: `0 0 1080 480`

### Purpose

Multi-scale view of a contact matrix showing three characteristic structures at once. Reader leaves with a mental model of "what's at which scale".

### Content

Main panel (800×800 approx): a ~10 Mb region contact matrix with three visible feature types:

- **Large-scale checkerboard**: plaid A/B compartment pattern visible across the full matrix.
- **Medium-scale triangles**: ~15 TAD blocks along the diagonal, each ~300 kb–1 Mb.
- **Small-scale focal spots**: 4–6 bright focal pixels (loops) near the corners of some TADs.

Three scale markers pointing to each feature with labels: "A/B compartments (Mb)", "TADs (~200 kb – 2 Mb)", "loops (kb–hundreds of kb)".

Below the matrix: a 1D track of the insulation score along the diagonal, with local minima marked as TAD boundaries using small downward-pointing triangles.

Colour-scale annotation and dimension bar.

### Style notes

- Compartment checkerboard: subtle tint variation (lighter cells = B, denser = A).
- TAD blocks: slightly brighter triangular regions along the diagonal.
- Loop dots: distinct `--accent-bright` or `--warning` dots.
- Insulation track: `--accent` line, 1.5px.
- Boundary markers: `--warning` triangles.

---

## Figure 11 — A/B compartments as first eigenvector

**File**: `diagrams/lecture-10/11-compartments-eigenvector.svg`
**Lecture anchor**: §6.1 Hi-C contact matrix as covariance
**ViewBox**: `0 0 1080 460`

### Purpose

The EE-story figure. Show: O/E matrix → correlation matrix → first eigenvector. Reader leaves recognising that compartmentalisation is PCA.

### Content

Three panels side-by-side:

**(a) O/E matrix.** A 60×60 heatmap with checkerboard pattern visible. Label "observed / expected contact matrix · distance-detrended".

**(b) Correlation matrix.** Same size, computed as correlation across columns of (a). Checkerboard much sharper; clear bipartite block structure. Label "correlation matrix (rowwise)".

**(c) First eigenvector.** A 1D track of the first eigenvector plotted against bin index. Positive regions shaded in `--accent`, negative in `--warning`. Below, a binarised track showing A (accent) vs B (warning) compartment assignment.

Arrows between panels showing the computational flow.

### Style notes

- Heatmaps in panels (a) and (b): diverging colour scale centred at zero — negative in cool blue, positive in warm red.
- Eigenvector track (c): filled area above and below zero with diverging colours.
- Arrows with operators labelled: `corr()` between (a) and (b); `eigen()` between (b) and (c).

---

## Figure 12 — TAD calling as change-point detection

**File**: `diagrams/lecture-10/12-tad-insulation.svg`
**Lecture anchor**: §6.2 TAD calling as block-diagonal detection
**ViewBox**: `0 0 1080 480`

### Purpose

Show the reduction from a 2D matrix to a 1D change-point problem. Reader leaves recognising the classical statistical analog.

### Content

Three stacked panels, shared x-axis (bin index):

**(a) Contact matrix** — rectangular window showing a stretch of the matrix with 4 triangular TAD blocks clearly visible. X-axis is bin index.

**(b) Insulation score track** — 1D signal, the insulation score at each bin. Local minima correspond to TAD boundaries. Several local minima clearly visible; boundary detection threshold marked as a horizontal dashed line; minima below threshold marked with arrows.

**(c) Called TADs** — coloured intervals showing the resulting TAD calls, each a different accent-family hue.

Annotation banner at the bottom: "2D → 1D reduction → change-point detection = TAD boundary detection".

### Style notes

- Panel (a): standard heatmap.
- Panel (b): line plot, filled under, with threshold line in `--warning` dashed.
- Panel (c): coloured intervals in 4 accent-family hues.
- Arrows pointing from 2D → 1D → called TADs.
