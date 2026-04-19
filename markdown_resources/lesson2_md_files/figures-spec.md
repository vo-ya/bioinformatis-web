# Lecture 2 — Figures Specification

> **Scope**: Static diagrams for Lecture 2 (Read Alignment).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG unless explicitly marked otherwise. No photographs in Lecture 2 — the content is algorithmic.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering, per `diagram-style-guide.md` §11.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Every figure is **line-art first**: stroke-only by default, fill only where the fill earns the ink.
- Base colors (`--base-a`, `--base-t`, `--base-g`, `--base-c`) are used *only* when actual nucleotide characters are shown. In conceptual diagrams where the alphabet is generic (like the BWT rotation matrix or the FM index diagrams), use the neutral line palette. The genetic meaning takes priority over decorative variety.
- Read and reference strings use **JetBrains Mono**; all other labels use **Inter**.
- Arrows follow the shared `<marker id="arrow-accent">` pattern from `diagram-style-guide.md` §7.

## Figure Budget

Ten figures for a ~3h 50min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Alignment problem overview | Part 1 | Custom SVG |
| 2 | Brute-force sliding window | Part 2 | Custom SVG |
| 3 | Suffix array structure | Part 2 | Custom SVG |
| 4 | Hash k-mer index | Part 2 | Custom SVG |
| 5 | BWT rotation matrix | Part 3 | Custom SVG |
| 6 | BWT LF-mapping invariant | Part 3 | Custom SVG |
| 7 | FM-index checkpoints tradeoff | Part 3 | Custom SVG |
| 8 | Smith-Waterman matrix with traceback | Part 4 | Custom SVG |
| 9 | CIGAR as run-length encoding | Part 4 | Custom SVG |
| 10 | Seed-and-extend pipeline | Part 5 | Custom SVG |

---

## Figure 1 — Alignment Problem Overview

**File**: `diagrams/lecture-02/01-alignment-problem-overview.svg`
**Lecture anchor**: §1.2 Formal statement and scale
**ViewBox**: `0 0 860 360`

### Purpose

Set the vocabulary for the entire lecture. Establish what a "reference", a "read", and an "alignment" look like on the page, so subsequent diagrams can refer to these without re-teaching them.

### Content

A long horizontal band labeled **REFERENCE** (with an uppercase kicker label `REFERENCE · ~3 × 10⁹ bp`). Above it, about eight short colored bars of varying length, each labeled **READ**, scattered at their aligned positions. Three of them carry small annotations:

- one labeled "exact match" pointing straight down to the reference
- one labeled "1 mismatch" with a small red-muted tick showing where the base differs
- one labeled "2 bp insertion" with a small gap notch

Below the reference, a short kicker `SCALE · ~10⁹ reads per run`, in small uppercase Inter.

### Style notes

- Reference is a single horizontal bar, `--fg` stroke, `--bg-muted` fill, 2.5-unit stroke.
- Reads are narrow horizontal bars in `--accent`, 1.5-unit stroke, no fill.
- The three annotations use a thin dashed line (`stroke-dasharray="4 3"`, `--fg-muted`) and a short label in Inter size 11.

### Annotations

- Small monospace coordinate labels at each end of the reference: `0` and `~3 Gb`.
- Bracket on a representative read showing its length with a label `~150 bp`.

---

## Figure 2 — Brute-Force Sliding Window

**File**: `diagrams/lecture-02/02-brute-force-sliding.svg`
**Lecture anchor**: §2.1 Brute force and its complexity
**ViewBox**: `0 0 860 300`

### Purpose

Visualize the naive algorithm's fundamental move — the read slides one position at a time, and at each position every character is compared. Makes the O(|R|·|r|) cost geometric rather than symbolic.

### Content

Three horizontally offset rows, each showing the read aligned to a different position of the reference. Between rows, a small right-pointing arrow labeled `shift +1`. Above the read in each row, a short sequence of tiny vertical ticks (one per character compared). A mismatch is drawn as a small red-muted cross; a match as a small `--fg` dot.

Below the three rows, a bracket labels the total comparison count for those three positions, with an ellipsis and `× (|R| − |r| + 1) positions total`.

### Style notes

- The read monospace text uses JetBrains Mono, weight 500, size 13.
- Compare-ticks are 1-unit strokes, 8 units tall.
- Match dots: filled `--fg`, radius 2. Mismatch crosses: `--error` (muted red), 1.5-unit stroke.
- No color on the reference or the read bases themselves — this is a conceptual figure, not a per-base coloring.

---

## Figure 3 — Suffix Array Structure

**File**: `diagrams/lecture-02/03-suffix-array-structure.svg`
**Lecture anchor**: §2.2 Suffix arrays and binary search
**ViewBox**: `0 0 720 480`

### Purpose

Show the suffix array as a single object: three linked columns — index, SA value (position in R), and the suffix itself. The student should see the sorted suffixes as what justifies binary search.

### Content

Two side-by-side panels separated by a thin vertical rule.

**Left panel** — the input string `BANANA$` rendered with character index labels above it (`0 1 2 3 4 5 6`).

**Right panel** — a table of seven rows:

```
 i   SA[i]   suffix
 0     6     $
 1     5     A$
 2     3     ANA$
 3     1     ANANA$
 4     0     BANANA$
 5     4     NA$
 6     2     NANA$
```

On rows 2 and 3 (the `ANA$` and `ANANA$` rows), a soft accent-bg band highlights both rows, with a bracket on the right side labeled "query `ANA` → range [2, 3]".

### Style notes

- Table is drawn with horizontal rules only, no vertical rules.
- Column headers in uppercase Inter kicker style (`I`, `SA[I]`, `SUFFIX`).
- The suffixes use JetBrains Mono.
- The highlighted band uses `--accent-bg` fill and an accent-color bracket (2-unit stroke).

---

## Figure 4 — Hash k-mer Index

**File**: `diagrams/lecture-02/04-hash-kmer-index.svg`
**Lecture anchor**: §2.3 Hash maps and k-mer indices
**ViewBox**: `0 0 860 380`

### Purpose

Make the k-mer → positions mapping concrete so the space cost is visible. The student should see why 4ᵏ entries plus growing position lists are the problem.

### Content

Two columns.

**Left column**: a small reference string, e.g., `ACGTACGTAC`, with a sliding k-mer window (k = 3) shown at position 0, with an arrow to the right panel.

**Right column**: a hash-table schematic — a vertical stack of buckets, each labeled with a 3-mer (`AAA`, `AAC`, ..., `TTT`, abbreviated with an ellipsis in the middle). From a few buckets — `ACG`, `CGT`, `GTA`, `TAC` — short horizontal chains extend to the right, each chain terminated by a small position label (e.g., `ACG → [0, 4]`).

At the top of the hash-table column, a kicker label `4ᵏ buckets` and at the bottom a small annotation `|R| total positions across all chains`.

### Style notes

- Buckets are uniform rectangles, `--border` stroke, no fill.
- Position chains use 1.5-unit strokes with a small filled circle at each end.
- Ellipsis in the bucket list is three `--fg-subtle` dots, centered in the bucket.

---

## Figure 5 — BWT Rotation Matrix

**File**: `diagrams/lecture-02/05-bwt-rotation-matrix.svg`
**Lecture anchor**: §3.1 The Burrows-Wheeler Transform
**ViewBox**: `0 0 720 520`

### Purpose

Present the canonical BWT construction figure in house style. This is the figure most students will scroll back to during homework. It must be the definitive version.

### Content

Use the example `ACAACG$` (length 7).

**Top section** — the seven unsorted cyclic rotations, one per line, JetBrains Mono.

**Arrow down** — labeled "lex. sort" in accent-color Inter kicker.

**Bottom section** — the seven sorted rotations:

```
$ACAACG
AACG$AC
ACAACG$
ACG$ACA
CAACG$A
CG$ACAA
G$ACAAC
```

Above the sorted block, a small accent-color `F` label over the first column; above the last column, a small accent-color `L` label. Below the sorted block, a boxed rendering of `BWT = L = "GC$AAAC"`.

### Style notes

- The `$` character is rendered in `--fg-subtle` to visually distinguish it from the base alphabet.
- F and L column labels use the uppercase kicker style (Inter 600, letter-spacing 0.08em, `--accent`).
- The BWT output box uses `--bg-inset` fill with a 1-unit `--border` stroke.

---

## Figure 6 — BWT LF-Mapping Invariant

**File**: `diagrams/lecture-02/06-bwt-lf-mapping.svg`
**Lecture anchor**: §3.2 Inverting the BWT: the LF mapping
**ViewBox**: `0 0 860 440`

### Purpose

Visualize the LF property — the i-th occurrence of character `c` in F corresponds to the i-th occurrence of `c` in L. This is the invariant that makes backward search work; it deserves its own figure, separate from the construction figure.

### Content

Continuing from Figure 5, show only the F and L columns side by side, drawn wider than in the previous figure. Draw curved connectors (in `--accent`) between:

- the A's in L (rows 2, 4, 6) and the A's in F (rows 1, 2, 3)
- the C's in L (row 1) and the C's in F (row 4, say)
- a G pair, an `$` pair, etc.

Each connector has a small numeric label (1st, 2nd, 3rd) at both endpoints, showing the i-th-occurrence correspondence.

Below the figure, a short boxed annotation:

> The i-th occurrence of c in F is the i-th occurrence of c in L.

### Style notes

- Connectors are thin (1.5-unit), `--accent` stroke, with soft curves (Bézier). The curvature is functional here — it keeps connectors from overlapping.
- Per-character color is still *not* used; this is not a nucleotide semantic figure. Use only the neutral palette.
- The annotation box uses `--bg-muted` fill.

---

## Figure 7 — FM-Index Checkpoints Tradeoff

**File**: `diagrams/lecture-02/07-fm-checkpoints-tradeoff.svg`
**Lecture anchor**: §3.4 Checkpoints: the space–time tradeoff made explicit
**ViewBox**: `0 0 860 360`

### Purpose

Make the classic space–time tradeoff visible. Two schematic rank-array layouts on the left; a small memory-vs-time plot on the right.

### Content

**Left half** — two stacked schematics:

1. **Dense rank arrays**: a long horizontal strip with `|R|` small tick marks, each one holding a filled dot. Label: `rank stored at every position — O(|R|·|Σ|) memory, O(1) query`.
2. **Sparse rank arrays (checkpointed)**: the same strip, but only every 32nd (or 128th) position is filled. Between checkpoints, a bracket labeled "scan on the fly". Label: `rank stored every d positions — memory ÷ d, query O(d)`.

**Right half** — a simple axis plot: x-axis labeled "memory (bytes/base)", y-axis labeled "query cost (ns)". Two points plotted: one dense (bottom-right), one sparse (top-left), with a dashed curve connecting them indicating the continuous knob.

### Style notes

- The strip schematic is drawn with a 2.5-unit horizontal `--fg` stroke and 1-unit tick marks.
- Filled checkpoint dots use `--accent`, radius 3.
- The tradeoff curve uses `stroke-dasharray="4 3"` in `--fg-muted`.
- Axes are labeled with Inter numerical labels, size 10, `--fg-subtle`.

---

## Figure 8 — Smith-Waterman Matrix with Traceback

**File**: `diagrams/lecture-02/08-smith-waterman-matrix.svg`
**Lecture anchor**: §4.3 Smith-Waterman
**ViewBox**: `0 0 860 540`

### Purpose

The definitive static version of the Smith-Waterman matrix for the two example strings used throughout §4. The artifact (Artifact #5) lets the student play; this figure is the reference they can scroll back to.

### Content

Example strings: reference `GATTACA`, query `GCATGCA` (or similar short pair with a non-trivial traceback). Scoring: match +2, mismatch -1, gap -2.

An 8 × 8 matrix (including the 0-row and 0-column) with every cell filled. The global maximum cell is outlined in `--accent`, 2.5-unit stroke. From it, the traceback path is drawn as a sequence of arrows along cell edges back to a zero — diagonal, up, and left arrows distinguishable by shape.

Below the matrix:

- The aligned strings on two lines (JetBrains Mono, gaps shown as `-`).
- The resulting CIGAR string.

### Style notes

- Cells are a square grid, 1-unit `--border` stroke.
- Cell numbers in JetBrains Mono, size 12, `--fg`.
- Traceback arrows use the `arrow-accent` marker. They are drawn *alongside* the cell boundary, not through cell centers, to avoid occluding the numbers.
- Maximum cell is **outlined**, not filled — the outline plus the converging traceback is enough signal without a color wash.
- The alignment block at the bottom is within a `--bg-inset` rectangle.

---

## Figure 9 — CIGAR as Run-Length Encoding

**File**: `diagrams/lecture-02/09-cigar-rle-diagram.svg`
**Lecture anchor**: §4.4 CIGAR strings
**ViewBox**: `0 0 860 280`

### Purpose

Drive home the EE framing: CIGAR is RLE of an edit-operation sequence. Show the raw traceback symbol sequence on top, the RLE beneath, and the final CIGAR string underneath.

### Content

Three horizontal rows, aligned vertically:

**Row 1** — raw symbol sequence, e.g., `= = = = X = = D D = = = I = = =`, each symbol in a small equal-width box.

**Row 2** — the same sequence with horizontal brackets grouping runs of identical symbols (`4 =`, `1 X`, `2 =`, `2 D`, `3 =`, `1 I`, `3 =`), each labeled above the bracket with the run length.

**Row 3** — the final CIGAR string: `4=1X2=2D3=1I3=`, in JetBrains Mono, size 16, centered.

To the right of the whole figure, a small inset block titled **`EE FRAMING`** (kicker label) reading:

> CIGAR : traceback = run-length encoding : raster line

### Style notes

- Row 1 boxes are uniform width, `--border` stroke, no fill.
- Row 2 brackets are 1.5-unit `--accent` strokes with small tick returns at each end.
- Row 3 is drawn inside a slightly larger `--bg-inset` rectangle with 2-unit `--border-strong` outline.
- The EE framing inset uses `--bg-muted` fill and an accent kicker label.

---

## Figure 10 — Seed-and-Extend Pipeline

**File**: `diagrams/lecture-02/10-seed-and-extend-pipeline.svg`
**Lecture anchor**: §5.1 Two-step alignment as a design pattern
**ViewBox**: `0 0 860 420`

### Purpose

Single diagram that compresses the whole lecture into one picture: a read enters, FM index finds seeds fast, candidate windows get Smith-Watermanned, one alignment wins. Used as the summary figure for Part 5 and as the visual anchor at the top of the wrap-up.

### Content

Horizontal flow, left-to-right, four stages:

1. **Read** — a single horizontal read bar, 50 bp, labeled. JetBrains Mono sequence if there's room, else just the label.
2. **FM index seed lookup** — a small schematic box with the BWT label and the `[sp, ep]` interval, output arrow to a set of "candidate positions" on a reference track.
3. **Reference track** with 3–4 seed-hit positions marked, each with a ±30 bp extension window shown as a shaded `--accent-bg` region.
4. **Smith-Waterman extension** — inside each shaded window, a tiny 3×3 matrix glyph; only the best window has its matrix drawn larger, with a mini traceback.
5. **Output** — a single alignment with its CIGAR string.

Arrows between stages use `arrow-accent` markers; stage boxes have accent-color kicker labels above them: `SEED`, `CANDIDATES`, `EXTEND`, `ALIGNMENT`.

Beneath the diagram, a three-row comparison strip in small text:

```
BWA-MEM    : FM-index variable-length seeds (SMEMs) + banded SW
Bowtie2    : FM-index fixed-length seeds + SSE-accelerated SW
minimap2   : minimizer hash seeds + chaining + banded SW
```

### Style notes

- Stage boxes are uniform height, `--border` stroke, `--bg-muted` fill.
- The "winning" extension window uses `--accent-bg` fill plus a 2.5-unit `--accent` outline — both, to avoid color-alone signaling.
- The comparison strip at the bottom is JetBrains Mono, size 12, `--fg-muted`, in a `--bg-inset` rectangle spanning the figure width.

---

## Cross-Figure Consistency Notes

- **Reference colors**: the reference is always drawn with the same fill (`--bg-muted`) and stroke (`--fg`) across figures — Figures 1, 2, and 10 all use it. A student scrolling should recognize the reference as the reference.
- **FM/BWT neutral palette**: Figures 5, 6, 7 are all "string-as-abstract-alphabet" and stay in the neutral line palette. Figures that show actual nucleotides (1, 2, 10) use the `--base-*` palette.
- **Accent color is scarce**: accent color is reserved for the element the figure is teaching (the highlighted SA range in Figure 3, the LF connectors in Figure 6, the winning extension window in Figure 10). If every element were accented, none would be.
- **Monospace for sequences, always**: any DNA string, BWT string, CIGAR string, or suffix is JetBrains Mono. Never Inter.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All ten figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient, drop shadow, glow, or 3D effect.
- [ ] All sequences are in JetBrains Mono; all labels are in Inter or Source Serif 4.
- [ ] Base colors appear only where actual nucleotides are shown.
- [ ] Accent color appears only on the teaching-point element of each figure.
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
- [ ] No figure duplicates an artifact's dynamic behavior — figures are the static reference; artifacts are the playground.
