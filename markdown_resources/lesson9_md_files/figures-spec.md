# Lecture 9 — Figures Specification

> **Scope**: Static diagrams for Lecture 9 (ChIP-seq, ATAC-seq, and Peak Calling).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-09.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. Content is algorithm- and chemistry-heavy; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales cleanly up to 1200 px.
- Monospace (JetBrains Mono) for: nucleotide letters in sequences, gene symbols, numerical values (p-values, scores, fragment sizes), PWM weights. Inter for UI labels.
- DNA-strand colours use the locked base palette (`--base-a`, `--base-t`, `--base-g`, `--base-c`) **only** when literal bases are drawn; any other DNA cartoons use `--fg` / `--fg-muted`.
- Arrows follow the shared `<marker id="arrow-accent">` and `<marker id="arrow-muted">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture, slightly heavier than prior lectures to cover the detection + motif + sequence-model trio. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | What ChIP-seq and ATAC-seq measure | Part 1 | Custom SVG |
| 2 | ChIP-seq workflow | Part 1 | Custom SVG |
| 3 | ATAC-seq Tn5 chemistry | Part 2 | Custom SVG |
| 4 | Fragment-length distribution — nucleosome laddering | Part 2 | Custom SVG |
| 5 | Read pile-up at a TF binding site | Part 3 | Custom SVG |
| 6 | MACS2 local-Poisson detection | Part 3 | Custom SVG |
| 7 | Narrow vs broad peaks | Part 3 | Custom SVG |
| 8 | Differential accessibility — MA plot | Part 4 | Custom SVG |
| 9 | TF motif as PWM and sequence logo | Part 5 | Custom SVG |
| 10 | Motif scanning as matched filtering | Part 5 | Custom SVG |
| 11 | ATAC footprint at a TF motif | Part 5 | Custom SVG |
| 12 | Enformer architecture sketch | Part 6 | Custom SVG |

---

## Figure 1 — What ChIP-seq and ATAC-seq measure

**File**: `diagrams/lecture-09/01-chip-vs-atac-overview.svg`
**Lecture anchor**: §1.2 ChIP-seq — where proteins bind DNA
**ViewBox**: `0 0 1080 500`

### Purpose

Side-by-side comparison of the two assays. Reader should leave able to say what each assay measures, what its output is, and why the two are complementary.

### Content

Two panels at the top with a shared genomic cartoon below.

**Left panel — ChIP-seq.** Cartoon icons showing: crosslinked cells → sonicated chromatin fragments → antibody pulling down a specific TF/histone-marked fragment → the resulting fragment DNA being sequenced. Arrow leads to a read pile-up visualisation showing sharp narrow coverage at bound sites.

**Right panel — ATAC-seq.** Cartoon: Tn5 transposase icons inserting adapters into open regions of chromatin → short fragments output → sequenced. Arrow leads to a read pile-up visualisation showing coverage at open regulatory regions.

**Bottom band — shared cartoon**: a ~5 kb stretch of genome with drawn elements:
- Nucleosomes as beige-grey oval bumps at regular intervals.
- Two TFBSs where specific TFs (coloured squares) are bound.
- An open/accessible stretch between nucleosomes.

Aligned below, two coverage tracks (ChIP in accent, ATAC in amber) showing where each assay produces signal: ChIP spikes at the TF-bound positions; ATAC spikes at the open region and falls within the nucleosome-protected regions.

### Style notes

- Panel boxes: `--border` 1px, rounded 6px.
- TF icons: coloured squares (accent, amber).
- Nucleosomes: `--bg-muted` ovals.
- Coverage tracks: filled area charts under a baseline.
- Arrows: `--accent` 1.5 with standard marker.

---

## Figure 2 — ChIP-seq workflow

**File**: `diagrams/lecture-09/02-chip-seq-workflow.svg`
**Lecture anchor**: §1.2 ChIP-seq
**ViewBox**: `0 0 1080 360`

### Purpose

Four-step workflow making the chemistry concrete. Reader leaves able to name each step.

### Content

Horizontal sequence of 4 schematic panels, connected by `--accent` arrows:

1. **Crosslink.** A cell outline with chromatin inside and a "HCHO" label with a small molecule structure; little covalent bonds drawn between a TF (coloured square) and the DNA it contacts.
2. **Sonicate.** Cross-linked chromatin being broken into fragments by wavy sound-wave icons; labelled "~300 bp fragments".
3. **Immunoprecipitate.** A Y-shaped antibody capturing one TF-bound fragment from a mix; others shown discarded.
4. **Reverse crosslinks + sequence.** Labelled "pure DNA → Illumina" with a small sequencer icon; output = set of reads.

Bottom annotation line: "~30M reads → align → peak calling".

### Style notes

- Each panel in a rounded rectangle, `--bg-muted` fill.
- DNA strands in `--fg` curves.
- TF coloured squares in `--accent`.
- Antibody Y-shape in `--warning`.
- Arrows `--accent` 1.5 with standard marker.

---

## Figure 3 — ATAC-seq Tn5 chemistry

**File**: `diagrams/lecture-09/03-atac-tn5.svg`
**Lecture anchor**: §1.3 ATAC-seq and DNase-seq
**ViewBox**: `0 0 1080 420`

### Purpose

Show how Tn5 transposase accesses only open chromatin and produces fragments of specific sizes based on nucleosome spacing. Reader leaves able to explain the nucleosome-ladder outcome.

### Content

**Top** — a ~1.5 kb chromatin cartoon: 6 nucleosomes (beige-grey ovals each with DNA wrapped as curved lines on top) spaced along a horizontal DNA line. Between nucleosomes: short "linker" stretches of open DNA. One segment has a larger open region (nucleosome-free promoter or enhancer).

**Middle** — Tn5 transposase enzymes drawn as small teal blobs approaching the DNA. Arrows show insertion sites falling **only** in the open linker regions (not through nucleosomes). Mark inserted sequencing adapters as small coloured bars.

**Bottom** — a list of fragment classes produced:
- **Sub-nucleosomal** (< 100 bp) — both Tn5 cuts in the same open stretch, zero nucleosomes spanned.
- **Mono-nucleosomal** (~150–180 bp) — cuts flank a single nucleosome.
- **Di-nucleosomal** (~300 bp) — cuts span two adjacent nucleosomes.
- **Tri-nucleosomal** (~500 bp) — three nucleosomes, much rarer.

Each fragment class rendered as a small DNA segment with its length annotation.

### Style notes

- Nucleosomes: `--bg-muted` ovals with `--border-strong` outline.
- DNA: `--fg` 2px curves.
- Tn5 enzymes: `--base-g` (teal-green) blobs.
- Adapters: short `--accent` bars.
- Fragment-class labels: Inter 10 with JetBrains Mono for sizes.

---

## Figure 4 — Fragment-length distribution, nucleosome laddering

**File**: `diagrams/lecture-09/04-fragment-length.svg`
**Lecture anchor**: §2.2 Fragment-length distributions
**ViewBox**: `0 0 960 440`

### Purpose

The canonical ATAC-seq QC plot. Reader leaves able to identify the ladder peaks and read a real fragment-size distribution.

### Content

Large x-axis: fragment length (bp), range 0 to ~800, log-scale.
Y-axis: read density (linear).

**Primary curve**: ATAC-seq fragment distribution in `--accent` — clear multi-modal shape with peaks at:
- ~50 bp (sub-nucleosomal, tallest)
- ~180 bp (mono-nucleosomal)
- ~330 bp (di-nucleosomal)
- ~500 bp (tri-nucleosomal, faint)

Annotate each peak with a label pointing to the peak and its biological meaning.

**Overlay curve**: a lighter-weight dashed distribution in `--warning` representing a ChIP-seq sonicated library — a single smooth mode centred at ~300 bp with no sub-structure. Label it "ChIP — sonicated, no laddering" to contrast.

Legend at top-right. Annotation at bottom: "laddering is the ATAC QC signature — absence = bad experiment".

### Style notes

- Axis labels in JetBrains Mono 10 `--fg-muted`.
- ATAC curve: `--accent` 2px solid with fills `--accent-bg` 0.5 opacity.
- ChIP curve: `--warning` 1.5px dashed.
- Peak labels: small arrows with `--fg-muted` 0.5px lines + Inter 10 text.

---

## Figure 5 — Read pile-up at a TF binding site

**File**: `diagrams/lecture-09/05-read-pileup.svg`
**Lecture anchor**: §3.1 The detection problem
**ViewBox**: `0 0 1080 360`

### Purpose

Genome-browser-style view of what peak calling needs to detect. Reader leaves recognising the shape of the detection task.

### Content

**Top band** — a 5 kb genomic coordinate axis (in kb) with tick marks.

**Middle band** — the coverage track: a per-base height plot over the 5 kb window. Baseline coverage of ~3 reads/base across most of the window; one central ~300 bp region shows a spike rising to ~25 reads/base; another smaller rise at ~2 kb reaching ~8 reads/base (ambiguous signal); a flat tail elsewhere.

**Callouts** on the main peak:
- "peak summit" label with arrow pointing to the highest point.
- "peak boundaries" brackets marking the ~200 bp region above threshold.
- "local background" shaded region showing the baseline around the peak.

**Bottom band** — annotation track showing that the main peak falls on a promoter region (a rectangle labelled "TSS of GENE_X"), demonstrating the connection between signal and biology.

### Style notes

- Coverage track filled area in `--accent` over a `--bg-muted` plot background.
- Annotations: `--warning` for the peak summit pointer; `--fg-muted` brackets for peak boundaries; `--accent-bg` shaded region for local background.
- Axis in JetBrains Mono 9.

---

## Figure 6 — MACS2 local-Poisson detection

**File**: `diagrams/lecture-09/06-macs2-local-poisson.svg`
**Lecture anchor**: §3.2 MACS2 algorithm
**ViewBox**: `0 0 1080 460`

### Purpose

The algorithmic core of peak calling. Reader leaves able to state how the local-adaptive $\lambda$ is computed and why it works.

### Content

Three columns:

**Column 1 — Coverage track with candidate window.** A section of coverage track with a small box marking "candidate window (200 bp, count c = 24)".

**Column 2 — Local windows.** Same coverage track but now with nested larger windows drawn: 1 kb, 5 kb, 10 kb surrounding the candidate. Each window labelled with its estimated per-200bp $\lambda$:
- $\lambda_{\text{bg}}$ = 0.5 (genome-wide)
- $\lambda_{1000}$ = 2.1
- $\lambda_{5000}$ = 3.8
- $\lambda_{10000}$ = 2.9
- $\lambda_{\text{input}}$ = 1.2

Bottom text: "$\lambda_{\text{local}}$ = max(…) = 3.8".

**Column 3 — Poisson null and test.** A Poisson($\lambda = 3.8$) PMF plotted as a histogram; mark the candidate count $c = 24$ on the far right tail. Text: "P(X ≥ 24 | λ = 3.8) = 2.4 × 10⁻¹⁴ → reject H₀".

Between columns, arrows with `→` marking the algorithmic flow.

### Style notes

- Poisson PMF histogram: `--accent` bars.
- Rejection region: `--warning` shading under the tail.
- Local-$\lambda$ labels in JetBrains Mono; equations in Source Serif italic where needed.
- Box annotations on the nested windows: `--border-strong` 0.5px dashed.

---

## Figure 7 — Narrow vs broad peaks

**File**: `diagrams/lecture-09/07-narrow-vs-broad.svg`
**Lecture anchor**: §3.4 Narrow vs broad peak modes
**ViewBox**: `0 0 1080 460`

### Purpose

Show the three characteristic signal shapes side-by-side. Reader leaves able to pick the right MACS2 mode per mark type.

### Content

Three stacked coverage tracks over the same 50 kb coordinate axis:

**(a) CTCF ChIP (TF) — narrow peaks.** ~5 sharp peaks spaced irregularly, each ~200 bp wide, tall (~30× baseline). Label each peak with a small vertical line and "TFBS".

**(b) H3K4me3 (active promoter mark) — broader peaks.** ~3 peaks, each ~1 kb wide, more gradual shoulders. Label each with "TSS".

**(c) H3K27me3 (silenced domain) — broad domain peaks.** Two wide plateaus each spanning 10–20 kb, not resolved into individual peaks. Label as "silenced domain".

Under each track: a small annotation tag indicating the appropriate caller mode — "--narrow (default)" for (a) and (b), "--broad" for (c).

### Style notes

- Track heights scaled so all three fit comfortably.
- Different accent hues per track: (a) cobalt, (b) teal, (c) amber.
- Mode-annotation tags in a small pill badge style.

---

## Figure 8 — Differential accessibility, MA plot

**File**: `diagrams/lecture-09/08-differential-accessibility.svg`
**Lecture anchor**: §4.1 The DB / DA problem
**ViewBox**: `0 0 900 440`

### Purpose

Show DA results in the same visual language as L6's differential-expression MA plot. Reader leaves recognising the callback.

### Content

MA scatter plot:
- X-axis: mean log₂ count across samples (0 to 16).
- Y-axis: log₂ fold change between conditions (−5 to +5), with zero marked.

**Cloud of points** distributed around y=0 at moderate x, showing natural scatter.

**Significant peaks**:
- Cluster of red dots at positive y (gained accessibility) in the x range ~8–12.
- Cluster of blue dots at negative y (lost accessibility) similar x range.

Threshold lines at y = ±1 (2-fold change) as faint horizontal dashed lines.

**Highlighted peaks** with callouts:
- One red-cluster peak labelled "ENH1 near GENE_A (up 3.2×)".
- One blue-cluster peak labelled "PROM2 near GENE_B (down 2.1×)".

Inset legend: red = up in treated; blue = down in treated; grey = not significant.

### Style notes

- Scatter dots: grey for non-significant, red (`--error`) for up, blue (`--accent`) for down.
- Threshold lines: `--fg-muted` 0.5px dashed.
- Callout annotations: small connecting lines + Inter 9 labels.

---

## Figure 9 — TF motif as PWM and sequence logo

**File**: `diagrams/lecture-09/09-pwm-motif.svg`
**Lecture anchor**: §5.1 TF binding motifs and the PWM
**ViewBox**: `0 0 1080 420`

### Purpose

The PWM — the central representation of sequence specificity. Reader leaves able to read both the matrix form and the logo form, and translate between them.

### Content

**Left — PWM numerical matrix.** A 4 × 19 grid representing a CTCF motif (approximate). Rows labelled A, C, G, T (each base in its base-palette colour). Columns numbered 1 to 19. Cell contents: log-odds weights (numbers roughly between −2 and +3), coloured by value — positive in `--accent`, negative in `--fg-subtle`, zero-ish in `--bg-muted`.

**Right — Sequence logo.** Standard IC-height logo: at each position, a stack of A/C/G/T letters with heights proportional to their bits of information; letter height encodes frequency, letter identity shown by the character itself in base-palette colours. Tall core "CCCTC" / "GCCCCCTC" pattern in the middle (positions ~8–14); variable flanks.

**Between them** — an arrow and the text: "PWM rendered as a logo — same information, two views".

**Bottom annotation line**: "each column = one base position · column entropy = conservation · IC = max entropy − observed entropy".

### Style notes

- PWM cells: small rectangles (about 20 × 20 px each), coloured per weight.
- Numbers inside cells in JetBrains Mono 9.
- Logo letters: 4 base-palette colours, sized by information content.
- Matrix labels in Inter 10.

---

## Figure 10 — Motif scanning as matched filtering

**File**: `diagrams/lecture-09/10-matched-filter.svg`
**Lecture anchor**: §5.3 Motif scanning as a matched filter
**ViewBox**: `0 0 1080 460`

### Purpose

The EE framing — PWM scan = correlation of template against one-hot sequence. Reader leaves recognising the identity.

### Content

Three horizontally aligned bands:

**Top — the 4-channel template (PWM).** 4 rows × L columns (L ≈ 12), each cell a coloured square by weight. Label "PWM (4 channels × L columns)".

**Middle — the 4-channel signal (one-hot DNA).** 4 rows × N columns (N ≈ 60), each column has exactly one cell "hot" (coloured, based on which base is present); other cells pale. Label "sequence as 4-channel one-hot" and annotate a few positions with the corresponding nucleotide letter.

**Bottom — the scan output.** A 1D score track showing the sliding-window inner-product output. Mostly flat-low; a clear spike at one position where the sequence matches the PWM. Annotate: "hit: score = 18.3, p < 10⁻⁴".

**Between template and signal**: arrows indicating sliding window — dashed box around the template shown at three candidate positions along the signal, with small thin arrows from template-positions to output-values.

**Bottom caption band**: "inner product at each position = Σ_i w_{s_i, i} = PWM score = matched-filter output on 4-channel signal".

### Style notes

- Template cells: PWM colour coding from Figure 9.
- Signal cells: pale grey cells; hot cells in base-palette colour.
- Output track: filled area in `--accent` under a baseline.
- Dashed-box sliding window: `--warning` 1.5px dashed.

---

## Figure 11 — ATAC footprint at a TF motif

**File**: `diagrams/lecture-09/11-atac-footprint.svg`
**Lecture anchor**: §5.4 ATAC footprinting
**ViewBox**: `0 0 1080 460`

### Purpose

Show the footprint as an ensemble phenomenon — noisy individually, clean in aggregate. Reader leaves with the √N averaging intuition.

### Content

**Main (central) panel — aggregated footprint.** A curve showing Tn5 cut density (y) as a function of distance from motif centre (x, −200 to +200 bp). Shape: elevated flanking shoulders (~2× baseline) and a clear dip at the motif centre (~0.3× baseline). Width of dip: ~20 bp. Labelled "aggregated over 5,000 CTCF motif instances".

**Side panel (left) — individual noisy traces.** 5 small stacked sub-plots, each showing the noisy cut-count trace at a single motif instance. All look essentially random; no clear dip visible in any one of them.

**Arrow between the two panels**: labelled "sum over N sites → SNR × √N".

**Bottom annotation band**: "motif-centred aggregation recovers the footprint as coverage-trough signal — the TF's 'binding shadow'".

### Style notes

- Aggregated curve: `--accent` 2.5 with `--accent-bg` fill under.
- Individual noisy traces: `--fg-muted` 1px lines, pale.
- Arrow between panels: `--accent` 2px with standard marker.

---

## Figure 12 — Enformer architecture sketch

**File**: `diagrams/lecture-09/12-enformer-architecture.svg`
**Lecture anchor**: §6.1 Enformer and its siblings
**ViewBox**: `0 0 1080 400`

### Purpose

Show the CNN + transformer hybrid architecture at a conceptual level. Reader leaves able to point at the components and name them.

### Content

Left-to-right flow:

1. **Input**. A long DNA strand labelled "100 kb DNA sequence (200K bases × 4 channels)". Show a short highlighted stretch with base letters visible.

2. **CNN block** — 3 stacked grey rectangles labelled "conv + ReLU", with downsampling indicated by arrow showing resolution going from 1 bp → 128 bp bins. Output label: "local sequence features".

3. **Transformer block** — 3 stacked coloured boxes labelled "transformer layer (self-attention)", connected by arrows. Output label: "long-range contextualised features".

4. **Output heads** — branches into multiple parallel prediction heads, each labelled:
   - "CAGE (5000 tracks, TSS activity)"
   - "DNase/ATAC (~700 tracks)"
   - "TF ChIP (~4000 tracks)"
   - "Histone marks (~2000 tracks)"

Output representation: a per-128-bp-bin signal prediction for each track.

**Caption band at bottom**: "Enformer · Avsec et al. 2021 · predicts epigenomic signal from sequence · forward pointer to Lecture 16".

### Style notes

- Input DNA: base-palette letters against `--bg-muted`.
- CNN rectangles: `--bg-muted` with `--border-strong`.
- Transformer boxes: `--accent-bg` with `--accent` outline.
- Output heads: coloured pills in distinct accent hues.
- Arrows: `--accent` 1.5 standard.
