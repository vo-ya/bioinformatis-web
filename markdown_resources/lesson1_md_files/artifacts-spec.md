# Lecture 1 — Interactive Artifacts Specification

> **Purpose**: Detailed build specs for the six interactive HTML artifacts embedded in Lecture 1.
> Each artifact is a single self-contained `.html` file. See `website-spec.md` for the shared design system and the iframe embedding contract.

---

## Global Conventions

### Tech stack

- **Pure HTML + CSS + vanilla JavaScript.** No build step, no framework.
- **Charting** (where needed): Chart.js loaded from a CDN. Plotly or D3 only if genuinely required.
- **Math rendering** (where needed): KaTeX via CDN (`katex.min.js` + `katex.min.css` + `auto-render.min.js`).
- **No external network calls** beyond CDN library fetches. No analytics.

### Shared files

Every artifact imports:

```html
<link rel="stylesheet" href="_shared/artifact-theme.css">
<script defer src="_shared/resize.js"></script>
```

- `artifact-theme.css` — CSS custom properties matching the site design system (colors, typography, spacing). See `website-spec.md` §4.2–4.3.
- `resize.js` — posts a `{type: 'artifact-resize', height}` message to the parent window whenever the artifact's height changes (`ResizeObserver` on `document.body`). The parent lecture page listens and sets `iframe.style.height` accordingly. This removes the need to hardcode iframe heights.

### Design rules (all artifacts)

1. Use **CSS variables from the design system only**. No hex colors inline.
2. DNA/RNA bases always use `--base-a / t / g / c / u / n`. Establish the color association in artifact #1 and never deviate.
3. Monospace (`--font-mono`) for any sequence, FASTQ line, or quality string. Sans (`--font-sans`) for UI chrome.
4. Every artifact has a visible `<h1>` or `<h2>` title at the top — it's part of the artifact, not the embedding page, so it works standalone.
5. Every artifact includes a short "What this shows" caption (1–2 sentences) beneath the title.
6. Controls grouped in a clearly-delimited panel; outputs in a clearly-delimited panel. Don't mix.
7. Reasonable defaults — the artifact should be interesting the moment it loads, with no user input required.
8. Target minimum width: 600px. Below that, gracefully degrade (some artifacts will need horizontal scroll; that's OK).
9. **No keyboard traps.** All controls reachable by Tab. Focus visible.
10. **No browser storage APIs** (no `localStorage`, no `sessionStorage`) — state lives in memory.

### Naming

- Files: `NN-name-kebab.html` (`01-dna-explorer.html`, …).
- Each artifact has a semantic `<title>` matching its display title.
- Root container: `<main class="artifact" data-artifact="dna-explorer">`.

---

## Artifact #1 — DNA Explorer

**File**: `01-dna-explorer.html`
**Lecture placement**: §1.3
**Teaching goal**: Make abstract "A pairs with T" tactile. Introduce the shared base-color palette used in every subsequent artifact.
**Estimated class time**: 5 min

### UI layout

```
┌─────────────────────────────────────────────────────────┐
│  DNA Explorer                                            │
│  Type or paste a DNA sequence to see its properties.    │
├─────────────────────────────────────────────────────────┤
│  Input:                                                  │
│  [ GATTACAGATTACAGATTACA__________________ ] [Random]   │
│                                                          │
│  Length: 21    GC content: 33.3%    Valid: ✓            │
├─────────────────────────────────────────────────────────┤
│  5'   G A T T A C A G A T T A C A G A T T A C A   3'  │
│       | | | | | | | | | | | | | | | | | | | | |        │
│  3'   C T A A T G T C T A A T G T C T A A T G T   5'  │
│                                                          │
│  Complement:         CTAATGTCTAATGTCTAATGT              │
│  Reverse:            ACATTAGACATTAGACATTAG              │
│  Reverse complement: TGTAATCTGTAATCTGTAATC              │
├─────────────────────────────────────────────────────────┤
│  Base composition bar: [AAAA|TTTTTT|GGG|CCCC]           │
│  A: 8 (38%)  T: 7 (33%)  G: 3 (14%)  C: 3 (14%)        │
└─────────────────────────────────────────────────────────┘
```

### Behavior

- Text input accepts A/C/G/T/N (case-insensitive, normalize to uppercase, strip whitespace).
- Invalid characters highlighted in red inline; show an error below the input.
- **Random** button generates a random 30-base sequence.
- As the user types, **live update** all panels.
- Each base rendered in its signature color (CSS var). Pairing bars between strands visible.
- GC content = `(G+C) / (A+C+G+T) × 100`, shown to 1 decimal.
- Base composition bar is a horizontal stacked bar with counts + percentages.

### Implementation notes

- No libraries needed.
- Use a `<div>` with `display: flex; gap: 2px; font-family: var(--font-mono);` for each strand; one child span per base with class `base-A / T / G / C / N`.
- For the complement strand, render bottom-aligned; for reverse-complement, show as a plain sequence block.

### Edge cases

- Empty input: show placeholder state with example text grayed out.
- Sequences longer than the viewport: wrap in a horizontally scrollable container; maintain alignment between top and bottom strands.
- Very long sequences (>500 bp): cap the double-strand rendering at 200 bp and show a message; continue to compute composition stats for the full sequence.

---

## Artifact #2 — Central Dogma Translator

**File**: `02-central-dogma.html`
**Lecture placement**: §2.2
**Teaching goal**: The "one conceptual artifact" of the lecture — DNA → mRNA → protein, codon by codon, with genetic code redundancy made visible.
**Estimated class time**: 7 min
**Complexity**: Medium-high. This is also the **design-system proof-of-concept**; build it first.

### UI layout

```
┌──────────────────────────────────────────────────────────┐
│  Central Dogma Translator                                 │
│  Watch DNA transcribe to mRNA and translate to protein.  │
├──────────────────────────────────────────────────────────┤
│  DNA (template):  [ ATGGCTTTCGTTAAATAG________ ] [Rand]  │
│                                                           │
│  [▶ Play]  [⏸ Pause]  [⏭ Step]  [↺ Reset]                │
│  Speed: [slow ─────●──── fast]   Frame: ● 1 ○ 2 ○ 3      │
├──────────────────────────────────────────────────────────┤
│  DNA:   5' A T G G C T T T C G T T A A A T A G 3'       │
│            │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │ │           │
│  mRNA:  5' A U G G C U U U C G U U A A A U A G 3'       │
│           ╰─────╯ ╰─────╯ ╰─────╯ ╰─────╯ ╰─────╯         │
│            AUG     GCU     UUC     GUU     …            │
├──────────────────────────────────────────────────────────┤
│  Protein:    M       A       F       V       K    STOP  │
│            (Met)  (Ala)  (Phe)  (Val)  (Lys)            │
├──────────────────────────────────────────────────────────┤
│  Codon table (click a codon to see synonymous codons):   │
│  [ 4×16 compact codon-table grid, colored by amino acid ]│
└──────────────────────────────────────────────────────────┘
```

### Behavior

- User types DNA sequence. mRNA is derived by replacing T→U on the **sense** strand (we render the sense strand as input for simplicity; note in caption that biologically, transcription reads the template strand).
- **Play** animates: character by character, highlight DNA base, show it "becoming" mRNA, then group into codons, then produce amino acids. Adjustable speed.
- **Step** advances by one codon.
- **Frame** selector: shift reading frame by 0, 1, or 2 bases. Recompute accordingly.
- **Stop codon** (UAA, UAG, UGA) ends translation with a visible "STOP" marker.
- **Start codon** (AUG) highlighted specially the first time it appears in the chosen frame (for the "translation starts here" teaching moment).
- **Codon table** at the bottom:
  - 4×16 compact grid OR the traditional 4×4×4 layout — pick one, annotate clearly.
  - Each codon cell colored by its amino acid (group amino acids into ~5 chemical categories: hydrophobic, polar, acidic, basic, special — use distinct colors that aren't in the base palette).
  - **Clicking a codon** highlights all synonymous codons (same amino acid) to show degeneracy.
  - Current codon being translated (in the main animation) pulses in the table.

### Implementation notes

- Single codon table as a JS constant. Include 3-letter and 1-letter amino acid abbreviations.
- Amino acid categories for coloring:
  - Hydrophobic: A, V, L, I, M, F, W, P
  - Polar: S, T, N, Q, Y, C
  - Acidic: D, E
  - Basic: K, R, H
  - Special: G, (STOP)
- Animation via `requestAnimationFrame` + a state machine (`{step: 'dna-highlight' | 'transcribe' | 'codon' | 'translate', index: i}`).
- Reserve space so layout doesn't shift as amino acids appear — pre-allocate slots.

### Edge cases

- Sequence length not divisible by 3: show the trailing 1–2 bases with a "partial codon" label, don't translate them.
- No start codon in the frame: still translate from position 0; annotate.
- Stop codon before the end of the sequence: annotate "rest of the sequence is untranslated 3' region."

---

## Artifact #3 — Sequencing Cost Explorer

**File**: `03-cost-explorer.html`
**Lecture placement**: §3.5
**Teaching goal**: Visualize the cost collapse; let students see Moore's Law fall behind; connect immediately to the 23andMe / Dante Labs discussion.
**Estimated class time**: 5 min
**Complexity**: Low-medium. Uses Chart.js.

### UI layout

```
┌────────────────────────────────────────────────────────────┐
│  Sequencing Cost Explorer                                   │
│  NHGRI's genome cost vs Moore's Law, 2001–present.          │
├────────────────────────────────────────────────────────────┤
│  Axis: [●] Log  [ ] Linear    Show: [✓] Cost per genome    │
│                                     [✓] Cost per Mb        │
│                                     [✓] Moore's Law        │
│                                                             │
│  [ Big line chart, log-y by default, 2001 → today ]        │
│  [ Annotations at key events: HGP complete, NGS arrival,   │
│    $1000 genome, etc. ]                                    │
├────────────────────────────────────────────────────────────┤
│  Family sequencing calculator                               │
│  Year: [────●────── 2015]   People: [4]                    │
│                                                             │
│  In 2015, sequencing 4 genomes would have cost ~$16,000.   │
│  Today (2024): ~$400.                                       │
└────────────────────────────────────────────────────────────┘
```

### Behavior

- Line chart with two data series (cost per genome, cost per Mb) plus a dashed Moore's Law reference line (halves every 24 months, anchored to the 2001 genome cost for visual comparison).
- Log-y toggle (default ON).
- Hover: tooltip shows `{year, cost per genome, cost per Mb, dominant technology}`.
- Annotations on the chart for: 2003 (HGP complete), 2008 (NGS mainstream — curve break), 2014 ($1000 genome), 2022 ($200 genome).
- **Family calculator**: year slider + people spinner → multiplies cost-per-genome at that year × count. Shows a natural-language line comparing to today's cost.

### Data

Embed as a JSON constant in the file. Approximate NHGRI values (verify against latest NHGRI cost data when building):

```js
const costData = [
  { year: 2001, costPerGenome: 100_000_000, costPerMb: 5292.39, tech: "Sanger (HGP)" },
  { year: 2003, costPerGenome:  70_000_000, costPerMb: 3898.64, tech: "Sanger" },
  { year: 2006, costPerGenome:  20_000_000, costPerMb:  651.81, tech: "Sanger / early NGS" },
  { year: 2008, costPerGenome:   1_000_000, costPerMb:   15.03, tech: "Illumina GA, 454, SOLiD" },
  { year: 2010, costPerGenome:      50_000, costPerMb:    0.35, tech: "HiSeq 2000" },
  { year: 2012, costPerGenome:       7_666, costPerMb:    0.09, tech: "HiSeq 2500" },
  { year: 2014, costPerGenome:       1_000, costPerMb:    0.04, tech: "HiSeq X Ten" },
  { year: 2019, costPerGenome:         942, costPerMb:    0.01, tech: "NovaSeq 6000" },
  { year: 2022, costPerGenome:         525, costPerMb:   0.006, tech: "NovaSeq, DNBSEQ" },
  { year: 2024, costPerGenome:         200, costPerMb:   0.003, tech: "NovaSeq X, Ultima" },
];
```

> **Note**: Verify figures at build time against https://www.genome.gov/about-genomics/fact-sheets/Sequencing-Human-Genome-cost — the above are approximate.

### Implementation notes

- Chart.js, one `<canvas>`. Use `type: 'line'`, logarithmic y-axis.
- Annotations via the `chartjs-plugin-annotation` plugin (also CDN-loadable).
- Moore's Law line: start at 2001 cost, halve every 2 years: `moore(y) = cost_2001 / 2^((y - 2001)/2)`.

### Edge cases

- Users might drag year slider to future years — cap at the latest data year.

---

## Artifact #4 — Illumina Base Caller

**File**: `04-illumina-basecaller.html`
**Lecture placement**: §4.4
**Teaching goal**: Show that base-calling is a signal-detection problem on image data. Make the Illumina SBS cycle concrete.
**Estimated class time**: 8 min
**Complexity**: High. Along with #5, this is the most involved artifact.

### UI layout

```
┌────────────────────────────────────────────────────────────┐
│  Illumina Base Caller                                       │
│  Watch clusters on a simulated flow-cell tile emit          │
│  fluorescence, cycle by cycle, and get called as bases.     │
├────────────────────────────────────────────────────────────┤
│  Template sequence: [ GATTACAGATTACA____ ] [Rand]          │
│  Clusters: [ ──●── 8 ]   Cycles to run: [ ──●── 14 ]       │
│  Noise σ: [ ──●── 0.15 ]   Phasing: [ ──●── 0.03 ]         │
│  Chemistry: [●] 4-channel  [ ] 2-channel                   │
│                                                             │
│  [▶ Play]  [⏸]  [⏭ Step]  [↺ Reset]                       │
├────────────────────────────────────────────────────────────┤
│  Flow-cell tile (current cycle: 7)                         │
│  ┌──────────────────────────────────────┐                  │
│  │  🔴  🟢  🔵  🟡  🔴  🟢  🔵  🟡       │  ← clusters    │
│  │  .. 4-color fluorescent spots ..      │                  │
│  └──────────────────────────────────────┘                  │
│                                                             │
│  Per-cluster intensity bars (A/C/G/T), current cycle:      │
│  Cluster 1: [▇▃▁▂]  → called: A  Q=37                      │
│  Cluster 2: [▁▇▁▁]  → called: C  Q=40                      │
│  ...                                                        │
├────────────────────────────────────────────────────────────┤
│  Basecall output so far:                                   │
│  Cluster 1:  GATTAC(A)____    true: GATTACA                │
│  Cluster 2:  GATTAC(A)____    true: GATTACA                │
│  Cluster 3:  GATTAC(T!)____   ← error, noise roll          │
│  ...                                                        │
└────────────────────────────────────────────────────────────┘
```

### Behavior

- Simulated clusters, each reading the same template (or user can pick different templates per cluster — start simple, same for all).
- Each cycle:
  1. Determine the "true" base about to be incorporated.
  2. Generate 4 intensity values: true channel gets `1.0 + N(0, σ)`, others get `0.0 + N(0, σ)`. Add a small cross-contamination term and a phasing term (probability `p_phasing` that the cluster read ahead or lagged by one position).
  3. Call the base as `argmax(intensities)`.
  4. Compute a Phred quality from the intensity ratio (simple monotonic mapping; doesn't need to be real).
- **Fluorescence visualization**: render clusters as circles with fill color mixed from the 4 channel intensities. In 4-channel mode: A=red, C=green, G=blue, T=yellow (historical Illumina convention — or pick and document). In 2-channel mode: show how "G = dark" and "A = both."
- **Intensity bars** update per cycle per cluster.
- **Running basecall string** accumulates per cluster, with errors highlighted in red.
- **Controls**: play/pause/step/reset, noise slider, phasing slider, cluster count, chemistry toggle.

### Pedagogical moments

- Crank noise up to σ=0.4: watch error rate climb. Teaches signal-to-noise intuition.
- Crank phasing up: watch reads drift out of register. Teaches why late-cycle quality drops.
- Switch to 2-channel: watch the "dark G" problem surface — uninformative cycles at the end of runs.

### Implementation notes

- Canvas-based rendering for the flow cell (SVG would work too; canvas is fine for ~100 clusters).
- Intensity bar chart: small Chart.js instances OR hand-rolled `<div>` bars with widths.
- State machine: `{cycle: n, clusters: [{template, reads, qualities}]}`. Step function advances all clusters one cycle.
- Gaussian noise: Box-Muller transform. Fine.

### Edge cases

- User enters a template shorter than `cycles to run`: stop that cluster at template length; show "—" for remaining cycles.
- Cluster count > 16: render as a grid rather than a single row.

---

## Artifact #5 — Nanopore Squiggle Decoder

**File**: `05-nanopore-squiggle.html`
**Lecture placement**: §5.3
**Teaching goal**: Show that nanopore sequencing is a time-series-to-sequence problem. Make k-mer context in the pore visceral.
**Estimated class time**: 6 min
**Complexity**: High. Alongside #4, the most involved artifact.

### UI layout

```
┌────────────────────────────────────────────────────────────┐
│  Nanopore Squiggle Decoder                                  │
│  DNA threads through a pore; ionic current drops as bases   │
│  block it. A basecaller turns current into sequence.        │
├────────────────────────────────────────────────────────────┤
│  Template: [ GATTACAGATTACAGATTACA ] [Rand]                │
│  [▶ Play]  [⏸]  [⏭]  [↺]                                  │
│  Speed: [slow ──●── fast]   Noise: [──●── 5 pA]            │
├────────────────────────────────────────────────────────────┤
│  Pore cartoon (5 bases in the constriction):               │
│                                                             │
│      ══════════════════════════════  ← membrane            │
│            │  A G T T A  │  ← bases in pore                │
│            └─────────────┘                                  │
│      ══════════════════════════════                        │
│          ↑ base entering       ↓ base exiting              │
│     GATTACAGATT            │  AGATTACA                     │
│                                                             │
├────────────────────────────────────────────────────────────┤
│  Ionic current trace (scrolling window, pA vs time):       │
│  [ ~~~~~\______/-----\___/-----\_____ ] ← live squiggle    │
│                                                             │
│  Basecaller output (streaming): GATTA...                   │
│  True sequence (for comparison): GATTA...                  │
└────────────────────────────────────────────────────────────┘
```

### Behavior

- User enters a template. Play simulates the DNA translocating through the pore.
- At each time step:
  - A **5-base window** of the template is "in" the pore.
  - The expected current = a lookup from a (simulated) k-mer → current table. Generate this table deterministically from the 5-mer string (e.g., hash to a value in a reasonable pA range ~60–120).
  - Emit current = expected + Gaussian noise.
  - Move the window forward one base after some variable dwell time.
- **Scrolling squiggle plot**: shows last ~5 seconds of current. Levels are visibly stepped.
- **Basecaller** (educational, not a real ML model): a simple rule — classify each level segment by which k-mer its mean is closest to, emit the center base of the matched k-mer. Demonstrate it works with low noise and degrades with high noise.
- **Pore cartoon** is a schematic, not physically accurate — just enough to make the "5 bases at once" point stick.

### Pedagogical moments

- Raise the noise slider: watch the basecaller make mistakes. Critical teaching moment: real basecallers use deep learning precisely because this problem is hard.
- Pause and point at a single current level and say "this one level is determined by 5 bases — if any of them changes, the level shifts." That's the whole intuition.
- Highlight homopolymers (runs like AAAAA): even with zero noise, the current doesn't change while the homopolymer slides through — the basecaller literally can't count how many A's there are. That's why homopolymer indels are the classic ONT error mode.

### Implementation notes

- Canvas for the scrolling current trace (smoother than SVG for constantly-updating lines).
- K-mer-to-current map: `const kmerCurrent = hashToRange(kmer, 60, 120)` — stable pseudo-random but deterministic so the same template always produces the same squiggle.
- Pore cartoon: SVG, static structure with `<text>` elements for the current k-mer, animated position.

### Edge cases

- Template shorter than 5: disable play, show a prompt.
- Extreme noise: basecall may degrade to nonsense — that's the point; annotate it.

---

## Artifact #6 — FASTQ Inspector + Phred Calculator

**File**: `06-fastq-inspector.html`
**Lecture placement**: §6.2–6.3
**Teaching goal**: Demystify the 4-line FASTQ record. Make the Q ↔ P_error ↔ ASCII triad concrete.
**Estimated class time**: 8 min
**Complexity**: Low-medium.

### UI layout — two panels stacked

**Panel A — FASTQ Inspector**

```
┌────────────────────────────────────────────────────────────┐
│  FASTQ Inspector                                            │
│  Paste a FASTQ record. Each quality char is decoded.        │
├────────────────────────────────────────────────────────────┤
│  Encoding: [●] Phred+33 (modern)  [ ] Phred+64 (legacy)    │
│                                                             │
│  [ multi-line textarea, 4 lines prefilled with example ]   │
│                                                             │
│  [Load example]  [Clear]  [Random]                         │
├────────────────────────────────────────────────────────────┤
│  Parsed:                                                    │
│    ID:          HWI-D00360:5:H814YADXX:1:1101:1230:2106   │
│    Description: 1:N:0:ACAGTG                               │
│    Length:      27                                          │
│    Mean Q:      36.2    Min Q: 33   Max Q: 40             │
│                                                             │
│  Per-base view (hover for details):                        │
│   G   A   T   C   C   T   A   C   T   G   ...             │
│   B   B   B   F   F   F   F   F   F   F   ...             │
│  Q=33 Q=33 Q=33 Q=37 Q=37 ...                              │
│                                                             │
│  Per-base error probability (bar chart):                   │
│  [ small bar per position, height = P_error, log y-axis ]  │
└────────────────────────────────────────────────────────────┘
```

**Panel B — Phred Calculator**

```
┌────────────────────────────────────────────────────────────┐
│  Phred ↔ Error ↔ ASCII                                      │
├────────────────────────────────────────────────────────────┤
│  Change any field; others update.                           │
│                                                             │
│  Q (Phred score):     [ 30 ]                                │
│  P (error probability): [ 0.001 ]  (1 in 1,000)             │
│  ASCII char (Phred+33): [ ? ] (value 63)                    │
│                                                             │
│  Equation: Q = −10 · log₁₀(P)                              │
│                                                             │
│  Intuition:  Q30 = 30 dB of confidence, 99.9% accurate.    │
└────────────────────────────────────────────────────────────┘
```

### Behavior

**Inspector**:
- Parse pasted text as a FASTQ record (expect 4 lines; forgive a trailing newline; show validation errors inline).
- Validate: line 1 starts with `@`, line 3 starts with `+`, line 2 length == line 4 length.
- Per-base view: three stacked rows (base, ASCII char, Q). Hover a column → popover with `{position, base, char, ASCII code, Q, P_error, 1-in-N}`.
- Bar chart of per-base P_error (log y-axis).
- Summary stats: length, mean/min/max Q.
- Encoding toggle switches between Phred+33 and Phred+64 — quality scores recompute.

**Calculator**:
- Three bound fields: Q, P, ASCII char. Edit any → others update.
- Show the equation with KaTeX.
- Short intuition caption referencing dB.

### Implementation notes

- Chart.js for the per-base error bar chart.
- KaTeX for the equation.
- All computation is O(N) and trivial — no perf concerns.
- Example FASTQ records to rotate through (hardcode 3–4):
  - Short clean Illumina read (like the one in the lecture).
  - A read with a bad tail (show quality drop-off).
  - A PacBio HiFi read (very high Q).

### Edge cases

- Mismatched seq/qual lengths: highlight and refuse to parse.
- Non-printable ASCII in quality string: show ASCII value anyway, mark as suspicious.
- Phred+64 decoding of a Phred+33 string will produce negative Q values; catch and warn.

---

## Build Order Recommendation

1. **#2 — Central Dogma Translator**: build first. It's visually distinctive and exercises most of the design-system patterns (colors, monospace, controls, animation, codon table). Approving this approves the aesthetic for the rest.
2. **#6 — FASTQ Inspector**: simple and high-utility; good second build.
3. **#1 — DNA Explorer**: quick win, reinforces base palette.
4. **#3 — Cost Explorer**: Chart.js-centric, moderate effort.
5. **#4 — Illumina Base Caller** and **#5 — Nanopore Squiggle**: both complex, both simulation-heavy. Save for last.

## Testing checklist (for each artifact)

- [ ] Opens standalone in browser with no console errors.
- [ ] Embeds cleanly in the lecture page iframe.
- [ ] `resize.js` reports height correctly; iframe grows/shrinks to fit.
- [ ] All controls reachable by keyboard; focus visible.
- [ ] Renders sensibly at 600px, 900px, 1400px widths.
- [ ] Works in Chrome, Firefox, Safari (latest two versions each).
- [ ] No external requests except documented CDN libraries.
- [ ] Design tokens (colors, fonts) all sourced from CSS vars — no inline hex values.
