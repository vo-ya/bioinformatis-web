# Lecture 12 — Figures Specification

> **Scope**: Static diagrams for Lecture 12 (Population Genetics Fundamentals).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-12.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG; content is population-model-diagram-heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Plots use standard x/y-axis conventions; use Source Serif italic for math (p, q, t, N_e, r²).
- Monospace (JetBrains Mono) for: allele labels (A, a, B, b), numerical values, equations.
- Inter for UI labels and annotations.
- Ancestry / compartment / population colour palette: use the same 6 accent-family hues across Figures 1, 10, 11 (African cobalt, European amber, East Asian red, South Asian green, Admixed violet, Oceanian teal).
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Allele frequencies across a population | Part 1 | Custom SVG |
| 2 | Hardy-Weinberg genotype curves | Part 1 | Custom SVG |
| 3 | Wright-Fisher model schematic | Part 2 | Custom SVG |
| 4 | Genetic drift trajectories | Part 2 | Custom SVG |
| 5 | Selection coefficients and allele trajectories | Part 2 | Custom SVG |
| 6 | LD matrix along a chromosome | Part 3 | Custom SVG |
| 7 | LD decay with distance | Part 3 | Custom SVG |
| 8 | A coalescent tree | Part 4 | Custom SVG |
| 9 | Coalescent ↔ Wright-Fisher duality | Part 4 | Custom SVG |
| 10 | PSMC output — N_e(t) | Part 5 | Custom SVG |
| 11 | ADMIXTURE bar plot | Part 5 | Custom SVG |
| 12 | Selection scan — Manhattan plot | Part 6 | Custom SVG |

---

## Figure 1 — Allele frequencies across a population

**File**: `diagrams/lecture-12/01-allele-frequency.svg`
**Lecture anchor**: §1.1 From variants to frequencies
**ViewBox**: `0 0 1080 400`

### Purpose

The foundational diagram for counting alleles. Reader leaves able to convert a cohort of diploid individuals into allele and genotype frequencies.

### Content

**Left — cohort view.** 50 diploid individuals rendered as a 10×5 grid of "individual" icons. Each individual is two coloured boxes stacked (blue = A allele, red = a allele). Distribution: 22 AA (both boxes blue), 20 Aa (one of each), 8 aa (both red).

**Middle — counting.** A vertical tally panel:
- "Total chromosomes: 100 (50 × 2)"
- "A-alleles: 65"
- "a-alleles: 35"
- "p = #A / 2N = 0.65"
- "q = #a / 2N = 0.35"

**Right — genotype-frequency bar chart.** Three bars: AA count = 22 (44%), Aa = 20 (40%), aa = 8 (16%).

### Style notes

- Individual icons: two small squares (10×10 px each) stacked.
- A allele: `--accent` cobalt.
- a allele: `--error` red.
- Heterozygote: clearly two-colour.
- Bar chart: stacked rectangles, base-palette colours, value labels on top.

---

## Figure 2 — Hardy-Weinberg genotype curves

**File**: `diagrams/lecture-12/02-hardy-weinberg.svg`
**Lecture anchor**: §1.2 Hardy-Weinberg equilibrium
**ViewBox**: `0 0 960 440`

### Purpose

The canonical HWE curves. Reader leaves able to read off genotype frequency from allele frequency and recognise the pattern.

### Content

Single plot area. X-axis: allele frequency $p$ from 0 to 1 (minor-axis label shows $q = 1-p$ decreasing from 1 to 0). Y-axis: genotype frequency from 0 to 1.

Three curves:
- $p^2$ (AA genotype): parabola from (0, 0) to (1, 1). Cobalt. Filled under in pale cobalt.
- $q^2$ (aa): parabola from (0, 1) to (1, 0). Red.
- $2pq$ (Aa): symmetric hump, maximum 0.5 at $p = 0.5$. Green.

Annotations:
- Peak of $2pq$ at $(0.5, 0.5)$ marked and labeled.
- At $p = 0$ and $p = 1$: fixation annotations.
- Side caption: "at any $p$, the three frequencies sum to 1".

### Style notes

- Curve fills at 0.3 opacity; stroke 2 px solid.
- Grid lines every 0.25 on both axes.
- Axis labels in Source Serif italic; values in JetBrains Mono.

---

## Figure 3 — Wright-Fisher model schematic

**File**: `diagrams/lecture-12/03-wright-fisher.svg`
**Lecture anchor**: §2.1 The Wright-Fisher model
**ViewBox**: `0 0 1080 440`

### Purpose

Show the generative sampling rule visually. Reader leaves with a concrete picture of binomial sampling each generation.

### Content

Three generations stacked vertically. Each generation is a horizontal row of $2N = 10$ chromosome slots, rendered as coloured cells (blue A / red a).

**Generation 0**: 5 A, 5 a (p = 0.5). Arranged 10 cells across.

**Generation 1**: 4 A, 6 a (p = 0.4). Each cell has a thin arrow pointing to its "parent" cell in generation 0 (random sampling with replacement — some parents have 0 children, others 2).

**Generation 2**: 3 A, 7 a (p = 0.3). Same arrow structure to generation 1.

Side annotations:
- "each child chromosome samples a parent with probability 1/(2N) per chromosome"
- "p_t+1 | p_t ~ Binomial(2N, p_t) / (2N)"
- Stochastic drift visible: frequencies wander even with no selection.

### Style notes

- Cells: small rectangles.
- Arrows: `--fg-muted` thin (0.6 px) with tiny markers.
- Some "parent" cells have multiple children (thick arrow style); some have none (no outgoing arrow).
- Grid to emphasise generations.

---

## Figure 4 — Genetic drift trajectories

**File**: `diagrams/lecture-12/04-drift-trajectories.svg`
**Lecture anchor**: §2.2 Genetic drift
**ViewBox**: `0 0 1080 460`

### Purpose

Show drift as a function of population size. Reader leaves with the intuition that small populations drift fast and large ones drift slowly.

### Content

Two side-by-side plots, same axes: X = generations (0 to 500), Y = allele frequency (0 to 1).

**Left panel — Small population (N = 50).** 8 coloured trajectories starting at $p_0 = 0.5$. Several hit 0 (fixation of a) or 1 (fixation of A) within 500 generations. Visibly volatile. Annotate "N = 50, ~half fix by gen 200".

**Right panel — Large population (N = 5000).** 8 trajectories starting at $p_0 = 0.5$. All stay clustered near 0.5 across 500 generations. Annotate "N = 5000, all frequencies within (0.42, 0.58)".

### Style notes

- Trajectories: each a different accent-family hue, thin (1 px) to allow visual overlap.
- Horizontal dashed lines at $p = 0$ and $p = 1$ marking absorbing boundaries.
- Panel outlines: `--border-strong`.

---

## Figure 5 — Selection coefficients and allele trajectories

**File**: `diagrams/lecture-12/05-selection.svg`
**Lecture anchor**: §2.3 Selection
**ViewBox**: `0 0 1080 440`

### Purpose

Contrast drift alone, weak selection, and strong selection. Reader leaves with the intuition of the critical $2Ns$ scaling.

### Content

Single plot area. X = generations (0 to 500). Y = allele frequency of the favoured allele (0 to 1).

Three sets of trajectories:

- **$s = 0$ (neutral)** in grey: 5 runs, zigzag drift around 0.5.
- **$s = 0.01$ (weak)** in blue: 5 runs, each a noisy upward trend, some fix by gen 500.
- **$s = 0.1$ (strong)** in red: 5 runs, all fix cleanly within ~100 generations, close to the deterministic curve.

Overlay the deterministic theoretical curve (no noise) for $s = 0.1$ as a black dashed line.

Side annotations:
- "s = 0: drift only"
- "s = 0.01: slow selection + drift"
- "s = 0.1: selection dominates"
- Box stating "critical scaling is $2Ns$, not $s$ or $N$ alone".

### Style notes

- Colour-coded trajectories per regime; label-legend bottom-right.
- Deterministic curve: `--fg` 1.5 px dashed.

---

## Figure 6 — LD matrix along a chromosome

**File**: `diagrams/lecture-12/06-ld-matrix.svg`
**Lecture anchor**: §3.3 LD decay and haplotype blocks
**ViewBox**: `0 0 1080 480`

### Purpose

The canonical LD-heatmap view. Reader leaves recognising block-diagonal structure and recombination hotspots.

### Content

Large square heatmap: 100 × 100 SNPs. Cell $(i, j)$ coloured by $r^2$ between SNP $i$ and SNP $j$.

Visual features:
- Strong diagonal (each SNP in perfect LD with itself).
- Triangular blocks along the diagonal, 10–20 SNPs wide, showing high-LD regions.
- Sharp white (low LD) boundaries between blocks — recombination hotspots.
- Off-diagonal LD mostly low (pale).

Below the matrix, a 1D genomic-position track showing:
- SNP positions as tick marks.
- Recombination hotspot positions marked with vertical dashed lines, matching the LD-block boundaries.

Annotations on the right:
- "Block 1: 25 SNPs in high LD"
- "Recombination hotspot at position 30 breaks block 1 from block 2"

### Style notes

- Heatmap cells: gradient from pale (`#f6f4ee`) to deep red (`#c4342c`) via intermediate amber.
- Diagonal explicitly boxed in `--border-strong`.
- Hotspot markers: `--warning` 1 px dashed vertical.

---

## Figure 7 — LD decay with distance

**File**: `diagrams/lecture-12/07-ld-decay.svg`
**Lecture anchor**: §3.3 LD decay and haplotype blocks
**ViewBox**: `0 0 960 420`

### Purpose

Show the exponential LD-decay curve. Reader leaves able to read off typical LD scales for humans.

### Content

Single plot. X-axis: genomic distance in kb (log scale, 0.1 to 1000). Y-axis: average $r^2$ (0 to 1).

**Decay curve**: starts at $r^2 = 1$ at distance 0 bp, decays smoothly to ~0.05 at 1 Mb. Fit shape: $r^2 \approx \exp(-\lambda d)$ with $\lambda$ chosen so half-max sits at ~30 kb.

Annotations:
- Point marker at 30 kb / $r^2 = 0.5$, labelled "half-max".
- Point marker at 100 kb / $r^2 = 0.2$.
- Point marker at 1 Mb / $r^2 = 0.05$.
- Shaded band around the curve for "typical range across populations (10th-90th percentile)".

Side annotation: "human 1000G cohorts typically show $\lambda \approx 1/30$ per kb".

### Style notes

- Main curve: `--accent` 2 px.
- Shaded band: `--accent-bg` at 0.5 opacity.
- Annotation dots: small filled circles.

---

## Figure 8 — A coalescent tree

**File**: `diagrams/lecture-12/08-coalescent-tree.svg`
**Lecture anchor**: §4.2 The standard coalescent
**ViewBox**: `0 0 1080 460`

### Purpose

The canonical coalescent-tree diagram. Reader leaves able to interpret waiting times and the asymmetry between early and late intervals.

### Content

A binary-tree diagram with:
- Time flowing bottom (present) to top (past).
- 8 sampled chromosomes at the bottom, labelled (S1, S2, …, S8).
- Pairs of lineages coalescing at different times. The specific realisation:
  - T_8 = short (S1+S2 coalesce near the bottom).
  - T_7 = short.
  - T_6 = medium.
  - T_5, T_4 = longer.
  - T_3 = very long.
  - T_2 = longest (the final two lineages take the most time).
- All waiting-time intervals labelled to the left with the expected values: $T_k \sim \text{Exp}(\binom{k}{2}/(2N))$, $\mathbb{E}[T_k] = 2N/\binom{k}{2}$.

Right panel: a sliver showing the ratio of waiting times — the T_2 bar is larger than the sum of T_3..T_8.

Top annotation: "MRCA" at the root node, with TMRCA ≈ $4N(1 - 1/n)$ ≈ $4N \cdot 7/8 = 3.5N$ labelled.

### Style notes

- Tree branches: `--fg` 2 px.
- Sample leaves: coloured circles (one per sample).
- Coalescent nodes: filled circles.
- Waiting-time labels in JetBrains Mono.

---

## Figure 9 — Coalescent ↔ Wright-Fisher duality

**File**: `diagrams/lecture-12/09-coalescent-wf-duality.svg`
**Lecture anchor**: §4.3 From coalescent trees to summary statistics
**ViewBox**: `0 0 1080 460`

### Purpose

Show that forward WF and backward coalescent are two views of the same underlying process. Reader leaves understanding the duality.

### Content

Two panels side-by-side.

**Left panel — Forward Wright-Fisher.** A 100-generation WF realisation, with $2N = 20$ chromosomes per generation. Each generation is a horizontal row of 20 small coloured circles. Lines connect each chromosome to its parent (random from previous generation). The visual is a "spaghetti" of parent-child lines from top (past) to bottom (present).

Highlight: 5 specific present-day chromosomes (bottom row) and trace their lineages backward through the spaghetti; those traced lineages are in bold colour. All other lineages are faded.

**Right panel — Backward coalescent of the same 5 samples.** Just the bold lineages from the left panel, but rearranged into a clean coalescent tree. Shows coalescent events at the same generation times as on the left.

Arrow between panels labelled "same process, two views — the coalescent keeps only the lineages reaching the present".

### Style notes

- Faded lineages in left panel: `--fg-subtle` 0.4 opacity.
- Traced lineages: distinct accent-family hues with 1.5 px stroke.
- Coalescent tree on right panel: same colours, tree-renderer layout.

---

## Figure 10 — PSMC output, N_e(t) for human populations

**File**: `diagrams/lecture-12/10-psmc-ne.svg`
**Lecture anchor**: §5.1 PSMC — from a single genome
**ViewBox**: `0 0 1080 440`

### Purpose

The canonical PSMC result plot. Reader leaves able to interpret ancient bottlenecks and population expansions.

### Content

Log-log plot. X-axis: years before present (10³ to 10⁷). Y-axis: effective population size $N_e$ (10³ to 10⁶).

Three population trajectories (step functions, as PSMC produces):
- **African** (cobalt): stays at $N_e \approx 20{,}000$ for most of history; slight dip around 100 kya; slight increase in recent past.
- **European** (amber): tracks African until ~70 kya; dramatic bottleneck to $N_e \approx 2{,}000$ around 50 kya (out-of-Africa bottleneck); recovery in last 20 kya.
- **East Asian** (red): similar to European — tracks African, then bottleneck around 50 kya, slightly stronger, recovery more recent.

Annotations:
- Shaded region labelled "out-of-Africa bottleneck ~50 kya".
- Shaded region labelled "recent expansion (agriculture?) last 10 kya".
- Note near the left edge: "PSMC loses resolution in the most recent past (hasn't accumulated enough coalescence signal yet)".

### Style notes

- Step functions rendered as connected horizontal-vertical segments.
- Log axes with gridlines at each decade.
- Shaded event regions in `--accent-bg` with text labels.

---

## Figure 11 — ADMIXTURE bar plot

**File**: `diagrams/lecture-12/11-admixture.svg`
**Lecture anchor**: §5.3 Admixture inference
**ViewBox**: `0 0 1080 400`

### Purpose

The canonical admixture visualisation. Reader leaves recognising the bar-chart format and understanding individual-level ancestry mixtures.

### Content

Stacked bar chart. X-axis: 100 individuals (sorted by continental origin). Y-axis: ancestry proportion (0 to 1, summing to 1 per individual).

4 ancestry components, each a colour:
- Component 1 (African, cobalt)
- Component 2 (European, amber)
- Component 3 (East Asian, red)
- Component 4 (South Asian, green)

Visual structure:
- Left 25 individuals: nearly 100% component 1 (African).
- Next 25: nearly 100% component 2 (European).
- Next 25: nearly 100% component 3 (East Asian).
- Right 25: nearly 100% component 4 (South Asian).
- A cluster in the middle-right: admixed individuals with mixtures (e.g. half European, half East Asian).

Below: a row of text labels showing group origins. Above: heading "K=4 ancestry components fit to 1000 Genomes subset".

### Style notes

- Bars: each individual is a thin vertical stack (~10 px wide).
- Colour scheme matches other L12 figures.
- Group divider lines between populations.

---

## Figure 12 — Selection scan, Manhattan plot

**File**: `diagrams/lecture-12/12-selection-scan.svg`
**Lecture anchor**: §6.1 Tajima's D and neutrality tests
**ViewBox**: `0 0 1080 440`

### Purpose

The canonical Manhattan plot. Reader leaves recognising the format and understanding that peaks are selection candidates.

### Content

X-axis: genomic position across 22 chromosomes (laid out as a single concatenated axis, with chromosome boundaries marked).

Y-axis: $-\log_{10}(p\text{-value})$ for a selection statistic (iHS), 0 to 30.

Dots: ~50,000 SNPs plotted. Most at low $-\log_{10}(p)$. Chromosomes alternate in colour (cobalt / grey).

Tall peaks (above y=15) clearly visible at specific loci. Annotations point to 4 known selection-target loci:
- **LCT** on chr 2 (lactase persistence; Northern European selection).
- **SLC24A5** on chr 15 (skin pigmentation; European selection).
- **ABCC11** on chr 16 (earwax type; East Asian).
- **EDAR** on chr 2 (hair thickness; East Asian).

A horizontal dashed line at $-\log_{10}(p) = 7.3$ labelled "genome-wide significance".

### Style notes

- Dots: small (1.5 px), coloured per chromosome.
- Labelled peaks: larger dots with arrow-annotations.
- Chromosome boundary lines: very thin grey vertical.
- Threshold line: `--warning` 1 px dashed.
