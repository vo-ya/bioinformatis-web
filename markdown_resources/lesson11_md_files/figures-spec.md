# Lecture 11 — Figures Specification

> **Scope**: Static diagrams for Lecture 11 (Long Reads and the Pangenome).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-11.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Base-palette colours (`--base-a/t/g/c`) only where literal nucleotide letters appear (Figures 3, 8).
- Graph elements use consistent visual vocabulary: nodes as rectangles with sequence labels; directed edges with arrow markers; paths as coloured traces overlaying the graph.
- Monospace (JetBrains Mono) for: DNA letters, GFA field values, node IDs, numerical values. Inter for UI labels. Source Serif for math in Figure 10 Viterbi update equation.
- Escape `&`, `<`, `>` as XML entities in text content (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | PacBio HiFi vs Oxford Nanopore in 2024 | Part 1 | Custom SVG |
| 2 | SVs short reads can't span | Part 2 | Custom SVG |
| 3 | Tandem repeats resolved by long reads | Part 2 | Custom SVG |
| 4 | Haplotype phasing — short vs long reads | Part 2 | Custom SVG |
| 5 | T2T assembly — before vs after | Part 2 | Custom SVG |
| 6 | Reference bias — concrete example | Part 3 | Custom SVG |
| 7 | Linear reference vs graph genome | Part 3 | Custom SVG |
| 8 | GFA format anatomy | Part 3 | Custom SVG |
| 9 | Seed-and-extend on a DAG | Part 4 | Custom SVG |
| 10 | Viterbi on a linear chain vs DAG | Part 5 | Custom SVG |
| 11 | vg toolchain workflow | Part 4 | Custom SVG |
| 12 | HPRC-v1 pangenome structure | Part 6 | Custom SVG |

---

## Figure 1 — PacBio HiFi vs Oxford Nanopore in 2024

**File**: `diagrams/lecture-11/01-hifi-vs-ont.svg`
**Lecture anchor**: §1.2 PacBio HiFi vs Oxford Nanopore in 2024
**ViewBox**: `0 0 1080 440`

### Purpose

Two-axis position of long-read platforms vs short reads. Reader leaves able to state the read-length × accuracy tradeoff and place each tech generation in the right quadrant.

### Content

Main scatter plot:

- X-axis: read length, log scale, 50 bp → 2 Mb. Label major ticks (100 bp, 1 kb, 10 kb, 100 kb, 1 Mb).
- Y-axis: per-base accuracy in Phred Q units, linear, Q10 → Q35. Label Q20, Q30 explicitly (with dashed horizontal lines).

Points plotted:

- **Illumina** (short-read short-accurate): 150 bp × Q30+. Cluster in upper-left.
- **PacBio HiFi 2018**: 15 kb × Q30. Labelled with year.
- **PacBio HiFi 2024** (Revio): 20 kb × Q32. Labelled.
- **ONT 2018** (R9.4): 10 kb × Q12. Labelled with year.
- **ONT 2022** (R10.4): 20 kb × Q18. Labelled.
- **ONT 2024** (R10.4.1 Q20+): 20 kb × Q22, with ultra-long sub-population extending to 1 Mb × Q22. Labelled.

An annotation arrow showing "reference-grade frontier" across the Q20 horizontal line — the line both long-read platforms now cross. Shade the "reference-grade long-read" region (accuracy ≥ Q20 AND length ≥ 10 kb) in accent-bg.

Legend placing each point's platform name with a small icon.

### Style notes

- Point markers: ONT dots in one hue family (e.g. `--base-t` family), PacBio in another (`--accent` family), Illumina grey.
- Frontier line: `--warning` 1.5px dashed.
- Shaded region: `--accent-bg` at 0.3 opacity.
- Axis labels in Inter 10; tick labels in JetBrains Mono 9.

---

## Figure 2 — SVs short reads can't span

**File**: `diagrams/lecture-11/02-sv-long-vs-short.svg`
**Lecture anchor**: §2.1 SVs resolved directly
**ViewBox**: `0 0 1080 440`

### Purpose

Show why short reads fail on mid-size SVs and long reads succeed by direct spanning. Reader leaves with an intuitive SV-resolution argument.

### Content

Three horizontal bands aligned to a shared genomic coordinate axis (0 kb, 5 kb, 10 kb labelled).

**Top — reference structure.** A ~10 kb reference segment with a 5 kb tandem duplication in the middle clearly marked. Annotations: "reference" and "duplicate copy".

**Middle — short-read coverage.** ~15 short reads (150 bp each) rendered as little rectangles above the reference. Key visual: reads stop at duplication boundaries; a few reads "split-map" (shown as two halves with a connecting line) hint at the event. A coverage track underneath shows a subtle coverage-doubling inside the duplicated region. Label: "inference from fragment geometry — fragile".

**Bottom — long-read coverage.** ~5 long reads (20 kb each), rendered as longer rectangles. At least two reads fully span the duplication. The reads directly contain both duplicated copies. Label: "direct spanning — robust".

### Style notes

- Short reads: `--fg-muted` small rectangles.
- Split-read halves: `--warning` with thin connecting line.
- Long reads: `--accent` larger rectangles with subtle gradient fill.
- Coverage track: filled area chart.
- Annotation arrows: `--accent` 1px with standard marker.

---

## Figure 3 — Tandem repeats resolved by long reads

**File**: `diagrams/lecture-11/03-tandem-repeats.svg`
**Lecture anchor**: §2.2 Spanning repeats
**ViewBox**: `0 0 1080 440`

### Purpose

Show the repeat-counting problem concretely — why short reads can't count repeat units when the repeat length exceeds the read length. Reader leaves understanding the HTT / CAG example as motivation.

### Content

Three bands:

**Top — the HTT locus.** A gene model of HTT with the CAG repeat region highlighted at exon 1. Below the gene, a zoom-in showing the repeat sequence as a series of `CAG` triplets in JetBrains Mono, highlighted. Label "CAG repeat × N" with pathogenic threshold (N ≥ 36) annotated.

**Middle — short-read ambiguity.** ~6 short reads (150 bp, covering ~50 CAG units each) rendered above a 180 bp repeat region (60 CAG units, just above the pathogenic threshold). Each short read contains some CAGs but cannot see the boundaries; the read count doesn't disambiguate between e.g. 50 and 70 units. Annotation: "which read came from which section? repeat count uncertain".

**Bottom — long-read counting.** A single 3 kb long read that spans the entire repeat region plus flanks. The CAG units are counted directly by reading off the molecule. Annotation: "repeat count: 60 units → pathogenic, above threshold of 36".

### Style notes

- Gene model: standard gene rendering, exons as rectangles, introns as lines.
- Repeat units: 3-bp coloured blocks in base-palette (C blue, A red, G amber), repeating pattern.
- Short reads: small rectangles, some highlighted to show ambiguous coverage.
- Long read: a single wide rectangle spanning the repeat + flanks.
- Clinical-threshold annotation: `--warning` horizontal dashed line.

---

## Figure 4 — Haplotype phasing, short vs long reads

**File**: `diagrams/lecture-11/04-phasing-short-vs-long.svg`
**Lecture anchor**: §2.3 Haplotype phasing directly from reads
**ViewBox**: `0 0 1080 440`

### Purpose

Show why short reads give ambiguous phase and long reads give deterministic phase. Reader leaves with the "link-via-molecule" intuition.

### Content

Three bands:

**Top — reference with 5 heterozygous SNPs.** A 20 kb reference segment. 5 heterozygous SNP positions marked: each site has ref / alt alleles (A/G, C/T, etc.). The true haplotypes are labelled: maternal carries {A, C, T, G, A}; paternal carries {G, T, C, A, T}.

**Middle — short-read phasing.** ~12 short reads (150 bp each), each covering 1–2 heterozygous sites at most. Annotations show that each individual read's 2-site combinations don't uniquely link the haplotypes. A "?" next to the inferred phase. Label: "statistical phasing from LD panel required".

**Bottom — long-read phasing.** 2 long reads (20 kb each), each spanning all 5 heterozygous sites. Read 1 shows {A, C, T, G, A} — maternal haplotype directly. Read 2 shows {G, T, C, A, T} — paternal. Label: "direct phase from a single molecule — no statistics needed".

A small side annotation: "short-read block: ~3 kb typical; long-read block: chromosome-scale".

### Style notes

- Heterozygous SNP positions: vertical ticks on the reference with coloured base labels.
- Short reads: small rectangles.
- Long reads: wide rectangles with allele labels at each SNP position.
- Ambiguity markers: `--warning` question marks.
- Haplotype colours: maternal `--accent`, paternal `--warning`.

---

## Figure 5 — T2T assembly, before vs after

**File**: `diagrams/lecture-11/05-t2t-before-after.svg`
**Lecture anchor**: §2.4 Telomere-to-telomere (T2T) assemblies
**ViewBox**: `0 0 1080 420`

### Purpose

Show concretely what the T2T-CHM13 assembly completed. Reader leaves with a sense of where the ~200 Mb addition sits (centromeres, acrocentric arms, segmental duplications).

### Content

Two chromosome bars, same chromosome (use chr 9 as the canonical example — has a visible centromere and acrocentric p-arm).

**Top — GRCh38.** Chromosome 9 drawn with:
- Left (p-arm): truncated/gapped before the centromere. Gap region labelled "N-masked centromere".
- Centre: centromere shown as a small hatched region.
- Right (q-arm): complete.
- Annotations: "unplaced contigs" (15 Mb), "centromeric gap" (5 Mb), "heterochromatin gap" (20 Mb).

**Bottom — T2T-CHM13.** Same chromosome, now:
- p-arm extends fully to the telomere with the acrocentric short arm resolved (labelled "rDNA arrays, alpha-satellite").
- Centromere fully sequenced (no gaps, labelled "alpha-satellite HOR, fully resolved").
- q-arm the same as GRCh38 (mostly — some added segmental duplications).
- Annotations showing the new ~15% of chr 9 that was added.

Between the two chromosomes, a small "added" total: "+200 Mb genome-wide".

### Style notes

- Chromosome bars: `--fg-muted` outline, `--bg-muted` fill for euchromatin.
- Centromere: hatched pattern with `--accent` diagonals.
- Gaps: `--warning` dashed boundaries.
- Resolved regions in T2T: `--accent` fill with subtle opacity.
- Telomere caps: small circular end markers.

---

## Figure 6 — Reference bias, concrete example

**File**: `diagrams/lecture-11/06-reference-bias.svg`
**Lecture anchor**: §3.2 Reference bias — a concrete example
**ViewBox**: `0 0 1080 460`

### Purpose

Show reference-bias clinical impact concretely. Reader leaves understanding that GRCh38 linear alignment under-calls non-European variants.

### Content

Three bands.

**Top — locus with three common variants.** A 2 kb locus with three variant sites (V1, V2, V3), each labelled with its major-population carrier:

- V1: common in European populations (present in GRCh38).
- V2: common in African populations (absent from GRCh38).
- V3: common in East Asian populations (absent from GRCh38).

**Middle band — three populations × linear alignment result.** Three sub-panels per population showing reads from that population aligned to GRCh38. Each sub-panel renders per-variant call-quality:

- European reads: V1 called correctly (solid blue), V2/V3 called as ref (no-call, grey).
- African reads: V1 called as ref (miss), V2 missed (flagged red as "no-call / false ref"), V3 missed.
- East Asian reads: V1 called as ref, V2 missed, V3 missed (red flag).

**Bottom — same reads aligned to pangenome graph.** Three sub-panels again; all three variants are called correctly in every population (all green). Label: "pangenome graph includes all three variants as graph branches".

Annotations: "clinical impact: rare-disease variant commonly found in a given population can be missed entirely at the alignment step on GRCh38".

### Style notes

- Linear-alignment panels: `--bg-muted` background.
- Pangenome panels: `--accent-bg` background.
- Correct calls: `--success` green fill.
- Missed calls: `--error` red flag.
- No-call: `--fg-muted` grey.

---

## Figure 7 — Linear reference vs graph genome

**File**: `diagrams/lecture-11/07-linear-vs-graph.svg`
**Lecture anchor**: §3.3 The graph genome concept
**ViewBox**: `0 0 1080 420`

### Purpose

Side-by-side comparison of the two representational paradigms. Reader leaves able to sketch a pangenome graph given a list of variants.

### Content

**Left panel — Linear reference.** A horizontal 1 kb string with:
- A heterozygous SNP at position 200.
- A small indel at position 500 (insertion of 20 bp).
- A larger SV at position 800 (deletion of 100 bp in one allele).

Each variant shown as an annotation with its ref/alt call "A/G", "del", "20 bp ins" but the reference itself is still a single linear string. Annotations are off to the side.

**Right panel — Pangenome graph.** Same genomic region, but as a graph:
- Shared flanking sequences as long horizontal nodes.
- At the SNP position: two short nodes (one per allele) with two incoming and two outgoing edges.
- At the indel: a branch where one path goes through an insertion node, the other skips.
- At the SV: two parallel nodes (one with the deletable sequence, one without).
- Paths for two haplotypes drawn as coloured traces through the graph: `hap1` in `--accent`, `hap2` in `--warning`.

A mapping line connecting each variant in the linear panel to its corresponding branch in the graph panel.

### Style notes

- Graph nodes: rounded rectangles with node IDs (`N1`, `N2`, …).
- Edges: directed arrows with `--fg-muted` colour; paths highlighted by thick coloured overlays.
- Shared sequence: `--fg-muted` outline.
- Branch points emphasised with colour.

---

## Figure 8 — GFA format anatomy

**File**: `diagrams/lecture-11/08-gfa-anatomy.svg`
**Lecture anchor**: §3.4 GFA — Graphical Fragment Assembly format
**ViewBox**: `0 0 1080 440`

### Purpose

Tie GFA text format to its graph interpretation. Reader leaves able to read a GFA file by eye.

### Content

Two side-by-side panels with connecting lines.

**Left — GFA text.** A syntax-highlighted GFA file of about 10 lines:

```
H  VN:Z:1.0
S  1  ACGTAAGTTTG
S  2  ACGCAAGAATG
S  3  CGATCGATCGA
L  1  +  3  +  0M
L  2  +  3  +  0M
P  hap1  1+,3+  *
P  hap2  2+,3+  *
W  sample1  1  chr1  0  22  >1>3
```

Each line type coloured distinctly: `H` in `--fg-subtle`, `S` in `--accent`, `L` in `--warning`, `P` in `--base-g`, `W` in `--base-u`.

**Right — rendered graph.** The same graph visualised:
- Two parallel `S` nodes (1 and 2) at the top — "ACGTAAGTTTG" and "ACGCAAGAATG".
- One `S` node (3) downstream — "CGATCGATCGA".
- Two `L` edges from nodes 1 and 2 into node 3.
- Path `hap1` traced in green from node 1 → node 3.
- Path `hap2` traced in violet from node 2 → node 3.

Connecting dashed lines from each GFA line to its graph element.

### Style notes

- GFA text in monospace with line-type colour-coding.
- Graph on right: rounded-rectangle nodes with sequence labels.
- Path traces: thick coloured curves overlaying the edges.
- Connecting lines: `--fg-subtle` 0.5px dashed.

---

## Figure 9 — Seed-and-extend on a DAG

**File**: `diagrams/lecture-11/09-graph-seed-extend.svg`
**Lecture anchor**: §4.1 Seed-and-extend generalised to DAGs
**ViewBox**: `0 0 1080 440`

### Purpose

Show the generalisation of linear seed-and-extend to graph references. Reader leaves understanding how alignment branches at graph branches.

### Content

Three bands.

**Top — graph reference.** A DAG with ~8 nodes: three consecutive "trunk" nodes, then a branch into two parallel paths, then merge back into a shared tail. Node IDs visible.

**Middle — the read.** A 30 bp read with a 6-bp k-mer seed highlighted. The seed matches a k-mer in node 3 (the last trunk node before the branch).

**Bottom — extension.** Two extension candidates at the branch point:
- Path A: trunk → node 4 (upper branch) → shared tail.
- Path B: trunk → node 5 (lower branch) → shared tail.

Each candidate's running score displayed. Winner is path A (better score); trace highlighted in `--accent`.

Alignment output shown at the bottom: the read's aligned path as a sequence of (node, offset) tuples.

### Style notes

- Graph: nodes as rectangles; edges with arrow markers.
- Seed: highlighted k-mer cell in `--warning`.
- Extension candidates: parallel dashed traces in two different hues.
- Winning path: `--accent` thick overlay.

---

## Figure 10 — Viterbi on a linear chain vs DAG

**File**: `diagrams/lecture-11/10-viterbi-linear-vs-dag.svg`
**Lecture anchor**: §5.2 Generalising Viterbi to a DAG
**ViewBox**: `0 0 1080 460`

### Purpose

The core EE framing. Show that graph alignment and linear alignment are the same Viterbi algorithm with different state-graph structures.

### Content

Two side-by-side trellis panels.

**Left — linear-chain Viterbi.** A trellis:
- Y-axis: reference state (a chain of ~12 positions).
- X-axis: read time step (~12 columns).
- Vertical arrows between consecutive states show the single predecessor relation.
- One optimal path traced as a thick coloured line through the trellis.

Caption: "linear chain — single predecessor per state, O(N × T) DP".

**Right — graph (DAG) Viterbi.** Same trellis structure but:
- Y-axis: the states form a DAG (branching structure).
- Predecessor relations now include multiple incoming edges to each state (wherever the DAG has in-edges).
- Viterbi path now branches at the DAG branch points — the thick coloured path chooses one branch to go through.

Caption: "DAG — multiple predecessors, O(E × T) DP, same update equation".

**Bottom equation (shared)**: $V(s, t) = \max_{s' \in \text{pred}(s)} [V(s', t-1) + \text{match\_score}(s, x_t)]$

Annotation: "only $\text{pred}(s)$ changes between linear and graph cases".

### Style notes

- Trellis grid lines: `--fg-subtle` 0.3px.
- State nodes: small `--bg-muted` rounded rectangles.
- Predecessor edges: light grey arrows.
- Optimal path: `--accent` 2px thick.
- Equation in Source Serif italic at bottom.

---

## Figure 11 — vg toolchain workflow

**File**: `diagrams/lecture-11/11-vg-toolchain.svg`
**Lecture anchor**: §4.2 The graph-alignment toolchain
**ViewBox**: `0 0 1080 400`

### Purpose

Show the end-to-end vg pipeline from inputs to outputs. Reader leaves able to name each tool and its role.

### Content

Left-to-right pipeline diagram.

**Inputs column (left)**:
- Reference FASTA file icon.
- VCF of known variants, OR
- Collection of assembled genomes (FASTA files).

**Graph-build step**: `vg construct` OR `pggb` → produces `.vg` or `.gfa`.

**Alignment step**: `vg giraffe` (short reads) OR `GraphAligner` (long reads) → produces `.gam` (graph alignment format).

**Variant-call step**: `vg call` → produces VCF.

**Optional projection step**: `vg surject` → produces BAM (for downstream linear-reference tools).

Each arrow labelled with file formats. Each tool in a box with its name and one-line description.

### Style notes

- Boxes: rounded rectangles, `--bg-muted` fill, `--border-strong` outline.
- Arrows: `--accent` 1.5 with standard marker.
- File format labels: JetBrains Mono 9 on the arrows.
- Pipeline stages in a flow line.

---

## Figure 12 — HPRC-v1 pangenome structure

**File**: `diagrams/lecture-11/12-hprc-pangenome.svg`
**Lecture anchor**: §6.2 HPRC v1 — 47 haplotypes
**ViewBox**: `0 0 1080 440`

### Purpose

Show what the HPRC v1 pangenome actually contains. Reader leaves with a concrete sense of population coverage and graph structure.

### Content

**Left — population coverage.** A simple tree showing 6 ancestry groups with the number of haplotypes per group:

- African: 14 haplotypes
- European: 11
- East Asian: 7
- South Asian: 6
- Admixed American: 6
- Oceanian: 3

Tree drawn as a simple two-level dendrogram; leaf nodes labelled with the haplotype counts.

**Right — sample locus in the graph.** A simplified view of a ~2 kb region of the HPRC pangenome graph:
- A GRCh38-like linear "backbone" through the middle.
- Alternate paths branching off at 5 variant sites (3 SNPs, 1 indel, 1 SV).
- Each alternate path labelled with which haplotype groups carry it (e.g. "African ×11 haps", "European ×2 haps", etc.).

### Style notes

- Tree: `--fg-muted` lines with leaf labels in Inter 10.
- Graph: backbone in `--fg-muted` thick, alternate paths in ancestry-specific hues.
- Haplotype-count labels in JetBrains Mono 9.
- Some annotation numbers (total haplotypes = 47) displayed prominently.
