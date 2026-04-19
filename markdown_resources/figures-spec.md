# Lecture 1 — Figures Specification

> **Purpose**: Detailed build specs for static diagrams in Lecture 1. Analogous to `artifacts-spec.md` but for non-interactive figures.
> **Companion files**: `diagram-style-guide.md` (house style), `lecture-01.md` (content + placement), `website-spec.md` (embedding).

---

## Overview

Lecture 1 gets **12 custom SVG diagrams** and **up to 2 real photographs**. This adds visual rhythm without becoming illustrated-textbook busy. The goal is roughly one figure every 15–20 minutes of reading, interleaved with the six interactive artifacts so the page feels varied.

**Figure numbering** is per-lecture, starting at `01`. File names follow `diagrams/lecture-01/NN-name-kebab.svg`.

---

## Figure Placement Summary

| #  | Figure | Section | Type | Rough complexity |
|----|--------|---------|------|------------------|
| 01 | Cell organelles schematic | §1.2 | SVG | Medium |
| 02 | DNA double helix structure | §1.3 | SVG | Medium-high |
| 03 | Central dogma flow | §2.1 | SVG | Low-medium |
| 04 | Genetic code table (compact) | §2.2 | SVG | Medium |
| 05 | End-to-end bioinformatics pipeline | §2.3 | SVG | Low-medium |
| 06 | Sanger chain termination mechanism | §3.2 | SVG | Medium-high |
| 07 | PCR thermal cycle | §3.3 | SVG | Medium |
| 08 | Bridge amplification steps | §4.3 | SVG | High |
| 09 | Adapter-flanked fragment structure | §4.2 | SVG | Low |
| 10 | PacBio ZMW cross-section | §5.2 | SVG | Medium-high |
| 11 | Nanopore + squiggle schematic | §5.3 | SVG | High |
| 12 | FASTQ record anatomy | §6.1 | SVG | Medium (proof file exists) |
| P1 | NovaSeq device (photograph) | §4.6 | Photo | Source & crop only |
| P2 | MinION device (photograph) | §5.3 | Photo | Source & crop only |

---

## 01 — Cell Organelles Schematic

**File**: `diagrams/lecture-01/01-cell-organelles.svg`
**Placement**: §1.2 "The cell"
**Teaching goal**: Visual reference for organelles mentioned in the text. Viewer should be able to locate nucleus, ribosome, mitochondrion, ER/Golgi.

**Composition**:
- A schematic eukaryotic cell (not a real micrograph — a clean line-art simplification)
- Outer membrane as a primary 2.5-weight line
- Nucleus (labeled) with nuclear envelope showing pores
- 4–6 mitochondria scattered, with internal cristae line hatching
- Rough ER as a folded line network near the nucleus
- Golgi apparatus as stacked parallel curves
- A few ribosomes as small filled circles, some on the ER, some free in cytoplasm
- Plasma membrane lipid bilayer abstracted as a double line — don't draw individual phospholipid heads

**Style**:
- Line art, `--fg` default stroke
- Very subtle `--bg-muted` fill on the cytoplasm
- Labels in Inter 11–12, connected to structures with thin `--fg-muted` leader lines

**Do NOT**:
- Draw every organelle — stick to the ones named in the text
- Add color to individual organelles (this is a schematic, not an illustration)
- Include a centriole, peroxisome, lysosome, etc. — not teaching-relevant here

**ViewBox**: ~`0 0 720 480`

---

## 02 — DNA Double Helix Structure

**File**: `diagrams/lecture-01/02-dna-double-helix.svg`
**Placement**: §1.3 "DNA: structure and language"
**Teaching goal**: Show the helical structure, antiparallel strands, base pairing (A-T 2 bonds, G-C 3 bonds), and 5'→3' directionality.

**Composition**:
- **Left panel**: short segment of helix (6–8 base pairs), slightly stylized — two sugar-phosphate backbones as ribbons, base pairs as horizontal rungs
- **Right panel**: "flat ladder" representation of the same segment — base letters shown explicitly, hydrogen bonds shown as dots (2 dots for A-T, 3 for G-C)
- 5' and 3' labels on each strand end, showing antiparallel orientation
- Tiny arrow notation indicating directionality of each strand

**Style**:
- Helix backbones: 2.5-weight cobalt `--accent` curves (two intertwined)
- Base rungs: 1.5-weight `--fg` lines connecting them
- In the flat-ladder panel: base letters in JetBrains Mono 16, colored per base (`--base-a/t/g/c`)
- Hydrogen bond dots: small filled circles in `--fg-muted`
- Labels: Inter 11 for ends; Inter 12 for panel titles

**Do NOT**:
- Render the helix realistically (no 3D cylindrical shading, no atoms shown explicitly)
- Use more than two colors for backbones — they should clearly be "the two strands"

**ViewBox**: ~`0 0 860 380`

---

## 03 — Central Dogma Flow

**File**: `diagrams/lecture-01/03-central-dogma-flow.svg`
**Placement**: §2.1 "The central dogma"
**Teaching goal**: Visualize DNA → RNA → Protein with transcription, splicing, and translation labeled.

**Composition**:
- Three horizontal "track" bands, top to bottom: DNA, mRNA (with introns spliced out), protein
- DNA track: a line with exons as filled amber rectangles and introns as unfilled gaps
- Between DNA and mRNA track: an arrow labeled "Transcription"
- On the mRNA track: show the pre-mRNA (with introns) transitioning to mature mRNA (without), with a "Splicing" label
- Between mRNA and protein: an arrow labeled "Translation"
- Protein track: a string of colored circles representing amino acids, with the first few codons shown above them in JetBrains Mono
- Optional small "Reverse transcription" and "RNA replication" dashed arrows as exceptions, annotated

**Style**:
- Track bands with subtle `--bg-muted` fills
- DNA/mRNA in JetBrains Mono where sequences are shown
- Arrows: `--accent` cobalt, tapered, medium-weight
- Exception arrows: `--fg-muted`, dashed (`stroke-dasharray="4 3"`)

**ViewBox**: ~`0 0 860 380`

---

## 04 — Genetic Code Table (Compact)

**File**: `diagrams/lecture-01/04-genetic-code-table.svg`
**Placement**: §2.2 "The genetic code"
**Teaching goal**: Canonical reference. The reader should be able to look up any codon and see redundancy grouping.

**Composition**:
- Use a **4×4×4 grid** layout (traditional arrangement): first base row, second base column, third base sub-column
- Each cell: codon in JetBrains Mono + 1-letter amino acid in Inter bold + 3-letter abbreviation in Inter 10
- Amino acids colored by chemical category:
  - Hydrophobic (A, V, L, I, M, F, W, P): pale amber tint
  - Polar uncharged (S, T, N, Q, Y, C): pale green tint
  - Acidic (D, E): pale red tint
  - Basic (K, R, H): pale blue tint
  - Special (G): pale gray tint
  - Stop codons: filled `--fg` black
- Legend below the table showing the category colors
- Start codon AUG subtly accented with a cobalt border

**Style**:
- Cell backgrounds use very pale tints (8–10% opacity) so text remains readable
- Grid lines: `--border`, 1px
- Outer border: `--border-strong`, 1.5px
- Codon text: JetBrains Mono, weight 500
- Amino acid 1-letter: Inter, weight 700

**Do NOT**:
- Use strong saturated colors for cell fills — they're backgrounds, not the main information
- Include the IUPAC nomenclature in Greek — simplify to three-letter + one-letter codes

**ViewBox**: ~`0 0 720 480`

---

## 05 — End-to-End Bioinformatics Pipeline

**File**: `diagrams/lecture-01/05-bioinformatics-pipeline.svg`
**Placement**: §2.3 "End-to-end bioinformatics pipeline"
**Teaching goal**: The big-picture map. This figure is referenced multiple times throughout the lecture, so it must be readable at a glance.

**Composition**:
- Six boxes in a left-to-right flow: Sample → Library Prep → Sequencer → Raw Reads (FASTQ) → QC & Trimming → Alignment/Assembly → Variant Calling → Biological Interpretation
- Arrows between boxes
- Boxes are slightly different widths to fit their labels; all aligned on a baseline
- Below each box: a thin caption in `--fg-muted` listing a key detail (e.g., under "Sequencer": "Illumina / PacBio / ONT")
- Annotation band above the flow showing which boxes are covered in Lecture 1 (the first four) vs. "subsequent lectures" (the last three) — use a subtle horizontal bracket

**Style**:
- Boxes with `--bg-muted` fill, `--border-strong` 1.5px stroke, 6px rounded corners
- Primary-flow arrows: `--accent`, tapered
- The "covered in this lecture" bracket: `--accent` color, thin 1px, with a small tick on each end
- Labels inside boxes: Inter, weight 500, size 13

**ViewBox**: ~`0 0 900 280`

---

## 06 — Sanger Chain Termination Mechanism

**File**: `diagrams/lecture-01/06-sanger-termination.svg`
**Placement**: §3.2 "Sanger sequencing"
**Teaching goal**: Show how ddNTPs cause termination and produce a ladder of fragments, each ending in a different colored base.

**Composition**:
- **Top panel**: The biochemistry
  - A template strand, a primer, and a growing complementary strand
  - Incoming dNTPs shown as small squares near the polymerase
  - A ddNTP incorporated at the end — visually distinct (filled marker indicating the missing 3' OH)
  - Callout: "no 3'-OH → chain cannot extend"
- **Bottom panel**: The result
  - A gel/capillary ladder showing fragments of increasing length, each terminating in a colored base (A/C/G/T)
  - The sequence read off the bottom to top
- Arrow connecting top to bottom: "fluorescent ladder → sequence"

**Style**:
- Template and synthesized strand in JetBrains Mono where bases are labeled
- Base colors `--base-a/t/g/c`
- The ddNTP marker: a small cobalt asterisk or filled dot annotation
- Gel ladder: horizontal bands separated by a clean gray grid

**ViewBox**: ~`0 0 860 520`

---

## 07 — PCR Thermal Cycle

**File**: `diagrams/lecture-01/07-pcr-thermal-cycle.svg`
**Placement**: §3.3 "PCR"
**Teaching goal**: Show the three steps and the doubling-per-cycle exponential growth.

**Composition**:
- **Top**: a circular/cyclical diagram with three labeled arcs: Denature (95°C), Anneal (55–65°C), Extend (72°C). Each arc shows a cartoon of the DNA state (separated / primers bound / extended).
- **Bottom-left**: a mini temperature-vs-time plot (sawtooth) for 3 cycles, showing the temperature ramps
- **Bottom-right**: a "molecule count" panel — 1, 2, 4, 8, 16, …, 2³⁰ ≈ 10⁹ — illustrated as doubling icons or a simple exponential bar

**Style**:
- Circular flow: `--accent` curved arrows between stages, stages labeled in Inter 12
- DNA cartoons: the same helix motif from Figure 02 but very simplified (2–3 bp segments)
- Temperature plot axes in `--fg-muted` 10
- Doubling count in JetBrains Mono

**ViewBox**: ~`0 0 860 500`

---

## 08 — Bridge Amplification Steps

**File**: `diagrams/lecture-01/08-bridge-amplification.svg`
**Placement**: §4.3 "Cluster generation — bridge amplification"
**Teaching goal**: Show the stepwise bridge amplification that turns a single fragment into a cluster on the flow cell.

**Composition**:
- 5–6 panels in a horizontal strip, each showing a stage:
  1. Surface with oligos (P5 and P7 kinds, two colors)
  2. Library fragment hybridizing by its P5 end
  3. First extension: second strand synthesized covalently to surface oligo
  4. Bridge formation: free end bends over to a nearby P7 oligo
  5. Bridge amplification: extension produces two surface-attached strands
  6. After ~35 cycles: a cluster of ~1000 molecules on a single spot
- Each panel has a small caption with step number and name
- Final panel's "cluster" drawn as a small filled semi-transparent blob with dozens of fine line traces

**Style**:
- Surface: horizontal gray bar with small vertical "oligo" ticks
- P5 and P7 oligos in two distinct colors (use `--accent` cobalt for one and `--fg-muted` for the other — don't introduce a new color)
- Fragment strands: `--base-a` red and `--base-c` green respectively (complementary strands)
- Each panel has a thin `--border` rectangle around it
- Step number kicker: Inter 11, weight 600, `letter-spacing="0.08em"`, `--accent`

**ViewBox**: ~`0 0 980 360` (this one is wider to fit the panel strip)

---

## 09 — Adapter-Flanked Fragment Structure

**File**: `diagrams/lecture-01/09-adapter-structure.svg`
**Placement**: §4.2 "Library preparation"
**Teaching goal**: Replace the ASCII diagram in the markdown with a clean SVG version.

**Composition**:
- A single horizontal strand drawn as a wide pill shape
- Divided into 7 segments (labeled and colored differently):
  P5 adapter | i5 index | Read 1 primer | INSERT (the genomic fragment) | Read 2 primer | i7 index | P7 adapter
- The INSERT segment is notably wider than the others and uses a `--bg-muted` fill to distinguish it from the functional adapter segments
- 5' and 3' labels on the ends
- Segment labels below each segment, in Inter 11
- Optional: small arrows indicating where Read 1 and Read 2 sequencing primers begin extending from

**Style**:
- Functional adapter segments: each a distinct pale color (pastel versions of the accent and base palette; no saturated colors)
- INSERT segment: `--bg-muted` fill, emphasizing it's "the sample"
- All segments bordered by a thin `--border-strong` line
- Primer arrows: `--accent`, thin, tapered

**ViewBox**: ~`0 0 980 200`

---

## 10 — PacBio ZMW Cross-Section

**File**: `diagrams/lecture-01/10-pacbio-zmw.svg`
**Placement**: §5.2 "PacBio SMRT sequencing"
**Teaching goal**: Show the zero-mode waveguide concept — a sub-wavelength hole, evanescent illumination, polymerase at the bottom.

**Composition**:
- A cross-sectional view of the ZMW structure:
  - Glass substrate below (horizontal line with hatching)
  - Metal film above it with a sub-wavelength hole drilled through (the "ZMW")
  - Polymerase molecule fixed at the glass surface at the bottom of the hole
  - Incoming laser light from below (arrow with wavy line indicating wavelength)
  - **Evanescent field**: a graded fill (use a series of horizontal rectangles with decreasing opacity rather than a gradient — gradients are forbidden, but we can approximate) extending ~20–30 nm up the hole
  - **Outside the illuminated region**: noted as "dark zone"
- Nucleotide molecules shown in the bulk solution (dark zone) and one being incorporated at the polymerase (in the illuminated zone)
- Labels: "zero-mode waveguide", "evanescent field", "polymerase", "fluorescent dNTPs"
- Small inset: "cutoff wavelength physics" — a tiny note saying visible light (λ ~ 500nm) can't propagate through a ~70nm hole

**Style**:
- Metal film: dark gray `--fg-muted` fill with fine horizontal hatching
- Glass: `--bg-muted` fill
- Evanescent field: stepped rectangles of `--accent-bg` with decreasing opacity (1.0, 0.6, 0.3, 0.1)
- Polymerase: a simple rounded-rectangle abstraction labeled "pol"
- Nucleotides: small filled circles in base colors
- Light arrow: `--accent`

**ViewBox**: ~`0 0 720 520`

---

## 11 — Nanopore + Squiggle Schematic

**File**: `diagrams/lecture-01/11-nanopore-squiggle.svg`
**Placement**: §5.3 "Oxford Nanopore"
**Teaching goal**: Show the pore, the DNA threading through, and the ionic current trace it produces.

**Composition**:
- **Left side** (40% of width): Schematic cross-section of the nanopore
  - Membrane as a horizontal band with lipid-bilayer abstraction (two thin parallel lines)
  - Protein pore embedded in the membrane — drawn as a wider-than-constriction shape with an obvious narrow throat
  - DNA strand threading through, with 5 bases highlighted in the constriction ("5 bases in the pore simultaneously")
  - Motor protein sitting atop the pore (labeled)
  - Voltage indicator (+ above, − below) with an arrow showing ion flow
- **Right side** (60% of width): An ionic current trace (squiggle) scrolling left-to-right
  - Y-axis: current in pA (labeled ~60–120)
  - X-axis: time
  - Stepped current levels clearly visible (each level = a k-mer)
  - Below the trace: small labels showing which k-mer produced each level
  - A vertical line/marker connects the "current position" in the squiggle to the bases in the pore on the left side

**Style**:
- Pore protein: `--fg-muted` stroke, `--bg-muted` fill
- DNA bases in the pore: JetBrains Mono 12, colored per base
- Current trace: `--accent` cobalt line, 1.5-weight
- k-mer labels below trace: JetBrains Mono 10, `--fg-muted`
- Axis labels: Inter 10, `--fg-subtle`

**ViewBox**: ~`0 0 980 400`

---

## 12 — FASTQ Record Anatomy

**File**: `diagrams/lecture-01/12-fastq-anatomy.svg`
**Placement**: §6.1 "The 4-line record"
**Teaching goal**: Same as the proof file — annotated FASTQ record + Phred score decoding zoom panel.

**Status**: **Already built as proof file** (`../fastq-anatomy-proof.html`). Extract the `<svg>` block from the proof, save as standalone `.svg`, and place in `diagrams/lecture-01/12-fastq-anatomy.svg`.

---

## P1 — NovaSeq Device Photograph

**File**: `diagrams/lecture-01/photos/novaseq-device.jpg`
**Placement**: §4.6 "Platforms & use cases"
**Teaching goal**: One real photograph so students see what an industrial sequencer actually looks like.

**Source**: Look for a CC-licensed image on Wikimedia Commons ("Illumina NovaSeq") or request one from Illumina's press kit. If neither is suitable, use a MiSeq photograph instead (more commonly photographed).

**Crop / treatment**:
- Straighten and crop to 3:2 aspect
- Neutral background if possible; if the source shows a lab environment, keep only the device
- No color filtering — real photograph

**Caption**: "Illumina NovaSeq X. For scale, the instrument is roughly the size of a large chest freezer." — plus attribution.

---

## P2 — MinION Device Photograph

**File**: `diagrams/lecture-01/photos/minion-device.jpg`
**Placement**: §5.3 "Oxford Nanopore" — paired with the squiggle diagram (one shows the schematic, one shows the real thing)
**Teaching goal**: Emphasize the dramatic form-factor of nanopore — a full sequencer the size of a USB stick.

**Source**: ONT has press-kit photographs of the MinION. Wikimedia Commons also has some. Need CC-licensed or press-use terms.

**Crop / treatment**:
- Show the MinION with something for scale — ideally held in a hand, or beside a standard object
- Neutral background

**Caption**: "Oxford Nanopore MinION. A complete sequencer in a USB-powered device." — plus attribution.

---

## Build Order Recommendation

1. **12 (FASTQ anatomy)** — already done as proof, just extract.
2. **05 (Pipeline)** — simple, high-value, referenced repeatedly.
3. **09 (Adapter structure)** — simple, replaces the ASCII diagram.
4. **03 (Central dogma)** — simple-medium, pairs with interactive artifact #2.
5. **07 (PCR cycle)** — medium, visually distinctive.
6. **02 (Double helix)** — medium-high.
7. **01 (Cell organelles)** — medium.
8. **04 (Genetic code table)** — medium, detailed but systematic.
9. **06 (Sanger termination)** — medium-high.
10. **11 (Nanopore squiggle)** — high, pairs with interactive artifact #5.
11. **10 (ZMW cross-section)** — high.
12. **08 (Bridge amplification)** — highest complexity, 6 panels in sequence.
13. **P1, P2** — source and crop photos once all SVGs are done.

## Testing checklist (per figure)

- [ ] Opens standalone in browser, no external dependencies
- [ ] All colors from design-token palette only (no stray hex values)
- [ ] Only Inter / JetBrains Mono / Source Serif 4 in `font-family`
- [ ] `role="img"`, `<title>`, `<desc>` present
- [ ] Readable at ~720px display width
- [ ] Bases (if shown) use `--base-*` colors correctly and consistently
- [ ] Passes the diagram-style-guide §14 checklist
