# Lecture 4 — Figures Specification

> **Scope**: Static diagrams for Lecture 4 (Variant Calling).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-04.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. Variant calling is algorithmic; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Line-art first; fill only where the fill earns the ink.
- Base colors (`--base-a/t/g/c`) only where actual nucleotide characters appear. Otherwise use the neutral palette.
- DNA sequences, VCF field values, and codon labels use **JetBrains Mono**; all other labels use **Inter**.
- Arrows use the shared `<marker id="arrow-accent">` pattern.
- Escape `&` in text content as `&amp;` (see `diagram-style-guide.md` §6).

## Figure Budget

Eleven figures for a ~3h 30min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Variant taxonomy and size spectrum | Part 1 | Custom SVG |
| 2 | SNV consequences in coding regions | Part 1 | Custom SVG |
| 3 | Frameshift vs in-frame INDEL | Part 1 | Custom SVG |
| 4 | Anatomy of a samtools pileup | Part 2 | Custom SVG |
| 5 | Anatomy of a variant call at one position | Part 2 | Custom SVG |
| 6 | INDEL misalignment and local realignment | Part 3 | Custom SVG |
| 7 | Bayesian genotyping at a position | Part 3 | Custom SVG |
| 8 | The three SV detection signals | Part 4 | Custom SVG |
| 9 | The five canonical SV types | Part 4 | Custom SVG |
| 10 | Anatomy of a VCF file | Part 5 | Custom SVG |
| 11 | Variant annotation pipeline | Part 5 | Custom SVG |

---

## Figure 1 — Variant Taxonomy and Size Spectrum

**File**: `diagrams/lecture-04/01-variant-taxonomy.svg`
**Lecture anchor**: §1.2 Classification by size and structure
**ViewBox**: `0 0 900 360`

### Purpose

Establish the full size scale of genomic variation in one image, so SNVs, INDELs, SVs, and CNVs live in a shared frame from the first reference. Readers should leave knowing the approximate size range of each class and the approximate per-genome count.

### Content

A horizontal log-scale size axis from 1 bp to 10 Mb running across the middle of the figure (tick marks at 1 bp, 10 bp, 100 bp, 1 kb, 10 kb, 100 kb, 1 Mb, 10 Mb).

Four bins drawn as coloured horizontal bands above the axis, each labeled and sized to its typical length range:

1. **SNVs** — 1 bp. Single-position tick; per-genome count label "~3,000,000".
2. **INDELs** — 1 to 50 bp. Accent-coloured band covering that range; per-genome count "~500,000".
3. **SVs** — 50 bp to 10 Mb. Wider accent-coloured band; count "~25,000 per haplotype (long-read)".
4. **CNVs** — 1 kb to 10 Mb. Overlapping band shown at slightly different vertical position to convey that CNVs are a sub-class of SVs; count "~100 large events per genome".

Below the axis, annotation bubbles for each bin naming their primary functional consequences:

- SNV bin → "missense / nonsense / silent / splice"
- INDEL bin → "frameshift / in-frame"
- SV bin → "gene disruption / fusion / regulatory"
- CNV bin → "gene-dose changes"

### Style notes

- Log axis ticks in JetBrains Mono 10 `--fg-subtle`.
- Bins in soft `--accent-bg` fill with `--accent` outline.
- Per-genome counts in Inter 11 `--fg`, bold.
- Consequence annotations in Inter 10 `--fg-muted`.
- No gradient or glow; flat bands.

---

## Figure 2 — SNV Consequences in Coding Regions

**File**: `diagrams/lecture-04/02-snv-consequences.svg`
**Lecture anchor**: §1.3 SNVs and their consequences
**ViewBox**: `0 0 860 420`

### Purpose

Let the reader see at a glance how a single-base change produces four different consequence classes depending on which codon it lands in.

### Content

**Top band — reference.** A short mRNA segment shown as a sequence of three adjacent codons in JetBrains Mono, with each base coloured by the `--base-*` palette. Example: `ATG-CCT-GAG` (Met-Pro-Glu). Underneath the codons, the corresponding amino acids in one-letter code.

**Below — four parallel panels**, one per consequence class. Each panel shows the same reference codon with one specific base changed, and the resulting codon + amino acid:

1. **Silent** — e.g. `CCT` (Pro) → `CCC` (Pro). Base change highlighted in `--accent`, amino acid labeled unchanged, caption "silent: same amino acid."
2. **Missense** — e.g. `CCT` (Pro) → `CCG` (Ala? — or verify with spec; use an actually-missense example). Base change highlighted; old vs new amino acid shown side-by-side.
3. **Nonsense** — e.g. `CCT` (Pro) → `CAT`... (use an actually-nonsense example: e.g., `TGG` → `TAG` is Trp → STOP). Base change highlighted; amino acid label replaced with "STOP" and a small "×" to indicate truncation.
4. **Splice-region** — a codon at an exon-intron boundary, with the boundary drawn as a vertical dashed line. A SNV at the conserved GT donor or AG acceptor disrupts splicing; show the exon/intron boundary, mark the disrupted splice site, and caption "splice disruption: may skip exon or retain intron."

### Style notes

- Codons in JetBrains Mono 16, base-coloured.
- Amino acids in Inter 12 weight 500.
- SNV base highlighted with a soft `--accent-bg` fill behind the letter.
- STOP codon drawn with `--error` text color and a small "⊗" glyph after.
- Panel kickers in small-caps Inter 11 weight 600 `--accent`.

---

## Figure 3 — Frameshift vs In-Frame INDEL

**File**: `diagrams/lecture-04/03-indel-frameshift.svg`
**Lecture anchor**: §1.4 INDELs and frameshifts
**ViewBox**: `0 0 860 360`

### Purpose

Make the arithmetic reason for frameshifts visceral — a 3-base deletion spares the reading frame downstream; a 1-base deletion destroys it.

### Content

**Top band — reference.** A 24-base mRNA segment drawn in JetBrains Mono, with codon boundaries marked by small vertical ticks every 3 bases, and the resulting amino-acid sequence below (eight amino acids).

**Middle band — in-frame deletion (3 bp).** Same sequence with three adjacent bases removed. Codon boundaries re-drawn; the deletion region highlighted; the resulting amino-acid sequence shows one missing amino acid, with the rest unchanged.

**Bottom band — frameshift deletion (1 bp).** Same sequence with one base removed. Codon boundaries re-drawn from the deletion onward; the resulting amino-acid sequence is completely different from the deletion point forward. Mark a premature stop codon somewhere downstream with a clear "STOP" label and a note: "truncated protein, usually non-functional."

### Style notes

- Reference bases in JetBrains Mono 14, base-coloured.
- Codon ticks in `--fg-muted` at y-offset above bases.
- Deletion region drawn as a strikethrough on the missing bases plus a red "×" in `--error`.
- Amino acid letters in Inter 12 weight 500.
- Downstream scrambled region drawn with a soft `--warning` tint behind the bases.

---

## Figure 4 — Anatomy of a samtools Pileup

**File**: `diagrams/lecture-04/04-pileup-anatomy.svg`
**Lecture anchor**: §2.2 The pileup
**ViewBox**: `0 0 920 400`

### Purpose

Decompose a single pileup line into its column-by-column meaning. The reader should be able to read a real pileup after seeing this figure.

### Content

**Top band — a real-looking pileup line** in JetBrains Mono, with columns visually separated by small gaps:

```
chr1    1000    A    12    ,,..,A,..,aA,    IIIIII!IIIII
```

Each column has a leader line dropping to a labeled callout below it:

- Column 1 → "Chromosome (reference contig name)"
- Column 2 → "Position (1-based)"
- Column 3 → "Reference base"
- Column 4 → "Depth: 12 reads overlap this position"
- Column 5 → "Read bases: `.`,`,` = matches REF; uppercase/lowercase letters = non-REF (case encodes strand)"
- Column 6 → "Base qualities (Phred ASCII); one char per read in col 5"

**Below — an inset "special characters" legend.** A small panel showing the extra tokens:

- `^]` — "start of read, followed by mapping-quality character"
- `$` — "end of read"
- `+3AAA` — "insertion at next position of 3 bases AAA"
- `-2TT` — "deletion at next position of 2 bases TT"

### Style notes

- Pileup line background tint `--bg-muted`.
- Leader lines in `--fg-subtle` 1 px dashed.
- Callouts in Inter 11 weight 400 `--fg-muted`.
- Inset panel bordered in `--border-strong`, background `--bg-inset`.
- Non-ref bases in the example row coloured using `--base-*` tokens to be consistent with other figures.

---

## Figure 5 — Anatomy of a Variant Call at One Position

**File**: `diagrams/lecture-04/05-variant-record-anatomy.svg`
**Lecture anchor**: §2.3 The vocabulary of a variant record
**ViewBox**: `0 0 860 420`

### Purpose

Tie the pileup-level evidence to the populated VCF-level fields. Every field in the vocabulary (REF, ALT, DP, AD, VAF, GT, GQ, QUAL) should be visible in one place, with its arithmetic derivation.

### Content

**Top — pileup column snapshot.** A simulated "zoom-in" of reads piling up on one reference position. Draw 17 small horizontal read bars stacked vertically, each labeled with the base it carries at the focus position (12 "A" in `--base-a`, 5 "G" in `--base-g`). All at the same x-coordinate where the reference base is labeled below.

**Middle — derivation arrows.** From the pileup column, arrows drop into a table:

| Field | Value | How computed |
|---|---|---|
| REF | A | reference base |
| ALT | G | most-common non-REF base |
| DP | 17 | total reads |
| AD | 12, 5 | reads for REF, reads for ALT |
| VAF | 0.294 | AD[1] / DP = 5 / 17 |
| GT | 0/1 | heterozygous (VAF near 0.5) |
| GQ | 95 | Phred-scaled genotype confidence |
| QUAL | 420 | Phred-scaled variant-site quality |

**Bottom — the VCF-row equivalent.** A single VCF line rendered in JetBrains Mono showing how the above fields populate the record:

```
chr1  1000  .  A  G  420  PASS  DP=17;AF=0.294  GT:AD:DP:GQ  0/1:12,5:17:95
```

### Style notes

- Reads in `--accent` outline, base letter inside each coloured by `--base-*`.
- Derivation table with thin `--border` lines between cells, no heavy grid.
- VCF line at bottom on `--bg-inset` background.
- Connecting arrows from pileup to table, from table to VCF line.

---

## Figure 7 — Bayesian Genotyping at a Position

**File**: `diagrams/lecture-04/07-bayesian-genotyping.svg`
**Lecture anchor**: §3.2 Caller families — Bayesian
**ViewBox**: `0 0 860 440`

### Purpose

Make the Bayesian inference explicit: prior × likelihood = posterior, evaluated over the three diploid genotypes. The reader should leave knowing what the caller is actually computing.

### Content

**Top — a minimal pileup.** Shown as 10 reads: 6 REF "A", 4 ALT "G", all at Q30. Small tight figure, top-left.

**Middle — three parallel likelihood calculations**, one per genotype:

- `P(D | 0/0)` — probability of observing 4 "G" reads assuming true genotype is homozygous reference. Worked out as `(1−ε)^6 · (ε/3)^4` where ε = 10⁻³ (Q30). Compute the number; show it scientific notation.
- `P(D | 0/1)` — probability assuming heterozygous. `(0.5)^10` modified by quality. Show the number.
- `P(D | 1/1)` — probability assuming homozygous alternate. Very low, like `(ε/3)^6 · (1−ε)^4`.

Each likelihood shown as a small formula and a numerical value.

**Below — prior row.** A small bar chart of `P(G)` from a typical human prior: `P(0/0) ≈ 0.999`, `P(0/1) ≈ 7·10⁻⁴`, `P(1/1) ≈ 2·10⁻⁴`. Three bars in `--bg-muted` with `--accent` outline.

**Below — posterior row.** `P(G | D) ∝ P(D | G) · P(G)`. Three bars showing the normalised posterior; the MAP (maximum-a-posteriori) genotype highlighted with a thick `--accent` border.

**Side annotation bubble.** "Caller picks the MAP genotype. GQ = Phred-scaled ratio of MAP to runner-up."

### Style notes

- Formulas in JetBrains Mono 12 for consistency with other math.
- Bars in `--bg-muted` + `--accent` outline; MAP winner with thick `--accent` border and `--accent-bg` fill.
- Connecting arrows from likelihoods × prior → posterior in `--fg-muted` 1 px.

---

## Figure 6 — INDEL Misalignment and Local Realignment

**File**: `diagrams/lecture-04/06-indel-realignment.svg`
**Lecture anchor**: §3.1 Pre-calling pipeline
**ViewBox**: `0 0 860 420`

### Purpose

Show why INDEL realignment is a real, distinct problem: the same true 2-bp deletion gets expressed two different ways by the aligner, and only one of the two is correct.

### Content

**Top band — reference + reads, naive alignment.** A reference line with the label "REFERENCE" running horizontally. Six reads aligned below it, all carrying a true 2-bp deletion near the middle. In the naive alignment, most reads have been aligned as-if-ungapped — producing a cluster of 5–6 mismatched bases at the deletion site instead of a single gap. Mismatches coloured `--error`; the true deletion region drawn with a dashed outline above the reference.

Annotation: "without local realignment: deletion looks like a cluster of SNVs."

**Middle — arrow "local realignment" pointing down.**

**Bottom band — after local realignment.** Same reads, now aligned with gaps properly opened at the deletion site. All reads cleanly show a 2-bp gap in the same location; no mismatches in the region. Gap characters drawn as small horizontal dashes in `--fg-muted`.

Annotation: "with local realignment: clean 2-bp deletion, no false SNVs."

### Style notes

- Reference line in `--fg` 2.5px.
- Reads in `--accent` 1.5px outline.
- Mismatched bases highlighted `--error` bold; gap characters in `--fg-muted`.
- Before/after bands visually separated by a horizontal rule.

---

## Figure 9 — The Five Canonical SV Types

**File**: `diagrams/lecture-04/09-sv-types.svg`
**Lecture anchor**: §4.3 SV types and naming
**ViewBox**: `0 0 960 400`

### Purpose

One-stop reference for the standard SV vocabulary. Reader should come away with a mental image for each type.

### Content

**Five equal side-by-side panels**, each with a reference segment (or two, for translocation) on top and the SV-altered version below.

1. **Deletion (DEL).** Reference: `—A—B—C—D—`. Altered: `—A——D—` (B and C absent). Caption: "Segment removed."
2. **Insertion (INS).** Reference: `—A—B—C—`. Altered: `—A—X—B—C—` where X is a new inserted segment (shown in `--accent-bg`). Caption: "Segment added."
3. **Inversion (INV).** Reference: `—A—B—C—D—`. Altered: `—A—⟵C⟵B⟵—D—` (B-C reversed). Caption: "Segment flipped."
4. **Duplication (DUP, tandem).** Reference: `—A—B—C—`. Altered: `—A—B—B—C—` (B repeated). Caption: "Segment duplicated."
5. **Translocation (TRA / BND).** Two chromosomes shown:
   - Chromosome 1 reference: `—A—B—C—`. After event: `—A—Y—Z—`.
   - Chromosome 2 reference: `—X—Y—Z—`. After event: `—X—B—C—`.
   - Arrows showing the reciprocal swap. Caption: "Segments exchanged between chromosomes."

### Style notes

- Reference segments drawn as horizontal bars in `--bg-muted` with `--border-strong` outline and a small capital letter label for each segment.
- Altered segments in the same style but with affected regions highlighted in `--accent-bg` or `--warning-bg` to signal the change.
- Panel kickers in small-caps Inter 11 weight 600 `--accent`.
- Brief plain-language caption under each panel in Inter 10 `--fg-muted`.

---

## Figure 8 — The Three SV Detection Signals

**File**: `diagrams/lecture-04/08-sv-detection-signals.svg`
**Lecture anchor**: §4.1 Why SVs need a different algorithm
**ViewBox**: `0 0 900 480`

### Purpose

Show the three signal channels SV callers exploit — all reporting on the same 3 kb deletion event, but each providing different, complementary information.

### Content

A single reference region drawn horizontally across the top, with a true 3 kb deletion marked as a grey striped band between positions X and Y.

Three evidence tracks stacked below, each annotated with its caller-mechanism label:

- **Track 1: Discordant read pairs.** Five read pairs drawn as arcs connecting a left-read (on one side of the deletion) with its mate (on the far side), with insert-size annotations like "mapped insert = 3,400 bp · expected 400 bp → discordant." The mismatched arcs cluster around the breakpoints.
- **Track 2: Split reads.** Three single reads drawn, each half-mapping to the left flanking region and half-mapping to the right flanking region, with a clear visual break at the breakpoint. Annotation: "single-base-resolution breakpoint."
- **Track 3: Read depth.** A coverage-profile curve running under the reference showing flat coverage at 30× outside the deletion, dropping to ~0 between the breakpoints, then returning to 30× on the far side. Annotation: "depth drop = homozygous deletion; 0.5× drop = heterozygous."

### Style notes

- Reference bar: `--bg-muted` fill, `--fg-muted` stroke, deleted region overlaid with diagonal stripe pattern in `--warning-bg`.
- Track labels in small-caps Inter 11 weight 600 `--accent`.
- Read-pair arcs in `--accent` thin stroke; split reads with a visible break gap (small crossed-out segment in the middle).
- Depth curve filled under in `--bg-muted`, outline in `--fg`.

---

## Figure 10 — Anatomy of a VCF File

**File**: `diagrams/lecture-04/10-vcf-anatomy.svg`
**Lecture anchor**: §5.1 The VCF file format
**ViewBox**: `0 0 960 460`

### Purpose

Decompose a small multi-line VCF into its structural components: header meta-lines, column-header line, and data rows. Every field of a data row should be reachable via an annotated callout.

### Content

A rendered VCF fragment in JetBrains Mono 11 taking the main area of the figure:

```
##fileformat=VCFv4.3                                    ← meta: version
##reference=GRCh38                                      ← meta: reference
##INFO=<ID=DP,Number=1,Type=Integer,Description=...>    ← meta: INFO declaration
##FORMAT=<ID=GT,Number=1,Type=String,...>               ← meta: FORMAT declaration
#CHROM  POS    ID  REF  ALT  QUAL  FILTER  INFO      FORMAT       SAMPLE1
chr1    1000   .   A    G    420   PASS    DP=17     GT:AD:DP:GQ  0/1:12,5:17:95
chr1    1500   .   C    T    612   LowQual DP=8      GT:AD:DP:GQ  0/0:7,1:8:22
```

Labeled callouts pointing at each column header and each data row's fields. Special attention to:

- The header-line `##FORMAT=<...>` declaration and how it drives parsing of the per-sample colon-separated value.
- The difference between `.` (missing ID) and `rs12345` (dbSNP ID, not shown in this example but mentioned in a small note).
- The per-sample column parsing: `GT=0/1, AD=12,5, DP=17, GQ=95`.

### Style notes

- Meta-lines shaded `--bg-inset`.
- Column-header line bold, `--accent`.
- Data rows with alternating faint row backgrounds (`--bg-muted` every other row).
- Callouts with thin leader lines in `--fg-subtle`.
- Inset at corner labeled "VCF is indexable with tabix: `bgzip file.vcf && tabix -p vcf file.vcf.gz`" — one-line note in `--fg-muted`.

---

## Figure 11 — Variant Annotation Pipeline

**File**: `diagrams/lecture-04/11-annotation-pipeline.svg`
**Lecture anchor**: §5.3 Variant annotation
**ViewBox**: `0 0 1020 340`

### Purpose

Trace a single variant through the annotation pipeline, showing what each stage adds to the record.

### Content

A horizontal flow diagram: one variant enters on the left, with the bare-minimum VCF fields. At each stage, new annotation fields are added.

**Input VCF (left):** `chr7:140753336 A>T QUAL=520 DP=65 GT=0/1`

**Stage 1 — VEP / snpEff.** Adds gene and consequence annotation. Output adds `Gene=BRAF, Consequence=missense, ProteinChange=p.Val600Glu, Transcript=ENST00000288602`.

**Stage 2 — gnomAD join.** Adds population frequency. Output adds `gnomAD_AF=1.2e-5, gnomAD_AF_popmax=3e-5`.

**Stage 3 — ClinVar join.** Adds clinical significance. Output adds `ClinVar_Significance=Pathogenic, ClinVar_DiseaseName=melanoma, ClinVar_ID=VCV000013961`.

**Stage 4 — Functional prediction (REVEL, CADD, SpliceAI).** Adds numeric impact scores. Output adds `REVEL=0.932, CADD=26.5, SpliceAI=0.01`.

**Output VCF (right):** Same variant with all annotations present. Below the terminal box: a verdict banner reading "BRAF V600E · Pathogenic · Rare in population · High functional impact."

### Style notes

- Stages drawn as horizontal boxes with `--bg-muted` fill and `--border-strong` 1.5 px stroke.
- Annotations stack vertically inside each box, with new additions highlighted in `--accent`.
- Cobalt arrows between stages.
- Verdict banner at end drawn in a soft `--warning-bg` (since Pathogenic warrants attention), with `--warning` stroke.

---

## Cross-Figure Consistency Notes

- **SNV bases and codons** appear in Figures 2, 3 and should use the same `--base-*` palette.
- **Pileup-column representation** is shared between Figures 4 and 5 — use the same read-bar style.
- **SV diagrams** in Figures 8 and 9 must use the same segment-labeling convention (capital letter segments like `A`, `B`, `C`).
- **VCF rendering** in Figures 5, 10 should use identical JetBrains Mono sizing and column alignment.
- **Kicker labels** (uppercase, small-caps Inter 11 weight 600, letter-spacing 0.12em, colour `--accent`) are used for section headers inside every figure.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All eleven figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient, drop shadow, glow, or 3D effect.
- [ ] All sequences, codons, and VCF fields are in JetBrains Mono; all other labels in Inter.
- [ ] Base colors appear only where actual nucleotides are shown (Figures 2, 3, possibly 4 and 5).
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
- [ ] All `&`, `<`, `>` in text content are XML-escaped.
