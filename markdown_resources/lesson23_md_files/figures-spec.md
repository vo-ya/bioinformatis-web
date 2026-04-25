# Lecture 23 (proposed L18) — Figures Specification

> **Scope**: Static diagrams for Lecture 23 (Metagenomics and the Microbiome).
> **Companion files**: `lecture-23.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Phylum colours: Firmicutes amber, Bacteroidetes cobalt, Proteobacteria red, Actinobacteria green, Verrucomicrobia teal, others muted.
- Diversity-metric colours: alpha (within-sample) violet, beta (between-sample) cobalt.
- 16S vs shotgun colours: 16S amber, shotgun cobalt.
- Typography: Inter for chrome; JetBrains Mono for taxonomic names, sequences, abundances.

## Figure budget — 12 figures

| # | Title | Part |
|---|---|---|
| 1 | 16S vs shotgun metagenomics decision | Part 1 |
| 2 | DADA2 16S workflow | Part 2 |
| 3 | Shotgun metagenomics pipeline | Part 3 |
| 4 | Diversity metrics: alpha and beta | Part 4 |
| 5 | Microbiome-disease association case study | Part 5 |
| 6 | Microbiome analysis workflow | Part 6 |
| 7 | 16S rRNA gene structure (variable regions) | Part 2 |
| 8 | Compositional data problem | Part 4 |
| 9 | MAG recovery from contigs | Part 3 |
| 10 | UniFrac phylogenetic distance | Part 4 |
| 11 | Strain-level resolution: short vs long reads | Part 7 |
| 12 | Antibiotic resistance gene tracking | Part 7 |

---

## Figure 1 — 16S vs shotgun metagenomics decision

**File**: `diagrams/lecture-23/01-16s-vs-shotgun.svg`
**ViewBox**: `0 0 1200 600`

Two side-by-side decision panels: 16S (amber border) and shotgun (cobalt border). Each lists:

- Cost: 16S $30/sample, shotgun $300/sample.
- Resolution: 16S genus, shotgun species/strain.
- Taxa covered: 16S bacteria-only, shotgun all kingdoms.
- Functional content: 16S no, shotgun yes.
- Bioinformatics complexity: 16S simple, shotgun complex.

Bottom decision tree: question pyramid leading to choice. "Is functional content needed?" → if yes → shotgun. "Is budget tight?" → if yes → 16S. "Is sample biomass low?" → if yes → 16S to avoid host contamination.

---

## Figure 2 — DADA2 16S workflow

**File**: `diagrams/lecture-23/02-dada2.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart:

1. Raw paired-end reads (FASTQ).
2. Quality-trim + primer-remove.
3. Learn per-sample error rates.
4. Denoise → ASVs.
5. Merge paired ends.
6. Remove chimeras.
7. ASV table (samples × ASVs).
8. Taxonomy assignment (SILVA/GTDB).
9. Phylogenetic tree (FastTree).
10. Diversity + downstream analysis.

Side panel: schematic of error model — Illumina sequencer error PMF; how DADA2 distinguishes single-base sequencing-error from real biological variant.

---

## Figure 3 — Shotgun metagenomics pipeline

**File**: `diagrams/lecture-23/03-shotgun-pipeline.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom:

1. Raw shotgun reads (~5 Gb FASTQ).
2. QC + adapter trimming (fastp).
3. Host-DNA removal (KneadData / BBDuk).
4. Branch:
   a. Taxonomic classification: Kraken2 or MetaPhlAn → relative abundance.
   b. Functional profiling: HUMAnN → pathway abundance.
   c. Assembly: megahit/SPAdes → contigs → binning → MAGs.

Each stage annotated with typical runtime and tool choice. Side panel: MAG quality criteria (CheckM completeness > 90%, contamination < 5%).

---

## Figure 4 — Diversity metrics

**File**: `diagrams/lecture-23/04-diversity.svg`
**ViewBox**: `0 0 1200 720`

Two panels:

**Top — alpha diversity**:
- Two example sample composition pie charts: rich (10 species evenly distributed) vs poor (2 species, 80%+20%).
- Computed Shannon, Simpson, observed richness for each.
- Labels show: rich Shannon = 2.30, poor = 0.50.

**Bottom — beta diversity**:
- Distance matrix between 5 samples (heat map).
- PCoA scatter showing sample clustering.
- Two clusters visible: case (red) vs control (cobalt).
- Annotation: "PERMANOVA p < 0.001 between groups".

Side annotations: Shannon = $-\sum p_i \log p_i$; Bray-Curtis = compositional difference.

---

## Figure 5 — Microbiome-disease association case study

**File**: `diagrams/lecture-23/05-microbiome-disease.svg`
**ViewBox**: `0 0 1200 720`

Three-panel case study (immunotherapy response):

- Panel 1: PCoA of pre-treatment microbiomes; responders (red) and non-responders (grey) cluster differently.
- Panel 2: Differential abundance of *Akkermansia muciniphila* — significantly higher in responders (volcano plot).
- Panel 3: Kaplan-Meier survival curves split by *A. muciniphila* abundance — high vs low.

Bottom annotation: "Routy et al. 2018 — gut microbiome predicts anti-PD-1 response in NSCLC + RCC patients."

---

## Figure 6 — Microbiome analysis workflow

**File**: `diagrams/lecture-23/06-workflow.svg`
**ViewBox**: `0 0 1200 720`

A horizontal-banded flowchart showing the full analysis from sample to publication:

- Sample collection, storage, extraction.
- Sequencing.
- DADA2 / shotgun pipeline.
- Diversity + ordination.
- Differential abundance.
- Functional profiling.
- Visualisation + interpretation.

Each band annotated with reproducibility-critical steps (DB version pinning, primer specification, FDR control).

---

## Figure 7 — 16S rRNA gene structure

**File**: `diagrams/lecture-23/07-16s-structure.svg`
**ViewBox**: `0 0 1200 480`

Linear 16S gene drawn as horizontal bar:

- 9 alternating regions: V1 (variable, blue), C1 (conserved, grey), V2 (variable), ... up to V9.
- Total length ~1500 bp.
- Common primer pairs annotated: V3-V4 (~460 bp), V4 (~250 bp), full-length (1500 bp on PacBio).
- Annotation: "primers anneal to conserved regions, amplify variable region in between."

---

## Figure 8 — Compositional data problem

**File**: `diagrams/lecture-23/08-compositional.svg`
**ViewBox**: `0 0 1200 600`

Two scenarios shown side-by-side:

- Scenario 1 (true growth): species A doubles in absolute count; species B stays the same.
- Scenario 2 (technical artefact): both species' absolute counts unchanged, sequencing depth halved for sample 2.

Both produce the same relative-abundance pattern. Annotation: "relative abundance can change without absolute change. Compositional analysis (Aitchison transforms, ANCOM-BC) handles this; naive t-tests don't."

---

## Figure 9 — MAG recovery from contigs

**File**: `diagrams/lecture-23/09-mag-recovery.svg`
**ViewBox**: `0 0 1200 600`

Left panel: assembled metagenomic contigs (~1000 contigs, varied lengths) plotted as scatter:
- X-axis: tetranucleotide composition (PCA coordinate 1).
- Y-axis: read coverage depth.

Three clusters visible — three distinct genomes binned by composition + coverage.

Right panel: each cluster's CheckM quality:
- MAG 1: 95% complete, 1% contamination → high quality.
- MAG 2: 75%, 3% → medium.
- MAG 3: 50%, 8% → low (needs more data).

---

## Figure 10 — UniFrac phylogenetic distance

**File**: `diagrams/lecture-23/10-unifrac.svg`
**ViewBox**: `0 0 1200 540`

Phylogenetic tree of 8 ASVs from 2 samples (A, B):

- Branches coloured red if only sample A has descendants; cobalt if only B; black if both.
- Unweighted UniFrac = ratio of single-colour branch length to total branch length.
- Weighted UniFrac additionally factors in abundance per branch.

Side panel: comparison example showing same Bray-Curtis but different UniFrac (when shared taxa are phylogenetically close vs far).

---

## Figure 11 — Strain-level resolution

**File**: `diagrams/lecture-23/11-strain-resolution.svg`
**ViewBox**: `0 0 1200 540`

Two columns:

- Short-read 16S (V3-V4): assembled to ~genus level. Multiple species clustered together due to incomplete sequence info.
- Long-read 16S (PacBio HiFi, full 1500 bp): species-level resolution. Same sample, different inferred composition.

Annotation: "long reads resolve closely-related species that V3-V4 fragments cannot."

---

## Figure 12 — Antibiotic resistance gene tracking

**File**: `diagrams/lecture-23/12-arg-tracking.svg`
**ViewBox**: `0 0 1200 600`

A timeline of ARG abundance in a hospital sewage metagenome over 2 years:

- Y-axis: log abundance of ARG class.
- X-axis: time (months).
- 5 ARG classes plotted: tetracycline, beta-lactam, fluoroquinolone, vancomycin, carbapenem.
- Vertical line at month 8 marking institutional antibiotic policy change; visible inflection in some classes.

Side annotation: "shotgun metagenomics tracks horizontal gene transfer in real time across hospital ecosystems."
