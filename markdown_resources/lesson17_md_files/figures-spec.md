# Lecture 17 — Figures Specification

> **Scope**: Static diagrams for Lecture 17 (Clinical Genomics, Variant Interpretation, and Ethics).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-17.md`.

---

## 0. Conventions for This Lecture

- Figures are custom SVG; content is workflow / table / schematic heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Pathogenic colour: red `#c4342c`; likely pathogenic: orange `#b45309`; VUS: grey `#525252`; likely benign: muted cobalt `#6e88c8`; benign: cobalt `#1e3a8a`.
- Regulatory colours: FDA-cleared amber-dark; LDT amber-light; RUO grey; IVD cobalt.
- Ancestry palette consistent with Lectures 12 + 13: African cobalt, European amber, East Asian red, South Asian green, Admixed violet, Oceanian teal.
- Typography: Inter for UI labels; JetBrains Mono for HGVS variant names, gene names, star alleles, allele-frequency values.
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Research vs clinical pipeline comparison | Part 1 | Custom SVG |
| 2 | ACMG/AMP five-class distribution | Part 2 | Custom SVG |
| 3 | ACMG/AMP evidence code table | Part 2 | Custom SVG |
| 4 | ClinVar growth + classification distribution | Part 2 | Custom SVG |
| 5 | gnomAD allele frequency and PM2/BA1 thresholds | Part 3 | Custom SVG |
| 6 | REVEL / AlphaMissense / SpliceAI comparison | Part 3 | Custom SVG |
| 7 | Pharmacogenomics worked examples | Part 4 | Custom SVG |
| 8 | ACMG SF list v3.x categories | Part 5 | Custom SVG |
| 9 | FDA approval pathway | Part 6 | Custom SVG |
| 10 | GINA coverage map | Part 7 | Custom SVG |
| 11 | Clinical lab regulatory landscape | Part 1 | Custom SVG |
| 12 | Ancestry bias in ClinVar and predictors | Part 7 | Custom SVG |

---

## Figure 1 — Research vs clinical pipeline comparison

**File**: `diagrams/lecture-17/01-research-vs-clinical.svg`
**Lecture anchor**: §1.3 LDTs, IVDs, RUO
**ViewBox**: `0 0 1200 540`

### Purpose

Two parallel horizontal workflow tracks contrasting research and clinical pipelines step-by-step.

### Content

**Top band — Research pipeline**. Horizontal steps with arrows between:

- Sample → Sequence → Align → Variant call → Analysis → Manuscript.

Under each step small annotations: "1-off analyses", "tool version floats", "months TAT", "publication-grade".

**Bottom band — Clinical pipeline**. Mirror horizontal steps:

- Sample (with chain-of-custody icon) → Sequence (pinned reagent lot) → Align (pinned reference MD5) → Variant call (validated SOP) → Curation (named variant scientist) → Signed clinical report (physician signature icon) → EHR deposit.

Under each step, annotations: "SOP documented", "CLIA/CAP audited", "3–14 day TAT", "physician-signed".

**Central divider**. Vertical line with the caption: "Regulatory boundary — everything below is regulated software".

### Style notes

- Research track: `--fg-muted` boxes, thin borders.
- Clinical track: `--accent` cobalt boxes with amber borders (regulatory emphasis).
- Signature / chain-of-custody icons drawn inline with simple geometric primitives.

---

## Figure 2 — ACMG/AMP five-class distribution

**File**: `diagrams/lecture-17/02-acmg-classes.svg`
**Lecture anchor**: §2.2 The five-class output
**ViewBox**: `0 0 1200 520`

### Purpose

Show the ACMG/AMP five-class output, the probability bands behind each, and the typical population distribution of variants observed in clinical exomes.

### Content

**Top band — classification spectrum**. Horizontal bar divided into five segments labelled B, LB, VUS, LP, P. Colour per segment using the L17 conventions. Below each segment: probability range ("≤ 0.1% pathogenic", "0.1–10%", "10–90%", "90–99%", "≥ 99%") and an action label ("not reported", "not reported", "reported, no action", "reported, actionable", "reported, actionable").

**Bottom band — real-world frequency distribution**. A stacked histogram or a shaped distribution curve showing the proportion of variants in a typical clinical exome that fall into each class: ~85% B/LB (short bars), ~10% VUS (large bar), ~5% P/LP (small bar).

**Annotation band below the histogram**. "VUS reclassification is a continuous field-wide effort — today's VUS is tomorrow's P or LB."

### Style notes

- Five-segment spectrum with the canonical class-colour palette.
- Histogram bars: grey fill with colour-tinted outlines per class.

---

## Figure 3 — ACMG/AMP evidence code table

**File**: `diagrams/lecture-17/03-evidence-codes.svg`
**Lecture anchor**: §2.3 Evidence codes and their strengths
**ViewBox**: `0 0 1200 640`

### Purpose

Provide a complete, visually scannable reference for the 28 ACMG/AMP evidence codes with their strengths and brief descriptions.

### Content

**Two-column table** (pathogenic on the left, benign on the right).

**Left column — pathogenic evidence** (top to bottom, strongest to weakest):

- **PVS1** (Very Strong): LoF in gene where LoF is a known mechanism.
- **PS1** (Strong): same amino-acid change as known pathogenic.
- **PS2** (Strong): de novo with confirmed parentage.
- **PS3** (Strong): well-established functional studies show damaging.
- **PS4** (Strong): significantly enriched in cases vs controls.
- **PM1** (Moderate): hotspot / critical domain.
- **PM2** (Moderate): absent or ultra-rare in gnomAD.
- **PM3** (Moderate): in trans with pathogenic in recessive.
- **PM4** (Moderate): protein-length change.
- **PM5** (Moderate): novel missense where diff AA at same residue is pathogenic.
- **PM6** (Moderate): assumed de novo.
- **PP1** (Supporting): co-segregation.
- **PP2** (Supporting): missense in low-benign-missense gene.
- **PP3** (Supporting): multiple predictors agree damaging.
- **PP4** (Supporting): highly-specific phenotype.
- **PP5** (Supporting): reported pathogenic (deprecated 2018).

**Right column — benign evidence**:

- **BA1** (Stand-alone): ≥ 5% in gnomAD.
- **BS1** (Strong): frequency inconsistent with disease.
- **BS2** (Strong): observed in healthy adults.
- **BS3** (Strong): functional studies show no damage.
- **BS4** (Strong): lack of segregation.
- **BP1** (Supporting): missense in gene where LoF is mechanism.
- **BP2** (Supporting): in trans with P in healthy individual.
- **BP3** (Supporting): in-frame indel in repetitive region.
- **BP4** (Supporting): multiple predictors agree benign.
- **BP5** (Supporting): alternative diagnosis present.
- **BP6** (Supporting): reported benign (deprecated 2018).
- **BP7** (Supporting): silent, no predicted splice effect.

### Style notes

- Left column: red gradient from deep (PVS1) to light (PP).
- Right column: cobalt gradient from deep (BA1) to light (BP).
- Deprecated codes: strikethrough with "(deprecated)" tag.
- Code names in JetBrains Mono.

---

## Figure 4 — ClinVar growth + classification distribution

**File**: `diagrams/lecture-17/04-clinvar.svg`
**Lecture anchor**: §2.7 ClinVar as community truth
**ViewBox**: `0 0 1200 560`

### Purpose

Show both the growth of ClinVar over time and the current distribution of classifications across unique variants.

### Content

**Top half — growth curve**. Line plot. X-axis: year from 2013 (launch) to 2024. Y-axis (log): ClinVar submissions, 10⁴ to 10⁷. Curve rising from ~50k at launch to ~5M by 2024. Annotations at major milestones: "expert-panel submissions begin", "ClinGen harmonisation", "large-submitter policy tightened".

**Bottom half — classification distribution**. A donut chart with five slices:

- P + LP: ~15% (red + orange).
- VUS: ~40% (grey).
- LB + B: ~40% (cobalt + muted cobalt).
- Conflicting: ~5% (amber).

Labels on each slice with percentage.

### Style notes

- Growth curve: cobalt line with annotated milestones as diamond markers.
- Donut chart: standard per-class colour conventions.

---

## Figure 5 — gnomAD allele frequency and PM2/BA1 thresholds

**File**: `diagrams/lecture-17/05-gnomad-thresholds.svg`
**Lecture anchor**: §3.1 Population frequency: gnomAD
**ViewBox**: `0 0 1200 480`

### Purpose

Visualise the key allele-frequency thresholds used in ACMG/AMP classification, mapped to the actual distribution of variants.

### Content

**Horizontal frequency axis**, log-scale 10⁻⁷ to 10⁰ (singleton to 100%). Above the axis:

- **Shaded cobalt region**, < 1/200k: "PM2 applicable for dominant disease".
- **Shaded cobalt region**, < 1/20k: "PM2 applicable for recessive disease".
- **Shaded grey region**, 1/20k to 1/100: "neither PM2 nor BA1".
- **Shaded amber region**, ≥ 5%: "BA1 stand-alone benign".

Between the shaded regions, thin vertical threshold lines labelled with the exact frequencies.

**Below the axis** — two example variants plotted:

- A rare missense at 1/500k (ultra-rare; "PM2_Supporting applicable").
- A common SNP at 12% (above BA1; "BA1 stand-alone").

### Style notes

- Log-scale axis with major ticks at each decade.
- Shaded regions with semi-transparent fill.
- Example variants as labelled dots with arrows.

---

## Figure 6 — REVEL / AlphaMissense / SpliceAI comparison

**File**: `diagrams/lecture-17/06-predictors.svg`
**Lecture anchor**: §3.2 Functional impact; §3.3 SpliceAI
**ViewBox**: `0 0 1200 560`

### Purpose

Compare the three major in-silico predictors by input type, output range, and typical threshold for ACMG/AMP use.

### Content

**Three side-by-side panels**, one per predictor.

**Panel 1 — REVEL**:

- Input: missense variants.
- Architecture sketch: ensemble of 13 simpler predictors (SIFT, PolyPhen-2, MutationTaster, etc.) fed into a random-forest combiner.
- Score range: 0–1.
- Thresholds: ≥ 0.7 → PP3; ≤ 0.15 → BP4.

**Panel 2 — AlphaMissense**:

- Input: missense variants.
- Architecture sketch: AlphaFold2-derived; uses structural context; per-variant classifier head.
- Score range: 0–1 (with "likely pathogenic / ambiguous / likely benign" classes).
- Thresholds: ≥ 0.564 → PP3; ≤ 0.340 → BP4.

**Panel 3 — SpliceAI**:

- Input: any variant ± 10 kb of splice context.
- Architecture sketch: deep CNN on one-hot sequence.
- Score range: 0–1 per (acceptor gain / loss / donor gain / loss).
- Thresholds: ≥ 0.5 → splice effect supporting evidence; ≥ 0.8 → strong splice disruption.

**Bottom caption band**. "Multiple predictors must agree for PP3 or BP4; no single score is actionable alone."

### Style notes

- Panels visually similar but colour-coded per predictor.
- Architecture sketches simplified (box-and-arrow).
- Threshold values in monospace.

---

## Figure 7 — Pharmacogenomics worked examples

**File**: `diagrams/lecture-17/07-pgx-examples.svg`
**Lecture anchor**: §4.2 The canonical cases
**ViewBox**: `0 0 1200 640`

### Purpose

Four-panel grid covering the four canonical PGx examples, each with the full "drug / gene / variants / phenotype / action" chain.

### Content

**Panel 1 — CYP2D6 × codeine**:

- Drug: codeine.
- Gene: CYP2D6.
- Representative variants: *4 (LoF), *10 (reduced function), *1×N (gene duplication).
- Phenotypes: poor / intermediate / normal / rapid / ultra-rapid.
- Action: alternative opioid for PM or UM; standard for normal/intermediate.

**Panel 2 — Warfarin × VKORC1 + CYP2C9**:

- Drug: warfarin.
- Genes: VKORC1 (-1639G>A), CYP2C9 (*2, *3).
- Phenotype: warfarin sensitivity score.
- Action: genotype-guided dose nomogram (show a small example dose table).

**Panel 3 — Abacavir × HLA-B*57:01**:

- Drug: abacavir.
- Gene: HLA-B.
- Variant: HLA-B*57:01.
- Phenotype: HLA-B*57:01 positive / negative.
- Action: if positive, avoid abacavir (swap to alternative HIV regimen).

**Panel 4 — Thiopurines × TPMT / NUDT15**:

- Drugs: azathioprine, 6-mercaptopurine.
- Genes: TPMT (*2, *3A, *3B, *3C), NUDT15 (common in East Asians).
- Phenotypes: normal / intermediate / poor metaboliser.
- Action: dose reduce 10× for poor metabolisers.

### Style notes

- 2×2 grid of panels.
- Each panel with colour-coded drug at top, variant table in middle, action label at bottom.
- Star-allele labels in JetBrains Mono.

---

## Figure 8 — ACMG SF list v3.x categories

**File**: `diagrams/lecture-17/08-acmg-sf.svg`
**Lecture anchor**: §5.2 The ACMG SF list
**ViewBox**: `0 0 1200 560`

### Purpose

Visualise the composition of the ACMG Secondary Findings list by disease category.

### Content

**Left half — donut chart** of gene-count per category:

- Hereditary cancer syndromes (30 genes) — dominant slice.
- Cardiomyopathies + channelopathies (25 genes).
- Familial hypercholesterolemia (3 genes).
- Malignant hyperthermia (2 genes).
- Metabolic + other (21 genes).

**Right half — table** listing the top 10 SF genes with their associated diseases:

- BRCA1 / BRCA2 — hereditary breast / ovarian cancer.
- MLH1 / MSH2 / MSH6 / PMS2 — Lynch syndrome.
- TP53 — Li-Fraumeni syndrome.
- APC — familial adenomatous polyposis.
- LDLR — familial hypercholesterolemia.
- MYH7 / MYBPC3 — hypertrophic cardiomyopathy.
- RYR1 — malignant hyperthermia.
- KCNQ1 — long QT syndrome.
- NF1 — neurofibromatosis type 1.
- VHL — von Hippel-Lindau syndrome.

### Style notes

- Donut chart with segment colours per disease category.
- Gene names in monospace; disease names in Inter italic.

---

## Figure 9 — FDA approval pathway

**File**: `diagrams/lecture-17/09-fda-pathway.svg`
**Lecture anchor**: §6.1 FDA-cleared sequencing assays; §6.2 The 2023 FDA LDT rule
**ViewBox**: `0 0 1200 560`

### Purpose

Flowchart of possible regulatory paths for a genomic assay from conception to market.

### Content

**Starting node at the left**: "Assay concept (sequencing-based test)".

**Three branches**:

1. Upward branch: "RUO — Research Use Only". Annotated "non-clinical; label explicit". Ends with a "Research-only applications" terminal.
2. Middle branch: "LDT — Laboratory-Developed Test". Pre-2024: enforcement discretion (no FDA review). Post-2024: phased FDA oversight (4 stages shown as sub-boxes: adverse-event reporting → quality-system compliance → registration → premarket review). Terminal: "Clinical LDT in regulated lab".
3. Downward branch: "FDA-reviewed IVD". Splits into "510(k) — substantially equivalent" (faster) and "PMA — premarket approval" (rigorous; used for CDx approvals). Terminal: "Commercial IVD kit".

**Annotations** at each branch:

- RUO: "e.g. research reagents, academic tools".
- LDT: "e.g. most clinical exome labs pre-2024; transitioning 2024–2028".
- 510(k): "e.g. Oncomine Dx Target".
- PMA: "e.g. FoundationOne CDx, MSK-IMPACT authorisation".

### Style notes

- Top branch: `--fg-muted` grey (RUO).
- Middle branch: amber (LDT, with darker amber for post-2024 boxes).
- Bottom branch: cobalt (FDA-reviewed), with darker cobalt for PMA.

---

## Figure 10 — GINA coverage map

**File**: `diagrams/lecture-17/10-gina-coverage.svg`
**Lecture anchor**: §7.1 GINA and its limits
**ViewBox**: `0 0 1200 540`

### Purpose

Visualise the scope of GINA's protections and its significant carve-outs.

### Content

**Central Venn-like region** labelled "Covered by GINA":

- Health insurance (group + individual).
- Employment (large employers, ≥ 15 employees).

**Outer regions** labelled "NOT Covered by GINA":

- Life insurance.
- Disability insurance.
- Long-term care insurance.
- Small employers (< 15 employees).
- US military.

**State-level additions** shown as small inserts / annotations:

- California (state-level genetic nondiscrimination in insurance).
- Florida.
- Vermont.
- A handful of other states.

**Bottom caption band**. "Carve-outs reflect 2008 legislative compromises; bills to close them have been introduced periodically but not passed."

### Style notes

- Covered region: cobalt-filled Venn circle.
- Not-covered regions: red outlines.
- State-level additions: small green highlights over the red regions where state law fills the gap.

---

## Figure 11 — Clinical lab regulatory landscape

**File**: `diagrams/lecture-17/11-regulatory-landscape.svg`
**Lecture anchor**: §1.2 CLIA, CAP, and the lab ecosystem; §1.3 LDTs, IVDs, RUO
**ViewBox**: `0 0 1200 560`

### Purpose

Show the layered regulatory environment of a clinical lab — physical infrastructure, lab certification, test classification — and how they stack.

### Content

**Bottom layer — physical infrastructure**. Wide band. Items shown as labelled icons:

- Sequencer (FDA-cleared hardware, e.g. Illumina NovaSeq).
- Reagent kit (IVD-certified).
- Library-prep tubing.

**Middle layer — lab certification**. Three overlapping rings:

- CLIA (required, baseline).
- CAP accreditation (voluntary, expected).
- NYS CLEP (strictest).

Overlap region in the centre labelled "fully-certified clinical genomics lab".

**Top layer — test output classification**. Three tiles:

- RUO (research use only; greyed out).
- LDT (lab-developed test; amber; post-2024 regulated by FDA).
- IVD (FDA-reviewed commercial kit; cobalt).

**Arrows** from bottom → middle → top indicating the regulatory flow: physical tools + certified lab + classified output = defensible clinical result.

### Style notes

- Layered bands visually distinct.
- Ring overlaps with tint blending.
- 2024 FDA LDT change highlighted with a small amber-highlight marker on the LDT tile.

---

## Figure 12 — Ancestry bias in ClinVar and predictors

**File**: `diagrams/lecture-17/12-ancestry-bias.svg`
**Lecture anchor**: §7.2 Ancestry bias in clinical databases
**ViewBox**: `0 0 1200 640`

### Purpose

Quantify ancestry bias across three layers: ClinVar submissions, gnomAD composition, and predictor accuracy.

### Content

**Three stacked horizontal stacked-bar charts**, each spanning the same width.

**Top chart — ClinVar submissions by cohort ancestry**:

- European: ~70%.
- East Asian: ~10%.
- African: ~6%.
- South Asian: ~3%.
- Admixed American: ~2%.
- Unknown / mixed: ~9%.

Title: "ClinVar submissions (as of 2024)".

**Middle chart — gnomAD v4 composition**:

- European (non-Finnish): 36%.
- Finnish: 5%.
- East Asian: 10%.
- African: 20%.
- Admixed American: 13%.
- South Asian: 6%.
- Ashkenazi Jewish: 3%.
- Other: 7%.

Title: "gnomAD v4 ancestry composition (~800k exomes)".

**Bottom chart — predictor accuracy by ancestry**:

- European F1: ~0.88.
- Admixed American: ~0.83.
- East Asian: ~0.84.
- South Asian: ~0.82.
- African: ~0.78.

Title: "REVEL / AlphaMissense validation F1 by ancestry".

**Annotation band at bottom**. "Accuracy follows training-distribution coverage. Diverse data, not algorithm tweaks, closes the gap."

### Style notes

- Ancestry palette consistent with Lectures 12 and 13.
- Stacked bars in top two charts; grouped bars in bottom chart.
- F1 axis: 0.7–0.95 range.
