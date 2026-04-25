# Lecture 24 (proposed L21) — Figures Specification

> **Scope**: Static diagrams for Lecture 24 (CRISPR Screens and DepMap).
> **Companion files**: `lecture-24.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Screen-direction colours: positive selection red, dropout cobalt.
- Library colours: Brunello amber, TKOv3 green, custom violet.
- Cell-line lineages: hematopoietic red, breast cobalt, lung amber, colon green, glioma violet.
- Score-distribution colours: essential dropout cobalt, non-essential null grey, resistance positive red.
- Typography: Inter for chrome; JetBrains Mono for sgRNA sequences, gene symbols, scores.

## Figure budget — 12 figures

| # | Title | Part |
|---|---|---|
| 1 | Pooled CRISPR screen workflow | Part 1 |
| 2 | sgRNA library quality characteristics | Part 2 |
| 3 | MAGeCK analysis flow | Part 3 |
| 4 | DepMap dependency landscape | Part 4 |
| 5 | MAVE variant effect map | Part 5 |
| 6 | CRISPR screen analysis workflow | Part 6 |
| 7 | CRISPRko vs CRISPRi vs CRISPRa comparison | Part 1 |
| 8 | sgRNA on-target / off-target scoring | Part 2 |
| 9 | RRA rank aggregation in MAGeCK | Part 3 |
| 10 | Synthetic lethality discovery via DepMap | Part 4 |
| 11 | Pooled screen as compressed sensing | Part 2 |
| 12 | From screen to clinical drug | Part 7 |

---

## Figure 1 — Pooled CRISPR screen workflow

**File**: `diagrams/lecture-24/01-pooled-screen.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flow:

1. sgRNA library cloned into lentiviral plasmid pool.
2. Lentivirus production + transduction at MOI ~ 0.3.
3. Cells split into T0 (initial) and T-final (post-selection or growth).
4. Genomic DNA extraction from each timepoint.
5. PCR-amplify sgRNA region.
6. Sequence amplicons.
7. Count table (sgRNAs × samples).
8. MAGeCK / DrugZ analysis.
9. Gene-level fitness scores.

Side annotation: typical scale (~80k sgRNAs, 4 sgRNAs/gene, ~10⁸ cells per replicate, $20k cost).

---

## Figure 2 — sgRNA library quality characteristics

**File**: `diagrams/lecture-24/02-library-quality.svg`
**ViewBox**: `0 0 1200 600`

Two panels:

- Top: distribution of on-target activity scores across libraries (Brunello, TKOv3, GeCKO). Modern libraries skewed right (high activity); older libraries broader.
- Bottom: distribution of predicted off-target counts. Modern libraries skewed left (few off-targets).

Side panel: number of sgRNAs per gene by library.

---

## Figure 3 — MAGeCK analysis flow

**File**: `diagrams/lecture-24/03-mageck.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom:

1. Count tables (T0, T-final).
2. Median normalisation.
3. Per-sgRNA negative-binomial test → sgRNA p-values.
4. Per-gene RRA aggregation (4 sgRNAs → median rank).
5. Beta(1, k) significance test → gene p-value.
6. FDR correction (BH).
7. Gene rankings + significance.

Side panel: example output table — top 10 essential genes for a cancer cell line.

---

## Figure 4 — DepMap dependency landscape

**File**: `diagrams/lecture-24/04-depmap.svg`
**ViewBox**: `0 0 1200 720`

Heat map: ~30 cell lines (rows) × top ~50 genes (columns). Cells coloured by Chronos score.

- Pan-essential genes (essential everywhere): solid blue column.
- Pan-non-essential: solid white column.
- Lineage-selective: blue in some lineages, white in others (bands of selectivity).

Annotation: "lineage-selective dependencies = drug targets. Pan-essential = avoid (would kill normal cells too)."

---

## Figure 5 — MAVE variant effect map

**File**: `diagrams/lecture-24/05-mave-map.svg`
**ViewBox**: `0 0 1200 720`

A heat map: amino acid positions (x-axis, ~400 positions) × possible amino acid changes (y-axis, 20 AAs).

- Cells coloured by experimental fitness score (red = damaging, blue = neutral, green = beneficial).
- Wild-type residues shown as black dots.
- Functional sites (catalytic residues, binding sites) annotated as bands of red columns.

Side panel: ClinVar Pathogenic variants overlaid; concordance score with MAVE.

---

## Figure 6 — CRISPR screen analysis workflow

**File**: `diagrams/lecture-24/06-workflow.svg`
**ViewBox**: `0 0 1200 720`

Horizontal banded flowchart:

- Sample acquisition + library transduction.
- T0 + T-final sequencing.
- QC checks (MOI, recovery, replicate concordance, CEG dropout).
- MAGeCK / DrugZ / BAGEL2 (run all three).
- CN-bias correction (CRISPRcleanR).
- Gene-level rankings.
- Pathway enrichment.
- Hit triage (druggability, biological context).
- Validation arrayed screen.

Each stage with reproducibility-critical step annotations.

---

## Figure 7 — CRISPRko vs CRISPRi vs CRISPRa

**File**: `diagrams/lecture-24/07-modalities.svg`
**ViewBox**: `0 0 1200 600`

Three side-by-side panels:

- CRISPRko (Cas9 + sgRNA): cuts DNA → indels → frameshift → loss of function. Permanent, fatal for essential genes.
- CRISPRi (dCas9-KRAB + sgRNA at TSS): represses transcription. Reversible, useful for essential genes.
- CRISPRa (dCas9-VP64 + sgRNA at TSS): activates transcription. Gain-of-function screens.

Each panel: schematic, typical use cases, library size, advantages / drawbacks.

---

## Figure 8 — sgRNA on-target / off-target scoring

**File**: `diagrams/lecture-24/08-sgrna-scoring.svg`
**ViewBox**: `0 0 1200 540`

Two panels:

- On-target: sgRNA features (PAM context, GC content, position-specific nucleotides) → Doench Rule Set 2 score. Histogram of typical-library scores.
- Off-target: PAM-distal mismatches → CFD score. Heat map of mismatch tolerance per position.

Side annotation: "high-quality sgRNA = high on-target + low off-target. Library design balances both."

---

## Figure 9 — RRA rank aggregation in MAGeCK

**File**: `diagrams/lecture-24/09-rra.svg`
**ViewBox**: `0 0 1200 540`

A scatter plot: 4 sgRNAs per gene, x-axis = rank-percentile under null.

- Highly-essential gene: all 4 sgRNAs in top 1% (concordant strong signal).
- Modestly-essential: 2-3 sgRNAs in top 10%, 1-2 in middle.
- Non-essential: sgRNAs uniformly distributed.

RRA computes the minimum percentile under Beta(1, k) null; outputs gene p-value.

Annotation: "MAGeCK's modified RRA is robust to a few non-functional sgRNAs — uses median rank, not mean."

---

## Figure 10 — Synthetic lethality via DepMap

**File**: `diagrams/lecture-24/10-synthetic-lethality.svg`
**ViewBox**: `0 0 1200 600`

Two-panel analysis:

- Panel 1: scatter — cell lines with MTAP deletion (red, n=200) vs MTAP-intact (cobalt, n=800).
- Panel 2: PRMT5 essentiality (Chronos score) — significantly stronger in MTAP-deleted cells.

→ MTAP-deletion + PRMT5 dependency = synthetic-lethal interaction → PRMT5 drug development for MTAP-deleted cancer.

Annotation: "PRMT5 inhibitors (TANGO / GSK) are in clinical trials based on this DepMap finding."

---

## Figure 11 — Pooled screen as compressed sensing

**File**: `diagrams/lecture-24/11-compressed-sensing.svg`
**ViewBox**: `0 0 1200 600`

Mathematical schematic:

- Population with sgRNA abundance vector $\mathbf{x}$ (size $N$).
- Phenotype readout = linear measurement $y = A x + \epsilon$ where $A$ is the sgRNA → gene mapping matrix and $\epsilon$ is sequencing noise.
- Inverse problem: recover gene-level effect from sgRNA-level abundance changes.

Side annotation: "MAGeCK / DrugZ / BAGEL2 are statistical solvers for this inverse problem. Same family as compressed sensing in EE."

---

## Figure 12 — From screen to clinical drug

**File**: `diagrams/lecture-24/12-screen-to-drug.svg`
**ViewBox**: `0 0 1200 720`

Pipeline diagram:

1. CRISPR screen identifies dependency.
2. Hit validation (arrayed knockout, isogenic).
3. Druggability assessment (small-molecule pocket, surface accessible).
4. Hit-to-lead chemistry.
5. Pre-clinical validation (xenograft).
6. Phase I clinical trial (safety).
7. Phase II clinical trial (efficacy).
8. Phase III + FDA approval.

Time annotation: ~7-10 years from screen to FDA. Side panel: example real timelines (PARP inhibitors, PRMT5 inhibitors).
