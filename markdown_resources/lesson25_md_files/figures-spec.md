# Lecture 25 (proposed L20) — Figures Specification

> **Scope**: Static diagrams for Lecture 25 (Causal Inference and Mendelian Randomisation).
> **Companion files**: `lecture-25.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Causal-DAG colours: exposure cobalt, outcome red, instrument amber, confounder grey.
- IV assumptions colours: relevance green, independence cobalt, exclusion red (the strict).
- Pleiotropy colours: vertical (OK) green, horizontal (NOT OK) red.
- Typography: Inter for chrome; JetBrains Mono for SNP IDs, $\beta$ values, p-values.

## Figure budget — 12 figures

| # | Title | Part |
|---|---|---|
| 1 | Causation vs correlation: three explanations | Part 1 |
| 2 | 2SLS estimator illustrated | Part 2 |
| 3 | MR analysis flow | Part 3 |
| 4 | Mediation analysis decomposition | Part 4 |
| 5 | Drug-target MR pipeline | Part 6 |
| 6 | Causal inference toolkit map | Part 7 |
| 7 | Wald ratio: per-SNP estimator | Part 2 |
| 8 | MR-Egger plot: detecting pleiotropy | Part 3 |
| 9 | MR-PRESSO outlier detection | Part 3 |
| 10 | Steiger test: directionality | Part 3 |
| 11 | Vertical vs horizontal pleiotropy | Part 2 |
| 12 | Cis-MR for drug target validation | Part 6 |

---

## Figure 1 — Causation vs correlation: three explanations

**File**: `diagrams/lecture-25/01-three-explanations.svg`
**ViewBox**: `0 0 1200 540`

Three side-by-side causal-DAG panels:

- Panel 1 (causation): X → Y. Same correlation as observed.
- Panel 2 (reverse causation): Y → X. Same correlation but opposite direction.
- Panel 3 (confounding): C → X, C → Y. Same correlation but neither causes the other.

Bottom annotation: "the same observed association is consistent with all three. Distinguishing them requires extra information — instrumental variables, time-ordering, or experimental manipulation."

---

## Figure 2 — 2SLS estimator illustrated

**File**: `diagrams/lecture-25/02-2sls.svg`
**ViewBox**: `0 0 1200 720`

Two-panel mechanical illustration:

- Panel 1: stage-1 regression — scatter of $X$ vs $Z$, OLS line.
- Panel 2: stage-2 regression — scatter of $Y$ vs $\widehat{X}$ (the fitted X from stage 1), OLS line.

Side annotation:

- $\widehat{X}$ is the part of $X$ that's instrument-driven.
- This part is uncorrelated with confounders (under IV assumptions).
- Stage-2 slope = causal effect.

Bottom: "2SLS = projection (stage 1) + OLS (stage 2)."

---

## Figure 3 — MR analysis flow

**File**: `diagrams/lecture-25/03-mr-flow.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart:

1. Identify exposure GWAS hits (~100 SNPs).
2. Look up effect sizes in outcome GWAS.
3. Compute per-SNP Wald ratios.
4. IVW meta-analysis → primary estimate.
5. Sensitivity analyses: MR-Egger, weighted median, MR-PRESSO.
6. Convergence check across methods.
7. Steiger directionality.
8. Final causal effect with calibrated CI.

Side panel: example output for LDL → CAD.

---

## Figure 4 — Mediation analysis decomposition

**File**: `diagrams/lecture-25/04-mediation.svg`
**ViewBox**: `0 0 1200 600`

Causal DAG:

- $X$ → $Y$ (direct effect, NDE).
- $X$ → $M$ → $Y$ (indirect effect, NIE).

Numerical example:

- Total effect: 1.0.
- Direct effect (NDE): 0.7.
- Indirect via M (NIE): 0.3.

Bottom annotation: "for BMI → CAD via LDL: ~30% of the effect is mediated by LDL; ~70% is direct (or via other pathways)."

---

## Figure 5 — Drug-target MR pipeline

**File**: `diagrams/lecture-25/05-drug-target-mr.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom:

1. Identify drug target gene (e.g., PCSK9).
2. Find cis-eQTLs in tissue (GTEx).
3. Use cis-eQTLs as IVs for target's expression.
4. MR for "target inhibition → disease" (LDL, CAD, etc.).
5. PheWAS-style scan: target inhibition vs many phenotypes.
6. Causal effect on disease + side-effect risks.
7. Validation: experimental + clinical trials.

Side panel: PCSK9 example — predicted CAD reduction confirmed by alirocumab/evolocumab trials.

---

## Figure 6 — Causal inference toolkit map

**File**: `diagrams/lecture-25/06-toolkit.svg`
**ViewBox**: `0 0 1200 720`

A 2D map:

- Y-axis: data setting (observational, quasi-experimental, randomised).
- X-axis: method (descriptive correlation → mediation → MR → DiD → RD → RCT).

Each method placed by its quasi-experimental requirements; MR sits in the upper-quasi-experimental quadrant.

Side annotation: "Methods are not interchangeable; pick by data + assumptions."

---

## Figure 7 — Wald ratio: per-SNP estimator

**File**: `diagrams/lecture-25/07-wald.svg`
**ViewBox**: `0 0 1200 540`

For one SNP: scatter showing $X$ vs $Z$ (slope = $\beta_{ZX}$) and $Y$ vs $Z$ (slope = $\beta_{ZY}$).

Wald ratio = $\beta_{ZY} / \beta_{ZX}$. Visually drawn as the rise-over-rise of the two scatter relationships.

Side annotation: "if SNP moves X by 0.5 and Y by 0.3, the per-unit-X causal effect is 0.6."

---

## Figure 8 — MR-Egger plot: detecting pleiotropy

**File**: `diagrams/lecture-25/08-mr-egger.svg`
**ViewBox**: `0 0 1200 600`

Scatter plot:

- X-axis: $\beta_{ZX}$ (per-SNP exposure effect).
- Y-axis: $\beta_{ZY}$ (per-SNP outcome effect).
- 100 SNP points.
- IVW regression line passes through origin (slope = causal effect).
- MR-Egger regression has free intercept; intercept ≠ 0 indicates systematic pleiotropy.

Two scenarios:

- Left: Egger intercept ≈ 0 → no pleiotropy → IVW is unbiased.
- Right: Egger intercept ≠ 0 → systematic pleiotropy → trust Egger slope, not IVW.

---

## Figure 9 — MR-PRESSO outlier detection

**File**: `diagrams/lecture-25/09-mr-presso.svg`
**ViewBox**: `0 0 1200 540`

Same MR scatter as Figure 8, but:

- Most SNPs cluster around a consistent slope.
- 5 SNPs are visible outliers.
- MR-PRESSO test identifies these; computes "global heterogeneity" significance.
- After outlier removal, IVW estimate becomes consistent with weighted median.

Annotation: "outliers often correspond to known pleiotropic loci (HLA, ABO, lipid hubs); flag and exclude."

---

## Figure 10 — Steiger directionality test

**File**: `diagrams/lecture-25/10-steiger.svg`
**ViewBox**: `0 0 1200 540`

For a candidate IV: bar chart showing per-SNP "fraction of variance explained" in $X$ (cobalt) vs $Y$ (red).

- Most SNPs explain more in $X$ than $Y$ → consistent with $X \to Y$ direction.
- A few SNPs explain more in $Y$ than $X$ → suggests reverse causation or pleiotropy.

Annotation: "Steiger test = significance of the rank-based comparison. Robust direction inference."

---

## Figure 11 — Vertical vs horizontal pleiotropy

**File**: `diagrams/lecture-25/11-pleiotropy.svg`
**ViewBox**: `0 0 1200 540`

Two side-by-side DAGs:

- Vertical (OK): $Z$ → $X$ → $Y$. Wald ratio is unbiased.
- Horizontal (NOT OK): $Z$ → $X$ → $Y$ AND $Z$ → $Y$ via independent path. Exclusion restriction violated.

Side panel: example of each — vertical: SNP near LDL receptor → LDL → CAD; horizontal: SNP affecting both lipid metabolism and inflammation independently.

---

## Figure 12 — Cis-MR for drug target validation

**File**: `diagrams/lecture-25/12-cis-mr.svg`
**ViewBox**: `0 0 1200 720`

A genome track view of a target gene (e.g., PCSK9, 25 kb region):

- TSS marked.
- Gene body coloured cobalt.
- ±100 kb window highlighted.
- cis-eQTL SNPs marked as red ticks above the track; trans-eQTLs (genome-wide) in grey.

Side panel:

- Cis-MR uses only the red ticks → narrow IV set.
- Strict cis-IVs reduce horizontal pleiotropy risk → more interpretable causal effect.

Annotation: "drug-target MR uses cis-IVs for cleaner inference; the closer to the gene, the cleaner."
