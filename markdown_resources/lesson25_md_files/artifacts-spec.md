# Lecture 25 (proposed L20) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 25 (Causal Inference and Mendelian Randomisation).
> **Companion files**: `lecture-25.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-25/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- DAG colours: exposure cobalt, outcome red, instrument amber, confounder grey.
- IV-assumption colours: relevance green, independence cobalt, exclusion red.
- Pleiotropy colours: vertical green, horizontal red.
- Typography: Inter for chrome; JetBrains Mono for SNP IDs, $\beta$, p-values.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget — 7 interactive tools

| # | Title | Anchor |
|---|---|---|
| 1 | Causal vs correlational scenario explorer | §1 |
| 2 | 2SLS regression simulator | §2.3 |
| 3 | Two-sample MR with sensitivity panel | §3 |
| 4 | Mediation analysis decomposer | §4 |
| 5 | Pleiotropy detector | §3.3 |
| 6 | Drug-target MR explorer | §6 |
| 7 | IV-assumption diagnostic | §2.4 |

---

## Artifact #1 — Causal vs Correlational Scenario Explorer

**File**: `artifacts/lecture-25/01-causal-scenario.html`
**Anchor**: §1

### Teaching purpose

Generate synthetic data under each of three scenarios (causation, reverse causation, confounding); show that the correlation is identical, but the causal answer differs.

### UI layout

- Three tabs: causation | reverse | confounding.
- Each tab generates 1000 samples with controllable effect sizes.
- Output: scatter plot, Pearson correlation, "true causal effect" (known by simulation).
- Outcome banner: "all three scenarios produce r = 0.5; but only causation has a non-zero true X→Y effect."

### Target aha

The same observational pattern is consistent with three different causal stories.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — 2SLS Regression Simulator

**File**: `artifacts/lecture-25/02-2sls-simulator.html`
**Anchor**: §2.3

### Teaching purpose

Walk through stage-1 and stage-2 of 2SLS on synthetic data with a known true causal effect.

### UI layout

- Slider: confounding strength $\rho$ (0 to 1).
- Slider: instrument strength (F-statistic 5 to 100).
- Pre-loaded data simulating $X$, $Y$, $Z$, $C$ with confounding.
- Side-by-side panels: stage-1 fit, stage-2 fit.
- Computed estimate vs true effect.
- Outcome banner: "with strong instrument and modest confounding: 2SLS recovers true effect 1.0 (vs OLS biased to 1.4 by confounding)."

### Target aha

2SLS removes confounding bias; weak instruments inflate variance.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Two-Sample MR with Sensitivity Panel

**File**: `artifacts/lecture-25/03-mr-analysis.html`
**Anchor**: §3

### Teaching purpose

Run a complete two-sample MR analysis: IVW, MR-Egger, weighted median, MR-PRESSO. Show convergence vs divergence.

### UI layout

- Pre-loaded GWAS summary stats for ~50 SNPs (LDL → CAD example).
- Slider: pleiotropy injection level (perturb summary stats with horizontal pleiotropy).
- Output dashboard:
  - IVW estimate + 95% CI.
  - MR-Egger slope + intercept p-value.
  - Weighted median estimate.
  - MR-PRESSO global p-value + outlier list.
- Forest plot of all four methods.
- Outcome banner: "no pleiotropy: all converge on $\beta = 0.55$. With added pleiotropy: IVW biased but MR-Egger and weighted median recover unbiased estimate."

### Target aha

MR's robustness depends on pleiotropy level; sensitivity analyses reveal validity.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Mediation Analysis Decomposer

**File**: `artifacts/lecture-25/04-mediation.html`
**Anchor**: §4

### Teaching purpose

Compute total / direct / indirect effects via Baron-Kenny and counterfactual mediation.

### UI layout

- Sliders for path coefficients: $X \to M$, $M \to Y$, $X \to Y$ (direct).
- Computed total effect (NDE + NIE).
- Decomposition pie chart: NDE vs NIE.
- Counterfactual: "if we held M at baseline, what would the effect of X on Y be?"
- Outcome banner: "with these settings: total effect = 1.0; 70% direct; 30% via mediator."

### Target aha

Mediation decomposes total causal effect; choice of mediator matters for intervention design.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — Pleiotropy Detector

**File**: `artifacts/lecture-25/05-pleiotropy.html`
**Anchor**: §3.3

### Teaching purpose

Investigate horizontal pleiotropy in a synthetic MR. Apply MR-Egger and identify the systematic intercept.

### UI layout

- Pre-loaded scatter ($\beta_{ZX}$ vs $\beta_{ZY}$).
- Slider: pleiotropy magnitude (multiplied across all SNPs).
- IVW line + Egger line both drawn live.
- Egger intercept p-value displayed.
- Outcome banner: "pleiotropy 0: Egger intercept = 0.01, p = 0.4 → no pleiotropy. Pleiotropy 0.05: Egger intercept = 0.05, p = 0.001 → systematic pleiotropy detected."

### Target aha

MR-Egger detects systematic pleiotropy via the intercept; the slope is the pleiotropy-corrected causal effect.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — Drug-Target MR Explorer

**File**: `artifacts/lecture-25/06-drug-target-mr.html`
**Anchor**: §6

### Teaching purpose

Pick a drug target; artifact runs cis-MR using simulated cis-eQTLs.

### UI layout

- Drug target dropdown: PCSK9, IL6R, HMGCR, IL23R, etc.
- For each target: list of cis-eQTLs (SNP, effect size on target expression).
- MR runs against several outcome phenotypes (CAD, AMD, IBD, T2D, etc.).
- PheWAS-style heat map: target × phenotype causal-effect estimates.
- Outcome banner: "PCSK9 inhibition: significant LDL reduction → CAD reduction (β = 0.5). No effect on AMD or IBD. Confirms cardiovascular indication."

### Target aha

Cis-MR pre-validates drug indications and identifies side-effect risks before clinical trials.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — IV-Assumption Diagnostic

**File**: `artifacts/lecture-25/07-iv-diagnostic.html`
**Anchor**: §2.4

### Teaching purpose

Run an IV diagnostic battery on a candidate instrument set: F-statistic, exclusion-restriction proxies, balance checks.

### UI layout

- Pre-loaded IV set (~10 SNPs for an exposure).
- Diagnostic tests:
  - F-statistic per SNP and combined.
  - Population-stratification balance (by ancestry covariates).
  - Pleiotropy check via PheWAS.
  - Steiger directionality.
- Pass/fail per assumption.
- Outcome banner: "instrument set passes relevance (F = 47 > 10), independence (p > 0.4), and Steiger (p < 0.001). Exclusion restriction passes proxy tests."

### Target aha

Three IV assumptions; each has diagnostic checks; pre-screen instruments before MR.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- DAG colours consistent across #1, #2, and lecture figures.
- Pleiotropy colour scheme consistent across #3, #5 and lecture figures.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
