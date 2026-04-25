# Lecture 22 (proposed L10) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 22 (Network Biology and Pathway Analysis).
> **Companion files**: `lecture-22.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-22/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Network types: PPI cobalt, GRN red, metabolic green, signalling violet, disease teal.
- Enrichment colours: significant red, marginal amber, non-significant grey.
- Module colours: rotating palette (cobalt, amber, violet, teal, green, red).
- Typography: Inter for chrome; JetBrains Mono for gene names, p-values, scores.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget

Seven interactive tools.

| # | Title | Anchor |
|---|---|---|
| 1 | ORA enrichment calculator | §3.2 |
| 2 | GSEA running-enrichment explorer | §3.4 |
| 3 | Network propagation walker | §4.2 |
| 4 | Spectral clustering visualizer | §5.3 |
| 5 | Pathway database query | §2 |
| 6 | Drug-target prediction with RWR | §6.1 |
| 7 | Module detection (Louvain) | §5.4 |

---

## Artifact #1 — ORA Enrichment Calculator

**File**: `artifacts/lecture-22/01-ora-calculator.html`
**Anchor**: §3.2

### Teaching purpose

Compute a hypergeometric / Fisher's-exact p-value for over-representation of a gene set in a gene list. Adjust the universe; see how p-value changes.

### UI layout

- Three numeric inputs: gene-list size $|L|$, set size $|G|$, observed overlap $a$.
- Universe size slider $|U|$ (typically $\sim$ 20,000 for human).
- Output: hypergeometric p-value, log-odds, fold enrichment, FDR (with Bonferroni-adjusted estimate over 1000 hypotheses).
- Visualisation: 2x2 table updates live; bar chart of expected vs observed.
- Outcome banner: "$|L| = 200, |G| = 50, a = 12, |U| = 20000$ → p = 1e-7, fold enrichment = 24x. Highly enriched."

### Target aha

Hypergeometric: bigger overlap given size of $L$ and $G$ → smaller p-value. Universe size matters; gene-set background doesn't fix bias.

### Acceptance criteria

- Default opens with realistic numbers and significance.
- Slider on universe size moves p-value visibly.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — GSEA Running-Enrichment Explorer

**File**: `artifacts/lecture-22/02-gsea-explorer.html`
**Anchor**: §3.4

### Teaching purpose

Walk along a synthetic ranked gene list. See running enrichment score build up. Compare gene sets with strong, weak, no enrichment.

### UI layout

- Pre-loaded ranked gene list (~500 genes) with three example gene sets:
  - Strong (clustered at top of ranking).
  - Mixed (spread across ranking).
  - No enrichment (uniform).
- Per-set running ES curve plotted live.
- Slider: position cursor along the ranking — see ES value at that position.
- Permutation histogram showing the null distribution.
- Outcome banner: "set A: ES = 0.6, p = 0.001 — strong. set B: ES = 0.2, p = 0.5 — non-significant."

### Target aha

GSEA detects coherent shifts; non-clustered gene sets don't accumulate ES.

### Acceptance criteria

- Three gene sets demonstrate clear contrast.
- Permutation null computed in-browser.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Network Propagation Walker

**File**: `artifacts/lecture-22/03-rwr-walker.html`
**Anchor**: §4.2

### Teaching purpose

Run RWR on a small synthetic PPI. Watch heat propagate from seeds to network neighbours. Tune restart probability $r$.

### UI layout

- ~30-node network displayed as SVG with force-directed layout.
- Slider: restart probability $r$ (0 to 1).
- Step button: advance one iteration of RWR.
- Run-to-convergence button.
- Each node's RWR score visible as colour intensity.
- Outcome banner: "with $r = 0.5$ from 3 seeds: convergence in 12 iterations. Top scoring non-seed: gene XYZ."

### Target aha

Propagation spreads signal across network; restart probability controls the spread radius.

### Acceptance criteria

- Default seed configuration shows clear propagation pattern.
- Slider on $r$ visibly changes spread.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Spectral Clustering Visualizer

**File**: `artifacts/lecture-22/04-spectral-cluster.html`
**Anchor**: §5.3

### Teaching purpose

On a small network with known communities, compute the Laplacian, eigendecompose, project onto eigenvectors, k-means cluster.

### UI layout

- A 30-node network with 3 communities.
- Three tabs: "Network" | "Laplacian" | "Eigenspace".
- Network tab: standard force-directed layout.
- Laplacian tab: heat map of $L = D - A$.
- Eigenspace tab: 2D scatter on $\mathbf{v}_2, \mathbf{v}_3$; node colour by inferred cluster.
- k-slider for number of clusters.
- Outcome banner: "k=3: NMI vs ground truth = 0.95. k=5: NMI = 0.78 (over-segmented)."

### Target aha

Eigenvectors at low frequencies separate communities; k-means in eigenspace recovers them.

### Acceptance criteria

- Three-community network with clear ground truth.
- Eigendecomposition computed in JS (small matrices ~30×30).
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — Pathway Database Query

**File**: `artifacts/lecture-22/05-pathway-query.html`
**Anchor**: §2

### Teaching purpose

Pick a gene; see what pathways it's in across KEGG, Reactome, GO. Compare granularities.

### UI layout

- Gene dropdown (50 examples: TP53, BRCA1, MYC, EGFR, ...).
- Side-by-side panels for KEGG, Reactome, GO BP, GO MF, MSigDB Hallmark.
- Each panel: pathway count, list of top pathways the gene appears in.
- Highlight overlap: pathways present in multiple databases.
- Outcome banner: "TP53 in KEGG pathways: 14. Reactome: 87. GO BP terms: 198. The granularity differs — pick the database that matches your question."

### Target aha

Gene-set sources differ in granularity; same gene appears in different "pathways" across databases.

### Acceptance criteria

- 50 genes pre-loaded with realistic pathway annotations.
- Cross-database overlap visualised.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — Drug-Target Prediction with RWR

**File**: `artifacts/lecture-22/06-drug-target-rwr.html`
**Anchor**: §6.1

### Teaching purpose

Build a heterogeneous drug-target-disease network. Run RWR from one drug. Predict novel targets.

### UI layout

- Heterogeneous network with drugs, targets, diseases.
- Drug dropdown: pick a drug (~10 presets).
- RWR runs from drug's known targets.
- Output: top 10 predicted novel targets ranked by RWR score.
- Held-out validation: 5 of the predicted "novel" are actually known but excluded for the demo; user can toggle to compare.
- Outcome banner: "for drug imatinib, top novel target is FGR (RWR score 0.34); validated literature reports FGR as a known imatinib target."

### Target aha

RWR predicts plausible novel drug-target associations from network structure alone.

### Acceptance criteria

- Realistic toy network (~50 drugs, ~100 targets, ~30 diseases).
- Held-out validation works.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — Module Detection with Louvain

**File**: `artifacts/lecture-22/07-louvain.html`
**Anchor**: §5.4

### Teaching purpose

Run Louvain step by step on a small network. Watch communities emerge from local optimisation.

### UI layout

- ~50-node network with 4 ground-truth communities.
- Step-forward button advances one Louvain iteration (one node moved or one aggregation).
- Modularity Q displayed; updates after each step.
- Run-to-convergence button.
- Final partition coloured by community.
- Outcome banner: "Louvain converged in 8 iterations to Q = 0.87 with 4 communities matching ground truth."

### Target aha

Louvain greedily improves modularity; converges fast in practice; produces interpretable communities.

### Acceptance criteria

- Step controls work; modularity tracker updates correctly.
- Final partition visible.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Network-type colours consistent across #3, #6, #7 and lecture figures.
- Enrichment colours consistent across #1, #2 and lecture figures.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
