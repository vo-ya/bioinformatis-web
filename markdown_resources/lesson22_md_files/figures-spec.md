# Lecture 22 (proposed L10) — Figures Specification

> **Scope**: Static diagrams for Lecture 22 (Network Biology and Pathway Analysis).
> **Companion files**: `lecture-22.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Network types: PPI cobalt, GRN red, metabolic green, signalling violet, disease teal.
- Enrichment colours: significant red, marginal amber, non-significant grey.
- Module colours: stable rotation per module (cobalt, amber, violet, teal, green, red).
- Typography: Inter for chrome; JetBrains Mono for gene names, p-values, scores.

## Figure budget

Twelve figures.

| # | Title | Part |
|---|---|---|
| 1 | Biological network taxonomy | Part 1 |
| 2 | Pathway database landscape | Part 2 |
| 3 | ORA vs GSEA | Part 3 |
| 4 | Network propagation as low-pass filtering | Part 4 |
| 5 | Spectral clustering of a network | Part 5 |
| 6 | Network analysis workflow | Part 7 |
| 7 | KEGG glycolysis pathway map | Part 2 |
| 8 | GO DAG hierarchy | Part 2 |
| 9 | RWR in a heterogeneous drug-target network | Part 6 |
| 10 | GSEA running enrichment statistic | Part 3 |
| 11 | Disease-gene module on PPI | Part 6 |
| 12 | Hub-centrality bias visualisation | Part 7 |

---

## Figure 1 — Biological network taxonomy

**File**: `diagrams/lecture-22/01-network-types.svg`
**ViewBox**: `0 0 1200 720`

Five panels, one per network type:

- PPI (cobalt): nodes proteins, edges physical interactions, scale ~20k nodes / ~600k edges in human.
- GRN (red): nodes TFs/genes, edges regulatory; scale ~1500 TFs / 100k+ predicted edges.
- Metabolic (green): nodes metabolites, edges enzymes; KEGG-scale.
- Signalling (violet): nodes signalling proteins, edges activations; Reactome-curated.
- Disease (teal): nodes diseases, edges shared genes / comorbidity.

Each panel shows: typical visual style, scale, example database.

---

## Figure 2 — Pathway database landscape

**File**: `diagrams/lecture-22/02-pathway-databases.svg`
**ViewBox**: `0 0 1200 720`

A grid:

- KEGG: ~550 pathways, manually curated, metabolic + signalling + disease.
- Reactome: ~2,500 pathways, hierarchical.
- GO BP: ~30,000 terms, DAG structure.
- GO MF: ~12,000 terms.
- GO CC: ~4,000 terms.
- MSigDB: ~30,000 gene sets, aggregator.

Bottom panel: example query "What pathway is TP53 in?" with each database's heterogeneous answer.

Side annotation: "Different databases give different answers; modern analyses run multiple in parallel."

---

## Figure 3 — ORA vs GSEA

**File**: `diagrams/lecture-22/03-ora-vs-gsea.svg`
**ViewBox**: `0 0 1200 540`

Left panel — ORA:

- 2x2 contingency table: in $L$ × in $G$.
- Numbers: $a = 12$, $b = 38$, $c = 100$, $d = 850$.
- Hypergeometric p-value: $p = 0.0023$.

Right panel — GSEA:

- Ranked gene list along x-axis (left = most up-regulated, right = most down-regulated).
- Tick marks below for "in $G$" hits, distributed near the top of the ranking.
- Running enrichment-score curve overlaid; max ES = 0.6 marked.
- Permutation histogram showing null distribution; observed |ES| in the 99th percentile.
- p-value: $p = 0.001$.

Bottom annotation: "ORA throws away rank; GSEA uses the full list. GSEA detects coherent shifts ORA misses."

---

## Figure 4 — Network propagation as low-pass filtering

**File**: `diagrams/lecture-22/04-network-propagation.svg`
**ViewBox**: `0 0 1200 720`

Top: small PPI network (~30 nodes) with 3 seed nodes (red intensity 1) and the rest (white, 0).

Middle: after 5 RWR iterations — direct neighbours have intensity ~0.4; second-degree neighbours ~0.1.

Bottom: after 50 iterations (steady state) — exponential decay with graph distance from seeds.

Side panel: graph-signal-processing analogy. Show:

- Original signal (impulse at seed nodes).
- After "low-pass filter" (graph diffusion): smoothed across edges.
- Compare to image Gaussian smoothing.

Annotation: "Propagation = low-pass filter in the graph Fourier basis."

---

## Figure 5 — Spectral clustering of a network

**File**: `diagrams/lecture-22/05-spectral-clustering.svg`
**ViewBox**: `0 0 1200 720`

Top: a 3-community network (~30 nodes) with three visually-apparent dense clusters connected sparsely.

Middle: 2D scatter plot — eigenvectors $\mathbf{v}_2$ vs $\mathbf{v}_3$ of the Laplacian. Each node's coordinates are its components in the two slowest non-zero eigenvectors. Three clean clusters visible.

Bottom: k-means partition in eigenspace recovers the correct community assignment.

Side panel: comparison to Louvain — same partition, different algorithm.

---

## Figure 6 — Network analysis workflow

**File**: `diagrams/lecture-22/06-workflow.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart:

1. Gene list (DEGs / GWAS hits / cancer drivers).
2. Network selection (STRING / BioGRID / Reactome).
3. Filter (confidence > 0.7).
4. Run propagation / community detection.
5. Pathway enrichment.
6. Visualise (Cytoscape).
7. Biological hypothesis.

Each step annotated with tool choice and runtime. Side panel: variant for GWAS-driven workflow (NetWAS, DOMINO).

---

## Figure 7 — KEGG glycolysis pathway map

**File**: `diagrams/lecture-22/07-kegg-glycolysis.svg`
**ViewBox**: `0 0 1200 720`

Stylised KEGG-style pathway map of glycolysis:

- Nodes: metabolites (glucose, glucose-6-phosphate, fructose-6-phosphate, ..., pyruvate).
- Edges: enzymes (hexokinase, phosphofructokinase, ..., pyruvate kinase).
- Branches: pentose-phosphate-pathway entry, lactate output, mitochondrial entry.

Annotation: "KEGG's iconic ball-and-stick diagrams are the standard for human-readable pathway visualisation."

---

## Figure 8 — GO DAG hierarchy

**File**: `diagrams/lecture-22/08-go-dag.svg`
**ViewBox**: `0 0 1200 720`

A subset of the Biological Process DAG:

- Root: "biological process".
- Several deep paths: ... → "cellular process" → "cell cycle" → "cell cycle phase transition" → "G1/S transition".
- ... → "metabolic process" → "carbohydrate metabolism" → "glycolysis".

Some terms have multiple parents (DAG, not tree). Visualise this with crossing edges.

Side panel: number of GO terms at each level of the hierarchy.

---

## Figure 9 — RWR in a heterogeneous drug-target network

**File**: `diagrams/lecture-22/09-drug-target-rwr.svg`
**ViewBox**: `0 0 1200 720`

A small heterogeneous graph:

- Nodes: drugs (squares), proteins (circles), diseases (triangles).
- Edges: drug-target (red), protein-protein (cobalt), drug-disease (amber), gene-disease (teal).

Run RWR from a known drug:

- Iteration 0: only drug seed has intensity 1.
- Iteration 50: heat spreads to direct targets, then to PPI neighbours, then to other drug-disease relationships.

Annotation: "novel drug-target candidates rank high in the RWR steady state."

---

## Figure 10 — GSEA running enrichment statistic

**File**: `diagrams/lecture-22/10-gsea-running.svg`
**ViewBox**: `0 0 1200 540`

A ranked gene list with ~5000 genes along x-axis.

- Genes in target gene set marked as black tick marks below the axis.
- Running enrichment score curve overlaid (signed; positive when accumulating gene-set hits, negative when accumulating non-hits).
- Max |ES| highlighted.
- Side panel: leading-edge subset (the genes contributing to ES at the max).

Annotation: "the running statistic = signed Mann-Whitney rank-sum walk = KS-style maximum deviation."

---

## Figure 11 — Disease-gene module on PPI

**File**: `diagrams/lecture-22/11-disease-module.svg`
**ViewBox**: `0 0 1200 720`

A subgraph extracted from STRING, showing:

- Known cardiomyopathy genes (highlighted red).
- Their high-confidence (STRING > 0.9) neighbours.
- Module structure visible — most known genes cluster together; a few outlier genes connect via hubs.

Side panel: pathway enrichment of the module — top 5 pathways (sarcomere assembly, contraction, etc.).

Annotation: "disease genes form a connected module on the PPI; novel candidates rank high by network proximity."

---

## Figure 12 — Hub-centrality bias visualisation

**File**: `diagrams/lecture-22/12-hub-bias.svg`
**ViewBox**: `0 0 1200 540`

Two histograms:

- Left: degree distribution of TP53 (a hub: degree ~250).
- Right: degree distribution of a representative non-hub protein (degree ~5).

Below: P-value of "TP53 enriched in this gene list" computed naively (p < 10⁻¹⁰) vs after hub-centrality correction (p ~ 0.05).

Annotation: "Hubs appear in every analysis. Always correct for hub centrality before claiming biological significance."
