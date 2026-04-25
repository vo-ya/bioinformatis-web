# Lecture 22 (proposed L10) — Network Biology and Pathway Analysis

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology assumed
> **File**: to be rendered as `lectures/lecture-22.html` (provisional name; renumber to `lecture-10.html` when curriculum is reordered)

> **Proposed placement**: insert after L6 (differential expression, becomes new L9 in reordered curriculum). The natural arc: bulk RNA-seq → DE analysis → "you have a gene list, what does it mean?" → networks and pathways. Network biology gives DEGs (and GWAS hits, methylation outliers, single-cell markers — all from upstream lectures) somewhere to land. The lecture frames networks as graphs and pathway enrichment as a hypothesis-testing problem on an enrichment statistic; this connects directly back to the multiple-testing framework of L6 and forward to single-cell trajectory (L8) and GWAS prioritisation (L13).

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Distinguish the major biological-network types: protein-protein interaction (PPI), gene regulatory (GRN), metabolic, signalling, and disease networks; recognise the data sources behind each.
2. Apply pathway enrichment via over-representation analysis (ORA, hypergeometric / Fisher's exact test) and gene-set enrichment analysis (GSEA, signed Mann-Whitney / Kolmogorov-Smirnov on a ranked list).
3. Compute network propagation (random-walk-with-restart, heat diffusion) and explain it as a Laplacian smoother on the network.
4. Detect modules / communities via spectral clustering and modularity maximisation.
5. Interpret pathway databases (KEGG, Reactome, MSigDB, GO) and their hierarchical / DAG structure.
6. Apply network-based methods to drug-target prediction, disease-gene prioritisation, and pathway-context interpretation.
7. Frame network biology as graph signal processing: the network Laplacian is a spectral operator; propagation = low-pass filter; module detection = clustering on graph spectrum.

---

## Part 1 — Why Networks (≈25 min)

### 1.1 The gene-list problem (≈5 min)

By the time you get to this lecture, you've generated several gene lists in this course:

- DEGs from RNA-seq differential expression (L6).
- TF-bound regions from ChIP-seq (L9).
- Cell-state-marker genes from scRNA-seq clustering (L7).
- GWAS-significant SNPs and their nearest genes (L13).
- Tumour-driver mutations (L18).

In every case, the question **"what does this gene list mean?"** is the next step. Networks and pathways are the framework for answering it.

### 1.2 Network types in biology (≈8 min)

Several distinct network types, each with different data sources and uses:

**Protein-protein interaction (PPI) networks**:

- Nodes: proteins.
- Edges: physical interactions detected by affinity purification mass spec (AP-MS), yeast two-hybrid (Y2H), proximity labelling (BioID, APEX), structural prediction (AlphaFold-multimer).
- Major databases: **STRING** (broadest, includes predicted interactions), **BioGRID** (curated experimental), **IntAct**, **HuRI** (human reference interactome).
- Use: identify protein complexes, infer function by guilt-by-association.

**Gene regulatory networks (GRNs)**:

- Nodes: genes (or TFs and their targets).
- Edges: TF → gene transcriptional regulation.
- Sources: ChIP-seq motif enrichment, perturbation experiments (CRISPRi screens), ATAC-seq + chromatin annotation.
- Use: predict expression cascades; identify master regulators.

**Metabolic networks**:

- Nodes: metabolites; edges: enzymes (or vice versa).
- Sources: KEGG, MetaCyc, Recon3D.
- Use: pathway analysis, flux balance analysis, drug-side-effect prediction.

**Signalling networks**:

- Nodes: signalling proteins (kinases, GPCRs, etc.).
- Edges: phosphorylation, allosteric activation, complex formation.
- Sources: PhosphoSitePlus, OmniPath, Reactome.
- Use: drug-mechanism inference, pathway perturbation modelling.

**Disease networks**:

- Nodes: diseases.
- Edges: shared genes, comorbidity, drug-target overlap.
- Sources: DisGeNET, OpenTargets, OMIM.
- Use: drug repurposing, comorbidity prediction.

**FIGURE — Figure #1: Biological network taxonomy** → `diagrams/lecture-22/01-network-types.svg`
*Five panels, one per network type. Each shows: typical node and edge appearance, scale (number of nodes / edges in human), example database, and example application. Side panel: key cross-network workflows (PPI + GRN → "regulome", PPI + metabolic → "interactome").*

### 1.3 Network properties at a glance (≈5 min)

Real biological networks share common topology features:

- **Scale-free degree distribution**: few hub nodes (high degree); many low-degree nodes. Power-law $P(k) \sim k^{-\gamma}$.
- **Small-world**: average path length grows logarithmically with network size.
- **Modular**: dense clusters (functional modules) connected sparsely.

These are not unique to biology — social networks, the Internet, citation networks all show similar patterns. Network science as a field exploits this commonality.

### 1.4 Why graphs are the natural representation (≈4 min)

Genes and proteins don't function in isolation. Their behaviour depends on:

- Which other proteins they interact with.
- Which pathways they're embedded in.
- How perturbations propagate through their neighbourhood.

A list misses all of this; a graph captures it. Once you have the graph, the toolkit of graph algorithms — shortest paths, centrality, community detection, propagation — applies directly.

### 1.5 The deep dive (≈3 min)

> **EE framing — networks as graph signal processing**: A biological network is a **graph** $G = (V, E)$ with adjacency matrix $A$ and Laplacian $L = D - A$ (where $D$ is the degree matrix). A "signal" on the graph is a function $f: V \to \mathbb{R}$ (e.g., expression levels). The graph Fourier transform decomposes $f$ in the eigenbasis of $L$: low-frequency components are smooth across edges; high-frequency components vary rapidly. Network propagation (Part 4) is a **low-pass filter**: it smooths the input signal across the network. Spectral clustering (Part 5) operates in this same eigenbasis. The whole field of network biology can be re-cast as graph signal processing — a relatively young area in EE that bioinformatics implicitly invented in the 2000s.

---

## Part 2 — Pathway Databases and the Gene-Set Universe (≈25 min)

### 2.1 What's a pathway (≈4 min)

A **pathway** is a curated list of genes / proteins / metabolites that work together to accomplish a biological function:

- Glycolysis: ~20 enzymes converting glucose to pyruvate.
- p53 signalling: TP53 + ATM + MDM2 + CDKN1A + ~50 downstream targets.
- Wnt-β-catenin: ligands, receptors, intracellular cascade, transcriptional output.

Pathways are the unit of human-readable biological interpretation.

### 2.2 KEGG (≈5 min)

**KEGG (Kyoto Encyclopedia of Genes and Genomes)**, established 1995, ~550 manually-curated pathways.

- Metabolic pathway maps (the iconic "ball-and-stick" diagrams).
- Signalling, regulatory, disease pathways.
- Cross-species: pathways are abstract; genes mapped per species.

KEGG pathway IDs (e.g., hsa04110 = "Cell cycle - human") are the de facto standard for pathway-level annotation.

### 2.3 Reactome (≈5 min)

**Reactome**, established 2003, ~2,500 pathways (more granular than KEGG):

- Hierarchical: pathways nest into super-pathways and sub-pathways.
- Reaction-level detail: each pathway is a sequence of biochemical reactions.
- Strong coverage of human signalling and immune pathways.

Reactome's hierarchy means you can run analysis at the granularity that fits — e.g., "DNA repair" (broad) or "Translesion synthesis" (specific).

### 2.4 Gene Ontology (GO) (≈5 min)

**GO**, established 1998, ~50,000 terms organised as three DAGs:

- **Biological Process (BP)**: e.g., "cell cycle progression", "immune response".
- **Molecular Function (MF)**: e.g., "kinase activity", "DNA binding".
- **Cellular Component (CC)**: e.g., "nucleus", "plasma membrane".

GO is the most widely-used annotation system in the world. Each gene is annotated with multiple GO terms, often dozens.

The DAG structure means terms have parent-child hierarchy: "cell cycle phase transition" ⊂ "cell cycle" ⊂ "cellular process" ⊂ "biological process".

### 2.5 MSigDB and other compendia (≈3 min)

**MSigDB** (Molecular Signatures Database, Broad Institute): aggregates many pathway / gene-set sources. ~30,000 gene sets across:

- Hallmark gene sets (50 well-curated functional collections).
- KEGG, Reactome, BioCarta translations.
- Curated from literature.
- Cancer-specific signatures.
- Cell-type marker genes.

For typical RNA-seq DE analysis, MSigDB is the first place to test enrichment.

### 2.6 The deep dive (≈3 min)

The **gene-set universe** has rich structure: genes belong to many overlapping sets; sets have hierarchical relationships; some are well-curated, others statistical aggregates. Modern analyses must account for this overlap (similar terms often produce correlated p-values).

**FIGURE — Figure #2: Pathway database landscape** → `diagrams/lecture-22/02-pathway-databases.svg*
*A grid showing each major database (KEGG, Reactome, GO BP/MF/CC, MSigDB) with: scale, granularity, hierarchy structure, typical use case. Bottom panel: example query — "What pathway is TP53 in?" — and the heterogeneous answers from each database. Side annotation: "Different databases give different answers; modern analyses run multiple in parallel."*

---

## Part 3 — Pathway Enrichment Analysis (≈45 min)

### 3.1 The over-representation problem (≈4 min)

Setup: you have

- $L$ = your gene list (e.g., DEGs at FDR < 0.05).
- $G$ = a gene set (e.g., a KEGG pathway).
- $U$ = the universe (typically all genes tested in the experiment).

Question: are genes from $G$ over-represented in $L$ relative to expectation under the null hypothesis "$L$ is a random subset of $U$"?

### 3.2 Hypergeometric / Fisher's exact test (≈8 min)

The natural null: $L$ is sampled uniformly without replacement from $U$. Under this null, the count of $L \cap G$ follows a **hypergeometric distribution**.

Fisher's exact test (mathematically equivalent in this setting):

| | In $G$ | Not in $G$ | Total |
|---|---|---|---|
| In $L$ | $a$ | $b$ | $a + b$ |
| Not in $L$ | $c$ | $d$ | $c + d$ |
| Total | $a + c$ | $b + d$ | $N$ |

P-value = probability of observing $\geq a$ overlap by chance. Computed exactly via the hypergeometric CDF.

This is **ORA (Over-Representation Analysis)** — the simplest pathway-enrichment method.

### 3.3 The multiple-testing problem (≈5 min)

If you test 1000 pathways, you'll get many "significant" by chance. Apply Benjamini-Hochberg FDR (from L6) to the per-pathway p-values. Cap reported pathways at FDR < 0.05.

For overlapping gene sets (e.g., parent and child GO terms), the p-values are correlated, so BH is conservative but valid.

### 3.4 GSEA (≈12 min)

**Gene Set Enrichment Analysis** (Subramanian et al. 2005) addresses ORA's biggest weakness: ORA throws away the rank information.

ORA: only counts genes above the significance threshold.
GSEA: uses the full ranked list of all tested genes (e.g., ranked by log fold change or t-statistic).

Algorithm:

1. Rank all $N$ genes by some statistic (e.g., signed log fold change).
2. For each gene set $G$:
   - Walk down the ranked list.
   - At each gene in $G$, increment a running enrichment score (ES).
   - At each gene not in $G$, decrement.
   - Track the maximum |ES| reached — this is the **enrichment score** for $G$.
3. Compute statistical significance via permutation: shuffle gene labels and re-compute ES; the p-value is the fraction of permutations with $|ES'| \geq |ES|$.

GSEA detects coherent shifts in pathway expression even when individual genes don't pass the strict significance threshold. Far more sensitive than ORA on real data.

### 3.5 Modern variants (≈6 min)

- **fgsea** (Sergushichev 2016): a fast preranked GSEA.
- **camera** (Wu & Smyth 2012): variance-inflation correction for inter-gene correlation.
- **roast / mroast** (Wu et al. 2010): rotation-based gene-set tests.
- **SAFE / SCREENER / CERNO**: alternative running-statistic methods.

For typical RNA-seq workflows: run **fgsea** with MSigDB Hallmark + KEGG + Reactome, FDR < 0.05.

### 3.6 The signed vs unsigned distinction (≈4 min)

If the ranked statistic is signed (e.g., log fold change), GSEA detects **directional** enrichment (pathway up vs down). If unsigned (absolute value), GSEA detects **deregulation** (pathway perturbed in either direction).

Choice depends on the question:

- "Is glycolysis up in tumour vs normal?" → signed.
- "Is glycolysis perturbed (up or down) by drug X?" → unsigned.

### 3.7 The deep dive (≈6 min)

> **EE framing — GSEA as a Mann-Whitney / Kolmogorov-Smirnov test on signed rankings**: GSEA's running enrichment score is the difference between two CDFs: the CDF of "genes in $G$" and the CDF of "all genes" along the ranked axis. The maximum deviation is the **Kolmogorov-Smirnov statistic** for the two-sample test. The signed version is the **signed-rank Mann-Whitney**: for each gene in $G$, sum its rank; compare to expected sum under the null. This is well-known statistical machinery; GSEA's contribution was the running-statistic visualisation and the permutation-test calibration that's robust to gene-gene correlation. The Subramanian paper is one of the most-cited bioinformatics papers ever (~30,000 citations) — partly because the visualisation made the method instantly interpretable.

**FIGURE — Figure #3: ORA vs GSEA** → `diagrams/lecture-22/03-ora-vs-gsea.svg`
*Left panel: ORA — a 2x2 contingency table with the hypergeometric / Fisher's-exact p-value computation visualised. Right panel: GSEA — a ranked list of all genes; tick marks for genes in the gene set; running enrichment score curve overlaid; max ES highlighted. Bottom annotation: "ORA throws away rank; GSEA uses the full list. GSEA detects coherent shifts ORA misses."*

---

## Part 4 — Network Propagation (≈30 min)

### 4.1 The smoothing intuition (≈4 min)

If gene $g$ is connected to many DEGs in the PPI network, $g$ is likely involved even if its own expression didn't change. **Network propagation** formalises this: spread the DEG signal across the network so neighbours of DEGs accumulate score.

Used in:

- Disease-gene prioritisation: spread known-disease genes; rank novel candidates by accumulated score.
- Drug-target prediction: spread known drug targets; rank candidates.
- GWAS interpretation: spread GWAS hits to identify implicated pathways.
- Cancer driver discovery: spread mutations across the PPI to find context-driving modules.

### 4.2 Random walk with restart (RWR) (≈8 min)

The most common propagation algorithm:

$$\mathbf{p}^{(t+1)} = (1 - r) W \mathbf{p}^{(t)} + r \mathbf{p}^{(0)}$$

where $W$ is a normalised adjacency matrix (typically column-normalised to be a transition matrix), $r$ is the restart probability (typical 0.5), and $\mathbf{p}^{(0)}$ is the initial signal (e.g., 1 for known-disease genes, 0 for everything else).

At each step:

- With probability $1 - r$, walk to a neighbour.
- With probability $r$, restart at $\mathbf{p}^{(0)}$.

Iterate to convergence (~50 iterations). The stationary distribution $\mathbf{p}^{(\infty)}$ gives each node a score reflecting its proximity to the seeds in the network.

Closed-form: $\mathbf{p}^{(\infty)} = r (I - (1-r) W)^{-1} \mathbf{p}^{(0)}$.

### 4.3 Heat diffusion (≈4 min)

An alternative formulation: model the network signal as a heat distribution that diffuses according to the graph Laplacian:

$$\frac{d\mathbf{p}}{dt} = -L \mathbf{p}$$

with $L$ the normalised Laplacian. Solution: $\mathbf{p}(t) = e^{-Lt} \mathbf{p}(0)$.

The continuous heat kernel $e^{-Lt}$ is the analog of RWR. Both are **low-pass filters** in the graph Fourier basis.

### 4.4 What propagation does, formally (≈6 min)

In the eigenbasis of $L$ (with eigenvalues $\lambda_i$ and eigenvectors $\mathbf{v}_i$):

$$\mathbf{p}(t) = \sum_i e^{-\lambda_i t} \langle \mathbf{p}(0), \mathbf{v}_i \rangle \mathbf{v}_i$$

High-eigenvalue components decay fast; low-eigenvalue components persist. The result: the propagated signal is **smooth** across the graph (low-frequency).

This is the same low-pass filter logic as Gaussian smoothing in image processing, applied to a graph instead of a grid.

### 4.5 PRINCE, HotNet, and the workflow (≈4 min)

**PRINCE** (Vanunu et al. 2010): RWR for disease-gene prioritisation; the canonical implementation.

**HotNet2** (Leiserson et al. 2015): heat-diffusion-based discovery of "hot subnetworks" enriched in cancer mutations.

**NetWAS**: GWAS hit prioritisation via network propagation.

Practical: load network → set seeds → run propagation (50 iterations or matrix inversion) → rank nodes by score → interpret top hits.

### 4.6 The deep dive (≈4 min)

> **EE framing — network propagation as Laplacian smoothing**: RWR and heat diffusion are both **low-pass filters in the graph Fourier domain**. The graph Laplacian $L$ is the spectral operator: smooth signals (slowly varying across edges) live in low-eigenvalue eigenvectors; rough signals (rapidly varying) live in high-eigenvalue ones. Diffusing a signal for time $t$ multiplies its eigenvector decomposition by $e^{-\lambda_i t}$; the high-frequency components are damped. This is exactly the spatial smoothing intuition — propagation enforces "neighbours should have similar scores" in the same way image smoothing enforces "adjacent pixels should be similar". The connection to graph signal processing is direct.

**FIGURE — Figure #4: Network propagation as low-pass filtering** → `diagrams/lecture-22/04-network-propagation.svg`
*Top: a small PPI network with a few seed nodes (red, intensity 1) and the rest (white, intensity 0). Middle: after 5 iterations of RWR, intensity has spread to direct neighbours. Bottom: after 50 iterations, the signal has reached steady state — high near seeds, low far away. Side panel: a graph-signal-processing analogy showing the same effect as Gaussian smoothing on an image grid.*

---

## Part 5 — Module Detection and Communities (≈25 min)

### 5.1 The module hypothesis (≈4 min)

A **module** is a tightly-connected subgraph corresponding to a functional unit (a protein complex, a pathway, a cell-type-specific regulatory module). Detecting modules from a network is **community detection**.

The biological assumption: genes in the same module have related functions and are more likely to be co-regulated, co-expressed, and disease-associated.

### 5.2 Modularity (≈5 min)

The Newman modularity score (Newman 2004) quantifies how strongly partition $\{C_1, ..., C_k\}$ corresponds to communities:

$$Q = \frac{1}{2m} \sum_{i,j} \left[ A_{ij} - \frac{k_i k_j}{2m} \right] \delta(c_i, c_j)$$

where $m$ = total edges, $k_i$ = degree of node $i$, $\delta(c_i, c_j) = 1$ if $i$ and $j$ are in the same community, 0 otherwise.

Maximising $Q$ produces the optimal community partition. The Louvain algorithm and Leiden algorithm (Traag et al. 2019) are the standard fast heuristics.

### 5.3 Spectral clustering (≈8 min)

An alternative: eigendecompose the graph Laplacian $L$, then cluster on the smallest non-zero eigenvectors.

Algorithm:

1. Compute $L = D - A$.
2. Find $k$ smallest non-zero eigenvalues with eigenvectors $\mathbf{v}_1, ..., \mathbf{v}_k$.
3. Form an $n \times k$ matrix $V$.
4. Cluster rows of $V$ with k-means.

The eigenvectors at low frequencies highlight **cuts** of the graph that separate dense communities. Spectral clustering is mathematically equivalent to graph-cut optimisation and produces similar results to modularity maximisation.

### 5.4 Louvain and Leiden (≈4 min)

**Louvain** (Blondel et al. 2008): iterative modularity maximisation:

- Initialise each node as its own community.
- For each node, move to the neighbouring community that increases $Q$ most.
- Aggregate communities; repeat at the meta-level.

**Leiden** (Traag et al. 2019): improves Louvain by guaranteeing connectivity within communities (Louvain occasionally produces disconnected modules). Used in scanpy / Seurat for single-cell clustering (L7).

### 5.5 The deep dive (≈4 min)

> **EE framing — community detection as graph-spectrum clustering**: Both modularity maximisation (Louvain / Leiden) and spectral clustering operate on the **graph spectrum**. The smallest-eigenvalue eigenvectors of $L$ encode the slowest-varying signals on the graph; these are the natural "split directions" — the analogs of Fourier components that highlight the graph's natural cuts. The connection between modularity $Q$ and spectral clustering was made rigorous in the early 2010s. They give **similar results in practice** for biological networks; the difference is mostly speed (Louvain/Leiden scale better) and theoretical guarantees (spectral has a tighter cut-quality bound).

**FIGURE — Figure #5: Spectral clustering of a network** → `diagrams/lecture-22/05-spectral-clustering.svg`
*Top: a small network with three visually-apparent dense clusters connected sparsely. Middle: 2-D projection onto the second and third smallest Laplacian eigenvectors; the three clusters separate cleanly. Bottom: k-means partition in the eigenspace returns the correct community assignment. Side panel: comparison to Louvain, which produces the same partition.*

---

## Part 6 — Applications: Drugs, Disease, and Beyond (≈20 min)

### 6.1 Drug-target prediction (≈5 min)

Network-based drug-target prediction:

- Build a heterogeneous network: drugs, targets, diseases, side effects.
- Use random-walk-with-restart from drugs of known efficacy or known disease genes.
- Predict new drug-target relationships from high RWR scores between drug nodes and protein nodes.

Methods: **GBA (Guilt By Association)**, **NeoDTI** (deep-learning extension), **DrugBank-based RWR**.

In practice: network methods are useful for **drug repurposing** (finding new uses for existing drugs) and **off-target prediction** (anticipating safety issues before clinical trials).

### 6.2 Disease-gene prioritisation (≈5 min)

Given a chromosomal region (e.g., from L13 GWAS or L17 clinical exome), several genes lie in the candidate region. Which is most likely the disease gene?

Network-based approach:

- Take a known set of seed disease genes for the same disease.
- Run RWR from the seeds.
- Among candidate genes in the region, the one with the highest RWR score is the most likely culprit.

Tools: **PRINCE** (Vanunu 2010), **DOMINO** (Levi 2021), **GeneMania** (Warde-Farley 2010). These are routine in clinical genomics for novel disease-gene discovery.

### 6.3 Pathway-context interpretation of DEGs (≈4 min)

When you have 200 DEGs from an RNA-seq experiment, ORA and GSEA tell you which pathways are enriched. Network methods take this further:

- Project DEGs onto the PPI network.
- Run propagation to find connected modules among DEGs.
- Identify "core" hubs that connect multiple DEGs.
- Hypothesise causal mechanism: hub gene's perturbation explains the cascade.

This is the workflow behind tools like **OmicsNet**, **NetworkAnalyst**, **STRING-driven analysis**.

### 6.4 Cell-cell communication (≈3 min)

A growing application: from single-cell data (L8), infer cell-cell signalling networks based on ligand-receptor expression patterns and downstream effector activation.

Tools: **CellChat** (Jin 2021), **CellPhoneDB**. The signalling network is inferred per cell-type pair; integrative analysis identifies dominant communication pathways in tissue.

### 6.5 The 2024 frontier (≈3 min)

- **Graph neural networks** (GNNs, L16) are starting to replace classical RWR for some tasks. They learn task-specific node embeddings instead of using fixed propagation.
- **Multi-modal networks**: integrating PPI + GRN + drug-target + GWAS into one heterogeneous graph; GNNs learn unified embeddings.
- **Network medicine**: applying network biology to clinical decision support; still mostly research, slowly entering clinical practice.

**EMBED — Artifact #6: Drug-target prediction with RWR** → `artifacts/lecture-22/06-drug-target-rwr.html`
*Pick a drug; artifact runs RWR from its known targets and ranks novel candidates. Compare to ground truth.*

---

## Part 7 — Tools, Visualisation, and Workflows (≈15 min)

### 7.1 Visualisation tools (≈4 min)

- **Cytoscape**: the desktop standard for network visualisation. Plugins for layout, statistics, integration.
- **Gephi**: another popular desktop tool.
- **igraph** (R / Python): scriptable analysis.
- **networkx** (Python): the de facto standard for programmatic network analysis.
- **DG-graphs** (R Bioconductor): for biological-network-specific manipulation.

For a typical RNA-seq paper: igraph / networkx for analysis; Cytoscape for figure generation.

### 7.2 The standard workflow (≈4 min)

1. **Define seeds**: DEGs / GWAS hits / cancer drivers / etc.
2. **Pick network**: STRING (broad), BioGRID (curated), Reactome (pathway-based).
3. **Filter network**: keep highest-confidence interactions (STRING confidence > 0.7).
4. **Run propagation / community detection**.
5. **Visualise**: subgraph induced by seeds + their high-score neighbours.
6. **Annotate**: pathway enrichment on each module.

### 7.3 Common pitfalls (≈4 min)

- **Network bias**: the network has more edges among well-studied proteins; ranks may reflect study-bias rather than biology.
- **False edges**: high-throughput interaction screens have false-positive rates 30-50%. STRING's score thresholds matter.
- **Tissue specificity**: PPI networks are species-aggregate, not tissue-specific. Tissue-specific filtering improves relevance.
- **Hub bias**: hubs (e.g., TP53) appear in every analysis as significant; correct for hub centrality before interpreting.

### 7.4 Hands-on exercise (≈3 min)

For a typical analysis:

1. Get DEGs from one of the L6 exercises (or use a public list).
2. Upload to STRING; threshold confidence to 0.7.
3. Apply MCL clustering or Louvain.
4. For each module, run pathway enrichment via STRING's built-in tool.
5. Identify the top 1-2 functional modules.

This takes ~15 minutes and is the standard "first network analysis" of any DE result.

**FIGURE — Figure #6: Network analysis workflow** → `diagrams/lecture-22/06-workflow.svg`
*Top-to-bottom flowchart: gene list → network selection → propagation / clustering → pathway enrichment → visualisation → biological hypothesis. Each step annotated with tool choices and runtime estimates. Side panel: variant for GWAS-driven workflow (NetWAS / DOMINO).*

---

## Wrap-up (≈10 min)

### What you should take away

- **Networks are the natural representation** for biology beyond gene lists. Five major types: PPI, GRN, metabolic, signalling, disease.
- **Pathway enrichment** comes in two flavours: ORA (hypergeometric on a thresholded gene list) and GSEA (Mann-Whitney / KS on a ranked list). GSEA is more sensitive for coherent shifts.
- **Network propagation** = RWR or heat diffusion. Mathematically a low-pass filter in the graph spectrum. Used for disease-gene prioritisation, drug-target prediction, GWAS interpretation.
- **Module detection** = modularity maximisation (Louvain / Leiden) or spectral clustering. Both operate on the graph spectrum.
- **Pathway databases** (KEGG, Reactome, GO, MSigDB) differ in granularity and coverage; modern analyses run multiple in parallel.
- **EE framings**: networks as graphs with Laplacian; propagation as low-pass filter; module detection as spectral clustering; GSEA as KS test on ranked statistic.

### Next lecture

Single-cell RNA-seq fundamentals (existing L7, becomes new L11). The Leiden clustering used in scanpy / Seurat is the same Leiden algorithm you saw here.

### Homework

1. Take a list of 200 DEGs (use any public RNA-seq DE result). Run ORA against MSigDB Hallmark, KEGG, Reactome via clusterProfiler or g:Profiler. Report the top 5 enriched pathways at FDR < 0.05.
2. Run GSEA on the same data using the full ranked list. Compare to ORA: how many pathways agree? Where does GSEA find pathways ORA misses?
3. Build a small (~50-node) PPI subgraph from STRING for a gene family of your choice. Run RWR from one seed gene; rank all nodes by RWR score. Compare top 10 to direct neighbours of the seed.
4. Implement Louvain community detection in Python on a synthetic network with known communities. Compare to ground truth using normalised mutual information.
5. For one disease of interest, get the DisGeNET disease genes. Run propagation in the STRING network. Identify novel candidates (high RWR score, not yet associated with the disease).

### Recommended reading

- Subramanian, A., Tamayo, P., Mootha, V. K., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102, 15545–15550.
- Khatri, P., Sirota, M., & Butte, A. J. (2012). Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology* 8, e1002375.
- Vanunu, O., Magger, O., Ruppin, E., et al. (2010). Associating genes and protein complexes with disease via network propagation. *PLoS Computational Biology* 6, e1000641. (PRINCE.)
- Leiserson, M. D., Vandin, F., Wu, H. T., et al. (2015). Pan-cancer network analysis identifies combinations of rare somatic mutations across pathways and protein complexes. *Nature Genetics* 47, 106–114. (HotNet2.)
- Newman, M. E. (2006). Modularity and community structure in networks. *PNAS* 103, 8577–8582.
- Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9, 5233.
- Cowen, L., Ideker, T., Raphael, B. J., & Sharan, R. (2017). Network propagation: a universal amplifier of genetic associations. *Nature Reviews Genetics* 18, 551–562.
- Barabási, A. L., Gulbahce, N., & Loscalzo, J. (2011). Network medicine: a network-based approach to human disease. *Nature Reviews Genetics* 12, 56–68.
- STRING: <https://string-db.org/>
- KEGG: <https://www.kegg.jp/>
- Reactome: <https://reactome.org/>
- MSigDB: <https://www.gsea-msigdb.org/gsea/msigdb/>
- Gene Ontology: <http://geneontology.org/>
- Cytoscape: <https://cytoscape.org/>
- networkx: <https://networkx.org/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Why Networks                                       | 25 min | 0:25 |
| Part 2 — Pathway Databases and the Gene-Set Universe         | 25 min | 0:50 |
| Part 3 — Pathway Enrichment Analysis                         | 45 min | 1:35 |
| Part 4 — Network Propagation                                  | 30 min | 2:05 |
| Part 5 — Module Detection and Communities                     | 25 min | 2:30 |
| Part 6 — Applications: Drugs, Disease, and Beyond              | 20 min | 2:50 |
| Part 7 — Tools, Visualisation, and Workflows                   | 15 min | 3:05 |
| Wrap-up                                                          | 10 min | 3:15 |

**Total:** ~3h 15min of content.
