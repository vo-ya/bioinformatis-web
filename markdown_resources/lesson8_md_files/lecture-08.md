# Lecture 8 — Advanced Single-Cell: Trajectories, Integration, Multi-Modal

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-08.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Explain what pseudotime is, why it can be estimated from a single snapshot of cells, and state the assumption that makes it work.
2. Compare three trajectory-inference approaches (Monocle3, Slingshot, PAGA) and pick the right one for linear vs branching vs disconnected topologies.
3. Describe the RNA velocity dynamical model and derive the velocity estimator from spliced/unspliced ratios as a state-space problem.
4. Explain why scRNA-seq datasets from different labs/batches cannot be naively concatenated, and compare Harmony (linear correction in PCA space) to scVI (VAE on counts with NB likelihood).
5. State how scVI's variational ELBO decomposes, relate it to iterative soft-decision decoding, and connect it to the Salmon EM thread from Lecture 5.
6. Describe CITE-seq and scATAC-seq as additional modalities on the same cells, and name a joint-embedding method (WNN, totalVI, MOFA+).
7. Contrast Visium, MERFISH, and Xenium spatial platforms on resolution, throughput, and panel size; describe the mixed-pixel problem and one deconvolution approach.
8. Describe ligand-receptor inference (CellChat, NicheNet) and name two reasons the resulting communication scores are interpretive rather than mechanistic.

---

## Part 1 — Pseudotime and Trajectory Inference (≈30 min)

### 1.1 Why pseudotime from a static snapshot (≈5 min)

scRNA-seq captures a single snapshot in time — the cells are dissociated, lysed, sequenced, and gone. You cannot follow an individual cell over time. So how do we talk about "developmental progression" or "differentiation trajectories" from a snapshot?

The answer is pseudotime. The assumption: at any instant, a population of differentiating cells contains cells at every point along the trajectory — some just starting, some partway through, some almost done. If you assume the trajectory is continuous and that the population is asynchronous, the snapshot is effectively a time-series where every sample is taken at the same wall-clock instant but at different trajectory positions.

**Pseudotime** is the scalar coordinate along the trajectory. A cell's pseudotime value is "how far along the trajectory it is," in arbitrary units from 0 (start) to 1 (end) — or from "start cell" to "leaf cell" in graph units.

> **Intuition box**: Think of a crowded escalator viewed from a photograph. You cannot follow any individual person up the escalator — but the escalator is long enough that every step has a person on it, and you can read the whole progression from the single photograph. Pseudotime treats the cell population the same way: a dense sample of an unobserved process, with cells at every stage simultaneously.

### 1.2 Monocle, Slingshot, and PAGA (≈15 min)

Three dominant pseudotime methods in 2024, each with different assumptions:

**Monocle3** (Trapnell et al. 2014, Cao et al. 2019). Builds a principal-graph skeleton (via reverse-graph-embedding or `SimplePPT` algorithm) through the UMAP, then projects every cell to its nearest skeleton node. Pseudotime is graph distance from a designated "root" node. Handles branches natively (the skeleton is a tree / DAG). Best for datasets where the trajectory topology is unknown a priori.

**Slingshot** (Street et al. 2018). Runs clustering first (e.g. Leiden from Lecture 7), then connects clusters into a minimum-spanning tree, and fits a smooth principal curve through each lineage. Pseudotime is arc-length along the principal curve. Best when clusters are distinct and the trajectory topology is roughly tree-shaped; gives smoother, more interpretable curves than Monocle3.

**PAGA** (Wolf et al. 2019). Builds a graph at the cluster level — nodes are clusters, edges are weighted by the fraction of inter-cluster kNN edges in the cell graph. Does not attempt per-cell pseudotime by default; instead gives a **topology summary** showing which clusters are connected. Best as a first-pass sanity check before committing to a specific trajectory method — if PAGA shows three disconnected components, a tree-based trajectory method is wrong.

**FIGURE — Figure #1: Pseudotime — from snapshot to ordering** → `diagrams/lecture-08/01-pseudotime-intuition.svg`
*Top: a 2D UMAP showing cells sampled at a single instant, with a curved "true trajectory" drawn through them. Bottom: the same cells re-coloured by pseudotime (0 → 1 gradient), with a principal curve overlaid. Caption: one snapshot, continuous ordering.*

**FIGURE — Figure #2: Three trajectory methods side-by-side** → `diagrams/lecture-08/02-trajectory-methods.svg`
*Three small UMAPs showing the same branching dataset analyzed by (a) Monocle3 principal graph, (b) Slingshot principal curves, (c) PAGA abstracted cluster graph. Annotate topology assumptions per method.*

**EMBED — Artifact #1: Pseudotime Ordering Demo** → `artifacts/lecture-08/01-pseudotime-ordering.html`
*A synthetic branching trajectory with cells coloured by true time. Fit a principal curve; compute pseudotime. Compare inferred ordering to true time via Kendall τ.*

> **Historical pointer**: The original Monocle paper (Trapnell et al. 2014, *Nature Biotechnology*) was among the first scRNA-seq methods papers to cross 10,000 citations. The insight — that a snapshot population contains a developmental timeline — was non-obvious at the time. Monocle2 reformulated it using reverse-graph-embedding; Monocle3 (2019) scaled to 1M+ cells and added the principal-graph skeleton approach. Pseudotime is now canonical; the original Trapnell 2014 paper framed a full new sub-field in half a page.

### 1.3 Branching, roots, and evaluation (≈10 min)

Three practical issues that every trajectory analysis has to handle:

**Choosing a root.** Pseudotime orderings are direction-free by default — the method gives you a curve, but the curve could run either way. You pick the root either biologically (a cluster known to be the progenitor, e.g. HSCs for hematopoiesis) or algorithmically (the cluster with highest expression of early markers, or with the lowest RNA velocity flux — Part 2).

**Branching.** At a lineage decision point, cells commit to one of two fates. A trajectory method has to detect the branch point and assign each cell to one of the downstream branches (or neither, if the cell is at the branch point itself). Slingshot handles this well; Monocle3 handles it natively via the principal graph.

**Evaluation.** Unlike clustering, trajectory inference is rarely evaluated against a gold standard, because true developmental ordering is almost never known experimentally. When it *is* (e.g. time-course experiments where cells were harvested at known time points), evaluation uses rank correlations (Kendall τ, Spearman ρ) between inferred pseudotime and true time. Published benchmarks (Saelens et al. 2019) test all methods against 100+ simulated + real trajectories.

> **Warning box**: Pseudotime assumes the trajectory is **continuous** and that cells are sampled **densely** along it. If the underlying process has discrete state transitions (cell-cycle phases G1 / S / G2M as independent states; a sudden shock response triggering an all-or-nothing switch), fitting a continuous pseudotime will produce smooth trajectories that bear no biological meaning. Always check that the dataset actually has the shape the method expects — PAGA is the cheap sanity check.

---

## Part 2 — RNA Velocity (≈40 min)

### 2.1 Spliced vs unspliced reads (≈8 min)

In scRNA-seq, a fraction of reads come not from mature mRNA but from pre-mRNA — transcripts that still contain their introns because splicing hasn't finished yet. Cell Ranger (Lecture 7 §2.2) marks these separately: reads aligning fully within exons are **spliced**, reads overlapping intron-exon boundaries are **unspliced**. Typical scRNA-seq run: ~80% spliced, ~20% unspliced (depending on library prep and cell cycle state).

La Manno et al. 2018 noticed something useful: the ratio of unspliced to spliced transcripts per gene carries information about whether that gene's expression is increasing or decreasing at the time of sampling.

Intuition: unspliced transcripts are newer (just transcribed), spliced are older. For a gene being *turned up* right now, newly-transcribed pre-mRNA arrives faster than old mRNA is degraded — unspliced/spliced ratio is high. For a gene being *turned down*, transcription has stopped but old mature mRNA lingers — unspliced/spliced ratio is low.

This ratio lets you estimate whether each gene is being induced, at steady state, or repressed. And by summing across all genes, you get a **velocity vector** for each cell — a direction in gene space pointing toward the cell's likely next state.

**FIGURE — Figure #3: Spliced and unspliced reads** → `diagrams/lecture-08/03-spliced-unspliced.svg`
*Top: a gene model showing introns and exons. Middle: three scenarios — fully-spliced mature mRNA (all reads in exons), unspliced pre-mRNA (reads spanning intron boundaries), partially-spliced intermediate. Bottom: scRNA-seq read distribution across the gene showing both types.*

### 2.2 The dynamical model and state-space estimator (≈15 min)

The RNA velocity model (La Manno et al. 2018, Bergen et al. 2020) treats each gene as a dynamical system with two state variables per cell:

- $u$: unspliced transcript concentration (pre-mRNA)
- $s$: spliced transcript concentration (mature mRNA)

The continuous-time dynamics:

$$\frac{du}{dt} = \alpha(t) - \beta u$$

$$\frac{ds}{dt} = \beta u - \gamma s$$

where $\alpha(t)$ is the (time-varying) transcription rate, $\beta$ is the splicing rate (pre-mRNA → mature), and $\gamma$ is the degradation rate (mature → gone). Most implementations assume $\beta = 1$ (normalized).

At steady state ($du/dt = 0$, $ds/dt = 0$): $u = \alpha/\beta$ and $s = \beta u / \gamma = \alpha / \gamma$. The steady-state ratio is $u/s = \gamma / \beta$ — constant per gene.

The **velocity estimator**: for each gene, fit the steady-state ratio $\gamma$ from the cells that appear to be at steady state (upper-right and lower-left in the u-s plane — high in both or low in both). Then for each cell, compute:

$$\text{velocity}_s = \beta u - \gamma s$$

Positive velocity means $s$ is increasing (gene being induced); negative means $s$ is decreasing (gene being repressed). Aggregate across genes → a velocity vector per cell in gene space.

Project velocities into the UMAP / PCA embedding to get the familiar "arrow per cell" plot — arrows point from current state toward likely next state.

**FIGURE — Figure #4: RNA velocity dynamical model** → `diagrams/lecture-08/04-velocity-model.svg`
*Left: the two-compartment kinetic diagram (α → u → s → ∅) with rate labels. Middle: a u-s phase diagram showing steady-state line and off-diagonal cells (inducing / repressing). Right: a 2D UMAP with velocity arrows overlaid on cells.*

**EMBED — Artifact #2: RNA Velocity Vector Field** → `artifacts/lecture-08/02-rna-velocity.html`
*Synthetic cells on a 2D embedding; each cell has simulated (u, s) per gene. Fit γ per gene; compute per-cell velocities; project arrows onto the embedding. Adjustable noise and gene count sliders.*

> **EE framing — state-space estimator**: RNA velocity is a literal state-space estimator. The cell's state is $(u, s)$ per gene; the dynamics are linear ODEs with known structure but unknown rate parameters; the data is a noisy observation of $(u, s)$ at one time point. The estimator fits the steady-state relation (a linear regression on the $u$–$s$ plane) to recover $\gamma$, then evaluates the residual $\beta u - \gamma s$ as the per-cell velocity. This is exactly the form you'd see in Kalman-filter-adjacent texts: identify the state, model the dynamics, fit rate constants from steady-state behaviour, then use the fitted model to propagate state forward. The genomics-specific twist is that the "time points" are different cells rather than different instants — a population-level rather than time-series-level estimator.

### 2.3 velocyto and scVelo (≈10 min)

Two dominant implementations:

**velocyto** (La Manno et al. 2018). The original. Implements the **steady-state model** described above — fits $\gamma$ per gene from the steady-state ratio, computes velocities, projects to embedding. Conservative: assumes the population contains both steady-state and transient cells; the fit anchors on the steady-state subset. Works well for typical lineage data.

**scVelo** (Bergen et al. 2020). Extends with two additional models:

- **stochastic model**: adds second-moment information (variance across cells) to the fit, improving estimates in low-count regimes.
- **dynamical model**: fits the full time-dependent $\alpha(t)$ per gene, allowing for induction / repression kinetics (not just the steady-state assumption). Uses EM-like iteration to fit $\alpha$, $\beta$, $\gamma$, and per-cell latent times jointly.

The dynamical model is strictly more expressive but also strictly slower, and it only pays off on datasets with clear temporal structure (embryonic development, directed differentiation). For steady-ish tissue atlases, the steady-state or stochastic model is usually enough.

scVelo is the 2024 default. Integrates with Scanpy/AnnData directly; most single-cell velocity papers use scVelo's dynamical model.

> **Historical pointer**: La Manno et al.'s 2018 *Nature* paper was one of two big scRNA-seq method papers that year (the other was the original SCENIC paper). The velocity idea had existed conceptually for years — "newly-transcribed pre-mRNA marks which genes are turning on" — but nobody had the cell counts and UMI precision to use it until 10x chemistry arrived. The moment the data existed, the method followed within 12 months.

### 2.4 Caveats and failure modes (≈7 min)

Velocity has well-documented failure modes. The most important:

**Failure when assumptions break.** The steady-state model assumes that $\alpha$ is constant or piecewise-constant and that cells at the "upper tail" of the u-s plot are at steady state. If the true system has rapid, oscillatory, or bursty transcription (which many systems do), $\gamma$ is mis-estimated and the velocities point in wrong directions.

**Direction reversals in cell-cycle genes.** Cell-cycle genes oscillate. Their velocity arrows can point backward at certain cell-cycle phases, contaminating trajectories that aren't cell-cycle-related. Standard practice: regress out cell cycle before velocity analysis, or include only non-cell-cycle HVGs.

**Projection artefacts.** Velocities are computed in gene space (thousands of dimensions) but projected to UMAP (2D). The projection is via cosine similarity with nearby cells' displacements. When the UMAP is not trustworthy (small local neighborhoods distorted), the projected arrows can be misleading.

Bergen, Soldatov, et al. 2023 (*Nature Biotechnology*) documented these failure modes systematically and proposed diagnostic tests. The field's current posture: velocity is a useful hypothesis generator, but not a direct measurement of cell state evolution. Always cross-check against lineage-tracing data when you have it.

> **Discussion prompt**: You compute RNA velocity on a dataset of lung-regeneration-after-injury. The arrows point from alveolar type 2 cells toward alveolar type 1 cells, matching the known regeneration trajectory. A collaborator on the paper says: "great, this proves AT2 cells differentiate into AT1 cells during regeneration." What's the methodological objection? (Velocity shows direction of gene-expression change under the model assumption of linear induction/repression dynamics. It does not prove lineage — that would require cell-tracing experiments, e.g. genetic lineage labels. Velocity is correlative; it is consistent with the AT2 → AT1 hypothesis but does not by itself establish it.)

> **Warning box**: RNA velocity confidence drops in tissues with low cell turnover (liver, kidney) because unspliced/spliced ratios become dominated by noise when transcription is at slow steady state. Use scVelo's stochastic model (which accounts for ratio variance) for these datasets; or use direct lineage-labeling approaches (Zman-seq, metabolic labeling with 4sU, SLAM-seq) when velocity clearly isn't working.

---

## Part 3 — Batch Integration (≈45 min)

### 3.1 The batch effect problem (≈10 min)

Every scRNA-seq experiment has a **batch**: a specific capture run on a specific day, with a specific reagent lot, by a specific operator, on a specific sequencer. Batches systematically differ in ways the biology does not care about but the data reflects: small shifts in mean expression, differences in capture efficiency, variations in dropout rate.

When you want to compare or combine multiple batches — which you almost always do, because any interesting study has technical replicates, or pools samples from different patients, or integrates with public atlases — naive concatenation fails. Cells from the same biological cell type in different batches cluster *by batch*, not by biology. The UMAP of the concatenated data shows two sets of clusters that are actually just replicates.

**FIGURE — Figure #5: Batch effect before and after integration** → `diagrams/lecture-08/05-batch-integration.svg`
*Left: concatenated UMAP coloured by batch — clusters separate by batch. Right: same data post-integration — clusters separate by cell type instead, with batches mixed within each cluster.*

Formally: the cell × gene count matrix $X$ can be decomposed as

$$X = f(\text{biology}, \text{batch}, \text{technical noise})$$

Integration is the task of estimating the biology component while marginalizing out the batch component.

> **EE framing — source separation**: Batch integration is source separation. Two (or more) signal sources — the shared biological variation and the batch-specific nuisance — are mixed additively (or more generally, nonlinearly) in the observed counts. The goal is to recover the biological component. This is structurally the same problem as blind source separation in audio (ICA decomposing mixed microphones into speakers), multi-channel interference rejection in communications (separating the desired signal from co-channel interference), or anything involving "signal + nuisance, both unlabeled." The genomics-specific flavour is that the batch labels are known, which makes the problem *supervised* source separation — easier than the fully blind case.

### 3.2 Harmony — linear correction in PCA space (≈10 min)

**Harmony** (Korsunsky et al. 2019) is the simple, fast baseline. It operates entirely in PCA space and does the correction iteratively:

1. Run PCA on the concatenated data (as in Lecture 7).
2. Soft-cluster the PC embedding into K groups (fuzzy k-means).
3. For each cluster, compute the per-batch centroid. Compute a correction vector for each cell: move it toward the cross-batch cluster centroid, weighted by its cluster membership.
4. Apply corrections. Re-cluster on corrected PCs. Repeat until convergence (usually 10-50 iterations).

The result is a batch-corrected PC matrix, same shape as the input, ready for UMAP + clustering as in Lecture 7. Harmony does not touch the original counts — only the PC embedding.

Pros: fast (scales linearly with cells); simple algorithm; easy to understand; implemented in both Seurat and Scanpy. Cons: linear correction in a linear space — if the true batch effect is nonlinear, Harmony cannot fix it. Tends to over-correct: subtle biological differences between batch-specific populations can get erased along with the batch noise.

Harmony is a reasonable first pass for most datasets; it handles "cosmetic" batch differences well and runs on a laptop. For anything harder, reach for scVI.

### 3.3 scVI — variational autoencoders on counts (≈20 min)

**scVI** (Lopez et al. 2018, Gayoso et al. 2022) is a deep generative model trained per-dataset. Architecturally it's a **variational autoencoder (VAE)** with a **negative-binomial reconstruction likelihood** — a direct callback to Lecture 6's NB count model.

The generative model:

$$z_i \sim \mathcal{N}(0, I)$$

$$\rho_i = f_\theta(z_i, s_i)$$

$$x_{ij} \sim \text{NB}(\rho_{ij} \cdot l_i, \phi_j)$$

where $z_i$ is a low-dimensional latent representation of cell $i$ (typically 10-dimensional), $s_i$ is the batch label (one-hot), $\rho_i$ is the decoded mean expression profile, $l_i$ is a per-cell library-size latent (learned), and $\phi_j$ is a per-gene dispersion parameter. The encoder $q_\phi(z | x, s)$ is a neural network producing a Gaussian posterior; the decoder $f_\theta(z, s)$ is another neural network producing $\rho$.

Training maximizes the evidence lower bound (ELBO):

$$\text{ELBO} = \mathbb{E}_{q_\phi}[\log p(x | z, s)] - \text{KL}(q_\phi(z|x,s) \| p(z))$$

The first term is the reconstruction likelihood under the NB model; the second is a regularizer pulling the posterior toward the unit-Gaussian prior.

**The batch-integration mechanic**: $s_i$ is fed as input to both encoder and decoder. The encoder sees $(x, s)$ and produces $z$; the decoder sees $(z, s)$ and reconstructs $x$. At inference time, the latent $z$ is what you use for clustering and visualization. Because the decoder already has access to $s$, the latent $z$ does *not* need to encode batch identity to reconstruct $x$ — the KL regularizer explicitly penalizes any batch information that leaks into $z$. The latent $z$ is the biology-only representation.

**FIGURE — Figure #6: scVI VAE architecture** → `diagrams/lecture-08/06-scvi-architecture.svg`
*Encoder path: $x$ + batch label $s$ → latent $z$. Decoder path: $z$ + $s$ → mean expression $\rho$ → NB likelihood on raw counts. Annotated with the ELBO terms.*

**EMBED — Artifact #3: Harmony vs scVI Batch Integration** → `artifacts/lecture-08/03-batch-integration.html`
*Two synthetic batches of the same 3 cell types, with batch-specific shifts. Apply naive concatenation, Harmony, and an scVI-style encoder side-by-side. Compare how each recovers the cell-type structure.*

> **EE framing — variational EM**: scVI is variational EM on counts. The E-step is the encoder (compute posterior over latents given data); the M-step is the decoder-plus-likelihood update (maximize reconstruction given latents). This is directly analogous to the Salmon/Kallisto EM from Lecture 5 — read-to-transcript soft assignment was the E-step, transcript-abundance re-estimation was the M-step. The key architectural upgrade: Salmon's E-step had a closed-form update; scVI's latent $z$ has no closed form, so we approximate the posterior with a neural network (amortized inference) and optimize its parameters by gradient descent. The spiritual continuity is precise — both methods are iterative soft-assignment decoders, and the Lecture 5 thread on iterative soft-decision decoding carries directly into modern deep-learning genomics.

> **Historical pointer**: The VAE architecture (Kingma & Welling 2013) came from the deep-learning side of the street; applying it to single-cell counts with a domain-appropriate likelihood (NB rather than Gaussian) was the Lopez et al. 2018 contribution. scVI, totalVI (multimodal), scANVI (semi-supervised), and the whole scvi-tools family now covers most single-cell modelling needs. It's a rare case where an EE-trained student will find the underlying architecture entirely familiar from general ML courses — the domain adaptation is small, the architecture itself is standard VAE.

### 3.4 Evaluating integration quality (≈5 min)

Integration is a tradeoff. Over-correct, and you erase biological differences (different responder cell types from the same cell class get fused). Under-correct, and batches still separate. Three metrics evaluate the tradeoff:

**LISI — local inverse Simpson index** (Korsunsky 2019). For each cell, count how many batches are represented in its k-nearest neighbours. Average across cells. High LISI means batches are well-mixed locally; low LISI means batches still separate. Evaluated as `iLISI` (batch LISI — want it high) vs `cLISI` (cell-type LISI — want it low).

**kBET — k-nearest-neighbour batch-effect test** (Büttner et al. 2019). For each cell's k-nearest neighbourhood, test whether the batch composition matches the global batch composition (chi-squared test). Reports the fraction of neighbourhoods that "pass" the test.

**ASW — average silhouette width**. Standard clustering metric, but computed separately using batch labels (want it low: batches don't separate) and cell-type labels (want it high: cell types do separate).

The scib benchmark (Luecken et al. 2022) combines all three plus others into a single integration score and ranks methods on diverse datasets. The 2022 benchmark verdict: no single method wins across all scenarios; scVI and scANVI lead on large, complex atlases; Harmony leads on small, simple batches.

---

## Part 4 — Multi-modal Single-Cell (≈30 min)

### 4.1 CITE-seq — RNA plus surface protein (≈10 min)

**CITE-seq** (Stoeckius et al. 2017) measures RNA and surface-protein abundance on the same cells. Before the cells enter the 10x chip, they are incubated with a panel of antibodies, each carrying a DNA tag with its own unique barcode (the antibody-derived tag, ADT). When the cell is captured and lysed inside its droplet, the antibody tags release and are reverse-transcribed alongside the mRNA. A single sequencing run produces:

- The standard mRNA count matrix (Lecture 7).
- A per-cell protein matrix: ~100–300 surface proteins, measured as ADT counts.

CITE-seq gives you the same cells measured two ways. Protein markers like CD4, CD8, CD19 — which are the canonical cell-type indicators in flow cytometry — can now be measured directly per cell alongside transcript expression.

Why it matters: protein expression does not track transcript expression perfectly. CD4 mRNA and CD4 protein are moderately correlated but not identical — a cell with high CD4 mRNA may not yet have CD4 protein on its surface. For cell-type annotation, protein is often the gold standard; CITE-seq lets you validate transcript-based annotations directly.

**FIGURE — Figure #7: CITE-seq schematic and joint matrices** → `diagrams/lecture-08/07-cite-seq.svg`
*Top: a cell with antibody tags bound to surface proteins; the tags are RT'd and sequenced together with mRNA. Bottom: two count matrices — the RNA matrix (cells × genes) and the protein matrix (cells × antibodies) — sharing the same cell barcodes.*

### 4.2 scATAC-seq — open chromatin per cell (≈10 min)

**scATAC-seq** (Buenrostro et al. 2015, Cusanovich et al. 2015) measures chromatin accessibility per cell. The assay uses the Tn5 transposase — introduced in bulk ATAC-seq (Lecture 9) — adapted for single-cell chemistry. Tn5 inserts sequencing adapters preferentially into open (nucleosome-free) regions of DNA; reads pile up at regulatory elements and active gene promoters.

The output format is different from scRNA-seq. Instead of a cells × genes count matrix, scATAC-seq produces a cells × peaks matrix where each peak is a region of the genome (typically a 500 bp window). Entries are binary or low count — a cell either has a Tn5 insertion in a given peak or it doesn't. Peaks are extremely sparse per cell (~10,000 fragments per cell distributed over ~200,000 peaks).

The sparsity is more severe than scRNA-seq. Most peaks are 0 in most cells; dimensionality reduction typically uses **LSI** (latent semantic indexing, a TF-IDF-weighted SVD borrowed from text retrieval) rather than PCA. UMAP + Leiden follows as in Lecture 7.

Multi-omic chemistry (10x Genomics Multiome) measures RNA and ATAC on the same nucleus in the same droplet. Each cell produces a gene expression vector and a chromatin accessibility vector — the two modalities are directly linked.

### 4.3 Joint embeddings — WNN, totalVI, MOFA+ (≈10 min)

Given two (or more) modalities per cell, how do you produce a single embedding that uses both? Three approaches:

**Weighted Nearest Neighbours (WNN)** (Hao et al. 2021, Seurat v4). Run PCA on each modality separately → get two kNN graphs. Weight each cell's two graphs by the relative signal-to-noise of each modality at that cell (compute this locally). Fuse into a single WNN graph; cluster and UMAP on the fused graph. Simple, no per-dataset training, works in Seurat directly.

**totalVI** (Gayoso et al. 2021). Extends scVI to handle RNA + protein jointly. Two decoders (one NB for RNA, one for protein with a specialized likelihood) share a single latent $z$. The latent encodes both modalities simultaneously. Trained by VAE with joint ELBO. Handles batch integration natively.

**MOFA+** (Argelaguet et al. 2020). Factor analysis on multi-modal data. Each cell is represented as a linear combination of factors; factors can be modality-specific or shared. Unsupervised discovery of joint structure. Fast, interpretable (factors map back to genes/proteins), but less powerful than VAE approaches for high-dimensional data.

**EMBED — Artifact #4: CITE-seq Joint Embedding** → `artifacts/lecture-08/04-cite-seq-joint.html`
*Synthetic CITE-seq data: 3 cell types, each with distinct RNA and protein profiles but with one cell type where protein and RNA disagree (post-transcriptional regulation). Compute RNA-only UMAP, protein-only UMAP, and a WNN joint embedding. See which modality catches the discrepancy.*

> **EE framing — joint embedding**: WNN is sensor fusion. Two sensors (RNA and protein) measure the same target (cell state) with different noise characteristics; the fused estimate weights each sensor by its local reliability. This is precisely the structure of Kalman-filter sensor fusion, multi-sensor radar tracking, and any "combine two measurement channels" problem. scVI / totalVI do the same thing in a learned latent space rather than via weighted kNN — the architectural analog is an autoencoder with two input heads sharing a bottleneck, trained to reconstruct both modalities from the shared latent.

> **Intuition box**: For most cells, RNA and protein agree — the cell expressing CD19 mRNA also shows CD19 protein on its surface. The interesting cells are the ones where they *disagree*: a transitioning cell that has upregulated CD19 mRNA but hasn't yet accumulated protein, or a cell that has lost mRNA but retains surface protein from past expression. Joint embeddings let you find these cells as outliers from the main "agreement" axis.

---

## Part 5 — Spatial Transcriptomics (≈40 min)

### 5.1 Platforms and resolution (≈10 min)

Spatial transcriptomics measures gene expression *with the cells' positions preserved*. Three major platform families in 2024, trading off panel size, resolution, and throughput:

**Visium (10x Genomics)** — sequencing-based. A tissue section is laid on a capture slide covered in ~5000 spots, each ~55 μm in diameter and ~100 μm centre-to-centre. Each spot captures mRNA from 1–20 cells (depending on tissue density). All genes are measured (~20,000), but cells are *not* individually resolved — each spot is a **mixed pixel** of the cells under it.

**MERFISH / seqFISH / HybISS** — imaging-based with multiplexed FISH. Pre-designed panels of 100–500 targeted genes. Multiple rounds of fluorescent-probe hybridisation with combinatorial labelling yield a per-cell gene-expression vector at sub-cellular resolution. Throughput is limited by microscopy speed; a typical experiment covers ~500,000 cells per slide.

**Xenium (10x Genomics)** — imaging-based, commercial. A preset panel of ~400 genes per run; sub-cellular resolution; fully automated. A direct competitor to MERFISH for commercial labs.

Key tradeoff: sequencing-based (Visium) measures every gene but loses cellular resolution; imaging-based (MERFISH, Xenium) achieves cellular resolution but measures only a panel. The right platform depends on whether you need an unbiased survey or a targeted interrogation.

**FIGURE — Figure #8: Spatial platforms — Visium vs MERFISH vs Xenium** → `diagrams/lecture-08/08-spatial-platforms.svg`
*Three small schematics: (a) Visium slide with capture spots and mixed cells under each, (b) MERFISH cells imaged with fluorescent probes binding mRNAs, (c) Xenium — similar to MERFISH but commercial instrument. Comparison table underneath: resolution, panel size, throughput.*

### 5.2 The mixed-pixel problem and deconvolution (≈15 min)

Visium's core challenge: each spot contains multiple cells, so the spot's measured counts are a weighted sum of contributions from each cell type present. If you want per-cell-type expression, you must **deconvolve** the mixed signal — estimate the per-cell-type contribution at each spot.

Formally: each spot $i$ has observed counts $x_i$. Each cell type $k$ has a reference profile $\mu_k$ (learned from scRNA-seq). The spot's signal is a mixture:

$$x_i \approx \sum_k w_{ik} \mu_k$$

where $w_{ik}$ is the proportion of spot $i$'s content attributable to cell type $k$. Solve for $\mathbf{w}_i$ at each spot (constrained: non-negative, sums to 1).

Two dominant deconvolution methods:

- **RCTD** (Cable et al. 2021). Maximum-likelihood fit under a Poisson count model. Uses a scRNA-seq atlas of the same tissue as the cell-type reference. Two modes: "singlet" mode assumes one cell per spot (appropriate for high-resolution platforms); "doublet" mode allows two cell types per spot (for Visium).
- **cell2location** (Kleshchevnikov et al. 2022). Bayesian model with spatial priors. Produces posterior distributions over cell-type proportions with explicit uncertainty. Slower but more principled for downstream spatial analysis.

**FIGURE — Figure #9: Visium spot deconvolution** → `diagrams/lecture-08/09-spot-deconvolution.svg`
*Left: a Visium spot containing 5 cells of 3 different types overlaid on a tissue image. Right: the deconvolution output — a bar chart of estimated proportions for 5 cell types at that spot. Cells under the spot are coloured; bar heights match the true proportions.*

**EMBED — Artifact #5: Visium Spot Deconvolution** → `artifacts/lecture-08/05-spot-deconvolution.html`
*Simulated Visium spots with 3 cell types. Given reference scRNA-seq profiles, fit per-spot cell-type proportions via constrained least squares. Compare inferred proportions to ground truth.*

> **EE framing — mixed-pixel unmixing**: Spatial deconvolution is the remote-sensing **mixed-pixel problem**. A multispectral satellite pixel observes a weighted sum of contributions from multiple land-cover classes beneath it (soil, vegetation, water, urban); the unmixing algorithm solves for per-class fractional abundances under spectral endmember constraints. The spatial-transcriptomics version replaces "land cover" with "cell type," "spectral endmember" with "cell-type expression profile," and "pixel" with "Visium spot." The mathematical structure is identical — constrained linear unmixing $x \approx M w$ with $w \geq 0, \sum w = 1$ — and the same algorithmic toolbox (non-negative least squares, spectral unmixing under Bayesian priors) transfers directly. Gene-expression "endmembers" come from the reference scRNA-seq atlas rather than from lab spectral libraries; everything else is the same.

### 5.3 Spatial niches and neighbourhood analysis (≈10 min)

Once cells (or spots) have both expression profiles and spatial coordinates, new analyses become possible:

**Niche detection.** Group cells or spots by their spatial neighbourhood composition — find regions where a specific combination of cell types co-occur. E.g., tumor border niches (tumor + immune cells + fibroblasts), specific crypt-villus neighbourhoods in gut epithelium.

**Spatial differential expression.** For a given cell type, test whether its expression profile depends on its spatial context. A macrophage near a tumor may express different genes than a macrophage in healthy tissue — both labeled "macrophage" in scRNA-seq, but spatially distinguished.

**Gradients and axes.** Many tissues have axes along which gene expression varies continuously — cortical layers, liver zonation (pericentral vs periportal), intestinal crypt-villus. Spatial analysis extracts these axes directly from the data rather than requiring them to be annotated a priori.

Tools: **Squidpy** (Palla et al. 2022) is the Scanpy-integrated spatial analysis toolkit; **Giotto** (Dries et al. 2021) is the R equivalent; **BANKSY** (Singhal et al. 2024) does spatially-aware clustering.

### 5.4 What's still hard about spatial (≈5 min)

Spatial transcriptomics is new enough to have limitations worth naming:

- **Resolution-vs-coverage tradeoff is real.** No current platform gives sub-cellular resolution *and* whole-transcriptome coverage *and* high throughput. Study design forces a choice.
- **3D is nascent.** Most spatial experiments are 2D sections. Reconstructing 3D structure from serial sections is open problem; native 3D (light-sheet microscopy combined with barcoding) exists but is experimental.
- **Single-cell-plus-spatial requires either high-resolution imaging or probabilistic assignment.** Visium alone doesn't give single cells; the integration step matters.
- **Batch effects across slides.** Analogous to scRNA-seq batch effects (Part 3) but with additional spatial autocorrelation complications. Active area of method development.

> **Warning box**: A "spatial" analysis is only as good as its spatial registration. Slides from different experiments, different days, or different microscopes need alignment — and the coordinate systems don't align trivially. Always check that your "spatially varying gene" is not just a slide-to-slide offset in probe efficiency. The pre-registration pipeline (SpaceRanger for Visium; Xenium Onboard Analysis; MERFISH's registration software) is critical to trust.

---

## Part 6 — Cell-Cell Communication (≈15 min)

### 6.1 Ligand-receptor inference (≈8 min)

Cells communicate via secreted ligands (cytokines, growth factors, chemokines) that bind receptors on other cells. Given a scRNA-seq dataset with labeled cell types, we can ask: which cell types are *sending* signals to which others, and via what ligand-receptor pairs?

The inference is indirect — scRNA-seq measures transcripts, not ligand secretion or receptor binding. The usual proxy: if cell type A highly expresses ligand $L$ and cell type B highly expresses the cognate receptor $R$, infer that A can signal to B via the $(L, R)$ axis. Score each cell-type-pair × ligand-receptor-pair by the product (or geometric mean) of ligand and receptor expression.

Two dominant tools:

**CellChat** (Jin et al. 2021). Uses a curated database of ~2000 ligand-receptor pairs (CellChatDB) including known complexes (e.g., IL2 binds a 3-subunit receptor; the score uses all three). Permutation testing for significance. Outputs ranked pairs and network visualisations.

**NicheNet** (Browaeys et al. 2020). Goes further: models ligand → receptor → downstream target gene chains. A ligand is inferred to act if its downstream target genes are differentially expressed in the receiver cell population. This tests for *functional* signalling, not just co-expression.

**FIGURE — Figure #10: Ligand-receptor communication network** → `diagrams/lecture-08/10-ligand-receptor.svg`
*A directed graph: nodes are cell types, edges are ligand-receptor signalling axes. Edge thickness = total communication score; edge colour = ligand family. Example: T cells → macrophages via IFN-γ; macrophages → endothelial via VEGF.*

**EMBED — Artifact #6: Ligand-Receptor Network Browser** → `artifacts/lecture-08/06-ligand-receptor.html`
*A synthetic 5-cell-type dataset with known ligand and receptor expression. Compute per-cell-type-pair communication scores across a toy L-R database. Visualise the network; filter by ligand family; threshold edges to declutter.*

### 6.2 Caveats and what communication scores mean (≈7 min)

Ligand-receptor inference is one of the most interpretation-heavy steps in single-cell analysis. The scores look quantitative but they rest on several not-always-true assumptions:

**Expression ≠ secretion.** A cell expressing a cytokine transcript may or may not actually secrete protein. Many cytokines are post-translationally regulated (requires cleavage, requires stimulus, requires correct secretion machinery). High transcript expression in scRNA-seq does not guarantee high protein output.

**Receptor expression ≠ active receptor.** Similarly, cells can express a receptor transcript without surface expression (stored in vesicles, not trafficked), or with surface expression but inactive (desensitized, internalized).

**Spatial information is usually ignored.** Ligands diffuse with short range (tens of microns for most cytokines); two cell types that co-express a ligand-receptor pair are not communicating if they're in different tissue compartments. Pure-scRNA-seq-based tools (CellChat, NicheNet) ignore space; spatial-aware tools (COMMOT, Giotto's cell-cell-interaction module) incorporate tissue coordinates when available.

**The statistics are shallow.** Permutation tests in CellChat test "is this pair's score higher than chance?" — but the null distribution does not account for many biological sources of correlation. Published communication networks often have many false positives.

> **Discussion prompt**: You run CellChat on a tumor dataset and get 300 significant ligand-receptor pairs across 15 cell types. A wet-lab collaborator asks which one to follow up on experimentally. What filters do you apply before picking one? (Possibilities include: restrict to pairs where the ligand is secreted/extracellular — not all hits are; require downstream target-gene corroboration via NicheNet; check for known biology — a pair involving IL-2 or IFN-γ between immune cell types is much better-supported than a novel axis; check spatial co-location if possible; prioritize druggable receptors for translational value. The point of the exercise: communication inference produces hypotheses, not answers.)

> **Intuition box**: Communication scores are like a dating app's compatibility score — computed from the data both people put into their profile. They tell you whether two people *might* be a good match in principle. Whether they actually ever meet, speak, or form a connection requires evidence the app cannot see. The same discount applies: scRNA-seq ligand-receptor inference tells you which communication axes are plausible given expression, not which are occurring in vivo.

---

## Wrap-up (≈10 min)

### What you should take away

- **Pseudotime turns a snapshot into an ordering.** The assumption is that a differentiating population is asynchronous — cells exist at every trajectory position. The method is principal-curve / principal-graph fitting to the embedding.
- **RNA velocity is a state-space estimator.** Spliced and unspliced transcript ratios, under a linear-ODE kinetic model, give a per-cell velocity vector in gene space. Project to 2D for the familiar arrow plot. Treat as hypothesis-generating, not lineage-proving.
- **Batch integration is supervised source separation.** Harmony is the fast linear baseline in PCA space; scVI is a VAE with NB likelihood that learns a biology-only latent. scVI's variational EM is the modern incarnation of the Salmon EM thread — iterative soft-assignment through a learned encoder.
- **Multi-modal single-cell is sensor fusion.** CITE-seq adds protein; scATAC adds chromatin accessibility; totalVI / WNN / MOFA+ combine modalities into joint embeddings.
- **Spatial transcriptomics preserves position.** Visium measures all genes with mixed-pixel spots; MERFISH/Xenium measure gene panels at sub-cellular resolution. Spot deconvolution is mixed-pixel unmixing — the same math as remote-sensing land-cover analysis.
- **Ligand-receptor inference is hypothesis generation.** Expression-based communication scores have known blind spots (secretion, activity, diffusion range). Treat outputs as a ranked hypothesis list to triage with wet-lab or spatial data.

### Next lecture

Epigenomics. Where proteins bind DNA (ChIP-seq), where chromatin is open (ATAC-seq at bulk scale). Peak calling as CFAR-style detection. Motif analysis as matched filtering on a 4-channel signal. A forward pointer to Enformer — the deep-learning model that predicts regulatory landscapes end-to-end from DNA sequence.

### Homework

1. Run scVelo's dynamical model on the pancreas endocrinogenesis tutorial dataset. Report: which cell types show the clearest velocity arrows? Which have noisy or inconsistent velocity? Discuss why.
2. Take the 10x PBMC 3k dataset from Lecture 7 and the 10x PBMC 10k dataset (separate runs). Concatenate naively; compute UMAP. Then apply Harmony; compute UMAP again. Report LISI scores before and after. Now apply scVI; compare all three.
3. Use the Seurat CITE-seq tutorial (cord blood mononuclear cells, 10 proteins). Compute RNA-only UMAP, protein-only UMAP, and WNN joint UMAP. Find one cell population where the WNN embedding reveals structure that RNA-only missed.
4. Download a Visium brain section dataset. Apply RCTD deconvolution against a reference scRNA-seq brain atlas. Map the inferred cortical-layer cell-type proportions to the spatial coordinates — do you see the expected cortical-layer stratification?
5. Run CellChat on any tutorial PBMC or tumor dataset. Pick the top 10 ligand-receptor pairs by combined score. For each: is the ligand known to be secreted? Is the receptor known to be on the cell surface? How many survive this filter?

### Recommended reading

- Trapnell, C., Cacchiarelli, D., Grimsby, J., et al. (2014). The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nature Biotechnology* 32, 381–386. (The original Monocle paper.)
- Saelens, W., Cannoodt, R., Todorov, H., &amp; Saeys, Y. (2019). A comparison of single-cell trajectory inference methods. *Nature Biotechnology* 37, 547–554.
- La Manno, G., Soldatov, R., Zeisel, A., et al. (2018). RNA velocity of single cells. *Nature* 560, 494–498.
- Bergen, V., Lange, M., Peidli, S., Wolf, F. A., &amp; Theis, F. J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology* 38, 1408–1414. (The scVelo paper.)
- Bergen, V., Soldatov, R. A., Kharchenko, P. V., &amp; Theis, F. J. (2023). RNA velocity — current challenges and future perspectives. *Molecular Systems Biology* 19, e11799.
- Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., &amp; Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods* 15, 1053–1058. (The scVI paper.)
- Korsunsky, I., Millard, N., Fan, J., et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* 16, 1289–1296.
- Luecken, M. D., Büttner, M., Chaichoompu, K., et al. (2022). Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* 19, 41–50. (The scib benchmark paper.)
- Stoeckius, M., Hafemeister, C., Stephenson, W., et al. (2017). Simultaneous epitope and transcriptome measurement in single cells. *Nature Methods* 14, 865–868. (The CITE-seq paper.)
- Hao, Y., Hao, S., Andersen-Nissen, E., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell* 184, 3573–3587. (The Seurat v4 / WNN paper.)
- Cable, D. M., Murray, E., Zou, L. S., et al. (2021). Robust decomposition of cell type mixtures in spatial transcriptomics. *Nature Biotechnology* 40, 517–526. (The RCTD paper.)
- Jin, S., Guerrero-Juarez, C. F., Zhang, L., et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications* 12, 1088.
- scvi-tools documentation: <https://docs.scvi-tools.org/>
- scVelo tutorials: <https://scvelo.readthedocs.io/>
- Squidpy tutorials: <https://squidpy.readthedocs.io/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Pseudotime and Trajectory Inference        | 30&nbsp;min | 0:30 |
| Part 2 — RNA Velocity                                | 40&nbsp;min | 1:10 |
| Part 3 — Batch Integration                           | 45&nbsp;min | 1:55 |
| Part 4 — Multi-modal Single-Cell                     | 30&nbsp;min | 2:25 |
| Part 5 — Spatial Transcriptomics                     | 40&nbsp;min | 3:05 |
| Part 6 — Cell-Cell Communication                     | 15&nbsp;min | 3:20 |
| Wrap-up                                               | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
