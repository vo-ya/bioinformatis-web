# Lecture 16 — ML in Genomics: Architectures, Pitfalls, Frontiers

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-16.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Map the major genomic prediction problems to their dominant ML architectures: pileups → CNN, long-range regulation → dilated convolution + transformer, count matrices → VAE, molecular graphs → GNN, protein 3D → equivariant transformer.
2. Explain *why* each architecture matches its problem — translate each architecture's inductive bias into the structural property of the data it exploits.
3. Describe where labelled data in genomics comes from, why it's chronically limited, and how self-supervised pretraining compensates.
4. Diagnose data-leakage patterns in genomics: sequence homology, family relatedness, cross-study batch effects. Design train/test splits (chromosome-wise, species-wise, family-wise) that avoid them.
5. Walk through DNA language models (DNABERT, Nucleotide Transformer, HyenaDNA, Evo) and give an honest assessment of what pretraining-on-DNA does and does not buy.
6. Assess cell foundation models (scGPT, Geneformer, scFoundation) — what they claim, what the evidence actually shows, and what "zero-shot" means at the single-cell scale.
7. Predict what kinds of genomics ML advances are likely to make the biggest difference in the next 3 years and why.

---

## Part 1 — Architectural Pattern Survey (≈45 min)

### 1.1 Why now, a synthesis lecture (≈4 min)

You've already met the main architectures in this course:

- **L4 — DeepVariant**: CNN classifier on read pileups.
- **L6 — DESeq2 / Limma**: GLM-based differential-expression testing (pre-deep-learning, included for contrast).
- **L8 — scVI / totalVI**: VAE for single-cell count matrices.
- **L9 — Enformer / Borzoi**: dilated convolution + transformer for long-range regulatory sequence.
- **L13 — SuSiE / LDpred**: Bayesian sparse regression for fine-mapping + PRS.
- **L15 — AlphaFold2**: equivariant transformer + iterative refinement for protein structure.

Individually, each lecture treats its method as "the right approach for this problem". This lecture pulls the pattern out: for each problem, the architecture's **inductive bias** matches the data's **structural property**. Architecture-problem matching is the core design skill in genomics ML.

### 1.2 Pileups → Convolutional Neural Network (≈8 min)

**Problem** (Lecture 4): variant calling. Given a pileup of reads at a candidate position, classify the genotype (ref/ref, ref/alt, alt/alt, no-variant).

**Input representation**: a small 2D image. Rows = reads; columns = positions around the candidate; channels = (base, base-quality, strand, mapping-quality, variant-support, ...). Typical size: ~100 rows × ~200 columns × ~6 channels.

**Architecture**: 2D CNN. Convolutional layers extract local patterns (a few bases across a few reads); pooling reduces to a representation that's increasingly context-aware; dense layers produce the final classification.

**Why CNN works**:

- **Translation invariance**: the exact column of the candidate doesn't matter; what matters is the local read support around it. Convolutional weights are shared across positions, so the model is invariant to "where in the window" by construction.
- **Local feature extraction**: mismatches, deletions, and strand-bias patterns are local features in the 2D image. Small convolutional kernels (3 × 3) detect them efficiently.
- **Compositional hierarchy**: early layers detect individual base mismatches; later layers detect concerted patterns (many reads with the same mismatch = strong variant evidence).

DeepVariant (Poplin et al. 2018) demonstrated that a CNN on pileup images could match or exceed the hand-engineered HaplotypeCaller at variant-call precision/recall. The inductive bias of the CNN — translation-invariant local feature extraction — matches the structure of pileup evidence exactly.

**FIGURE — Figure #1: DeepVariant pileup-CNN** → `diagrams/lecture-16/01-deepvariant-cnn.svg`
*Left: pileup image with reads as rows, positions as columns, bases colour-coded. Middle: CNN architecture — stack of Conv → ReLU → MaxPool layers; final dense → softmax over 4 genotype classes. Right: output probability vector. Annotation at bottom: "translation-invariant local-feature extraction matches the pileup's structure: mismatches are local patterns; position-within-window doesn't matter".*

### 1.3 Long-range sequence regulation → Dilated Conv + Transformer (≈9 min)

**Problem** (Lecture 9): predict regulatory output (TF binding, chromatin accessibility, gene expression) from 100 kb of surrounding genomic DNA.

**Input representation**: a one-hot encoded DNA sequence, ~100,000 × 4.

**Architectures**:

- **Stage 1 (local)**: stacked dilated convolutions. Each layer's receptive field grows exponentially with dilation rate. After 10 layers with dilations 1, 2, 4, ..., 512, the receptive field is ~2 kb per position. Extracts local motif-scale features.
- **Stage 2 (long-range)**: transformer layers with multi-head self-attention operating on the dilated-conv output. Attention heads learn to route information from distal regulatory elements (enhancers ~100 kb away) to the gene body.

The combination (dilated conv → transformer) is the **Enformer** / **Borzoi** architecture.

**Why this works**:

- **Hierarchy matches scale**: local features (motifs ~10 bp) → intermediate features (enhancer blocks ~1 kb) → long-range interactions (enhancer-to-gene ~100 kb). Dilated conv gets efficiency for the first two scales; attention gets the third.
- **Attention is the right tool for sparse long-range interactions**: only a handful of distal enhancers actually regulate any given gene. Attention's per-query routing (high weight to a few keys, low weight to everything else) captures this sparsity naturally.
- **Receptive-field management**: full self-attention at 100 kb resolution would be $O(L^2) = O(10^{10})$ — infeasible. The dilated-conv stage down-samples first (factor ~100), bringing the sequence length into transformer-feasible range.

Alternative: **HyenaDNA** (Nguyen et al. 2023) replaces attention with **long convolutions** parameterised via implicit operators; scales to >1M bp contexts with linear cost. Same overall idea — extract long-range dependencies — different computational machinery.

**FIGURE — Figure #2: Enformer architecture** → `diagrams/lecture-16/02-enformer.svg`
*Left: 100 kb DNA sequence one-hot encoded. Middle: stack of dilated convolutions (dilations 1, 2, 4, ..., 512) expanding the receptive field. Middle-right: multi-head self-attention over the down-sampled sequence (each position attends to all others). Right: output heads predicting chromatin accessibility, CAGE (gene expression), ChIP-seq per cell type. Annotation: "dilated conv = efficient local-to-mid scale; attention = sparse long-range routing".*

### 1.4 Count matrices → Variational Autoencoder (≈8 min)

**Problem** (Lecture 8): model single-cell count data. Observations: cells × genes count matrix with high sparsity, batch effects, and technical dropout.

**Input representation**: per-cell count vector over ~20k genes.

**Architecture**: Variational Autoencoder (VAE).

- **Encoder**: cell's count vector → low-dimensional latent $\mathbf{z}$ (~10–30 dims) via a neural network.
- **Decoder**: $\mathbf{z}$ + batch covariates → parameters of a count likelihood (typically negative binomial).
- **Loss**: negative ELBO = reconstruction (NB log-likelihood of counts given parameters) + KL regularisation (match $q(\mathbf{z}|\text{cell})$ to prior $p(\mathbf{z}) = \mathcal{N}(0, I)$).

**scVI** (Lopez et al. 2018) is the canonical implementation. **totalVI** extends to multi-modal (RNA + protein) via a shared latent. **MultiVI** extends to RNA + ATAC.

**Why VAE works**:

- **Dimensionality reduction via learned prior**: a generative model with a Gaussian latent naturally produces a compact representation. PCA does this too but linearly; VAE's nonlinear encoder captures cell-state manifolds.
- **Right likelihood** matters. Gene counts are over-dispersed; using Gaussian reconstruction (as vanilla VAE would) grossly misfits. Negative binomial (or zero-inflated NB) matches the empirical distribution.
- **Batch correction** via conditioning: decoder takes batch as input, so the latent $\mathbf{z}$ learns a batch-invariant representation. Analogous to conditioning on confounders in a regression.
- **Probabilistic framework**: the ELBO gives a likelihood; uncertainty in the latent is calibrated (crudely); downstream Bayesian analyses remain valid.

> **EE framing — VAE with negative binomial likelihood as a statistically calibrated auto-encoder**: A vanilla auto-encoder minimises reconstruction MSE, which implicitly assumes Gaussian noise. Count data don't follow Gaussian distributions — they follow NB, Poisson, or ZINB. Switching the reconstruction likelihood from Gaussian MSE to NB log-likelihood is the right domain adaptation. The same pattern appears in signal processing when the measurement noise model is wrong (assume Gaussian when it's Poisson → estimation bias); using the correct likelihood fixes it. scVI's core trick is **not** "we used deep learning on sc-RNA"; it's "we matched the noise model to the data."

### 1.5 Molecular graphs → Graph Neural Network (≈6 min)

**Problem**: predict molecular properties (solubility, binding affinity, toxicity) from molecular structure. Often cheminformatics-adjacent; in genomics, relevant for protein-ligand interaction, metabolomics.

**Input representation**: a graph. Nodes = atoms (with features: element, hybridisation, formal charge); edges = bonds (with features: bond type, aromaticity). Typical size: 10–100 atoms.

**Architecture**: Graph Neural Network (message-passing).

- At each layer, each node aggregates messages from its neighbours (summing or attending).
- Node features update; after several layers, each node's representation encodes a neighbourhood of increasing radius.
- Final graph-level readout (sum/mean/attention pool) produces a single vector per molecule.

Notable GNN variants: **GCN** (Kipf &amp; Welling 2017), **GraphSAGE** (Hamilton et al. 2017), **GIN** (Xu et al. 2019), **message-passing neural networks** (Gilmer et al. 2017).

**Why GNN works**:

- **Permutation invariance**: atoms in a molecule have no canonical order; the GNN is invariant to reordering the nodes.
- **Locality**: most chemical properties depend on local neighbourhoods (functional groups, rings). Message-passing in a few layers captures a ~3-bond radius, matching what chemists intuit.
- **Explicit graph structure**: edges encode bonds; atoms that are far apart in sequence but close in graph distance are handled correctly.

### 1.6 Protein 3D → Equivariant Transformer (≈6 min)

**Problem** (Lecture 15): predict 3D protein structure from sequence.

**Input representation**: MSA + pair representation (L × L).

**Architecture**: Evoformer (axial attention + triangle updates) + structure module (Invariant Point Attention, SE(3)-equivariant).

**Why this works**:

- **Axial attention on MSA**: 2D data (sequences × residues) — full 2D attention is prohibitive. Axial attention is quadratic in one axis at a time.
- **Triangle updates**: respect geometric consistency of pair representations — if (i,j) and (j,k) are close, (i,k) should be constrained. Inductive bias from physics.
- **SE(3) equivariance in structure module**: rotate input, output rotates the same. No augmentation by rotation needed during training; much more sample-efficient.

**FIGURE — Figure #3: Architecture-to-problem map** → `diagrams/lecture-16/03-architecture-map.svg`
*A five-column table / grid. Columns: Problem / Input representation / Architecture / Key inductive bias / Example model. Rows: variant calling (pileup CNN / translation invariance / DeepVariant); long-range regulation (dilated conv + transformer / multi-scale + sparse long-range / Enformer); count matrix (VAE with NB / NB likelihood + latent prior / scVI); molecular graph (GNN / permutation invariance + locality / GIN); protein structure (equivariant transformer / axial attention + SE(3) / AlphaFold2). Each row colour-tagged.*

### 1.7 Common threads (≈4 min)

Looking across all five:

- Each architecture **commits** to a structural prior. That prior must match the data.
- "Commit to a prior" is the essence of **inductive bias**. Stronger prior → more sample-efficient on matching data, but worse on mismatched data.
- Architectures don't need to be exotic to work. Most genomics ML is CNN, transformer, VAE, GNN, or combinations. Innovation comes from matching the right variant to the right problem.

**EMBED — Artifact #1: Architecture-to-Problem Matcher** → `artifacts/lecture-16/01-architecture-matcher.html`
*A matching game. Presents a genomic problem (e.g. "predict enhancer activity from 100 kb sequence"). Student picks from CNN / dilated conv + transformer / VAE / GNN / equivariant transformer. Immediate feedback explains why each choice does or doesn't fit. Target aha: the choice is usually obvious once you identify the data's structural property.*

---

## Part 2 — Inductive Biases: Why Does Each Work? (≈35 min)

### 2.1 What "inductive bias" really means (≈6 min)

An **inductive bias** is any assumption an algorithm makes that isn't justified by the training data alone. Every learning algorithm has them, whether explicit or hidden.

- **Linear regression**: assumes the response is a linear function of the predictors.
- **k-NN**: assumes nearby points in feature space have similar labels.
- **Random forest**: assumes response decomposes into axis-aligned piecewise-constant regions.
- **CNN**: assumes features are local + translation-invariant.
- **Transformer (with positional embeddings)**: assumes contextual interactions matter, but is relatively *weak* on locality or translation.
- **Equivariant network**: assumes invariance to a specific group (SO(3), SE(3), permutation).

Without an inductive bias, learning needs exponentially more data to generalise (the **no-free-lunch theorem** makes this precise).

The genomics ML design problem is: **pick an architecture whose inductive bias matches your problem's structure.** Mismatch is the #1 cause of genomics ML failures.

> **Intuition box**: Inductive bias is what the architect brings to the table — a prior about what "should work." If you know the feature you want is local and position-independent (a motif), baking translation invariance into your network (CNN) is worth many thousands of training examples. If the feature is inherently rotation-invariant (3D protein structure), baking SO(3) equivariance into your network buys you sample-efficient training. If the data has no useful prior structure, a fully general transformer + mountains of data is your remaining option. Picking the right inductive bias is why genomics ML isn't just "more parameters, more data" — it's choosing the prior that fits the problem.

### 2.2 CNN vs MLP vs Transformer on pileup-like data (≈8 min)

**Scenario**: you're classifying pileup images.

- **MLP** (fully connected): treats each pixel independently. No locality, no translation invariance, no weight sharing. Must learn each feature separately for each position.
- **CNN**: local + translation-invariant + weight-sharing. Dramatically more sample-efficient.
- **Transformer**: attends globally, no locality bias. More data needed to learn that only local patterns matter.

On the DeepVariant task, at comparable training budget:

- CNN: ~99.5% F1 (the DeepVariant-reported number).
- Transformer with absolute positional embeddings: ~98.0% F1 with 10× more training data (numbers illustrative; the pattern is real).
- Transformer *with 2D positional embeddings + local-attention bias*: ~99.3% F1 — **restoring the CNN's inductive bias by other means**, the gap narrows.
- MLP: ~95% F1 — underfits systematic position/translation structure.

**Lesson**: inductive bias isn't "old vs new"; it's "right vs wrong for this data". When transformers beat CNNs (image classification at scale, protein structure), it's because they eventually amass enough data to learn the equivalent local priors themselves. On limited-data genomics tasks, the explicit inductive bias of CNNs often still wins.

**FIGURE — Figure #10: CNN vs Transformer vs MLP at matched training budget** → `diagrams/lecture-16/10-architecture-comparison.svg`
*A line plot. X-axis: training-set size (log, from 10³ to 10⁷). Y-axis: validation F1 on a synthetic pileup-like task, 0.5–1.0. Three curves: MLP (lowest throughout), CNN (strong from the start, early plateau near 0.995), plain Transformer with absolute positional embeddings (steep climb at large data, crosses CNN around 10⁶ samples, pulls ahead at 10⁷). A fourth dotted curve: Transformer with 2D positional + local-attention inductive bias (matches CNN at 10⁵ samples; pulls ahead at 10⁶). Annotation: "the right inductive bias buys 10×–100× data efficiency".*

**EMBED — Artifact #5: Inductive Bias Explorer** → `artifacts/lecture-16/05-inductive-bias.html`
*Pit three architectures (MLP, CNN, Transformer) against a synthetic "pileup" classification task at user-controlled training-set size. Show validation accuracy for each under a matched training budget. Student sees: at small data, CNN's translation-invariance prior wins; at huge data, architectures converge. Target aha: architecture choice is data-regime-dependent; "transformer always wins" is a large-data mirage.*

### 2.3 Why NB likelihood wins over MSE for count data (≈6 min)

Counts from single-cell RNA-seq: mean much less than variance (over-dispersed), often zero-inflated.

If you use **MSE reconstruction** in a VAE:

- Implicitly assumes Gaussian $P(x|\hat{x}) \propto \exp(-(x-\hat{x})^2)$.
- Heavily penalises large deviations. Over-dispersion → constant large deviations → model drags predictions toward a mean that can't be right.
- Predicts negative counts; NaN under log transforms; generally wrong everywhere.

If you use **Poisson reconstruction**:

- $P(x|\hat{x}) \propto \hat{x}^x e^{-\hat{x}} / x!$. Appropriate for integer counts.
- Assumes mean = variance. Single-cell data have variance >> mean. Still under-fits.

If you use **Negative Binomial reconstruction**:

- $P(x|\hat{x}, r) \propto \binom{x+r-1}{x} p^x (1-p)^r$ with mean $\hat{x}$ and variance $\hat{x} + \hat{x}^2 / r$. Over-dispersion captured.
- Matches empirical single-cell count distributions.

If you use **Zero-Inflated Negative Binomial (ZINB)**:

- Mixture: $\pi \cdot \delta_0 + (1 - \pi) \cdot \text{NB}(\hat{x}, r)$. Dropout inflation modelled.
- Slightly better on 10× Chromium-era data where technical dropout is real.

The right likelihood buys orders of magnitude in statistical efficiency. This is a crucial cross-domain EE lesson: **the noise model is the problem-matching decision, not the architecture depth.**

> **EE framing — inductive bias as the signal-processing equivalent of a matched filter**: A matched filter for detection is optimal when the signal template and the noise covariance are known exactly. If you get either wrong — wrong signal shape, wrong noise model — performance collapses. An ML architecture's inductive bias is the generalised version: "template" = architectural symmetries / features the network can efficiently represent; "noise model" = the likelihood chosen for training. Matched filter for a rectangular pulse in white Gaussian noise is the CNN with MSE loss on denoising; matched filter for counts is NB-likelihood VAE; matched filter for 3D rotation-invariant features is equivariant network. Architecture design = matched-filter design at scale.

### 2.4 When to break the rule: transformers eating CNNs (≈5 min)

ViT (Dosovitskiy et al. 2020) showed that transformers beat CNNs on ImageNet-scale data. Why?

- With enough data (100M+ images), transformers learn the equivalent of local priors themselves.
- Transformers' weakness (no built-in locality) becomes a strength (can learn any prior from data).
- Their inductive bias is weaker and more general; at sufficient data scale, this is better.

The same dynamic applies in genomics:

- Protein structure (AlphaFold2): transformers + specific inductive biases (axial, triangle, equivariance) + massive data (150M+ sequences) → SOTA.
- Single-cell foundation models (scGPT, Geneformer): transformers on millions of cells. Work OK but the jury is still out on whether they beat specialised VAEs at matched data.

When does the transformer win? When you have **a lot** of pretraining data (cross-task, cross-species), and the downstream tasks aren't so specialised that a specialised prior helps. For most genomics tasks, the answer today is "not yet — keep using specialised priors."

### 2.5 Equivariant networks: when physics gives you the prior (≈5 min)

If your problem has a known symmetry — rotation of 3D coordinates, permutation of items in a set, translation of a signal — you can **bake the symmetry into the network**.

- **Translation equivariance**: standard CNN property.
- **Rotation equivariance in 2D**: group-convolutional networks (Cohen &amp; Welling 2016), steerable filters.
- **SE(3) equivariance in 3D**: E(3)-equivariant networks (Geiger et al. 2022), NequIP, SchNet, DimeNet. Used in AlphaFold2's IPA.
- **Permutation equivariance / invariance**: GNNs for graphs, Deep Sets (Zaheer et al. 2017) for unordered sets.

Baking symmetry in > learning it via data augmentation. Rule of thumb: factor of 10×–100× data efficiency.

### 2.6 Summary table (≈5 min)

**FIGURE — Figure #4: Inductive-bias table** → `diagrams/lecture-16/04-inductive-biases.svg`
*A three-column table. Column 1: Architecture (CNN / dilated conv+transformer / VAE+NB / GNN / equivariant transformer / fully-connected / transformer base). Column 2: Inductive bias (translation invariance + locality / multi-scale + sparse long-range / negative binomial likelihood + latent prior / permutation invariance + local aggregation / SO(3)-equivariance / none / no locality). Column 3: Data regime where it wins (small-to-medium, position matters; long-range structured; count data with batch effects; molecular graphs; 3D geometry; very large generic / data-rich without structure). Each row colour-coded with the architecture's signature colour.*

---

## Part 3 — Training Data Is the Bottleneck (≈25 min)

### 3.1 Where labels come from in genomics (≈8 min)

A recurring issue: architecture choices are constrained by **how much labelled data** you have. In genomics, this is almost always much less than comparable ML domains.

Sources of labels in genomics, roughly ranked by abundance:

- **Unlabelled sequences** (DNA, protein). **Billions** of entries (UniRef, NCBI, SRA). Labels = "this is a valid biological sequence". Usable for self-supervised pretraining (Part 3.3).
- **Genotype-phenotype associations** from GWAS. Hundreds of thousands of SNPs × hundreds of traits. Noisy; requires biobank-scale cohorts.
- **Transcriptomic data** (bulk RNA-seq). Millions of samples (SRA, GEO). Labels = tissue, condition, phenotype. Noisy annotations; batch effects.
- **Single-cell expression data**. ~100M cells publicly available (cellxgene, HCA). Labels = cluster / cell type — often imperfect; re-clustering changes them.
- **Protein structures** (PDB). ~180k experimental; 200M predicted (AF DB). Predicted structures usable as labels with careful confidence filtering.
- **Gene annotations** (Ensembl, RefSeq). Tens of thousands per genome. Labels = transcript structure; mostly-curated.
- **Variant pathogenicity** (ClinVar). Hundreds of thousands of variants. Labels = pathogenic / VUS / benign; domain-expert-curated; slow to accumulate.
- **Functional experimental data** (deep mutational scans, massively parallel reporter assays). ~1k genes; hundreds of thousands of mutations per study. Gold-standard; tiny.
- **Clinical outcomes**. Access-controlled (dbGaP, EGA). Incomparable across studies.

This is an **inverted pyramid** vs generic ML:

- Text: trillions of tokens, clear labels (next word).
- Images: hundreds of millions with weak labels (hashtags).
- Genomics variant pathogenicity: hundreds of thousands with expert labels.
- Functional outcome: single-digit thousands with experimental labels.

**FIGURE — Figure #5: Labelled-data pyramid in genomics** → `diagrams/lecture-16/05-data-pyramid.svg`
*Inverted pyramid / funnel shape. Top (widest): unlabelled sequences (10¹⁰ entries). Below: bulk RNA-seq (10⁶ samples). Below: GWAS associations (10⁵ SNP × trait pairs). Below: single-cell (10⁸ cells). Below: ClinVar variants (10⁶ with pathogenicity calls). Bottom (narrowest): MPRA / DMS experiments (10⁴ variants per gene). Each tier colour-coded; annotations indicate typical label quality.*

### 3.2 Why this shape matters (≈5 min)

Supervised deep learning wants millions of labelled examples. Most genomics tasks have hundreds to thousands of labels that meet the quality bar. Consequences:

- **Foundation-model strategies dominate**: pretrain on the huge unlabelled sequence corpus, fine-tune on the small labelled target task.
- **Self-supervised objectives matter more** than in generic ML. Masked-language-model on DNA / protein sequences is the workhorse.
- **Evaluation is perilous**: with small test sets, benchmarks are noisy; differences within a few % are often not significant.
- **Transfer between tasks** is the only way many tasks get solved.

### 3.3 Self-supervised pretraining (≈7 min)

**Self-supervised pretraining**: take unlabelled data, invent a label (from the data itself), train.

Two dominant objectives:

- **Masked language modelling** (MLM). Mask ~15% of tokens; predict them from context. BERT (2018) popularised this for text; DNABERT (2021) applied it to DNA k-mers; ESM (2021) applied it to protein sequences.
- **Next-token / causal prediction**. Predict token $t+1$ given tokens $1..t$. GPT for text; HyenaDNA, Evo for DNA; ProGen, ESM-IF for proteins.

After pretraining, the model has learned **contextual representations** — per-position vectors that encode which amino acid / nucleotide goes where, what's statistically expected in the local context, etc.

Downstream uses:

- **Feature extraction**: freeze the pretrained network; use its activations as features for a small downstream classifier. Works well when downstream labels are scarce.
- **Fine-tuning**: continue training the pretrained network on the labelled downstream task. Better performance when downstream labels are moderate.
- **In-context / prompting**: for truly foundation-model-style networks, supply context and ask for the answer without further training. Less established in genomics than in NLP.

> **Intuition box**: Self-supervised pretraining is the "read the whole library before answering any specific question" strategy. Instead of teaching the network exactly what you want from a small labelled set, you let it read billions of unlabelled examples with invented puzzles (guess the masked token, continue the next one). In doing so, it builds a rich representation of "what kinds of sequences exist, how they behave, what's normal vs surprising". When you later hand it a small supervised task, it doesn't start from scratch — it starts from "I already know this language; you're asking me to label one aspect of it". The supervised fine-tuning needs ~100× fewer examples than training from scratch would. This is the single most important trick that made genomics ML work with small labelled datasets.

**EMBED — Artifact #6: Self-Supervised Pretraining Simulator** → `artifacts/lecture-16/06-ssl-pretraining.html`
*Simulate a small transformer learning on a synthetic sequence corpus. Toggle pretraining on/off; set downstream label budget. Watch validation accuracy as downstream training proceeds. Target aha: pretraining gives a large and persistent head start — with 1000 labelled examples, pretrained matches scratch-trained at 100,000 labelled examples.*

**FIGURE — Figure #6: Self-supervised pretraining flow** → `diagrams/lecture-16/06-ssl-flow.svg`
*Left: massive unlabelled sequence corpus (UniProt ~650M proteins). Middle: pretraining — mask some tokens, predict them; iterate over millions of sequences. Right: the pretrained model is frozen / fine-tuned on a small labelled downstream task (e.g. "predict which variants are pathogenic"; label count ~50k). Annotation: "scale in unlabelled pretraining bought ~10× more downstream efficiency than training from scratch".*

> **EE framing — self-supervised pretraining as representation learning from unlabelled data**: Self-supervised pretraining is the generalised answer to "your labels are expensive but your unlabelled data is abundant". You invent a pretext task (mask some, predict), optimise it across unlabelled data, and learn a representation that's generally useful for downstream tasks. The analogy in signal processing: unsupervised feature learning (k-means, sparse coding, PCA) from unlabelled data before supervised classification. Modern variant: contrastive learning (SimCLR), which extends MLM's "predict what was hidden" idea to "predict whether these two augmentations are of the same example".

### 3.4 Weak supervision and noisy labels (≈5 min)

When labels are scarce, use noisy ones:

- **Distant supervision**: infer labels from auxiliary data (e.g. gene expression inferred from bulk tissue samples; protein function labels inferred from GO annotations, which are themselves partially noisy).
- **Heuristic labels**: rule-based labelling of large corpora (e.g. variants annotated as "likely pathogenic" by a rule rather than experts).
- **Multi-task learning**: train on related tasks with more labels, transfer representations.
- **Active learning**: prioritise requesting expert labels for the samples where the model is most uncertain.

All of these are standard ML methods, particularly valuable in genomics because label scarcity is structural.

---

## Part 4 — The Data Leakage Problem (≈30 min)

### 4.1 Why genomics leaks worse than almost any domain (≈6 min)

**Data leakage**: the phenomenon where train and test sets share information they shouldn't, producing falsely-optimistic performance. Standard ML advice ("shuffle and split 80/20") is catastrophic in genomics because genomic data is **highly correlated**.

Sources of correlation:

- **Sequence homology**: protein sequences from related species share 60–95% identity. A random 80/20 split puts closely-related sequences on both sides. The test set "answers" are trivial to predict if the training set has a homologue.
- **Genomic locality**: nearby positions on a chromosome are in LD. A per-SNP random split breaks LD structure but puts correlated SNPs on both sides. Models can memorise the correlation and look like they've learned something.
- **Family relatedness**: cohort individuals are often distant relatives. Sib-pairs and cousin-pairs share segments of DNA. Per-individual random splits leak via relatedness.
- **Batch effects across studies**: samples processed together have systematic technical signatures. Random splits across studies put correlated batch artefacts on both sides; the model learns "which lab?" rather than "which phenotype?".

Consequence: **published numbers are often dramatically optimistic**. Papers reporting 99% test accuracy on naïvely-split data frequently drop to 60–80% on properly-split data.

> **EE framing — data leakage as violated independence in correlated samples**: Every train/test evaluation assumes the samples are drawn i.i.d. from the same distribution, and that train / test are independent draws. In genomics, neither assumption holds. Neighbouring SNPs are autocorrelated along a chromosome; homologous sequences are autocorrelated across species; samples from the same lab are autocorrelated via batch effects. This is the same statistical issue as **spatial autocorrelation** in geostatistics or **temporal autocorrelation** in time-series evaluation: a random hold-out from correlated data gives you no meaningful test because the test points are statistically tethered to the training points. The fix is the same as in spatial / time-series ML: **blocked cross-validation** — hold out whole blocks (chromosomes, families, studies, time windows) that share correlation structure internally but are independent of other blocks. Get this wrong and your reported accuracy measures how well the model memorised the training distribution, not how well it generalises.

**FIGURE — Figure #7: Data leakage — random vs chromosome split** → `diagrams/lecture-16/07-leakage-splits.svg`
*Two panels side-by-side. Left: "Random split" — a chromosome bar with train (blue) and test (amber) positions interleaved at random. Annotation: "adjacent SNPs are in LD; test info is trivially predictable from train". Right: "Chromosome split" — chromosome 1..21 = train; chromosome 22 = test. Spatially separated; no LD bleed. Annotation: "evaluates generalisation to unseen regions, not memorisation of LD structure".*

### 4.2 Homology-based splits for protein tasks (≈6 min)

**CAFA** (Critical Assessment of Functional Annotation, 2010+): the community benchmark for protein-function prediction. Uses strict time-based splits: train = sequences known by date T; test = sequences added by date T + 2 years. Avoids sequence leakage via time.

**Sequence identity clustering**: cluster all sequences by CD-HIT or MMseqs2 at 30–50% identity; split cluster-wise. Each test cluster contains only sequences that share <50% identity with any training cluster. This is the **MCMD standard** for protein task splits.

Common mistakes:

- Splitting by sequence ID without clustering: homologous sequences on both sides.
- Splitting by UniProt accession randomly: same thing.
- Splitting by species: better but not enough if species share homologues.

Rule of thumb: **if a random baseline scores 60% on your test set, you have leakage, not a model.**

### 4.3 Chromosome splits for sequence regulation (≈5 min)

For Enformer-style long-range regulation prediction:

- Random per-position split: catastrophic — adjacent positions are in LD and share regulatory context.
- Per-window split (e.g. 100 kb windows): still leaky if test windows are within the Enformer receptive field of training windows.
- **Chromosome split**: hold out entire chromosomes for test. Enformer uses chr8 + chr9 for test. This is the community standard.
- **Whole-species split**: for cross-species generalisation studies, train on human, test on mouse. Stricter; rarely done for foundation-level benchmarks but more honest.

### 4.4 Cryptic relatedness in biobank cohorts (≈5 min)

UK Biobank has ~500k individuals. By chance, thousands of pairs are cousins-or-closer. Simple random splits leak via these pairs:

- A case in training + a second-cousin control in test = the model learns to distinguish them by shared haplotype patterns.
- **Kinship-aware splitting**: compute the kinship matrix; ensure no train-test pair has kinship > 0.05.

GWAS evaluation pipelines always implement this. ML pipelines often don't, resulting in optimism.

> **Intuition box**: The right mental image for genomics data leakage is "adjacent pixels in an image". If you hand out a photo to a model pixel-by-pixel — train on 80% of random pixels, test on the other 20% — the model barely needs to look at anything; it just averages the neighbouring training pixels. You haven't tested anything. Genomic data is structurally identical: SNPs along a chromosome are the "adjacent pixels"; related individuals are "adjacent patients"; samples from one lab run are "adjacent batches". Random hold-outs on any of these produce an exam you've already seen the answers to. Leakage-aware splits hold out *whole regions* of this correlation structure — a chromosome, a family, a lab — so the test genuinely probes generalisation.

### 4.5 Batch effects across studies (≈4 min)

Single-cell data: a given cell-type classifier trained on some studies and tested on others consistently over-performs vs. a classifier trained on most-studies-but-held-out-one-study. Reason: the model learns "this came from Study X". Fixes:

- **Leave-one-study-out** evaluation: each held-out study evaluates cross-study generalisation.
- **Batch correction before split**: remove known batch effects; remaining signal is biology.
- **Benchmarks that enforce this**: Luecken et al. 2022 (benchmarking single-cell integration) explicitly uses leave-one-study-out.

### 4.6 Real-world case studies (≈4 min)

Two widely-cited cases of published work affected by leakage:

- **Protein-function prediction benchmarks in the 2010s**: naive sequence splits produced papers claiming 95%+ accuracy; CAFA-style time-based splits showed ~40–60% — the difference was leakage.
- **Early cancer-driver prediction** (2014–2018 era): some tools trained and tested on overlapping known-driver databases, claiming 99% accuracy. Independent benchmarks dropped them to 70%.
- **Enformer's initial-release test set** (2021): passed chromosome splits; but a 2023 reanalysis showed within-chromosome regulatory correlation meant even chromosome splits leaked to some degree. Tightened splits exist.

**FIGURE — Figure #11: Leakage case studies — random vs corrected accuracy** → `diagrams/lecture-16/11-leakage-cases.svg`
*A grouped bar chart with three case studies on the x-axis. For each case study, two bars: "as reported" (tall) and "after tightened split" (much shorter). Case 1: "Protein function prediction, ~2015 era": 95% → 45%. Case 2: "Cancer driver prediction, ~2016 era": 99% → 70%. Case 3: "Enformer 2021 chromosome-split": 0.82 Pearson → 0.74 Pearson. Annotation below the plot: "reproducing any published genomics-ML paper should start by auditing its split; the corrected number is usually the honest one".*

> **Warning box**: When you see a genomics ML paper with accuracy that seems too good to be true, your first question should be: **how was the train/test split done?** If the answer is "random 80/20" or isn't explained in detail, discount the number heavily. The paper's method may still be sound — but the reported metric is likely inflated. Reproducing the paper with a leakage-aware split is a standard sanity check, and it frequently changes the headline number.

> **Discussion prompt**: You're reviewing a genomics ML paper claiming state-of-the-art on a variant-pathogenicity classifier. The paper's train/test split is: "ClinVar dataset randomly split 80/20 stratified by gene". Identify three potential leakage sources. Suggest a better splitting strategy. What would you predict the paper's true held-out accuracy to be if the splits were tightened? (Hint: leakage via variant neighbours in the same gene; leakage via related genes with shared structural features; leakage via the training corpus containing variants highly correlated with test variants in the same family — think about sequence / structural / functional proximity.)

**EMBED — Artifact #2: Data Leakage Demonstrator** → `artifacts/lecture-16/02-data-leakage.html`
*Simulate a simple variant-effect-prediction task with controllable sequence correlation between nearby positions. Show test performance under: random-position split, chromosome split, clustered-by-homology split. Target aha: the accuracy gap can be 30+ percentage points — "random split" is a well-intentioned disaster on correlated genomics data.*

---

## Part 5 — DNA Language Models (≈25 min)

### 5.1 What DNA LMs are (≈5 min)

A **DNA language model (DNA-LM)** is a large transformer (or transformer-variant) pretrained on genomic sequence, with an objective similar to masked-language-modelling in NLP:

- Chop the genome into overlapping windows of some length (500 bp – 1 Mb, depending on the model).
- Tokenise (k-mer, byte-pair-encoded, or per-nucleotide).
- Mask some tokens; predict the masked ones from context. (Or: predict next token causally.)
- Pretrain for millions of steps on billions of bases.

Downstream use: fine-tune or feature-extract for specific tasks (promoter detection, TF binding, variant effect, regulatory element classification).

Key entries (2021–2024):

- **DNABERT** (Ji et al. 2021). BERT-style; k-mer tokenisation; 512-token context. ~100M parameters.
- **Nucleotide Transformer** (Dalla-Torre et al. 2023, InstaDeep). Larger (~500M–2.5B parameters); multi-species pretraining; context 1000 tokens.
- **HyenaDNA** (Nguyen et al. 2023, Stanford). Replaces attention with Hyena's long convolutions. Scales to 1 million base context. Sub-quadratic cost.
- **Evo** (Nguyen et al. 2024, Arc Institute). Hyena-style, trained on millions of genomes across prokaryotes, eukaryotes. Generates biologically-plausible full genomes. Context up to ~131k.

**FIGURE — Figure #8: DNA-LM landscape** → `diagrams/lecture-16/08-dna-lm-timeline.svg`
*Timeline 2021–2024 with key DNA LMs plotted by year (x-axis) and context length in bases (y-axis, log scale). DNABERT at 512 bp, 2021. Nucleotide Transformer at 12 kb, 2023. HyenaDNA at 1 Mb, 2023. Evo at 131 kb, 2024. Each point with model parameter count. Annotation: "context length is the main axis of competition — it determines what scale of regulation is learnable".*

### 5.2 What DNA pretraining gives you (≈5 min)

The honest account of what DNA LMs buy:

- **Motif discovery** at no cost. Trained models can be probed to recover known TF binding motifs (Alipanahi et al. 2015 showed this is possible; DNA LMs do it naturally).
- **Moderate improvements on regulatory prediction**. Downstream tasks (promoter detection, enhancer prediction, splice site prediction) see 1–5% absolute improvements over specialised architectures. Not dramatic; real.
- **Cross-species transfer**. Models pretrained across species generalise better on held-out species than species-specific models.
- **In-silico mutagenesis**: differential scoring of reference vs mutant sequences. Faster than SHAP-style attribution on specialised models.

### 5.3 What DNA pretraining does NOT give you (≈5 min)

Honest caveats:

- **No uniform large wins over specialised methods**. Enformer (dedicated architecture for 100 kb regulatory prediction) still beats DNA LMs on the enformer-style tasks. DNA LMs are improving but haven't collapsed the field.
- **Generative quality is suspect**. DNA LMs generate plausible-looking sequences; unclear whether generation respects subtle biology (e.g. whether a generated "enhancer" actually would function).
- **Long-range dependency isn't free**. HyenaDNA / Evo have long context in name, but empirically, attention-like interactions between distal elements are still weakly learned.
- **Species bias**. Most DNA LMs are dominated by human + a few model organisms. Performance on bacteria / non-model organisms lags.
- **Evaluation-set contamination**. DNA LMs' evaluation tasks (BENDBenchmark 2023) may overlap sequences the LMs saw during pretraining. Audit carefully.

### 5.4 Why DNA is NOT like text (≈5 min)

A critical EE point: DNA is **not** just text in a 4-letter alphabet. Differences that matter:

- **Orders of magnitude more "text"**. Human genome ~3 billion bp; prokaryotic pan-genomes are effectively unbounded. Pretraining corpora are tens of TB.
- **Very low information density per token** vs English words. Most DNA is intergenic / non-coding / repetitive. Effective information per position is much lower.
- **Functional structure is multi-scale and sparse**. A promoter motif is 10 bp; an enhancer is 300 bp; a TAD is 1 Mb. Attention over 1 Mb to find 10 bp signals is a needle-in-a-haystack problem.
- **Frame-shift and reverse-complement symmetries**. The signal on the + strand is the reverse-complement of the − strand. Language models usually need explicit data augmentation; text-style pretraining misses this.
- **Mutation statistics differ from language editing statistics**. DNA variants are SNP / indel / SV; text edits are character substitutions, deletions, swaps. Pretraining objectives designed for text may not align with biology.

> **EE framing — DNA LMs as language models on a 4-letter alphabet with very different statistics**: A naive analogy with NLP foundation models predicts large gains from scale. In practice, the gains have been modest. Reasons: the signal-to-noise ratio per base is much lower than per word; functional elements are much sparser than in natural language; tokenisation choices (k-mer vs per-nucleotide vs BPE) materially affect what's learnable. These aren't reasons to give up — they're reasons to be careful with the analogy. DNA is its own modality, with its own statistics, and the field is still working out which pretraining objectives work best.

### 5.5 Current state of the art (≈5 min)

As of 2025:

- For very long-range prediction (>100 kb): Enformer-class specialised models + limited DNA-LM contributions.
- For mid-range (10 kb) prediction: DNA LMs starting to match specialised methods.
- For short-range (<1 kb) motif discovery: DNA LMs are competitive.
- For generation (designed sequences): Evo and other generative DNA LMs are the frontier; experimental validation still ongoing.
- For cross-species transfer: DNA LMs clearly win.

The honest assessment: DNA LMs are a useful tool in the genomics ML toolkit but haven't transformed the field the way ESM transformed protein prediction or GPT transformed NLP.

**EMBED — Artifact #3: DNA-LM Tokeniser Explorer** → `artifacts/lecture-16/03-dna-lm-tokeniser.html`
*Take a preset DNA sequence, tokenise it three ways (per-nucleotide, 3-mer, 6-mer, BPE-style). Show vocabulary size, typical token count for a fixed sequence, and how the receptive field of a 512-token window changes. Target aha: tokenisation choice determines what scale of feature the model can attend over — per-nucleotide is finest but slowest; k-mer is coarser but efficient.*

---

## Part 6 — Foundation Models for Cells (≈25 min)

### 6.1 The pitch (≈5 min)

**scGPT** (Cui et al. 2024), **Geneformer** (Theodoris et al. 2023), **scFoundation** (Hao et al. 2024), **UCE / Universal Cell Embedding** (Rosen et al. 2024). These are transformer-based models pretrained on **millions of single-cell transcriptomes** (tens to hundreds of millions of cells).

The pitch: "scGPT / Geneformer is to single-cell what GPT is to text." A foundation model whose pretrained representations transfer to **any** downstream single-cell task: cell-type classification, perturbation response, gene-gene interaction, cross-species alignment.

Each model takes a gene-expression vector, represents it somehow (ranked gene list, or value-binned, or tokenised) as a sequence, and runs transformer blocks. Different tokenisations:

- **Geneformer** (2023): rank-value encoding — for each cell, order genes by expression; use the rank as the token.
- **scGPT** (2024): binned expression values + gene identity as parallel token channels.
- **scFoundation** (2024): per-gene real-valued expression via a custom embedding.
- **UCE** (2024): cross-species shared embedding.

### 6.2 What works (≈5 min)

Strongest claims (and evidence):

- **Zero-shot cell-type classification**: classify a cell's type without any task-specific training. Some models achieve 70–80% accuracy on held-out studies. Decent but not replacing specialised classifiers.
- **Perturbation response prediction**: given a cell's baseline expression + a knock-out, predict the response. scGPT shows moderate success.
- **Gene-gene interaction discovery**: extract which genes' representations are "connected" through attention. Recovers known interactions.
- **Cross-species transfer via shared embedding**: UCE shows mouse-human cell-type alignment at improved accuracy.

### 6.3 What's oversold (≈6 min)

Honest critique:

- **Evaluation is immature**. Many papers report zero-shot accuracy on cell types that were present in pretraining (under different names). Proper leave-one-study-out evaluations give much lower numbers.
- **Specialised methods still win on specialised tasks**. scVI for integration; Harmony for batch correction; CellTypist for cell-type classification. Foundation models are not yet the default.
- **Cost is high**. Pretraining takes thousands of GPU-hours; fine-tuning takes hundreds; downstream utility is not yet proportional.
- **Biological interpretability is limited**. The "learned representations" are not straightforwardly interpretable as gene-regulatory networks or pathways.
- **The "GPT analogy" is loose**. Text has semantic meaning per token; gene-expression vectors encode a cell's state, which is a different kind of object. The tokenisation choices (rank, binned value, etc.) are still being figured out.

### 6.4 Where these models might matter (≈5 min)

Honest optimistic view:

- **Cross-species alignment** at scale: UCE's direction is promising, and matters for translating findings between model organisms and humans.
- **Patient-level predictions from single-cell data**: combining sc-data from multiple patients into a shared embedding space may enable disease-state classification the way radiomics enabled imaging-based diagnosis.
- **Perturbation screens at virtual scale**: pretrained representations + a small amount of experimental data could predict CRISPR screen outcomes across cell types, avoiding some experimental cost.
- **Multi-modal integration**: foundation models that cover scRNA + scATAC + protein are the likely future; current-gen foundation models are a stepping stone.

### 6.5 How to read a foundation-model paper honestly (≈4 min)

Questions to ask:

1. **What was the pretraining corpus?** Which cells, from which studies, at what technology? Any held-out studies for evaluation?
2. **What's the evaluation protocol?** Leave-one-study-out or leave-one-cell-out? The former is hard; the latter is trivial.
3. **What's the baseline?** Compared to scVI, or to "identity mapping"? The right baseline is the current best specialised tool, not the trivial baseline.
4. **What's the downstream task?** The model might be great at cell-type classification but poor at perturbation prediction. Conflating these produces misleading headline numbers.
5. **How much fine-tuning is involved?** "Zero-shot" with heavy pre-run embedding lookup is effectively fine-tuning.

**FIGURE — Figure #9: Cell foundation model architecture** → `diagrams/lecture-16/09-cell-foundation-model.svg`
*A schematic of scGPT / Geneformer style. Input: a cell's gene-expression vector (tokenised as rank-value or binned pairs). Transformer stack (12–24 blocks). Output: cell-level embedding for any downstream task. Right-side annotation box: pretraining corpus (~30M cells); downstream tasks (cell-type classification, perturbation, cross-species). Note: "specialised methods still win on specialised downstream tasks".*

> **Warning box**: Foundation models for cells are marketed as transformative but, as of 2025, the evidence is mixed. A responsible evaluation compares them on the *exact* benchmark and splitting protocol used by the specialised method being compared to. Many of the "scGPT beats X" headlines don't survive that comparison. Use foundation-model cell embeddings as one tool, not as a replacement for specialised methods.

**EMBED — Artifact #4: Cell Foundation Model Zero-Shot Explorer** → `artifacts/lecture-16/04-cell-foundation.html`
*Preset single-cell embeddings from a scGPT-like foundation model. Student selects a cell type and evaluation protocol (leave-one-cell-out vs leave-one-study-out). See the accuracy gap between the two protocols. Target aha: foundation-model zero-shot accuracy is a function of evaluation protocol, not a single number.*

---

## Part 7 — What's Coming (≈15 min)

### 7.1 Multimodal foundation models (≈4 min)

The likely next wave: foundation models that cover **more than one modality** in a single latent space.

- Protein sequence + structure + function (ESM-3 already, with richer downstream): one model that handles all three.
- DNA sequence + chromatin state + RNA expression: pretraining on cross-modality aligned data (e.g. matching DNase-seq + RNA-seq per cell).
- scRNA + scATAC + protein: true multi-modal cell foundation models.

The engineering recipe is clear (separate encoders per modality → shared transformer → per-modality heads); the data is the bottleneck. Datasets with matched multi-modal measurements are small.

**FIGURE — Figure #12: Multimodal foundation-model landscape** → `diagrams/lecture-16/12-multimodal.svg`
*A 2D axes. X-axis: modalities covered (single → 2 → 3+). Y-axis: parameter count (log, 100M → 100B). Dots for existing models: ESM-2 (protein seq only, ~15B), AlphaFold3 (protein + nucleic + ligand, moderate params), DNABERT (DNA only, 100M), Evo (DNA only but cross-species, ~7B), scGPT / Geneformer (cell + RNA, 10–50M), ESM-3 (seq + structure + function, ~100B). Shaded "frontier" region in the top-right (high-param, multi-modal). Annotation: "the frontier moves toward multi-modal — data availability, not model scale, is the current bottleneck".*

### 7.2 Protein-DNA co-design (≈4 min)

Combining AlphaFold (protein) + DNA-LM (nucleic acid) to design **both a protein and its target DNA / RNA sequence** jointly:

- Design a TF + its binding site.
- Design an RNA-binding protein + its bound RNA.
- Design a CRISPR variant + its guide RNA preference.

RFDiffusion-All-Atom (2024) and extensions are starting to do this. AlphaFold3 handles the structure-prediction side.

### 7.3 Generative models for designed biology (≈3 min)

Beyond prediction: generate novel biology.

- **Designed proteins** (L15): RFDiffusion, ProteinMPNN. Experimental validation is catching up.
- **Designed DNA / RNA**: Evo and successors can generate plausible genome-length sequences. Whether they function biologically is mostly untested.
- **Designed cells / whole organisms**: speculative; current tech far from this.

The dual-use question (L15 §6.4) becomes more acute with each generation.

### 7.4 Convergence with mainstream AI (≈4 min)

Historically, genomics ML was a niche that adapted techniques from mainstream ML with 2–3 year lag. The gap is closing:

- Protein structure (AlphaFold) has become mainstream ML.
- Transformer adoption in genomics has been fast.
- DNA LMs benefit from all of NLP's progress.
- Scaling laws, efficient training, inference optimisation — all cross-applicable.

The emerging scenario: genomics ML = mainstream ML + biology-specific inductive biases + biology-specific data challenges. Not a separate field.

> **Historical pointer**: The DeepVariant paper (Poplin et al. 2018) is a good marker for when genomics ML seriously entered the mainstream. Before 2018, most bioinformatics tools were hand-engineered statistical methods (GATK's HaplotypeCaller is essentially an HMM + hand-crafted features). DeepVariant showed that a plain 2D CNN trained on pileup images could match or beat a decade of careful hand-engineering — and importantly, retrain for a new sequencing technology in a few GPU-days. The 2018–2021 window saw this pattern repeat across the course's domains: scVI for single-cell (2018), DeepSEA / Enformer for regulation (2015 / 2021), AlphaFold for structure (2018 / 2021). The underlying driver was the same: enough labelled / pretraining data had accumulated that the right deep-learning architecture beat hand-engineered baselines.

**EMBED — Artifact #7: Genomics ML Pitfall Checklist** → `artifacts/lecture-16/07-pitfall-checklist.html`
*An interactive checklist where the student describes a proposed genomics ML project (architecture, data source, split strategy, evaluation metric, baselines). The artifact flags likely pitfalls: "this split looks leaky", "this baseline is too weak", "this loss function mismatches the data's noise model". Target aha: most genomics ML design bugs are recognisable from a checklist — the mistakes are predictable, not exotic.*

---

## Wrap-up (≈10 min)

### What you should take away

- **Architecture choice is inductive-bias matching**. For each genomics problem, the structural property of the data (translation invariance, long-range sparsity, count over-dispersion, rotation symmetry) should dictate the architecture.
- **Right likelihood > deep architecture**. scVI's NB likelihood is what matters; the VAE scaffold is scaffolding. Signal-processing analogy: matched-filter correctness > filter complexity.
- **Data is the bottleneck**. Unlabelled sequences are abundant; labelled-outcome data is scarce. Self-supervised pretraining + fine-tuning is the dominant pattern.
- **Data leakage is catastrophic in genomics and widely under-acknowledged**. Random splits leak via homology, LD, relatedness, batch. Leakage-aware splits (chromosome-wise, family-wise, time-based) are mandatory for honest evaluation.
- **DNA language models are useful but not transformative**. They buy moderate gains on regulatory prediction; they don't yet beat specialised architectures at the tasks those were built for.
- **Foundation models for cells are in their infancy**. The marketing outpaces the evidence; zero-shot evaluation is still maturing; specialised methods remain competitive.
- **Inductive bias design is THE skill**. Architecture choice + likelihood + split design + pretraining strategy together define the result.

### Next lecture

Clinical genomics, variant interpretation, and ethics. ACMG/AMP variant classification; the evidence ecosystem (population frequency, functional predictors, splice predictors, ClinVar); pharmacogenomics; incidental findings; FDA / regulatory landscape; ancestry bias and data sovereignty.

### Homework

1. Pick a genomics ML paper from 2022+ with a reported 95%+ accuracy. Reproduce the data split; design a tighter split (chromosome / family / time-based); re-evaluate. Report the accuracy delta.
2. Take a DNA-LM (DNABERT or Nucleotide Transformer via HuggingFace) and fine-tune it on a small downstream task (e.g. promoter classification from DeepProm data). Compare against a CNN trained from scratch on the same data. Which wins, by how much?
3. Train a simple VAE (PyTorch, 50 lines) on a sc-RNA dataset with MSE, Poisson, and NB likelihoods. Compare reconstruction quality. Which produces more biologically-coherent clusters on UMAP?
4. Read the scGPT paper and critically audit its benchmark. Does the evaluation use leave-one-study-out? How are "zero-shot" results framed? What baselines are compared? Write a one-page honest assessment.
5. Design an architecture + inductive bias for a novel problem: predicting chromatin 3D contacts at 10 kb resolution from 1 Mb of surrounding sequence. Justify your choice of architecture from first principles.

### Recommended reading

- Poplin, R., Chang, P., Alexander, D., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology* 36, 983–987. (DeepVariant.)
- Avsec, Ž., Agarwal, V., Visentin, D., et al. (2021). Effective gene expression prediction from sequence by integrating long-range interactions. *Nature Methods* 18, 1196–1203. (Enformer.)
- Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., &amp; Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods* 15, 1053–1058. (scVI.)
- Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
- Dosovitskiy, A. et al. (2020). An image is worth 16x16 words: Transformers for image recognition at scale. *ICLR 2021*.
- Ji, Y. et al. (2021). DNABERT: pre-trained Bidirectional Encoder Representations from Transformers model for DNA-language in genome. *Bioinformatics* 37, 2112–2120.
- Dalla-Torre, H. et al. (2023). The Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *bioRxiv*.
- Nguyen, E. et al. (2023). HyenaDNA: long-range genomic sequence modeling at single nucleotide resolution. *NeurIPS 2023*.
- Nguyen, E. et al. (2024). Sequence modeling and design from molecular to genome scale with Evo. *Science* 386.
- Theodoris, C. V. et al. (2023). Transfer learning enables predictions in network biology. *Nature* 618, 616–624. (Geneformer.)
- Cui, H. et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* 21, 1470–1480.
- Rosen, Y. et al. (2024). Universal Cell Embeddings: A foundation model for cell biology. *bioRxiv*.
- Luecken, M. D. et al. (2022). Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* 19, 41–50.
- Friston, K. J. et al. (2022). A free energy principle for biological self-organization. (Representation-learning perspective; not genomics-specific but relevant.)
- Teichmann, S. A., &amp; Efremova, M. (2020). Method of the year 2019: single-cell multimodal omics. *Nature Methods* 17, 1 (Editorial).
- Papers With Code — Genomics: <https://paperswithcode.com/area/biology>
- OpenFold: <https://github.com/aqlaboratory/openfold>
- HuggingFace Hub (DNA models): <https://huggingface.co/models?search=dna>
- scGPT repository: <https://github.com/bowang-lab/scGPT>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Architectural Pattern Survey                          | 45&nbsp;min | 0:45 |
| Part 2 — Inductive Biases: Why Does Each Work?                   | 35&nbsp;min | 1:20 |
| Part 3 — Training Data Is the Bottleneck                          | 25&nbsp;min | 1:45 |
| Part 4 — The Data Leakage Problem                                  | 30&nbsp;min | 2:15 |
| Part 5 — DNA Language Models                                       | 25&nbsp;min | 2:40 |
| Part 6 — Foundation Models for Cells                                | 25&nbsp;min | 3:05 |
| Part 7 — What's Coming                                              | 15&nbsp;min | 3:20 |
| Wrap-up                                                              | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
