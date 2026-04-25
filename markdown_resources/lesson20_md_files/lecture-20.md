# Lecture 20 (proposed L4) — Multiple Sequence Alignment, Phylogenetics, and Comparative Genomics

> **Duration**: ≈3h 55min content
> **Audience**: EE undergraduates / graduates, minimal biology assumed
> **File**: to be rendered as `lectures/lecture-20.html` (provisional name; renumber to `lecture-04.html` when curriculum is reordered)

> **Proposed placement**: insert after the new BLAST lecture, becoming the new L4. The natural arc is L2 (pairwise alignment) → L3 (BLAST: one-vs-many search) → L4 (MSA: many sequences, families, evolution). MSA is the input to profile HMMs (new L7), AlphaFold (existing L15→L19), and almost every comparative-genomics workflow. Phylogenetics gives the vocabulary for evolution that recurs in cancer clonal trees (L18→L22), microbial evolution, and population genetics (L12→L16).

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Formulate the multiple-sequence-alignment (MSA) problem; explain why exact dynamic programming is intractable beyond ~5 sequences and why progressive alignment (Clustal, MUSCLE, MAFFT) is the practical compromise.
2. Walk through progressive alignment with a guide tree; explain Feng-Doolittle "once a gap, always a gap" and its consequence (early errors freeze).
3. Apply iterative refinement (MUSCLE) and the Carrillo-Lipman bound; understand why MAFFT's FFT-based approach scales.
4. Construct phylogenetic trees by parsimony, by distance methods (UPGMA, Neighbor-Joining), and by maximum likelihood; describe the molecular-clock assumption and its violation.
5. Compute and interpret dN / dS ratios; identify positive vs purifying selection at the gene level.
6. Describe synteny, ortholog identification (RBH, OrthoFinder), and paralog distinction; navigate UCSC Genome Browser and Ensembl Compara.
7. Frame phylogenetics and MSA in EE terms: trees as graphical models; MSA as constrained DP; molecular clock as a Poisson process; dN/dS as a likelihood-ratio test.

---

## Part 1 — From Pairwise to Multiple (≈25 min)

### 1.1 Why MSA matters (≈5 min)

Pairwise alignment (L2) handles two sequences. But biology routinely deals with **families**: a gene family across 50 species, all isoforms of a protein, or all members of a Pfam domain. To compare these, we need multiple sequences aligned in **a single coherent column structure** where homologous residues line up.

Concrete uses:

- **Phylogenetics**: trees are inferred from MSAs.
- **Profile HMMs / Pfam**: column-by-column residue conservation gives a profile.
- **AlphaFold**: MSA columns provide coevolution signals.
- **Conserved-region detection**: identifies functionally-important residues.
- **Primer design**: degenerate primers cover all family members.

### 1.2 The MSA problem (≈6 min)

Given $N$ sequences of average length $L$, find the alignment with maximum **sum-of-pairs** (SP) score:

$$\text{SP} = \sum_{i < j} \text{score}(\text{aln}_i, \text{aln}_j)$$

Each pair's contribution is the standard pairwise score; the goal is jointly optimal placement of gaps so all pairs simultaneously score well.

### 1.3 Why exact DP is impossible (≈6 min)

Pairwise DP runs in $O(L^2)$ time. The natural generalization to $N$ sequences runs in $O(L^N)$ time and space — the **Carrillo-Lipman tensor**. For $N = 5$ and $L = 300$, that's $\sim 10^{12}$ states. For $N = 50$ (a small family), it's $10^{125}$. **NP-hard in general** (Wang & Jiang 1994).

So all practical MSA tools are **heuristic**.

### 1.4 The Carrillo-Lipman bound (≈4 min)

Carrillo & Lipman (1988): a partial bound that lets exact DP run for small $N$ and short $L$. The idea: prune branches of the DP tensor that can't beat the current best. Useful for $N \leq 5$ or so. Not used at scale; mentioned because it tells you exact MSA isn't entirely impossible — just unscalable.

### 1.5 The progressive-alignment idea (≈4 min)

Feng & Doolittle (1987): instead of optimizing all sequences jointly, align them **in tree order**:

1. Build a guide tree from rough pairwise distances.
2. Align the two closest sequences first.
3. At each internal node, align the **profile** built from previously-aligned sequences with the next sequence.
4. End with one alignment of all $N$ sequences.

This is suboptimal: errors made early are frozen ("once a gap, always a gap"). But it's tractable (close to $O(N^2 L^2)$) and works well in practice.

**FIGURE — Figure #1: Pairwise → MSA generalization** → `diagrams/lecture-20/01-pairwise-to-msa.svg`
*Three panels. Left: pairwise alignment of two sequences — DP matrix with diagonal traceback. Middle: 5-sequence DP tensor in 5D — visualised as nested cubes; total state count exploding. Right: progressive alignment with a guide tree — pairwise step then profile-vs-sequence then profile-vs-profile, showing how complexity is decomposed.*

---

## Part 2 — Progressive Alignment in Practice (≈40 min)

### 2.1 The Clustal family (≈5 min)

**ClustalW** (1994), **ClustalX** (1997), **Clustal Omega** (2011). The progenitors of progressive alignment, still widely used.

- Compute all pairwise distances (typically by k-tuple-counting, fast but rough).
- Build a guide tree (UPGMA or NJ).
- Progressive alignment along the tree.

Limitation: the guide tree has to be built before any alignment is done, but the distance estimates from raw sequences can be inaccurate at distant homology.

### 2.2 MUSCLE: iterative refinement (≈10 min)

**MUSCLE** (Edgar 2004) added two key innovations:

1. **Better guide tree**: MUSCLE iterates — build a quick alignment, derive distances from the alignment, build a new tree, re-align.
2. **Iterative refinement**: after the initial alignment, **bipartition** the tree at a random internal edge; re-align the two halves separately, then combine. If the SP score improves, accept; else reject. Repeat ~20 iterations.

Iterative refinement breaks the "once a gap, always a gap" lock-in by giving early errors a chance to be corrected later.

### 2.3 MAFFT: FFT-based scoring (≈10 min)

**MAFFT** (Katoh et al. 2002, 2013) is the modern default for large MSAs:

- For pairwise comparisons, uses **fast Fourier transform (FFT)** on physicochemical-property representations of residues. Approximates the distance computation in $O(L \log L)$ rather than $O(L^2)$.
- Multiple algorithm choices: FFT-NS-1 (fast), FFT-NS-2 (balanced), L-INS-i (most accurate, uses local pairwise alignments as anchors), G-INS-i (global).
- Scales to thousands of sequences.

For typical workflows: MAFFT L-INS-i is the gold standard; FFT-NS-2 the practical default.

### 2.4 The MSA quality landscape (≈5 min)

Benchmarks (BAliBASE, PREFAB, OXBench, HOMSTRAD) compare MSA tools:

- **MAFFT L-INS-i** and **PROBCONS**: highest accuracy, slower.
- **MUSCLE**, **Clustal Omega**: balanced.
- **Kalign**: fast for large datasets.

Differences are most visible at distant homology (~25% identity). Closely-related sequences are well-aligned by any tool.

### 2.5 The deep dive (≈5 min)

> **EE framing — MSA as constrained tensor DP**: The exact MSA problem is **dynamic programming on a $D$-dimensional tensor** where $D = N$. Each cell stores the best alignment up to that prefix combination. The transition is a max over $2^N - 1$ neighbours (each subset of sequences advancing or staying). Exponential in $N$ both in time and memory. Heuristic algorithms work by **decomposing the tensor**: progressive alignment computes 2D slices, then 2D-vs-1D slices, sequentially. Iterative refinement is essentially **coordinate-descent optimisation** on the SP-score landscape, restarting from a partition of the current solution. The shape of the solution space — heavily local-minima-ridden, gap-extension non-convex — explains why no closed-form algorithm exists and why empirical benchmarking is the gold standard.

### 2.6 Worked example (≈5 min)

**FIGURE — Figure #2: Progressive alignment step-by-step** → `diagrams/lecture-20/02-progressive-msa.svg`
*Show 5 toy globin sequences with ~30 residues each. Step 1: pairwise alignment of the two closest (sea cucumber vs lamprey). Step 2: profile from step 1 aligned with chicken. Step 3: profile from step 2 aligned with mouse. Step 4: profile from step 3 aligned with human. Final 5-row alignment with conserved positions highlighted (functional residues like the histidine that binds heme).*

---

## Part 3 — Phylogenetic Tree Construction (≈40 min)

### 3.1 What's a phylogenetic tree (≈4 min)

A **phylogenetic tree** is a rooted (or unrooted) bifurcating tree where:

- **Leaves** = observed sequences (taxa).
- **Internal nodes** = inferred ancestors.
- **Edges** = lineages with associated lengths (proportional to evolutionary distance).
- **Root** = most-recent common ancestor (MRCA).

Building such trees from sequence data is the core of evolutionary inference.

### 3.2 The data: distance, character, and likelihood (≈4 min)

Three flavours of phylogenetics:

- **Distance methods**: compute a pairwise distance matrix from the MSA, then tree-build from distances. UPGMA, Neighbor-Joining.
- **Character / parsimony methods**: count substitutions; minimise total changes on the tree. Maximum parsimony.
- **Likelihood methods**: probabilistic substitution model, find tree maximising P(data | tree).

Modern phylogenetics is dominated by likelihood and Bayesian methods (the latter outside scope here).

### 3.3 UPGMA (≈5 min)

**UPGMA (Unweighted Pair Group Method with Arithmetic Mean)**: hierarchical clustering on a distance matrix.

1. Find the closest pair of taxa.
2. Merge them into a new internal node; new branch length = distance / 2.
3. Update the distance matrix.
4. Repeat until one node remains.

**Assumption**: molecular clock — substitution rates are equal across all lineages. When this fails (it usually does), UPGMA produces incorrect topologies. Use it only as a quick first pass.

### 3.4 Neighbor-Joining (≈10 min)

**Neighbor-Joining** (Saitou & Nei 1987): a distance method that doesn't assume the molecular clock.

For each pair $(i, j)$, compute:
$$Q(i, j) = (n - 2) d(i, j) - \sum_k d(i, k) - \sum_k d(j, k)$$

The pair minimising $Q$ is the next pair to join. NJ produces an additive tree if input distances are additive; otherwise it's a robust heuristic.

NJ is fast ($O(n^3)$) and accurate for moderate datasets. Standard for "quick tree from MSA" tasks.

**FIGURE — Figure #3: Neighbor-Joining algorithm step-by-step** → `diagrams/lecture-20/03-nj-tree.svg`
*Six-taxon example. Input: pairwise distance matrix. Iteration 1: compute Q-matrix; identify the smallest Q entry; merge corresponding taxa; update distance matrix. Iterations 2-3: same process on shrinking matrix. Final tree displayed. Annotation: "Q-matrix correction removes molecular-clock assumption — branches can have different lengths."*

### 3.5 Maximum parsimony (≈5 min)

Count the minimum number of substitutions required to explain the data on the tree.

- Brute force: enumerate all $(2n-3)!! / 2^{n-1}$ unrooted topologies; for each, compute parsimony score using Fitch's or Sankoff's algorithm.
- Heuristics: nearest-neighbour interchange, subtree-pruning-and-regrafting.

**Limitation**: long-branch attraction. Two distantly-related sequences with high substitution rates appear closer to each other than to their true relatives. Parsimony is biased; ML and Bayesian methods are unbiased.

### 3.6 Maximum likelihood (≈8 min)

Fix a substitution model (Jukes-Cantor, K2P, F81, HKY, GTR; protein analogs JTT, WAG, LG); then for each candidate tree compute:

$$L(\text{tree}) = P(\text{data} \mid \text{tree}, \text{model})$$

Use Felsenstein's pruning algorithm: dynamic-programming over the tree to compute leaf likelihoods bottom-up.

Search the tree space with hill-climbing (RAxML, IQ-TREE, FastTree). Modern likelihood-tree software produces gold-standard trees in minutes for thousands of sequences.

### 3.7 Branch support: bootstrapping (≈4 min)

How confident are we in each branch? **Bootstrap**: resample MSA columns with replacement → rebuild tree → repeat 100-1000 times → record what fraction of bootstrap trees contain each branch of the original tree. Branches with bootstrap support $\geq 70$ are typically considered reliable.

> **EE framing — phylogenetic inference as graphical-model inference**: A phylogenetic tree with substitution model is exactly a **probabilistic graphical model**. Internal nodes are latent random variables (ancestral sequences); leaves are observed. The substitution model defines transition probabilities $P(\text{child} \mid \text{parent}, t)$ on each edge with edge length $t$. Computing $P(\text{data} \mid \text{tree})$ is **belief propagation** on the tree (Felsenstein's pruning is the bottom-up forward pass). MCMC over trees is **Metropolis-Hastings on the joint tree-and-parameter space**. The whole field is a beautiful instance of probabilistic ML applied long before "graphical models" were named — Felsenstein's 1981 paper predates Pearl's 1988 *Probabilistic Reasoning* by seven years.

---

## Part 4 — Molecular Clock and Dating (≈25 min)

### 4.1 The clock concept (≈4 min)

The **molecular clock hypothesis** (Zuckerkandl & Pauling 1962): substitutions accumulate at a roughly constant rate per unit time. If true, **distance ∝ time**: branch lengths are time elapsed since divergence.

Calibrated against fossil records or known divergence times, branch lengths translate to **absolute dates**.

### 4.2 The clock as a Poisson process (≈6 min)

Substitutions on a lineage are well-modelled as a **Poisson process**: substitution count $X(t) \sim \text{Pois}(\mu t)$ where $\mu$ is the per-site rate. Variance equals mean — overdispersion suggests rate variation across sites.

Per-site rate variation: model as Gamma distribution → site-specific rates → Felsenstein's "+G" correction (gamma-distributed rates, often discretised to 4-8 categories).

### 4.3 Clock violations and the relaxed clock (≈4 min)

The strict clock fails: substitution rates vary across lineages (e.g., rodents > primates), across genes, and across sites. **Relaxed clock** models (BEAST, MrBayes) allow rates to vary along the tree, with rate priors (uncorrelated lognormal, autocorrelated, etc.).

### 4.4 Calibration and dating (≈4 min)

To convert relative branch lengths to absolute times, anchor with **calibration points**: dated fossils or known divergence events. Combined with the relaxed clock, this gives **divergence-time estimation** with credible intervals.

Example: estimating when SARS-CoV-2 emerged from a bat coronavirus ancestor ($\sim$ Oct/Nov 2019, with confidence intervals from 2015 to 2020).

### 4.5 Practical tools (≈4 min)

- **MEGA**: friendly GUI for distance and parsimony trees + simple ML.
- **RAxML / IQ-TREE / FastTree**: production-grade ML.
- **MrBayes / BEAST**: Bayesian + relaxed clock.
- **PhyML**: classic ML.

For a typical paper: align with MAFFT, build tree with IQ-TREE, bootstrap 1000×, visualise with FigTree or iTOL.

### 4.6 Dating worked example (≈3 min)

**FIGURE — Figure #4: Molecular clock dating** → `diagrams/lecture-20/04-clock-dating.svg`
*Phylogenetic tree of 5 mammalian species (human, chimpanzee, mouse, rat, dog). Branch lengths in substitutions-per-site. Right side annotated with calibration points (chimpanzee-human split 6 Mya from fossils). Mapped to absolute time using a relaxed clock. Bottom annotation: "branch lengths are time when the clock holds; the relaxed clock is what we use in practice."*

---

## Part 5 — dN/dS and Selection at Sequence Level (≈25 min)

### 5.1 The codon-level view (≈4 min)

The genetic code is degenerate: most amino acids have multiple codons. **Synonymous substitutions** (codon change without amino-acid change) are evolutionarily near-neutral. **Non-synonymous substitutions** (codon change *with* amino-acid change) are subject to selection.

dS = synonymous substitutions per synonymous site.
dN = non-synonymous substitutions per non-synonymous site.

The **dN / dS ratio** $\omega = dN / dS$ measures selection:

- $\omega = 1$: neutral evolution.
- $\omega < 1$: purifying (negative) selection — most genes most of the time.
- $\omega > 1$: positive (Darwinian) selection — adaptive evolution.

### 5.2 Computing dN and dS (≈6 min)

Algorithms:

- **Nei-Gojobori (1986)**: simple counting; assumes equal mutation rates across codons.
- **Yang & Nielsen (2000)**: maximum-likelihood-based; codon-substitution model with unequal rates.
- **Codeml** (PAML package): standard ML implementation.

Practical: align coding sequences (always nucleotide alignment, codon-aware), run codeml, get per-gene $\omega$.

### 5.3 Detecting selection: the ratio test (≈6 min)

Two models:

- **Null**: $\omega$ is the same across all sites and lineages (clock-like).
- **Alternative**: $\omega$ varies — some sites have $\omega > 1$ (positive selection); others $\omega \ll 1$ (purifying).

Likelihood-ratio test: $2(\log L_{\text{alt}} - \log L_{\text{null}}) \sim \chi^2_{df}$. If significant → some sites under positive selection.

Variants:

- **Branch-site test**: positive selection on a specific branch of the tree.
- **MEME**: episodic positive selection detection (HyPhy).
- **PAML M0/M3/M7/M8 site models**: increasingly refined per-site $\omega$ distributions.

### 5.4 Examples (≈5 min)

**Genes under positive selection** (omega > 1 at specific sites):

- **Influenza HA**: surface protein evolves to evade antibodies. Annual updates to vaccine targets reflect this.
- **MHC**: peptide-binding groove residues are highly variable.
- **Sperm-egg interaction proteins**: Bindin, ZP3 — diverging fast across species.

**Genes under strong purifying selection** ($\omega \ll 1$):

- **Histones**: virtually identical across all eukaryotes.
- **Ribosomal RNA**: similar story.

### 5.5 The deep dive (≈4 min)

> **EE framing — dN/dS as a likelihood ratio test**: dN/dS is the **likelihood-ratio statistic** for the test "does selection deviate from neutral?". The null model has a fixed neutral substitution rate; the alternative allows site-specific deviations. The test statistic $\Lambda = 2 \log (L_{\text{alt}} / L_{\text{null}})$ is asymptotically $\chi^2$-distributed under the null. This is the same Neyman-Pearson hypothesis-testing framework you've seen for variant calling, peak detection, and regulatory-element identification — applied to molecular evolution. The substitution-process estimation is the "sufficient statistic" extraction; the LRT thresholds it.

**FIGURE — Figure #5: dN/dS distribution under selection regimes** → `diagrams/lecture-20/05-dnds.svg`
*Three panels. Panel 1: ω distribution under purifying selection (most sites clustered near 0). Panel 2: ω distribution under positive selection (long tail above 1). Panel 3: example real-world distributions for HA (broad with positive tail), histones (narrow at 0), housekeeping genes (centered around 0.1). Bottom annotation: "ω = dN/dS gives quantitative selection detection at the gene and site level."*

---

## Part 6 — Comparative Genomics and Synteny (≈25 min)

### 6.1 Genome-scale comparison (≈4 min)

Beyond single-gene comparisons, **comparative genomics** asks: how are *whole genomes* related?

Resources:

- **UCSC Genome Browser**: aligned genomes side-by-side.
- **Ensembl Compara**: pairwise genome alignments + synteny + ortholog assignments.
- **OrthoDB / OMA**: pre-computed ortholog clusters across all sequenced species.

### 6.2 Synteny (≈6 min)

**Synteny** = preserved gene order / co-linearity between genomes. Conserved synteny blocks are a hallmark of common ancestry.

Detection: align two genomes; identify regions where homologous genes appear in the same order. Disrupted synteny = chromosomal rearrangement, inversion, translocation since divergence.

**Synteny visualisation**: dot plots (genome 1 on x-axis, genome 2 on y-axis, dots at homolog pairs). Diagonal stripes = collinear blocks. Off-diagonals = rearrangements.

**FIGURE — Figure #6: Synteny dot plot** → `diagrams/lecture-20/06-synteny.svg`
*Two-genome dot plot: human chromosome 17 (x-axis) vs mouse chromosome 11 (y-axis). Diagonal stripes show conserved synteny blocks; off-diagonal points show rearranged regions. Side annotations call out specific gene clusters (e.g., HOXB cluster, conserved across mammals). Top inset: homolog-density per chromosome pair across full human–mouse comparison.*

### 6.3 Orthologs vs paralogs (≈5 min)

- **Orthologs**: genes in different species sharing a common ancestor through speciation. "Same gene in different species."
- **Paralogs**: genes in the same species (or different species) arising from gene duplication. "Different gene from the same family."

Distinguishing them is hard for gene families with both speciation and duplication events.

Tools:

- **RBH (Reciprocal Best Hits)**: simple, ~80% accurate. (Touched in L19 / new L3.)
- **OrthoFinder**: clustering-based; handles complex cases better.
- **OMA / OrthoDB**: pre-computed databases.

### 6.4 Conservation tracks (≈5 min)

UCSC's **PhyloP** and **PhastCons** tracks: per-base scores measuring **evolutionary conservation** across vertebrates (or other clades). High conservation = likely functional sites under purifying selection.

Used for:

- Variant prioritisation in clinical genomics.
- Functional element discovery in non-coding DNA (regulatory regions, miRNA target sites).
- Comparative regulatory genomics.

### 6.5 Comparative genomics for regulatory inference (≈5 min)

If a non-coding region is highly conserved across mammals, it's likely a regulatory element. The technique **phylogenetic footprinting** uses cross-species conservation to identify transcription factor binding sites without ChIP-seq.

Modern workflows combine:

- ChIP-seq (L9 / new L13) for empirical binding.
- Conservation tracks for evolutionary support.
- Motif-scanning (HOMER, FIMO) for sequence specificity.

The intersection is the highest-confidence regulatory inventory.

**EMBED — Artifact #6: Synteny Browser** → `artifacts/lecture-20/06-synteny-browser.html`
*Pick two genomes; artifact shows dot plot + summary stats. Aha: diagonal stripes are collinear blocks; off-diagonal points are rearrangements.*

---

## Part 6.5 — RNA Secondary Structure and Non-Coding RNAs (≈25 min)

### 6.5.1 Why RNA structure matters (≈4 min)

Most of this lecture has assumed sequences are interpreted at the residue level — letters compared, columns aligned, trees built. But RNA molecules **fold** into 3D structures stabilised by base-pairing, and their function depends on that structure as much as on their sequence. Examples:

- **tRNA** — folds into the cloverleaf required for ribosome recognition.
- **rRNA** — the ribosome itself is a structured RNA-protein machine.
- **Riboswitches** — bacterial mRNA elements that change conformation in response to ligand binding, regulating downstream translation.
- **miRNAs** — ~22-nt regulatory RNAs whose precursor stem-loop structure drives Dicer processing.
- **lncRNAs** — long non-coding RNAs (XIST, HOTAIR, NEAT1) function via specific structural domains.

For non-coding RNA family detection, alignment alone is not enough — you need **structure-aware alignment**.

### 6.5.2 RNA folding basics (≈6 min)

Watson-Crick base pairs (A·U, G·C) and the wobble pair (G·U) drive RNA secondary-structure formation. The space of foldings is constrained by:

- **No pseudoknots** in classical RNA folding (most algorithms ignore them; pseudoknots make folding NP-hard).
- **Minimum free energy (MFE)**: the most stable fold under thermodynamic energy parameters (Turner free energies).
- **Pair probabilities**: under the Boltzmann ensemble, each pair has a probability — useful for reliability assessment.

The MFE structure is computed by **Zuker's algorithm** (1989): O(n³) dynamic programming on the sequence with the energy parameters. Implementation: **mfold**, **RNAfold** (ViennaRNA package).

### 6.5.3 ViennaRNA and RNAfold (≈4 min)

**ViennaRNA** (Lorenz et al. 2011) is the standard RNA-folding suite:

- `RNAfold`: minimum-free-energy structure + base-pair probabilities.
- `RNAcofold`: heterodimer folding (e.g., miRNA + target).
- `RNAplfold`: local folding for long sequences.
- `RNAalifold`: consensus structure across an alignment.

Output uses **dot-bracket notation**: `((((....))))` represents a 4-bp stem-loop with 4-nt loop.

### 6.5.4 Profile RNA models — Infernal and Rfam (≈6 min)

For non-coding RNA families, the analog of profile HMMs (Lecture 7 / new L7) is **profile stochastic context-free grammars (SCFGs)**. SCFGs handle base-pairing in the same way HMMs handle linear states — just with the added expressive power of context-free grammars to describe nested structure.

**Infernal** (Eddy 2002, 2013) is the SCFG equivalent of HMMER:

- `cmbuild`: from structure-annotated MSA → covariance model (CM).
- `cmsearch`: scan a sequence database for matches to the CM.
- Sensitive at much greater sequence divergence than profile HMMs because the CM uses structural constraints.

**Rfam** (Kalvari et al. 2021) is the Pfam analog: ~4000 curated families of non-coding RNAs (tRNAs, rRNAs, snoRNAs, miRNAs, lncRNAs, riboswitches), each with seed alignment, consensus structure, and Infernal CM.

For ncRNA discovery in a new genome: run `cmsearch` against Rfam — it's the standard first pass.

### 6.5.5 RNA tertiary structure prediction (≈3 min)

Beyond secondary structure, predicting full 3D coordinates of an RNA is much harder. Tools:

- **RNAComposer**, **3dRNA** — fragment-assembly-based.
- **AlphaFold2-RNA / RoseTTAFold-NA** (2022-2023) — deep-learning extensions.
- **AlphaFold 3** (2024) — handles RNA + DNA + protein + ligand jointly.

Accuracy is good for small structured RNAs (tRNA, riboswitches) but limited for large flexible lncRNAs.

### 6.5.6 The deep dive (≈2 min)

> **EE framing — RNA folding as constrained context-free parsing**: The key complexity in RNA secondary-structure prediction comes from the **nested base-pairing** structure: each pair forms a properly-balanced bracket. Stochastic context-free grammars (SCFGs) handle this exactly via the inside-outside algorithm (CYK-style parsing) — the analog of the forward-backward algorithm for HMMs but extended to handle constituency / nesting. Pseudoknots (overlapping pairs) require context-sensitive grammars and become NP-hard, paralleling the jump from regular to context-sensitive languages in formal grammar hierarchy. The whole field is a clean instance of probabilistic-grammar inference applied biologically.

**FIGURE — Figure #13: RNA secondary structure and SCFG model** → `diagrams/lecture-20/13-rna-structure.svg`
*Top: tRNA cloverleaf with anticodon and acceptor stem highlighted; dot-bracket annotation underneath. Middle: corresponding SCFG production rules in profile-CM format (match-pair, match-singlet, insert, delete states). Bottom: Infernal cmsearch output showing a Rfam family hit with E-value and structure conservation.*

---

## Part 7 — Tools and Practice (≈20 min)

### 7.1 The standard workflow (≈5 min)

For a typical phylogenetics paper:

1. **Collect**: gather homologous sequences (BLAST, OrthoFinder, manual curation).
2. **Align**: MAFFT L-INS-i (or MUSCLE).
3. **Trim**: Gblocks or trimAl removes poorly-aligned columns.
4. **Tree**: IQ-TREE with model selection + 1000 bootstraps.
5. **Visualise**: iTOL or FigTree.
6. **Selection**: PAML / HyPhy if dN/dS is the question.
7. **Comparative**: UCSC Browser if genome context matters.

### 7.2 Common pitfalls (≈4 min)

- **Misaligned columns**: a few badly-aligned positions will throw off the tree. Trim aggressively.
- **Long-branch attraction**: distantly-related fast-evolving sequences appear close. Add intermediates if possible; use ML.
- **Saturated distances**: at very high divergence, multiple substitutions per site cause distance underestimation. Use codon-aware methods for protein-coding alignments.
- **Missing data**: indels and gaps must be handled correctly (either removed or modelled).
- **Wrong root**: for unrooted trees, picking a root requires an outgroup; placing the root wrongly inverts ancestor relationships.

### 7.3 The deep dive: tree of life (≈3 min)

The **Tree of Life** (Hug et al. 2016, *Nature Microbiology*) is a ~3000-tip ML tree built from 16S rRNA + concatenated marker genes. It re-classified the bacteria-archaea-eukarya split, with archaea now placed as sister to eukaryotes (the "two domains" hypothesis). Comparative genomics underlies this entire reclassification.

### 7.4 Getting hands-on (≈4 min)

**Quick self-paced:**

- Run MAFFT on a small protein family from UniProt (use the BLAST exercise from L19 to seed it).
- Build a tree with FastTree (one-line UNIX command).
- Run codeml on a coding alignment for $\omega$.
- Visualise the result on iTOL.

You can do all of this in an afternoon with no installation (NCBI's Web Phylogeny Toolkit covers the basics).

### 7.5 The 2024 frontier (≈4 min)

Recent developments at the frontier:

- **AlphaFold MSAs**: AlphaFold uses MSAs as input; the 200M-structure AlphaFold Database is a comparative-genomics product.
- **PhyloPGM**: deep learning on phylogenetic graph embeddings.
- **Phylogenetic foundation models**: pre-trained on 100k+ trees for transfer to novel inference tasks (e.g., Zhukov et al. 2024).
- **Bayesian phylogenetics at scale**: efficient MCMC with GPU acceleration (e.g., BEAGLE library).

These are mostly research-frontier; the bread-and-butter pipeline (MAFFT → IQ-TREE → bootstrap) is unchanged.

**FIGURE — Figure #7: Standard phylogenetics workflow** → `diagrams/lecture-20/07-workflow.svg`
*Top-to-bottom flowchart: collect homologs → MAFFT alignment → trimAl trimming → IQ-TREE model selection + ML tree + bootstrap → visualise + downstream selection / dating. Side panel: alternative branches (BLAST seeds, OrthoFinder ortholog input, PAML for selection, BEAST for dating).*

---

## Wrap-up (≈10 min)

### What you should take away

- **MSA is intractable in exact form**; progressive alignment + iterative refinement (Clustal, MUSCLE, MAFFT) is the practical compromise.
- **Phylogenetic trees** can be inferred from distances (NJ), characters (parsimony), or likelihood (ML). ML is the gold standard.
- **The molecular clock** is a Poisson-process model; relaxed-clock variants handle real-world rate variation.
- **dN/dS** is a likelihood-ratio test for selection; $\omega > 1$ = positive selection, $\omega < 1$ = purifying.
- **Comparative genomics** uses synteny, orthologs, and conservation tracks to align genome organisation across species.
- **EE framings**: trees as graphical models; MSA as constrained tensor DP; molecular clock as a Poisson process; dN/dS as a likelihood-ratio test.

### Next lecture

Genome assembly (existing L3, becomes new L5): from reads to contigs to scaffolds. The MSA you just learned is what's used to validate assembly correctness via cross-sample comparison.

### Homework

1. Take 5 protein sequences from the same family (e.g., haemoglobin orthologs). Run MAFFT L-INS-i. Identify conserved positions. Trim to a 50-residue core domain.
2. Build NJ and ML trees on the same dataset; compare topologies. Bootstrap 100×; report which branches are well-supported.
3. Compute dN/dS for the alignment from #1 using codeml. Test the M0 vs M2a (positive selection) models. Report ω and p-value.
4. From UCSC Browser, visualise the human-mouse synteny for a 1 Mb region of human chromosome 17. Identify a clear rearrangement; describe its inferred type (inversion, translocation).
5. Pick a recently-published phylogeny (any 2023+ Nature paper). Recreate Figure 1 from its data. Document deviations.

### Recommended reading

- Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Research* 32, 1792–1797.
- Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution* 30, 772–780.
- Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4, 406–425.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution* 17, 368–376.
- Yang, Z. (2007). PAML 4: phylogenetic analysis by maximum likelihood. *Molecular Biology and Evolution* 24, 1586–1591.
- Hug, L. A., Baker, B. J., Anantharaman, K., et al. (2016). A new view of the tree of life. *Nature Microbiology* 1, 16048.
- Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Molecular Biology and Evolution* 32, 268–274.
- UCSC Genome Browser: <https://genome.ucsc.edu/>
- Ensembl Compara: <https://www.ensembl.org/info/genome/compara/index.html>
- iTOL: <https://itol.embl.de/>
- IQ-TREE: <http://www.iqtree.org/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — From Pairwise to Multiple                  | 25 min | 0:25 |
| Part 2 — Progressive Alignment in Practice           | 40 min | 1:05 |
| Part 3 — Phylogenetic Tree Construction              | 40 min | 1:45 |
| Part 4 — Molecular Clock and Dating                  | 25 min | 2:10 |
| Part 5 — dN/dS and Selection at Sequence Level       | 25 min | 2:35 |
| Part 6 — Comparative Genomics and Synteny            | 25 min | 3:00 |
| Part 6.5 — RNA Secondary Structure and ncRNAs         | 25 min | 3:25 |
| Part 7 — Tools and Practice                           | 20 min | 3:45 |
| Wrap-up                                                | 10 min | 3:55 |

**Total:** ~3h 55min of content.
