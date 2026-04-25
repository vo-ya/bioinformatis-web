# Lecture 15 — Protein Structure Prediction (AlphaFold-era)

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-15.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Describe the four levels of protein structure (primary, secondary, tertiary, quaternary) and explain why the tertiary fold is the principal predictive target.
2. Summarise the classical era of structure prediction (homology modelling, threading, CASP) and explain why it plateaued prior to 2018.
3. Explain how multiple-sequence alignment (MSA) covariation encodes 3D contact information; describe DCA / EVfold as inverse-covariance estimators.
4. Walk through the AlphaFold2 architecture: Evoformer (axial attention over MSA + pair representation), structure module (invariant point attention, equivariant geometry), recycling.
5. Interpret AlphaFold's confidence metrics (pLDDT, PAE) and describe failure modes (orphan proteins, disorder, multi-state conformations).
6. Describe AlphaFold3's diffusion-based generation, multi-chain / ligand / nucleic-acid handling, and how single-sequence methods (ESMFold, ESM-3) differ.
7. Explain inverse folding (ProteinMPNN) and diffusion-based protein design (RFDiffusion) as inverse problems with structural priors, and place them in the de-novo design pipeline.
8. Describe the AlphaFold Database (~200M predicted structures) and what the field's practical workflow has become post-2021.

---

## Part 1 — Protein Structure Basics (≈25 min)

### 1.1 Why structure matters (≈5 min)

Lectures 1–14 have been about sequence: DNA, RNA, reads, variants, populations. Proteins are the functional molecules encoded by genes, and for almost every question about protein function — what does this enzyme do, how does this drug bind, why does this mutation cause disease — **3D structure** is the relevant representation, not sequence.

Proteins are polymers of amino acids. Twenty standard amino acids (plus a few rare exceptions) make up the chemical alphabet. A typical protein is 50–2000 residues long. The functional behaviour — catalysis, binding, signalling, transport — emerges from the three-dimensional arrangement of those residues in space.

A single residue by itself is inert. The chain folds, residues come into contact, those contacts form binding pockets, active sites, and interfaces. **Structure determines function.** Without structure, sequence is a flat book — readable but not interpretable.

**FIGURE — Figure #1: Four levels of protein structure** → `diagrams/lecture-15/01-structure-hierarchy.svg`
*A four-panel stack. Panel 1: primary structure — a linear chain of amino-acid letters (N-M-Q-V-L-...). Panel 2: secondary structure — the same chain with α-helix spirals and β-strand arrows highlighted. Panel 3: tertiary structure — the same chain folded into a compact 3D globule, helices and strands in context. Panel 4: quaternary structure — four subunits of haemoglobin assembled with the haem groups bound. Annotations: "primary = sequence", "secondary = local H-bond patterns", "tertiary = global 3D fold", "quaternary = multi-chain assembly".*

### 1.2 Primary, secondary, tertiary, quaternary (≈8 min)

**Primary structure**: the linear amino-acid sequence, N-terminus to C-terminus. A string in a 20-letter alphabet (plus two rare: selenocysteine U, pyrrolysine O). Encoded directly by the gene: codon → amino acid.

**Secondary structure**: local backbone conformation, driven by hydrogen bonds between backbone carbonyl (C=O) and amide (N-H) groups. Two dominant motifs:

- **α-helix**: right-handed spiral, 3.6 residues per turn. Backbone H-bonds within the helix. Stable, compact, common.
- **β-strand**: extended conformation; H-bonds to another strand (parallel or antiparallel) to form a **β-sheet**.

A typical globular protein is 30–50% helix + strand, the rest loops and turns. Secondary structure is predictable from primary sequence at ~85% accuracy (classical: PSIPRED; modern: essentially perfect from AlphaFold's internal representations).

**Tertiary structure**: the full 3D arrangement of a single chain. This is what "solving a structure" usually means. Parameterised by the backbone torsion angles (φ, ψ, ω) plus side-chain torsions (χ₁, χ₂, …). Given these, atomic coordinates are fully determined.

Representations:

- **Cartesian coordinates** (x, y, z per atom). ~8 × N atoms per residue. Standard PDB format.
- **Internal coordinates** (torsions + bond lengths + angles). More compact; AlphaFold outputs these and converts to Cartesian.
- **Contact map**: a binary / continuous N × N matrix where entry (i, j) = distance between residues i and j. Dominant representation for contact-prediction methods.

**Quaternary structure**: assembly of multiple chains (subunits). Haemoglobin is 4 chains; the ribosome is 80 chains + RNA. Predicting quaternary requires predicting both the individual chain structures and their interfaces.

### 1.3 The Anfinsen experiment (≈6 min)

**Christian Anfinsen (1961, Nobel Prize 1972)** showed that a small protein (RNase A, 124 residues) could be completely denatured (unfolded), then refolded in vitro back to its native, fully-active structure, without any help. The instructions for folding are entirely in the sequence.

The implication: **the folded structure is a minimum of a free-energy landscape over the conformational space**, and sequence uniquely determines that minimum. Computationally, this says the structure-prediction problem is well-defined: given a sequence, find the energy-minimising 3D arrangement.

In practice, finding that minimum is vastly harder than the Anfinsen formulation suggests. The conformational space is combinatorially enormous (Levinthal's paradox: if a protein sampled conformations randomly, it would take longer than the age of the universe to find the native state). Real proteins fold on millisecond-to-second timescales because the landscape is funnelled toward the native state — not sampled randomly.

> **Historical pointer**: Anfinsen's RNase A experiment used chemical denaturation with urea and reduction of disulfide bonds with β-mercaptoethanol. After removing both the denaturant and the reducing agent, the enzyme refolded spontaneously, with full catalytic activity restored. The paper that won the Nobel was Anfinsen (1973) *Science* 181: 223–230, "Principles that govern the folding of protein chains". The result established that sequence → structure is a well-defined computational problem; the 50 years from that paper until AlphaFold2 were essentially an attempt to find the algorithm.

### 1.4 Domains, folds, and the fold-space universe (≈6 min)

Most proteins are built from **domains** — semi-independent structural units that can appear in different proteins as mix-and-match modules. A domain is typically 50–250 residues and folds independently.

**Folds** are repeated architectural patterns: specific arrangements of secondary-structure elements. Classical fold classifications:

- **SCOP / SCOPe** (Structural Classification Of Proteins). Hierarchical: Class → Fold → Superfamily → Family.
- **CATH** (Class, Architecture, Topology, Homology). Similar hierarchy, different methodology.

Estimated number of distinct folds in nature: ~1000–2000. Before AlphaFold, most of them had been observed experimentally; new folds were rare. After AlphaFold applied to metagenomic sequences, a modest number of new folds have been identified.

**Domain composition is modular**: eukaryotic proteins frequently have 2–6 domains, each contributing an independently-foldable piece. A typical kinase has a catalytic domain + regulatory domain + targeting domain. Predicting the whole protein = predicting each domain + predicting their spatial relationship.

> **Intuition box**: Think of proteins like Lego structures assembled from standard brick shapes (folds) connected by flexible joints (linkers). The fold universe is the catalogue of standard brick shapes nature has evolved; there aren't infinitely many, because backbone physics + residue chemistry constrain what can stably fold. Domain architecture = which bricks this protein uses, in what order. AlphaFold's job is first to predict each brick's shape from its sequence, then to predict how the bricks sit relative to each other. The first sub-problem is what makes prediction feasible; the second is where most remaining error lives.

---

## Part 2 — The Classical Era (≈20 min)

### 2.1 Homology modelling (≈6 min)

For proteins with a sequence-similar experimental structure ("template") in the PDB, **homology modelling** is straightforward:

1. **Find a template**: BLAST / HMMER search against PDB sequences.
2. **Align**: compute the sequence alignment between query and template.
3. **Copy backbone**: for aligned residues, copy the template's backbone coordinates.
4. **Model side chains**: place the query's side chains using a rotamer library (SCWRL, Rosetta side-chain packer).
5. **Refine**: energy-minimise the result.

Tools: **Modeller** (Sali 1993+), **SWISS-MODEL** (Schwede et al. 2003+). SWISS-MODEL runs the full pipeline as a web server.

Quality depends on sequence identity to template:

- **>50% identity**: model is "close" (backbone RMSD ~1 Å), usable for most purposes.
- **30–50% identity**: "twilight zone"; models have structural errors but functional inference often still works.
- **<30%**: unreliable; structural details can be completely wrong.

Classical limit: for proteins without a close homologue in the PDB (at least 25–30% identity), homology modelling fails. For orphan proteins (found only in a single species), it's unusable.

### 2.2 Threading and ab initio (≈7 min)

**Threading** (Bowie, Lüthy &amp; Eisenberg 1991; **I-TASSER** = Zhang 2008+) generalises homology modelling: for each known fold in the PDB, score how well the query sequence "fits" the fold, using sequence-structure compatibility metrics. Pick the best-scoring fold, build the model.

Good when the query shares a fold with a known protein but diverges at the sequence level (analogous folds with different sequences — "convergent evolution of fold" or "remote homology").

**Ab initio** ("from the beginning"): predict structure without any template. Classical approaches:

- **Fragment assembly** (**Rosetta**, Baker 1997+). Represent short (9-residue) local fragments from the PDB as structural building blocks. Assemble them with Monte Carlo sampling; score with a learned energy function. Rosetta's fragment library is the canonical example.
- **Molecular dynamics** (Duan &amp; Kollman 1998). Simulate the folding trajectory atom-by-atom with a physics-based force field. In 1998 this worked only for very small proteins (villin headpiece, 35 residues) after months of compute.

Both approaches produced moderate success for small proteins but plateaued on anything non-trivial.

### 2.3 CASP: the community benchmark (≈5 min)

**CASP** (Critical Assessment of Structure Prediction) is the biennial blind-prediction experiment, run since 1994 by Moult et al. Participants receive sequences of soon-to-be-released experimental structures; they submit predictions; after the experimental structures are released, predictions are scored.

CASP tracks:

- **Template-based modelling (TBM)**: homology and threading.
- **Free modelling (FM)**: no homologue; ab initio only.
- **Quaternary (multimers)**, **refinement**, **contact prediction**, and others.

Primary metric: **GDT_TS** (Global Distance Test, Total Score). Scores ≤ 100; roughly, 50 is "recognisable fold", 70 is "good prediction", 90+ is "experimental quality".

History, abbreviated:

- **CASP1 (1994)–CASP8 (2008)**: slow, incremental progress. Homology models gradually improve; FM plateau.
- **CASP9–CASP12 (2010–2016)**: Rosetta dominates FM; Zhang-lab's I-TASSER leads TBM. FM best scores ~30–40 GDT_TS.
- **CASP13 (2018)**: DeepMind's **AlphaFold1** enters. FM scores jump to ~60 GDT_TS — an unprecedented single-edition leap, attributable to contact prediction (Part 3).
- **CASP14 (2020)**: **AlphaFold2** reaches ~90 GDT_TS on FM targets — effectively experimental quality. The field pivots.
- **CASP15–CASP16 (2022–2024)**: AlphaFold2 and successors dominate; methods compete within a narrower band. New targets: protein-protein complexes, protein-ligand docking, disorder.

**FIGURE — Figure #2: CASP GDT_TS progression** → `diagrams/lecture-15/02-casp-progression.svg`
*A timeline plot from CASP1 (1994) to CASP16 (2024). X-axis: CASP edition. Y-axis: median best GDT_TS on Free-Modelling targets. Dots and a line showing slow progress 1994–2016, a jump at CASP13 (AlphaFold1), a massive jump at CASP14 (AlphaFold2 at ~90), subsequent editions staying at the high plateau. Annotations point out Rosetta era, contact-prediction era, AlphaFold2 inflection.*

### 2.4 Why the classical era plateaued (≈2 min)

Two structural reasons the classical era couldn't break through:

- **Template-dependence**: without a PDB homologue, there was no structural signal to transfer.
- **Energy functions were imperfect**: Rosetta's energy function is an approximation to the real physics; its minimum doesn't always correspond to the native state. Improving the energy function incrementally didn't close the gap.

The breakthrough came from elsewhere: **the information was already in the sequences**, if you knew how to read it across many sequences at once.

---

## Part 3 — MSA-Based Methods and Coevolution (≈30 min)

### 3.1 Why coevolving residues are spatially close (≈8 min)

Given a protein family (e.g. hundreds of homologues of a kinase from different species), align them into a **multiple sequence alignment (MSA)**. Each column is one residue position; each row is one sequence.

Observation: some pairs of columns vary together. If column *i* mutates from L to V, column *j* frequently also mutates (in the same sequences). These residues are **coevolving**.

Why? The physical explanation: residues *i* and *j* are in 3D **contact** in the folded structure. A mutation at one changes the local chemistry; the contacting residue must co-mutate to maintain stability. Evolution across many species sweeps many such double-mutations across the family.

Formalised: the 3D contact structure leaves an imprint on the evolutionary record. Reading that imprint — which column pairs co-mutate — recovers the contact map.

**FIGURE — Figure #3: Coevolution and 3D contacts** → `diagrams/lecture-15/03-coevolution.svg`
*Top: a schematic MSA — ~15 sequences rows, ~40 columns. Two columns highlighted (i=12, j=34). In rows where column 12 is L, column 34 is V; in rows where 12 is I, 34 is T. Dashed line connecting the two columns, labelled "coevolving pair". Middle: the folded structure of this protein with residues 12 and 34 highlighted as neighbours in 3D, shown as contact. Bottom: inference arrow: "co-mutation pattern in MSA → 3D contact".*

### 3.2 Pairwise mutual information and its flaws (≈6 min)

Naive approach: for each column pair (i, j), compute **mutual information (MI)** between the two columns treated as categorical distributions. High MI = the columns carry information about each other.

$$\text{MI}(i, j) = \sum_{a, b} P(a, b) \log \frac{P(a, b)}{P(a) P(b)}$$

where $P(a)$ is the marginal probability of amino acid $a$ at position $i$, and $P(a, b)$ is the joint probability.

Problem: **indirect correlations**. If columns i–j and j–k are both in contact, then columns i–k will also have high MI even if they aren't in direct contact. The signal is spread along the chain of correlations.

Mutual information alone produces spurious "contacts" between residues that are merely both correlated with a common third residue. Classical MI-based contact prediction was unreliable because of this.

### 3.3 Direct Coupling Analysis (DCA) (≈8 min)

The insight that unlocked classical coevolution methods: **inverse covariance**. The partial correlation between two variables, controlling for all others, isolates direct statistical dependencies.

Define a Potts-model-like probability distribution over MSA sequences:

$$P(\sigma) \propto \exp\left( \sum_i h_i(\sigma_i) + \sum_{i < j} J_{ij}(\sigma_i, \sigma_j) \right)$$

where $\sigma = (\sigma_1, \ldots, \sigma_L)$ is a sequence, $h_i$ are per-position biases (conservation), and $J_{ij}$ are pairwise couplings. Pairs with large $|J_{ij}|$ are **directly** coupled — indirect correlations have been "explained away" by the model.

Two classical algorithms to estimate $J_{ij}$:

- **mfDCA** (mean-field DCA; Morcos et al. 2011). Close-form approximation. Fast.
- **plmDCA** (pseudolikelihood; Ekeberg et al. 2013). More accurate; the standard until the deep-learning era.

Given estimated $J_{ij}$, the **contact score** for pair (i, j) is a norm over the 20 × 20 coupling matrix:

$$C_{ij} = \| J_{ij} \|_F = \sqrt{\sum_{a, b} J_{ij}(a, b)^2}$$

Large $C_{ij}$ = predicted contact.

**EVfold** (Marks et al. 2011) was the first well-known DCA-based contact predictor. On proteins with large MSAs (>1000 diverse sequences), EVfold's top-scored contacts were ~80% correct — far better than MI-based methods and enough to seed ab-initio folding (Rosetta with DCA contacts as restraints achieved solid FM predictions at CASP10–12).

> **EE framing — MSA covariation as inverse covariance estimation**: Reading co-mutations from an MSA is exactly the **graphical-model / inverse-covariance-estimation** problem. The MSA is a matrix of categorical samples from a joint distribution; the true contact graph is the inverse-covariance sparsity pattern. Naive MI = covariance (direct + indirect). DCA = inverse covariance (direct only). Graphical lasso (Friedman et al. 2008) is the Gaussian version of the same problem; plmDCA is the discrete-variable version. The mathematical structure is identical, and the field rediscovered it — the DCA papers cite statistical-physics literature (Ising models, inverse statistical mechanics) rather than ML statistics, but the problems are the same.

**FIGURE — Figure #4: DCA contact prediction** → `diagrams/lecture-15/04-dca-contacts.svg`
*Three-panel workflow. Left: raw MI matrix for an MSA — noisy, many off-diagonal blobs. Middle: the DCA / plmDCA inverse-covariance matrix — most of the noise suppressed, diagonal band plus isolated sharp contacts. Right: the native contact map from the experimental PDB structure — binary, showing true contacts. Annotation: "DCA recovers the sparse true contacts by inverting the covariance; MI alone mixes direct and indirect correlations".*

### 3.4 MSA depth matters (≈4 min)

DCA's accuracy depends critically on **MSA depth** — how many diverse sequences are in the family.

- **< 50 sequences**: DCA unreliable; estimation variance too high.
- **100–1000 sequences**: moderate accuracy.
- **> 1000 sequences**: excellent, top-100 contacts ~80% correct.
- **> 10000 sequences** (metagenomic MSAs for abundant families): near-perfect.

Orphan proteins — those with no detectable homologues — are invisible to DCA. About 20% of human proteins are in this category or close to it, and DCA + Rosetta was uselessly bad on them.

This is a key limit that AlphaFold2 partially overcomes: shallower MSAs (10–50 sequences) still produce usable predictions via deep-learning-based contact refinement.

### 3.5 From contacts to full structures (≈4 min)

A predicted contact map doesn't give you coordinates directly. To go from contacts to 3D structure:

1. Treat predicted contacts as soft constraints: $d_{ij} < 8$ Å.
2. Run fragment assembly (Rosetta) or distance-geometry embedding (CNS) subject to these constraints.
3. Pick the best-scoring resulting model.

DCA + Rosetta was the CASP12 state-of-the-art. It worked well for proteins with deep MSAs; it failed on shallow MSAs or complex topologies.

AlphaFold1 (CASP13) replaced the Rosetta step with a neural-network-based distance-distribution predictor, plus a gradient-descent folder. This was incremental but significant. AlphaFold2 threw out the two-step architecture entirely in favour of end-to-end learning (Part 4).

**EMBED — Artifact #1: MSA Coevolution and Contact Prediction** → `artifacts/lecture-15/01-coevolution-contacts.html`
*Build a toy protein family (user sets contact topology: 1D chain, 2D sheet, α-helical bundle). Simulate an MSA by running a Metropolis chain on the Potts model with those contacts. Compute MI and DCA contact scores. Plot against ground truth. Target aha: DCA recovers the true contacts even when MI is hopelessly tangled — inverse covariance is the right estimator.*

---

## Part 4 — AlphaFold2 Architecture (≈60 min)

This is the single largest content block in the lecture. AlphaFold2 transformed structure prediction; its architecture is also an EE-readable tour of modern ML engineering for biological data.

### 4.1 The input representation (≈7 min)

AlphaFold2 (Jumper et al. 2021, *Nature*) takes a protein sequence and produces a 3D structure plus per-residue confidence. Its input is not just the sequence; it's the sequence plus **its MSA** (computed via search against large sequence databases: UniRef, MGnify, BFD).

Two data objects are maintained throughout the network:

- **MSA representation** $\mathbf{M} \in \mathbb{R}^{s \times r \times c}$ — $s$ sequences in the MSA × $r$ residues × $c$ channels (initially 256, later expanded).
- **Pair representation** $\mathbf{Z} \in \mathbb{R}^{r \times r \times c}$ — $r \times r$ residue-pair matrix × $c$ channels. Think of this as a "learned contact map with rich features per pair".

The network's job is to iteratively refine both representations, with each updating the other, until they encode enough geometric information for the structure module to produce coordinates.

### 4.2 Overall architecture (≈8 min)

Three major components:

1. **Evoformer** — ~48 blocks of axial attention and triangle updates that refine (MSA, pair) jointly. Most of the compute. Section 4.3.
2. **Structure module** — 8 iterations of invariant point attention (IPA) plus a parameterised geometry updater that produces 3D coordinates. Section 4.5.
3. **Recycling** — the output (updated pair repr + structure) is fed back to Evoformer as input; repeat 3–4 times. Section 4.6.

Plus: confidence heads (pLDDT, PAE, pTM) computed from the final representations; template embeddings if PDB templates exist.

**FIGURE — Figure #5: AlphaFold2 architecture overview** → `diagrams/lecture-15/05-af2-overview.svg`
*Top-to-bottom flow. Input: protein sequence + MSA + (optional) templates. Box 1: Evoformer (48 blocks), maintaining MSA + pair representations. Box 2: Structure module (8 iterations), consuming final pair representation + single representation, producing coordinates. Feedback arrow (recycling) from structure + pair back to Evoformer input. Output: 3D structure with per-residue pLDDT confidence and N×N PAE matrix. Parameter-count annotation: "~93M trainable parameters".*

### 4.3 The Evoformer (≈15 min)

The Evoformer is 48 repeated blocks, each updating the MSA and pair representations. Each block has ~six sub-operations:

1. **Row-wise MSA attention with pair bias**: each row (one MSA sequence) attended to itself with residue-position as the sequence axis. The pair representation $\mathbf{Z}$ biases the attention map. This is how pair information propagates into the MSA.
2. **Column-wise MSA attention**: each column (one residue position across species) attended to itself along the species axis. Different from row attention: informs each residue about its homologous versions.
3. **MSA transition** (feedforward).
4. **Outer product mean → pair update**: the MSA is summarised into an outer-product-mean statistic over species and fed into the pair representation. This propagates MSA-derived information into the pair.
5. **Triangle multiplicative update**: the pair representation updates itself using triangle geometry — if i–k and j–k are both "contact-like", then i–j should be "contact-like" too. Exploits the triangle inequality implicitly.
6. **Triangle attention**: attention-based version of the above. Each pair (i, j) attends over all pairs sharing a residue (i, k) or (k, j).

The architecture is carefully designed around the pair representation having a **geometric interpretation** — pair-updates respect the fact that (i, j) + (j, k) → (i, k) should be consistent.

**Axial attention**: rather than full O(N²) attention over the MSA (which would be O(s² r²) — infeasible at s = 512, r = 512), row attention operates over one axis (residues, fixed sequence), column attention over the other (species, fixed residue). Each is O(s × r²) or O(r × s²); the combination is quadratic in one axis at a time, not their product. This is the "axial attention" pattern popularised in image transformers (Ho et al. 2019) and adapted here.

> **EE framing — Evoformer as axial attention on MSA + pair**: The Evoformer is a two-axis transformer. One axis is residues (tokens along the chain); the other is species (tokens down the MSA). Full 2D attention would be prohibitive for large MSA × long sequence combinations. Axial attention alternates row-attention and column-attention to approximate full 2D attention at quadratic cost in one axis, not the product. This is the same design used in image transformers (Axial-DeepLab, Ho 2019), adapted to the MSA × residue 2D grid. Triangle updates in the pair representation are analogous to a **graph neural network** layer where edges update themselves based on two-step neighbours — imposing a geometric prior that pair features should respect the triangle inequality.

**FIGURE — Figure #6: Evoformer block** → `diagrams/lecture-15/06-evoformer-block.svg`
*A single Evoformer block as a DAG of operations. Nodes: MSA repr (top-left), pair repr (top-right). Arrows: row-wise MSA attention (with pair bias), column-wise MSA attention, outer product mean (MSA → pair), triangle multiplicative update, triangle self-attention. Each operation labelled with complexity (e.g. O(s × r²) for row attention). Output: updated MSA + updated pair. Annotation at right: "block repeats 48 times".*

### 4.4 Why does the Evoformer work? (≈4 min)

An intuitive reading: the Evoformer's main job is to build a high-quality pair representation whose learned attention patterns encode not just "who contacts whom" but **distances and orientations** — enough geometric information for a downstream 3D embedder.

The MSA is the input signal: millions of evolutionary examples of which mutations co-occur. The pair representation is the output signal: a learned, richly-featured N × N matrix where pair (i, j) encodes "what's the 3D relationship between residues i and j". The Evoformer's 48 blocks iteratively refine this.

Why 48 blocks? Empirical. Fewer → worse. More → diminishing returns, higher memory. The depth was tuned for CASP14 performance.

> **Intuition box**: The MSA is a field of evolutionary gossip about which residue pairs are in contact. Raw MSA = all the gossip mixed together. The Evoformer's job is to distill that gossip into a clean pair representation where each (i, j) entry has a well-defined answer to "are you two residues in contact, at what distance, in what orientation?" The structure module then reads off 3D coordinates from the pair representation's answer. Each Evoformer block is a round of cross-referencing — MSA informs pair, pair informs MSA, triangles update pairs — after 48 rounds, the answer is clean enough to produce coordinates.

### 4.5 The structure module + IPA (≈14 min)

After the Evoformer produces a final (MSA, pair) representation, the **structure module** converts them into 3D coordinates.

Architecture: 8 iterations of the same block. At each iteration:

1. Take current backbone frames $\mathbf{T}_i \in \text{SE}(3)$ for each residue (an $\text{SE}(3)$ frame = 3D position + orientation).
2. Run **Invariant Point Attention (IPA)** — attention where keys / queries / values live in 3D space relative to each residue's frame.
3. Update frames: each residue gets translated and rotated based on the IPA output.
4. Compute per-residue "points" (α-carbon and side-chain torsions).

The structure is **iteratively refined**: at each iteration, every residue knows its current 3D position and the 3D positions of all others (via IPA); it then updates its position based on local geometry signals.

**Invariant Point Attention (IPA)**:

- Attention keys and values are **points** expressed in each residue's own 3D frame.
- The attention weight between residue *i* and residue *j* depends on the distance between *i*'s query-points (expressed in its frame) and *j*'s key-points (expressed in *j*'s frame).
- When computing the distance, both points are re-expressed in a common frame — the attention is **equivariant** to global rotations and translations of the entire structure.

Equivariance means: if you rotate the input by any $R \in \text{SO}(3)$, IPA's outputs rotate by the same $R$. This respects the physics: the protein's energy, stability, and contacts don't depend on its orientation in the lab.

> **EE framing — IPA as SE(3)-equivariant networks**: The structure module's IPA is an **$\text{SE}(3)$-equivariant operation** — it commutes with the group of 3D rotations and translations applied to the input. The design approach is the same as in e3nn / NequIP (Geiger et al. 2022) and in SE(3)-Transformers (Fuchs et al. 2020): treat per-residue features as tensors with well-defined rotation behaviour (scalars, vectors, higher-order), and construct operations that preserve the rotation class. The structure module's output — 3D coordinates — must be equivariant under input rotation, so its computation must be too. Without equivariance, training would need to learn rotational invariance from data (data augmentation by rotation), which is expensive and incomplete.

Frames live in $\text{SE}(3) = \text{SO}(3) \ltimes \mathbb{R}^3$. Each residue's frame is parameterised by a 3D rotation + translation. Updates to the frame are composed via quaternion multiplication + translation addition.

**FIGURE — Figure #7: Structure module + IPA** → `diagrams/lecture-15/07-structure-module.svg`
*Diagram showing a cartoon of the structure module's iterative refinement. Left: initial residue positions (all at origin, random orientations). Middle: after IPA iteration 1 — residues moved to approximate positions. Right: after IPA iteration 8 — fully-folded structure. Below: a detail of IPA — query points in residue i's frame; key points in residue j's frame; both mapped to a common reference for distance computation; attention weight. Annotation: "SE(3)-equivariant by construction — rotate input, output rotates the same way".*

### 4.6 Recycling and confidence (≈8 min)

**Recycling**: after the structure module produces coordinates, feed the final pair representation **and** the current predicted structure back into the Evoformer as an input, then rerun. Typically 3–4 recycling iterations.

Why: Evoformer can be too aggressive in an early pass, reaching a weird intermediate state. Recycling lets the network iteratively correct, with each recycling cycle benefiting from the structure information produced in the previous.

Training was particularly careful about recycling — early-layer outputs during training used 0–3 recycles with randomised stochasticity, so the model generalises to any recycling count at inference time.

**Confidence heads**:

- **pLDDT** (predicted local distance difference test). Per-residue, 0–100. pLDDT = 100 means "I'm very confident this residue's position is correct within 1 Å"; 50 means "disordered or uncertain"; below 50 = probably wrong / disordered.
- **PAE** (Predicted Aligned Error). Per-pair, in Angstroms. PAE$_{ij}$ = expected error in residue *j*'s position if the structure were aligned to residue *i*'s frame. Low PAE = "j's position is well-defined relative to i". Large PAE = "I don't know where *j* is relative to *i*".
- **pTM** (predicted TM-score). A scalar, 0–1, summarising overall predicted correctness.

PAE is particularly useful for interpreting multi-domain proteins: residues within a domain have low inter-PAE (domain structure confident); residues across domains often have high inter-PAE (domain orientations uncertain).

> **Warning box**: pLDDT is **calibrated** in the sense that "pLDDT = 80" means "on average, these residues are correct within a few Å." But it measures **local** accuracy only. A protein where every residue is pLDDT = 90 can still have the two halves flipped if the inter-domain PAE is high. Always look at both. AlphaFold is not overconfident in local structure; it can be overconfident in domain arrangement.

**FIGURE — Figure #8: pLDDT and PAE confidence metrics** → `diagrams/lecture-15/08-plddt-pae.svg`
*Three-panel figure. Left: a cartoon of a two-domain protein, rainbow-coloured by pLDDT — most residues high pLDDT (blue) except a flexible linker (red/orange). Middle: the PAE matrix for the same protein — low PAE within each domain (two dark blocks on the diagonal), high PAE between domains (off-diagonal bright regions). Right: example interpretation text. "Individual domains well-folded; inter-domain orientation uncertain; use PAE for functional inference".*

### 4.7 Training regime (≈4 min)

AlphaFold2 was trained on:

- **PDB** (~170k experimental structures at training time).
- **UniRef90 + BFD** (hundreds of millions of sequences for MSA generation).
- **Self-distillation**: predict structures for UniRef sequences, use high-confidence predictions as additional training data. This dramatically expanded the effective training set.

Training ran for ~11 days on 128 TPUv3 cores. Inference on a single protein ~1–5 minutes on a single GPU (depending on sequence length).

**EMBED — Artifact #2: Evoformer Axial-Attention Visualiser** → `artifacts/lecture-15/02-evoformer-attention.html`
*Interactive toy. Shows a small MSA (10 × 30) and pair representation (30 × 30). Student can step through row-wise attention, column-wise attention, outer-product mean, and triangle update operations. Each step shows the matrices before and after. Target aha: the axial attention + triangle update pattern iteratively builds a clean pair map from the raw MSA.*

**EMBED — Artifact #3: IPA / Structure Module Walkthrough** → `artifacts/lecture-15/03-ipa-walkthrough.html`
*Visualise one residue's SE(3) frame being updated by IPA. Show the frame as axes + position. Show query points in the frame; key points in a neighbour's frame; distance computation; weight. Step through 8 iterations; watch residues converge to their folded positions. Target aha: IPA is attention in 3D, equivariant by construction — rotating the input rotates the output the same way.*

**EMBED — Artifact #4: pLDDT and PAE Confidence Interpreter** → `artifacts/lecture-15/04-plddt-pae.html`
*Load a preset AlphaFold output (well-folded single domain / multi-domain with uncertain linker / disordered region / orphan protein with shallow MSA). Show the 3D structure coloured by pLDDT, the PAE matrix, and a natural-language interpretation. Target aha: confidence metrics tell you which parts to trust and which to ignore — using the structure without reading them is malpractice.*

---

## Part 5 — AlphaFold3 and Successors (≈20 min)

### 5.1 AlphaFold3 architecture (≈7 min)

**AlphaFold3** (Abramson et al. 2024, *Nature*). Same basic goal, different architecture. Key changes:

1. **Diffusion-based structure module**. Replaces IPA / iterative frame updates with a diffusion model: generate coordinates from Gaussian noise by iteratively denoising, conditioned on the Evoformer's pair representation. More flexible; handles multi-chain and ligand positioning natively.
2. **Multi-chain + ligand + nucleic-acid support**. The input can be a mixture: protein chains, nucleic-acid chains, small molecule ligands, ions. Output includes all components in their predicted bound state.
3. **Simpler Evoformer-like backbone (the "Pairformer")**. Still uses axial attention on pair + single repr, but drops the MSA axis once early representations are formed. Faster.

Use cases AlphaFold3 opens up:

- **Protein-drug docking**. Predict the complex of a protein with a small-molecule ligand.
- **Nucleic-acid interactions**. RNA-binding proteins with their bound RNA; DNA-binding proteins with DNA.
- **Multi-chain assemblies**. Antibody-antigen complexes; multi-subunit enzymes.

Accuracy on the benchmarks released with the paper: ~60–80% of the multi-chain / nucleic-acid tasks are solved at experimental quality, depending on the sub-task.

**FIGURE — Figure #12: AlphaFold3 diffusion structure module** → `diagrams/lecture-15/12-af3-diffusion.svg`
*Top: the forward diffusion process — start with a clean folded structure (left), add Gaussian noise across T timesteps until the coordinates are pure noise (right). Bottom: the learned reverse process — start from noise, iteratively denoise using a network conditioned on the Pairformer's pair representation, end with clean coordinates. A right-side inset shows multi-component input (protein chain + DNA + small-molecule ligand) being denoised together in one pass. Annotation: "diffusion handles variable-composition inputs and multi-modal output distributions — the AF2 IPA module cannot".*

> **EE framing — diffusion as iterative denoising**: The AlphaFold3 structure module is a **denoising diffusion probabilistic model (DDPM)** applied to 3D coordinates. Forward process: add Gaussian noise to true coordinates over $T$ timesteps until the coordinates are pure noise. Learned reverse process: a neural network predicts the noise or the clean signal at each timestep, trained across many proteins + noise levels. At inference: start from Gaussian noise, denoise iteratively conditioned on the pair representation, end with clean coordinates. Same framework as image diffusion (Ho et al. 2020, DDPM), score-matching (Song et al. 2021), and molecular conformation generation (DiffDock, 2022). The flexibility advantage over IPA: diffusion handles generation of variable-size outputs (add or remove a ligand), multi-modal distributions (multiple conformations), and conditioning on heterogeneous inputs (protein + ligand).

### 5.2 RoseTTAFold and the open-source landscape (≈6 min)

**RoseTTAFold** (Baek et al. 2021). Baker lab's independent implementation, open-source. Similar overall architecture but different specific design choices (three-track architecture: sequence, pair, 3D). Slightly below AlphaFold2 accuracy but fully available with weights. Became the workhorse for labs that couldn't run AlphaFold2 (which was weights-only for non-commercial use before 2022).

**RoseTTAFold2** (2023) and **RoseTTAFold-All-Atom** (2024) are iterative upgrades, adding multi-chain and ligand support.

**OpenFold** (Ahdritz et al. 2022). Open-source reimplementation of AlphaFold2 in PyTorch. Identical architecture, re-trained from scratch. Released with full training code and weights. Powers a great deal of downstream work.

**ColabFold** (Mirdita et al. 2022). A usability wrapper around AlphaFold2 and RoseTTAFold with fast MSA generation (replaces AF2's slow search with MMseqs2-based search). Runs in a Colab notebook in minutes. For most labs, ColabFold is what "running AlphaFold" actually means.

### 5.3 Single-sequence methods (ESMFold) (≈5 min)

**ESMFold** (Lin et al. 2023). Uses a large **protein language model** (ESM-2, a transformer trained on 65M sequences via masked-residue prediction) as the feature extractor, replacing AlphaFold2's Evoformer + MSA search. Takes a single protein sequence, passes it through ESM-2 → projects to pair representation → structure module → coordinates.

**Crucial property**: no MSA needed. Runs ~60× faster than AlphaFold2 on the same hardware because the expensive MSA search is replaced by a (more expensive but amortised) single-sequence pass through ESM-2.

Accuracy: lower than AlphaFold2 when MSAs are deep; **comparable or better** than AlphaFold2 when MSAs are shallow or absent (orphan proteins, metagenomic proteins, synthetic / designed proteins). ESMFold has been used to predict ~600M metagenomic proteins as part of the ESM Metagenomic Atlas.

**ESM-3** (2024). The latest in the ESM line. Handles proteins, structures, and function annotations in a single transformer. Can condition on partial input (e.g. "design a sequence that has this active site"). The frontier of protein foundation models.

### 5.4 The field's current state (≈2 min)

As of 2025, the practical workflow is:

- **Default**: ColabFold (AlphaFold2 variant) for most proteins.
- **Fast / orphan proteins**: ESMFold.
- **Multi-chain / ligand**: AlphaFold3 (via the EMBL-EBI webserver) or RoseTTAFold-All-Atom.
- **Custom / research**: OpenFold with fine-tuning or feature probing.

Three years after AlphaFold2's CASP14 debut, the field is broadly saturated at ~90 GDT_TS for single-domain fold-prediction; remaining open problems are multi-chain, multi-state (conformational ensembles), and disorder (discussed in Part 7).

---

## Part 6 — Inverse Folding and Protein Design (≈25 min)

### 6.1 The inverse problem (≈5 min)

Structure prediction: **sequence → structure**. Forward direction.

**Inverse folding**: **structure → sequence**. Given a desired 3D fold, what amino-acid sequence would produce it?

The inverse problem is ill-posed in a specific sense: many sequences can produce the same structure (the sequence-to-structure map is many-to-one). The interesting direction is: produce sequences that fold *stably and at high accuracy* into the target structure.

Applications:

- **Stabilisation**: given a natural protein's native structure, find a sequence that folds to the same structure but is more stable (higher melting temperature, lower aggregation).
- **Functional engineering**: given a de-novo designed structure (from Part 6.3), find a sequence that folds to it.
- **Redesign**: given a protein with a defect, find a near-native sequence that fixes it.

### 6.2 ProteinMPNN (≈8 min)

**ProteinMPNN** (Dauparas et al. 2022, *Science*). Baker-lab inverse-folding network:

- Input: **backbone-only** 3D structure (N, Cα, C, O atoms per residue; no side chains, no sequence).
- Output: a probability distribution over 20 amino acids per residue; sampled sequences.

Architecture: message-passing graph neural network. Nodes = residues with their backbone coordinates. Edges = k-nearest-neighbours in 3D. Message-passing iterations exchange geometric information between neighbours. Final per-residue softmax over 20 amino acids. Autoregressive decoder: fill in residues one by one, conditioning on already-decoded ones.

Performance: sequences designed by ProteinMPNN fold to their target structure with much higher success rate than Rosetta-designed sequences. In experimental tests (tested in the Baker lab), ~50% of ProteinMPNN designs express solubly and fold correctly, vs ~10% for Rosetta designs.

**FIGURE — Figure #9: ProteinMPNN inverse folding** → `diagrams/lecture-15/09-proteinmpnn.svg`
*Left: a target 3D backbone structure (N, Cα, C, O atoms per residue; side chains absent). Middle: the ProteinMPNN GNN — nodes are residues with their backbone frames; edges to k-nearest-neighbours; message-passing arrows. Right: per-residue amino-acid probability distributions; sampled sequence. Annotation: "~50% of sampled sequences fold correctly to the target — vastly above random or Rosetta baselines".*

> **EE framing — inverse folding as an inverse problem with strong structural priors**: Inverse folding is a classical ill-posed inverse problem — the forward map (sequence → structure) is many-to-one, so the inverse has infinite solutions. ProteinMPNN regularises by incorporating **strong structural priors**: it's trained on millions of known backbone-sequence pairs from the PDB, learning which sequences are consistent with which backbones. The autoregressive decoder provides **compositional regularisation** — residue types are not sampled independently; each depends on its neighbours. This is analogous to image super-resolution (ill-posed; prior = natural-image statistics, learned by CNN) or compressed sensing (ill-posed; prior = sparsity). Good performance depends on having a strong, correct prior.

### 6.3 De novo protein design (≈8 min)

The dream: design proteins from scratch — with structures and functions **never seen in nature**.

**Rosetta-era design** (Baker lab 2000+):

- Start with an idealised topology (e.g. a 4-helix bundle; a particular $\beta$-barrel).
- Build a backbone matching the topology.
- Use Rosetta's sequence-design module to pack amino acids that stabilise it.
- Express in *E. coli*; test folding; iterate.

Key results: fully de-novo miniproteins with novel folds (Koga et al. 2012); de-novo enzymes (some catalytic); de-novo binders.

**RFDiffusion** (Watson et al. 2023, *Nature*). Uses a diffusion model trained on the PDB to generate novel protein backbones. Input: specifications like "a 4-helix bundle", "binds this ligand here", "matches this active site". Output: a backbone structure satisfying the constraints.

Workflow:

1. **RFDiffusion**: specify constraints, sample a backbone.
2. **ProteinMPNN**: generate amino-acid sequence for that backbone.
3. **AlphaFold2 / ESMFold validation**: check that the designed sequence folds back to the target backbone (self-consistency).
4. **Experimental validation**: synthesise, express, characterise.

Success rate: in a 2023 paper from the Baker lab, ~10–40% of RFDiffusion + ProteinMPNN designs expressed, folded, and exhibited the desired function (e.g. binding a small-molecule target). Without AlphaFold / ProteinMPNN, equivalent success rates were <1%.

> **Intuition box**: The design loop is a two-player relay. RFDiffusion is the sculptor — it hallucinates a plausible 3D shape meeting your constraints, ignoring whether any sequence can fold to it. ProteinMPNN is the scriptwriter — given that shape, it writes a sequence that, according to everything it's seen in nature, ought to fold to that shape. AlphaFold plays the critic — it folds the sequence and checks whether the result matches the sculptor's shape. If the three agree, you have a candidate to synthesise. The critic's role is crucial: before AlphaFold, there was no fast way to check whether a designed sequence actually folded as intended, so designs mostly failed at the wet-lab step. Adding a reliable folding critic inside the loop is what flipped design success rates from <1% to 10–40%.

**FIGURE — Figure #10: RFDiffusion protein design pipeline** → `diagrams/lecture-15/10-rfdiffusion.svg`
*Horizontal workflow. Panel 1: constraint specification (e.g. "binds Zn²⁺ at this position"). Panel 2: RFDiffusion-generated backbone (novel fold). Panel 3: ProteinMPNN sequence design. Panel 4: AlphaFold2 validation — does the designed sequence fold back to the target backbone? Panel 5: experimental validation. Annotation: "10–40% of designed proteins express + fold + function; previously <1% with Rosetta alone".*

### 6.4 Safety and ethics (≈4 min)

Protein design is a **dual-use** technology:

- **Beneficial uses**: therapeutic binders, novel enzymes for green chemistry, malaria vaccines, COVID-era receptor-binding-domain mimics, enzymes for plastic degradation.
- **Concerning uses**: designed toxins, designed pathogen proteins.

The field has started grappling with biosecurity. The Baker lab has adopted a screening policy (don't design certain classes of proteins); Anthropic's ASL-3 and similar thresholds explicitly reference bioweapon-uplift capabilities. Sharing model weights is now a deliberate policy decision, not an automatic open-source default.

As EE students with ML backgrounds, you may find yourself building these models. Know the dual-use issues.

> **Discussion prompt**: You've just trained a new RFDiffusion-like model that can design proteins with better binding affinity than the 2023 state-of-the-art. Before releasing weights + code, what questions should you ask? (Starting points: biosecurity risk assessment — what known bioweapon-adjacent proteins can your model improve vs baseline? Responsible release — weights public, restricted, or private? Dataset provenance — did training data contain curated dual-use sequences? Monitoring — can you detect misuse post-release?) No right answer; the goal is to learn to ask.

### 6.5 The design-predict-validate loop (≈0 min — covered in §6.3 above)

**EMBED — Artifact #5: ProteinMPNN Inverse-Folding Demo** → `artifacts/lecture-15/05-proteinmpnn.html`
*Pick a target backbone (preset mini-proteins: 4-helix bundle, $\beta$-barrel, random native). Run a simplified ProteinMPNN (toy message-passing + softmax). Sample 5 sequences. Check sequence diversity and predicted foldability (a toy self-consistency score). Target aha: many sequences fold to the same backbone; ProteinMPNN finds the manifold quickly.*

---

## Part 7 — What This Changes (≈20 min)

### 7.1 The AlphaFold Database (≈6 min)

**AlphaFold Database** (Varadi et al. 2022, EMBL-EBI hosted). Predicted structures for ~200M proteins covering almost all of UniProt. Public, freely accessible. Each structure tagged with pLDDT colouring and PAE matrix.

Before 2021, the PDB had ~180k experimentally-determined structures. After the AlphaFold DB, the ratio of "available structural information" to "unexplored sequence" flipped — from <0.1% of proteins having any structure to >99%.

Users of the AlphaFold DB:

- **Drug discovery**. Every protein target now has a predicted structure to inspect. Druggability / pocket prediction, structure-based virtual screening, all feasible on targets that previously had no structure.
- **Functional annotation**. Structure → function inference (detecting fold similarity to known proteins of known function). Pansa et al. 2024 identified ~10% of previously-unannotated proteins via structural homology.
- **Evolutionary studies**. Comparing predicted structures across species reveals structural conservation not visible from sequence.
- **Enzyme engineering**. Given a predicted structure with a specific pocket, design mutations to alter substrate specificity.

The DB is updated regularly; v2 (2024) adds cryptic-pocket predictions and ligand-binding probability scores.

**FIGURE — Figure #11: AlphaFold Database scale** → `diagrams/lecture-15/11-af-database.svg`
*A comparison bar chart. Left bar: PDB experimental structures as of 2021 (~180k). Right bar: AlphaFold DB predicted structures (~200M). Log-scale x-axis. Below the bars: a small grid showing pLDDT-quality buckets: >90 (very high, 40%), 70–90 (confident, 30%), 50–70 (low, 20%), <50 (disordered / orphan, 10%). Annotation: "1000× more protein structural information than the combined history of X-ray crystallography".*

### 7.2 Drug discovery (≈5 min)

Pre-AlphaFold, **structure-based drug design** required an experimentally-determined structure of the target. X-ray crystallography or cryo-EM, years of work, many failed targets. Most drug targets had no structure.

Post-AlphaFold: every human protein has a predicted structure. Drug discovery pipelines now include AlphaFold predictions as a routine input step.

Caveats:

- AlphaFold predicts the **apo** (unbound) structure. The drug-bound (holo) structure can differ, sometimes significantly. For kinases, AlphaFold predicts the active conformation; many drugs bind the inactive conformation.
- **Conformational ensembles**: proteins don't have a single structure; they fluctuate. Drug binding often exploits transient states AlphaFold doesn't capture.
- **Allosteric sites**: pockets that only form in specific conformational states are poorly predicted.

Nonetheless, the pharmaceutical industry has absorbed AlphaFold deeply. Every major pharma now has an AlphaFold-based workflow for target validation.

### 7.3 What's still hard (≈5 min)

Three open frontiers:

- **Conformational dynamics**. A protein isn't one structure; it's an ensemble. AlphaFold predicts a dominant state but not the distribution. Molecular dynamics simulations complement but are expensive. Recent work (AF2-multimer with MSA subsampling, AlphaFlow, Distributional Graphormer) addresses this partially.
- **Intrinsically disordered regions**. ~30% of eukaryotic proteins have significant disordered content. AlphaFold predicts these as "low pLDDT" — it knows they're uncertain but doesn't give you a distribution over their states. These are often functional (transcription factors, signalling motifs).
- **Multi-state proteins**. Kinases, transporters, ribosomes cycle through distinct functional states. AlphaFold predicts one state (usually the most abundant in the training data). For multi-state function, predicting all states is required.

> **Warning box**: AlphaFold predictions are **not experimental structures**. For publication, regulatory, or clinical work that requires real structural data — drug binding modes, mutagenesis experiments targeting specific atoms, structure-function causation claims — the prediction is a starting hypothesis, not a definitive answer. Always cross-check with experiment when the stakes are real.

### 7.4 Single-cell / multi-modal extensions (≈2 min)

Emerging work extends structure prediction to:

- **Protein-RNA complexes**: RoseTTAFold-Nucleic, AlphaFold3.
- **Protein-DNA complexes**: AlphaFold3, RoseTTAFold2.
- **Large assemblies**: the ribosome, the nuclear pore complex (hundreds of chains). AlphaFold3 handles 5–10 chains well; 100+ is still research.
- **Membrane proteins**: AlphaFold was trained on soluble proteins; membrane protein prediction is strong but benefits from membrane-aware fine-tuning (AlphaFold-multimer variants).

### 7.5 AI safety and the regulatory landscape (≈2 min)

As of 2024–2025, **biosecurity screening** is becoming norm for major model releases. Meta's ESM-3 underwent biosecurity review before release; Anthropic's ASL-3 explicitly references bioweapon-uplift risks; the US Executive Order on AI (2023) requires major model developers to notify the government of frontier capabilities.

The field is at an interesting moment where the open-science norm (Part 7 of L14) collides with dual-use concerns. How this resolves affects the pace of future progress.

**EMBED — Artifact #6: AlphaFold Database Explorer** → `artifacts/lecture-15/06-afdb-explorer.html`
*Preset catalogue of ~20 proteins (well-folded single-domain, multi-domain, orphan, disordered, engineered de-novo). For each, show the 3D structure coloured by pLDDT, the PAE matrix, and a metadata panel with species / length / UniProt ID / function. Filter by pLDDT / length / disorder fraction. Target aha: the database is searchable by confidence — you can focus on high-quality regions and ignore uncertain ones.*

**EMBED — Artifact #7: Protein Structure Level Explorer** → `artifacts/lecture-15/07-structure-levels.html`
*Interactive viewer showing a small protein at each structural level: linear sequence, secondary-structure-highlighted sequence, tertiary 3D backbone, quaternary assembly (if applicable). Toggle between levels; rotate the 3D view; click residues to highlight across views. Target aha: a protein isn't "a string of letters" or "a blob" — it's a hierarchy, and each level informs different analyses.*

---

## Wrap-up (≈10 min)

### What you should take away

- **Protein structure is four levels**. Primary = sequence; secondary = local H-bond patterns; tertiary = 3D fold; quaternary = multi-chain assembly. Tertiary is what "structure prediction" usually means.
- **Anfinsen established that sequence → structure is a well-defined problem**. The 50-year journey to AlphaFold2 was finding the algorithm.
- **The classical era (homology modelling, threading, Rosetta) plateaued** because it lacked a way to exploit the rich evolutionary information in MSA coevolution.
- **DCA unlocked coevolution-based contact prediction** by computing inverse covariances of amino-acid co-occurrence. EE framing: it's graphical-model estimation on categorical data.
- **AlphaFold2's Evoformer + structure module** is the breakthrough architecture. Evoformer refines (MSA, pair) iteratively via axial attention + triangle updates; structure module converts to 3D via $\text{SE}(3)$-equivariant Invariant Point Attention.
- **pLDDT and PAE are calibrated confidence metrics**. Use them — a structure without reading its confidence is dangerous.
- **AlphaFold3, ESMFold, RoseTTAFold** extend the toolbox: multi-chain + ligand (AF3), single-sequence / speed (ESMFold), open-source + multi-modal (RFAA).
- **ProteinMPNN and RFDiffusion** enable de-novo design. Inverse folding = inverse problem with strong structural prior. Design loop = RFDiffusion (backbone) → ProteinMPNN (sequence) → AF2 (validation) → experiment.
- **The AlphaFold Database (~200M predicted structures)** has flipped the field. Drug discovery, functional annotation, and evolutionary studies now assume a predicted structure.
- **What's still hard**: dynamics, disorder, multi-state function, large assemblies, dual-use / biosecurity.

### Next lecture

ML in genomics: architectures, pitfalls, frontiers. A synthesis across previous lectures — pileups → CNN (DeepVariant); long-range sequence regulation → dilated conv + transformer (Enformer, Borzoi); count matrices → VAE (scVI); protein 3D → equivariant transformer (AlphaFold); molecular graphs → GNN. Plus the shared pitfalls: out-of-distribution generalisation, calibrated confidence, data leakage at biological scale.

### Homework

1. Install ColabFold locally (or use the web notebook). Predict structures for three of your favourite proteins (you choose): one with a close PDB homologue, one with a moderate MSA, one orphan. Compare pLDDT distributions. Which does best and why?
2. From the AlphaFold Database, fetch predicted structures for ten human "dark proteins" (uncharacterised, lacking PDB homologues). What fraction have high-confidence (>70 pLDDT) predictions? Pick one and use foldseek or Dali to find structural homologues; does it suggest a function?
3. Use ProteinMPNN (available through ColabDesign or Baker lab GitHub) to generate 10 sequences for a given backbone. Evaluate the sequence recovery rate (fraction of positions matching the original native) and AlphaFold's self-consistency score (does AF2 fold the new sequence back to the same backbone?).
4. Read the AlphaFold2 paper (Jumper et al. 2021) Methods section. For one Evoformer sub-operation (your choice), write a one-page explanation in your own words, including the tensor shapes at input and output.
5. Pick an AlphaFold-predicted structure with a large low-pLDDT region. Is it disordered (check DisProt or IUPred) or is it genuinely uncertain? How would you distinguish these experimentally?

### Recommended reading

- Anfinsen, C. B. (1973). Principles that govern the folding of protein chains. *Science* 181, 223–230.
- Jumper, J., Evans, R., Pritzel, A., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
- Abramson, J., Adler, J., Dunger, J., et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630, 493–500.
- Baek, M., DiMaio, F., Anishchenko, I., et al. (2021). Accurate prediction of protein structures and interactions using a three-track neural network. *Science* 373, 871–876. (RoseTTAFold.)
- Lin, Z., Akin, H., Rao, R., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure. *Science* 379, 1123–1130. (ESMFold.)
- Dauparas, J., Anishchenko, I., Bennett, N., et al. (2022). Robust deep learning based protein sequence design using ProteinMPNN. *Science* 378, 49–56.
- Watson, J. L., Juergens, D., Bennett, N. R., et al. (2023). De novo design of protein structure and function with RFdiffusion. *Nature* 620, 1089–1100.
- Marks, D. S., Colwell, L. J., Sheridan, R., et al. (2011). Protein 3D structure computed from evolutionary sequence variation. *PLoS ONE* 6, e28766. (EVfold.)
- Morcos, F., Pagnani, A., Lunt, B., et al. (2011). Direct-coupling analysis of residue coevolution captures native contacts across many protein families. *PNAS* 108, E1293–E1301. (DCA.)
- Varadi, M., Anyango, S., Deshpande, M., et al. (2022). AlphaFold Protein Structure Database. *Nucleic Acids Research* 50, D439–D444.
- Mirdita, M., Schütze, K., Moriwaki, Y., et al. (2022). ColabFold: making protein folding accessible to all. *Nature Methods* 19, 679–682.
- Baker, D. (2019). What has de novo protein design taught us about protein folding and biophysics? *Protein Science* 28, 678–683.
- AlphaFold database: <https://alphafold.ebi.ac.uk/>
- ESM Metagenomic Atlas: <https://esmatlas.com/>
- ColabFold: <https://github.com/sokrypton/ColabFold>
- RFDiffusion: <https://github.com/RosettaCommons/RFdiffusion>
- ProteinMPNN: <https://github.com/dauparas/ProteinMPNN>
- Foldseek (structural search): <https://search.foldseek.com/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Protein Structure Basics                               | 25&nbsp;min | 0:25 |
| Part 2 — The Classical Era                                       | 20&nbsp;min | 0:45 |
| Part 3 — MSA-Based Methods and Coevolution                       | 30&nbsp;min | 1:15 |
| Part 4 — AlphaFold2 Architecture                                  | 60&nbsp;min | 2:15 |
| Part 5 — AlphaFold3 and Successors                                | 20&nbsp;min | 2:35 |
| Part 6 — Inverse Folding and Protein Design                        | 25&nbsp;min | 3:00 |
| Part 7 — What This Changes                                        | 20&nbsp;min | 3:20 |
| Wrap-up                                                             | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
