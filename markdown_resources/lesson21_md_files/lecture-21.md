# Lecture 21 (proposed L7) — HMMs, Profile HMMs, and Gene Finding

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology assumed
> **File**: to be rendered as `lectures/lecture-21.html` (provisional name; renumber to `lecture-07.html` when curriculum is reordered)

> **Proposed placement**: insert after L4 (variant calling, becomes new L6 in reordered curriculum), so this becomes the new L7 between variant calling and bulk RNA-seq. The placement is deliberate: gene finding turns an assembled, variant-called genome into a usable annotation; profile HMMs (HMMER, Pfam) operationalise the MSA → profile pipeline introduced in the new L4; and the Viterbi-on-DAG idea recurs in the long-reads + pangenome lecture (L11→L15). Many later lectures (RNA-seq, ChIP-seq peak calling, chromatin segmentation via ChromHMM) make HMM-class assumptions; this lecture is where students see the foundations.

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Define a Hidden Markov Model: states, transitions, emissions, initial distribution; explain why the Markov assumption holds and how it breaks.
2. Solve the three canonical HMM problems: forward / backward (likelihood), Viterbi (decoding), Baum-Welch (training); recognise EM and belief propagation under the hood.
3. Construct profile HMMs from MSAs (HMMER, Pfam): match, insert, delete states; the 3-state-per-column architecture.
4. Apply profile HMMs to remote-homology detection; interpret HMMER output (sequence and domain E-values).
5. Walk through eukaryotic gene finding (GENSCAN, AUGUSTUS): exon / intron / splice / UTR states; combination of HMM with sequence priors.
6. Connect HMMs to chromatin segmentation (ChromHMM), splice-aware alignment (HISAT2 graph extensions), and the Viterbi-on-DAG of L11 / new L15.
7. Recognise modern alternatives: CRFs, RNNs, transformers — and where HMMs still win.

---

## Part 1 — HMM Foundations (≈30 min)

### 1.1 Why Markov chains aren't enough (≈4 min)

A **Markov chain** on observations: $P(x_t \mid x_{t-1}, x_{t-2}, ..., x_1) = P(x_t \mid x_{t-1})$. The next observation depends only on the current one.

For DNA, this is too weak: a CpG island's "C is followed by G" pattern requires distinguishing two regimes — inside-CpG-island vs outside — that a single Markov chain on letters can't represent.

The fix: introduce a **latent state** that captures "which regime are we in?". Now:

- Latent state $z_t$ at each position, with a Markov transition $P(z_t \mid z_{t-1})$.
- Observation $x_t$ depends on the state via emission $P(x_t \mid z_t)$.

This is a **Hidden Markov Model**.

### 1.2 The HMM definition (≈5 min)

Formally:

- States $\{S_1, S_2, ..., S_K\}$.
- Initial distribution $\pi_i = P(z_1 = S_i)$.
- Transition probabilities $a_{ij} = P(z_t = S_j \mid z_{t-1} = S_i)$.
- Emission probabilities $b_i(x) = P(x_t = x \mid z_t = S_i)$.
- Sequence of observations $X = x_1, ..., x_T$.
- Sequence of latent states $Z = z_1, ..., z_T$.

The model is parameterised by $\theta = (\pi, a, b)$.

Joint probability:
$$P(X, Z \mid \theta) = \pi_{z_1} \prod_{t=2}^{T} a_{z_{t-1}, z_t} \prod_{t=1}^{T} b_{z_t}(x_t)$$

### 1.3 The three canonical problems (≈5 min)

Given an HMM, three things you might want:

1. **Likelihood**: $P(X \mid \theta)$. How probable is this observation? Solved by **forward algorithm**.
2. **Decoding**: argmax$_Z$ $P(Z \mid X, \theta)$. What's the most probable state sequence? Solved by **Viterbi algorithm**.
3. **Training**: argmax$_\theta$ $P(X \mid \theta)$. Given observations only, learn the parameters? Solved by **Baum-Welch (EM)**.

Each has $O(T K^2)$ time complexity — efficient.

### 1.4 The CpG island example (≈8 min)

A canonical first example. CpG dinucleotides are rare in mammalian DNA (because methylated C tends to mutate to T). But CpG **islands** — short genomic regions with elevated CpG density — are functional regulatory features (often near promoters).

Two-state HMM:

- State 1: "CpG island" (high CpG frequency).
- State 2: "background" (low CpG frequency).

Transitions: small probability of switching states ("islands are rare"), high self-transition ("once inside an island, you stay for a while").

Emissions: 16 dinucleotide probabilities per state (CG, CA, CT, ..., GG).

Train on labelled data → run Viterbi on a new genome → get the CpG-island annotation.

This pattern — **two-or-more state HMM + Viterbi on the genome → segmentation** — recurs throughout genomics: ChromHMM, gene finders, peak callers (with HMM extensions), variant callers handling CNV.

### 1.5 The deep dive (≈8 min)

> **EE framing — HMMs as state-space models on discrete symbols**: An HMM is a **discrete-state-space dynamical system** with stochastic transitions and stochastic observations — exactly the setting of Kalman filtering with a discrete state. Forward-backward is the discrete-time-discrete-state analog of the Kalman filter / smoother. Viterbi is the **MAP decoder** for the most probable trajectory. Baum-Welch is **EM** for the latent state and parameters jointly. The HMM was published in 1966; Kalman 1960. The two communities developed independently for decades before the convergence in the 1990s. For an EE student, the punchline is: every HMM technique you'll use here has a state-space-model analog you've seen in signal processing — and the algorithms transfer.

**FIGURE — Figure #1: HMM graph and joint-probability factorisation** → `diagrams/lecture-21/01-hmm-graph.svg`
*Top: a 4-state HMM as a directed graph: states as nodes, transitions as labelled edges; emission distributions as histogram annotations. Bottom: the unrolled-in-time graphical model for $T$ observations: states $z_1, z_2, ..., z_T$ in a horizontal chain; observations $x_t$ hanging below each state. Annotations on conditional independence: $z_t$ depends only on $z_{t-1}$ given the chain.*

---

## Part 2 — The Three HMM Algorithms (≈45 min)

### 2.1 Forward algorithm (≈10 min)

Goal: compute $P(X \mid \theta) = \sum_Z P(X, Z \mid \theta)$.

Naive: sum over $K^T$ state paths — exponential.

Forward DP: define $\alpha_t(i) = P(x_1, ..., x_t, z_t = S_i \mid \theta)$. Then:

- Initialise: $\alpha_1(i) = \pi_i b_i(x_1)$.
- Recurse: $\alpha_t(j) = b_j(x_t) \sum_i \alpha_{t-1}(i) a_{ij}$.
- Terminate: $P(X \mid \theta) = \sum_i \alpha_T(i)$.

$O(T K^2)$ time, $O(T K)$ space.

### 2.2 Backward algorithm (≈6 min)

Symmetric: $\beta_t(i) = P(x_{t+1}, ..., x_T \mid z_t = S_i, \theta)$.

- Initialise: $\beta_T(i) = 1$.
- Recurse: $\beta_t(i) = \sum_j a_{ij} b_j(x_{t+1}) \beta_{t+1}(j)$.

Combined with forward: $P(z_t = S_i \mid X, \theta) = \alpha_t(i) \beta_t(i) / P(X)$.

This is **per-position posterior** — useful for confidence (Posterior Decoding).

### 2.3 Viterbi algorithm (≈10 min)

Goal: argmax$_Z P(Z, X \mid \theta)$. Most probable state sequence.

Define $\delta_t(i) = \max_{z_1, ..., z_{t-1}} P(z_1, ..., z_{t-1}, z_t = S_i, x_1, ..., x_t \mid \theta)$.

- Initialise: $\delta_1(i) = \pi_i b_i(x_1)$.
- Recurse: $\delta_t(j) = b_j(x_t) \max_i [\delta_{t-1}(i) a_{ij}]$.
- Backtrack: store argmax pointers; trace from $\arg\max_i \delta_T(i)$.

The forward algorithm is "sum over paths"; Viterbi is "max over paths". Both run in $O(T K^2)$.

In log space (to avoid underflow):
$$\log \delta_t(j) = \log b_j(x_t) + \max_i [\log \delta_{t-1}(i) + \log a_{ij}]$$

A pure Viterbi computation looks identical to log-domain pairwise alignment, which is no accident — pairwise alignment IS a Viterbi computation on a trivial 1-state HMM.

**FIGURE — Figure #2: Viterbi trellis with backtrace** → `diagrams/lecture-21/02-viterbi-trellis.svg`
*A 4-state HMM unrolled over $T = 8$ time steps. Trellis lattice with one node per (state, time) cell. Edges show transitions; node colour intensity proportional to $\log \delta_t(i)$. The Viterbi backtrace path drawn in bold red from end to start. Annotation: "Viterbi = max-plus DP; same shape as pairwise alignment."*

### 2.4 Baum-Welch (EM training) (≈10 min)

Given observations $X$ but no labels $Z$, learn $\theta$.

EM iterates:

- **E-step**: compute posterior expectations $\gamma_t(i) = P(z_t = S_i \mid X, \theta)$ and $\xi_t(i, j) = P(z_t = S_i, z_{t+1} = S_j \mid X, \theta)$ using forward-backward.
- **M-step**: re-estimate parameters:
  $$\hat{\pi}_i = \gamma_1(i)$$
  $$\hat{a}_{ij} = \frac{\sum_t \xi_t(i, j)}{\sum_t \gamma_t(i)}$$
  $$\hat{b}_i(x) = \frac{\sum_{t: x_t = x} \gamma_t(i)}{\sum_t \gamma_t(i)}$$

Converges to a local optimum. Initialisation matters; restart from multiple random points.

Convergence is monotonic in $P(X \mid \theta^{(k)})$ — that's the EM guarantee.

### 2.5 The connection to belief propagation (≈4 min)

Forward-backward IS belief propagation on the chain-structured graphical model that defines the HMM. The forward pass computes "messages from past"; the backward pass computes "messages from future"; their product gives the per-position posterior.

The connection is exact: forward-backward is BP on a chain. Felsenstein's pruning (L4-new) is BP on a tree. The HMM and the phylogenetic tree differ only in graph structure — algorithmically they're the same family.

### 2.6 Practical tips (≈5 min)

- **Log space everything**: probabilities multiply across long sequences, underflowing immediately. Convert to log-prob.
- **Pseudocounts**: avoid zero probabilities in $a$ and $b$ from rare states (Laplace smoothing or Dirichlet prior).
- **Multiple Baum-Welch restarts**: EM is local. 5–10 random initialisations is standard.
- **Hyperparameter $K$ (state count)**: cross-validation or BIC/AIC. BIC is the safer default for HMM model selection.

---

## Part 3 — Profile HMMs and HMMER (≈35 min)

### 3.1 The profile-HMM idea (≈4 min)

A **profile HMM** is an HMM whose state structure encodes an MSA:

- One **match (M) state** per column of the MSA.
- One **insert (I) state** per inter-column gap.
- One **delete (D) state** per column allowing the column to be skipped.

Match-state emission distributions reflect the column's amino-acid composition. Transition probabilities reflect the gap statistics of the input MSA.

Result: a probabilistic model of the protein family that scores any new sequence by Viterbi alignment to the model.

### 3.2 The Krogh-Mian-Mitchison architecture (≈5 min)

The standard profile-HMM architecture (Krogh, Mian, Sjölander, Haussler 1994):

```
       +---I0           +---I1           +---I2
       |   |            |   |            |   |
B ---> M1 --> M2 ----> M3 --> M4 ----> M5 ---- E
       |   ^            |   ^            |
       D1            D2            D3
```

States:

- **M (match)**: position-specific emission of an amino acid.
- **I (insert)**: position-specific insertion of one or more amino acids.
- **D (delete)**: position skipped (no emission).
- **B (begin)**: start state.
- **E (end)**: end state.

Transitions are tightly constrained: from $M_k$ you can go to $M_{k+1}$, $I_k$, or $D_{k+1}$. From $I_k$, you can go to $M_{k+1}$, $I_k$ (self-loop), or $D_{k+1}$. From $D_k$, $M_{k+1}$, $D_{k+1}$, or $I_k$.

This 3-state-per-column structure — match / insert / delete — is the **canonical** profile-HMM topology.

### 3.3 HMMER (≈8 min)

**HMMER** (Eddy 1998, 2011): the dominant profile-HMM software.

- `hmmbuild`: from MSA → profile HMM.
- `hmmsearch`: query a profile HMM against a sequence database.
- `hmmscan`: query a sequence against a database of profile HMMs.

The `nhmmer` variant handles nucleotide sequences (e.g., for non-coding RNA detection).

Output:

- **Sequence E-value**: probability the entire sequence is from a random null model.
- **Domain E-value**: probability for the specific domain alignment.
- **Bit score**: log-odds of model vs null.
- **Bias correction**: adjusts for compositionally-biased queries (similar to BLAST's composition correction, L19-new).

### 3.4 The Pfam database (≈5 min)

**Pfam** (Bateman et al. 2002, Mistry et al. 2021, now hosted at InterPro): a curated database of ~20,000 protein families. Each family has:

- A seed alignment (manually curated MSA of representative sequences).
- A full alignment (HMMER auto-aligned to all UniProt sequences with $E < 10^{-3}$).
- A profile HMM built from the seed.
- Metadata: function, GO terms, structural classification, literature.

Pfam profile HMMs are the **default** for protein-domain annotation. When you assemble a new genome and translate ORFs, the first protein-level annotation is "scan against Pfam".

### 3.5 Worked example (≈5 min)

**FIGURE — Figure #3: Profile HMM as MSA → 3-state-per-column model** → `diagrams/lecture-21/03-profile-hmm.svg`
*Top: a 5-row × 8-column MSA. Bottom: the corresponding profile HMM with 8 match-state columns. Each match state shows its emission distribution as a small bar chart; insert states drawn as self-loops above; delete states drawn as bypass arrows below. Side annotation: "the model encodes both the column-wise composition and the gap statistics."*

### 3.6 The deep dive (≈8 min)

> **EE framing — profile HMMs as position-specific filters**: A profile HMM is a **position-specific matched filter** with explicit gap states. Compared to a fixed-template matched filter (which assumes the signal length is known), the profile HMM's insert and delete states allow **length variability**. The Viterbi alignment of a query to a profile HMM is the analog of dynamic-time-warping (DTW) in speech recognition: aligning a variable-rate query to a fixed template, allowing local stretches and skips. Indeed, DTW and profile HMM alignment share the same DP machinery; they were developed in parallel in the 1970s-80s and unified in the 1990s.

The HMMER bit-score $S$ has the standard log-likelihood-ratio interpretation; Karlin-Altschul-style E-value statistics apply (with HMM-specific calibration).

---

## Part 4 — Gene Finding (≈40 min)

### 4.1 The eukaryotic gene-finding problem (≈4 min)

Given an assembled eukaryotic genome (~3 Gbp for human), find:

- Gene boundaries (start and stop).
- Exon-intron structure.
- Coding (CDS) regions.
- 5'/3' untranslated regions (UTRs).
- Promoter regions (less commonly).

This is harder than prokaryotic gene finding (which has near-uniform coding density and no introns) because eukaryotic genes are sparse, fragmented, and overlap with regulatory elements.

### 4.2 GENSCAN (≈10 min)

**GENSCAN** (Burge & Karlin 1997): the canonical statistical gene finder.

Architecture: a giant HMM with state classes:

- **Intergenic** (single state).
- **Promoter / 5'-UTR** states.
- **Initial exon** (with start codon emission).
- **Internal exon** states.
- **Final exon** (with stop codon emission).
- **Single-exon gene** state.
- **Intron** states (with various phase tracking — coding frame can be 0, 1, or 2 across an intron).
- **3'-UTR** states.
- **Polyadenylation signal** state.

Each state has length distributions (geometric for introns, gamma for exons). Emission probabilities reflect codon-usage statistics, splice-site consensus sequences, and CpG-island context.

Viterbi decoding on this HMM produces a complete gene structure annotation.

### 4.3 Splice site detection (≈6 min)

A critical sub-problem: identify donor and acceptor splice sites within introns.

- **Donor (5' splice site)**: GT consensus at intron start.
- **Acceptor (3' splice site)**: AG consensus at intron end, plus a polypyrimidine tract upstream.

These are detected by **position-specific weight matrices (PWMs)** on a small window around each candidate site. The PWM is exactly a profile HMM with no inserts or deletes — the simplest case.

Modern tools layer **deep learning** (SpliceAI, L17-new = clinical genomics) on top of PWM-style models for higher precision.

### 4.4 AUGUSTUS (≈5 min)

**AUGUSTUS** (Stanke et al. 2003-2008): the modern eukaryotic gene finder. Improvements over GENSCAN:

- **Conditional Random Field (CRF) version**: more flexible than HMM; conditions on observations rather than generative.
- **Hint integration**: incorporates external evidence (RNA-seq alignments, protein-protein hits, conservation tracks).
- **Iterative training**: bootstraps from initial gene set on the new genome.

For a new mammalian genome assembly, AUGUSTUS + RNA-seq hints is the standard.

### 4.5 RNA-seq-driven annotation (≈5 min)

A modern approach: **train the gene finder on the genome's own RNA-seq data**.

Pipeline:

1. Map RNA-seq reads to the genome (HISAT2, STAR — L5).
2. Assemble transcripts (StringTie, Cufflinks).
3. Use the assembled transcripts as training set for AUGUSTUS / BRAKER.
4. Run gene finder; integrate gene model predictions with transcript evidence.

Output: a gene set that's both algorithmically predicted and experimentally supported.

### 4.6 Prokaryotic gene finding (≈5 min)

Simpler than eukaryotic: most prokaryotic genes are continuous (no introns), high gene density (~1 gene / kb), and well-conserved start codons.

Tools: **GeneMark**, **Glimmer**, **Prodigal**. Use shorter HMMs with codon-frequency emission. Run in ~minutes on a typical bacterial genome (5 Mb).

For metagenomics (mixed organisms), MetaGeneMark and Prodigal-meta handle community-level gene finding.

### 4.7 Worked example (≈5 min)

**FIGURE — Figure #4: GENSCAN HMM architecture** → `diagrams/lecture-21/04-genscan.svg`
*A high-level state diagram of GENSCAN. State classes (intergenic, promoter, exon-initial, exon-internal, exon-final, intron-phase-0/1/2, UTR-5/3) drawn as boxes; canonical transitions as arrows; length distributions noted (geometric for intron, gamma for exon). Annotation: "Viterbi on this HMM segments a 100 Mb chromosome into a list of genes in minutes."*

---

## Part 5 — HMM Applications Beyond Gene Finding (≈25 min)

### 5.1 ChromHMM (≈8 min)

**ChromHMM** (Ernst & Kellis 2012): epigenome segmentation. Uses a multi-mark HMM to integrate ChIP-seq tracks for multiple histone modifications and call distinct chromatin states.

States (typical 15-state model):

- Active promoter (high H3K4me3, H3K27ac).
- Active enhancer (H3K4me1, H3K27ac).
- Strong enhancer (additional H3K27ac).
- Weak enhancer.
- Insulator (CTCF-bound).
- Strong transcription (H3K36me3 in gene body).
- Polycomb-repressed (H3K27me3).
- Heterochromatin (H3K9me3).
- ... etc.

Each state has emission probabilities for each histone mark; the HMM is trained on observed mark co-occurrences.

Practical: download ChromHMM, run on your own ChIP-seq data, get a chromatin-state segmentation. Used in ENCODE, Roadmap Epigenomics — connects directly to L9 / new L13 (ChIP-seq).

### 5.2 Splice-aware alignment in HISAT2 (≈4 min)

**HISAT2** uses an HMM extension to handle splicing:

- Reads can span an intron.
- The aligner treats intron states as "skip-no-emission" states in a graph-aware HMM.
- Splice junction quality is scored by junction-likelihood.

This is the engine behind RNA-seq alignment in modern pipelines (L5).

### 5.3 Pair-HMMs in variant calling (≈4 min)

**HaplotypeCaller** (GATK, L4): the per-read alignment likelihood is computed via a **pair-HMM** — a 3-state HMM for matches, insertions, deletions in pairwise alignment.

Pair-HMMs handle indel-aware alignment with explicit gap-open / gap-extend transitions; they're the basis of all modern haplotype-based variant callers (DeepVariant uses CNN on pair-HMM-aligned pileups).

### 5.4 Connection to Viterbi-on-DAG (≈5 min)

Lecture 11 / new L15 (long reads + pangenome) introduces **Viterbi on a DAG** — alignment of a long read to a pangenome graph. The same machinery: states = positions in the graph; transitions = edges. The HMM background you have here makes Viterbi-on-DAG immediately interpretable.

### 5.5 What HMMs can't do (≈4 min)

HMMs assume:

- **Markov property**: state depends only on previous state. Long-range dependencies (e.g., trans-splicing across a 100 kb intron) violate this.
- **Conditional independence**: emissions are independent given state. Real dependencies (e.g., codon-context correlations) are imperfectly captured.
- **Discrete state space**: continuous variables (e.g., expression levels) need extension to GMM-HMM or autoregressive HMM.

For genomic data, these are mostly fine. For natural language and music, **transformers replaced HMMs** because long-range dependencies dominate. Genomics is in the middle: HMMs still dominate gene finding and chromatin segmentation; transformers are gaining ground on regulatory prediction (Enformer, L9).

---

## Part 6 — Modern Alternatives and the HMM Niche (≈20 min)

### 6.1 CRFs (≈5 min)

**Conditional Random Fields** generalize HMMs by removing the generative assumption:

- HMM models $P(X, Z)$.
- CRF models $P(Z \mid X)$ directly.

CRFs allow rich feature functions on observations without the conditional independence constraint. AUGUSTUS uses CRF for gene finding; ChromHMM has CRF variants (Segway).

### 6.2 Deep learning replacements (≈5 min)

For sequence classification tasks, RNNs and transformers have largely replaced HMMs:

- **SpliceAI**: deep CNN on splice-site detection.
- **DeepBind / DeepSEA**: CNN on transcription-factor binding.
- **Enformer / Borzoi**: dilated CNN + transformer on regulatory sequence (L9, L16-new).

These models have higher accuracy but are **less interpretable** than HMMs and require much more training data.

### 6.3 The HMM niche (≈5 min)

HMMs still win for:

- **Small training data**: when you only have a handful of labelled examples (e.g., a new pathogen's gene structure), HMMs train robustly.
- **Interpretable annotation**: when each state needs a human-readable label.
- **Sparse-state inference**: when most positions are in one of a few well-defined states (chromatin segmentation, gene segmentation).
- **Deterministic decoding**: when you need a single most-probable assignment, not a probability distribution.

### 6.4 The hybrid future (≈5 min)

Most modern bioinformatics pipelines use **HMMs + deep learning hybrids**:

- HMM segments first; deep learning refines per-state predictions.
- Or: deep learning extracts features; HMM provides a calibrated state assignment.

Examples: BRAKER (AUGUSTUS HMM + deep learning hints), DeepHMM (deep emissions on HMM topology), HMMRATAC (HMM segmentation on ATAC-seq with deep features).

**FIGURE — Figure #5: HMM ↔ deep learning landscape** → `diagrams/lecture-21/05-hmm-vs-dl.svg`
*A 2D scatter plot. X-axis: training-set size (log scale, $10^2$–$10^9$ examples). Y-axis: prediction accuracy. Two curves: HMM (rises early, plateaus) and deep learning (rises late, exceeds HMM). The crossover is at ~$10^4$–$10^5$ examples for typical genomic-classification tasks. Annotation: "HMM regime: limited data + interpretability. DL regime: abundant data + flexibility."*

---

## Part 7 — Worked End-to-End Example (≈20 min)

### 7.1 Annotating a new bacterial genome (≈8 min)

A typical 2024 workflow:

1. **Assemble** reads to contigs (L3 / new L5).
2. **Run Prodigal** for prokaryotic gene finding. Get ~5,000 ORF predictions.
3. **Translate** ORFs to protein sequences.
4. **HMMER scan against Pfam** for domain assignment. Get ~80% of ORFs annotated.
5. **DIAMOND BLAST against UniProt** for the rest. Get function transfer for ~90% of ORFs.
6. **Map to KEGG pathways** for metabolic context.
7. **Final gene set** with ~95% functional annotation.

Total runtime: hours on a laptop.

### 7.2 Annotating a new mammalian genome (≈6 min)

For a 3 Gbp mammalian genome:

1. **Assemble** with long reads (L11 / new L15).
2. **Repeat-mask** with RepeatMasker (~40% of genome).
3. **Train AUGUSTUS** on a related species' gene set as starting point.
4. **Run BRAKER** with own RNA-seq data → trained eukaryotic gene model.
5. **Run AUGUSTUS** on the masked genome → ~20,000 gene predictions.
6. **Validate** with HMMER/Pfam for domain annotation; cross-check with InterProScan.
7. **Assign function** via ortholog identification (OrthoFinder, comparing to UniProt).

Total runtime: days on a server.

### 7.3 The 2024 frontier (≈3 min)

- **DeepBind / DeepSEA** for regulatory annotation (replace pure-HMM regulatory state inference with deep models).
- **AlphaFold-driven annotation** for protein function (predict structure → infer function from structural homologs).
- **Metagenomic gene callers** (MetaProdigal, Prokka) for community-level annotation.

The HMM-driven core (Pfam, AUGUSTUS, ChromHMM) remains the foundation; deep learning is layered on top.

### 7.4 Hands-on (≈3 min)

To do the full workflow yourself:

- Assemble 1 bacterial genome with SPAdes (~30 min on laptop).
- Run Prodigal: 1 min.
- Run hmmer search against Pfam: ~5 min.
- All open-source; tutorials available at the Pfam, HMMER, Prodigal websites.

**EMBED — Artifact #6: End-to-End Annotation Walkthrough** → `artifacts/lecture-21/06-annotation-walkthrough.html`
*Step through assembly → ORF call → Pfam → BLAST function transfer for a mock 1 kb contig. Aha: each step builds on the last; the final annotation is the synthesis.*

---

## Wrap-up (≈10 min)

### What you should take away

- **HMMs are state-space models on discrete observations.** Three algorithms: forward (likelihood), Viterbi (decoding), Baum-Welch (training). $O(T K^2)$ time complexity each.
- **Profile HMMs operationalise MSA → probabilistic family.** HMMER + Pfam are the workhorses for protein-domain annotation.
- **Gene finding is HMM-driven.** GENSCAN architecture: state classes for intergenic, promoter, exon-initial / internal / final, intron-phase 0/1/2, UTR. AUGUSTUS adds CRF + RNA-seq hints.
- **HMM applications beyond gene finding**: ChromHMM (chromatin segmentation), pair-HMMs (variant calling), splice-aware alignment, Viterbi-on-DAG (pangenome alignment).
- **Modern alternatives** (CRFs, RNNs, transformers) have replaced HMMs in regulatory prediction but not in segmentation tasks. The HMM niche: interpretable annotation, small training sets, sparse states.
- **EE framings**: HMMs as state-space models; forward-backward as belief propagation; Viterbi as MAP decoding; profile HMMs as position-specific filters with gap states; Baum-Welch as EM.

### Next lecture

Bulk RNA-seq (existing L5, becomes new L8): from reads to expression counts. The splice-aware alignment in HISAT2 you saw mentioned here is what does the heavy lifting.

### Homework

1. Build a 2-state HMM for CpG islands (states: island vs background). Run Viterbi on a synthetic 10 kb sequence; identify the islands. Compute precision-recall against ground truth.
2. Run HMMER on a protein sequence of your choice against Pfam (web service or local). Report all domain hits with $E < 10^{-3}$. For each, note the alignment quality.
3. Use AUGUSTUS to predict genes on a 100 kb genomic region of your choice (any species AUGUSTUS supports). Compare predictions to existing UCSC annotations; report concordance.
4. Implement Baum-Welch training in Python for a 3-state HMM. Train on synthetic data; verify the recovered transitions and emissions match the generating process. (~50 lines of code.)
5. Pick a chromatin state from a ChromHMM-annotated genome (e.g., ENCODE H1 ESC). Compute the genome-wide proportion of bases in that state. Compare to the published numbers.

### Recommended reading

- Eddy, S. R. (1998). Profile hidden Markov models. *Bioinformatics* 14, 755–763.
- Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology* 7, e1002195.
- Krogh, A., Brown, M., Mian, I. S., et al. (1994). Hidden Markov models in computational biology: applications to protein modeling. *Journal of Molecular Biology* 235, 1501–1531.
- Burge, C., & Karlin, S. (1997). Prediction of complete gene structures in human genomic DNA. *Journal of Molecular Biology* 268, 78–94.
- Stanke, M., & Waack, S. (2003). Gene prediction with a hidden Markov model and a new intron submodel. *Bioinformatics* 19, ii215–ii225.
- Ernst, J., & Kellis, M. (2012). ChromHMM: automating chromatin-state discovery and characterization. *Nature Methods* 9, 215–216.
- Mistry, J., Chuguransky, S., Williams, L., et al. (2021). Pfam: the protein families database in 2021. *Nucleic Acids Research* 49, D412–D419.
- Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. *Proceedings of the IEEE* 77, 257–286. (Classic tutorial; still excellent.)
- HMMER: <http://hmmer.org/>
- Pfam (now at InterPro): <https://www.ebi.ac.uk/interpro/>
- AUGUSTUS: <https://bioinf.uni-greifswald.de/augustus/>
- ChromHMM: <http://compbio.mit.edu/ChromHMM/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — HMM Foundations                                 | 30 min | 0:30 |
| Part 2 — The Three HMM Algorithms                          | 45 min | 1:15 |
| Part 3 — Profile HMMs and HMMER                            | 35 min | 1:50 |
| Part 4 — Gene Finding                                       | 40 min | 2:30 |
| Part 5 — HMM Applications Beyond Gene Finding              | 25 min | 2:55 |
| Part 6 — Modern Alternatives and the HMM Niche             | 20 min | 3:15 |
| Part 7 — Worked End-to-End Example                          | 20 min | 3:35 |
| Wrap-up                                                      | 10 min | 3:45 |

**Total:** ~3h 45min of content.
