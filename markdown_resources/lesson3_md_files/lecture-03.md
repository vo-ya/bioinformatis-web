# Lecture 3 — DNA Sequence Assembly: From Reads to Genomes

> **Duration**: ≈205 min content (~3h 25min)
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-03.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Distinguish de novo assembly from resequencing and explain when each is the right problem to solve.
2. Describe the six stages of a typical assembly pipeline from raw reads to polished scaffolds.
3. Reason about coverage — including its Poisson distribution and the consequences of low-coverage regions — for a given sequencing run.
4. Explain how k-mer error correction works and why it's a prerequisite for graph construction.
5. Construct a small de Bruijn graph by hand from short reads and read off the contigs.
6. Compare de Bruijn graphs and overlap graphs on the axes of complexity, runtime, and read-length suitability.
7. Identify the three main topological signatures of errors in an assembly graph — tips, bubbles, and tangles — and describe how assemblers resolve each.
8. Interpret N50, NG50, and related metrics, and state their known failure modes.

---

## Part 1 — The Assembly Problem (≈30 min)

### 1.1 Why assemble? (≈6 min)

Lecture 2 ended with a pile of reads aligned against a reference genome. That assumes a reference exists — and for humans, mice, *E. coli*, and a few hundred other well-studied organisms, one does. For everything else — a soil microbe nobody has sequenced before, a tumor with structural variants not in any reference, a novel viral isolate from a patient swab — there is no reference to align to. You have to reconstruct the genome from the reads themselves.

That reconstruction is **assembly**. Given millions of short reads drawn (with errors) from an unknown genome, the assembler outputs its best guess at the original sequence: a set of long strings called **contigs** that together cover as much of the genome as the data allows, with gaps wherever the data does not.

> **Intuition box**: Assembly is a jigsaw puzzle. You have millions of pieces, each about 150 pixels wide, from a picture that is three billion pixels long. The pieces are sampled randomly from overlapping regions of the picture, so each region appears on roughly thirty pieces. About one in a hundred ink strokes on each piece is wrong. There is no picture on the box. Reconstruct the picture.

### 1.2 De novo vs resequencing (≈8 min)

Two problems sound similar and are often confused.

**Resequencing** is what Lecture 2 covered. You have a reference, the sample is expected to be nearly identical to it, you align the reads to the reference, and the output is a list of differences (variants). The hard part is fast alignment; the biology is known. This is how 23andMe, clinical diagnostic panels, and most of human population genetics work.

**De novo assembly** is what this lecture covers. You do not have a reference. Or the sample differs from the reference so much that aligning would miss most of what's interesting — a cancer with large structural rearrangements, a new bacterial species, a plant with a recent whole-genome duplication. The hard part is reconstruction itself; there is no template to match against.

> **EE framing**: De novo assembly is blind reconstruction of a signal from samples — no pilot, no training sequence, no channel estimate. Resequencing is matched-filter detection against a known template. The problems are the same only when you lack the template, which is the one case that makes de novo hard. Every technique in this lecture either exploits the algebraic structure of the sample overlaps (de Bruijn graphs) or searches directly over possible assemblies (overlap graphs), because without a template there is no shortcut.

**FIGURE — Figure #1: De novo vs resequencing** → `diagrams/lecture-03/01-denovo-vs-resequencing.svg`
*Two panels side by side: the resequencing panel shows short reads aligning to a pre-drawn reference; the de novo panel shows the same reads without a reference, assembling themselves into a contig through pairwise overlaps.*

### 1.3 The typical assembly pipeline (≈10 min)

An assembly run is rarely a single program. It's a pipeline of six stages, and production assemblers like SPAdes, Flye, and Canu are orchestrators that run roughly this sequence:

1. **Input QC.** Adapter trimming, length filtering, per-read quality trimming. Same tools as Lecture 1 — `fastp`, `trim_galore`.
2. **Error correction.** Clean up systematic read errors before they poison the graph. Tools: `BFC`, `Lighter`, SPAdes's built-in BayesHammer. §2.3.
3. **Graph construction.** Build the de Bruijn graph (for short reads) or the overlap graph (for long reads). §3.
4. **Graph cleanup.** Remove tips, collapse bubbles, simplify tangles. §4.2.
5. **Contig construction.** Walk the remaining graph to extract confident linear sequences. §4.1.
6. **Scaffolding.** Use paired-end or long-range information to link contigs across gaps. §4.4.
7. (Optional) **Polishing.** Align the reads back to the assembly and correct any residual errors — usually with a different tool class, `Pilon` for short reads or `Medaka` / `Racon` for long reads.

**FIGURE — Figure #2: The assembly pipeline** → `diagrams/lecture-03/02-assembly-pipeline.svg`
*A seven-box left-to-right flow: Reads → QC → Error correction → Graph construction → Graph cleanup → Contigs → Scaffolds, with Polishing shown as an optional branch.*

Every stage can fail independently, and most assembly post-mortems are about figuring out which stage let a problem through. Parts 2 through 5 cover these stages in order.

### 1.4 Reads and read types for assembly (≈6 min)

Not all reads are equally useful for assembly. The relevant axes are length, accuracy, and whether the reads come in pairs.

- **Short accurate reads** (Illumina, 150–300 bp, ≤0.5% per-base error). Cheap, deep, low error. Bad at repeats longer than 300 bp. Dominant for bacterial and small-eukaryote assembly.
- **Long noisy reads** (PacBio CLR, Oxford Nanopore raw, 10–100 kb, 5–15% raw error). Span most repeats. Noise corrects out with consensus (§2.3). Essential for mammalian-genome assembly.
- **Long accurate reads** (PacBio HiFi, 10–25 kb, <0.1% error). The best of both worlds, at a price premium. The current gold standard for de novo mammalian assembly.
- **Paired-end / mate-pair reads**. Two reads from opposite ends of a single fragment; the distance between them is roughly known. Paired-end inserts are 300–600 bp, mate-pair ("jumping library") inserts are 2–40 kb. Paired information helps bridge short-repeat gaps; mate-pair information helps scaffold.

> **Intuition box**: Short reads are small puzzle pieces with crisp, reliable edges. Long reads are bigger pieces that cover more of the picture at once, but their edges are blurry. The cheapest modern assemblies mix both: short reads for accuracy, long reads for span.

**FIGURE — Figure #3: Read types for assembly** → `diagrams/lecture-03/03-read-types.svg`
*Four horizontal tracks on the same genomic scale bar showing typical read lengths and accuracies: Illumina (short, tight), PacBio HiFi (medium-long, tight), ONT (very long, fuzzy), and a paired-end pair linked by a dashed "insert" arc.*

---

## Part 2 — Data Before Assembly (≈30 min)

### 2.1 Coverage and average coverage (≈10 min)

Before any assembly algorithm runs, the raw reads are characterised by one summary number more than any other: **coverage**.

Coverage at a position is the number of reads that span that position. Averaged across the genome, it's the **average coverage** or **depth**:

```
average coverage = (N reads × read length) / genome length
```

For a 5 Mb bacterial genome sequenced with 100,000 Illumina 150 bp reads, average coverage is 100,000 × 150 / 5,000,000 = 3×. For a typical WGS run that value is 30× to 60×.

> **EE framing**: Coverage is Poisson sampling. If you throw N reads uniformly over an L-base genome, the depth at any single position follows a Poisson distribution with mean λ = N·|r|/L — same distribution that photon-counting detectors, spike-counting neuroscientists, and telephone-switching-network engineers have been working with for a century. At λ = 30, the probability of zero coverage at a random position is e⁻³⁰ ≈ 9 × 10⁻¹⁴ — mathematically vanishing. In practice it is not, because GC bias, PCR dropout, and repetitive regions generate real zeros that the Poisson model doesn't capture.

**FIGURE — Figure #4: Coverage distribution** → `diagrams/lecture-03/04-coverage-distribution.svg`
*Top: a horizontal reference with read pileup drawn above, coverage depth at each position plotted as a stem below. Bottom: a histogram of per-position coverage overlaid with the theoretical Poisson curve at the observed mean.*

**EMBED — Artifact #1: Coverage Simulator** → `artifacts/lecture-03/01-coverage-simulator.html`
*Sprinkle N reads of length L uniformly along a reference bar, watch the coverage depth fill in, and compare the resulting histogram to the Poisson distribution at the given mean.*

The practical rule of thumb: **average coverage is a weak summary.** It tells you nothing about uniformity. An assembly with 30× average can be unassemblable if the coverage is distributed as 0× across 10% of the genome and 33× across the rest. Always look at the coverage histogram, not just the mean.

> **Discussion prompt**: Your assembly of a new bacterial isolate shows 30× average coverage but includes several kilobase-long zero-coverage gaps. Is this a sequencing problem, a mapping problem, or a biology problem? (Could be any: GC-biased PCR dropout during library prep, regions too repetitive to map uniquely, or a real absence from your sample that's present in the reference. You can't tell without inspecting the gap locations.)

### 2.2 Uniformity and coverage dropouts (≈6 min)

The Poisson model assumes reads are drawn uniformly at random. They are not.

Three sources of non-uniformity matter in practice:

- **GC bias.** PCR amplification during library prep is more efficient on moderate-GC fragments than on very high-GC or very low-GC regions. Extreme-GC regions end up under-represented by a factor of two to ten. This is the dominant cause of short-read assembly gaps in real genomes.
- **Repeat collapse.** Reads from repetitive elements map equivocally to multiple positions. Some pipelines weight them fractionally, some drop them, some assign them randomly — and the resulting coverage at repeat positions can look artificially low.
- **Strand bias and read-end effects.** The first and last few bases of a read have elevated error rates and can be trimmed aggressively, which concentrates coverage in the middle of reads and leaves the outer edges thin.

The cleanup for GC bias is a second PCR-free library prep (or a PCR-free chemistry like 10X Genomics or PCR-free Illumina kits). The cleanup for repeats is long reads. There is no cleanup for strand bias — you tolerate it.

### 2.3 Error correction (≈14 min)

Raw short reads have per-base error rates of roughly 0.1–1%. A 150 bp read therefore has on average 0.15 to 1.5 errors. Long raw reads from ONT can have 5–15% error. Feed these reads directly into a de Bruijn graph and the graph will be littered with spurious edges — each error creates a wrong k-mer that branches off the true path.

Before graph construction, every production assembler **corrects** the reads. The standard technique is k-mer spectrum analysis.

**How k-mer error correction works, in three steps:**

1. Count every k-mer in every read. For a 5 Mb genome with 30× coverage and k = 21, you get a k-mer frequency histogram with two populations: a peak around frequency 30 (true k-mers, each seen about 30 times because of coverage) and a long tail at frequencies 1 and 2 (erroneous k-mers, each a unique typo produced by a single read error).

2. Set a threshold. K-mers with frequency below the threshold (say, 3) are presumed to be errors. K-mers above it are presumed to be true.

3. For every read, walk across the read character by character. When you hit an error k-mer, try single-base substitutions and look for the one that produces a true k-mer. Replace the offending base and move on.

**FIGURE — Figure #5: k-mer spectrum** → `diagrams/lecture-03/05-kmer-spectrum.svg`
*A log-y histogram of k-mer frequency: tall peak at frequency 1 labeled "errors", larger peak around the coverage depth labeled "true k-mers", horizontal line showing the error-threshold cutoff separating them.*

**EMBED — Artifact #2: k-mer Error Correction** → `artifacts/lecture-03/02-kmer-error-correction.html`
*Input a set of short reads with a controlled error rate; see the k-mer spectrum build up with a clear error peak; set a threshold and watch erroneous bases get corrected to the nearest true k-mer.*

> **EE framing**: K-mer error correction is signal-noise separation in the k-mer frequency domain. True k-mers sit in a Poisson peak around the coverage depth; erroneous k-mers sit in a long tail near frequency 1. The two populations are well-separated at 30× coverage; they overlap at 5× coverage, which is why low-coverage assembly is hard. Setting the threshold is the same problem as setting a detection threshold on a matched filter — you pick the false-alarm / miss-detection tradeoff that's right for your downstream noise tolerance.

Modern assemblers (SPAdes, Canu, Flye) do this automatically with BayesHammer or similar Bayesian corrector. You don't typically run it as a separate step; you benefit from it implicitly.

---

## Part 3 — Graph Theory for Assembly (≈50 min)

### 3.1 A detour through graph theory (≈7 min)

Before de Bruijn graphs, a minimum viable vocabulary for the rest of the lecture.

A **graph** is a set of **nodes** connected by **edges**. Edges can be **directed** (A → B) or **undirected** (A – B). Graphs can have **cycles** (paths that return to their starting node) or be **acyclic**. Nodes can have **labels** or **attributes**. Edges can have **weights**.

Three observations matter for what follows:

- Graphs are everywhere in EE: circuit topology, control-flow diagrams of a program, the state-transition diagram of a finite-state machine, the Markov chain of a channel, the trellis of a Viterbi decoder.
- An **Eulerian path** on a graph visits every edge exactly once. A **Hamiltonian path** visits every node exactly once. These are different problems. Eulerian paths are easy (polynomial time). Hamiltonian paths are NP-hard. This distinction determines the entire design of assembly algorithms.
- The **degree** of a node is the number of edges touching it. For directed graphs, you have **in-degree** and **out-degree** separately. Graph algorithms usually look at degree as their first signal.

> **Intuition box**: A graph is a bookkeeping system for "A is connected to B". In circuits A is a node and B is another, and the edge is a wire with some impedance. In assembly A is a short substring (a k-mer minus its last base) and B is another short substring, and the edge is the single k-mer that overlaps them. The same algorithms work on both.

### 3.2 De Bruijn graphs (≈14 min)

Here is the idea that made modern short-read assembly feasible.

Given a set of reads and a parameter k, the **de Bruijn graph** has:

- One **node** per distinct (k−1)-mer that appears as a prefix or suffix of some k-mer in the reads.
- One **directed edge** per distinct k-mer in the reads, connecting the node labeled with its first (k−1) characters to the node labeled with its last (k−1) characters.

That's it. Two-line definition. Implementations vary in how they store it, but the concept is this simple.

Example. Reads: `ACGTCC`, `CGTCCA`, `GTCCAT`. With k = 4, the k-mers are ACGT, CGTC, GTCC, TCCA, CCAT. The (k−1)-mers — nodes — are: ACG, CGT, GTC, TCC, CCA, CAT. Edges:

```
ACG --ACGT--> CGT --CGTC--> GTC --GTCC--> TCC --TCCA--> CCA --CCAT--> CAT
```

One path from ACG to CAT, walking through every edge exactly once. Reading off the node labels, keeping only the last character of each: `ACGTCCAT`. That's the genome the reads came from.

**FIGURE — Figure #6: De Bruijn graph** → `diagrams/lecture-03/06-debruijn-graph.svg`
*Three input reads at the top, k-mer decomposition in the middle showing each k-mer split into its (k−1)-mer prefix and suffix, and the assembled graph at the bottom with nodes as circles and edges labeled with the k-mers that produced them.*

**EMBED — Artifact #3: De Bruijn Graph Builder** → `artifacts/lecture-03/03-debruijn-builder.html`
*Input a set of short reads and a value of k; watch the graph build itself node by node; hover an edge to see which read and which offset it came from.*

> **EE framing**: A de Bruijn graph is the state-transition diagram of a shift register whose state is the current (k-1)-mer. Every time you slide the window forward by one base, the oldest letter shifts out and the new letter shifts in — that's the register's transition. The edges of the graph are exactly the transitions the shift register performs as it processes a string. Which means assembly is deriving a string from observations of its shift-register state transitions — exactly the problem a Viterbi decoder solves against a convolutional-code trellis. The algebraic structure is the same; only the cost function and the target representation differ.

> **Historical pointer**: The graph is named after Nicolaas de Bruijn (1946), who studied sequences with special combinatorial properties. Idury and Waterman applied it to DNA assembly in 1995. Pevzner made it the foundation of a practical assembler in 2001, and the Velvet assembler (Zerbino & Birney, 2008) turned the approach into something that ran on bacterial genomes overnight on a desktop. By 2010 every serious short-read assembler was a de Bruijn graph variant.

### 3.3 Solving de Bruijn graphs: Eulerian paths (≈12 min)

Once you have the graph, the assembly problem reduces to: **find a path that walks every edge exactly once.** This is called an Eulerian path.

Euler (the Euler) proved in 1735 — studying whether you could walk all seven bridges of Königsberg without crossing any of them twice — that an Eulerian path exists in a connected graph if and only if every node has even degree, or exactly two nodes have odd degree (and those two are the start and end of the path). For directed graphs: in-degree equals out-degree everywhere, or exactly two nodes have imbalanced degree by one (one source, one sink).

The construction, Hierholzer's algorithm, is linear-time: start at the source, walk edges until you return to your starting point (forming a cycle), then walk any remaining unvisited edges by finding an intermediate node on your current path that has unvisited edges and splicing a new cycle in there.

```
walk the graph, keeping track of unvisited edges
if stuck mid-path:
    find the most recent node with unvisited edges
    recurse into the subgraph from that node
    splice the returned subpath into the main path
```

This is O(|E|) — linear in the number of edges. For a bacterial genome at 30× coverage and k = 25, |E| ≈ 5 × 10⁶ · 30 = 1.5 × 10⁸, which is a few seconds on modern hardware.

**EMBED — Artifact #4: Eulerian Path Finder** → `artifacts/lecture-03/04-eulerian-path.html`
*Step through Hierholzer's algorithm on a small de Bruijn graph; see unvisited edges shrink to zero as the path is built; watch how the algorithm handles branching.*

> **EE framing**: Euler's 1735 result is the reason assembly is polynomial. If we had to find a Hamiltonian path (visit every node exactly once) instead, we'd be in NP-hard territory — and for a long time, people did approach assembly as a Hamiltonian-path problem over overlap graphs, which is why pre-2001 assemblers were much slower. The whole de Bruijn shift is a reformulation that turns assembly from an intractable problem into a tractable one.

### 3.4 Overlap graphs (≈10 min)

The older, more intuitive way to formalize assembly.

An **overlap graph** has:

- One **node** per read.
- One **edge** from read *a* to read *b* if the suffix of *a* matches the prefix of *b* for at least some threshold length *w*.

To assemble, find a Hamiltonian path in this graph — a path that visits every read exactly once. Concatenate the reads along the path (with overlaps merged) and you have the genome.

This is the **Overlap-Layout-Consensus (OLC)** approach used by Celera (the software that assembled the human genome for Celera Genomics in 2001) and still used by modern long-read assemblers like Canu, Flye, and miniasm.

**FIGURE — Figure #7: Overlap graph** → `diagrams/lecture-03/07-overlap-graph.svg`
*Five short reads listed on the left, a pairwise overlap table in the middle, and the resulting directed overlap graph on the right with nodes labeled by read index and edges labeled by overlap length.*

The catch is Hamiltonian: NP-hard in general. In practice, long-read assemblers exploit the fact that the graph is usually almost-linear (each read has exactly one or two high-confidence overlap neighbors), so a greedy heuristic does well. When it doesn't — at repeats — the graph diverges and the assembler reports the ambiguity rather than guessing.

### 3.5 Overlap vs de Bruijn: when to use which (≈7 min)

The choice is almost never free. Read length forces it.

| Attribute | De Bruijn graph | Overlap graph |
|---|---|---|
| Node is a | (k−1)-mer | read |
| Edge is a | k-mer (single overlap between two k-mer contexts) | read-to-read suffix-prefix overlap |
| Construction | O(total read length) | O(reads²) in the naive case, O(reads · log reads) with indexing |
| Search | Eulerian path (polynomial) | Hamiltonian path (NP-hard; heuristic in practice) |
| Memory | Scales with graph size, not read count | Scales with read count |
| Best for | Short accurate reads (Illumina) | Long reads (PacBio, ONT) |

> **Warning box**: The question "de Bruijn or overlap?" rarely has a free answer. For reads under 500 bp you effectively have to use de Bruijn — overlap graphs would be quadratic in the number of reads. For reads over 5 kb you effectively have to use overlap — de Bruijn with k ≈ 1000 has no redundancy to leverage. Hybrid assemblers like SPAdes-hybrid and Unicycler use both: de Bruijn to build a clean short-read backbone, then long-read overlaps to bridge repeats and gaps.

---

## Part 4 — Contig Construction (≈50 min)

### 4.1 From graph to contigs (≈12 min)

The graph is not the assembly. The assembly is a set of strings — **contigs** — each representing a region of the genome the assembler is confident in.

A contig is produced by walking a **non-branching path** through the graph: a maximal path where every internal node has exactly one in-edge and one out-edge. The walk terminates at a branch point (where the node has multiple out-edges or multiple in-edges) or at a dead end. Each non-branching path becomes one contig.

> **Intuition box**: A contig is a stretch of the genome where there is no ambiguity — one and only one way to proceed through the graph from start to end. Every time the graph forks, you break the contig. Forks are where the data stops being confident.

Why break at forks? Because forks are exactly the places where more than one reconstruction is consistent with the reads. The assembler could guess, but guessing hides the uncertainty. Producing two separate contigs at a fork is the honest output.

For a clean bacterial genome with 30× coverage and no repeats longer than the read length, the de Bruijn graph is close to a single linear path — one contig, the whole genome. Real genomes are not like that. Real genomes have repeats, errors, and structural variation, and the graph has forks.

### 4.2 Graph topology: tips, bubbles, and tangles (≈13 min)

Three patterns of "not-a-linear-path" show up in every assembly graph. Learn to recognize them.

- **Tips** are short dead-end branches a few edges long. They come from single sequencing errors near the end of a read — the error creates a unique k-mer that branches off the true path and then terminates because no other read extends it. Cleanup: remove any tip shorter than a threshold (typically 2·k).
- **Bubbles** are parallel paths between the same pair of branch nodes. They come from SNPs (one read is heterozygous, one isn't), from repeats of exactly the wrong length, or from small indels. Cleanup: identify the two paths, compare their coverage, keep the higher-coverage one, collapse the bubble.
- **Tangles** are densely interconnected subgraphs with many in- and out-edges. They come from repeats — every copy of a repeat collapses onto the same path in the graph, and the reads coming into the repeat can exit through any of the copies. Cleanup: no clean cleanup. Report as an unresolved region; use long reads or paired-end information to try to walk it.

**FIGURE — Figure #8: Graph topology signatures** → `diagrams/lecture-03/08-topology-signatures.svg`
*Three side-by-side panels showing a tip (a short dead-end branch labeled "single-read error"), a bubble (two parallel paths between shared endpoints labeled "SNP / small variant"), and a tangle (a dense multi-path subgraph labeled "repeat").*

**EMBED — Artifact #5: Assembly Topology Inspector** → `artifacts/lecture-03/05-topology-inspector.html`
*Load a canned small-genome assembly graph with tips, bubbles, and tangles highlighted. Apply each cleanup operation and watch the graph simplify toward linear segments.*

> **EE framing**: Topology in the graph is the signal pattern of different error types. A sequencing error shows up as a tip one or two edges long. A biological SNP shows up as a bubble exactly two paths wide. A repeat of length R and copy number N shows up as a tangle with N entry and N exit edges. Each has a distinct signature, and the assembler's cleanup stage is essentially pattern matching followed by a fixed remediation — matched filter, then take the action the match implies.

### 4.3 Repeats and ambiguity (≈13 min)

Repeats are the fundamental obstacle to assembly. Everything else is engineering; repeats are physics.

A **repeat** is a sequence that appears in the genome more than once, at different positions. Examples in the human genome: Alu elements (300 bp, ~1 million copies, ~10% of the genome), LINEs (6 kb, ~500,000 copies, ~20% of the genome), centromeric satellite repeats (arrays of short motifs stretching for megabases), segmental duplications (large blocks of 10 kb or more with >90% identity).

Here's why they break assembly. Take two identical copies of a 500 bp repeat, located at positions X and Y in the genome. For any k ≤ 500, every k-mer in the repeat appears in both copies. In the de Bruijn graph, the two copies **collapse into a single path** — the graph has one path for the repeat, with the two upstream-of-the-repeat regions converging into it and the two downstream-of-the-repeat regions diverging out of it. The assembler can tell that the repeat is there, but it can't tell which upstream belongs to which downstream.

**FIGURE — Figure #9: Repeat collapse** → `diagrams/lecture-03/09-repeat-collapse.svg`
*Top: a schematic genome with two identical repeat copies (labeled R) at positions X and Y, surrounded by unique flanking sequence. Bottom: the de Bruijn graph showing the two flanking regions converging into a single R node and diverging out of it, with the two possible reconstructions marked as ambiguous.*

> **Warning box**: **Repeats are the fundamental limit of assembly.** Any repeat longer than the read length cannot be resolved from the reads alone — the reads don't span it, so there is no information to distinguish the copies. Long reads help (if the read spans the repeat, the repeat is resolved). Mate-pair reads help (if the insert spans the repeat, scaffolding can bridge it). But no algorithm resolves all repeats. Published genomes with "complete" in the title usually have gaps inside centromeric or telomeric repeat arrays.

Assemblers handle repeats in three ways, in increasing order of honesty:

1. **Collapse and report.** The repeat appears once in the assembly with elevated coverage. This is what short-read-only assemblers do for long repeats. The copy number is lost.
2. **Fork and report both.** The assembler emits two contigs at the repeat boundary and marks them as alternate paths. This requires reasoning about coverage.
3. **Span with long reads.** A read that is longer than the repeat unambiguously resolves the copy. This is why long-read assembly is necessary for any mammalian genome.

### 4.4 Gaps and scaffolding (≈12 min)

After contig construction and graph cleanup, the assembler has a set of contigs. These contigs are ordered and oriented with respect to the genome **only within the limits of the read data**. Unresolved repeats, coverage gaps, and low-information regions all break contig contiguity.

**Scaffolding** is the post-step of ordering and orienting contigs across the unresolved regions, using information that was ignored during graph construction.

The two workhorses:

- **Paired-end reads.** A pair of reads from the same ~500 bp fragment. If one read maps to contig A and the other to contig B, you know A and B are within ~500 bp of each other and you can order them.
- **Mate-pair libraries.** Same idea, larger inserts — 2 kb to 40 kb. The long inserts jump over short repeats that paired-ends can't bridge.

Modern equivalents include Hi-C (chromosomal-conformation capture, which links contigs that are near each other in 3D even across megabases) and optical maps (Bionano) which provide long-range order information orthogonal to the reads.

The output of scaffolding is a **scaffold**: a sequence consisting of contigs linked by gaps. Gaps are filled with N characters whose length is the estimated gap size from the linking evidence. A scaffold of length 2 Mb might be 1.8 Mb of real sequence and 0.2 Mb of N's distributed across a dozen gaps.

**FIGURE — Figure #10: Scaffolding** → `diagrams/lecture-03/10-scaffolding-mate-pairs.svg`
*Three contigs drawn as horizontal bars at the top; below, several mate-pair reads with arcs linking one end in contig A to the other end in contig B, another set linking B to C; the resulting scaffold at the bottom shows contigs A, B, C placed in order with Ns between them.*

> **Historical pointer**: Mate-pair ("jumping") libraries were introduced by the Celera Genomics human genome assembly (Myers, Sutton et al., 2000) precisely to solve the repeat problem. The whole-genome shotgun approach without mate-pairs would have produced a million fragments no-one could order. With 2-, 10-, and 50-kb insert libraries layered together, Celera delivered an assembly competitive with the HGP's clone-by-clone approach at one-fifth the cost.

> **Discussion prompt**: Why don't assemblers just paper over every gap by inserting some guess for what's in it, so the output is a continuous sequence? (Answers: because the gap N's carry no signal — they are a placeholder that downstream tools can recognize and skip over. A guess would silently corrupt variant calling, expression quantification, and any other analysis that maps onto the assembly. Gaps are honest about what the data doesn't resolve. Hiding gaps would produce nicer-looking assemblies that give wrong answers.)

---

## Part 5 — Whole Genome Shotgun and Quality Control (≈35 min)

### 5.1 De novo whole-genome shotgun assembly (≈12 min)

Putting all of Parts 1–4 together: here's how a bacterial genome gets assembled end-to-end.

1. Extract DNA; fragment it to 300–500 bp; build a sequencing library; sequence on Illumina to get ~30× of 150 bp reads. Cost: roughly $100 for a small bacterial genome in 2024.
2. Run `fastp` on the reads for adapter trimming and quality filtering. Drop reads shorter than 50 bp.
3. Run SPAdes with `--isolate` (or `--careful`) mode. SPAdes internally runs BayesHammer error correction, builds de Bruijn graphs at several k values, merges them, cleans tips and bubbles, emits contigs.
4. (Optional) Sequence a second library of mate-pairs or do a PacBio HiFi run for scaffolding.
5. Run a scaffolder (SSPACE, BESST, or SPAdes's own) to order and orient contigs.
6. Run `Pilon` or `Racon` to polish any remaining single-base errors.
7. Run `QUAST` (§5.2) to compute metrics.

Total compute for a 5 Mb bacterial genome: roughly 30 minutes on a modern laptop. The same pipeline on a 3 Gb mammalian genome requires 100× the compute, multi-day runtime, and long reads throughout; typical modern runs use PacBio HiFi with the HiCanu or hifiasm assembler.

> **Historical pointer**: The shotgun-versus-clone-based debate that drove the 1998–2001 human-genome race is over: shotgun won. The hierarchical-clone approach used by the public HGP produced a slightly higher-quality reference but was ten times more expensive. Celera's shotgun approach — the technique described in this section — became the universal standard for every genome sequenced since, from microbes to mammals.

### 5.2 Assembly metrics: N50, NGx, consistency, coverage (≈15 min)

How do you tell a good assembly from a bad one?

**N50** is the dominant summary. Definition: sort contigs by length, largest first; walk down the list adding up lengths; N50 is the length of the contig at which the running total first exceeds half of the total assembly length. Equivalently: the length L such that contigs of length ≥ L cover at least 50% of the total assembly.

> **EE framing**: N50 is the median of the contig-length distribution weighted by length. The usual median would be the "middle contig by count"; N50 is the "middle contig by base." That weighting is what makes the metric biologically meaningful — you care more about where half the genome sits, not about where the middle contig out of ten thousand sits.

**NG50** is the same idea but using expected genome size as the denominator instead of total assembly length. Useful because you can compare NG50 across different assemblers and different assemblies of the same genome.

**Nx for other x.** N90 is the length at which 90% of the assembly is covered. N10 is the length at which 10% of the assembly is covered. The full curve of Nx as x varies from 0 to 100 is more informative than any single N50.

**FIGURE — Figure #11: N50 metric** → `diagrams/lecture-03/11-n50-metric.svg`
*A sorted-by-length horizontal bar chart of contigs on the left; a cumulative-length curve on the right showing the running total crossing 50% of the genome size at the N50 contig; N50 value labeled on both.*

**EMBED — Artifact #6: Assembly Metrics Calculator** → `artifacts/lecture-03/06-n50-calculator.html`
*Enter a list of contig lengths. See N50, N90, NG50 computed; plot the Nx curve; compare two assemblies side by side to see how different length distributions produce the same N50 with different shapes.*

> **Warning box**: N50 is easy to game. Aggressive joining through uncertain regions increases N50 while decreasing accuracy. Always report N50 alongside **total length**, **number of contigs**, **largest contig**, and **genome fraction covered**. A tight 5-number summary tells you more than any single metric.

**Other metrics to know:**

- **Assembly consistency** is the fraction of the assembly that aligns correctly to a reference (when one exists for comparison). Requires a trusted reference.
- **Genome coverage** or **completeness** is the fraction of the expected genome that the assembly covers. For a bacterial genome of 5 Mb with total assembly 4.8 Mb, completeness is 96%. Tools like BUSCO use expected single-copy orthologs to estimate completeness without needing the full reference.
- **Misassembly count** is the number of places the assembly disagrees with a reference in ways not explainable by real variation — inverted segments, translocated segments, spurious junctions. Computed by tools like QUAST.

### 5.3 BAM format for aligned reads (≈8 min)

Once you have an assembly, the first thing you do is align the reads back to it — to polish, to estimate coverage, to detect residual errors, to check for sample contamination. Those alignments are stored in **BAM** format.

BAM is the binary, compressed version of SAM (Sequence Alignment/Map format). One BAM file stores the alignments of an entire sequencing run. Typical structure:

- **Header** (marked with `@` lines if you open it as SAM text): reference sequences and their lengths, program that produced the file, read group information.
- **Alignment records**: one per read, with 11 mandatory fields — read name, flags (paired? unmapped? secondary?), reference name, position, MAPQ, CIGAR, mate reference, mate position, insert size, sequence, quality — plus optional tags.

A BAM file has the same alignment as a SAM file but stored in BGZF-compressed binary — roughly 3–4× smaller than gzipped SAM, and indexable so you can seek to any genomic region in constant time.

For assembly work, the most relevant BAM tools are:

- `samtools view` to read and filter.
- `samtools depth` to compute per-position coverage.
- `samtools flagstat` to count reads in each alignment category.
- `samtools sort` and `samtools index` to prepare a BAM for fast region access.
- `bcftools mpileup` or `DeepVariant` to call variants from the pileup.

The most common newcomer mistake:

> **Warning box**: BAM files from an assembler pipeline are not interchangeable with BAM files from a resequencing pipeline. An assembly BAM has reads aligned to the **assembly itself** (the contigs), not to a reference genome — so the "reference" column lists contig names like `NODE_1_length_12345`, not chromosome names like `chr3`. Tools that assume reference-aligned BAMs will silently produce nonsense output on an assembly BAM. Always confirm what the reference was before using a BAM downstream.

---

## Wrap-up (≈10 min)

### What you should take away

- **De novo assembly is template-free reconstruction.** When no reference is available, assembly is the only way to get a sequence, and the problem is fundamentally harder than resequencing.
- **The de Bruijn graph is the central algorithmic idea of short-read assembly.** It turns reconstruction into an Eulerian-path problem, which is polynomial-time and scales to mammalian genomes. Every major short-read assembler since 2008 uses it.
- **Repeats are the limit.** No algorithm fully resolves repeats longer than the read length. Long reads, paired-end reads, and Hi-C all help; none of them solve it. Assembly-with-gaps is the honest output.
- **Coverage is a Poisson process with systematic non-uniformity.** Mean is a weak summary; the distribution — and especially its left tail at low coverage — tells you where the assembly is going to fail.
- **N50 is median contig length, weighted by length.** Useful but gameable. Always pair it with total length, contig count, and genome fraction.

### Next lecture

Variant calling: given aligned reads on a reference (from Lecture 2) or an assembly (this lecture), identify where the sample differs. SNPs, indels, structural variants. The statistical machinery of Bayesian genotype callers.

### Homework

1. Download a small bacterial genome (*E. coli* K-12, ~4.6 Mb). Simulate 30× Illumina-style reads at 150 bp with `wgsim` (default error rate).
2. Assemble the simulated reads with SPAdes in isolate mode. Report: number of contigs, largest contig, N50, NG50, total assembly length. Compare to the true genome length.
3. By hand, construct the de Bruijn graph (k = 4) from these three reads: `ACAGACGT`, `CAGACGTA`, `AGACGTAC`. Draw the graph. Find the Eulerian path. Read off the contig.
4. For a hypothetical genome with contig lengths [4500, 3200, 2800, 2100, 1500, 900, 700, 600, 400, 300], compute N50, NG50 (assume genome size = 18 kb), and plot the Nx curve.
5. Explain in one paragraph why a perfect tandem repeat of exactly (k − 1) base pairs will NOT collapse in a de Bruijn graph with parameter k. (Hint: think about what happens to the k-mer at the boundary.)

### Recommended reading

- Compeau, P. E. C., Pevzner, P. A., & Tesler, G. (2011). How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology* 29, 987–991.
- Idury, R. M., & Waterman, M. S. (1995). A new algorithm for DNA sequence assembly. *Journal of Computational Biology* 2, 291–306.
- Myers, E. W., Sutton, G. G., Delcher, A. L., et al. (2000). A whole-genome assembly of Drosophila. *Science* 287, 2196–2204. (The Celera paper.)
- Bankevich, A., Nurk, S., Antipov, D., et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology* 19, 455–477.
- Koren, S., Walenz, B. P., Berlin, K., et al. (2017). Canu: scalable and accurate long-read assembly. *Genome Research* 27, 722–736.
- SAM/BAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf

---

## Appendix — Speaker cues and timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — The Assembly Problem | 30 min | 0:30 |
| Part 2 — Data Before Assembly | 30 min | 1:00 |
| Part 3 — Graph Theory for Assembly | 50 min | 1:50 |
| Part 4 — Contig Construction | 50 min | 2:40 |
| Part 5 — WGS and Quality Control | 35 min | 3:15 |
| Wrap-up | 10 min | 3:25 |

**Total**: ~3h 25min.
