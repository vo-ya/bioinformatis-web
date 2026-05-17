#import "../theme/book-theme.typ": *

= #idx("BLAST")BLAST and Sequence Search Statistics <ch:blast>

#matters[
  More than any other tool, BLAST is what working bioinformaticians
  actually run. The NCBI web service alone serves a few million queries
  per day; every newly assembled genome on the planet is annotated by
  first asking BLAST what each predicted protein looks most like. The
  reason a 1990 piece of software survived three rewrites of the field's
  hardware and two generations of replacements is not nostalgia. BLAST
  combines two ideas — heuristic #idx("seed-and-extend")seed-and-extend search and the
  Karlin–Altschul extreme-value statistics — into a tool that gives up
  exact optimality and recovers an honest #idx("E-value")E-value in exchange. That
  trade is the design move worth internalising. Once you understand it,
  the modern replacements (#idx("DIAMOND")DIAMOND, #idx("MMseqs2")MMseqs2) are variations on the same
  theme, and the application patterns — function transfer, #idx("ortholog")ortholog
  identification, off-target prediction — fall out as direct
  consequences of how the score and the E-value are defined.
]

You will spend a meaningful fraction of your bioinformatics career
running database searches. The shape of the problem is simple and
recurring: you have one sequence — a predicted gene from a newly
assembled microbial genome, a peptide from a mass-spectrometry hit, a
candidate effector cloned out of an interesting plant — and you want to
know what it is. The honest answer almost never comes from inspecting
the sequence directly. It comes from finding the closest known
relatives in a database of a few hundred million annotated proteins
and reading the answer off their labels. That operation is the
database search, and BLAST is what people mean when they say "I
BLAST-ed it."

This chapter walks the algorithm and the statistics in equal measure.
The first half is the engineering: how the seed-and-extend cascade
takes a query that would cost four months of CPU time to align exactly
against the database and reduces it to thirty seconds without losing
any of the hits you would actually have cared about. The second half
is the statistics: how Karlin and Altschul's 1990 result turned the
output of a heuristic into a calibrated false-discovery measurement,
and how the #idx("bit score")bit score and the E-value separate _what_ you found from
_how surprised you should be_. The last sections cover #idx("PSI-BLAST")PSI-BLAST's
iterated profile refinement, the modern fast successors, and the
application patterns that make BLAST a verb.

Read-alignment in Chapter 2 set up one half of the picture: many short
queries against one fixed reference. Database search is the dual: one
or a few long queries against many unknown references. The data
structure flips, the right metric changes, and the statistics
genuinely matter for the first time. Everything else follows.


== From Read Mapping to Database Search <sec:read-vs-search>

Read mapping and database search look superficially similar — both
are "find approximately matching substrings" problems — and they are
solved by completely different algorithm classes for reasons worth
naming explicitly. The differences are scale, metric, and statistics,
in that order.

In read mapping (Chapter 2), the geometry is: hundreds of millions of
queries, each one short and clean, against a single reference genome
of a few billion bases. The reference is fixed across the whole run;
the cost-effective move is to build an FM-index of the reference once
and reuse it for every read. The metric is edit distance — Hamming
plus a handful of small indels — because the queries are conspecific
and the only differences from the reference are sequencing errors and
real variants. Statistical significance is implicit: a 150-base read
that matches a unique reference location within two mismatches is
unambiguously placed, and the mapping quality MAPQ encodes the
remaining ambiguity.

Database search inverts each of those three setups. The query side is
small: usually one sequence at a time, or a small batch, on the order
of hundreds to thousands of amino acids each. The database side is
enormous and structurally heterogeneous: NCBI's `nr` collection holds
roughly $10^9$ amino-acid residues across hundreds of millions of
protein records pulled from every sequenced organism. The metric is
not edit distance. Distantly related proteins routinely share less
than 25 % amino-acid identity while folding into nearly the same 3-D
structure and performing the same biochemistry; an edit-distance score
gives them nothing. You need a *substitution-aware* score that
recognises that an aspartate and a glutamate (both small,
negatively-charged, often interchangeable) are almost as good a match
as two aspartates, while an aspartate and a tryptophan are not. And
because the database contains many sequences that share no
evolutionary relationship with the query, statistical significance
becomes the central question: out of the millions of database records
your search just walked past, _which one_ scored highly enough that
you should believe the match is not random?

#note[
  The number to keep in mind is six orders of magnitude. A naive
  exact Smith–Waterman alignment of a 1000-residue query against
  $10^{11}$ database residues is roughly $10^{14}$ cell updates. At
  $10^9$ updates per second on a modern core that is a hundred
  thousand seconds — about a day — per single query. The web service
  does that in three seconds. Almost all of the speedup is the
  heuristic that the rest of this chapter walks.
]

The algorithmic consequence is that exhaustive #idx("dynamic programming")dynamic programming is
out. BLAST's design move — the move that fixed everything for
thirty-five years and counting — was to observe that the vast majority
of database sequences share no meaningful similarity with any given
query, and that the small fraction that _do_ share similarity always
contain at least one short region where the local similarity is strong
enough to detect by hashing. So the algorithm runs a cheap filter
first to find candidate regions, and the expensive full alignment runs
only on the survivors. Three stages — *seed, extend, score* — are the
backbone. Seeds are short word-level matches found by hashing.
Extension grows seeds into local alignments under a substitution
matrix and an #idx("X-drop")X-drop stopping rule. The final score is converted to a
bit score and then to an E-value that tells you how surprised to be.

#figure(
  image("../../diagrams/lecture-19/01-blast-pipeline.svg", width: 92%),
  caption: [
    The BLAST pipeline. Five stages cascade from the raw query through
    hashing-based seed detection, #idx("two-hit")two-hit gating, gapless then gapped
    extension, and an E-value-based significance filter. Each stage
    is cheap enough to apply to the whole survivor set of the previous
    stage; only the expensive Smith–Waterman runs on the few hundred
    candidates left after gating.
  ],
) <fig:ch19-pipeline>

The detection-theoretic reading of this pipeline is worth pausing on
because it sets up everything that follows. Stage one — the
hash-based seed scan — is a coarse #idx("matched filter")matched filter applied to every
$k$-mer in the database. Stage two — the X-drop ungapped extension —
is a finer matched filter applied to the small subset of positions
that survived the seed scan. Stage three — the banded gapped
alignment — is the full Smith–Waterman #idx("likelihood ratio test")likelihood ratio test, run
only on the survivors of the ungapped extension. The final E-value is
a multiple-testing correction: it is the expected number of false
positives in a database of size $N$ at the given bit-score threshold.
A reader who has internalised Neyman–Pearson detection from a signal
processing class can read every BLAST output as the output of a
two-stage detector cascade with a built-in Bonferroni-flavoured
correction.

#note[
  BLAST gives up exact optimality and recovers statistical
  interpretability. The trade is one-directional: you cannot in
  general recover the Smith–Waterman-optimal alignment from a BLAST
  hit (the heuristic may have missed the global optimum's seed). But
  you can read off, in calibrated units, _how unlikely it would be
  to see this alignment by chance_ — and for almost every actual
  bioinformatics task that is the question you care about.
]


== Substitution Matrices and Bit Scores <sec:matrices>

The first half of BLAST's substantive content lives in its scoring
matrix. Aligning two #idx("nucleotide")nucleotide sequences is a relatively trivial
business: A matches A, G matches G, and mismatches and gaps incur fixed
penalties. Aligning two amino-acid sequences is not, because the twenty
amino acids are not exchangeable. Some pairs (D ↔ E, K ↔ R, I ↔ L)
substitute for each other constantly across evolution because their
physico-chemical properties are nearly identical; other pairs (W ↔ A,
C ↔ E) are essentially never seen in homologous positions. A useful
score has to encode this asymmetry.

=== Why Protein Dominates Distant Search

For the rest of this chapter, *protein BLAST (`blastp`)* is the default
program in mind, and the lecture's logic generalises directly to
translated-nucleotide searches (`blastx`, `tblastn`, `tblastx`) where
the database, the query, or both are translated through the genetic
code before scoring. Nucleotide-only searches (`blastn`) exist and
have their place — finding near-identical genomic regions, or
amplicon clustering at the #idx("OTU")OTU level — but for any homology question
older than tens of millions of years, the protein search wins.

Two reasons. The first is the #idx("genetic code")genetic code's degeneracy: the third
position of most codons can change without changing the encoded amino
acid, so synonymous mutations accumulate as noise on the nucleotide
track while the protein track stays still. Aligning at the protein
level discards that noise for free. The second reason is the
substitution structure itself: at the protein level, conservative
substitutions (D ↔ E) carry positive evidence of homology even where
the two sequences differ; at the nucleotide level, every difference is
penalised the same way. The practical consequence is that you can
detect homology at 25 % protein identity that is invisible at the
nucleotide level beneath the synonymous-substitution floor.

=== #idx("PAM")PAM Matrices (Dayhoff, 1978)

Margaret Dayhoff's *Point Accepted Mutation* matrices, published in
the _Atlas of Protein Sequence and Structure_ in 1978, were the first
substantial empirical substitution model. Dayhoff's group collected
about 1500 closely-related protein pairs (sequences over 85 % identical
across about 70 protein families) and tabulated the observed
substitutions at homologous positions. From the resulting counts she
built a single-step substitution probability matrix that represents
the changes expected after 1 % of residues have been replaced — one
*PAM unit* of evolutionary distance.

The cleverness was the extrapolation. To model distant evolution,
Dayhoff raised the PAM1 matrix to integer powers. PAM250 is PAM1 raised
to the 250th power, giving the substitution probabilities expected after 250 % expected
substitutions per residue — far more than 100 %, because the same
position is allowed to mutate multiple times. The matrix exponentiation
is the move that makes PAM workable for distant searches even though
the underlying data was collected on close pairs only.

Each PAM matrix is converted to a log-odds score before use:

$ s_(i j) = log_2 frac(p_(i j), f_i dot f_j) $

where $p_(i j)$ is the observed joint substitution frequency between
residues $i$ and $j$ and $f_i, f_j$ are their background frequencies.
The score is positive when residues $i$ and $j$ co-occur in homologous
positions more often than chance — that is, when the substitution
carries evidence of common ancestry. PAM30 and PAM70 are tuned for
short, closely-related searches (e.g., peptide matching). PAM250 is
the classical default for distant comparisons.

=== #idx("BLOSUM")BLOSUM Matrices (Henikoff & Henikoff, 1992)

The BLOSUM family supersedes PAM for almost every modern application,
and the reason is the data. Steven and Jorja Henikoff in 1992
constructed their matrices directly from the *BLOCKS* database — a
hand-curated collection of ungapped, conserved blocks across hundreds
of distantly related protein families. Where Dayhoff's input was
closely related pairs and her extrapolation made the matrix work at
distance, the Henikoffs' input _was_ distant alignments, so no
extrapolation was needed.

The naming convention encodes the input cutoff. BLOSUM$n$ is built
from blocks in which sequences sharing more than $n$ % identity have
been clustered down to a single representative — so #idx("BLOSUM62")BLOSUM62 sees
substitutions across pairs that are at most 62 % identical. Smaller
numbers correspond to more distant pairs in the training data:
BLOSUM45 for very distant homology, BLOSUM90 for close homology,
BLOSUM62 for the broad middle that captures most realistic searches.
NCBI BLAST's default is BLOSUM62 and has been since the mid-1990s.

#figure(
  image("../../diagrams/lecture-19/06-blosum62.svg", width: 95%),
  caption: [
    The BLOSUM62 #idx("substitution matrix")substitution matrix. Self-substitution scores along
    the diagonal are highest; chemically similar pairs (D ↔ E, K ↔ R,
    I ↔ L) score positive; chemically distant pairs (W ↔ A) score
    sharply negative. The row-column ordering groups amino acids by
    chemistry so the positive blocks cluster visibly.
  ],
) <fig:blosum62>

A useful piece of historical accounting: BLOSUM62 _was_ in fact built
on a known bug in the original 1992 implementation, which produced
slightly different entries from what the published algorithm
specified. Styczynski and colleagues noticed in 2008. By that point
BLOSUM62 had been in production for sixteen years, was tuned around in
every benchmark, and an order of magnitude of empirical work would
have been invalidated by quietly fixing it. So nobody did. The
BLOSUM62 that ships with NCBI BLAST today is the buggy 1992 version,
and that is the matrix every benchmark in modern use is calibrated
against. The lesson is the same one Chapter 14 named in passing: in
established pipelines, _stability_ of the numerical answer beats
correctness of the recipe.

#figure(
  image("../../diagrams/lecture-19/07-pam-blosum.svg", width: 95%),
  caption: [
    Choosing a substitution matrix by the expected homology range.
    BLOSUM62 spans the broad middle and is the universal default;
    PAM250 and BLOSUM45 are the distant-homology choices; PAM30 and
    BLOSUM90 are tuned for very close homologs, including peptide
    matching.
  ],
) <fig:pam-blosum>

#tip[
  In practice you almost never have a reason to change away from
  BLOSUM62. The exceptions are short peptides (BLOSUM90 or PAM30
  gives sharper discrimination on a small number of residues) and
  searches deliberately targeting deep homology where you have
  already established that BLOSUM45 finds hits BLOSUM62 misses. If
  you find yourself fiddling with the matrix to make a borderline
  hit significant, the bit score is telling you the hit is
  borderline. The matrix is not the variable to adjust.
]

=== Bit Scores

A raw BLAST score $S$ — the sum of substitution-matrix entries along
the alignment, minus gap penalties — depends on which matrix was
used. A score of 100 under BLOSUM62 is not the same evidence as a
score of 100 under PAM250, because the matrices have different
underlying scale factors. To compare scores across runs and matrices,
BLAST converts the raw score to a *bit score*:

$ S' = frac(lambda S - log K, log 2) $

where $lambda$ and $K$ are matrix-and-gap-penalty-specific constants
that come out of the Karlin–Altschul derivation in
@sec:karlin-altschul. The $lambda$ rescales the raw score to a
natural log-odds unit; the $log K$ subtracts a normalisation; the
divide-by-$log 2$ converts natural logs to bits. The end product is
a matrix-independent score with a clean log-likelihood-ratio
interpretation: *a bit score of $S'$ means the alignment is
$2^(S')$ times more likely under the homology model than under the
random null.*

Bit scores are additive in evidence. Each extra bit doubles the
likelihood ratio. A bit score of 30 says the alignment is a billion
times more likely under homology than under random. That sounds
overwhelming until you remember that the database may contain a
billion places where chance alone could produce that score — which is
exactly what the E-value, in @sec:karlin-altschul, will fix.

=== Gap Penalties

Real homologous alignments are almost never gap-free. Insertions and
deletions accumulate alongside substitutions over evolutionary time,
and any practical scoring scheme has to handle them. BLAST uses an
*#idx("affine gap")affine gap penalty*

$ "gap"(n) = -G_o - (n - 1) dot G_e $

where $G_o$ is the gap-opening cost, $G_e$ is the gap-extension cost
per additional residue, and $n$ is the gap length. The affine shape
captures the empirical observation that gap openings are evolutionarily
expensive (typically requiring a structural rearrangement) while
extending an existing gap by one more residue is cheap (a small
local #idx("indel")indel slides easily along an existing loop). NCBI BLAST
defaults for BLOSUM62 are $G_o = 11$, $G_e = 1$.

#note[
  The affine model is the simplest of a family. Linear gap penalties
  ($G_e$ only, $G_o = 0$) over-penalise long gaps; log-affine and
  double-affine variants exist for tasks with very high indel
  densities (e.g., signal peptides). For protein homology search,
  affine with the BLAST defaults is what every benchmark is built on.
]


== The Seed-and-Extend Cascade <sec:seed-extend>

The pipeline in @fig:ch19-pipeline is implemented as a strict cascade where
each stage filters the survivors of the previous one. The early
stages are cheap and have high false-positive rates; the late stages
are expensive and only see candidates that have already cleared the
filters. The discipline of throwing away as much as possible as early
as possible is what makes the search tractable.

=== Word Lists and Neighbourhoods

The cascade starts by decomposing the query into all overlapping
$k$-mers. The default $k$ for BLASTP is 3; for BLASTN it is 11 (the
nucleotide alphabet is smaller, so longer words are needed to keep
the false-positive rate at the seed stage under control). A
300-residue protein query has 298 overlapping 3-mers.

For each query 3-mer, BLAST does something a casual reader might
miss: it expands the 3-mer into a *neighbourhood* of related 3-mers
whose substitution-matrix score with the query 3-mer is at or above a
threshold $T$. So the query 3-mer `KVL` does not just generate a
single hash lookup. Under BLOSUM62 at $T = 11$ it generates roughly
50 hash lookups — `KVL` itself, plus chemically similar variants
like `KVI`, `RVL`, `KVM`, `RVI`, and so on, each of which would still
produce a positive substitution score against the original. The
neighbourhood is the move that gives BLAST sensitivity to inexact
matches at the seed stage. Without it, a query 3-mer would only
seed against an _identical_ database 3-mer, and the substitution
structure that justifies BLOSUM62 in the first place would never
get used at the cheap filter.

The neighbourhood size is the sensitivity-speed knob. Lower $T$ makes
the neighbourhood larger — more seed hits, more candidates to extend,
higher recall but slower search. Higher $T$ shrinks the neighbourhood
toward the identical-only case and runs faster at the cost of recall
on distant homologs. Default $T$ values are tuned for the typical
BLOSUM62 query and rarely changed.

#figure(
  image("../../diagrams/lecture-19/02-seed-extend.svg", width: 95%),
  caption: [
    Seed–extend in one panel. The query 3-mer `KVL` hits the database
    subject at position 47; its substitution-matrix neighbourhood
    (`KVI`, `RVL`, `KVM`) also produces hits in the same neighbourhood.
    Ungapped extension grows the alignment in both directions under
    the X-drop rule, producing a final HSP bounded by where the
    running score first drops $X$ below its current maximum.
  ],
) <fig:seed-extend>

=== The Two-Hit Method

Early BLAST versions (Altschul et al., 1990) extended every seed match.
The seed-stage false-positive rate was high enough that extension
dominated the runtime, and most extensions led nowhere. The
breakthrough in BLAST 2.0 (Altschul et al., 1997, _Nucleic Acids
Research_) was the *two-hit method*: a database location triggers
extension only if it contains _two_ non-overlapping seed hits within
distance $A$ on the same database sequence, on the same diagonal of
the alignment matrix.

The intuition is statistical. A genuine homologous region typically
contains several short conserved patches close together — multiple
3-mers above threshold within a few dozen residues of each other.
Random database positions, by contrast, contain at most one
above-threshold seed by chance; the probability of two close seeds
under the null is much lower than the probability of one. Gating on
two close seeds drops extension cost by roughly an order of magnitude
while preserving sensitivity on real homologs, which still produce
the close-seed pattern. The default spacing constraint is $A = 40$
residues for BLASTP.

#figure(
  image("../../diagrams/lecture-19/08-two-hit.svg", width: 95%),
  caption: [
    The two-hit gating method. Extension fires only when two seed
    hits on the same diagonal lie within $A = 40$ residues of each
    other. Isolated seeds — overwhelmingly false positives at the
    BLAST sensitivity setting — are filtered out without running
    extension. The two-hit constraint cuts extension cost roughly
    tenfold while preserving sensitivity on true homologs.
  ],
) <fig:two-hit>

=== Ungapped Extension and the X-Drop Rule

Each two-hit-paired seed is extended in both directions from the
seed centre, summing substitution-matrix scores without admitting
gaps. The extension continues until the running score drops more than
$X$ below the maximum score reached so far on this extension. When
that happens, the extension stops, the current alignment is trimmed
back to the maximum-score boundary, and the result is recorded as a
*High-scoring Segment Pair (HSP)*: a maximal ungapped local
alignment.

The X-drop rule is the matched-filter analog of "if the running
cross-correlation has fallen this far below the peak, you have left
the peak; stop." Setting $X$ too low truncates real HSPs in the
middle of a homologous region where one or two consecutive
mismatches happen to push the score down temporarily. Setting $X$
too high lets the extension wander into non-homologous flanking
sequence and inflates the HSP boundary noisily. The BLAST defaults
($X = 22$ raw score units for the ungapped step under BLOSUM62) are
calibrated to catch real homologs with a typical bit-score
#idx("trajectory")trajectory while rejecting noise.

#figure(
  image("../../diagrams/lecture-19/09-x-drop.svg", width: 92%),
  caption: [
    Three score trajectories during ungapped extension. The
    true-homology trajectory climbs steadily and never triggers
    X-drop; the mid-strength trajectory builds a moderate HSP and
    then triggers X-drop on the noisy flank; the noise trajectory
    triggers X-drop immediately. The X-drop rule is the matched-filter
    decision to stop when the running score has fallen far enough
    below its peak that further extension cannot recover.
  ],
) <fig:xdrop>

=== Gapped Extension

HSPs that survive the ungapped extension and clear a per-HSP score
threshold are passed to the *gapped extension* stage. This is the
only place in the pipeline where full Smith–Waterman dynamic
programming runs, and it is run in a narrow band around the
ungapped HSP rather than across the whole sequence pair. The bandwidth
($\pm 10$ residues from the HSP diagonal is typical) is wide enough
to capture realistic insertions and deletions and narrow enough that
the dynamic programming cost is linear in the HSP length rather than
quadratic in the sequence lengths. Modern BLAST runs gapped extension
only on the single best HSP per query-database pair.

=== Composition-Based Score Adjustment

The last refinement in the extension pipeline is *composition-based
statistics*, introduced by Schäffer and colleagues in 2001. Real
protein databases contain many sequences with locally unusual amino-acid
composition: signal peptides rich in hydrophobic residues, low-complexity
linkers full of glycine and serine, repetitive motifs in structural
proteins. The BLOSUM62 matrix was calibrated on average protein
composition; when applied to a glycine-rich query against a
glycine-rich subject, the matrix systematically overestimates the
score because the background frequencies $f_a$ in the log-odds
denominator no longer match the local composition.

Composition-based statistics reweights the substitution matrix on a
per-alignment basis to use the empirical composition of the two
sequences in place of the universal background frequencies. The
effect is to bring the bit score back down to the value it would have
had under matched composition. The correction is on by default in
modern BLAST and is the reason the output line reads "`Method:
Composition-based stats`" — it is the algorithm telling you the bit
score you see has already been adjusted.

#warn[
  Composition correction does not save you from low-complexity hits
  that survive the filter. The standard practice is to mask
  low-complexity regions of the query with the SEG filter (Wootton &
  Federhen, 1993) before search; modern BLAST does this by default.
  If you turn the filter off — there are occasional reasons to —
  expect to wade through long lists of high-scoring nonsense against
  signal peptides and proline-rich linkers. The bit score is real;
  what the bit score is _evidence of_ is not.
]


== Karlin–Altschul Statistics and E-values <sec:karlin-altschul>

The first half of BLAST is the algorithm. The second half — the part
that turns BLAST from a heuristic into a measurement instrument — is
the statistics. The single most important question a search must
answer is: given that I just scored an HSP at bit score $S'$, how
surprised should I be? The answer comes from a result Karlin and
Altschul derived in 1990 and refined in 1993, predating the BLAST
software itself.

=== The Null Model

Imagine running BLAST against a *random database* — a database of
the same total size $n$ as the real one, but with residues drawn
independently from the background distribution $f_a$. What is the
distribution of HSP bit scores you would see purely by chance?

If the answer were "Gaussian with some mean and variance," BLAST
would be in serious trouble. The Gaussian tail decays slowly enough
that, at the bit scores actual homologs produce (50 to 200 bits),
the random-database tail would still contain many false positives.
But the Karlin–Altschul theorem shows that the answer is much
better-behaved. The maximum local-alignment score between two random
sequences of length $m$ and $n$ follows an *extreme-value (Gumbel)
distribution*, whose density at the right tail decays as
$exp(-lambda S)$ — exponentially in the raw score, with a matrix-
and-gap-specific rate parameter $lambda$. The expected number of
HSPs at bit score $S'$ or greater is

$ E = K dot m dot n dot e^(-lambda S) $

with $K$ another matrix-specific constant. Converting to bit scores
via $S' = (lambda S - log K) / log 2$, this collapses to the form
worth memorising:

$ E = m dot n dot 2^(-S') $

The expected-by-chance count of HSPs at bit score $S'$ falls
exponentially in $S'$ and grows linearly in the product of query
length and database size. That single equation is the entire
Karlin–Altschul payload, and it is what makes BLAST a measurement
instrument rather than a heuristic.

#note[
  The Gumbel distribution arises here for the same reason it arises in
  any extreme-value problem: the maximum of many independent
  identically distributed random variables, suitably normalised,
  converges in distribution to a Gumbel. In BLAST's setting, the
  random variables are HSP scores at many independent diagonals of
  the alignment matrix, and the maximum is the score of the best HSP.
  Karlin and Altschul's contribution was to extend the classical
  extreme-value theory to handle dependent sequence positions and
  arbitrary substitution matrices.
]

=== Reading the E-value

Let the database be `nr` ($n approx 10^{11}$ residues), the query
length $m = 100$ residues. Then

$ E approx 10^(13) dot 2^(-S') $

so a bit score of 30 gives $E approx 10^{13} dot 10^(-9) approx 10^4$ —
ten thousand HSPs of this strength expected by chance, plainly not
significant. A bit score of 50 gives $E approx 10^{13} dot 10^(-15)
approx 10^(-2)$ — borderline, worth a second look. A bit score of 80
gives $E approx 10^(-11)$ — overwhelming. Each extra bit halves the
expected by-chance count. The exponential dependence is what makes a
ten-bit difference matter so much: between $S' = 40$ and $S' = 50$
the E-value changes by a factor of a thousand.

The classic significance thresholds people quote:

- $E < 10^(-3)$: probably significant; worth investigating.
- $E < 10^(-10)$: almost certainly homologous; safe to act on.
- $E < 10^(-50)$: bulletproof; transitive function transfer is justified.

These are not statistical theorems — they are community-empirical
conventions calibrated against the fact that a typical search returns
a manageable number of false positives at $E < 10^(-3)$ on a typical
database.

#figure(
  image("../../diagrams/lecture-19/03-evalue-curve.svg", width: 95%),
  caption: [
    The Karlin–Altschul map from bit score to E-value. The left panel
    shows the Gumbel tail of random-database HSP scores; the right
    panel plots E-value against bit score for three database sizes
    spanning six orders of magnitude. Doubling the database doubles
    the E-value at any fixed bit score; ten additional bits drop the
    E-value by a factor of a thousand.
  ],
) <fig:evalue>

#figure(
  image("../figures/ch19/f1-detection-cascade.svg", width: 95%),
  caption: [
    BLAST as a two-stage Neyman–Pearson detector with built-in
    multiple-testing correction. The seed scan plays the role of the
    coarse matched filter; the ungapped extension plays the role of
    the fine matched filter; the bit-score threshold is the
    likelihood-ratio decision boundary; the E-value is the expected
    false-positive count at that decision boundary against a database
    of size $N$.
  ],
) <fig:detection>

=== P-values, E-values, and the Distinction

The relationship is $P = 1 - exp(-E)$. For the small-$E$ regime that
actually matters, $P approx E$ and you can read either off the BLAST
output line interchangeably. They start to differ when $E$ grows
toward and past one: $P$ is bounded above by 1 because it is a
probability, while $E$ can exceed 1 because it is an _expected count_.
A BLAST hit reporting $E = 5$ means "you expect to see this many
HSPs at this bit score by chance against a database of this size"
— which is more informative than the equivalent $P = 0.993$, because
the count form makes it obvious that the hit is unimpressive without
having to reason about how close 0.993 is to 1. BLAST reports
E-values for exactly this reason.

=== Database-Size Effects, and Why They Matter

The same alignment with the same bit score produces _different_
E-values depending on which database you searched against. Search
your query against the SwissProt curated subset ($n approx 2 dot 10^8$
residues) and the E-value reflects multiple-testing against a few
hundred thousand records. Search the same query against `nr` ($n
approx 10^(11)$ residues) and the E-value is roughly five hundred
times larger for the same alignment. The bit score is invariant; the
E-value is not.

#figure(
  image("../../diagrams/lecture-19/10-database-scaling.svg", width: 95%),
  caption: [
    E-value scaling with database size at fixed bit score $S' = 50$.
    The E-value grows linearly with database size; the bit score does
    not change. Same alignment, different statistical interpretation
    depending on the haystack you searched it against.
  ],
) <fig:db-scaling>

The practical consequence is that the "right" E-value depends on
what you are trying to claim. If you are asking "does this alignment
exist as the best match anywhere in the universe of known proteins?"
then your database is the whole biosphere and your E-value is the
`nr`-sized one. If you are asking "does this alignment exist among
the manually-curated reference proteome of a single species?" then
your database is SwissProt-sized and your E-value reflects that
smaller multiple-testing burden. Some workflows run the search
against a small curated database to get a manageable hit list, then
re-test the surviving hits against the full database to get an
honest measure of significance.

#tip[
  When reading someone else's BLAST results, the first three things
  to check are which database they searched, which matrix they used,
  and whether composition correction was on. The bit score and E-value
  are meaningful only relative to those three choices. A reported
  $E = 10^(-30)$ against `swissprot` is not the same evidence as an
  $E = 10^(-30)$ against `nr`. The headline number does not stand
  alone.
]

=== Anatomy of a BLAST Hit

The BLAST output is dense and worth learning to read. A typical hit
record looks like this:

```
> tr|A0A024R7T9|A0A024R7T9_HUMAN  Some Protein
Length=247
Score = 187 bits (475),  Expect = 4e-49, Method: Composition-based stats.
Identities = 95/120 (79%),  Positives = 102/120 (85%),  Gaps = 0/120 (0%)
```

The fields, in the order they appear:

- *Score*: bit score (187), with the raw score in parentheses (475).
  The bit score is matrix-independent and is what you read for
  cross-comparison; the raw score is matrix-dependent and is reported
  for historical reasons.
- *Expect*: the E-value (4 × 10⁻⁴⁹). The headline statistical claim.
- *Method*: which score-correction scheme was used. "Composition-based
  stats" is the modern default.
- *Identities*: the number and proportion of positions where the two
  sequences agree exactly.
- *Positives*: the number and proportion of positions where the
  substitution-matrix score is positive — that is, conservative
  substitutions plus identities. Always at least as large as
  Identities.
- *Gaps*: the number and proportion of alignment positions where
  one sequence is gapped.

#figure(
  image("../../diagrams/lecture-19/12-blast-output.svg", width: 95%),
  caption: [
    Anatomy of a BLAST hit record. Each numeric field carries a
    specific meaning; the bit score is the matrix-independent
    headline, the E-value is the calibrated false-discovery
    measurement, and the identities-versus-positives split reveals
    how much of the alignment is conservative substitution.
  ],
) <fig:output>

#figure(
  image("../figures/ch19/f2-bit-score-to-evalue.svg", width: 95%),
  caption: [
    A reference card for the bit-score-to-E-value conversion at three
    canonical database sizes. The same alignment can be safely
    significant against a curated database and merely suggestive
    against `nr`; the bit score is the invariant the table is
    indexed on.
  ],
) <fig:lookup>

A useful diagnostic discipline: when the Identities figure is low
(say, 25 %) but the Positives figure is high (say, 60 %), the
alignment is finding a homolog whose surface residues have drifted
extensively but whose structural and functional residues remain
conservatively substituted. That is what distant homology looks
like at the protein level. When Identities and Positives are both
high (say, 90 %+) the alignment is finding a close #idx("paralog")paralog or
ortholog. When Identities are very low and Positives are barely
above them, the alignment is borderline and worth inspecting by
eye before believing.


== PSI-BLAST and Profile-Based Search <sec:psi-blast>

Single-pass BLASTP detects homologs reliably down to about 25 %
identity. Below that — in the so-called _twilight zone_, named by
Rost (1999) for the region where pairwise identity is no longer a
reliable predictor of homology — the sensitivity falls off. But
biologically important relationships often live in the twilight zone:
distantly related kinases share only 15-20 % identity across their
catalytic domains, #idx("transcription factor")#idx("transcription")transcription factor families diverge to single-
digit pairwise identities while preserving the same DNA-binding
fold. A single-query, single-matrix BLAST will miss many of them.

*PSI-BLAST* — Position-Specific Iterated BLAST, Altschul et al.,
_NAR_, 1997 — is the iterated profile-refinement procedure that
extends sensitivity into the twilight zone. The construction is
straightforward in outline:

1. *Iteration 0*: run standard BLASTP with the query against the
   database. Collect all hits with $E < E_("include")$ (default
   $0.005$).
2. *Profile construction*: build a #idx("multiple sequence alignment")multiple sequence alignment from
   the surviving hits and from it construct a *Position-Specific
   Scoring Matrix (PSSM)* — a $20 times L$ matrix giving a per-position
   substitution score, where $L$ is the query length.
3. *Iteration $gt.eq 1$*: run BLAST again, but use the PSSM in place
   of BLOSUM62 for scoring.
4. Iterate until convergence (no new hits cross the inclusion
   threshold) or until a maximum iteration count (default 5).

Each iteration refines the PSSM by incorporating the patterns of
substitution seen across the hit set. Conserved positions accumulate
sharp scoring profiles — high score for the conserved residue,
strongly negative for everything else. Variable positions develop
flat profiles that score near zero. After a few iterations, the
PSSM is no longer the query's blunt BLOSUM62 view of the world but
a family-aware scoring scheme that lets distant homologs score above
the inclusion threshold and feed back into the next iteration's
profile.

#figure(
  image("../../diagrams/lecture-19/04-psi-blast.svg", width: 92%),
  caption: [
    PSI-BLAST as iterated profile refinement. Iteration 0 runs
    standard BLAST; subsequent iterations build a PSSM from the
    current hit set and use it as the scoring matrix for the next
    search. Each iteration recruits more distant homologs as the
    profile sharpens at conserved positions. The procedure caps at
    5 iterations to limit drift.
  ],
) <fig:psi-blast>

=== Why PSI-BLAST Works

The procedure is *expectation–maximisation in disguise*. The E-step is
the BLAST search: given the current PSSM, identify which database
sequences are in the homology set. The M-step is the PSSM rebuild:
given the homology set, recompute the best position-specific scoring
matrix. Iterating the two converges (when it converges) on a fixed
point at which the PSSM and the hit set are mutually consistent — the
PSSM scores the hit set members above threshold, and no new database
sequence falls above threshold under the PSSM.

The PSSM construction itself is straightforward. For each column $i$
of the multiple alignment, compute the observed frequency $q_(i, a)$
of amino acid $a$, and define

$ M_(i, a) = log frac(q_(i, a), f_a) $

where $f_a$ is the background frequency. Dirichlet pseudo-counts are
added to avoid zero-frequency artefacts at columns where some amino
acid never appears in the current hit set — a column where every
sequence happens to show D should still give Glu (E) a small but
non-zero positive score, because D ↔ E substitutions are conservative
even if the current sample did not exhibit one.

=== Where PSI-BLAST Fails

The EM-style structure brings the EM-style failure mode:
*model drift*. If a single non-homologous sequence — a chance hit
that crossed the inclusion threshold by random alignment — gets
added to the profile in iteration 1, its residue frequencies bias
the PSSM in iteration 2. In iteration 3, the drifted PSSM admits
more sequences that match the contamination's pattern rather than
the original query's. By iteration 5 the profile may have wandered
into a different protein family entirely.

The practical mitigations are operational rather than algorithmic.
Tighten the inclusion threshold to $E < 10^(-4)$ or so for
distantly related searches; inspect the alignment between iterations
and manually exclude obvious non-homologs from the profile; cap the
iteration count at 3 rather than 5; treat the final hit list with
appropriate skepticism if any of the early-iteration profiles looked
mixed.

#warn[
  PSI-BLAST drift is silent. The E-values that drifted iterations
  produce are computed against the (drifted) profile, so the
  per-hit statistics still look fine from inside the procedure. The
  only way to catch drift is to look at the multiple alignment that
  produces the PSSM and verify that the sequences in it are
  biologically sensible. Trust the procedure's headline only as far
  as you have audited the profile.
]

=== Connection to #idx("HMMER")HMMER

The next algorithmic step up from PSI-BLAST is the *#idx("profile HMM")profile #idx("HMM")HMM*
implemented in HMMER (Eddy, 1995 and ongoing) and used to power
#idx("Pfam")Pfam's domain classifications. A profile HMM is what you get when
you replace PSI-BLAST's position-specific scoring matrix with an
explicit *#idx("hidden Markov model")hidden Markov model* whose states are aligned columns,
with separate insertion and deletion states between them, and whose
emission and transition probabilities are fit by #idx("maximum likelihood")maximum likelihood
on the multiple alignment. The expressive power is strictly greater
than PSI-BLAST's PSSM — HMMs model insertion and deletion as
position-specific events with their own transition probabilities,
which a flat PSSM cannot do.

PSI-BLAST is to HMMER what $k$-means is to a Gaussian mixture model:
the same modelling intent, simpler machinery, faster on small
problems, less expressive on hard ones. The two are complementary in
practice. PSI-BLAST is the fast iterative refinement for moderate
homology; HMMER (and the closely related HHblits) is the
heavier-weight tool for deep homology where you have already invested
in building or downloading a profile HMM for the family. Chapter 20
returns to multiple alignment and profile models as the basis for
phylogenetic and structural inference.


== Modern Fast Alternatives <sec:modern>

Classical BLAST was designed for the hardware of the early 1990s —
MHz-class CPUs, single-digit megabytes of RAM, gigabyte-class disks.
The algorithm scales linearly in the database size $n$, but with a
constant factor that, on hardware from 2025 running queries against
a 2025-sized database, makes a `nr`-scale BLASTP search take minutes
per query at best. The protein-sequence universe has grown about
five orders of magnitude since 1990. The hardware has grown about
five orders of magnitude. The algorithm has not.

The modern successor tools — DIAMOND, MMseqs2, USEARCH — preserve the
seed-extend-score architecture but exploit one or more modern
hardware affordances: large RAM (allowing the entire database index
to live in memory), SIMD vectorisation (allowing the
substitution-matrix inner loop to process 16 or 32 amino acids per
instruction), multi-core parallelism (allowing the seed-scan to scale
linearly with thread count), and cache-aware data layout (allowing
the database to stream through L2 cache without thrashing). The
result is a 100- to 1000-fold speed-up at sensitivity that is
typically within a few percent of classical BLAST on biologically
relevant tasks.

=== DIAMOND

*DIAMOND* (Buchfink, Xie, & Huson, _Nature Methods_, 2015) is the
de-facto BLASTP replacement at metagenomic scale. The key engineering
moves are *double-indexed seeds* — the algorithm indexes both the
query (or query batch) and the database in advance, then matches
index entries against each other rather than scanning one against the
other — and *spaced seeds*, which allow certain non-contiguous
patterns of matching positions in the seed (e.g., a pattern of `11011`
that requires identity at positions 1, 2, 4, 5 and tolerates anything
at position 3). Spaced seeds were introduced by PatternHunter (Ma,
Tromp, & Li, 2002); they raise sensitivity per seed by allowing more
informative seeds without raising the false-positive rate. DIAMOND
also uses *block-based cache-friendly batch processing* so that
chunks of the database stream through the CPU caches without random
access.

The net effect is a 100- to 1000-fold speedup over BLASTP at roughly
95 % sensitivity on metagenomic-scale tasks. DIAMOND is the default
search engine in essentially every metagenomic annotation pipeline
written after about 2016. For high-throughput protein-to-database
search where the user can tolerate missing a few of the most
distant homologs, it is the right default.

=== MMseqs2

*MMseqs2* (Steinegger & Söding, _Nature Biotechnology_, 2017) takes
a different approach to the same target: deeper sensitivity (down to
the ~10 % identity range) at speeds competitive with DIAMOND. The
core moves are a *cascaded $k$-mer prefilter* — a series of
increasingly stringent $k$-mer-based filters that progressively
narrow the candidate set — followed by SIMD-vectorised gapless and
gapped extension. MMseqs2 in cluster mode handles UniRef-scale
clustering in hours where pre-MMseqs tools required weeks; it is the
engine behind several major reference cluster databases (UniRef30,
ColabFoldDB).

=== USEARCH and UCLUST

*USEARCH* (Edgar, _Bioinformatics_, 2010) and its open-source variant
UCLUST are the speed champions at the high-identity end. USEARCH
sorts database sequences by length, precomputes $k$-mer profiles per
sequence, and uses these as cheap pre-filters. The tool is
extraordinarily fast (orders of magnitude beyond BLAST) at clustering
and search where the targets share 70 % or more identity — exactly
the regime that dominates #idx("microbiome")microbiome #idx("16S rRNA")16S rRNA workflows, amplicon
clustering, and chimera detection. It is less effective at distant
homology and is rarely the right choice when the question is "find
the closest UniProt protein to my newly discovered effector."

=== The Speed–Sensitivity Pareto

@fig:pareto sketches where each tool sits in the speed–sensitivity
plane. The frontier is real: at fixed sensitivity, MMseqs2 is the
fastest; at fixed throughput, MMseqs2 and DIAMOND are nearly
indistinguishable until the deep-homology regime where MMseqs2 pulls
ahead. USEARCH dominates the high-identity corner; BLASTP holds the
high-sensitivity corner at the cost of throughput. The right choice
depends on the workload, not on which tool has the highest brand
recognition.

#figure(
  image("../../diagrams/lecture-19/05-pareto.svg", width: 95%),
  caption: [
    The speed–sensitivity Pareto for protein database search. MMseqs2
    sits on the frontier across most distance regimes; DIAMOND in its
    sensitive mode is essentially indistinguishable on close to
    moderate homology; classical BLASTP holds the deep-sensitivity
    corner at the cost of throughput; USEARCH dominates the
    high-identity, high-throughput corner that 16S microbiome
    pipelines occupy.
  ],
) <fig:pareto>

#figure(
  image("../figures/ch19/f4-tool-decision-tree.svg", width: 95%),
  caption: [
    A decision tree for choosing among BLASTP, DIAMOND, MMseqs2, and
    USEARCH. The question to answer first is what identity range you
    expect; the second is whether throughput or sensitivity dominates
    your constraint; the third is whether the workflow has a curated
    standard that you should not deviate from for reproducibility.
  ],
) <fig:decision>

#note[
  BLAST is still the canonical workhorse despite being slower because
  three factors keep it there. Its E-value calibration is the
  community reference standard — every benchmark in the field is
  ultimately compared to BLAST's significance numbers, and the
  modern tools justify themselves relative to BLAST's sensitivity.
  Its output format is the lingua franca that downstream tools
  expect to parse. And for any single-query, manual search at the
  NCBI web interface, the latency is acceptable and the workflow is
  the one every collaborator already knows. The high-throughput
  alternatives win on pipeline-scale work; the web interface remains
  BLAST.
]


== Application Patterns <sec:applications>

The previous sections taught you the algorithm and the statistics.
This one walks the handful of recurring workflows that consume
BLAST output. Each is built from the same primitives — search,
filter by E-value, transfer some property from the hit to the query —
but the choice of database, the choice of E-value threshold, and the
follow-up interpretation differ by task.

=== Genome Annotation

When a new genome is assembled (Chapter 3), the first-pass functional
annotation runs every predicted protein-coding gene through a
database search and labels each gene with the closest hit's name,
function, and domain assignments. The database is typically UniProt
(or a clade-appropriate subset) plus Pfam (Chapter 20) for domain
hits. The search engine is BLAST or, increasingly, DIAMOND for
throughput.

NCBI's PGAP (Prokaryotic Genome Annotation Pipeline), Ensembl's
gene-build pipeline, MAKER for eukaryotic annotation, and #idx("BRAKER")BRAKER for
RNA-seq-assisted prokaryotic annotation all use BLAST or DIAMOND
under the hood as the protein-search step. The pattern is universal
enough that the field treats "annotated by sequence similarity" as a
defining feature of any modern genome project.

=== Reciprocal Best Hits for Ortholog Identification

*Orthologs* — genes in different species that descend from a single
ancestral gene by speciation, with no intervening duplication — are
the comparative-genomics primitive most often asked of a BLAST
output. The simplest and most widely used method is *reciprocal best
hits (RBH)*:

1. BLAST gene $A$ from species 1 against species 2. Call the best
   hit $B$.
2. BLAST gene $B$ from species 2 against species 1. If the best hit
   is $A$, then $(A, B)$ is a reciprocal best-hit pair.
3. RBH pairs are predicted orthologs, with caveats around recent
   gene duplications.

#figure(
  image("../../diagrams/lecture-19/11-rbh.svg", width: 92%),
  caption: [
    Reciprocal best hits for ortholog identification. Bidirectional
    best hits across species predict orthology with about 80 %
    accuracy; the failure mode is recent paralog duplication, where
    the "best hit" in one direction is the wrong paralog. Modern
    extensions (#idx("OrthoFinder")OrthoFinder, OrthoMCL) handle paralog confounds by
    clustering.
  ],
) <fig:rbh>

RBH recovers about 80 % of true orthologs on benchmark datasets and
fails predictably on recent gene duplications: if species 1 has a
recently duplicated $A_1, A_2$ and species 2 has the unduplicated
ancestor $B$, then both $A_1$ and $A_2$ may best-hit $B$, but $B$'s
best hit is only one of them — the other is co-orthologous to $B$
through the duplication, and the RBH procedure misses the
relationship. *OrthoMCL* (Li, Stoeckert, & Roos, 2003) and
*OrthoFinder* (Emms & Kelly, 2015 and 2019) extend RBH by clustering
all-versus-all BLAST hits with Markov clustering or by inferring
duplication-aware ortho-groups; both are the modern replacements when
the RBH approximation is not good enough.

=== Protein Function Transfer

Once a BLAST hit clears a high E-value threshold against a curated
database, the dominant downstream move is *transitive function
transfer*: take the curated functional annotation of the hit and
attach it to the query. The workflow:

1. Query protein $X$ against UniProt with `blastp` or `diamond`.
2. The top hit is some $Y$ with $E approx 10^(-50)$ and 60 % identity
   across 80 % of the query length.
3. $Y$ has a curated #idx("Gene Ontology")Gene Ontology annotation, a Pfam domain call,
   an enzyme-commission number, perhaps a structure in the PDB.
4. Transfer those annotations to $X$, with confidence proportional
   to the identity, the alignment #idx("coverage")coverage, and the E-value.

The empirical reliability of function transfer is sharply
identity-dependent. Above 50 % identity over most of the protein
length, function is preserved across about 90 % of hits — the
biochemistry, the substrate specificity, the cellular role are
usually the same. Below 30 % identity, function transfer is risky:
the proteins likely share fold and possibly active-site
architecture, but specific function (which substrate, which
inhibitor, which protein partner) often diverges. The 25-30 %
identity band is the twilight zone where transfer is permissible only
with explicit caveats and ideally with corroborating evidence from
structure or experiment.

#figure(
  image("../figures/ch19/f3-function-transfer-pyramid.svg", width: 95%),
  caption: [
    Function-transfer confidence tiers by sequence identity. Above
    50 % identity and 70 % coverage, transitive annotation is reliable
    for specific function; the 30-50 % band carries fold-level but not
    specific-function confidence; below 30 % the transfer is
    permissible only with corroborating evidence from structure,
    experiment, or profile-based search.
  ],
) <fig:transfer>

#warn[
  Function transfer is a chain of inferences-from-inferences. The
  hit's annotation may itself have been transferred from another
  hit, which in turn was transferred from an experimental
  determination several jumps back. Each hop accumulates uncertainty
  the per-hit E-value does not capture. UniProt's evidence codes
  (`EXP` for experimentally determined, `ISS` for inferred by
  sequence similarity) flag the provenance; treat `ISS`-annotated
  hits below 40 % identity with the skepticism the chain warrants.
]

=== Horizontal Gene Transfer Detection

Most genes in a bacterial genome are vertically inherited and produce
BLAST hits whose top matches are in closely related species. *Genes
acquired by horizontal transfer* — phage, plasmid, or environmental
#idx("DNA")DNA picked up by a microbe and integrated into its genome — show an
*anomalous BLAST signature*: a gene whose closest hits lie in a
distant phylum rather than in sibling species. The detection
discipline is simple in principle: for each gene, compare the
taxonomic distribution of high-scoring BLAST hits against the
expected distribution under vertical inheritance, and flag the
anomalies.

The pattern is the workhorse for tracking the spread of antibiotic-
resistance genes across pathogen species, identifying eukaryote-to-
bacterium transfers in symbionts, and reconstructing phage-host
integration histories. The false-positive rate is non-trivial — gene
loss in some lineages and rapid evolution in others can mimic the
HGT signature — so the convention is to combine BLAST-based
anomaly detection with phylogenetic analysis of the candidate gene
and #idx("synteny")synteny analysis of the surrounding genomic context.

=== Drug Target and Off-Target Prediction

When a drug binds to protein $A$ with known affinity, BLAST against
the human proteome reveals which other proteins are closely enough
related to $A$ that they may also bind the drug — the *off-target
candidates*. The kinase family is the canonical example: human
kinases share roughly 30-50 % identity in their catalytic domains,
so even highly selective kinase inhibitors hit five to ten off-target
kinases at clinically relevant concentrations. Pharmaceutical
selectivity profiling routinely starts with a BLAST off-target list,
then narrows it by experimental screening.

The same pattern appears in antibiotic discovery (target a bacterial
enzyme, BLAST to confirm no human homolog), in herbicide design
(target a plant enzyme, BLAST to confirm no animal homolog), and in
toxicology (a candidate drug-toxic-effect protein, BLAST to find
related proteins in non-target organisms). BLAST does not
predict whether a specific small molecule binds; it predicts which
proteins are similar enough that the same binding mode is
geometrically plausible. The full off-target prediction also requires
#idx("docking")docking, structural alignment, and ideally experimental binding
assays — but the BLAST off-target list is the universal first pass.


== Summary <sec:ch19-summary>

- *BLAST is heuristic search with statistical guarantees.* The
  seed-and-extend cascade gives up exact Smith–Waterman optimality
  and in exchange recovers a calibrated E-value: the expected number
  of hits this strong by chance against a random database of this
  size.

- *The pipeline has five stages.* Build a word list with substitution-
  matrix neighbourhoods, scan the database by hash for seeds, gate
  on two close seeds, ungap-extend with X-drop, gap-extend with
  banded Smith–Waterman. Each stage is cheap enough to apply to the
  survivors of the previous one.

- *Bit scores are log-likelihood ratios.* The conversion
  $S' = (lambda S - log K) / log 2$ rescales any raw score to a
  matrix-independent additive unit. A bit score of $S'$ means the
  alignment is $2^(S')$ times more likely under homology than under
  random — additive in evidence, comparable across runs.

- *Karlin–Altschul gives the E-value.* Random-database HSP scores
  follow an extreme-value distribution; the expected count above bit
  score $S'$ is $E = m n 2^(-S')$. Linear in query length and
  database size, exponential in bit score. Doubling the database
  doubles the E-value at fixed bit score; ten extra bits drop the
  E-value by a factor of a thousand.

- *PSI-BLAST is EM-style profile refinement.* Build a PSSM from the
  current hit set, re-search with the PSSM, iterate to convergence.
  Sensitivity extends into the twilight zone (~10-20 % identity) at
  the cost of drift risk; cap at three iterations in practice and
  audit the profile between rounds.

- *Modern alternatives buy 100-1000× speed.* DIAMOND and MMseqs2
  preserve the BLAST architecture and exploit large RAM, SIMD,
  spaced seeds, and cache-friendly indexing. They lose a few percent
  sensitivity on the deepest homologs and are the right default for
  any pipeline-scale work.

- *Applications follow from the primitives.* Genome annotation is
  bulk BLAST against UniProt. Ortholog identification is reciprocal
  best hits. Function transfer is annotation propagation along
  high-identity hits. HGT detection is anomalous-taxonomy hit lists.
  Off-target prediction is BLAST against the host proteome. Every
  workflow is a particular E-value threshold and a particular
  downstream interpretation of the hit.

- *Read the bit score and the E-value together.* The bit score is
  the invariant evidence; the E-value is the multiple-testing-
  corrected interpretation of that evidence against a particular
  database. Both numbers belong on the page; neither alone is
  sufficient.


== Exercises <sec:ch19-exercises>

#strong[1.] #emph[Bit score to E-value, by hand.] A query of length
100 amino acids is searched against a database of $10^9$ amino-acid
residues. An HSP comes back with a bit score of 45. Compute the
expected E-value. Compute the bit score that would give $E = 10^(-3)$
on the same query and database. Now hold the bit score fixed at 45
and recompute the E-value for a database of $10^{11}$ residues — how
much does the significance change?

#strong[2.] #emph[Read a BLAST hit.] Given the hit record below,
identify (a) the bit score, (b) the raw score, (c) the E-value, (d)
the percentage of identical positions, (e) the percentage of
positively scoring positions, and (f) whether composition correction
was applied. Comment on whether this hit looks like a close paralog,
a distant homolog, or something borderline.

```
> sp|P12345|EXAMPLE_HUMAN  Example protein
Length=312
Score = 142 bits (358),  Expect = 2e-37, Method: Composition-based stats.
Identities = 84/240 (35%),  Positives = 138/240 (57%),  Gaps = 8/240 (3%)
```

#strong[3.] #emph[Substitution matrix arithmetic.] Look up the
BLOSUM62 entries for the pairs (D, D), (D, E), (D, K), and (W, A).
Convert each to a log-odds interpretation: how much more (or less)
likely are these pairs to co-occur in homologous positions than under
the random background? Repeat for PAM250 and comment on which pairs
differ most between the two matrices.

#strong[4.] #emph[Implement the X-drop rule.] Write a Python or
MATLAB function that takes two equal-length sequences and a
substitution matrix and performs ungapped extension from a given
seed position in both directions under an X-drop rule with
configurable $X$. Test it on a synthetic homologous region with
embedded noise. Plot the running score and mark where X-drop
triggers. Vary $X$ from 5 to 50 and observe the effect on the HSP
boundary.

#strong[5.] #emph[Reciprocal best hits.] Take the predicted protein
sets of two closely related bacterial species from NCBI (e.g.,
two _Escherichia coli_ strains, or _E. coli_ K-12 and _Salmonella
enterica_ Typhimurium). Run all-versus-all BLASTP. Identify
reciprocal best hits at $E < 10^(-10)$. Report the size of the RBH
set as a fraction of each genome's gene count. Sample five RBH
pairs and inspect them for evidence of paralog confounds (the
"best hit" in each direction is the wrong paralog).

#strong[6.] #emph[PSI-BLAST drift.] Pick a kinase catalytic domain
(about 250 amino acids) from a curated source such as Pfam's
`Pkinase` family. Run PSI-BLAST against `nr` for five iterations
with the default inclusion threshold. Track (a) the number of hits
per iteration, (b) whether new sequences entering the hit list at
each iteration are kinases by their UniProt annotation, and (c)
whether you can detect drift — non-kinase domain hits creeping in
late. Repeat with a tightened inclusion threshold of $E < 10^(-4)$
and comment on the difference.

#strong[7.] #emph[Database-size effect.] Choose a single query
protein. Run BLASTP against SwissProt, then against UniRef90, then
against `nr`. For the same top hit (it should be the same hit in
each search), report the bit score and the E-value from each
database. Verify by hand that the E-values scale linearly with the
quoted database sizes.

#strong[8.] #emph[(Open-ended.)] Pick a published tool or method
that uses BLAST in a non-obvious way — examples might be HGT
detection in a specific clade, phylogenetic placement of
metagenomic reads, or an antibiotic-resistance surveillance system.
Describe in one paragraph how the tool uses the BLAST output: what
database it searches against, what E-value threshold it applies,
what it does with the hit-versus-no-hit decision, and how it
handles the inevitable false-positive rate at scale.


== Further Reading <sec:ch19-further-reading>

- *Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman,
  D. J.* (1990). "Basic local alignment search tool." _Journal of
  Molecular Biology_ 215: 403–410. The original BLAST paper.
  Surprisingly readable; the seed-and-extend logic and the
  Karlin–Altschul application are both laid out clearly.

- *Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J.,
  Zhang, Z., Miller, W., & Lipman, D. J.* (1997). "#idx("gapped BLAST")Gapped BLAST and
  PSI-BLAST: a new generation of protein database search programs."
  _Nucleic Acids Research_ 25: 3389–3402. The BLAST 2.0 paper,
  introducing gapped extension, two-hit seeding, and the PSI-BLAST
  iterative profile algorithm.

- *Karlin, S., & Altschul, S. F.* (1990). "Methods for assessing the
  statistical significance of molecular sequence features by using
  general scoring schemes." _PNAS_ 87: 2264–2268. The Gumbel
  distribution result; the entire statistical foundation of BLAST is
  derived in these few pages.

- *Henikoff, S., & Henikoff, J. G.* (1992). "Amino acid substitution
  matrices from protein blocks." _PNAS_ 89: 10915–10919. The BLOSUM
  paper. Sets out the BLOCKS-derived methodology and the resulting
  matrices that displaced PAM as the community default.

- *Buchfink, B., Xie, C., & Huson, D. H.* (2015). "Fast and sensitive
  protein alignment using DIAMOND." _Nature Methods_ 12: 59–60.
  The DIAMOND paper. Short; the technical detail lives in the
  supplement.

- *Steinegger, M., & Söding, J.* (2017). "MMseqs2 enables sensitive
  protein sequence searching for the analysis of massive data sets."
  _Nature Biotechnology_ 35: 1026–1028. The MMseqs2 paper.
  The clearest current statement of how to push deep-homology
  sensitivity at metagenomic scale.

- *Rost, B.* (1999). "Twilight zone of protein sequence alignments."
  _Protein Engineering_ 12: 85–94. The empirical mapping between
  pairwise identity and reliable homology inference, and the
  source of the "twilight zone" terminology.

- *Pearson, W. R.* (2013). "An introduction to sequence similarity
  ('homology') searching." _Current Protocols in Bioinformatics_
  Chapter 3, Unit 3.1. A working reference for the entire space; if
  you read one chapter-length overview to back up this one, this is it.
