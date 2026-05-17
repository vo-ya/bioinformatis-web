#import "../theme/book-theme.typ": *

= HMMs, Profile HMMs, and Gene Finding <ch:hmms>

#matters[
  Hidden Markov Models are the workhorse of segmentation in genomics.
  Roughly every protein domain assignment in UniProt, every gene in a
  newly assembled bacterial genome, every chromatin-state call in the
  ENCODE atlas, and every indel-aware variant likelihood inside #idx("GATK")GATK's
  #idx("HaplotypeCaller")HaplotypeCaller is — at the bottom of the stack — a #idx("Viterbi")Viterbi or a
  forward pass on an #idx("HMM")HMM. The architecture was designed in the late
  1960s for speech recognition, was adapted to biology in the early
  1990s, and refuses to retire. Deep neural networks have displaced
  HMMs in some tasks and not in others, and the distinction is
  predictable enough to be a design principle. An EE student who has
  seen Kalman filtering already knows half of this chapter; the other
  half is the careful adaptation of the same machinery to discrete
  symbol streams over an alphabet of four bases or twenty amino acids.
]

A genome, viewed as an information channel, is a discrete observation
stream punctuated by long stretches that share a hidden regime. Inside
a gene the #idx("codon")codon statistics look one way; inside an #idx("intron")intron, another;
inside an intergenic spacer, a third. A #idx("CpG island")#idx("CpG")CpG island has a different
dinucleotide composition than the bulk genome that surrounds it. A
#idx("promoter")promoter is statistically distinguishable from the open reading frame
it precedes. None of these labels — gene, intron, island, promoter —
appears in the FASTA file. They are latent regimes that the sequence
emits clues about. The task of recovering them, position by position,
is the segmentation problem that defines half of all annotation work
in genomics.

The #idx("hidden Markov model")Hidden Markov Model is the cleanest probabilistic statement of that
task. It posits a small unobservable state that controls the local
emission distribution and a Markov transition between states. Given a
genome, three things you might want come for free in $O(T K^2)$ time:
the likelihood of the observation under the model, the most probable
state sequence consistent with it, and — given only observations and
no labels — the maximum-likelihood parameters that generated it. Those
three computations are the forward, the Viterbi, and the #idx("Baum-Welch")Baum-Welch
algorithms. Everything else in this chapter is a particular choice of
state space.

The historical arc is worth keeping in mind as you read. Andrey
Markov defined the chains that bear his name in 1906, in an analysis
of vowel and consonant transitions in Pushkin's _Eugene Onegin_.
Leonard Baum and Ted Petrie at the Institute for Defense Analyses
extended the chain to a hidden-state variant in 1966 and gave the
EM-style training algorithm its first complete derivation in 1972.
Andrew Viterbi's max-product decoder appeared in 1967, in a paper on
convolutional codes; the same algorithm was rediscovered for speech
recognition a decade later. Lawrence Rabiner's 1989 tutorial in the
_Proceedings of the IEEE_ canonicalised the notation. The transfer to
biology — profile HMMs for protein families (Krogh, Mian, Sjölander,
and Haussler in 1994) and gene finding (Burge and Karlin's #idx("GENSCAN")GENSCAN in
1997) — followed within a decade. The #idx("HMMER")HMMER software package (Sean
Eddy, 1998 and 2011) packaged the protein-family inference into a
tool that has been the default first pass in genome annotation for
twenty-five years.

The chapter is laid out in the order most readers will need it. The
first section sets up the model and the three canonical problems.
Sections 2 and 3 develop the forward, backward, Viterbi, and Baum-Welch
recurrences with the numerical-stability and initialisation tricks that
working implementations require. Section 4 specialises the model to
protein-family inference (profile HMMs) and walks through HMMER and
#idx("Pfam")Pfam. Section 5 specialises it to gene finding and walks through GENSCAN
and #idx("AUGUSTUS")AUGUSTUS. Section 6 surveys the other places HMMs hide inside
the tools you already use — #idx("ChromHMM")ChromHMM, #idx("HISAT2")HISAT2's splice-aware aligner,
the pair-HMM inside HaplotypeCaller. Section 7 closes with the design
question that decides whether HMMs or a deep neural network is the
right tool for a new problem.


== From Markov Chains to Hidden State <sec:foundations>

A *Markov chain* over a sequence of observations $x_1, x_2, dots, x_T$
makes one structural assumption: the next observation depends only on
the current one, not on any earlier history. Formally,
$ P(x_t | x_(t-1), x_(t-2), dots, x_1) = P(x_t | x_(t-1)). $
This is a strong but not entirely useless assumption for #idx("DNA")DNA. A first-order
Markov chain on the four bases can capture, for example, the elevated
frequency of "CG" in CpG islands or the depletion of "TA" in coding
regions. Its limitation is that it has exactly one transition matrix.
There is no room in the model to say that the "CG" frequency is
elevated _here_ but background _there_.

The fix is to introduce a *latent state* $z_t$ at each position that
captures "which regime are we in?" and to let the observation be drawn
from a state-specific distribution:
- The latent variable $z_t in {S_1, dots, S_K}$ evolves as a Markov
  chain itself, with transition probabilities $a_(i j) = P(z_t = S_j | z_(t-1) = S_i)$.
- The observation $x_t$ is generated by the state-specific emission
  distribution $b_i(x) = P(x_t = x | z_t = S_i)$.

That gives the *Hidden Markov Model* of Baum and Petrie. The full set of
parameters is $theta = (pi, a, b)$, where $pi_i = P(z_1 = S_i)$ is the
initial-state distribution, $a$ is the $K times K$ transition matrix,
and $b$ is a $K$-by-alphabet matrix of emission probabilities. The
joint distribution over observations and states factorises into a
chain of transitions multiplied by a chain of emissions:
$ P(X, Z | theta) = pi_(z_1) product_(t=2)^T a_(z_(t-1), z_t) product_(t=1)^T b_(z_t)(x_t). $
Every algorithm in this chapter exploits that factorisation.

#figure(
  image("../../diagrams/lecture-21/01-hmm-graph.svg", width: 95%),
  caption: [
    A four-state HMM (top) and the same model unrolled in time (bottom).
    The unrolled graph makes the conditional independences explicit:
    each $x_t$ depends only on its own $z_t$, and each $z_t$ depends only
    on $z_(t-1)$. Algorithms exploit that chain structure to avoid the
    exponential blow-up of marginalising over all $K^T$ paths.
  ],
) <fig:hmm-graph>

The unrolled-in-time picture is the one to keep in your head. The
hidden states form a horizontal Markov chain; the observations hang
below them. The joint factorises into Markov transitions and
emissions, and the conditional independences in the graph determine
which messages each algorithm needs to pass.

=== Three Canonical Problems

Given the model, three questions recur. They are the same three
questions Rabiner's 1989 tutorial begins with, and they recur because
every concrete biological task reduces to one of them.

The first is the *likelihood*: given parameters $theta$ and an observed
sequence $X$, what is $P(X | theta)$? This is the quantity you compare
across alternative models — a CpG-island HMM against a background-only
chain, two competing gene-finder parameterisations, the trained model
against a permuted baseline. The naive computation sums over all
$K^T$ state paths and is exponential in sequence length. The forward
algorithm collapses it to $O(T K^2)$ by #idx("dynamic programming")dynamic programming.

The second is *decoding*: what is the most probable state sequence
$arg max_Z P(Z | X, theta)$? This is the segmentation question — given
my genome and my trained gene-finder HMM, where are the exons? The
Viterbi algorithm answers it in $O(T K^2)$ with max-product DP,
storing back-pointers to reconstruct the optimal path.

The third is *training*: given only observations, learn $theta$. The
Baum-Welch algorithm — expectation-maximisation specialised to HMMs
— iterates between using the current parameters to compute posterior
expectations of state occupancies and using those expectations to
re-estimate the parameters. It is the #idx("EM algorithm")EM algorithm of Dempster, Laird,
and Rubin (1977), specialised to a chain-structured latent variable,
and like all EM it converges monotonically to a _local_ optimum of the
likelihood.

#note[
  An HMM is a *discrete-state-space dynamical system* with stochastic
  transitions and stochastic observations — the discrete-symbol analog
  of the Kalman filter. #idx("forward-backward")Forward-backward is the
  discrete-time-discrete-state version of the Kalman filter and
  smoother. Viterbi is the maximum-a-posteriori decoder of the most
  probable #idx("trajectory")trajectory. Baum-Welch is EM for the latent state and
  parameters jointly. Baum and Petrie published the HMM in 1966; Kalman
  published the filter in 1960. The two communities developed
  independently for decades before the convergence in the 1990s. For
  an EE student, the punchline is that every HMM technique you will see
  here has a state-space-model analog you have already met in signal
  processing — and the algorithms transfer in both directions.
]

#figure(
  image("../figures/ch21/f1-state-space-analog.svg", width: 95%),
  caption: [
    The state-space-model family tree. The Kalman filter assumes a
    continuous latent state with linear Gaussian dynamics; the HMM
    swaps the continuous state for a finite alphabet and the Gaussian
    emissions for arbitrary discrete distributions. Forward-backward
    on the HMM is the structural analog of the Kalman filter-smoother
    pair; Viterbi corresponds to the maximum-a-posteriori trajectory
    decoder.
  ],
) <fig:state-space>

=== The CpG-Island Example

The canonical introductory example is CpG-island detection. The
cytosine in a CpG dinucleotide is the substrate for DNA
methyltransferase; methylated cytosine spontaneously deaminates to
thymine over evolutionary time, with the result that the human genome
contains only about 20% of the CpGs expected under independence of
adjacent bases. CpG _islands_ — short regions (typically 200–2000 bp)
of unmethylated, CpG-rich sequence — survive because they sit at
promoters and are kept unmethylated by transcription-factor binding.
About 70% of human gene promoters overlap a CpG island.

A two-state HMM segments the genome into island and non-island
regions:
- State 1 (island) emits the 16 dinucleotides with elevated CpG
  frequency — concretely, $P("CG"|"island") approx 0.05$ vs
  $P("CG"|"background") approx 0.01$.
- State 2 (background) emits the 16 dinucleotides with depleted CpG.
- Transition probabilities reflect the rarity and stickiness of
  islands: $P("island"|"island") approx 0.99$, $P("background"|"island") approx 0.01$.

Train the model on labelled regions, run Viterbi on a new #idx("chromosome")chromosome,
read off an island annotation. The pattern that this example
establishes — _few-state HMM plus Viterbi decoding equals genome
segmentation_ — recurs throughout the rest of the chapter. ChromHMM
uses 15-state versions of it. Gene finders use 30-state versions of
it. The state space changes; the decoding machinery does not.


== The Forward, Backward, and Viterbi Recurrences <sec:dp>

The three core algorithms share the same dynamic-programming spine.
Each defines a quantity at position $t$ and state $i$, gives a recurrence
that builds position $t+1$ from position $t$, and terminates with a sum
or a max over the final column of the table. The bookkeeping is
mechanical; the structure rewards careful indexing.

=== The Forward Algorithm

Define the forward variable
$ alpha_t(i) = P(x_1, x_2, dots, x_t, z_t = S_i | theta), $
the joint probability of having observed the prefix $x_1 ... x_t$ and
landed in state $S_i$ at time $t$. The recurrence is
- *Initialise:* $alpha_1(i) = pi_i  b_i(x_1)$.
- *Recurse:* $alpha_t(j) = b_j(x_t) sum_(i=1)^K alpha_(t-1)(i) a_(i j)$.
- *Terminate:* $P(X | theta) = sum_(i=1)^K alpha_T(i)$.

The recursion says: the probability of being in state $j$ at time $t$
having emitted everything so far is the emission probability $b_j(x_t)$
times the sum, over all states $i$ you could have been in at time
$t-1$, of the forward probability at $(t-1, i)$ times the transition
$a_(i j)$. Time complexity is $O(T K^2)$: $K$ recursions per cell,
$T K$ cells.

=== The Backward Algorithm

Symmetric to forward, the backward variable
$ beta_t(i) = P(x_(t+1), dots, x_T | z_t = S_i, theta) $
holds the probability of the _suffix_ given the current state. The
recurrence runs right-to-left:
- *Initialise:* $beta_T(i) = 1$ for all $i$.
- *Recurse:* $beta_t(i) = sum_j a_(i j) b_j(x_(t+1)) beta_(t+1)(j)$.

Multiplied together, $alpha$ and $beta$ give the per-position posterior:
$ P(z_t = S_i | X, theta) = (alpha_t(i)  beta_t(i)) / P(X | theta). $
This *posterior decoding* is what you want when you care about
per-position confidence — say, the probability that each base of a
gene lies inside an #idx("exon")exon rather than committing to the single most
probable exon-intron parse. Posterior decoding can produce state
sequences that are themselves never traversed in a single Viterbi path
because individually optimal positions need not be globally consistent;
that is a feature for many downstream uses and a bug for others.

#figure(
  image("../../diagrams/lecture-21/06-forward-backward.svg", width: 95%),
  caption: [
    Forward and backward messages on an HMM chain. The forward pass
    accumulates evidence from the left ($alpha_t$); the backward pass
    accumulates evidence from the right ($beta_t$). Their product at any
    node, normalised by the data likelihood, gives the per-position
    posterior $gamma_t(i)$. Mechanically this is the sum-product
    belief-propagation algorithm specialised to a chain graph.
  ],
) <fig:forward-backward>

=== The Viterbi Algorithm

Where the forward pass _sums_ over paths into each cell, Viterbi
_maximises_ over them. Define
$ delta_t(i) = max_(z_1, dots, z_(t-1))  P(z_1, dots, z_(t-1), z_t = S_i, x_1, dots, x_t | theta). $
The recurrence is structurally identical to forward, with one operator
swap:
- *Initialise:* $delta_1(i) = pi_i  b_i(x_1)$.
- *Recurse:* $delta_t(j) = b_j(x_t)  max_i [delta_(t-1)(i)  a_(i j)]$.
- *Back-pointer:* store $psi_t(j) = arg max_i  [delta_(t-1)(i)  a_(i j)]$.
- *Terminate:* read off $z^*_T = arg max_i delta_T(i)$ and trace back
  through $psi$ to recover the optimal path.

The DP table has $T K$ cells, each takes $O(K)$ work, so Viterbi too is
$O(T K^2)$.

The Viterbi recurrence is a *max-plus* DP in the log domain. Taking
logs:
$ log delta_t(j) = log b_j(x_t) + max_i [ log delta_(t-1)(i) + log a_(i j) ]. $
This is structurally identical to the log-domain pairwise alignment
recurrence from Chapter 2: pairwise alignment is Viterbi on a trivial
one-state-per-residue HMM with explicit gap states, and the gap
penalties are negative log transition probabilities. Once you have
implemented one of these algorithms, the other is a relabelling.

#figure(
  image("../../diagrams/lecture-21/02-viterbi-trellis.svg", width: 95%),
  caption: [
    Viterbi trellis with the backtrace highlighted. Each cell stores
    a maximum log-probability and a back-pointer to the predecessor
    cell that achieved it; the optimal state sequence is recovered by
    following back-pointers from the final column. The max-plus DP on
    the lattice is structurally the same machinery as the log-domain
    #idx("Smith-Waterman")Smith-Waterman recurrence of Chapter 2.
  ],
) <fig:viterbi-trellis>

=== Log Space, Underflow, and the Log-Sum-Exp Trick

A naive implementation of forward or Viterbi in floating-point
probabilities will underflow within a few hundred steps. Probabilities
in the range $10^(-4)$ to $10^(-2)$ multiply together; after a thousand
steps you are below the smallest representable double-precision number,
and your DP table fills with zeros. The fix is to work in log
probabilities throughout.

Viterbi in log space is mechanical because max distributes over addition:
$ log delta_t(j) = log b_j(x_t) + max_i [log delta_(t-1)(i) + log a_(i j)]. $
Forward is awkward because the forward recurrence sums products and
$log sum$ is not $sum log$. The standard trick is the
*log-sum-exp* identity:
$ log sum_i e^(L_i) = L^* + log sum_i e^(L_i - L^*), quad L^* = max_i L_i. $
Subtracting the maximum before exponentiating keeps every term in
$[0, 1]$ and avoids underflow; the additive constant cancels in the
final $log$. Every working HMM library does this. If you implement one
by hand and skip the trick, the symptom is that your log-likelihood
becomes $-infinity$ a few hundred bases into the sequence.

#figure(
  image("../figures/ch21/f2-log-space-viterbi.svg", width: 95%),
  caption: [
    Probability vs log-probability on a 200-position trellis. The
    left panel plots Viterbi cell values in native floating-point: the
    quantities fall below the smallest representable double around
    position 120. The right panel plots the same DP in the log domain
    with the log-sum-exp trick applied to forward sums: the values
    decay linearly with sequence length and never underflow.
  ],
) <fig:log-space>

#warn[
  Three numerical pitfalls trip nearly every new implementation. The
  first is probability underflow, which you fix by working in log
  space. The second is zero emission or transition probabilities, which
  freeze in maximum-likelihood training because the corresponding
  paths can never recover non-zero weight; smooth with pseudocounts or
  a Dirichlet prior. The third is the off-by-one in the back-pointer
  index. The remedies are old enough that any modern HMM toolkit
  (`hmmlearn`, `pomegranate`, `pyhsmm`) handles them transparently;
  resist the temptation to roll your own except for pedagogical purposes.
]


== Training: Baum-Welch and the Local-Optimum Problem <sec:training>

Given a labelled training set — pairs of observations and their true
state sequences — fitting an HMM is a single-pass exercise in counting.
Maximum-likelihood transition probabilities are observed transition
frequencies; maximum-likelihood emissions are observed emission
frequencies; the initial distribution is the observed first-state
frequency. Add a small pseudocount to every entry to avoid zeros, and
the result is a usable model in one matrix multiplication.

The problem genomics keeps running into is that the state sequences
are precisely what you do _not_ have. You have a megabase of DNA and a
guess about how many regimes it should be partitioned into. The
fitting problem becomes "given only the observations, maximise
$P(X | theta)$ over $theta$." This is the unsupervised training
problem, and it does not have a closed form. The classical solution
is the *Baum-Welch* algorithm, which is EM specialised to the HMM.

The iteration alternates between two steps. In the *E-step*, you treat
the current parameters as fixed and use forward-backward to compute,
for every position and every state, the posterior probability of
occupying that state and the posterior probability of every transition:
$ gamma_t(i) = P(z_t = S_i | X, theta^((k))) = (alpha_t(i)  beta_t(i)) / P(X | theta^((k))), $
$ xi_t(i, j) = P(z_t = S_i, z_(t+1) = S_j | X, theta^((k))) = (alpha_t(i) a_(i j) b_j(x_(t+1)) beta_(t+1)(j)) / P(X | theta^((k))). $
In the *M-step*, you treat these soft assignments as if they were
hard counts and re-estimate:
$ hat(pi)_i = gamma_1(i), $
$ hat(a)_(i j) = (sum_(t=1)^(T-1) xi_t(i, j)) / (sum_(t=1)^(T-1) gamma_t(i)), $
$ hat(b)_i (x) = (sum_(t: x_t = x) gamma_t(i)) / (sum_(t=1)^T gamma_t(i)). $
Iterate until the log-likelihood plateaus.

The theoretical guarantee Baum and Welch proved in 1972 is that
$log P(X | theta^((k+1))) >= log P(X | theta^((k)))$ at every iteration.
The likelihood never goes down. What the guarantee does _not_ say is
that the optimum it converges to is the global one. EM is local;
multiple random restarts are essential. A bad initialisation can leave
the algorithm stuck in a basin that does not match the data-generating
process at all — two states that should differ collapse onto each
other, or a third state captures noise that should have been pooled
with the dominant state.

#figure(
  image("../../diagrams/lecture-21/07-baum-welch.svg", width: 95%),
  caption: [
    Baum-Welch from five random initialisations on a synthetic
    three-state HMM. Each curve plots $log P(X | theta^((k)))$ versus
    iteration $k$. Three runs find the global optimum; two converge
    monotonically to inferior local optima. Multiple restarts and
    informed initialisation (e.g., from a brief $k$-means on
    observation windows) are standard practice.
  ],
) <fig:baum-welch>

#tip[
  In practice, five to ten random initialisations is the working
  default. A common trick is to seed initialisations from a quick
  unsupervised pass: $k$-means on overlapping windows of the
  observation sequence, or a Gaussian mixture on aggregate emission
  statistics. The seeded run usually beats fully random ones by a
  comfortable margin, and the cost is negligible.
]

=== Choosing the State Count

The number of states $K$ is a hyperparameter, not a learned quantity,
and Baum-Welch on too many states will happily overfit. The standard
tools for picking $K$ are penalised likelihood criteria:
$ "BIC"(K) = -2 log P(X | hat(theta)_K) + p_K log N, $
$ "AIC"(K) = -2 log P(X | hat(theta)_K) + 2 p_K, $
where $p_K$ is the parameter count (which grows roughly as $K^2$ from
the transition matrix plus $K times |"alphabet"|$ from emissions) and
$N$ is the sample size (typically the observation length). BIC penalises
extra parameters more aggressively than AIC and tends to recover the
right $K$ on simulated data with a known truth. AIC tends to err on
the side of too many states. Either way, you fit several values of $K$,
compute the criterion, and pick the minimum.

#figure(
  image("../figures/ch21/f3-state-count-bic.svg", width: 92%),
  caption: [
    Model-selection curves for HMM state count on a synthetic dataset
    generated from a true 5-state HMM. Log-likelihood (top) increases
    monotonically with $K$ as the model gets richer; BIC (middle)
    penalises extra parameters and reaches its minimum at the true
    $K = 5$; AIC (bottom) penalises more weakly and selects $K = 7$.
    BIC is the safer default for HMM model #idx("selection")selection.
  ],
) <fig:bic>

A more biology-aware alternative is to choose $K$ from a prior on the
number of regimes the genome contains. ChromHMM authors typically scan
$K$ from 8 to 25 and pick by visual inspection of the emission matrix:
states that collapse onto each other or that lack a clean biological
interpretation get pruned. The objective is not the lowest possible
BIC; it is the smallest model in which every state is interpretable.

=== Forward-Backward as Belief Propagation

There is one observation that ties the algorithms of this section
together with the algorithms of Chapter 20 on phylogenetics. The
forward-backward recurrence is exactly the *sum-product belief
propagation* algorithm specialised to a chain-structured graphical
model. The forward pass sends messages left-to-right; the backward
pass sends messages right-to-left; the product at each node is the
marginal posterior. #idx("Felsenstein")Felsenstein's pruning algorithm — the workhorse
of tree-based phylogenetic likelihood from Chapter 20 — is the same
sum-product BP on a _tree_ instead of a chain. Viterbi is the
*max-product* variant. The HMM and the phylogenetic likelihood differ
only in the underlying graph; the algorithm is the same.


== Profile HMMs and Protein-Family Annotation <sec:profile-hmms>

The first genomics application that gave HMMs their dominant place in
the toolbox was protein-family annotation. Anders Krogh, Saira Mian,
Kimmen Sjölander, and David Haussler published the profile-HMM
architecture in the _Journal of Molecular Biology_ in 1994. The idea
is to take a #idx("multiple sequence alignment")multiple sequence alignment of a known protein family —
the kind of #idx("MSA")MSA Chapter 19 covered — and turn it into a probabilistic
model that scores any new sequence against the family with full
position-specific resolution and explicit handling of insertions and
deletions.

=== The 3-State-Per-Column Architecture

The #idx("profile HMM")profile HMM has three kinds of states, one set per MSA column:
- A *match state* $M_k$ emits the amino acid at column $k$ with the
  column-specific frequency distribution.
- An *insert state* $I_k$ emits one or more amino acids _between_
  columns $k$ and $k+1$, allowing the query to be longer than the
  consensus. The insert state has a self-loop, so it can emit a run of
  inserted residues.
- A *delete state* $D_k$ silently skips column $k$, allowing the
  query to be shorter.

Transitions are tightly constrained. From $M_k$ you can go to $M_(k+1)$
(continue the match), $I_k$ (open an insertion), or $D_(k+1)$ (skip
the next match column). From $I_k$ you can go to $I_k$ (extend the
insertion), $M_(k+1)$, or $D_(k+1)$. From $D_k$ you can go to $M_(k+1)$,
$D_(k+1)$, or $I_k$. Two terminal states $B$ (begin) and $E$ (end)
flank the model.

#figure(
  image("../../diagrams/lecture-21/03-profile-hmm.svg", width: 95%),
  caption: [
    A profile HMM compiled from a five-row, eight-column MSA. Each
    match column has its own emission distribution; insert states sit
    between columns with self-loops for run extensions; delete states
    bypass match columns. Transition probabilities are estimated from
    the gap statistics of the input MSA. Viterbi alignment of any
    new sequence to this model gives a position-specific score and a
    sequence-to-model alignment.
  ],
) <fig:profile-hmm>

Match-state emissions are estimated by counting (with pseudocounts)
amino acids in the corresponding MSA column. Insert-state emissions
default to the background amino-acid distribution, since the inserted
residues are by construction not aligned to anything. Transition
probabilities are estimated from the observed transition counts in the
input MSA, again with pseudocounts.

The Viterbi alignment of a query sequence to the model gives a
position-specific score that is both a quality of fit and an
alignment. The forward score, which sums over all alignments rather
than picking the best one, is what HMMER's _local_ alignment mode
uses for sensitivity in the deep-homology regime, where no single
alignment dominates.

=== Position-Specific Filters with Gap States

The EE-flavoured reading of a profile HMM is as a *position-specific
#idx("matched filter")matched filter with explicit gap states*. A fixed-template matched
filter from radar or sonar correlates the input with a known signal
shape; the profile HMM does the same with two extensions. First, every
position of the template has its own emission distribution rather than
a single amplitude, so a position that varies across the family
(amino acid 50 of every known kinase is one of L, V, I, or M) gets a
spread distribution while a position that is invariant (the conserved
catalytic lysine) gets a near-Kronecker delta. Second, the gap states
allow length variability — the query can be longer or shorter than
the template within bounded slack.

Dynamic time warping in speech recognition is the immediate cousin.
DTW aligns a variable-rate speech utterance to a fixed phoneme template,
allowing local stretches and skips, and uses essentially the same
$O(T K)$ DP. The two communities developed the machinery in parallel
through the 1970s and 1980s and recognised the equivalence in the
1990s. If you have implemented DTW in a homework, you have implemented
the central recurrence of a profile HMM.

=== HMMER and Pfam

The dominant profile-HMM software is *HMMER*, written and maintained
by Sean Eddy since 1998. Four commands cover almost everything you
will want to do with profile HMMs in practice:
- `hmmbuild` compiles a profile HMM from an MSA.
- `hmmsearch` queries a profile HMM against a sequence database.
- `hmmscan` queries a sequence against a database of profile HMMs.
- `nhmmer` is the #idx("nucleotide")nucleotide variant, used for non-coding RNAs and
  regulatory #idx("motif")motif scanning.

HMMER's output report has four numbers per hit worth keeping straight.
The *sequence #idx("E-value")E-value* is the expected number of background sequences
that would score at least as well as the query, computed under a null
model of independent residues drawn from the background amino-acid
distribution. The *domain E-value* is the same statistic computed per
domain hit within the sequence. The *#idx("bit score")bit score* is the log-odds of the
model against the null in base-2 log; bit scores translate to E-values
through the database size. The *bias correction* downweights matches to
queries with unusual amino-acid composition (e.g., long proline runs or
highly hydrophobic stretches), in the same spirit as the #idx("BLAST")BLAST
compositional adjustment from Chapter 19. A sequence E-value below
$10^(-3)$ is the typical reporting threshold; below $10^(-10)$ is
unambiguous homology.

The companion database is *Pfam*, originally published by Bateman and
colleagues in 2002 and now hosted at the European Bioinformatics
Institute as part of InterPro. Pfam contains roughly twenty thousand
curated protein families, each described by:
- A *seed alignment* of representative sequences, manually curated.
- A *full alignment* of all UniProt sequences that match the model at
  $E < 10^(-3)$, computed automatically.
- A *profile HMM* compiled from the seed.
- Annotation: family name, function, GO terms, structural classification,
  reference papers.

When you assemble a new genome, predict open reading frames, and
translate them, Pfam scan is the default first-pass functional
annotation. About 80% of bacterial ORFs and 70% of eukaryotic ORFs get
at least one significant Pfam hit; the remainder go to BLAST or to
structural prediction.

#figure(
  image("../../diagrams/lecture-21/11-pfam-hit.svg", width: 95%),
  caption: [
    A 600-residue tyrosine kinase with three Pfam domains: an SH2
    binding module, the catalytic kinase domain, and a C-terminal PH
    domain. Each domain is an independent profile-HMM hit; the linker
    regions match no Pfam model. The combination of hits is the
    protein's functional architecture.
  ],
) <fig:pfam-hit>

The combination of HMMER and Pfam has been the default protein-family
annotation toolchain for twenty-five years. It is one of the cleanest
examples of HMMs out-performing more flexible alternatives because the
biological signal — position-specific amino-acid preferences plus
controlled gap structure — exactly matches the #idx("inductive bias")inductive bias the
profile-HMM imposes.


== Gene Finding: HMMs as Chromosome-Scale Tape Readers <sec:gene-finding>

The second canonical application of HMMs in genomics is gene finding.
The problem statement is straightforward: given an assembled
chromosome, annotate it with gene boundaries, exon-intron structure,
coding regions in the correct reading frame, untranslated regions,
and (sometimes) promoter regions. The hard cases are eukaryotic
genomes, where genes are sparse and fragmented — a typical human
gene occupies tens of thousands of base pairs of genomic DNA but
codes for a protein from only one or two thousand of them, with the
rest occupied by introns.

The eukaryotic gene-finding problem is naturally cast as HMM decoding
on a chromosome-length tape. The state space distinguishes
intergenic, untranslated, exon, intron, and splice-site regions; the
emissions encode the codon usage, splice-site sequences, and
CpG-island context that distinguish each regime; Viterbi on the
chromosome produces a full gene annotation in a single pass.

=== GENSCAN

The canonical statistical gene finder is *GENSCAN*, published by
Christopher Burge and Samuel Karlin in 1997 in the _Journal of
Molecular Biology_. GENSCAN is a generalised HMM — generalised
because it relaxes the geometric length distribution that a pure HMM
imposes on each state — with state classes for every major piece of
eukaryotic gene structure:
- *Intergenic* sequence with background statistics.
- *Promoter* and *5'-UTR* states.
- *Initial exon* with a start-codon emission at the boundary.
- *Internal exon* states.
- *Final exon* with stop-codon emission.
- *Single-exon gene* state for the small fraction of intronless genes.
- *Intron* states stratified by *phase* (0, 1, or 2): the codon
  reading frame can be split across an intron in three ways, and each
  phase has different splice-site statistics.
- *3'-UTR* and *polyadenylation* signal states.
- Mirrored states for the reverse strand, since genes can be encoded
  on either strand.

Each state has its own length distribution — geometric for introns
(which can be arbitrarily long), gamma-distributed for exons (which
have a strong peak around 130 bp), pointwise for splice signals.
Emission distributions encode codon-usage statistics, splice-site
consensus sequences, and the CpG-island enrichment near #idx("transcription")transcription
start sites. Viterbi decoding on this multi-class HMM segments a
hundred-megabase chromosome into a list of genes in minutes.

#figure(
  image("../../diagrams/lecture-21/04-genscan.svg", width: 95%),
  caption: [
    GENSCAN state architecture. Boxes are state classes; arrows are
    permitted transitions. Each class has its own length distribution
    (geometric for introns, gamma for exons) and its own emission
    model (codon-usage statistics, splice-site consensus, polyA signal).
    Viterbi on the full chromosome produces the gene annotation in a
    single pass.
  ],
) <fig:genscan>

The trained GENSCAN model — released in 1997, frozen ever since —
is still in heavy use. It is also still the canonical pedagogical
example of "HMM as chromosome-scale tape reader," and it is the
ancestor of every modern eukaryotic gene finder.

=== Splice-Site Detection

A sub-problem worth singling out is splice-site detection. The donor
(5') splice site at the start of an intron has the consensus dinucleotide
GT at intron positions +1 and +2, embedded in a roughly nine-base
context around it. The acceptor (3') splice site at the end of an
intron has the consensus AG at positions $-2$ and $-1$, with a
polypyrimidine tract (a run of C and T residues) extending 10–20
bases upstream.

These signals are detected by *position-specific weight matrices* — the
simplest case of a profile HMM, with no insert or delete states. A #idx("PWM")PWM
is a $L times 4$ matrix that scores each position of a fixed-length
window with the log-odds of seeing each base under the splice-site
model versus the background:
$ "score"(s) = sum_(i=1)^L log (P(s_i | "site", i)) / (P(s_i | "background")). $
The model parameters are estimated by counting bases at each position
in a training set of annotated splice sites, with pseudocounts.

#figure(
  image("../../diagrams/lecture-21/08-splice-pwm.svg", width: 95%),
  caption: [
    Donor and acceptor splice-site PWMs displayed as sequence logos.
    Bar height in bits at each position is the information content
    relative to background; tall stacks at GT (donor +1/+2) and AG
    (acceptor $-2$/$-1$) reflect near-invariant conservation.
    The polypyrimidine tract upstream of the acceptor is visible as
    elevated C and T preferences over a 10–20 base region.
  ],
) <fig:splice-pwm>

Modern splice-site prediction layers a deep convolutional network
(*SpliceAI*, Jaganathan et al. 2019) on top of the PWM-style models.
The #idx("CNN")CNN reads a kilobase or two of flanking context and outperforms
PWM-only models on the precision-recall frontier by integrating
long-range sequence features. The PWM remains the backbone inside
classical gene finders because of its interpretability and because
it gets most of the signal almost for free.

=== AUGUSTUS, #idx("BRAKER")BRAKER, and the #idx("RNA-seq")RNA-seq Era

The modern successor to GENSCAN is *AUGUSTUS*, written by Mario Stanke
and colleagues from 2003 onward. AUGUSTUS introduces three improvements
over GENSCAN that have made it the de facto standard.

First, AUGUSTUS supports a *conditional random field* mode in addition
to the classical HMM mode. A CRF (Lafferty, McCallum, Pereira 2001)
models $P(Z | X)$ directly rather than the joint $P(X, Z)$. The
discriminative formulation allows rich feature functions on the
observation sequence without the conditional-independence constraint
of an HMM's emissions — the score at each position can incorporate
arbitrary features of the surrounding genomic context, including
features that overlap or are correlated. The trade-off is that you
lose the ability to generate sequences from the model, which you
were not using in gene finding anyway.

Second, AUGUSTUS supports *hint integration*. External evidence —
RNA-seq read alignments, protein homologs, conservation tracks across
species — is folded into the scoring as additional features. A hint
that an RNA-seq junction read crosses a position elevates the
likelihood of the splice-site states at that position. The hint
mechanism is what makes AUGUSTUS responsive to experimental data
rather than purely sequence-driven.

Third, AUGUSTUS supports *iterative self-training*. Starting from a
related-species model, you predict an initial gene set on the new
genome, use the most confidently predicted genes as training data,
retrain, predict again, and iterate. The pipeline that wraps this
loop is *BRAKER* (Brůna et al. 2021), which combines AUGUSTUS, the
#idx("GeneMark")GeneMark family for the bootstrap, and the user's own RNA-seq for
hints. For a new mammalian or plant genome assembled today, BRAKER
plus AUGUSTUS plus RNA-seq is the standard annotation pipeline.

=== Prokaryotic Gene Finding

Bacterial and archaeal gene finding is substantially easier than
eukaryotic. There are no introns, gene density is high (typically one
gene per kilobase), and start codons are well-conserved. The
canonical tools — *GeneMark* (Borodovsky and McIninch 1993), *Glimmer*
(Salzberg et al. 1998), and *Prodigal* (Hyatt et al. 2010) — are all
HMMs with codon-frequency emissions and Markov-chain transition
models. They run in minutes on a typical five-megabase bacterial
genome and recover roughly 99% of annotated genes on benchmark
organisms.

The metagenomic case (mixed organisms in a single sample) is harder
because the codon-usage statistics differ between organisms, and the
training set cannot be assumed to be drawn from a single genome.
MetaGeneMark and Prodigal-meta train short-window codon models that
generalise across taxa, with a small accuracy penalty relative to
single-organism models.


== HMMs Hiding Inside Other Tools <sec:applications>

If gene finding and Pfam are the visible applications of HMMs, several
of the most heavily used tools in genomics use HMMs as silent
internal components. This section walks three of them.

=== ChromHMM and Chromatin-State Segmentation

*ChromHMM* (Jason Ernst and Manolis Kellis, _Nature Methods_ 2012) is
the standard tool for segmenting a genome into #idx("chromatin")chromatin states given
multi-track #idx("ChIP-seq")ChIP-seq data. The input is a per-genomic-bin binary or
count vector across several #idx("histone")histone modifications (typically nine to
twelve marks: H3K4me1, H3K4me3, H3K27ac, H3K36me3, H3K27me3, H3K9me3,
H3K9ac, H4K20me1, CTCF, and a few others). The model is a
multi-emission HMM where each chromatin state has a vector of
mark-specific emission probabilities, and the state space typically
covers fifteen or so biologically interpretable regimes.

The state names — _active promoter_, _strong enhancer_, _weak enhancer_,
_polycomb-repressed_, _heterochromatin_, _quiescent_ — are not
hard-coded. They are read off the emission matrix _after_ Baum-Welch
training: a state with high H3K4me3 and H3K27ac and low H3K27me3 is
labelled active promoter because that is what the histone-mark
combinatorics say it should be. The interpretability of the emission
matrix is what makes ChromHMM so popular: the user sees, for each of
the 15 states, the probability of observing each mark, and can read
off a biological meaning.

#figure(
  image("../../diagrams/lecture-21/09-chromhmm-emission.svg", width: 95%),
  caption: [
    The emission matrix of a 15-state ChromHMM model. Rows are states
    (active promoter through quiescent); columns are histone marks.
    Each cell is the probability of observing the mark given the state.
    The structure of the matrix is what defines the biological label
    of each state; after training, the user reads off the labels from
    co-occurrence patterns.
  ],
) <fig:chromhmm>

The ENCODE and Roadmap Epigenomics projects ran ChromHMM across
hundreds of cell types and tissues, producing the chromatin-state
maps that are the standard reference for regulatory genomics today.
The pipeline is straightforward: align ChIP-seq from your marks of
interest (Chapter 9), bin into 200-bp windows, threshold, run
ChromHMM with $K = 15$ or so, inspect the emission matrix, label
the states, project the annotation onto the genome.

=== Pair-HMMs Inside Variant Callers

GATK's *HaplotypeCaller* (Chapter 4) does indel-aware #idx("variant calling")variant calling
by re-aligning reads to candidate haplotypes with a three-state
*pair-HMM*. The pair-HMM has match, insert, and delete states, with
gap-open and gap-extend transitions controlling the indel-friendly
score:
- From match $M$ you can go to $M$ (continue match, probability
  $1 - 2 delta - tau$), $I$ (open insertion, probability $delta$),
  $D$ (open deletion, probability $delta$), or $E$ (end with
  probability $tau$).
- From insert $I$ you can extend the insertion ($I arrow I$ with
  probability $epsilon$) or close it ($I arrow M$ with probability
  $1 - epsilon$).
- Symmetric for delete.

#figure(
  image("../../diagrams/lecture-21/10-pair-hmm.svg", width: 90%),
  caption: [
    The pair-HMM at the bottom of GATK HaplotypeCaller and #idx("DeepVariant")DeepVariant.
    Three states (match, insert, delete) with gap-open probability
    $delta$ and gap-extend probability $epsilon$. Forward on this model
    gives the likelihood of a read given a candidate #idx("haplotype")haplotype with
    indels handled coherently. Viterbi gives an alignment.
  ],
) <fig:pair-hmm>

The forward score on the pair-HMM is the read-likelihood term in the
HaplotypeCaller's local realignment, and it propagates into the
#idx("genotype likelihood")genotype likelihood. The pair-HMM is the same machinery as the
profile HMM, applied between two sequences instead of between a
sequence and a model. DeepVariant runs a convolutional network on
top of pair-HMM-aligned pileups; the HMM contributes the indel-aware
re-alignment, the CNN contributes the classification head.

=== Splice-Aware Alignment in HISAT2

HISAT2 (Kim et al. 2015), the workhorse RNA-seq aligner from Chapter 5,
treats #idx("RNA")RNA reads as observations from a graph-extended HMM that allows
#idx("spliced")spliced alignments across introns. The base aligner is the
#idx("Burrows-Wheeler transform")Burrows-Wheeler transform from Chapter 2; the HMM extension is a
graph layer that adds intron-skip edges between exonic positions and
scores them by a junction likelihood combining donor and acceptor PWMs
with an intron-length prior. The HMM is the smallest part of HISAT2 by
code volume but the part that turns a DNA aligner into an RNA aligner.

The pattern generalises. *Viterbi-on-DAG* alignment — the alignment of
a long read or a query sequence to a #idx("pangenome")pangenome graph (Chapter 11) —
uses the same Viterbi machinery, on a directed acyclic graph instead
of a chain. States become graph positions; transitions become graph
edges. The chain-structured DP of this chapter is the special case
of the path-structured DP that the long-reads chapter develops in full.

=== What HMMs Cannot Do

HMMs make two assumptions that are worth being explicit about, because
they bound where HMMs are useful.

The *Markov assumption* says the current state depends only on the
immediately previous state. Long-range dependencies — say,
trans-splicing across a hundred-kilobase intron, or codon-context
correlations that span four or five codons — violate the assumption.
A high-order HMM can encode some longer-range context by enlarging
the state space, but the parameter count grows multiplicatively and
quickly becomes unfittable.

The *conditional-independence assumption* says emissions at different
positions are independent given the state sequence. Real dependencies
— between codons in the same gene, between #idx("methylation")methylation marks at
adjacent CpGs, between histone modifications on nearby nucleosomes —
are imperfectly captured by an HMM with diagonal emissions. The CRF
extension partly addresses this by conditioning on observations
without modelling them generatively, but the underlying inductive
bias is still local.

For genomic segmentation tasks, these assumptions are usually fine.
The state regimes are sticky (introns are long, chromatin domains
are large) and the emission signal is local. For tasks where
long-range dependencies dominate — long-range regulatory
interaction, alternative-splicing decisions that depend on
distant exonic enhancers — HMMs underperform deep architectures
that can route information further. The next section turns that into
a design rule.


== The HMM Niche, and How to Recognise It <sec:niche>

The history of the last decade is largely the history of deep neural
networks displacing HMMs from tasks where they were once dominant.
Splice-site prediction was a PWM-and-HMM problem until SpliceAI; TF
binding was a PWM problem until #idx("DeepBind")DeepBind; regulatory expression
prediction from sequence was a kernel-methods problem until #idx("Basenji")Basenji
and #idx("Enformer")Enformer (Chapter 16). The pattern suggests an obvious question:
should HMMs survive at all, or are they a transitional technology that
will be fully replaced as compute and data continue to grow?

The answer, after several years of empirical comparisons, is that
HMMs survive in a clearly defined niche and that the boundaries of
the niche are worth understanding as a design rule rather than as
nostalgia.

=== Where HMMs Still Win

Four properties together define the HMM regime.

*Small training data.* A handful of labelled examples — a new
pathogen's gene structure curated by a single researcher, a curated
seed alignment of forty proteins for a newly identified family — is
enough to train a robust HMM. The same data trains a deep network
poorly. The reason is the strong prior the HMM bakes in: the
state-space structure, the Markov dependence, the per-position
emission distribution are all hand-specified. A neural network with
millions of parameters has to learn those structural priors from
data, which it cannot do at a sample size of forty.

*Interpretable annotation.* Every HMM state has a name and an emission
distribution; the segmentation output assigns one state per position.
A user reading the output knows which positions were called intron and
why. The corresponding deep model produces a class probability per
position with no further explanation, and recovering interpretations
requires post hoc tools (attribution methods, #idx("attention")attention visualisation)
that are themselves contested. For applications that get audited —
clinical gene annotation, regulatory-region prediction in a research
publication — the interpretability is decisive.

*Sparse-state inference.* Most positions in a chromatin segmentation,
a gene annotation, or a domain annotation belong to one of a small
number of well-defined states; transitions between them are rare and
sticky. The HMM's prior — that adjacent positions usually share state
— exactly matches the data structure. Deep models can learn this
structure but pay a sample-complexity penalty for not assuming it.

*Deterministic decoding.* Viterbi returns a single optimal state
sequence — useful when the downstream pipeline wants a discrete
annotation, not a probability distribution. The deep alternative
returns soft probabilities that the user must threshold, often
producing fragmented or inconsistent annotations.

#figure(
  image("../../diagrams/lecture-21/05-hmm-vs-dl.svg", width: 92%),
  caption: [
    Accuracy as a function of training set size, schematic. The HMM
    curve rises early and saturates around $10^4$ labelled examples;
    the deep-learning curve starts lower but exceeds the HMM around
    $10^5$ and continues to improve with scale. The crossover region
    is where most genomic-classification tasks sit, and where the
    choice between HMM and deep model is a real design decision.
  ],
) <fig:hmm-vs-dl>

=== Where #idx("deep learning")Deep Learning Has Won

Three task types have moved decisively to deep learning.
*Splice-site detection* with long-range context: SpliceAI's 10 kb
receptive field captures sequence features that the local PWM cannot
see, and the precision-recall improvement is substantial.
*Transcription-factor binding* prediction from DNA: DeepBind, DeepSEA,
and the subsequent Enformer family rediscovered motifs as
first-layer filters and beat PWMs by exploiting context.
*Quantitative regulatory prediction* — CAGE expression, chromatin
accessibility, ChIP-seq intensity — from long-range sequence:
Enformer's two-hundred-kilobase receptive field is well beyond
anything an HMM can represent.

What these tasks share is *abundant labelled training data* (hundreds
of thousands of labelled examples available from ENCODE, GTEx,
TCGA) and *long-range dependencies* (the relevant features are
spread across thousands of bases). Both properties favour deep
networks: the data is plentiful enough to train the architectural
prior from scratch, and the dependencies are too long-range for the
HMM's chain structure to capture without enlarging the state space
absurdly.

#figure(
  image("../figures/ch21/f4-hmm-niche-quadrant.svg", width: 95%),
  caption: [
    A two-axis design diagram. The horizontal axis is the amount of
    labelled training data; the vertical axis is the strength of the
    downstream interpretability requirement. The four quadrants suggest
    different defaults: HMMs dominate the small-data /
    high-interpretability quadrant, deep models dominate the
    large-data / low-interpretability quadrant, and hybrid pipelines
    populate the other two.
  ],
) <fig:niche>

=== Hybrid Pipelines

The boundary between the two regimes is not sharp; most modern
production pipelines combine an HMM backbone with deep-learning
components. *BRAKER* combines AUGUSTUS's HMM with deep-learning
splice-site hints. *HMMRATAC* runs an HMM on #idx("ATAC-seq")ATAC-seq pileups with
deep-feature emissions. *DeepHMM* and friends keep the HMM topology
and replace the local emission models with shallow neural networks
that score richer features of the observation window.

The common pattern in these hybrids is that the HMM contributes the
*coherence prior* — adjacent positions belong to the same state — and
the deep network contributes the *local discrimination* — what does
this window look like? Each component does what it is best at. The
HMM's chain structure provides the smoothing that prevents the deep
classifier from issuing one-base flips between exon and intron calls;
the deep classifier's expressivity recovers signal that a simple
emission distribution would miss.

The design rule, condensed, is the following. If you have fewer than
about ten thousand labelled examples, the task admits a small,
biologically interpretable state space, and you care about producing
a deterministic annotation that an auditor can read, the HMM is the
right tool. If you have hundreds of thousands of examples, the
relevant signal is long-range, and the downstream consumer is another
machine that takes soft probabilities, the deep network is. If you
are in between — and most production genomics pipelines are — the
hybrid is.


== End-to-End: Annotating a New Genome <sec:pipeline>

The chapter so far has been algorithmic. This section walks an
end-to-end annotation pipeline of the kind a working bioinformatician
runs the first time a new assembled genome lands on disk, to give the
algorithms their concrete context.

=== A New Bacterial Genome

A typical 2025 bacterial annotation pipeline on a five-megabase
assembled genome runs as follows.
1. *Assemble* the reads to contigs (#idx("SPAdes")SPAdes for short reads, Flye or
   hifiasm for long reads). Output: a multi-contig FASTA file of
   chromosomal and plasmid sequences.
2. *Run Prodigal* for ORF prediction. The HMM-based prokaryotic gene
   finder produces about five thousand ORF predictions on a typical
   bacterial genome in a few minutes.
3. *Translate* the ORFs to protein sequences using the standard
   #idx("genetic code")genetic code (or the appropriate variant for archaea or
   organelle-targeted predictions).
4. *Run `hmmscan` against Pfam.* Each ORF is scored against the
   ~20,000 Pfam profile HMMs; hits with $E < 10^(-3)$ are recorded.
   Approximately 80% of ORFs receive at least one significant Pfam
   hit. Runtime is roughly one CPU-hour per thousand ORFs on a single
   modern core; the standard accelerated mode with `--cpu 8` brings
   that down to minutes.
5. *Run #idx("DIAMOND")DIAMOND BLAST against UniProt* for ORFs without strong Pfam
   hits. Adds function transfer for an additional 10% of ORFs.
6. *Map to #idx("KEGG")KEGG pathways* and #idx("Gene Ontology")Gene Ontology terms for metabolic and
   functional context.
7. *Output* a GFF3 annotation file plus a per-gene functional
   description.

The composite pipeline annotates 95% of the gene content of a typical
bacterial genome with definite function and another 4% with a
plausible function ("hypothetical protein, similar to ..."). Total
runtime on a modern laptop is under two hours.

#figure(
  image("../figures/ch21/f5-pipeline-laptop.svg", width: 95%),
  caption: [
    A complete bacterial-genome annotation pipeline with the
    laptop-scale runtimes that prevail in 2025. Each box has an
    input, an output, a tool, and an indicative runtime. Two of the
    boxes — Prodigal and `hmmscan` — are HMM-driven. The pipeline
    composes cleanly because each tool emits a well-defined data
    type that the next consumes.
  ],
) <fig:ch21-pipeline>

=== A New Mammalian Genome

A mammalian genome — three gigabases, on the order of twenty thousand
genes with strong intron structure — runs a related but heavier
pipeline.
1. *Assemble* with long reads plus #idx("Hi-C")Hi-C scaffolding (Chapter 11).
2. *Repeat-mask* with RepeatMasker, masking the roughly 40% of the
   genome that is repetitive sequence (LINEs, SINEs, LTRs, satellite
   DNA). Repeats are masked because gene finders trained on
   non-repetitive sequence over-call genes inside them.
3. *Bootstrap-train AUGUSTUS* from the gene set of a closely related
   species. A mouse genome bootstraps from rat; a primate bootstraps
   from human.
4. *Run BRAKER* with the user's own RNA-seq data as hints. The
   pipeline iteratively trains AUGUSTUS and a GeneMark variant on
   high-confidence predictions, refining the model with each pass.
5. *Final AUGUSTUS run* on the masked, repeat-masked, hint-guided
   genome. Output: roughly twenty thousand gene predictions for a
   typical mammalian genome.
6. *Validate* with `hmmscan` against Pfam for domain annotation.
   Most predicted proteins receive at least one Pfam hit; the
   remainder go to InterProScan or to BLAST against UniProt.
7. *Orthology assignment* with #idx("OrthoFinder")OrthoFinder or BUSCO against the
   nearest sequenced relatives, transferring gene names and
   functional annotation by orthology.

Total runtime is on the order of days on a server. The headline tools
in the pipeline — RepeatMasker, AUGUSTUS / BRAKER, HMMER —
are all HMM-based at their cores.

=== The 2025 Frontier

The frontier is hybrid. Several pieces of the classical HMM pipeline
have been augmented or partially replaced.

*SpliceAI* augments the splice-site PWMs inside AUGUSTUS with a deep
CNN that reads a 10 kb window of context. The integrated model
outperforms either component alone.

*AlphaFold-driven annotation* is now standard for the 10–20% of
predicted proteins that have no Pfam hit. The #idx("AlphaFold")AlphaFold structure
prediction is fed to Foldseek, which finds structural homologs in
the PDB and the AlphaFold database; function is inferred from
structural orthology rather than sequence orthology. This rescues a
substantial fraction of the "hypothetical protein" tail.

*Metagenomic gene callers* — MetaProdigal, Prokka — handle community
samples with mixed organisms. The HMM-driven core is the same; the
training corpus is broadened.

The architecture has not changed. HMMs still drive Pfam scanning,
gene prediction, chromatin segmentation, and indel-aware alignment.
Deep learning is layered on top, mostly to augment splice-site
scoring, motif scoring, and structural-orthology assignment. The
working bioinformatician's stack of 2025 is the HMM stack of 2005
with a layer of neural networks at the top and bottom; the middle
remains HMM-driven.


== Summary <sec:ch21-summary>

- Hidden Markov Models are *state-space dynamical systems on discrete
  observations*. Parameters are an initial-state distribution, a
  transition matrix, and per-state emission distributions; the joint
  $P(X, Z | theta)$ factorises along the chain and admits
  $O(T K^2)$ inference.
- The three canonical problems — *likelihood* (forward), *decoding*
  (Viterbi), *training* (Baum-Welch) — share the same DP spine.
  Forward sums over paths; Viterbi maximises over them; Baum-Welch
  alternates EM E-steps (forward-backward to compute posterior
  occupancies) with M-steps (re-estimate parameters from soft counts).
- Practical implementations work in *log space* with the log-sum-exp
  trick to avoid underflow, use *pseudocounts* to avoid frozen-zero
  parameters, and run *multiple random restarts* to mitigate the
  local-optimum failure mode of EM.
- *Profile HMMs* (Krogh et al. 1994) turn an MSA into a probabilistic
  family model with one match-state-per-column, one insert state,
  and one delete state. HMMER and Pfam are the dominant
  implementations; Pfam annotates roughly 80% of ORFs in a typical
  new bacterial genome.
- *Gene finders* (GENSCAN, AUGUSTUS, BRAKER, Prodigal) are
  chromosome-scale HMMs with state classes for intergenic, exon,
  intron, UTR, and splice signals. Viterbi on a chromosome produces
  a complete annotation in minutes.
- HMMs hide inside many other tools: *ChromHMM* for chromatin
  segmentation, *HISAT2* for splice-aware alignment, *GATK
  HaplotypeCaller* and DeepVariant via *pair-HMMs* for indel-aware
  read likelihoods.
- The HMM niche is *small labelled data, interpretable
  segmentation, sparse-state inference, deterministic decoding*.
  Deep learning has displaced HMMs where data is abundant and
  long-range dependencies dominate; HMMs remain dominant in
  segmentation, family annotation, and gene structure.
- The EE-framing thread runs throughout: HMM as Kalman analog,
  forward-backward as belief propagation, Viterbi as MAP trajectory
  decoding, profile HMM as position-specific matched filter with
  gap states, Baum-Welch as EM.


== Exercises <sec:ch21-exercises>

#strong[1.] #emph[CpG-island HMM.] Implement a two-state HMM with
dinucleotide emissions. Train it (with pseudocounts) on a 100 kb
human genomic region with known CpG-island annotations. Run Viterbi
on a held-out region; report precision and recall against the
ground-truth annotations. Plot the per-position posterior probability
of the island state along the region and compare to the Viterbi
hard call.

#strong[2.] #emph[Forward-backward by hand.] For a three-state HMM
with the parameters
$pi = (0.5, 0.3, 0.2)$,
$a = mat(0.7, 0.2, 0.1; 0.1, 0.8, 0.1; 0.1, 0.1, 0.8)$, and emissions
that produce symbol "1" with probabilities $0.9$, $0.1$, $0.5$ from
the three states respectively, compute $alpha$, $beta$, and $gamma$
by hand for the observation sequence $X = (1, 1, 0)$. Verify that
$P(X | theta) = sum_i alpha_T(i)$ matches the value from the
backward pass.

#strong[3.] #emph[Log-space implementation.] Implement Viterbi for a
five-state HMM in two versions: one in native floating-point and one
in log space with the log-sum-exp trick. Run both on a synthetic
sequence of length 1000. Verify they produce the same decoded path
when the native version still has numerical resolution, and identify
the sequence length at which the native version begins to fail.

#strong[4.] #emph[Baum-Welch local optima.] Implement Baum-Welch
training for a three-state HMM. Generate synthetic data from a
known generative process. Run training from twenty random
initialisations and plot the converged log-likelihoods. How often
does Baum-Welch find the global optimum? How does the answer
change if you seed initialisations from a $k$-means clustering of
overlapping observation windows?

#strong[5.] #emph[Profile HMM bit score.] Build a profile HMM with
`hmmbuild` from a 10-row, 50-column MSA of your choice (from Pfam's
seed alignments). Run `hmmsearch` against UniProt-Swiss-Prot.
Report the highest-scoring hit, the score gap to the next hit, and
the alignment to the model. Interpret the bit score as a base-2
log-odds against a background null.

#strong[6.] #emph[GENSCAN state space.] Sketch a GENSCAN-style HMM
state diagram for a hypothetical organism with three intron phases
and no UTRs. Count the states, count the non-zero entries in the
transition matrix, and estimate how many parameters the model has.
Comment on which parameters can be tied or shared to reduce the
total.

#strong[7.] #emph[ChromHMM emission inspection.] Download a published
ChromHMM emission matrix (the 15-state Roadmap Epigenomics model is
the standard reference). For three states of your choice, write one
sentence explaining the biological identity from the emission
pattern. Identify any state that is hard to interpret; speculate
why.

#strong[8.] #emph[(Open-ended.)] Pick a 2023 or 2024 published method
that explicitly replaces a classical HMM-based tool with a deep
neural network (e.g., a splice-site predictor, a chromatin-state
caller, a protein-family classifier). Read the evaluation section.
Identify what the deep method gains, what it loses, and which of
the four properties of the HMM niche the authors implicitly or
explicitly trade off. Write one page; cite the paper.


== Further Reading <sec:ch21-further-reading>

- *Rabiner, L. R.* (1989). "A Tutorial on Hidden Markov Models and
  Selected Applications in Speech Recognition." _Proceedings of the
  IEEE_ 77: 257–286. The canonical pedagogical reference; still the
  best place to start.
- *Baum, L. E., Petrie, T., Soules, G., and Weiss, N.* (1970). "A
  Maximization Technique Occurring in the Statistical Analysis of
  Probabilistic Functions of Markov Chains." _Annals of Mathematical
  Statistics_ 41: 164–171. The original EM-for-HMMs paper, predating
  Dempster, Laird, and Rubin's general EM by seven years.
- *Krogh, A., Brown, M., Mian, I. S., Sjölander, K., and Haussler, D.*
  (1994). "Hidden Markov Models in Computational Biology: Applications
  to Protein Modelling." _Journal of Molecular Biology_ 235:
  1501–1531. The profile-HMM paper.
- *Eddy, S. R.* (2011). "Accelerated Profile HMM Searches." _PLoS
  Computational Biology_ 7: e1002195. The HMMER 3 design paper; the
  practical reference for understanding HMMER output.
- *Burge, C., and Karlin, S.* (1997). "Prediction of Complete Gene
  Structures in Human Genomic DNA." _Journal of Molecular Biology_
  268: 78–94. The GENSCAN paper.
- *Stanke, M., and Waack, S.* (2003). "Gene Prediction with a Hidden
  Markov Model and a New Intron Submodel." _Bioinformatics_ 19:
  ii215–ii225. AUGUSTUS; pair with Brůna et al. 2021 on BRAKER.
- *Ernst, J., and Kellis, M.* (2012). "ChromHMM: Automating
  Chromatin-State Discovery and Characterization." _Nature Methods_
  9: 215–216. The ChromHMM paper; read alongside the Roadmap
  Epigenomics 2015 segmentation atlas.
- *Durbin, R., Eddy, S. R., Krogh, A., and Mitchison, G.* (1998).
  _Biological Sequence Analysis: Probabilistic Models of Proteins
  and Nucleic Acids._ Cambridge University Press. Still the
  definitive textbook treatment.
