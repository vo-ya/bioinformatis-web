#import "../theme/book-theme.typ": *

= Advanced Single-Cell: Trajectories, Integration, Multi-Modal <ch:advanced-single-cell>

#matters[
  Single-cell #idx("RNA")RNA sequencing solves the problem Chapter 7 set up: it
  hands you a measurement of every transcript in every cell of a tissue.
  But the measurement is a snapshot of a process that is fundamentally
  about change — cells differentiate, batches drift, modalities disagree,
  pixels mix. The methods in this chapter are the second-generation
  single-cell toolkit. They take the count matrix you already know how to
  cluster and ask harder questions: which cells came before which?
  Which gene is turning on right now? How do you compare two experiments
  that were never meant to be compared? What can you say about cells
  whose spatial positions you still know? Each of these is a textbook
  signal-processing problem dressed in biological clothing. Reading them
  that way is the only reasonable approach.
]

A scRNA-seq experiment ends with a cells-by-genes count matrix and a #idx("UMAP")UMAP
that looks, if you squint, like the cells were drawn that way on purpose.
Chapter 7 covered the path that lands you there: alignment, #idx("droplet")droplet
demultiplexing, normalisation, principal components, nearest-neighbour
graph, #idx("Leiden")Leiden clustering, marker genes. By the end of that pipeline you
have an annotated atlas — _T cells over here, monocytes over there, a
few rare populations in the corners_. This is enough to publish a
descriptive paper. It is not enough to ask any of the questions that
actually motivate the experiment.

The questions that motivate the experiment are dynamic, comparative, and
multi-modal. _Where are these cells going next?_ _How do these two
patients' samples compare?_ _Which of these cells are receiving signals
from which others?_ _Which #idx("Visium")Visium spot is mostly T cells and which is
mostly fibroblasts?_ Each of these requires going beyond the static
snapshot. This chapter walks through the second-generation single-cell
toolkit that answers them.

Six topics structure the chapter. Section 8.1 introduces #idx("pseudotime")pseudotime —
the trick that recovers temporal ordering from a single sampling moment
— and the three algorithms (Monocle3, #idx("Slingshot")Slingshot, #idx("PAGA")PAGA) that implement
it. Section 8.2 unpacks #idx("RNA velocity")RNA velocity, which reads the direction of
gene-expression change off the ratio of #idx("unspliced")unspliced to #idx("spliced")spliced reads.
Section 8.3 handles #idx("batch integration")batch integration, the supervised source-separation
problem that lets you merge datasets from different labs, days, or
chemistries. Section 8.4 expands beyond RNA to surface protein
(#idx("CITE-seq")CITE-seq) and #idx("chromatin")chromatin accessibility (scATAC-seq), and to the joint
embeddings that fuse them. Section 8.5 turns to #idx("spatial transcriptomics")spatial transcriptomics:
the platforms, the mixed-pixel #idx("deconvolution")deconvolution problem, and the
niche-detection analyses they enable. Section 8.6 closes with
#idx("ligand-receptor")ligand-receptor inference — what cells appear to be talking to which,
and why those "appears" are doing a lot of work.

Throughout, the EE framings keep their force. RNA velocity is a
state-space estimator with linear ODE dynamics. Batch integration is
blind source separation with a known mixing label. #idx("scVI")scVI is variational
EM, the same engine as Chapter 5's #idx("Salmon")Salmon. Spot deconvolution is the
remote-sensing mixed-pixel problem. None of these analogies is a
gimmick; each maps onto methods you would recognise from a different
domain, and each lets you reason about failure modes from a familiar
starting point.

The chapter assumes you have read Chapter 7. The count matrix, the
nearest-neighbour graph, and the UMAP appear here as inputs, not as
things to be derived.


== Pseudotime and #idx("trajectory")Trajectory Inference <sec:pseudotime>

scRNA-seq captures a single moment. The cells were dissociated, lysed,
encapsulated in droplets, reverse-transcribed, and sequenced — and by
the time you have data, every cell is gone. There is no second time
point, no follow-up measurement, no way to ask what an individual cell
became. Yet half the papers in the field talk about "developmental
trajectories," "differentiation paths," "cellular state transitions."
How?

The answer is _pseudotime_. The premise is straightforward: at any
instant, a population of differentiating cells contains some cells that
just started, some that are partway through, and some that are nearly
done. If the trajectory is continuous and the population is asynchronous,
the snapshot is implicitly a time-series in which every sample was
drawn at the same wall-clock moment but from different positions along
the trajectory. *Pseudotime* is the scalar coordinate that orders cells
along that implicit trajectory — by convention, zero at the start and
one at the end, in arbitrary units that count "progression" rather than
elapsed time.

#figure(
  image("../figures/ch08/f1-pseudotime-snapshot.svg", width: 95%),
  caption: [
    A snapshot of an asynchronous differentiating population contains,
    implicitly, the whole trajectory. Pseudotime is the position along
    the inferred curve — one scalar per cell, recovered from one
    sampling moment.
  ],
) <fig:pseudotime>

#note[
  Pseudotime is an ergodicity argument. If the cells in a population
  visit every state on the differentiation path with roughly the same
  frequency at which any single cell would visit those states over
  time, then the population distribution at one moment equals the
  time-distribution of any single cell. Replacing time-averaging with
  population-averaging is the same move ergodic theory makes in
  statistical mechanics, and the same one that lets a Doppler radar
  read velocity from a static-spectrum measurement of a rotating
  target. The validity of the move depends entirely on whether the
  ergodic assumption actually holds — and pseudotime, like ergodic
  theory, fails quietly when it does not.
]

=== Three Methods, Three Topology Assumptions

Three pseudotime methods dominate the field in 2024, each with a
different stance on what the trajectory topology looks like.

*Monocle3* (Trapnell et al. 2014, Cao et al. 2019). Fits a
*principal-graph skeleton* through the UMAP embedding using a
reverse-graph-embedding algorithm (`SimplePPT`). The skeleton is a
sparse graph — typically a tree or DAG — whose nodes are landmark points
in the embedding and whose edges trace the underlying trajectory. Every
cell is then projected to its nearest skeleton node, and pseudotime is
the graph distance from a designated "root" node along the skeleton.
Monocle3 handles branches natively because the skeleton itself can
branch; it is the most general option when the trajectory topology is
unknown _a priori_.

*Slingshot* (Street et al. 2018). Runs clustering first, then connects
cluster centroids into a minimum-spanning tree, then fits a smooth
*principal curve* through each lineage in the tree. Pseudotime is
arc-length along the principal curve, in coordinates the curve itself
defines. The result is much smoother than Monocle3's per-cell graph
distance — the principal curve is a regression, not a discrete
projection — and easier to interpret biologically. Slingshot is the
right choice when the clusters are distinct and the trajectory topology
is roughly tree-shaped.

*PAGA — partition-based graph abstraction* (Wolf et al. 2019). Builds a
graph at the *cluster level*: nodes are clusters from the Leiden step,
edges are weighted by the fraction of inter-cluster kNN edges in the
cell-level graph. PAGA does not by default produce per-cell pseudotime;
it produces a _topology summary_ saying which clusters are connected to
which others. The value of PAGA is as a first-pass sanity check before
committing to a specific trajectory method. If PAGA shows three
disconnected components, fitting a single principal curve through all
of them with Slingshot is biologically wrong — the data is telling you
the cells are not all on one trajectory.

#figure(
  image("../../diagrams/lecture-08/02-trajectory-methods.svg", width: 95%),
  caption: [
    Three trajectory-inference approaches on the same branching
    dataset. Monocle3 fits a per-cell principal-graph skeleton,
    Slingshot fits smooth principal curves per lineage, PAGA summarises
    topology at the cluster level. The three answers are not strictly
    comparable — they are different abstractions of the same data.
  ],
) <fig:trajectory-methods>

=== Roots, Branches, and Evaluation

Three practical issues recur in every trajectory analysis.

*Choosing a root.* Pseudotime orderings are direction-free by default:
the method gives you a curve, but the curve could run either way. You
have to nominate a root cell or root cluster, and the rest of the
ordering propagates outward from there. Two routes to a root: pick one
biologically (the cluster known to be the progenitor — haematopoietic
stem cells for blood-cell development, neural progenitors for brain
development), or pick one algorithmically (the cluster with highest
expression of early-stage markers, or — and this is one of the strong
arguments for the next section — the cluster with the lowest aggregate
RNA velocity flux).

*Branching.* Differentiation paths are rarely linear. At every lineage
decision point, cells commit to one of two or more fates, and the
trajectory method has to detect the branch and assign each cell to one
of the downstream branches (or to neither, if the cell is itself at the
branch point). Slingshot handles branches by fitting one principal curve
per lineage and assigning cells via weights; Monocle3 handles them
natively because the principal graph can branch.

*Evaluation.* Unlike clustering, trajectory inference is rarely
evaluated against a gold-standard ordering because true developmental
time is almost never known experimentally. When it _is_ — typically in
time-course experiments where cells were harvested at known points —
evaluation uses rank correlations (Kendall $tau$, Spearman $rho$)
between inferred pseudotime and true elapsed time. The Saelens et al.
(2019) benchmark, published in _Nature Biotechnology_, evaluated
forty-five trajectory methods against more than a hundred simulated and
real trajectories; it is the canonical reference for which method to
reach for first.

#warn[
  Pseudotime assumes the trajectory is _continuous_ and that cells are
  sampled _densely_ along it. If the underlying biological process has
  discrete state transitions — the cell-cycle phases G1, S, and G2/M
  as effectively independent states, or a sudden shock response that
  switches a population on in an all-or-nothing manner — fitting a
  continuous pseudotime will produce a smooth curve that has no
  biological meaning at all. The curve will look fine. It will not be.
  Always run PAGA first as a cheap topology sanity check; if PAGA
  shows isolated components, the data is telling you that a continuous
  trajectory model is the wrong abstraction.
]

#note[
  The original #idx("Monocle")Monocle paper (Trapnell et al. 2014, _Nature
  Biotechnology_) was among the first scRNA-seq methods papers to cross
  ten-thousand citations. The non-obvious insight — that a snapshot
  population can be read as a developmental timeline — opened up a
  sub-field in half a page. Monocle2 reformulated it via
  reverse-graph-embedding; Monocle3 (2019) scaled the approach to
  datasets above one million cells. Pseudotime is now so canonical
  that the originating insight has been absorbed into the field's
  default vocabulary.
]


== RNA Velocity <sec:velocity>

Pseudotime orders cells but does not tell you which way they are going.
A trajectory inferred by Monocle3 or Slingshot is direction-free until
you nominate a root. *RNA velocity* — La Manno et al. (2018), refined
into the #idx("scVelo")scVelo framework by Bergen et al. (2020) — closes the gap. It
reads the _direction_ of expression change directly off the data, without
needing a temporal label, by exploiting a small piece of biology that
scRNA-seq leaves on the floor: the difference between #idx("mRNA")mRNA that has
been spliced and mRNA that has not.

=== Spliced and Unspliced Reads

A eukaryotic gene is transcribed as a *pre-mRNA* that contains both
exons and introns. The pre-mRNA is then *spliced* — the introns are cut
out, the exons joined — to produce the mature mRNA that the #idx("ribosome")ribosome
translates. Splicing is fast compared to #idx("transcription")transcription, but it is not
instantaneous, and at any given moment a fraction of a cell's
transcripts have not yet been spliced.

In a 10x scRNA-seq run, both forms get captured. Cell Ranger separates
them at the alignment step (Chapter 7): reads aligning fully within
exons are *spliced*; reads that span an intron-exon boundary, or fall
entirely inside an #idx("intron")intron, are *unspliced*. A typical scRNA-seq run
ends up roughly eighty per cent spliced and twenty per cent unspliced,
with the ratio depending on the gene, the library prep, and the cell-cycle
state of the population.

La Manno et al. (2018) noticed something useful about this ratio. The
unspliced fraction is younger — it was transcribed more recently than
the spliced fraction, because every spliced molecule used to be an
unspliced molecule. So the unspliced-to-spliced ratio per gene encodes
information about whether transcription is currently increasing,
holding steady, or decreasing for that gene. A gene that is being
_turned up_ right now will have a high unspliced-to-spliced ratio (new
pre-mRNA arriving faster than old mature mRNA is degraded). A gene
being _turned down_ will have a low ratio (transcription has stopped
but mature mRNA lingers). A gene at steady state sits at the
characteristic ratio its splicing and degradation rates dictate.

#figure(
  image("../../diagrams/lecture-08/03-spliced-unspliced.svg", width: 95%),
  caption: [
    Spliced and unspliced reads are separated by alignment: reads
    within exons are spliced, reads crossing intron-exon boundaries
    are unspliced. The ratio across genes encodes the direction in
    which transcription is moving.
  ],
) <fig:spliced-unspliced>

=== The Dynamical Model

Treat each gene as a two-compartment chemical system. Let $u$ be the
abundance of unspliced transcript per cell and $s$ the abundance of
spliced transcript per cell. The continuous-time dynamics are

$ (d u) / (d t) = alpha(t) - beta u $
$ (d s) / (d t) = beta u - gamma s $

where $alpha(t)$ is the (possibly time-varying) transcription rate,
$beta$ is the splicing rate (pre-mRNA $arrow$ mature), and $gamma$ is
the degradation rate (mature $arrow$ gone). Most implementations
normalise $beta = 1$ and absorb the units into $alpha$ and $gamma$.

At steady state both derivatives vanish, so $u = alpha slash beta$ and
$s = beta u slash gamma = alpha slash gamma$. The steady-state ratio
$u slash s = gamma slash beta$ is constant per gene — a single line
through the origin in the $(s, u)$ plane.

The *velocity estimator* exploits exactly this geometry. For each gene,
plot all cells in the $(s, u)$ plane. Cells at steady state — those
neither inducing nor repressing the gene — lie along the steady-state
line. Cells _above_ the line have more unspliced transcript than steady
state would predict, which means the gene is being induced ($s$ is
about to rise). Cells _below_ the line have less unspliced than steady
state, meaning the gene is being repressed. The vertical residual

$ v_s = beta u - gamma s $

is the rate of change of $s$ per cell — the *spliced-mRNA velocity*.
Positive $v_s$ means the gene is being turned up at this cell;
negative $v_s$ means it is being turned down. Aggregate $v_s$ across
all genes (or, in practice, across a few hundred highly variable genes)
and the per-cell vector defines the direction in gene space that the
cell is heading.

#figure(
  image("../figures/ch08/f2-velocity-phase-space.svg", width: 95%),
  caption: [
    The RNA-velocity geometry in the spliced-versus-unspliced plane
    for one gene. The steady-state line has slope $gamma slash beta$;
    cells above the line are inducing, cells below are repressing.
    The vertical residual $beta u - gamma s$ is the spliced-mRNA
    velocity; aggregated across genes, it gives a per-cell direction
    of change in gene space.
  ],
) <fig:velocity-phase>

To visualise velocity in the UMAP that everyone is already familiar
with, the per-cell gene-space vector is projected to two dimensions by
computing, for each cell, the cosine similarity between its velocity
vector and the displacement vectors connecting it to its nearest
neighbours. The result is the familiar arrow-per-cell flow field on
top of the UMAP — every cell has a short arrow pointing toward its
inferred next neighbourhood.

#figure(
  image("../../diagrams/lecture-08/04-velocity-model.svg", width: 95%),
  caption: [
    RNA velocity as a state-space estimator. A two-compartment
    transcription-splicing-degradation model produces a steady-state
    ratio that fits $gamma$ per gene; the residual $beta u - gamma s$
    is the per-cell velocity; projection to UMAP yields the
    arrow-per-cell flow field.
  ],
) <fig:velocity-model>

#note[
  RNA velocity is a literal state-space estimator. The state vector of
  each cell is $(u, s)$ per gene; the dynamics are linear ODEs with
  known structure but unknown rate parameters; the observation is a
  noisy single-time-point measurement of the state. The estimator fits
  the rate parameters from steady-state behaviour — a linear regression
  in the phase plane — then evaluates the model's prediction
  $beta u - gamma s$ as the velocity. Every step has a direct analogue
  in Kalman-filter-adjacent texts: identify the state, model the
  dynamics, fit rate constants from observed regimes, propagate the
  fitted model. The biology-specific twist is that the "samples at
  different times" are different _cells_ at the same time rather than
  the same cell at different times, which is exactly the ergodic move
  pseudotime is built on.
]

=== velocyto and scVelo

Two implementations dominate.

*velocyto* (La Manno et al. 2018). The original. Implements the
*steady-state model* described above: fit $gamma$ per gene from the
upper-right and lower-left of the $(s, u)$ scatter, where the cells
are likeliest to be at steady state, then compute residual velocities
for every cell. Conservative, fast, biologically interpretable. Works
well for typical lineage data with clear induction and repression
regimes.

*scVelo* (Bergen et al. 2020). Extends the framework with two
additional models. The *stochastic model* augments the fit with
second-moment information — the variance across cells, not just the
mean — and gives better $gamma$ estimates in low-count regimes. The
*dynamical model* drops the steady-state assumption entirely, fitting
the full time-dependent $alpha(t)$ per gene along with $beta$,
$gamma$, and a per-cell latent time, in an EM-like iteration. The
dynamical model is strictly more expressive than the steady-state model
and strictly slower; it pays off on datasets with clear temporal
structure (embryonic development, directed differentiation) and is
overkill on steady-ish tissue atlases.

scVelo is the 2024 default. It integrates directly with the #idx("Scanpy")Scanpy /
AnnData stack, ships with sane defaults, and the dynamical-model
results are what most current single-cell velocity papers report.

=== Caveats and the "Hypothesis Generator" Posture

RNA velocity has well-documented failure modes, and the field has
settled into a careful posture about what its output can and cannot
support.

*Failure when assumptions break.* The steady-state model assumes $alpha$
is roughly constant and that the cells in the upper-right of the
$(s, u)$ plot are at steady state. If the true biology has fast,
oscillatory, or bursty transcription — which many systems do — $gamma$
is mis-estimated and the velocities point in wrong directions.

*Direction reversals from cell-cycle genes.* Cell-cycle genes
oscillate. Their velocity arrows can point backward at certain
cell-cycle phases, contaminating trajectories that have nothing to do
with cell-cycle progression. Standard practice is to regress out
cell-cycle scores before velocity analysis, or to restrict velocity to
non-cell-cycle highly variable genes.

*Projection artefacts.* Velocities live in gene space (thousands of
dimensions); the arrows you see on the UMAP are 2D projections via
cosine similarity to nearby cells. When the local UMAP geometry is
distorted — small populations, mis-tuned neighbourhood sizes — the
projected arrows can be misleading even when the gene-space velocities
are correct.

Bergen, Soldatov, et al. (2023, _Molecular Systems Biology_) documented
these failure modes systematically and proposed diagnostic tests. The
field's current consensus: RNA velocity is a useful hypothesis
generator, not a direct measurement of cell-state evolution. A velocity
arrow pointing from cell type A to cell type B is consistent with the
hypothesis that A differentiates into B, but it does not prove that
hypothesis; lineage tracing (with genetic labels or metabolic labelling
methods like 4sU SLAM-seq) is the gold standard for actually
establishing fate.

#warn[
  RNA velocity confidence drops in tissues with low cell turnover —
  liver, kidney, mature endothelium — because the unspliced-to-spliced
  ratio becomes dominated by noise when transcription has reached slow
  steady state. The stochastic model in scVelo accounts for ratio
  variance and helps somewhat in this regime. For tissues where
  velocity clearly is not working, direct lineage-labelling approaches
  (metabolic labelling, Zman-seq) are the right tools instead.
]


== Batch Integration <sec:batch>

Every scRNA-seq experiment has a *batch*: a specific capture run on a
specific day, with a specific reagent lot, by a specific operator, on a
specific sequencer. Batches systematically differ from each other in
ways the biology never cared about — small shifts in mean expression,
differences in capture efficiency, variation in dropout rate. The
#idx("batch effect")batch effect is real, large, and reproducible.

The problem becomes acute the moment you want to combine batches —
which you almost always do. Any interesting study has technical
replicates. Any clinical study pools patients. Any atlas project
integrates dozens of public datasets generated by labs that have never
spoken to each other. Naive concatenation fails on every one of these:
the UMAP of the concatenated data shows cells clustering _by batch_,
not by biology. Two parallel sets of clusters appear where one merged
set of cell types should be.

#figure(
  image("../../diagrams/lecture-08/05-batch-integration.svg", width: 95%),
  caption: [
    Batch effect before and after integration. Naive concatenation
    produces parallel batch-specific clusters; integration recovers
    cell-type structure across the two batches.
  ],
) <fig:batch-integration>

Formally, the count matrix factors as

$ X = f("biology", "batch", "technical noise") $

and integration is the task of estimating the biology component while
marginalising out the batch component. The two are mixed
non-trivially — there is no clean additive decomposition in the count
domain — so the methods that work have to be more sophisticated than
"subtract the batch mean."

#note[
  Batch integration is supervised source separation. Two (or more)
  signal sources — biological variation and batch nuisance — are
  mixed in the observed counts. The goal is to recover the biological
  component. This is structurally the same problem as blind source
  separation in audio (ICA decomposing mixed microphones into
  speakers), multi-channel interference rejection in communications
  (separating a desired signal from co-channel interference), or any
  "signal plus nuisance, both unlabelled" problem. The genomics
  twist is that the batch labels are _known_ — which makes the
  problem supervised source separation, strictly easier than the
  fully blind case.
]

The pre-history of this problem is worth a paragraph. Microarray studies
in the 2000s faced exactly the same difficulty: a chip's mean and
variance depended on the laboratory and on the lot of reagents, and
two studies of the "same" disease produced incompatible expression
profiles. #idx("ComBat")Combat (Johnson, Li, Rabinovic 2007) introduced empirical
Bayes location-and-scale adjustment per gene per batch — a linear
mean-shift correction with #idx("shrinkage")shrinkage toward a global prior. Combat was
the standard for a decade. It fails on single-cell data because the
batch effect is no longer just a linear shift in mean: capture
efficiency differs, dropout rate differs, and the same gene can show
opposite-sign batch effects in different cell types. Single-cell
batch integration is Combat's harder cousin.

=== #idx("Harmony")Harmony — Linear Correction in #idx("PCA")PCA Space

*Harmony* (Korsunsky et al. 2019, _Nature Methods_) is the fast,
simple baseline. It operates entirely in PCA space and iterates a
correction:

1. Run PCA on the concatenated data, as in Chapter 7.
2. Soft-cluster the PC embedding into $K$ groups via fuzzy $k$-means.
3. For each cluster, compute the per-batch centroid. For each cell,
   compute a correction vector that moves it toward the cluster's
   cross-batch centroid, weighted by its membership in the cluster.
4. Apply the corrections. Re-cluster on the corrected PCs. Repeat to
   convergence — typically ten to fifty iterations.

The output is a batch-corrected PC matrix, the same shape as the input,
ready to feed back into the UMAP and clustering pipeline from Chapter 7.
Harmony never touches the original counts; it only adjusts the
projection.

*Pros.* Fast — scales linearly with cell count and runs on a laptop on
millions of cells. Simple algorithm with one well-tuned hyperparameter
($K$, the number of soft clusters, which defaults to a function of the
batch count). Implemented in both #idx("Seurat")Seurat and Scanpy.

*Cons.* Linear correction in a linear space: if the true batch effect
is non-linear, Harmony cannot fix it. The method is also prone to
*over-correction* — when batch-specific populations are subtly
biologically different (different patients with different disease
states, for example), Harmony will sometimes fuse them along with the
batch noise.

Harmony is the right first pass for most datasets. It handles cosmetic
batch differences well, runs cheaply, and gives a baseline you can
compare against. For anything harder, reach for scVI.

=== scVI — Variational Autoencoders on Counts

*scVI* (Lopez et al. 2018, _Nature Methods_; Gayoso et al. 2022 for the
modern toolkit) is the deep-generative answer to the same question.
Architecturally it is a *variational autoencoder* with a
*negative-binomial reconstruction likelihood* — a direct callback to
the NB count model that turns up in every modern bulk #idx("RNA-seq")RNA-seq
#idx("differential expression")differential expression tool (Chapter 6).

The generative model:

$ z_i ~ cal(N)(0, I) $
$ rho_i = f_theta(z_i, s_i) $
$ x_(i j) ~ "NB"(rho_(i j) dot l_i, phi_j) $

where $z_i$ is a low-dimensional latent representation of cell $i$
(typically ten-dimensional), $s_i$ is the one-hot batch label,
$rho_i$ is the decoded mean expression profile, $l_i$ is a per-cell
library-size latent (learned jointly), and $phi_j$ is a per-gene
#idx("dispersion")dispersion. The encoder $q_phi(z bar.v x, s)$ is a small neural
network mapping the observed counts and batch label to a Gaussian
posterior over $z$; the decoder $f_theta(z, s)$ is another network
mapping the latent and batch label back to expected expression.

Training maximises the evidence lower bound (ELBO):

$ "ELBO" = EE_(q_phi)[log p(x bar.v z, s)] - "KL"(q_phi(z bar.v x, s) || p(z)) $

The first term is the reconstruction likelihood under the
negative-binomial model. The second is a Kullback-Leibler regulariser
pulling the posterior toward the unit-Gaussian prior $p(z) = cal(N)(0, I)$.

The batch-integration mechanic is in where $s_i$ enters the
architecture. The batch label is fed as an input to _both_ the encoder
_and_ the decoder. The encoder sees $(x, s)$ and produces $z$. The
decoder sees $(z, s)$ and reconstructs $x$. At inference time, the
latent $z$ is what you use for clustering and visualisation. Because
the decoder already has direct access to $s$, the latent $z$ does
_not_ need to encode batch identity in order to reconstruct $x$ — the
KL regulariser explicitly penalises any batch information that leaks
into $z$. The latent is, by construction, a biology-only representation.

#figure(
  image("../figures/ch08/f3-scvi-elbo.svg", width: 95%),
  caption: [
    The scVI architecture annotated with its variational-EM
    correspondence. Encoder is the E-step (compute posterior over
    latents); decoder plus NB likelihood is the M-step. Batch label
    conditioning frees the latent from carrying batch information.
  ],
) <fig:scvi-elbo>

#figure(
  image("../../diagrams/lecture-08/06-scvi-architecture.svg", width: 95%),
  caption: [
    scVI as a graphical-model schematic. Inputs are counts $x$ and
    batch label $s$; the encoder produces a Gaussian posterior over
    the latent $z$; the decoder uses $z$ and $s$ to predict expected
    expression $rho$; the likelihood is negative-binomial with a
    library-size latent and a per-gene dispersion.
  ],
) <fig:scvi-arch>

#note[
  scVI is variational EM on counts. The E-step is the encoder
  (compute posterior over latents given data); the M-step is the
  decoder-plus-likelihood update (maximise reconstruction given
  latents). This is directly the structure of the Salmon and #idx("Kallisto")Kallisto
  EM from Chapter 5 — read-to-transcript soft assignment was the
  E-step, transcript-abundance re-estimation was the M-step. The
  architectural upgrade: Salmon's E-step had a closed-form update;
  scVI's posterior over $z$ has no closed form, so the encoder
  approximates it with a neural network (this is _amortized
  inference_) and the parameters are optimised by gradient descent.
  The spiritual continuity is precise: both methods are iterative
  soft-assignment decoders, and the EM thread runs straight from
  classical bioinformatics into modern deep-learning genomics.
]

scVI's ELBO admits a clean reading: a reconstruction term that rewards
the decoder for matching observed counts, and a regulariser that
prevents the latent from over-fitting. Variants extend the framework in
specific directions. *totalVI* (Gayoso et al. 2021) adds a second
decoder head for CITE-seq protein counts. *scANVI* (Xu et al. 2021)
adds semi-supervised cell-type labels. *MultiVI* handles RNA plus ATAC
on the same cells. The entire scvi-tools family is built on the same
VAE backbone and the same conditional-decoder trick for batch
integration.

=== Evaluating Integration

Integration is a tradeoff. Over-correct, and you erase real biological
differences (e.g. cell types that genuinely differ between disease and
control get fused). Under-correct, and the batches still separate.
Quantifying where on this trade-off a particular result sits requires
metrics.

*LISI — local inverse #idx("Simpson")Simpson index* (Korsunsky 2019). For each cell,
compute the inverse Simpson diversity of batch labels in its $k$-nearest
neighbours; average across cells. A high *iLISI* (batch LISI) means
batches are well-mixed locally; a low *cLISI* (cell-type LISI) means
cell types are still well-separated. You want iLISI up and cLISI down.

*kBET — k-nearest-neighbour batch-effect test* (Büttner et al. 2019).
For each cell's $k$-NN neighbourhood, run a chi-squared test on the
batch composition against the global batch composition. The output is
the fraction of neighbourhoods that pass the test — higher is better.

*ASW — average silhouette width*. The standard clustering metric,
computed separately on batch labels (want it low: batches don't
separate) and on cell-type labels (want it high: cell types do
separate).

The *scib* benchmark (Luecken et al. 2022, _Nature Methods_) combines
all of these plus additional graph-, cluster-, and trajectory-based
metrics into a single integration score and evaluates methods across
diverse datasets. The 2022 verdict was that no single method dominated
across all scenarios: scVI and scANVI led on large complex atlases,
Harmony led on small simple batches, and methods that look good in
simulation often look worse on real biology.

#tip[
  Always run the un-integrated UMAP first. If batches mix naturally
  on the un-integrated data, you don't need integration — and running
  integration will only erase structure. Integration is a tool for
  the cases where it is necessary, not a default step. The scib team
  reports that for at least a fifth of public datasets, no integration
  outperforms the simple baseline of running PCA on log-normalised
  counts and stopping.
]


== Multi-Modal Single-Cell <sec:multimodal>

scRNA-seq measures mRNA. mRNA is not the only thing that matters about
a cell. Cells also have proteins on their surface, chromatin in their
nuclei, #idx("methylation")methylation marks on their #idx("DNA")DNA, and an ATP economy that drives
all of it. Each of these is its own modality, each can be measured by
its own assay, and each tells you something the others do not. The
multi-modal-single-cell stack measures more than one of them on the
same cells in a single experiment.

=== CITE-seq — RNA Plus Surface Protein

*CITE-seq* (Stoeckius et al. 2017, _Nature Methods_) measures mRNA and
surface-protein abundance simultaneously on the same cell. The trick is
disarming: before the cells enter the 10x chip, they are incubated with
a panel of antibodies, each carrying a unique DNA tag — the
_antibody-derived tag_, ADT. The antibody binds its target protein on
the cell surface. When the cell is captured and lysed inside a droplet,
the antibody tags are released and reverse-transcribed alongside the
mRNA, sharing the same #idx("cell barcode")cell barcode. A single sequencing run therefore
produces two matrices over the same cell-barcode list:

- The standard mRNA count matrix from Chapter 7 — cells by ~20,000
  genes.
- A protein matrix — cells by ~100–300 surface proteins, measured as
  ADT counts.

CITE-seq panels typically include the classical flow-cytometry markers:
CD4, CD8, CD19, CD3, CD16, CD56, and the rest of the immunologist's
default panel. The point of running it is that protein expression does
not perfectly track transcript expression. A cell with high CD4 mRNA
may not yet have CD4 protein on its surface; a cell with surface CD4
may have downregulated its transcript long ago. For cell-type
annotation in tissues like blood, protein is closer to the
flow-cytometry gold standard than transcript is. CITE-seq lets you
validate transcript-based annotations directly against protein, on the
same cells, in the same experiment.

#figure(
  image("../../diagrams/lecture-08/07-cite-seq.svg", width: 95%),
  caption: [
    CITE-seq measures RNA and surface-protein abundance on the same
    cells in a single run. The output is two matrices sharing a
    cell-barcode index: cells-by-genes for RNA, cells-by-antibodies
    for protein.
  ],
) <fig:cite-seq>

=== scATAC-seq — Open Chromatin Per Cell

*scATAC-seq* (Buenrostro et al. 2015; Cusanovich et al. 2015) measures
chromatin accessibility per cell. The assay uses the #idx("Tn5")Tn5 #idx("transposase")transposase,
which inserts sequencing adapters preferentially into open
(nucleosome-free) regions of DNA. Where Tn5 inserts, sequencing reads
pile up; where chromatin is closed, no insertions, no reads. Reads
cluster at promoters, enhancers, and other regulatory elements where
the genome is actively accessible to transcription factors.

The output format is unlike scRNA-seq. Instead of a cells-by-genes
count matrix, scATAC-seq produces a cells-by-peaks binary or
low-count matrix, where each peak is a 500 bp region of the genome
called from the aggregated insertion profile. Each cell typically
yields about ten thousand fragments distributed over two hundred
thousand peaks — extreme sparsity, far worse than scRNA-seq.

That sparsity changes the analysis pipeline. PCA does not work well on
binary sparse matrices; dimensionality reduction for scATAC-seq
typically uses *latent semantic indexing* (LSI), a TF-IDF-weighted SVD
borrowed from text retrieval. UMAP and Leiden follow as in Chapter 7.

*Multiome* chemistry (#idx("10x Genomics")10x Genomics Multiome) measures RNA and ATAC on
the same #idx("nucleus")nucleus in the same droplet. Each cell yields a gene expression
vector and a chromatin accessibility vector that are directly paired —
no inference needed to link them. Multiome is the cleanest multimodal
experiment short of additional protein measurement, and is rapidly
becoming the default for atlas-scale studies.

=== Joint Embeddings — WNN, totalVI, MOFA+

Given two (or more) modalities per cell, how do you produce a single
embedding that uses both? Three approaches dominate.

*Weighted Nearest Neighbours* (Hao et al. 2021, _Cell_; the algorithm
introduced with Seurat v4). Run PCA on each modality separately, producing
two $k$-NN graphs over the same cells. For each cell, locally estimate
the relative information content of the two modalities by checking how
well one modality's neighbours predict the other modality's neighbours.
Weight the two graphs by these per-cell reliabilities, fuse them into a
single weighted graph, cluster and UMAP the fused graph. No per-dataset
training. Implemented directly in Seurat.

*totalVI* (Gayoso et al. 2021). Extends scVI to handle RNA plus
protein jointly. Two decoders — one negative-binomial for RNA, one with
a protein-specific likelihood (a mixture of negative-binomial foreground
and background components, because antibody binding is leaky) — share a
single latent $z$. The latent encodes both modalities simultaneously,
and batch integration is handled by the same conditional-decoder trick
scVI uses.

*MOFA+* (Argelaguet et al. 2020). Multi-omics factor analysis. Each
cell is represented as a linear combination of factors; factors can be
modality-specific (a factor that loads only on RNA) or shared (a
factor that loads on both RNA and protein). The decomposition is
unsupervised, fast, and interpretable — factor loadings map directly
back to genes and proteins — but less expressive than a VAE for
high-dimensional data.

#note[
  WNN is multi-sensor fusion. Two sensors (RNA and protein) measure
  the same underlying state (cell identity) with different noise
  characteristics; the fused estimate weights each sensor by its
  local reliability. This is structurally identical to
  Kalman-filter sensor fusion in robotics, multi-sensor radar
  tracking, and any "combine two measurement channels by their
  local SNR" pipeline. totalVI does the same in a learned latent
  space rather than a $k$-NN graph — an autoencoder with two
  input heads sharing a bottleneck, trained to reconstruct both
  modalities from a shared latent.
]

The interesting cells in a joint embedding are not the ones where the
modalities agree — those are the easy cases. They are the ones where
RNA and protein _disagree_: a transitioning cell that has upregulated
CD19 mRNA but not yet accumulated CD19 protein, or a cell that has lost
CD19 mRNA but still carries surface protein from past expression. Joint
embeddings make these cells visible as outliers from the
"RNA = protein" axis, and tagging them is often what the experiment
was for in the first place.


== Spatial Transcriptomics <sec:spatial>

Every method up to this point loses the cells' positions. The
dissociation step that lets you sort a tissue into a single-cell
suspension is also the step that throws away which cell sat next to
which. The downstream analyses can recover lineage from velocity and
identity from clustering, but spatial relationships — _these
fibroblasts are clustered around this tumour_, _this neuron is in
cortical layer four_ — are gone.

*Spatial transcriptomics* keeps the positions. The assay measures
gene expression while preserving each cell's, or each region's,
coordinates in the original tissue section. Three platform families
dominate in 2024, and they make different trade-offs between #idx("coverage")coverage,
resolution, and throughput.

=== Three Platforms, Two Modalities of Measurement

*Visium* (10x Genomics) is the sequencing-based platform. A tissue
section is laid on a capture slide patterned with about five thousand
*spots*, each fifty-five microns in diameter and one hundred microns
centre-to-centre. Each spot contains capture probes for poly-A mRNA,
and so each spot captures messenger RNA from whatever cells happen to
sit on top of it — anywhere from one to twenty cells, depending on
tissue density. The captured mRNA is then sequenced in bulk and assigned
to its spot of origin by the spot's positional barcode. Output: a spots
-by-genes count matrix with about twenty thousand genes per spot. Every
gene is measured, but cells are not individually resolved — each spot
is a *mixed pixel* of the cells underneath it.

*#idx("MERFISH")MERFISH* (Chen et al. 2015), *seqFISH* (Lubeck et al. 2014), and
their relatives are the imaging-based platforms. Pre-designed panels
of one hundred to five hundred targeted genes are visualised in the
tissue by multiplexed fluorescent in-situ hybridisation: each gene
gets a binary barcode across multiple imaging rounds, the barcode
identifies the transcript, and individual mRNA molecules are localised
to within a few hundred nanometres. The output is per-cell or even
sub-cellular spatial coordinates for each mRNA. Throughput is limited
by microscopy speed; a typical MERFISH experiment covers about half a
million cells per slide.

*Xenium* (10x Genomics) and CosMx (NanoString) are the commercial
imaging platforms. Both run pre-designed panels of three hundred to
five hundred genes, both achieve sub-cellular resolution, both are
fully automated. They are direct competitors to MERFISH for commercial
labs and are increasingly the workhorse for translational spatial
transcriptomics.

*Stereo-seq* (BGI) sits between Visium and the imaging platforms — a
sequencing-based platform with much smaller capture spots (down to
about 500 nm), buying near-cellular resolution at the cost of
substantial library-prep complexity.

The key trade-off is fundamental and unlikely to disappear soon:
*sequencing-based platforms measure every gene but lose cellular
resolution; imaging-based platforms achieve cellular resolution but
measure only a panel.* Which one is right depends on whether you need
an unbiased survey (Visium) or a targeted interrogation (MERFISH,
Xenium). For atlas-scale tissue mapping, the field's current default is
to run both — Visium for transcriptome-wide coverage at lower
resolution, then a targeted MERFISH/Xenium panel of three hundred
informative genes at sub-cellular resolution on a parallel section.

#figure(
  image("../../diagrams/lecture-08/08-spatial-platforms.svg", width: 95%),
  caption: [
    Three spatial-transcriptomics platforms compared. Visium
    sequences mixed-pixel spots with full transcriptome coverage;
    MERFISH and Xenium image individual cells with limited gene
    panels.
  ],
) <fig:spatial-platforms>

=== The Mixed-Pixel Problem

Visium's resolution challenge is the *mixed-pixel problem*: each spot
contains multiple cells of potentially different types, so the spot's
measured counts are a weighted sum of contributions from each cell type
present underneath it. If you want per-cell-type expression — and you
usually do — you must *deconvolve* the mixed signal.

Formally, each spot $i$ has observed counts $x_i$. Each cell type $k$
has a reference profile $mu_k$, learned from a scRNA-seq atlas of the
same tissue. The spot's signal is approximately a mixture:

$ x_i approx sum_k w_(i k) mu_k $

where $w_(i k)$ is the proportion of spot $i$'s content attributable to
cell type $k$, subject to non-negativity ($w_(i k) >= 0$) and
sum-to-one ($sum_k w_(i k) = 1$). Solving for $bold(w)_i$ at each
spot is constrained non-negative regression with a sum-to-one constraint
— the same problem the geosciences have been calling *spectral
unmixing* for forty years.

#figure(
  image("../figures/ch08/f4-spot-deconvolution.svg", width: 95%),
  caption: [
    Visium spot deconvolution as constrained non-negative unmixing.
    The observed spot signal is modelled as a weighted sum of
    cell-type reference profiles; non-negative least squares with a
    sum-to-one constraint recovers the proportions.
  ],
) <fig:spot-deconvolution>

#figure(
  image("../../diagrams/lecture-08/09-spot-deconvolution.svg", width: 95%),
  caption: [
    The same idea worked through for a single spot containing three
    cell types. The observed count vector is the weighted sum; the
    inferred proportions match the underlying mixture.
  ],
) <fig:spot-deconv-schematic>

Two deconvolution methods dominate practice.

*RCTD* (Cable et al. 2021, _Nature Biotechnology_) — _Robust Cell Type
Decomposition_. Maximum-likelihood fit under a Poisson count model
with a sparse mixture prior. Two modes: a "singlet" mode that assumes
exactly one cell per spot (appropriate for high-resolution platforms
like #idx("Slide-seq")Slide-seq), and a "#idx("doublet")doublet" mode that allows up to two cell types
per spot (appropriate for Visium's mixed-pixel regime). Uses a
scRNA-seq atlas of the same tissue as the cell-type reference.

*cell2location* (Kleshchevnikov et al. 2022) — a Bayesian hierarchical
model with explicit spatial priors. Produces posterior distributions
over cell-type proportions per spot, not just point estimates. Slower
than RCTD but better suited to downstream spatial analysis that wants
to propagate uncertainty.

#note[
  Spot deconvolution is the remote-sensing *mixed-pixel problem*.
  A multispectral satellite pixel observes a weighted sum of
  contributions from multiple land-cover classes underneath it
  (soil, vegetation, water, urban surfaces); a spectral-unmixing
  algorithm solves for per-class fractional abundances under
  spectral-endmember constraints. The spatial-transcriptomics
  version replaces "land cover" with "cell type", "spectral
  endmember" with "cell-type expression profile", and "pixel" with
  "Visium spot". The mathematical structure is identical —
  constrained linear unmixing $x approx M bold(w)$ with
  $bold(w) >= 0, sum bold(w) = 1$ — and the same algorithmic
  toolbox (non-negative least squares, spectral unmixing under
  Bayesian priors) transfers directly. The reference profiles come
  from a scRNA-seq atlas rather than a lab spectral library;
  everything else is the same.
]

=== Spatial Niches and Neighbourhood Analyses

Once cells (or spots) carry both expression profiles and spatial
coordinates, a new family of analyses becomes possible.

*Niche detection.* Group cells or spots by the cell-type composition of
their spatial neighbourhood. Tumour-border niches typically contain
tumour cells, infiltrating immune cells, and fibroblasts; intestinal
crypts contain a specific stem-cell-to-differentiated-cell gradient;
glomerular niches in kidney contain podocytes plus mesangial cells plus
endothelium. Each of these is a recurring spatial pattern that the data
makes visible.

*Spatial differential expression.* For a given cell type, test whether
its expression depends on its spatial context. A macrophage near a
tumour expresses a different set of genes than a macrophage in
adjacent healthy tissue — both labelled "macrophage" by the standard
clustering pipeline, but spatially distinguished, with implications
for the underlying biology.

*Gradients and zonation axes.* Many tissues have continuous axes along
which expression varies — cortical layers in brain, periportal-to-
pericentral zonation in liver, crypt-villus in intestinal epithelium.
Spatial analysis can extract these axes directly from the data rather
than requiring them to be annotated _a priori_.

*Squidpy* (Palla et al. 2022) is the Scanpy-integrated spatial-analysis
toolkit; *Giotto* (Dries et al. 2021) is the R equivalent; *BANKSY*
(Singhal et al. 2024) performs spatially aware clustering that explicitly
incorporates neighbourhood composition into the cell representation.

#warn[
  A spatial analysis is only as good as its spatial registration.
  Slides from different experiments, different days, or different
  microscopes have to be aligned, and coordinate systems do not
  align trivially. Before drawing conclusions from a "spatially
  varying gene", always check that the variation is not just a
  slide-to-slide offset in probe efficiency. The pre-registration
  pipelines (SpaceRanger for Visium, Xenium Onboard Analysis,
  MERFISH's image-registration software) are non-optional
  infrastructure.
]

=== What Spatial Is Still Bad At

Spatial transcriptomics is new enough to be honest about its
limitations.

The resolution-versus-coverage trade-off is real. No current platform
gives sub-cellular resolution _and_ whole-transcriptome coverage _and_
high throughput at the same time. Study design forces a choice; the
right choice depends on the question.

Three-dimensional spatial transcriptomics is nascent. Most experiments
are 2D sections, and reconstructing 3D structure from serial sections
is an active research problem. Native 3D approaches exist —
light-sheet microscopy combined with barcoding — but remain
experimental and low-throughput.

Single-cell-plus-spatial requires either high-resolution imaging from
the start (MERFISH, Xenium) or a probabilistic-assignment step on top of
Visium. Pure Visium does not give you single cells, and the deconvolution
step is itself a significant source of uncertainty.

Batch effects across slides are analogous to scRNA-seq batch effects,
with the added complication of spatial autocorrelation. This is an
active area of method development, and current tools mostly carry over
scVI-style conditional-decoder logic with a spatial prior bolted on.


== Cell-Cell Communication <sec:communication>

Cells communicate by secreting ligands — cytokines, growth factors,
chemokines — that bind receptors on other cells. Given a scRNA-seq
dataset with cells already annotated to cell types, the natural
question is: which cell types are _signalling_ to which others, and
through which ligand-receptor pairs? The methods that answer this are
the inferential equivalent of looking at a population of dating-app
users and asking who might match.

The inference is fundamentally indirect, because scRNA-seq measures
transcripts, not ligand secretion or receptor binding. The usual
proxy: if cell type A highly expresses ligand $L$ and cell type B
highly expresses the cognate receptor $R$, infer that A _can_ signal
to B via the $(L, R)$ axis. Score each (cell-type-pair $times$
ligand-receptor-pair) combination by some function of the ligand's
expression in A and the receptor's expression in B — typically the
product, the geometric mean, or a permutation-tested
co-expression score.

Two tools dominate.

*#idx("CellChat")CellChat* (Jin et al. 2021, _Nature Communications_). Uses a curated
database of about two thousand ligand-receptor pairs (CellChatDB)
that includes known multi-subunit complexes — IL-2, for instance,
binds a three-subunit receptor, and CellChat's score uses the
geometric mean of all three subunits' expression rather than treating
the pair as a single ligand-receptor link. Significance is assessed
by permutation testing. Outputs a ranked list of communication axes
plus network-style visualisations.

*#idx("NicheNet")NicheNet* (Browaeys et al. 2020, _Nature Methods_). Goes one step
further: rather than scoring only $L$-to-$R$ co-expression, NicheNet
models the full $L arrow R arrow$ target-gene chain. A ligand is
inferred to act on a receiver cell population if the downstream target
genes of that ligand's receptor are differentially expressed in the
receiver cells. This tests for _functional_ signalling — the receiver
actually responded — rather than just receptor presence.

#figure(
  image("../../diagrams/lecture-08/10-ligand-receptor.svg", width: 95%),
  caption: [
    A cell-cell communication network inferred from per-cell-type
    ligand and receptor expression. Nodes are cell types, edges
    are ligand-receptor signalling axes, edge thickness encodes
    aggregate score, edge colour encodes ligand family.
  ],
) <fig:ligand-receptor>

The catch is that the scores look quantitative but rest on several
assumptions that are routinely false.

*Expression is not secretion.* A cell expressing a cytokine transcript
may or may not actually secrete protein. Many cytokines are
post-translationally regulated: cleavage from a precursor, requires a
specific stimulus, requires correctly assembled secretion machinery.
High transcript expression does not guarantee high protein output.

*Receptor expression is not active receptor.* Cells can express a
receptor transcript without surface presentation (stored in vesicles,
not trafficked), or with surface expression but functionally inactive
(desensitised, internalised after a previous binding event).

*Spatial information is usually ignored.* Most cytokines diffuse only
tens of microns; two cell types that co-express a ligand-receptor pair
are not communicating if they live in different tissue compartments.
Pure-scRNA-seq tools (CellChat, NicheNet) ignore this. Spatial-aware
tools (COMMOT, Giotto's interaction module) incorporate coordinates
when available and substantially improve the false-positive rate.

*The statistics are shallow.* Permutation tests in CellChat ask "is
this pair's score higher than chance under a label-shuffling null?",
but the null does not account for many biological sources of
correlation in expression. Published communication networks contain a
non-trivial fraction of false positives.

#note[
  Communication scores are a dating app's compatibility score —
  computed from data both people put into their profile. They tell
  you whether two people _might_ be a good match in principle.
  Whether they actually meet, speak, or form a connection requires
  evidence the app cannot see. The same discount applies to
  scRNA-seq ligand-receptor inference: the output is a list of
  axes that are plausible given expression, not a list of axes
  that are actually active _in vivo_.
]

The right posture is to treat ligand-receptor scores as a *ranked
hypothesis list*, not a list of answers. Filters to apply before
following a hit experimentally: restrict to ligands known to be
secreted (many co-expression hits involve transmembrane proteins that
do not signal across cells); require downstream target-gene
corroboration via NicheNet; cross-check with spatial co-location if a
spatial dataset is available; bias toward axes with prior biological
support and toward druggable receptors for translational value. The
goal is to triage three hundred hits down to the three or four worth a
month of wet-lab follow-up.


== Summary <sec:summary>

- Pseudotime turns a static snapshot into an ordering by assuming the
  population is asynchronous and the trajectory is continuous. Monocle3,
  Slingshot, and PAGA implement three different topology assumptions;
  the right method depends on whether the trajectory is a single curve,
  a tree, or a disconnected graph.
- RNA velocity is a state-space estimator. Spliced and unspliced read
  ratios, under a linear-ODE kinetic model, give a per-cell velocity
  vector in gene space. Project to UMAP for the arrow-per-cell flow
  field. Treat as a hypothesis generator, not as lineage proof.
- Batch integration is supervised source separation. Harmony is the
  fast linear baseline in PCA space; scVI is a VAE with a
  negative-binomial likelihood that learns a biology-only latent.
  scVI's variational EM is the modern descendant of the Salmon EM —
  both methods are iterative soft-assignment decoders.
- Multi-modal single-cell is sensor fusion. CITE-seq adds surface
  protein; scATAC-seq and Multiome add chromatin accessibility;
  totalVI, WNN, and MOFA+ produce joint embeddings that surface the
  cells where modalities disagree.
- Spatial transcriptomics preserves position. Visium covers the whole
  transcriptome at mixed-pixel resolution; MERFISH and Xenium cover
  panels at sub-cellular resolution. Spot deconvolution is the
  remote-sensing mixed-pixel problem in a different costume.
- Ligand-receptor inference is hypothesis triage, not measurement.
  Co-expression of $L$ and $R$ across two cell types is necessary but
  far from sufficient for actual signalling. The output is a ranked
  list to filter, not a directory of who-talks-to-whom.


== Exercises <sec:exercises>

#strong[1.] #emph[Pseudotime topology check.]
You have a scRNA-seq dataset of peripheral blood from a single patient.
After clustering you observe three well-separated clusters
corresponding to T cells, B cells, and monocytes. A colleague proposes
running Slingshot to "find the differentiation trajectory" connecting
them. In one paragraph: explain why this is a bad idea, what PAGA
would show, and what the right interpretation of the three-cluster
topology is.

#strong[2.] #emph[Velocity from a #idx("pileup")pileup.]
At one gene you have a population of 200 cells with average unspliced
count $u = 4.0$ #idx("UMI")UMI and average spliced count $s = 12.0$ UMI. The
steady-state line you fit through the data has slope
$gamma slash beta = 0.5$. (a) Compute the velocity $v_s = beta u -
gamma s$ with $beta = 1$. (b) Is this gene being induced or repressed
on average? (c) For a single cell with $u = 5$ and $s = 8$, is the
gene being induced, repressed, or at steady state for that cell? Show
the arithmetic.

#strong[3.] #emph[ELBO decomposition.]
The scVI ELBO is
$"ELBO" = EE_(q_phi)[log p(x bar.v z, s)] - "KL"(q_phi(z bar.v x, s) || p(z))$.
(a) Name the role of each term and which inference step it controls.
(b) What happens to the latent $z$ if you remove the KL term entirely?
(c) What happens to reconstruction if you scale the KL term by a very
large coefficient ($beta >> 1$, as in $beta$-VAE)? Connect each
behaviour to a corresponding failure mode in batch integration.

#strong[4.] #emph[Batch metric arithmetic.]
A 100-cell sample is divided into two batches (50 cells each) and three
cell types. After integration, for each cell you count batch labels in
its 10 nearest neighbours and average across cells. The result is
$"iLISI" = 1.85$. (a) What does the maximum-achievable iLISI equal
under perfect mixing of the two batches? (b) What does the minimum
iLISI equal? (c) Is 1.85 a good or a bad number, and why?

#strong[5.] #emph[Spot deconvolution.]
A Visium spot contains an unknown mixture of three cell types: T cell,
B cell, and macrophage. The reference scRNA-seq atlas gives mean
expression profiles $mu_("T"), mu_("B"), mu_("Mac")$ for two marker
genes — CD3 and CD20 — as
$mu_("T") = (100, 5)$,
$mu_("B") = (5, 80)$,
$mu_("Mac") = (10, 10)$
(in normalised UMI units). The spot's observed counts are $x =
(50, 35)$. Set up the constrained non-negative least squares problem
and solve for $bold(w) = (w_T, w_B, w_("Mac"))$ subject to $w >= 0$
and $sum w = 1$. (Hint: this is a 2-equation 3-unknown system —
expect a small simplex of solutions; pick the one minimising the
residual norm.)

#strong[6.] #emph[Triaging communication hits.]
CellChat reports the top five ligand-receptor pairs in a tumour dataset:
(1) IL-2 — IL2Ra/b/g chain (CD4 T cell $arrow$ NK cell),
(2) VEGFA — VEGFR2 (tumour $arrow$ endothelial),
(3) MS4A1 — MS4A1 (B cell $arrow$ B cell),
(4) IFNG — IFNGR (CD8 T cell $arrow$ macrophage),
(5) HLA-A — KLRD1 (tumour $arrow$ NK cell).
For each pair, decide whether to keep or discard it on the
basis of secretion plausibility and known biology. Identify the one
hit that is almost certainly a false positive and explain why.

#strong[7.] #emph[(Open-ended.)]
Pick one of the methods covered in this chapter — Monocle3, scVelo,
Harmony, scVI, totalVI, RCTD, CellChat, NicheNet — and read its
primary publication. In one paragraph, describe a single empirical
result from the paper that surprised you, and explain in one sentence
why a reader of this chapter would or would not have predicted it.


== Further Reading <sec:further-reading>

- *Trapnell, C., Cacchiarelli, D., Grimsby, J., et al.* (2014). "The
  Dynamics and Regulators of Cell Fate Decisions Are Revealed by
  Pseudotemporal Ordering of Single Cells." _Nature Biotechnology_ 32:
  381–386. The original Monocle paper — half a page that opened the
  trajectory-inference sub-field.
- *Saelens, W., Cannoodt, R., Todorov, H., and Saeys, Y.* (2019). "A
  Comparison of Single-Cell Trajectory Inference Methods." _Nature
  Biotechnology_ 37: 547–554. The canonical benchmark; read this before
  picking a pseudotime method.
- *La Manno, G., Soldatov, R., Zeisel, A., et al.* (2018). "RNA
  Velocity of Single Cells." _Nature_ 560: 494–498. The paper that
  introduced velocity.
- *Bergen, V., Lange, M., Peidli, S., Wolf, F. A., and Theis, F. J.*
  (2020). "Generalizing RNA Velocity to Transient Cell States Through
  Dynamical Modeling." _Nature Biotechnology_ 38: 1408–1414. The scVelo
  paper, including the dynamical model.
- *Bergen, V., Soldatov, R. A., Kharchenko, P. V., and Theis, F. J.*
  (2023). "RNA Velocity — Current Challenges and Future Perspectives."
  _Molecular Systems Biology_ 19: e11799. The follow-up that
  systematically documents velocity's failure modes.
- *Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., and Yosef, N.*
  (2018). "Deep Generative Modeling for Single-Cell Transcriptomics."
  _Nature Methods_ 15: 1053–1058. The scVI paper.
- *Korsunsky, I., Millard, N., Fan, J., et al.* (2019). "Fast,
  Sensitive and Accurate Integration of Single-Cell Data with Harmony."
  _Nature Methods_ 16: 1289–1296.
- *Luecken, M. D., Büttner, M., Chaichoompu, K., et al.* (2022).
  "Benchmarking Atlas-Level Data Integration in Single-Cell Genomics."
  _Nature Methods_ 19: 41–50. The scib benchmark.
- *Stoeckius, M., Hafemeister, C., Stephenson, W., et al.* (2017).
  "Simultaneous Epitope and Transcriptome Measurement in Single Cells."
  _Nature Methods_ 14: 865–868. The CITE-seq paper.
- *Hao, Y., Hao, S., Andersen-Nissen, E., et al.* (2021). "Integrated
  Analysis of Multimodal Single-Cell Data." _Cell_ 184: 3573–3587. The
  Seurat v4 / WNN paper.
- *Cable, D. M., Murray, E., Zou, L. S., et al.* (2021). "Robust
  Decomposition of Cell Type Mixtures in Spatial Transcriptomics."
  _Nature Biotechnology_ 40: 517–526. The RCTD paper.
- *Jin, S., Guerrero-Juarez, C. F., Zhang, L., et al.* (2021).
  "Inference and Analysis of Cell-Cell Communication Using CellChat."
  _Nature Communications_ 12: 1088.
- *scvi-tools documentation.* `docs.scvi-tools.org`. The reference
  for scVI, totalVI, scANVI, and MultiVI.
