#import "../theme/book-theme.typ": *

= Mass-Spectrometry Proteomics and Metabolomics: Counting What Actually Does the Work <ch:proteomics-metabolomics>

#matters[
  RNA-seq counts the messages a cell sends to itself. The proteins those
  messages encode are what actually carry out the work — enzymes catalyse
  reactions, transcription factors regulate gene expression, receptors
  detect signals, antibodies fight infection — and the metabolites those
  proteins produce are what most clinical assays end up measuring. The
  correlation between mRNA and protein abundance across genes is around
  $r = 0.4$ to $0.6$, low enough that a transcriptomic measurement is at
  best a noisy proxy for the protein-level state. This chapter is about
  the instrument that fills the gap. A modern mass spectrometer can
  measure the abundance of thousands of proteins from a drop of plasma,
  detect post-translational modifications that no sequencing assay can
  see, and identify the small molecules that proteins push around — and
  the bioinformatics of doing so is, when you squint, an EE problem in
  disguise. The whole chapter is matched-filter detection on
  one-dimensional signals, with a few twists.
]

A protein is a polymer of twenty amino acids that folds into a
three-dimensional shape and does something useful — catalyse a reaction,
bind a substrate, transport an ion, transmit a signal. The previous
chapters of this book have built up most of the machinery you need to
think about proteins computationally: how their sequences are encoded
(Chapter 1), how they are expressed (Chapters 5 and 6), how they fold
(Chapter 15), how they participate in networks (Chapter 22). What we
have not yet covered is how anyone _measures_ them at scale. Sequencing
gives you the message; the protein is the thing the message describes,
and reading the protein out of a cell is a separate engineering
discipline.

The dominant instrument for that engineering discipline is the
*mass spectrometer*. Built from physics that J. J. Thomson worked out
in 1912, scaled by twenty years of ion-source and analyser engineering,
and tied to bioinformatics by an algorithm that Eng, McCormack, and
Yates published in 1994, the mass spectrometer turns a protein mixture
into a list of identified peptides plus their abundances. The list is
then a count matrix not unlike the one RNA-seq produces — proteins by
samples instead of genes by samples — and the downstream statistical
machinery looks superficially similar. The differences are where the
interesting bioinformatics lives. Proteins do not have nucleotide-style
sequences to grep against the genome; they have peptides that have to
be matched against a database of theoretical fragment masses. Their
abundance dynamic range in plasma spans ten orders of magnitude, where
no single assay can cover the whole range. Their modifications —
phosphorylation, acetylation, ubiquitination — are biologically
load-bearing and analytically combinatorial. And the metabolites they
produce live in chemical space that has no sequence at all to organise
the search.

This is the final chapter of the book, and a closing note at the end
returns to the larger arc — the way RNA, protein, metabolite, and
function fit together in multi-omics integration, and where the field
appears to be heading. The bulk of the chapter is the proteomics and
metabolomics pipeline itself: how the instrument works, how the spectra
get identified, how the abundances get quantified, what the resulting
matrix looks like, and what the matrix is good for. The structure of
the pipeline is the structure of the chapter.


== From Spectrum to Peptide <sec:ch27-spectrum-to-peptide>

The protocol that produces most of the world's proteomic data is called
*bottom-up*, *shotgun*, or *LC-MS/MS proteomics*, and it has a fixed
shape. A protein mixture is digested with a protease into peptides;
the peptides are separated by liquid chromatography; the separated
peptides are ionised, fragmented, and mass-analysed; the resulting
fragment spectra are matched against a database to identify the peptides
they came from; and the peptide identifications are aggregated back to
their parent proteins. Bottom-up is the dominant workflow because the
intermediate object — the peptide, between five and twenty-five amino
acids — is the right size for both chromatography and mass-spectrometric
fragmentation. The alternative *top-down* workflow analyses intact
proteins and skips digestion entirely; it scales to a few thousand
proteoforms rather than tens of thousands of peptides and remains a
specialty technique for problems where the connection between
modifications on the same protein matters.

#figure(
  image("../../diagrams/lecture-27/01-lc-ms-ms.svg", width: 95%),
  caption: [
    The LC-MS/MS workflow end-to-end. Cells or tissue are lysed,
    trypsin digests the proteins into peptides, reverse-phase
    chromatography separates the peptides, electrospray ionises them,
    a tandem mass spectrometer measures intact and fragment masses,
    a search engine identifies the peptides and rolls them up to a
    protein-by-sample abundance matrix.
  ],
) <fig:lc-ms-ms>

=== The Instrument

Step one is *lysis and digestion*. Cells, tissue, or plasma are
solubilised in a denaturing buffer; disulfide bonds are reduced and
alkylated to keep cysteine residues from re-forming bonds; *trypsin*,
a protease purified from bovine pancreas, is added. Trypsin cleaves
proteins after lysine (K) and arginine (R) residues, with the rules
that an immediately following proline (P) blocks cleavage and that
missed cleavages happen at low frequency. K and R appear roughly once
every ten amino acids in a typical eukaryotic proteome, so tryptic
peptides are biased toward the five-to-twenty-five-residue length that
the rest of the workflow likes. The choice of trypsin is not
incidental — it falls out of physics. Tryptic peptides carry a charged
C-terminus (basic K or R) plus a free N-terminus, which means most of
them carry two positive charges in electrospray and fragment cleanly
along the backbone.

Step two is *liquid chromatography*. The digested mixture is loaded onto
a *reverse-phase* column — a fused-silica capillary packed with
hydrophobic C18-coated silica particles — and eluted with a gradient
running from an aqueous mobile phase to an organic one over one to three
hours. Peptides bind by their hydrophobicity and elute at characteristic
retention times. The result is a one-dimensional separation that spreads
roughly fifty thousand peptides over a few thousand seconds, presenting
them to the mass spectrometer a few hundred at a time rather than all at
once. Chromatography is what makes the rest of the workflow tractable:
without it, the instrument would see fifty thousand peptides
simultaneously and have no way to distinguish them.

Step three is *ionisation*. Peptides eluting from the LC column are
sprayed through a fine-bore needle held at a high voltage. The
combination of voltage, gas flow, and evaporation produces aerosol
droplets that, as the solvent evaporates, become charged peptides in
the gas phase. This is *electrospray ionisation* (ESI), the technique
that John Fenn at Yale demonstrated for large biomolecules in 1989
and that earned him a share of the 2002 Nobel Prize in Chemistry. The
older alternative — *matrix-assisted laser desorption ionisation*
(MALDI), introduced independently by Karas and Hillenkamp in Germany
and Koichi Tanaka in Japan in 1988 — uses a UV laser pulse to vaporise
a co-crystallised matrix, carrying the peptides into the gas phase.
ESI dominates because it couples cleanly to LC; MALDI persists for
specialised applications like imaging mass spectrometry, where spatial
resolution on a tissue section is the priority.

#note[
  Before ESI and MALDI, mass spectrometry was a tool for chemists who
  worked with small molecules — molecular weights below a thousand
  daltons — because larger molecules could not be transferred to the
  gas phase without falling apart. The 1988–1989 ionisation revolution
  is what made proteomics possible. A modern protein-grade mass
  spectrometer is the same physics that J. J. Thomson built into his
  first parabola spectrograph in 1912; what changed was how to get
  large fragile molecules into the instrument intact.
]

Step four is *mass analysis*. The ionised peptides are pulled into a
vacuum chamber and their masses are measured. A modern instrument
typically does this in two stages, hence the "MS/MS" or *tandem mass
spectrometry* name. The *MS1 scan* measures intact peptide masses with
high accuracy — better than five parts per million on an Orbitrap, the
analyser Alexander Makarov developed in the late 1990s and Thermo
commercialised in 2005. Intact-peptide mass alone is not enough to
identify a peptide: thousands of tryptic peptides share the same
nominal mass to within an instrument's resolution. So a *peptide
precursor* of interest is selected from the MS1 scan, accelerated into
a collision cell, smashed against an inert gas (most commonly via
*higher-energy collisional dissociation*, HCD), and the resulting
fragment ions are measured in the *MS2 scan*. The pair — MS1 mass plus
MS2 fragmentation pattern — is what identifies the peptide.

=== The Fragmentation Spectrum

When a peptide collides with the bath gas, the bond most likely to
break is the *amide bond* in the peptide backbone — the C–N bond
between adjacent residues. Cleavage produces two complementary fragments.
The N-terminal fragment retains the original N-terminus and is called
a *b-ion*; the C-terminal fragment retains the original C-terminus and
is called a *y-ion*. A peptide of length $n$ produces $n - 1$ possible
cleavage sites, so the full fragmentation spectrum can in principle
contain $n - 1$ b-ions and $n - 1$ y-ions, with masses determined by
which amino acids sit on which side of each cleavage point.

#figure(
  image("../figures/ch27/f2-peptide-fragmentation.svg", width: 92%),
  caption: [
    Backbone fragmentation of the tryptic peptide $sans("LSDPYHR")$.
    Each amide-bond cleavage produces a b-ion (N-terminal) and a
    y-ion (C-terminal). The full b/y series uniquely determines the
    sequence; consecutive y-ion mass differences read off the amino
    acid sequence from C-terminus to N-terminus.
  ],
) <fig:fragmentation>

The b/y series has a beautiful property the database search engines
exploit: consecutive y-ions differ in mass by exactly the residue
mass of one amino acid. Read the mass differences between consecutive
y-ions and you have read the sequence, one residue at a time, from
the C-terminus inward. The same trick works for b-ions from the
N-terminus. In practice, the spectrum contains gaps (some cleavages
produce no detectable ion), additional peaks (neutral losses of water
or ammonia, immonium ions, multiply-charged fragments), and noise. A
real peptide identification therefore relies on a statistical score
rather than a strict mass-difference walk.

#figure(
  image("../../diagrams/lecture-27/08-fragmentation.svg", width: 95%),
  caption: [
    Predicted fragmentation pattern for a longer peptide. The b-ion
    and y-ion series are shown above and below the sequence; the
    bottom panel illustrates an idealised stick spectrum with b-ions
    in cobalt and y-ions in red.
  ],
) <fig:frag-spectrum>

=== Database Search

The identification problem, in the abstract, is matching an observed
spectrum to one of millions of theoretical spectra computed from a
proteome. The standard approach is *target-database search*. Given a
reference proteome — for human work, the UniProt SwissProt human
proteome contains about twenty thousand reviewed entries — the search
engine performs *in silico* tryptic digestion to generate a candidate
peptide list, computes the theoretical b/y fragment masses for each
candidate at every plausible charge state, and scores each observed
MS2 spectrum against every candidate whose precursor mass falls within
the instrument's tolerance. The candidate with the best score becomes
the *peptide-spectrum match* (PSM).

#figure(
  image("../figures/ch27/f1-ms-history.svg", width: 100%),
  caption: [
    A century of mass spectrometry, with instrument and ionisation
    milestones above the timeline and the software and methods that
    made the instrument biologically useful below it. The
    1988–1989 ionisation revolution (MALDI from Karas-Hillenkamp and
    Tanaka; ESI from Fenn) is what made protein-scale mass
    spectrometry possible; everything below the line that follows is
    bioinformatics built on top of those two ideas.
  ],
) <fig:history>

The lineage of database search tools is part of the field's history,
sketched in @fig:history. *SEQUEST*, written by Jimmy Eng in John
Yates's lab and published in 1994, was the first practical
implementation. It scored candidates by
a cross-correlation between observed and theoretical spectra, which
remains the conceptual core of every modern search engine. *Mascot*,
released commercially by Matrix Science (Perkins, Pappin, Creasy, and
Cottrell, 1999), introduced a probability-based scoring scheme and
dominated commercial proteomics for a decade. *MaxQuant*, released
free of charge by Jürgen Cox and Matthias Mann at the Max Planck
Institute in 2008, bundled a search engine (Andromeda), an MS1-based
quantification module (MaxLFQ), and an analysis suite (Perseus) into a
single workflow that became the de-facto community standard.
*MSFragger*, developed by Andy Kong, Felipe Leprevost, Dmitry
Avtonomov, Sangtae Mellacheruvu, and Alexey Nesvizhskii in 2017,
implemented a fragment-index data structure that made open searches —
allowing arbitrary unspecified modifications — practical for the first
time, and ships inside the FragPipe pipeline. The Proteome Discoverer
suite from Thermo wraps SEQUEST and several other engines into a
commercial interface.

#tip[
  When picking a search engine in 2024, the practical hierarchy is
  roughly: MaxQuant for DDA bulk proteomics, FragPipe/MSFragger for
  open searches and PTM-heavy applications, DIA-NN for DIA, Mascot
  if you have a vendor-supported pipeline. The differences between
  the top tools at default settings are smaller than the differences
  between any of them and a poor parameter choice. The single most
  expensive mistake is forgetting to enable carbamidomethylation as a
  fixed modification on cysteine — the standard sample prep alkylates
  cysteines, and missing it shifts every cysteine-containing peptide
  by 57 Da.
]

=== Peptide Identification as Constrained Matched Filtering

Strip the chemistry away and the identification problem has a familiar
EE shape. An MS2 spectrum is a one-dimensional signal — intensity
versus mass-to-charge ratio, sampled at some resolution determined by
the analyser. A candidate peptide produces a theoretical *template*:
a set of expected peak positions with predicted intensities. The
search engine computes a correlation or a probability score between
template and signal, and picks the candidate with the best score.

This is *constrained matched filtering*. The classical matched filter
in radar detection is optimal when both the signal template and the
noise covariance are known. In peptide identification, the templates
are the theoretical spectra of all tryptically possible peptides, and
the noise is everything else — chemical noise, electronic noise,
co-fragmented contaminants, instrument artifacts. The constraint —
that the candidate be a tryptic peptide from a known protein — shrinks
the search space from $20^n$ (all possible $n$-residue sequences) to
the few million tryptic peptides predictable from a reference proteome,
which is what makes the search computationally tractable. Without
the constraint, *de novo* sequencing from MS2 alone is possible but
much harder, and remains a specialty technique.

#note[
  Cross-correlation scoring in SEQUEST is exactly the discrete
  cross-correlation an EE student knows from signal processing:
  $X("score") = sum_(m) S(m) T(m)$ over the discretised spectrum
  $S(m)$ and theoretical template $T(m)$, normalised to remove the
  bias of common peaks across all templates. Mascot's probability
  score adds a likelihood model on top. Both reduce to "the better
  the template matches, the higher the score" — exactly the principle
  the matched filter encodes.
]


== Acquisition Modes: DDA, DIA, and Targeted <sec:ch27-acquisition>

The "MS/MS" half of LC-MS/MS hides a design choice the instrument has
to make every few seconds: from the dozens or hundreds of peptides
co-eluting at any moment, which ones do you fragment? The answer
defines the acquisition mode, and the three canonical answers —
*data-dependent acquisition* (DDA), *data-independent acquisition*
(DIA), and *targeted acquisition* (MRM/PRM) — trade off the same axes
in different ways.

=== Data-Dependent Acquisition

DDA is the classical mode and was the only practical option for two
decades. The instrument runs a fast MS1 scan, picks the top $N$ most
intense precursors — typically $N = 10$ to $20$ — and fragments each
in turn before going back to MS1. The selection rule is dynamic: every
few seconds the instrument re-examines the current MS1 landscape and
chooses what to fragment next. Most pipelines also apply a *dynamic
exclusion* window of fifteen to thirty seconds, so that a precursor
once fragmented is not re-picked until enough time has passed that it
is likely to have eluted.

DDA's strength is that the resulting MS2 spectra are mostly clean —
each spectrum contains fragment ions from a single peptide, give or
take some co-isolation contamination. Its weakness is *stochastic
coverage*. Because the instrument can only fragment a fixed number of
precursors per unit time, and because the precursor selection is
biased toward whatever happens to be most intense at that moment, the
set of peptides identified varies between technical replicates of the
same sample. Reproducibility of peptide identification across DDA
replicates is around 70%, which is fine for discovery work and a
problem for cohort-scale studies where missing values pile up.

=== Data-Independent Acquisition

DIA dispenses with the precursor-selection question entirely. Instead
of picking which peptides to fragment, the instrument *systematically*
fragments every precursor in a defined window of $m/z$ values, then
moves to the next window, then the next, cycling through the full
mass range every few seconds. A typical DIA method covers 400 to 1,200
$m/z$ in forty 20-Da windows; some methods use variable-width windows
calibrated to keep the precursor density approximately constant. The
SWATH-MS method, published by Ruedi Aebersold's group in 2012, was the
breakthrough that made DIA practical on Orbitrap and time-of-flight
platforms.

The consequence is that DIA acquires MS2 spectra deterministically:
every peptide in every window gets fragmented in every cycle. The
identifications are reproducible across replicates — typically 85% or
better — and the quantification CV drops from around 15% for DDA-LFQ
to around 5% for DIA-LFQ. The price is *chimeric spectra*: each
DIA MS2 spectrum contains fragment ions from every peptide that happened
to be in the 20-Da precursor window at the moment of acquisition, which
might be five, ten, or twenty peptides at a time. Deconvolving the
mixture is a sparse-recovery problem that the search engine has to
solve.

#figure(
  image("../../diagrams/lecture-27/02-dda-vs-dia.svg", width: 95%),
  caption: [
    DDA picks top-$N$ precursors and fragments them sequentially; DIA
    cycles through fixed $m/z$ windows and fragments every precursor
    in each window. The cost of DDA is stochastic coverage; the cost
    of DIA is chimeric spectra that need deconvolution.
  ],
) <fig:dda-vs-dia>

#figure(
  image("../figures/ch27/f3-dda-vs-dia-coverage.svg", width: 95%),
  caption: [
    Identification overlap across three technical replicates under DDA
    (top, top-$N = 15$) and DIA (bottom, forty 20-Da windows). DDA
    leaves a stochastic ring of peptides identified in one replicate
    but not the next; DIA's deterministic acquisition closes most of
    that ring.
  ],
) <fig:ch27-coverage>

The DIA software ecosystem matured around 2020. *Spectronaut*, a
commercial product from Biognosys, and *DIA-NN*, an open-source tool
from Vadim Demichev's group (Demichev, Messner, Müller, Tolić-Nørrelykke,
and Ralser, _Nature Methods_ 2020), now dominate. Both use neural
networks trained on observed DDA spectra to predict the fragmentation
patterns and retention times of theoretical peptides, producing
*in silico* spectral libraries that DIA searches against. The
combination — *library-free DIA* — was one of the most consequential
methodological shifts of the past decade and is the reason DIA is
now standard for cohort-scale proteomics.

=== Targeted Acquisition: MRM and PRM

When the question is not "what proteins are in this sample" but "how
much of these specific peptides is in this sample", a third mode
applies. *Multiple Reaction Monitoring* (MRM, sometimes called SRM)
runs on a triple-quadrupole mass spectrometer and monitors specific
precursor-to-fragment-ion transitions. The first quadrupole filters
to a single precursor mass; the third quadrupole filters to a single
fragment mass; only ions matching both pass through to the detector.
The selectivity is enormous and the quantitative dynamic range is the
best of any MS technique — five orders of magnitude — but the
throughput is low, because each peptide of interest costs an analyser
cycle. *Parallel Reaction Monitoring* (PRM) is the modern variant on
high-resolution instruments: a single high-resolution MS2 scan per
precursor of interest captures all fragments at once, with comparable
selectivity. PRM is the workhorse for clinical biomarker quantification
where absolute concentrations matter — measuring troponin, BNP, or a
half-dozen tumour markers in a clinical-research setting.

=== DIA as Compressed Sensing

The mathematical structure of DIA identification is worth a paragraph
of its own, because it is structurally identical to a problem that
EE students meet in signal processing. Each DIA MS2 spectrum
$bold(y) in bb(R)^M$ is approximately a linear combination of
theoretical peptide spectra: $bold(y) approx A bold(x)$, where the
columns of $A$ are the theoretical spectra of the candidate peptides
in the precursor window and $bold(x)$ is the vector of peptide
abundances. The matrix $A$ has many more columns than rows (more
candidate peptides than observed peaks), so the system is
underdetermined. Most of $bold(x)$ is zero — only a handful of
peptides actually co-fragment at any moment — so the recovery is a
sparse-recovery problem.

This is *compressed sensing*, as formalised by Donoho and by
Candès–Romberg–Tao in 2006. The solution is an
$ell_1$-regularised regression, a greedy selection algorithm, or a
neural network trained to perform the same task. DIA-NN does the
neural-network version. The pattern recurs throughout this book —
RNA-seq deconvolution of bulk samples into cell types
(Chapter 7), CRISPR-screen recovery of per-gene effects from pooled
phenotypes (Chapter 24), Hi-C contact recovery from sparse
observations (Chapter 10) — all are underdetermined linear systems
solved under sparsity priors. DIA is the proteomics instance of the
same family.

#note[
  The compressed-sensing framing also explains why DIA performance
  improves so dramatically with a good spectral library. The columns
  of $A$ are the theoretical templates; if the templates are wrong —
  missing peaks, mispredicted intensities — the recovery degrades.
  Neural-network spectral prediction (AlphaPeptDeep, Prosit) gives
  the library a much sharper template, which sharpens the sparse
  recovery downstream. The architecture follows the matched-filter
  rule of Chapter 16: get the template right and the noise model
  right, and the rest is bookkeeping.
]


== Identification, FDR, and Protein Inference <sec:ch27-identification>

Every peptide-spectrum match comes with a score. Distributing scores
across all candidate spectra produces two overlapping distributions —
correct matches at higher scores, incorrect matches at lower scores —
and the bioinformatician's job is to set a threshold that controls
the rate of false identifications without throwing away too many real
ones. The mechanism for doing this is *target-decoy estimation*, the
single most important statistical idea in modern proteomics.

=== The Target-Decoy Approach

The target-decoy approach, formalised by Elias and Gygi in 2007 but
foreshadowed by Moore, Young, and Lee in 2002, is conceptually simple.
You search not only against the real proteome (the *target* database)
but also against a *decoy* database constructed to be biologically
implausible — typically by reversing each protein sequence, sometimes
by shuffling residues with the K/R termini preserved. Any spectrum
that matches a decoy peptide is by construction a false match: the
decoy sequence does not exist in nature. If you count $D$ decoy hits
at some score threshold $T$ and $N$ target hits at the same threshold,
the *false discovery rate* (FDR) is estimated as $hat("FDR") = D / N$,
under the assumption that decoys mimic the distribution of incorrect
target matches.

The standard threshold across the field is *1% FDR at the peptide
level*. Every published proteomic identification list in the past
fifteen years has been filtered at this level or stricter; benchmark
challenges enforce it; reviewers demand it. The community converged on
1% because tighter thresholds lose too many real identifications and
looser ones admit too many false positives to be useful for downstream
analysis. There is nothing magic about 1%; the calibration is empirical.

#figure(
  image("../../diagrams/lecture-27/09-target-decoy.svg", width: 95%),
  caption: [
    Target and decoy score distributions overlap at low scores and
    separate at high scores. Sliding a threshold trades sensitivity
    for specificity; the 1% FDR cutoff is the convention.
  ],
) <fig:target-decoy>

#figure(
  image("../figures/ch27/f4-target-decoy-fdr.svg", width: 95%),
  caption: [
    Target and decoy score distributions on a simulated dataset, with
    the empirical FDR curve as a function of score threshold. The
    crossover between the two distributions is where the FDR is
    determined; moving the threshold above the crossover loses real
    identifications faster than it removes false ones.
  ],
) <fig:fdr-curve>

The original target-decoy estimator counts hits and divides; modern
implementations add a machine-learning step. *Percolator* (Käll, Canterbury,
Weston, Noble, and MacCoss, _Nature Methods_ 2007) trains a
semi-supervised classifier — originally an SVM, now usually a more
flexible model — on a set of features per PSM (score, mass accuracy,
precursor charge, number of matched ions, etc.), with targets as
positives and decoys as negatives. The classifier's posterior probability
becomes the new ranking score, and the FDR is re-estimated on the
ranked list. Percolator routinely recovers 30% to 50% more
identifications at the same nominal FDR than the underlying search
engine produces alone. Every modern pipeline either uses Percolator
directly or implements a variant.

#warn[
  Target-decoy FDR estimates assume the decoy distribution is a faithful
  mirror of the incorrect-target distribution. The assumption fails in
  predictable ways. *Concatenated* target-decoy searches (one combined
  database) can interact badly with multi-protein peptide ambiguity.
  *Reversal* preserves amino-acid composition but breaks all consecutive
  pairs that might bias scoring. *Shuffling* preserves composition but
  not amino-acid pairs at all. The community defaults are reasonable
  for typical datasets; pathological samples (heavily modified peptides,
  unusual proteomes, mass-tolerance mismatches) can produce FDR
  estimates off by factors of two or more. A well-engineered pipeline
  reports both peptide-level and protein-level FDR and treats the latter
  as the more conservative summary.
]

=== Protein Inference

Identifying peptides is half the job; aggregating them into proteins
is the other half. Most peptides come from one specific protein, but
some sequences are shared across multiple proteins in the proteome —
typical sources are paralogous gene families, splice isoforms of the
same gene, and short conserved domains. A shared peptide cannot, on
its own, tell you which member of the family it came from.

The standard solution is *parsimony*, sometimes called Occam's razor
for proteomics: the smallest set of proteins that can explain all
observed peptides is preferred. MaxQuant's *razor peptide* rule
assigns each shared peptide to the protein that has the largest
number of unique peptides. Other tools build a bipartite graph of
peptides and proteins and solve the protein-inference problem more
formally, but the practical result is the same: shared peptides
generally end up assigned to their most-supported parent protein.

Once proteins are assembled, a *protein-level FDR* is computed and
held below 1% in the same target-decoy framework. The protein-level
threshold is usually more stringent than the peptide-level one,
because a single false-positive peptide can flip an entire protein
group from absent to present. A typical published proteomics paper
reports both: 1% peptide FDR, ≤ 1% protein FDR, and a list of
proteins identified at that threshold.

#figure(
  image("../../diagrams/lecture-27/03-ms2-id.svg", width: 95%),
  caption: [
    Identification cascade: an MS2 spectrum with annotated b/y peaks
    feeds the search engine, which produces target and decoy score
    distributions, which yield an FDR-versus-threshold curve. The
    horizontal line at FDR = 0.01 picks the score threshold that
    keeps the field's standard 1% false-discovery rate.
  ],
) <fig:ms2-id>

=== Post-Translational Modifications

A protein is not just a sequence; it carries chemical decorations that
the genome does not encode. *Post-translational modifications* (PTMs)
add or remove functional groups at specific residues, with mass shifts
that are characteristic and detectable. The five most-studied PTMs in
proteomics, with their mass deltas, are:

- *Phosphorylation* on serine, threonine, or tyrosine: $+79.97$ Da
  (addition of $-"PO"_3"H"$). Adds a negative charge at physiological pH;
  regulates roughly thirty percent of the proteome.
- *Acetylation* on lysine or the N-terminus: $+42.01$ Da (addition of
  acetyl, $-"COCH"_3$). Neutralises positive charge; central to histone
  regulation and metabolic enzyme tuning.
- *Methylation* on lysine or arginine: $+14.02$ Da per methyl group.
  Permits "code" of mono-, di-, and tri-methylation states.
- *Ubiquitination* on lysine: leaves a glycine-glycine remnant after
  trypsin digestion, $+114.04$ Da on the modified lysine ("K-GG").
  Targets proteins for proteasomal degradation; also signals trafficking
  and DNA damage repair.
- *Oxidation* on methionine: $+15.99$ Da. Mostly an artifact of sample
  handling but biologically meaningful for some proteins.

#figure(
  image("../../diagrams/lecture-27/10-ptm-shifts.svg", width: 95%),
  caption: [
    Mass shifts for common PTMs, drawn on a peptide diagram. Each
    modification adds a fixed offset to the b/y fragment ions
    containing the modified residue; the offset is the same for every
    peptide carrying that PTM, which is what makes them detectable
    by an expanded search.
  ],
) <fig:ptm-shifts>

Detecting modified peptides is conceptually straightforward: add the
PTM masses to the search as *variable modifications*, and let the
search engine consider both the unmodified and modified versions of
each candidate. The price is *combinatorial explosion*. A peptide
with three potential phosphorylation sites has eight possible
modification states ($2^3$); enabling several variable PTMs at once
multiplies the candidate set by an additional factor per modification.
Practical tools cap the number of variable mods (typically four to
six) and the number of sites per peptide (typically three). MSFragger's
fragment-index data structure handles open searches more efficiently
than older engines, but the underlying combinatorial cost remains.

The deeper conceptual point is that *PTM searching is matched-filter
codebook expansion*. The unmodified peptide's template is the base
signal; each PTM adds a fixed offset to specific fragments, expanding
the dictionary of templates the engine has to consider. The expansion
is exactly the codebook-size-versus-detection-power tradeoff that
shows up in modulation theory: a larger codebook covers more possible
inputs but takes longer to search through and has higher per-input
false-detection rate. In proteomics, the practical answer is to
enrich for the modification before the LC-MS step, so that most of
the signal comes from modified peptides and the search engine does
not have to look for needles in haystacks.

=== Phosphoproteomics

The canonical PTM-enrichment workflow is *phosphoproteomics*. Tryptic
peptides are passed over an affinity resin — *titanium dioxide* (TiO₂)
beads or an *iron-IMAC* (immobilised metal affinity chromatography)
column — that selectively binds phosphate groups. Non-phosphorylated
peptides flow through; phosphorylated ones are retained, eluted at
high pH, and analysed by LC-MS/MS in the standard pipeline. A single
phosphoproteomics experiment routinely identifies twenty thousand
phosphosites across the proteome.

#figure(
  image("../../diagrams/lecture-27/12-phospho-enrich.svg", width: 95%),
  caption: [
    Phosphoproteomics enrichment. Tryptic peptides flow over a TiO₂
    or IMAC column; non-phosphorylated peptides elute first;
    phosphopeptides are retained and eluted at high pH. The
    enriched fraction is analysed by standard LC-MS/MS with phospho
    as a variable modification.
  ],
) <fig:phospho>

Phosphoproteomics intersects with kinase signalling networks
(Chapter 22) and with the systematic screens of Chapter 24. Identifying
the phosphosites is only the first step; deciding which kinase
phosphorylated which site, and inferring the upstream signalling
pathway, is the harder downstream problem. The 2020s have seen
multiple foundation-model-style efforts to predict kinase-substrate
relationships from sequence and phosphoproteomic data; the field's
ground truth remains *in vitro* kinase assays and CRISPR/RNAi
perturbations.


== Quantification <sec:ch27-quantification>

Identifying that a peptide is in a sample is one question; measuring
how _much_ of it is there is another. Quantification in proteomics
takes three main forms — *label-free*, *metabolic labelling*, and
*chemical labelling* — and each comes with a different
precision-versus-throughput trade-off.

=== Label-Free Quantification

The cheapest and most flexible approach is *label-free quantification*
(LFQ). Two samples to compare? Run them as two separate LC-MS
acquisitions and compare the resulting peptide intensities. The
canonical implementation is *MaxLFQ* (Cox et al., 2014): for each
peptide, integrate the area under its MS1 elution profile across the
LC gradient, normalise across samples using a robust median-of-ratios
estimator, and aggregate the per-peptide intensities into a per-protein
intensity using a maximum-likelihood model. The output is a
protein-by-sample matrix of log-scale intensities, ready for
differential abundance analysis. Open-source DIA tools (DIA-NN,
OpenSWATH) compute LFQ-style values from DIA data with the
deconvolution step replacing the precursor-selection step.

LFQ's strengths are flexibility (any sample, any cohort size, no
labelling step) and cost (just the LC-MS runs). Its weaknesses are
between-run variability (typically 10–15% coefficient of variation on
DDA, 5% on DIA), retention-time drift, and a tendency to leave missing
values when peptides fall below the MS1 detection threshold in some
samples.

=== Stable-Isotope Labelling: SILAC

The second approach is *metabolic labelling* — incorporating a
mass-distinguishable variant of an amino acid into the protein in vivo.
The canonical recipe is *Stable Isotope Labeling by Amino acids in
Cell culture* (SILAC), introduced by Shao-En Ong and Matthias Mann in
2002. Cells are grown in two media: one with normal ("light")
lysine and arginine, the other with heavy isotope-labelled
($""^13 "C"_6, ""^15 "N"_2$) lysine and ($""^13 "C"_6, ""^15 "N"_4$)
arginine.
After several cell-doubling generations the heavy amino acids are
fully incorporated into the proteome. The two cell populations are
mixed at 1:1, processed together through the rest of the workflow,
and the LC-MS/MS spectra show every tryptic peptide as a *pair* of
peaks — light and heavy — separated by a known mass difference. The
ratio of intensities within each pair is the per-peptide relative
abundance.

SILAC's strength is precision. Because the two samples are mixed
before any of the steps that introduce technical variability —
digestion, chromatography, ionisation, fragmentation — the heavy-to-light
ratio cancels most of those sources of noise. Typical CV is 5%, the
lowest of any quantification method short of targeted MRM. The
weakness is that it requires cell culture; clinical samples, primary
tissue, and most patient-derived material cannot be SILAC-labelled.
The *Super-SILAC* extension mixes a SILAC-labelled cell-line reference
into unlabelled clinical samples as a common reference, partially
recovering the benefits.

=== Chemical Labelling: TMT and iTRAQ

The third approach is *isobaric chemical labelling*, where every
sample is tagged with a chemical reagent that introduces a known
mass-shift pattern into MS2 fragmentation. *iTRAQ* (Thompson et al.,
_Analytical Chemistry_ 2003) and its successor *TMT* (Tandem Mass
Tags) attach a chemical label to every peptide's amine groups (N-terminus
and lysine side chains) at a fixed mass — the labels are *isobaric*,
meaning identical in total mass — but each label fragments in MS2 to
release a sample-specific *reporter ion* at a characteristic low
$m/z$. Six, ten, eleven, or sixteen samples can be labelled with
different tag variants, mixed, and run as a single LC-MS acquisition.
The MS2 spectra contain both peptide fragment ions (for
identification) and reporter ions (for quantification).

The reporter-ion intensities give the relative abundance of the
peptide across all multiplexed samples in one spectrum, which is
TMT's killer feature. Sixteen samples per run cuts instrument time
sixteen-fold relative to LFQ. The price is *ratio compression*: in
DDA mode, every co-isolated peptide in the precursor window contributes
its own reporter ions, so the measured ratios are diluted toward 1:1.
Modern MS3 acquisition (TMT-MS3, "SPS-MS3") mitigates the problem by
fragmenting MS2 fragment ions a second time before reading reporter
intensities; the instrument time cost is steep, but the quantitative
precision recovers.

#figure(
  image("../../diagrams/lecture-27/04-quantification.svg", width: 95%),
  caption: [
    Four quantification approaches at a glance. LFQ integrates the
    MS1 elution profile of each peptide; SILAC compares the
    intensities of heavy/light pairs; TMT reads reporter ions in MS2;
    PRM monitors specific precursor-to-fragment transitions. Each
    sits at a different point on the CV-versus-throughput plane.
  ],
) <fig:ch27-quant>

#figure(
  image("../figures/ch27/f5-quant-tradeoff.svg", width: 95%),
  caption: [
    Quantification methods on the CV-versus-multiplexing plane.
    SILAC has the lowest CV (mixing happens before any wet-lab
    variability) but is cell-culture only. TMT multiplexes up to 16
    samples per run but pays in ratio compression. LFQ is the
    most flexible and the noisiest. Targeted PRM is the most
    precise per peptide and the lowest throughput.
  ],
) <fig:quant-tradeoff>

The DIA quantification advantage is worth re-emphasising. Because
DIA's deterministic acquisition produces matched peptides across
samples by construction, DIA-LFQ delivers TMT-comparable CVs (around
5%) at lower per-run cost and without the ratio-compression problem.
For most cohort-scale or clinical proteomics studies starting in
2024, DIA-LFQ is now the default first choice.


== The Protein Matrix and Differential Abundance <sec:ch27-stats>

The output of any of the above quantification workflows is a
*protein-by-sample matrix*. Rows are proteins, columns are samples,
entries are log-scale relative abundances. A typical study produces
a matrix of about 5,000 to 8,000 proteins by 30 to 100 samples —
roughly the same shape as a bulk RNA-seq count matrix (Chapter 5),
which is no accident.

The statistical machinery is mostly borrowed from RNA-seq with a few
proteomics-specific adjustments. *limma* (Smyth, 2004), the linear
model with empirical-Bayes variance shrinkage that became the standard
for microarray and then for RNA-seq, transfers directly to proteomics
data once the intensities are log-transformed and normalised. Proteomics
adaptations like *DEqMS* (Zhu, Orre, Tran, Mermelekas, et al., 2020)
extend limma's shrinkage to account for the fact that protein-level
variance depends systematically on the number of peptides measured per
protein — proteins identified from many peptides have lower variance
than singletons. *ProDA* (Ahlmann-Eltze and Anders, 2020) replaces the
Gaussian noise model with a probabilistic dropout-aware model that
handles missing values without imputation.

=== The Missing-Value Problem

Proteomics data is much more missing than RNA-seq data. A typical
DDA-LFQ matrix has 30%–50% of entries missing; even DIA matrices have
5%–15% missing values, mostly clustered in the low-abundance proteins
where peptides occasionally fall below the detection threshold. The
mechanism behind a missing value matters: a value can be missing
because the protein is genuinely absent — *missing not at random*
(MNAR), with mass deltas systematically below the limit of detection —
or because the protein is present but happened not to be picked for
fragmentation in that sample — *missing at random* (MAR), independent
of the protein's true abundance.

#figure(
  image("../../diagrams/lecture-27/11-missing-values.svg", width: 95%),
  caption: [
    Missing-value patterns in a typical proteomics matrix. Low-abundance
    proteins show block-like MNAR patterns (whole columns missing for
    samples below the detection threshold); MAR misses are scattered.
    Imputation choices that assume one pattern when the other holds
    introduce systematic bias.
  ],
) <fig:missing>

The choice of imputation strategy depends on the assumed mechanism.
*Half-min imputation* — replacing missing values with half the
minimum observed value — assumes MNAR and is appropriate for
low-abundance proteins below detection. *K-nearest-neighbours
imputation* in protein-space assumes MAR and is appropriate when the
missing values are stochastic at random positions. *Multiple imputation*
under a probabilistic model averages across both possibilities.
*No imputation* — analysing the missing values explicitly under a
dropout model, as ProDA does, or by working in DIA matrices where
missingness is rare — is increasingly the preferred approach. The
single biggest mistake is to assume a uniform mechanism: if a fraction
of the missing values are MNAR and a fraction are MAR, neither
half-min nor KNN is correct, and the resulting differential-abundance
analysis can be biased in opposite directions for high-abundance
versus low-abundance proteins.

#warn[
  Several published proteomics analyses report differential proteins
  that turn out, on closer inspection, to be artifacts of the
  imputation step rather than real signal. Half-min imputation in
  particular can create apparent up-regulation in one group simply
  because the protein is below detection in the other; the
  half-minimum stand-in is then "increased" by however much real
  abundance the group with detected values shows. A standard
  sanity check before publishing: rerun the differential analysis
  with at least two imputation strategies (e.g., half-min and KNN)
  and only report the proteins that are differential under both.
]

=== Pathway Analysis

Once you have a list of differentially abundant proteins, the
downstream analysis tracks the RNA-seq pipeline of Chapter 22.
Over-representation analysis (ORA) tests whether known biological
pathways are enriched among the hits; gene-set enrichment analysis
(GSEA) generalises to the full ranked list; protein-protein interaction
networks from STRING propagate the signal across the interactome.
KEGG and Reactome are the standard pathway databases; Gene Ontology
the standard functional vocabulary. The proteomics-specific wrinkle is
that some PTM-aware analyses look not just at protein abundance but
at *per-modification* changes — phosphosite-level differential
analysis, for instance, where the question is which phosphosites
change occupancy rather than which proteins change total abundance.


== Applications <sec:ch27-applications>

The proteomics pipeline is general, but the questions it answers
are specific. Three application domains carry most of the weight.

=== Biomarker Discovery

The first is *biomarker discovery* — finding plasma or serum proteins
whose abundance changes with disease and can be measured in a clinical
assay. The plasma proteome spans roughly ten orders of magnitude in
abundance, from albumin at $10^(-4)$ molar down to amyloid-β at
$10^(-13)$ molar. No single assay covers the whole range. *Untargeted
LC-MS/MS* on plasma reaches roughly $10^(-9)$ molar; affinity-based
panels like *Olink* (proximity extension assay, about 5,000 proteins
per panel) and *SomaScan* (modified-aptamer pulldown, also about 5,000
proteins) reach lower; *Simoa* (single-molecule digital ELISA, Quanterix)
reaches single-molecule sensitivity for a few hundred proteins at a
time. Each technology occupies a band of the abundance landscape, and
study design starts with a choice of which band.

#figure(
  image("../../diagrams/lecture-27/05-biomarkers.svg", width: 95%),
  caption: [
    The plasma proteome spans ten orders of magnitude in abundance.
    Each assay technology covers a different band; LC-MS/MS is
    quantitative and unbiased but bottoms out around $10^(-9)$ M;
    affinity arrays extend lower at the cost of binder availability;
    Simoa reaches single molecules for a few hundred targets. Most
    clinically useful biomarkers sit between $10^(-9)$ and $10^(-12)$ M.
  ],
) <fig:ch27-biomarkers>

The list of clinically used protein biomarkers is short — perhaps fifty
proteins routinely measured in clinical chemistry — but each one is
load-bearing. Cardiac *troponin I and T* are the diagnostic gold
standard for myocardial infarction; *NT-proBNP* and *BNP* diagnose
heart failure; *apolipoprotein B* is a more discriminating
cardiovascular-risk marker than total cholesterol or LDL-C and is
the target of statin and PCSK9-inhibitor therapy; *PSA* screens for
prostate cancer; *CA-125, CA 19-9, AFP* are tumour markers in
established oncology workflows. The
2024 landscape is shifting under two pressures. Olink-based discovery
studies in biobank-scale cohorts (UK Biobank, FinnGen) are publishing
hundreds of plausible new biomarker candidates per year; the FDA's
biomarker qualification process is slow and few candidates ever make
it into clinical use. Proteomics has produced many discovery hits;
turning them into clinical-grade assays remains the bottleneck.

=== Chemoproteomics and Target Validation

The second application is *chemoproteomics* — identifying the proteins
that a small-molecule drug binds in cells. Three canonical workflows
dominate. *ABPP* (activity-based protein profiling, introduced by
Cravatt and colleagues in 1999) labels active enzymes via a chemical
warhead that covalently modifies the active site, then enriches and
identifies the labelled proteins by MS. *TPP* (thermal proteome
profiling, introduced by Savitski, Reinhard, Franken, et al., _Science_
2014) exploits the fact that drug-bound proteins are more thermally
stable than unbound ones; cells treated with a drug are heated through
a gradient, the soluble fraction is collected at each temperature, and
proteins that resist denaturation at higher temperatures when the drug
is present are inferred to be drug-bound. *KiNativ* and similar
kinase-focused platforms use ATP-competitive probes to identify the
kinases a kinase inhibitor engages in cells. The output of all three
is a list of *targets* and *off-targets*: the proteins the drug binds
intentionally and the proteins it binds unintentionally, both of
which determine clinical activity and toxicity.

#figure(
  image("../../diagrams/lecture-27/06-chemoproteomics.svg", width: 95%),
  caption: [
    Chemoproteomics in two paragraphs. A drug coupled to a chemical
    probe (or labelled by a covalent warhead) marks the proteins it
    binds in live cells; the marked proteins are pulled down or
    enriched and identified by LC-MS/MS. The output is a target list
    and, more importantly, an off-target list. The kinase-inhibitor
    example identifies roughly five expected targets and ten to
    fifteen off-targets in a typical experiment.
  ],
) <fig:chemoproteomics>

Chemoproteomics is the experimental complement to the AlphaFold-based
docking workflows of Chapter 15. A predicted binding pose is a
hypothesis; a chemoproteomic pulldown is the live-cell measurement
that can falsify it. For drug discovery in 2024 the two are
increasingly combined — predict candidates computationally, confirm
binding in cells, iterate.

=== Structural Proteomics and Single-Cell Frontier

A handful of MS techniques look at protein structure rather than just
abundance. *Cross-linking mass spectrometry* (XL-MS) uses a chemical
crosslinker to lock interacting residues together; MS then identifies
the crosslinked peptide pairs, which constrain the spatial layout of
protein complexes. *Hydrogen-deuterium exchange* (HDX-MS) measures
how readily backbone amide hydrogens exchange with deuterated solvent;
exchange rates report on local flexibility and ligand-binding-induced
conformational changes. *Native MS* keeps protein complexes intact
through the ESI source and measures their assembled masses directly.
All three complement AlphaFold and cryo-EM by giving conformational
or interaction-level constraints in solution rather than in a crystal.

*Single-cell proteomics* is the field's current frontier. SCOPE-MS
and nano-PiMMS workflows (Slavov lab, ~2017 onwards) push LC-MS/MS to
single-cell sensitivity by miniaturising the sample prep and using
TMT to multiplex many cells per run. As of 2024, a typical single-cell
proteomics experiment quantifies 1,000–3,000 proteins per cell with
substantial cell-to-cell missingness — still a gap of an order of
magnitude versus single-cell RNA-seq, but closing. Spatial proteomics
techniques like *Imaging Mass Cytometry* (IMC) and *CODEX* measure
40–60 protein markers per cell at tissue-section resolution, sitting
between the genomic-scale of MS and the targeted nature of antibody
panels.


== Metabolomics: Same Instrument, Sparser Reference <sec:ch27-metabolomics>

The metabolome — the small-molecule output of cellular metabolism —
sits one step further downstream from the genome than the proteome.
Where proteomics measures 20,000 protein-coding gene products,
metabolomics measures the 1,000 to 5,000 small molecules detectable in
human plasma: amino acids, organic acids, nucleotides, lipids, sugars,
hormones, drug metabolites, and microbiome-derived compounds. The
metabolome integrates upstream genetic, transcriptional, translational,
and enzymatic activity into a single readout; it is often the closest
measurement to clinical phenotype.

The instrument is the same LC-MS platform, but the chromatography and
ionisation choices are different. Polar compounds (sugars, amino acids,
nucleotides) need *hydrophilic interaction* (HILIC) chromatography
rather than reverse-phase; lipids and other non-polar compounds use
reverse-phase but with different mobile phases. Most metabolomic
runs acquire in both *positive* and *negative* ESI modes — basic
compounds ionise as cations, acidic compounds as anions — so a full
workflow has up to four LC-MS runs per sample.

=== The Identification Gap

The bioinformatics is where metabolomics and proteomics diverge most
sharply. A peptide has a sequence; a tryptic peptide can be enumerated
from the genome. A metabolite has a molecular formula and a chemical
structure but no sequence; the space of "all possible metabolites" is
vast and only partially catalogued.

The reference databases are the bottleneck. The *Human Metabolome
Database* (HMDB, Wishart and colleagues, since 2009) lists about
250,000 putative human metabolites, but only a fraction have
experimentally measured tandem-MS spectra. *METLIN* covers tandem-MS
spectra for about a million molecules, commercially. *MoNA* (the
MassBank of North America) and *GNPS* (Global Natural Products Social)
provide community-shared spectral collections. Even with all four,
*30%–50% of detected peaks in a typical untargeted metabolomic run
remain unidentified.* The unidentified peaks include real metabolites
whose spectra are not in any database, drug metabolites the human
liver has produced, microbiome-derived compounds, and exotic
fragmentation artifacts.

The *Metabolomics Standards Initiative* (MSI) confidence levels
formalise the spectrum of identification quality. *Level 1* is a
match against an authentic standard analysed on the same instrument —
mass, retention time, and MS/MS all agree. *Level 2* is a putative
annotation by spectral library match without the standard. *Level 3*
is a putative compound class (e.g., "a phosphatidylcholine") without
species-level resolution. *Level 4* is an unknown feature with mass
and retention time but no annotation. Most peaks in an untargeted run
are Level 3 or 4; the burden of upgrading a Level-3 peak to Level 1
is what keeps metabolomic biomarker validation slow.

#figure(
  image("../../diagrams/lecture-27/13-metabolomics.svg", width: 95%),
  caption: [
    A metabolomics LC-MS workflow. Polar metabolites are separated by
    HILIC; non-polar by reverse-phase; both are acquired in positive
    and negative ESI modes. Identification cascades through accurate
    mass, isotope pattern, retention time, and MS/MS spectral
    matching, with the unfortunate 30%–50% residue that remains
    unidentified.
  ],
) <fig:metabolomics>

A handful of computational tools narrow the identification gap.
*XCMS* (Smith, Want, O'Maille, Abagyan, and Siuzdak, 2006) is the
canonical open-source platform for untargeted metabolomic peak
detection and alignment. *MetaboAnalyst* (Pang, Chong, Li, and Xia,
2024) wraps statistical and pathway-enrichment analysis. *SIRIUS*
(Dührkop and colleagues, since 2015) uses fragmentation-tree
analysis plus a machine-learning module (CSI:FingerID) to predict
chemical structure from MS/MS spectra of unknown peaks — when it
works, it lifts a Level-4 peak to Level 2 with reasonable confidence.
The 2020s have seen graph-neural-network and transformer-based
extensions of this idea; the field has not yet converged on a single
winner.

#note[
  The "identification gap" in metabolomics is itself a matched-filter
  problem with a degraded codebook. In proteomics, the codebook is
  the tryptic-peptide enumeration of a known proteome; in
  metabolomics, the codebook is a sparse, vendor-fragmented,
  partially-empirical collection of chemical structures with their
  measured tandem-MS spectra. The gap will close as the codebook
  fills out — either by more experimental measurements of
  authentic standards, or by neural networks that predict the
  spectra of structures the databases do not yet contain. The Open
  Targets, MoNA, and GNPS communities are racing the long tail; the
  finish line is not in sight.
]

=== Targeted Metabolomics

For clinical applications, the alternative is *targeted metabolomics*:
a defined panel of 50 to 200 metabolites measured with isotope-labelled
internal standards and absolute calibration. Commercial platforms like
Biocrates AbsoluteIDQ and Metabolon's targeted panels deliver
clinical-grade quantification of a curated metabolite list. Targeted
metabolomics is what most published clinical-metabolomics association
studies actually use; untargeted runs are for discovery.

=== Lipidomics

A specialty branch of metabolomics focuses on *lipidomics* — the
~10,000 lipid species across phospholipids, sphingolipids, sterols,
fatty acids, and glycerolipids. Lipids are organised by class (a head
group) and species (variable acyl-chain lengths and saturations);
identification has a more structured grammar than the rest of the
metabolome. LipidSearch, LipidBlast, and LipidMatch are the standard
tools. Clinically relevant lipidomic markers include LDL particle
subclasses, ceramide profiles in diabetes and heart failure, and
drug-induced free-fatty-acid shifts.


== Tools, Pitfalls, and the 2024 Frontier <sec:ch27-tools>

A modern proteomics analysis in 2024 typically runs on one of three
software stacks. *MaxQuant* plus *Perseus* for DDA-LFQ is the most
mature pipeline and the default for many academic labs. *FragPipe*
(MSFragger + IonQuant + Philosopher) is faster and supports open
searches for PTM-heavy applications. *DIA-NN* or *Spectronaut* for
DIA workflows. Downstream statistical analysis runs in R with limma,
DEqMS, MSstats, or proteomics-specific extensions; pathway analysis
runs in R or Python through standard ORA and GSEA tools; visualisation
through Perseus, ggplot2, or matplotlib.

#figure(
  image("../../diagrams/lecture-27/07-workflow.svg", width: 95%),
  caption: [
    The end-to-end proteomics workflow, with reproducibility-critical
    decisions called out at each stage. Sample collection, digestion
    quality, search-parameter choices, FDR threshold, quantification
    method, normalisation, missing-value handling, and statistical
    model all leave fingerprints in the final hit list.
  ],
) <fig:ch27-workflow>

=== Public Repositories

The infrastructural backbone of the field is *ProteomeXchange* — a
consortium of repositories including *PRIDE* (Vizcaíno and colleagues
at EMBL-EBI, since 2003), *MassIVE* at UCSD, and *jPOST* in Japan.
Every published proteomics dataset is supposed to be deposited in
one of these archives with associated metadata and raw files. The
practice is uneven — older datasets are sparser than newer ones,
metadata quality varies — but ProteomeXchange contains hundreds of
thousands of datasets that constitute a public training corpus for
the next generation of proteomics ML. *MetaboLights* (EMBL-EBI) and
*Metabolomics Workbench* play the same role for metabolomics, less
comprehensively.

=== Common Pitfalls

The recurring failure modes of a proteomics pipeline are predictable
enough to list. *Sample preparation* introduces *keratin contamination*
(skin cells from sample handling) that produces a characteristic set
of false-positive peptide IDs and that contaminates downstream
analyses; the standard fix is to filter keratin peptides explicitly
after the search. *Incomplete tryptic digestion* leaves missed
cleavages and shifts the peptide population; allowing up to two
missed cleavages in the search recovers most of the affected peptides
but the search space grows.

*Acquisition* artifacts include mass-accuracy drift (recalibrate
periodically; tools detect and correct), carryover between samples
(blank runs between samples or sample randomisation), and
DDA stochastic gaps (move to DIA where reproducibility matters).

*Identification* errors stem from wrong PTM searches (every additional
variable mod expands the false-positive count proportionally; turn
on only those you have biological reason to expect), database mismatch
(make sure the search database matches the sample's species and
reasonable splice isoforms), and multi-protein peptide ambiguity (rely
on the search engine's parsimony rule, and look at protein groups
rather than single proteins in the output).

*Quantification* requires explicit attention. Missing values need a
defensible imputation strategy, or none at all (ProDA-style dropout
modelling). TMT runs need ratio compression awareness — interpret
fold changes conservatively, especially for highly-multiplexed runs.
Batch effects are real even within a single experiment — randomise
sample order within and across runs, model batch as a covariate in
the differential analysis.

#tip[
  The fastest sanity check on a new proteomics dataset is to plot
  per-sample protein-ID count, per-sample missing-value count, and
  inter-sample peptide-ID overlap before any analysis. Out-of-the-ordinary
  samples — three standard deviations below the mean ID count, a
  systematic missing-value pattern, an outlier in PCA — usually
  reflect technical problems (incomplete digestion, LC failure,
  detector miscalibration) that need to be resolved before
  differential analysis is meaningful.
]

=== Reproducibility

Proteomics reproducibility is intermediate between RNA-seq and
microbiome studies. Cohort-scale DIA workflows now achieve 80%–85%
protein-ID overlap across technical replicates and 60%–70% replication
of differential-abundance hits across independent cohorts.
DDA-LFQ studies sit lower, around 70% identification overlap and
50%–60% differential replication. Microbiome studies sit lower still;
RNA-seq, when the same expression-quantification pipeline is used, sits
higher, around 85% replication. The proteomics improvement curve over
the past five years — driven primarily by the DIA transition — has
been the largest among the omics platforms.

=== The 2024 Frontier

Three threads define the 2024 state of the art. *DIA with neural
spectral libraries* — DIA-NN's library-free mode and AlphaPeptDeep's
spectrum prediction — has all but eliminated the historical bottleneck
of needing a pre-acquired DDA library. *Single-cell proteomics* is
crossing from technical demonstration to early biological application;
the SCOPE2 and related workflows now quantify roughly 1,500 proteins
per cell with reasonable CVs. *Foundation models for proteomics* —
trained on the hundreds of thousands of MS2 spectra in PRIDE and
MassIVE — promise to predict spectra, retention times, and PTM
fragmentation across the entire proteome without species-specific
training. AlphaPeptDeep (Mann group, 2022) and the Casanovo *de novo*
sequencing model (Yilmaz and colleagues, 2022) are early examples of
the genre. None has yet displaced the classical search-engine pipeline
for routine identification, but the displacement is plausible within
the next instrument generation.


== Summary <sec:ch27-summary>

- *Mass-spectrometry proteomics measures the protein-level state of
  a cell or sample*, complementing RNA-seq's transcript-level
  measurement. mRNA–protein abundance correlations are around
  $r = 0.5$, low enough that the proteomic readout is not redundant
  with the transcriptomic one.
- *LC-MS/MS* in bottom-up shotgun mode is the dominant workflow.
  Trypsin digests proteins into peptides; reverse-phase LC separates
  them; ESI ionises them; MS1 measures intact masses; MS2 measures
  backbone-fragmentation b/y patterns; a search engine matches
  observed spectra against theoretical templates from a reference
  proteome.
- *Peptide identification is constrained matched filtering*. The b/y
  template is the matched filter; the tryptic constraint is the search
  space; cross-correlation or probability scoring is the decision rule.
  SEQUEST, Mascot, MaxQuant/Andromeda, MSFragger, and Percolator are
  the lineage.
- *Target-decoy estimation* is the field's universal FDR mechanism.
  A reversed or shuffled decoy database produces a controlled
  false-positive rate; the community standard is 1% peptide-level FDR.
- *DDA vs DIA* is the central acquisition trade-off. DDA picks
  top-$N$ precursors and produces clean but stochastic identification;
  DIA scans systematic $m/z$ windows and produces deterministic
  identification at the cost of chimeric spectra. Modern DIA tools
  (DIA-NN, Spectronaut) solve the deconvolution problem as a
  compressed-sensing recovery on a peptide-coverage matrix.
- *Quantification methods* sit on a CV-versus-multiplexing trade-off
  plane. SILAC (5% CV, cell-culture only), TMT (multiplexed, ratio
  compression), LFQ (cheap, 10–15% CV on DDA, 5% on DIA), and PRM
  (absolute, low-throughput) each occupy a different operating point.
- *Post-translational modifications* — phosphorylation, acetylation,
  methylation, ubiquitination — add fixed mass shifts that the search
  engine can find via expanded codebook searches. Phosphoproteomics,
  via TiO₂ or IMAC enrichment, routinely identifies twenty thousand
  phosphosites per experiment.
- *Missing values* are the proteomics-specific statistical wrinkle.
  MNAR (truly absent) and MAR (stochastically missed) need different
  treatments; imputation strategies and DIA workflows that minimise
  missingness are the practical responses.
- *Applications* include biomarker discovery (plasma proteome,
  Olink/SomaScan/Simoa/PRM), chemoproteomics (ABPP, TPP, KiNativ),
  structural proteomics (XL-MS, HDX-MS, native MS), and the
  early-frontier single-cell proteomics.
- *Metabolomics* uses the same LC-MS instrument with HILIC chromatography
  and dual-mode ESI to measure 1,000–5,000 small molecules; the
  identification gap (30%–50% of peaks unidentified in untargeted
  runs) is the field's structural limitation, addressed by HMDB,
  METLIN, MoNA, GNPS, and ML-driven structure prediction (SIRIUS,
  CSI:FingerID).


== A Closing Note: The Multi-Omics Future <sec:ch27-closing>

This is the last chapter of the book, and it is worth pausing on what
the twenty-seven-chapter arc adds up to. The course began with a
chemical polymer and an instrument that turns it into a text file
(Chapter 1). It walked the path from raw reads to genomes to variants
(Chapters 2–4), from RNA counts to expression dynamics (Chapters 5–6),
into the single-cell era (Chapters 7–8) and the regulatory and
epigenomic layers (Chapters 9–10), out to long reads and pangenomes
(Chapter 11), into population genetics and GWAS (Chapters 12–13),
through data engineering (Chapter 14) and protein structure (Chapter 15),
into the machine-learning architectures that have reshaped every prior
chapter (Chapter 16), and finally — across the latter half — into
clinical applications, regulatory frameworks, networks, causal
inference, ML in biology, and the integrative methods that pull the
pieces together. The thread that runs through all of it is the
recurring move every working bioinformatician makes: pick a biological
question, identify the measurement that has structural information
about that question, find the algorithm whose inductive bias matches
the measurement's shape, and validate against the messy reality of
biological variability. This chapter on mass spectrometry is the
last instrument we cover, but the move is the same one Chapter 1
introduced and every chapter in between has rehearsed.

The natural next step — beyond the scope of this book, but the place
the field is heading — is *multi-omics integration*. No single
measurement layer captures a biological system. RNA-seq (Chapter 5)
tells you what a cell intends to make; the proteome of this chapter
tells you what it actually made; the metabolome tells you what those
proteins did; single-cell methods (Chapter 7) reveal which cells are
which; ChIP-seq and ATAC-seq (Chapter 9) localise the regulators;
the GWAS catalogue (Chapter 13) anchors variants to phenotypes;
network methods (Chapter 22) interpret the layers as a system; causal
inference (Chapter 25) tries to ascribe direction to the resulting
edges. A serious clinical or basic-research question almost always
needs at least two of these layers, and the next decade of
bioinformatics is, in large part, about doing the integration well —
matched multi-modal foundation models (Chapter 16's frontier), causal
mediation analyses that span omics layers (Chapter 25), and clinical
biomarker panels that combine genetic risk, protein abundance, and
metabolite concentrations into a single risk model. The technical
pieces are largely available; what remains is the patient, careful
work of putting them together for specific biological systems.

Two trajectories are worth naming. The first is that *biology is becoming
a data problem in the way physics became a data problem in the
twentieth century*: not because the underlying systems are simpler than
chemistry but because the measurements have caught up with the
complexity. Whole-genome sequencing is now under \$200; deep proteomics on a
hundred-sample cohort is feasible in an academic lab; single-cell
multi-omics datasets that would have been a doctoral thesis in 2015
now ship as supplementary data. The bottleneck is no longer the
measurement; it is the modelling. Engineers — particularly EE-trained
engineers who think naturally in terms of signals, systems, models, and
control — are well placed to do that modelling. The fact that you have
read this far is itself evidence that the field is yours to shape.

The second trajectory is harder to predict but worth flagging. The
*foundation-model era of biology* — DNA language models, protein
language models, cell foundation models, multi-omics models — is two
or three years old as this book goes to press, and its endpoint is
uncertain. The DNA-language-model promise of "BERT for genomics" has
delivered useful tools but not, yet, a paradigm shift. AlphaFold, by
contrast, was a paradigm shift, and the next AlphaFold-shaped
breakthrough — perhaps in cell-state prediction, perhaps in causal
inference from observational omics, perhaps in something we have not
yet named — is the kind of event the field's history of the past
two decades suggests is increasingly likely. Whatever shape it takes,
the engineering discipline of the previous twenty-seven chapters will
remain the right framing: identify the structural property of the data,
match it to an architecture's inductive bias, validate against held-out
data, and never trust an evaluation protocol you have not personally
torn apart. The instrument changes; the engineering does not.


== Exercises <sec:ch27-exercises>

#strong[1.] #emph[Fragment-ion masses.] For the tryptic peptide
$sans("LSDPYHR")$, compute the theoretical b-ion and y-ion masses at
charge state +1 assuming monoisotopic amino-acid residue masses
(L = 113.084, S = 87.032, D = 115.027, P = 97.053, Y = 163.063,
H = 137.059, R = 156.101). Verify that consecutive y-ion mass
differences equal the residue mass of one amino acid each, then
identify which residue is read off at each step.

#strong[2.] #emph[Target-decoy FDR.] You have a peptide-spectrum-match
search with 50,000 target PSMs and 5,000 decoy PSMs. Compute the FDR
at the score threshold above which 8,000 targets and 80 decoys remain.
At what FDR would you set the threshold to recover 10,000 targets? Be
explicit about the assumption you are making about the decoy
distribution.

#strong[3.] #emph[Quantification choice.] You are designing a study
of plasma proteomic response to a drug across 60 patients (drug arm
and placebo arm, three time points each). Discuss the trade-offs of
DDA-LFQ, DIA-LFQ, TMT-16, and PRM for this study. Which method would
you pick, and why? Spell out the assumptions about CV, multiplexing,
batch effects, and cost.

#strong[4.] #emph[Missing-value imputation.] Take a small public
proteomics dataset from PRIDE (any DDA-LFQ MaxQuant output). Run
differential abundance under three imputation strategies: half-min,
KNN with $k = 5$, and no imputation (drop rows with > 50% missing).
Compare the overlap of significant proteins at FDR < 0.05 across the
three strategies. Comment on which strategy you would trust if the
missingness is genuinely MNAR for low-abundance proteins.

#strong[5.] #emph[PTM search expansion.] Estimate the multiplicative
factor by which the search space grows for a peptide of length 15 with
five potential phosphorylation sites and two potential acetylation
sites, all enabled as variable modifications, given a per-PTM
multiplicity cap of three sites per peptide. Discuss why open searches
(MSFragger) are computationally feasible despite this explosion.

#strong[6.] #emph[Metabolomics identification cascade.] You detect a
plasma feature at $m/z = 184.0734$ at retention time 6.2 minutes in
HILIC-positive mode. Walk through the identification cascade:
(a) what molecular formulas are consistent with this $m/z$ at 5 ppm
mass accuracy? (b) Pick one candidate (e.g., phosphocholine,
C5H15NO4P, $m/z = 184.0733$) and describe what additional data you
would need to confirm a Level-1 identification. (c) What is the
probability your detected feature is one of the 30%–50% that
remain unidentified, and what would you do about it?

#strong[7.] #emph[DIA as compressed sensing.] Sketch the matrix
algebra of DIA peptide identification. Define the candidate peptide
template matrix $A$, the observed-spectrum vector $bold(y)$, and the
peptide abundance vector $bold(x)$. Write down the regularised
optimisation problem the search engine is solving. Explain why
$ell_1$-style regularisation is appropriate and what the sparsity
assumption corresponds to biologically.

#strong[8.] #emph[(Open-ended.)] Pick a clinical biomarker — cardiac
troponin I, NT-proBNP, apolipoprotein B, prostate-specific antigen,
or another of your choice — and trace its history from discovery to
clinical implementation. What measurement technique was used in the
discovery study? What technique is used in the clinical assay today?
What proteomics or affinity-based panel could plausibly replace the
current assay, and at what cost-versus-precision trade-off?


== Further Reading <sec:ch27-further-reading>

- *Aebersold, R., and Mann, M.* (2003). "Mass spectrometry-based
  proteomics." _Nature_ 422: 198–207. The canonical review of the
  field as it crystallised; still the right starting point twenty
  years later.
- *Eng, J. K., McCormack, A. L., and Yates, J. R.* (1994). "An
  approach to correlate tandem mass spectral data of peptides with
  amino acid sequences in a protein database." _Journal of the
  American Society for Mass Spectrometry_ 5: 976–989. The SEQUEST
  paper; the algorithm that defined database search in proteomics.
- *Cox, J., and Mann, M.* (2008). "MaxQuant enables high
  peptide identification rates, individualised p.p.b.-range mass
  accuracies and proteome-wide protein quantification." _Nature
  Biotechnology_ 26: 1367–1372. The MaxQuant paper; pair with
  Cox et al. (2014) on MaxLFQ for the quantification side.
- *Elias, J. E., and Gygi, S. P.* (2007). "Target-decoy search
  strategy for increased confidence in large-scale protein
  identifications by mass spectrometry." _Nature Methods_ 4:
  207–214. The target-decoy formalisation.
- *Demichev, V., Messner, C. B., Müller, S. I.,
  Tolić-Nørrelykke, S. F., and Ralser, M.* (2020). "DIA-NN: neural
  networks and interference correction enable deep proteome coverage
  in high throughput." _Nature Methods_ 17: 41–44. The DIA-NN
  introduction; pair with Aebersold's SWATH-MS (2012) for the DIA
  lineage.
- *Savitski, M. M., Reinhard, F. B. M., Franken, H., Werner, T., et al.*
  (2014). "Tracking cancer drugs in living cells by thermal profiling
  of the proteome." _Science_ 346: 1255784. The TPP paper; the
  cleanest demonstration of in-cell drug-target identification by
  proteomics.
- *Dührkop, K., Shen, H., Meusel, M., Rousu, J., and Böcker, S.*
  (2015). "Searching molecular structure databases with tandem mass
  spectra using CSI:FingerID." _PNAS_ 112: 12580–12585. The
  computational side of the metabolomics identification gap; the
  current frontier of structure prediction from MS/MS.
- *Wishart, D. S., et al.* (2022). "HMDB 5.0: the Human Metabolome
  Database for 2022." _Nucleic Acids Research_ 50: D622–D631. The
  reference resource that anchors metabolomic identification.
- *Slavov, N.* (2020). "Single-cell protein analysis by mass
  spectrometry." _Current Opinion in Chemical Biology_ 60: 1–9. The
  current state and the near-term prospects of single-cell
  proteomics.
