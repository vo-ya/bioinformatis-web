#import "../theme/book-theme.typ": *

= Drug Discovery and Chemoinformatics: From SMILES to FDA <ch:drug-discovery>

#matters[
  A small-molecule drug is a chemical decision made under three
  simultaneous constraints. It must bind a particular protein hard
  enough to bend its behaviour. It must reach that protein, in a sick
  human, at a concentration that does the bending without poisoning
  the rest of the cell. And it must come out of a search through
  chemical space — a discrete set of around $10^60$ drug-like
  molecules — at a cost the developer can afford and on a timeline a
  patent allows. Bioinformatics meets that decision four times. It
  turns a target into a structure. It turns a library of compounds
  into ranked predictions. It turns assay-positive hits into leads.
  And, increasingly, it turns a pocket into a designed ligand without
  passing through screening at all. The methods are the chapter; the
  reason the chapter exists is that drug discovery is, on the
  computational side, the closest engineering analogue you will meet
  in this book to the constrained, multi-objective optimisation
  problems already familiar from VLSI placement and signal-processing
  filter design — with the inconvenient property that the cost
  function only fully reveals itself on humans, ten years and a
  billion dollars later.
]

A drug is a chemical entity that binds a specific protein, modulates
its activity in a clinically useful way, and has an acceptable
combination of safety, pharmacokinetics, and route of administration.
Most of what is on a pharmacy shelf is *small molecules*: organic
compounds below about 900 daltons, made by chemical synthesis,
swallowed as a pill. A growing fraction is *biologics* — antibodies,
peptides, oligonucleotides — that are larger, made in cells, and
delivered by injection. This chapter is about the small-molecule side.
The methods divide cleanly into five moves: turn a molecule into a
representation a computer can handle, search chemical space by
similarity, dock candidates into a target, filter for druglikeness,
and, when none of those quite work, design new molecules from the
pocket up. The last move is the one that has changed most since 2018
and the one that justifies the chapter sitting in the second half of
this book.

A reader who has worked through Chapters 15 and 16 will recognise the
shape of the problem. The molecule is a graph, the pocket is a 3-D
surface, the assay is a noisy observation of a free-energy difference,
and the architecture that handles each well is the one whose
inductive bias matches the geometric arrangement of the data —
permutation equivariance for the graph, SE(3) equivariance for the
pocket, an over-dispersed likelihood for the assay. Drug discovery is
where those ML pieces converge into a single multi-stage pipeline,
and where the failure modes from Chapter 16 — wrong representation,
wrong loss, leaky split, weak baseline — show up under different
names but with the same root cause.


== The Discovery Pipeline <sec:pipeline>

The path from "we want a drug against target $X$" to "an FDA-approved
pill" runs, in industry-standard form, through seven stages. Target
identification, where genetic, functional, and clinical evidence
converge on a protein worth modulating. Hit discovery, where some
combination of high-throughput screening, virtual screening, and
structure-based design produces hundreds to thousands of compounds
with sub-micromolar activity. Hit-to-lead, where those hits are
triaged for druggability, selectivity, and early ADMET signal.
Lead optimisation, where medicinal chemists work analogue series until
a small number of compounds clear the threshold for animal studies.
Pre-clinical pharmacokinetics and toxicology in rodents and a non-human
primate. Phase 1 safety in healthy volunteers. Phase 2 efficacy in
patients. Phase 3 confirmatory trials. And, if all of that holds, a
new-drug application to the FDA. The whole process takes ten to
fifteen years and costs, in 2024 dollars, one to three billion per
approved drug. The attrition is concentrated at the end: roughly nine
out of ten compounds that enter Phase 2 do not reach approval.

#figure(
  image("../figures/ch26/f1-attrition-funnel.svg", width: 95%),
  caption: [
    The drug-discovery attrition funnel. Roughly seven orders of
    magnitude separate the top of the virtual library from the single
    approved drug. Computational methods cheapen the upper bands; the
    biology of late-stage clinical attrition is unmoved by them.
  ],
) <fig:funnel>

@fig:funnel makes the asymmetry explicit. The funnel is steep but
roughly linear at the top, where computation is increasingly the
cheap option — a ten-million-compound virtual screen with Vina plus
a GNN rescorer now costs a few hundred thousand dollars of GPU time,
not the tens of millions a wet-lab high-throughput screen used to
require. It narrows abruptly at the IND-filing line, below which the
work is human biology and the costs are in clinical operations rather
than chemistry. Computation has done a lot for the top of the funnel
and almost nothing measurable for the bottom. The reason most "AI for
drug discovery" presentations focus on the top is that the bottom is
where ML's leverage is genuinely small.

Where genomics fits into the pipeline is everywhere except the
chemistry. *Target identification* uses GWAS to find disease-associated
loci, CRISPR screens to find essentialities, DepMap to find cancer
dependencies, and Mendelian randomisation to test whether a
biomarker is on the causal path or downstream of it. *Mechanism* uses
pathway analysis, structural biology, and increasingly AlphaFold-2
to characterise the target's biophysics. *Patient stratification*
uses clinical genomics to select trial populations where the drug is
most likely to work. *Biomarker development* uses bulk RNA-seq,
single-cell, and cancer genomics to find correlates of response. This
chapter walks the chemistry side — what the computational pipeline
does once a druggable target is in hand — but it would be misleading
to treat the chemistry in isolation. The genomic decisions taken
upstream determine whether the chemistry has any chance of producing
a drug at all.

#note[
  Druggability is the upstream question that constrains everything
  downstream. Of roughly twenty thousand human proteins, about three
  thousand are considered classically druggable — they are enzymes
  (especially kinases), G-protein-coupled receptors, ion channels,
  and nuclear receptors with well-formed binding pockets. The
  remainder — transcription factors, structural scaffolds,
  intrinsically disordered proteins — were until recently called
  "undruggable" for the simple reason that they presented no obvious
  small-molecule binding surface. Three things have eroded that
  category: targeted protein degraders such as PROTACs and molecular
  glues, which hijack the cell's own ubiquitin machinery rather than
  needing a deep pocket; AlphaFold-driven cryptic-pocket discovery,
  which finds transient binding sites that crystallography missed;
  and structure-based design of pocket-mimic peptidomimetics. The
  3000 / 20000 ratio has not yet doubled, but the upper bound on
  what is reachable is moving.
]


== Representations of Molecules <sec:representations>

Computers need a representation of a molecule, and there are three
sensible choices for one. *SMILES* is a string — the storage and
transport format. The *molecular graph* is the working representation
— what cheminformatics libraries operate on internally and what graph
neural networks consume directly. The *fingerprint* is a fixed-length
bit vector — the feature representation that classical ML methods
expect and that supports fast pairwise similarity. The three are not
alternatives; they are stages along a conversion chain, and each
stage discards or invents information in a way that matters for some
tasks and not for others.

#figure(
  image("../../diagrams/lecture-26/02-representations.svg", width: 95%),
  caption: [
    Three representations of aspirin. The SMILES string is
    canonicalised, transported, and looked up against databases. The
    molecular graph is what RDKit, OpenBabel, and any GNN see. The
    Morgan fingerprint is a 1024-bit hash of local atomic environments
    suitable as input to a classical learner.
  ],
) <fig:reps>

SMILES, the *Simplified Molecular-Input Line-Entry System*, was
introduced by David Weininger in 1988 at Daylight Chemical Information
Systems and has been the default text encoding for small molecules
ever since. The rules are compact. Atoms are written as element
symbols, with uppercase indicating aliphatic carbon, nitrogen, oxygen,
and so on, and lowercase indicating aromatic atoms within a Kekulé
ring. Bonds are explicit only when they are not single — a double
bond is `=`, a triple bond is `#`, and an aromatic bond inside a
lowercase ring is implicit. Branches are written in parentheses; rings
are closed by matching digit labels. Water is `O`, methane is `C`,
ethanol is `CCO`, benzene is `c1ccccc1`, aspirin is
`CC(=O)Oc1ccccc1C(=O)O`. The same molecule has more than one
syntactically valid SMILES, so a canonical form — defined by the
Daylight canonicalisation algorithm — is used whenever equality
testing matters. *Canonical SMILES* are the database keys for
PubChem, ChEMBL, and DrugBank; non-canonical forms are what users
type at a terminal.

InChI and InChIKey are the IUPAC alternatives. *InChI*, the
International Chemical Identifier, is similar to canonical SMILES in
intent but is layered — separate layers for the formula, the
connection table, the hydrogen layer, the charge, the stereochemistry,
and the isotopes — and includes explicit normalisation of tautomers
and protonation states. *InChIKey* is a 27-character hash of an
InChI, designed so that two identical structures produce the same
key regardless of how they were drawn. InChIKey is what database
back-ends use as the primary key when joining records across PubChem,
ChEMBL, and the literature. SMILES is what humans read; InChIKey is
what databases compare.

The *molecular graph* is the canonical internal representation.
Nodes carry per-atom features — element, formal charge, hybridisation,
aromaticity, implicit hydrogen count, ring membership, stereochemistry.
Edges carry per-bond features — single, double, triple, or aromatic;
in-ring or not; E or Z if double; cis or trans if relevant. Most
implementations store implicit hydrogens — they are computed from
the explicit valence of each heavy atom and the bonds drawn between
them — because explicit hydrogens roughly double the graph size with
no information gain. RDKit, the open-source standard at this layer,
parses SMILES into a graph in microseconds and exposes it through
a Python API that pharma cheminformatics teams have largely
standardised on since the late 2010s. A typical first-day RDKit
session walks the conversion: `mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")`,
then `Descriptors.MolWt(mol)` returns 180.16, `Descriptors.MolLogP(mol)`
returns about 1.2, `Descriptors.TPSA(mol)` returns 63.6, and
`Lipinski.NumRotatableBonds(mol)` returns 3. From the graph, every
downstream representation is computable.

#tip[
  The single most useful early-pipeline cheminformatics sanity check
  is to round-trip a SMILES through RDKit and compare the canonical
  form to the input. Differences expose tautomers, atom-order
  variants, and the occasional malformed string. The cost is one
  function call; the payoff is catching the bugs that bite a week
  later when the same compound appears twice in a screen, once by
  each of two near-canonical SMILES.
]

Property descriptors fall out of the graph trivially. Molecular weight
is a sum over atomic masses. The *octanol-water partition coefficient*,
logP, is estimated by Crippen's atomic contribution scheme:
break the molecule into atom-type fragments, look up each fragment's
logP increment in a table, and sum. *Topological polar surface area*
TPSA is a similar fragment-sum, intended as a proxy for the
membrane-permeability barrier. *Hydrogen-bond donor* and *acceptor*
counts are graph counts of OH and NH versus O and N. *Rotatable
bonds* are non-ring, non-terminal single bonds where rotation produces
distinct conformers. These descriptors are not exact physical
measurements — they are surrogate scores tuned on training sets of
drug-like compounds — but they are fast, deterministic, and accurate
enough to triage libraries of millions of molecules.

The fingerprint is the third stage. It is a fixed-length, typically
binary, vector that encodes the presence or absence of specific
structural features. Two families dominate. *MACCS keys*, from MDL's
Molecular ACCess System, are 166 hand-curated structural patterns —
"contains an aromatic ring," "has a nitrogen attached to an oxygen,"
"contains a halogen at a benzylic position." They are interpretable,
compact, and adequate for crude first-pass filtering. *Morgan
fingerprints*, also called *Extended Connectivity Fingerprints*
(ECFP), are the modern default. Rogers and Hahn formalised the
algorithm in 2010, building on Harry Morgan's 1965 atom-numbering
work for canonicalisation. Each atom gets an initial identifier
derived from its element, charge, degree, hydrogen count, and
aromaticity. At each round, every atom's identifier is replaced by a
hash of its previous identifier together with the sorted identifiers
of its neighbours. After two or three rounds — the typical "ECFP4"
uses radius 2 — every atom carries an identifier that encodes its
environment out to that radius. The union of all atom identifiers
across all rounds is hashed into a 1024- or 2048-bit vector.

#figure(
  image("../../diagrams/lecture-26/03-morgan.svg", width: 92%),
  caption: [
    The Morgan / ECFP algorithm. Atom identifiers are iteratively
    refined to include the neighbourhood out to radius $r$; the union
    of identifiers across rounds is hashed into a fixed-length bit
    vector.
  ],
) <fig:morgan>

The result has two useful properties. First, similar molecules produce
similar fingerprints because the atom-environment dictionary they
hash from is largely shared. Second, the hash is invariant to the
SMILES atom ordering — two valid SMILES for the same molecule yield
the same fingerprint. Variants exist: *count-based ECFP* uses the
multiplicity of each substructure rather than a binary bit; *Avalon*
adds path-based features; *MAP4* uses MinHash signatures. The choice
is a feature-engineering decision tied to the downstream task.

#note[
  The Morgan algorithm and Weisfeiler — Lehman graph isomorphism test
  are the same procedure. Both refine node identifiers by hashing in
  sorted neighbour identifiers; both run for a fixed number of
  rounds; both produce a canonical signature. The connection is not
  coincidental — modern graph neural networks (notably the GIN
  network of Xu et al., 2019) are designed to be at most as
  expressive as 1-WL, and ECFP is the cheminformatics precursor of
  the same idea. The chapter on ML in genomics (Chapter 16) places
  GNNs in the same architecture-data pairing family as Morgan
  fingerprints sit in here, and the connection runs through the
  Weisfeiler — Lehman trick.
]


== Similarity Search and the Tanimoto Coefficient <sec:similarity>

Once compounds are fingerprints, comparing them is set arithmetic.
The Tanimoto coefficient — known to statisticians as the Jaccard
index, after the Swiss botanist Paul Jaccard's 1901 paper on the
distribution of alpine flora — is the ratio of the size of the
intersection to the size of the union:

$ T(A, B) = abs(A inter B) / abs(A union B) $

For binary fingerprints with $a$ bits set in $A$, $b$ bits set in $B$,
and $c$ bits set in both, $T = c / (a + b - c)$. The coefficient
ranges from 0 (no shared bits) to 1 (identical fingerprints).
$T >= 0.85$ is the operational definition of "near-duplicate" used
for library deduplication. $T >= 0.6$ is the rule-of-thumb for "same
chemical series." $T approx 0.3$ is what random pairs of drug-like
compounds produce. The threshold between informative and noise
depends on the library — for a structurally homogeneous corpus of
kinase inhibitors, random pairs sit higher; for a diverse natural-product
collection, lower.

#figure(
  image("../figures/ch26/f2-tanimoto-distribution.svg", width: 95%),
  caption: [
    Tanimoto similarity distributions on a drug-like library. Random
    pairs peak at about 0.15 and decay sharply by 0.4; pairs from a
    single-target activity set are shifted to the right and peak
    around 0.45, with a long tail toward 0.7.
  ],
) <fig:tanimoto>

Tanimoto is the workhorse of compound triage. It clusters libraries
into chemical series — pick a clustering algorithm (Butina is the
cheminformatics default, single-linkage clustering with a Tanimoto
cut-off; or Taylor — Butina, which is greedier), set the threshold to
0.65 or so, and the output is a tractable number of representative
compounds for ordering. It deduplicates libraries — pairs above 0.95
are very probably the same molecule drawn differently, and a 0.85
cut-off catches near-duplicates that differ by a methyl group or a
counter-ion. It powers *similarity search* on databases — give me
all compounds in ChEMBL within Tanimoto 0.6 of this known active, and
the answer comes back in seconds because the fingerprint comparison
is a 32-word population count rather than a graph isomorphism check.
And it supports *hit expansion* — for a confirmed active, retrieve
the structurally closest compounds in a vendor catalogue, order them,
and test them. Hit expansion is how nearly every modern medicinal
chemistry program begins.

#tip[
  When Tanimoto is too crude — when two compounds share atoms but bind
  the target in different modes, or when two structurally diverse
  compounds bind the same way — switch to a *pharmacophore
  fingerprint*. Pharmacophores abstract over atoms: a hydroxyl at
  position $i$ and an amine at position $j$ might both register as
  "hydrogen-bond donor in this geometric pocket," even though the
  underlying atoms differ. RDKit's `Generate3DDistGeometryFeatures`
  and LigandScout are the standard tools. The cost is the 3-D
  conformer generation; the payoff is that bioisosteres — chemically
  distinct groups that play the same biophysical role — start
  clustering together.
]

The structural analogue from EE worth naming is *locality-sensitive
hashing*. An LSH function is one that maps similar inputs to similar
outputs with high probability and dissimilar inputs to dissimilar
outputs. Morgan fingerprints are an LSH on the molecular graph —
similar molecules produce similar bit patterns; the Tanimoto
coefficient on bit vectors is the dissimilarity score that recovers
the underlying similarity. The cheminformatics literature predates
the LSH formalism by a decade, but the math is the same and the
algorithmic tricks transfer in both directions. MinHash-based
fingerprints (MAP4, MHFP6) make the connection explicit — they are
LSH constructions adapted to molecular graphs and used when the bit
budget is tight or when distance-bounded nearest-neighbour queries
need to be exact.


== Virtual Screening and Docking <sec:docking>

Given a target protein and a library of candidates, *virtual
screening* ranks the library by predicted binding and concentrates
follow-up on the top thousand or so. The two halves of the problem are
geometry and energy. Geometry is *pose prediction*: where, in the
target's binding pocket, does each ligand sit, and in what conformation?
Energy is *scoring*: how strongly does the predicted pose bind, in
free-energy terms? The two are entangled — the best pose is by
definition the one that minimises free energy — and every docking
program is a particular algorithmic compromise between sampling poses
exhaustively enough to find the true minimum and scoring them quickly
enough to triage a million compounds in a manageable wall-clock time.

Target preparation comes first. A PDB structure, or an AlphaFold-2
prediction, is loaded; hydrogens are added (most experimental
structures lack them because X-ray crystallography rarely resolves
hydrogens); partial charges are assigned (Gasteiger for fast work,
AM1-BCC or RESP from a quantum-chemistry calculation when accuracy
matters); waters and ions are decided on case by case — some are
structural and should be kept, others are not; the binding pocket is
identified, either by inspection of a known co-crystal ligand or by a
pocket-finding tool such as fpocket or DoGSite. The search space is
defined as a rectangular box around the pocket, typically 20 to 30
ångströms on a side.

The most-used open-source docker is *AutoDock Vina* (Trott and Olson,
*Journal of Computational Chemistry*, 2010). Vina's energy function
is a sum of six empirical terms — steric attraction, steric repulsion,
hydrophobic contact, hydrogen-bonding, conformational entropy from
rotatable bonds, and an out-of-pocket penalty — calibrated to
reproduce a curated set of experimental binding affinities. The
sampling algorithm is iterated local search with a Markov-chain step:
random perturbation of position, orientation, and torsion angles,
followed by local minimisation. Nine poses are returned by default,
each with a predicted $Delta G$ in kcal/mol. On modern hardware Vina
docks one ligand in roughly one to ten seconds on a CPU and faster
still on its GPU port (AutoDock-GPU). A million-compound virtual
screen with Vina takes a few days on a small cluster — well within
the budget of a typical hit-finding program.

#figure(
  image("../../diagrams/lecture-26/04-vina-docking.svg", width: 95%),
  caption: [
    A Vina docking run. Many candidate poses are generated within a
    bounding box around the pocket; each is scored by the empirical
    energy function; the best pose is returned with its predicted
    binding free energy and key contacts annotated.
  ],
) <fig:vina>

Vina's scoring function is the part that frustrates everyone. The
six terms are accurate to about $plus.minus 2$ kcal/mol, which is
useful for rank-ordering but not for absolute affinity prediction —
two kcal/mol is two orders of magnitude on the binding-constant
scale. The remedies are layered. *Molecular Mechanics with
Generalised Born Surface Area* (MM-GBSA) re-scores the top hundred
poses with a more careful, slower energy function and typically
improves the ranking. *Free Energy Perturbation* (FEP+, Schrödinger's
implementation, or open-source variants) is slower still — minutes
to hours per compound — but gets within 1 kcal/mol of experiment on
favourable systems. *GNINA* (Ragoza et al., 2017) trains a 3-D
convolutional neural network on co-crystal structures and uses it as
a learned scoring function over Vina's poses; on benchmark sets it
outperforms Vina's empirical function by a meaningful margin.

The 2020s frontier is *learned end-to-end docking*. *DiffDock* (Corso
et al., 2023) treats pose prediction as a generative-diffusion problem
on the rigid-body pose plus torsion-angle coordinates — the network
learns to push noise toward the binding pose, conditioned on the
target structure. The score is taken as the model's confidence rather
than a physics-based energy. DiffDock generates poses at GPU speed
and, on the PDBBind benchmark, recovers the experimental pose within
2 Å for a higher fraction of cases than Vina. The caveat is that
DiffDock's confidence is not calibrated to free energy and that the
PDBBind benchmark has known structural overlaps with the training set;
honest evaluation under stricter splits closes much of the headline gap.
*RoseTTAFold All-Atom* (Krishna et al., 2024) and *AlphaFold-3*
(Abramson et al., 2024) extend the structure-prediction networks of
Chapter 15 to joint protein-ligand structure, predicting bound
complexes from the protein sequence and the ligand graph together.

#figure(
  image("../figures/ch26/f3-enrichment-curve.svg", width: 95%),
  caption: [
    Two views of the same virtual-screen output. The enrichment curve
    on the left compares Vina, an ML rescorer, and random ranking by
    fraction of actives recovered as a function of library depth. The
    EF1% bar chart on the right summarises the practically relevant
    quantity for a fixed-order budget.
  ],
) <fig:enrichment>

How well does a virtual screen work? The two standard metrics are
*ROC AUC* and *enrichment factor at the top $x$ %* of the ranked
output, written $"EF"_x$. AUC is the area under the receiver-operator
characteristic curve — actives recall on the $y$-axis, false-positive
rate on the $x$-axis, area integrated over the curve. AUC = 0.5 is
random; AUC = 1.0 is perfect. $"EF"_x$ is the fold-enrichment of actives
among the top $x$ % of compounds compared to random selection —
$"EF"_(1%) = 35$ means the top 1 % of ranked compounds contains 35
times more actives than a random 1 %. For typical kinase or protease
targets with well-characterised pockets, Vina alone reaches AUC =
0.70–0.80 and $"EF"_(1%) = 5"—"15$; an ML rescorer over Vina's top poses pushes
both numbers up. For poorly characterised targets, or for chemotypes
absent from the training data, both numbers collapse.

#warn[
  The DUD-E benchmark (Mysinger et al., 2012), which has been the
  standard ROC test set for docking since the early 2010s, has known
  decoy-bias problems. The decoys are chosen to be 2-D-dissimilar to
  the actives, which makes them easy to discriminate by simple
  features that have nothing to do with binding. Modern docking
  papers report numbers on DUD-E that look impressive and then fail
  on out-of-distribution targets. Treat any docking benchmark with
  the same suspicion you would treat a random 80/20 split in a
  genomics paper (see Chapter 16 on data leakage). The cleaner
  evaluation is *prospective* — design the screen, run it, order
  the top fifty compounds, measure them. Most published EFs are
  retrospective; the prospective numbers are roughly half.
]

The structural analogue from EE worth holding onto is *placement and
routing*. In VLSI design, the engineering problem is to place cells
into chip area subject to spacing constraints and route wires between
them while minimising delay, area, and power. In docking, the problem
is to place a ligand into a pocket subject to clash constraints and
contact requirements while minimising binding free energy. Both
problems are non-convex with many local minima, both are
NP-hard in general, both use iterated local search and simulated
annealing as the workhorse heuristics. The methods transfer; the
intuitions transfer; the failure modes transfer. A docking pose that
optimises the score but leaves a key polar group unsatisfied is the
ligand analogue of a placement that meets timing closure but exceeds
the power budget. Both look fine on paper. Both fail in fabrication.


== ADMET and Druglikeness <sec:admet>

Potency is necessary and not sufficient. A drug that binds its target
in vitro but does not reach the target in vivo, or that reaches the
target but is metabolised in fifteen minutes, or that reaches the
target and stays there but blocks a cardiac potassium channel as a
side effect, will not become a drug. The umbrella term for these
constraints is *ADMET*: Absorption, Distribution, Metabolism,
Excretion, and Toxicity. Absorption is how much of an oral dose
crosses the intestinal wall. Distribution is where in the body the
drug ends up — relevant for central-nervous-system targets, where the
blood-brain barrier filters out the polar majority of drug-like
compounds. Metabolism is hepatic clearance, dominated by the
cytochrome-P450 family of enzymes. Excretion is renal clearance and
half-life. Toxicity is everything that goes wrong — hepatotoxicity
from reactive metabolites, cardiotoxicity from hERG potassium-channel
blockade, mutagenicity, teratogenicity, and the long tail of
idiosyncratic adverse events.

The historical first attempt at quantifying druglikeness is
*Lipinski's Rule of Five* (Lipinski et al., 1997). A retrospective
study of compounds that had reached Phase 2 trials versus those that
had not found that the failures clustered at the extremes of four
simple physicochemical properties: molecular weight, the
octanol-water partition coefficient logP, the hydrogen-bond donor
count, and the hydrogen-bond acceptor count. Lipinski codified the
edges:

- MW $lt.eq$ 500 daltons,
- $log P lt.eq 5$,
- hydrogen-bond donors $lt.eq$ 5,
- hydrogen-bond acceptors $lt.eq$ 10.

A compound that violates two or more of these is "less likely" to be
orally bioavailable. The framing matters — the rules are *not* a
filter that excludes compounds; they are a triage heuristic that
flags ones worth examining more carefully. Roughly 10 % of currently
approved oral drugs violate at least two of the rules. The rules are
biased toward CNS-active small molecules, the chemotype Pfizer was
working in when Lipinski did the study, and they are weaker on
natural products, kinase inhibitors, and protease inhibitors —
classes that routinely violate one or two rules while being
perfectly drug-like in practice.

#figure(
  image("../../diagrams/lecture-26/09-lipinski.svg", width: 92%),
  caption: [
    The Lipinski Rule-of-Five property distributions across about
    2 000 approved oral drugs. Roughly 90 % satisfy all four rules;
    the violations cluster in macrocyclic natural products and a
    handful of recent kinase inhibitors.
  ],
) <fig:lipinski>

The single-number successor is *QED*, the *Quantitative Estimate of
Druglikeness* (Bickerton et al., 2012). QED combines eight properties
— MW, logP, HBD, HBA, polar surface area, rotatable bonds, aromatic
ring count, and a count of structural alerts — through a desirability
function that maps each property to a 0–1 score based on the
empirical distribution across approved drugs. The geometric mean of
the eight desirabilities is the QED. Drug-like compounds score 0.5
to 0.9; non-drug-like compounds score below 0.3. QED is what most
modern library-filtering pipelines use, often in addition to the
unmodified Lipinski rules.

Specific toxicity flags are a separate axis. *PAINS*
(*Pan-Assay Interference Compounds*, Baell and Holloway, 2010) is a
set of about 480 structural alerts that flag compounds known to
produce promiscuous false positives in screens — reactive Michael
acceptors, redox-cycling quinones, catechols, rhodanines. PAINS hits
should be triaged out before order, not because they are guaranteed
artefacts but because the conditional probability that an apparent
PAINS-hit is a real ligand is low and the cost of follow-up is
high. *hERG* is the human Ether-à-go-go-Related Gene channel; small
molecules that block hERG can prolong the cardiac QT interval and
trigger arrhythmia. hERG predictors range from simple QSAR models
trained on the published patch-clamp data to GNNs that incorporate
docked poses into the channel's cryo-EM structure. *AMES* is the
bacterial reverse-mutation assay; *AMES predictors* flag mutagenicity
risk. *Hepatotoxicity predictors* are the noisiest of the bunch
because the underlying mechanisms (reactive-metabolite formation,
mitochondrial dysfunction, immune-mediated injury) are diverse and
the assay endpoint (drug-induced liver injury, DILI) is rare,
multi-causal, and frequently misattributed.

#figure(
  image("../../diagrams/lecture-26/05-admet-radar.svg", width: 92%),
  caption: [
    An ADMET radar chart, the visual triage tool that has been the
    medicinal-chemistry standard since the early 2000s. Larger
    radius is better. Drug-like compounds approach a regular hexagon;
    candidates with one or two failing axes show as collapsed lobes.
  ],
) <fig:radar>

Modern ADMET prediction is increasingly ML-driven. *Chemprop*
(Yang et al., 2019) is a directed message-passing neural network
that takes a molecular graph and outputs a scalar property
prediction; it has been trained on most of the public ADMET endpoints
(*MoleculeNet*, Wu et al., 2018) and routinely outperforms classical
QSAR on cross-target benchmarks. *ADMETlab 2.0* and the *Therapeutics
Data Commons* (Huang et al., 2021) provide ensemble predictions
across roughly 50 ADMET endpoints through web interfaces with
known confidence intervals. The 2024 frontier is multi-task
foundation models — *ChemBERTa* (Chithrananda et al., 2020),
*MolBERT* (Fabian et al., 2020), *MoLFormer* (Ross et al., 2022) —
pretrained on hundreds of millions of unlabelled SMILES with masked
language modelling, then fine-tuned on the ADMET endpoints jointly.
The Chapter 16 lesson applies directly: pretraining buys an order of
magnitude in label efficiency when downstream data is scarce, which
is the case for almost every internal ADMET dataset.

#warn[
  ADMET predictors share a label-scarcity pathology with the
  genomics models of Chapter 16. Public ADMET datasets are small
  ($10^3$ to $10^4$ compounds per endpoint) and heavily biased
  toward whatever a few large pharma libraries are. Compounds in a
  novel scaffold class — exactly the class a generative model is
  trying to design into — sit far from the training distribution
  and produce predictions that are confidently wrong. Treat ADMET
  predictions on novel scaffolds as approximate prior probabilities
  rather than measurements. The first wet-lab data on a new
  chemotype usually rebuilds the predictor for that chemotype from
  scratch.
]


== Deep Learning for Drug Discovery <sec:deep-learning>

Graph neural networks were the first deep-learning architecture to
land in cheminformatics, because molecules are graphs and the
inductive bias is exact. *Message-passing neural networks* (Gilmer
et al., 2017) generalised the earlier convolutional-graph approaches
of Duvenaud, Kipf-Welling, and others into a single framework: at
each layer, every node aggregates a learned function of its
neighbours' features and bond features; after several layers, every
node carries a representation of its neighbourhood out to that
radius; a final readout function reduces all node representations
into a single graph-level vector that feeds an MLP for the prediction.
*Chemprop* (Yang et al., 2019) is a directed MPNN — messages run
along bonds rather than between atoms — and remains the open-source
reference implementation for molecular property prediction. The
inductive biases are clean: permutation equivariance because the
network's output is invariant to atom relabelling, locality because
the message-passing range is bounded, and bond-order awareness
because edges carry features.

#figure(
  image("../../diagrams/lecture-26/10-gnn-message-passing.svg", width: 92%),
  caption: [
    GNN message passing on a small molecule. Atom representations are
    iteratively refined by aggregating learned functions of neighbour
    atoms and bonds; the final readout produces a graph-level
    representation that feeds the property prediction head.
  ],
) <fig:gnn>

For *molecular property prediction* — predicting IC50, solubility,
logD, blood-brain-barrier permeability, hERG inhibition, and the rest
— Chemprop and its descendants compete with fingerprint-plus-random-forest
baselines and often win, but the margin is task-dependent and the
baselines are stronger than they look. The *MoleculeNet* benchmark
(Wu et al., 2018) standardised the comparison and produced a
sometimes uncomfortable conclusion: on tasks with few thousand
labels and well-defined chemistries, random-forest-on-ECFP is hard
to beat. On tasks with tens of thousands of labels and diverse
chemistries, GNNs pull ahead. On tasks with bespoke pretraining and
careful augmentation, foundation-model fine-tuning pulls further
ahead. The architecture-matters-less-than-the-data lesson from
Chapter 16 reappears.

The architectural variant worth knowing about is the *3-D equivariant*
GNN. Where the 2-D message-passing networks of Gilmer and Yang take
the molecular graph as input, 3-D equivariant networks (*SchNet*,
*DimeNet*, *NequIP*, *MACE*) take atomic coordinates and learn
features that are equivariant under rotation and translation of the
input. They are the cheminformatics analogue of AlphaFold-2's
structure module — same SE(3) group, same equivariance machinery,
same payoff in data efficiency on the right tasks. For property
predictions that depend on conformer geometry (binding affinity for a
specific pocket, partition coefficients in lipid bilayers) the 3-D
equivariant networks pull ahead of 2-D message-passing on benchmarks
that are honest about the conformer-generation step. For 2-D-only
tasks, the simpler 2-D networks win on compute and on data efficiency.

*Generative models* are the second class of architecture worth
covering, because they invert the inference problem. Where a property
predictor takes a molecule and outputs a score, a generative model
takes a target score and outputs a molecule. *REINVENT* (Olivecrona
et al., 2017; Loeffler et al., 2024) trains a recurrent neural network
to generate SMILES, then fine-tunes it with reinforcement learning
against a reward — predicted potency, predicted ADMET, scaffold
novelty. The output of a REINVENT run is roughly $10^5$ to $10^6$
generated SMILES, of which a few thousand have a non-zero chance of
synthesisability and a few hundred are worth a medicinal chemist's
attention. *JT-VAE* (Jin, Barzilay, Jaakkola, 2018) uses a
junction-tree variational autoencoder — molecules are decomposed
into fragments, fragments are arranged in a tree, the tree is the
latent code — to support property-conditioned generation in a
chemically valid latent space. *MolDiff* (Li et al., 2023) and *EDM*
(Hoogeboom et al., 2022) are diffusion models on molecular graphs and
3-D structures respectively, descendants of the image-diffusion
machinery now standard in generative AI.

#figure(
  image("../../diagrams/lecture-26/12-alphafold-design.svg", width: 95%),
  caption: [
    The AlphaFold-driven design pipeline. Target sequence to predicted
    3-D structure to detected pocket to pocket-conditioned generation
    to ADMET filter to docking to synthesis. Insilico Medicine's
    INS018_055 — the first AI-discovered drug to reach Phase 2 —
    came out of a version of this pipeline in roughly 18 months.
  ],
) <fig:alphafold-design>

The most consequential development of the 2020s is *structure-based
generative design*. The recipe combines AlphaFold-2 (or
AlphaFold-3, or RoseTTAFold) for target structure, a pocket-finder
for binding-site detection, and a pocket-conditioned generative model
— *Pocket2Mol* (Peng et al., 2022), *DiffSBDD* (Schneuing et al.,
2023), *TargetDiff* (Guan et al., 2023) — to produce molecules that
fit the pocket geometry. The output is a small set of designed
ligands; the validation passes through docking, ADMET prediction, and
synthesis. The economic claim is that this pipeline compresses the
target-to-IND timeline from five-plus years to roughly eighteen
months. *Insilico Medicine's INS018_055*, an inhibitor of TRAF2- and
NCK-interacting kinase for idiopathic pulmonary fibrosis, is the
worked example: AI-identified target, AI-designed compound,
AI-prioritised animal experiments, IND filing in roughly 30 months
total. The compound entered Phase 1 in 2023 and Phase 2a in 2024 —
the first molecule designed by an end-to-end ML pipeline to reach
patient testing.

#note[
  Whether Insilico's INS018_055 will become a drug is the wrong
  question to ask. The right question is whether the compression of
  target-to-IND from five years to 18 months is reproducible across
  programs and across companies. The 2020s answer is "partially."
  The pipeline works best when the target is a kinase or kinase-like
  enzyme with a well-defined pocket and reasonable existing
  literature — exactly the case where classical structure-based
  design also works best. The pipeline struggles on protein-protein
  interfaces, on intrinsically disordered targets, and on
  first-in-class biology where the training distribution is empty.
  The productivity gap in @fig:productivity remains the right
  context for evaluating any single AI-discovery success story.
]


== The Productivity Gap <sec:productivity-gap>

The cost of bringing a single new drug to market has risen, in
inflation-adjusted dollars, by roughly a factor of 80 since 1950.
*Eroom's law*, coined by Jack Scannell and colleagues in 2012
(*Nature Reviews Drug Discovery*) as a deliberate inversion of
Moore's law, names the pattern: in pharma, productivity per R&D
dollar has been halving every nine years for sixty years. *Hansen*
(1979) estimated the cost per new drug at about \$54 million in 1976
dollars, equivalent to roughly \$250 million today. *DiMasi et al.*
(2003 and 2016 updates) put it at \$800 million and then \$2.6
billion in subsequent revisions. *BCG and Deloitte* sector reports
since 2020 have placed the inflation-adjusted figure between \$2 and
\$3 billion, with much of the increase concentrated in larger
clinical trials, more demanding regulatory expectations, and a
shifting target portfolio toward harder-to-drug biology.

#figure(
  image("../figures/ch26/f5-productivity-gap.svg", width: 95%),
  caption: [
    Eroom's law in a single plot. The cobalt line tracks the
    inflation-adjusted cost per approved drug from the 1950s to the
    early 2020s, rising roughly two orders of magnitude. The amber
    dotted line tracks the inverse — approvals per inflation-adjusted
    billion dollars of R&D. The shaded region marks the post-2018
    AI-led era, where the trend has, at most, paused.
  ],
) <fig:productivity>

The diagnosis is contested. Scannell's "four reasons" — the better-than-the-Beatles
problem (each new drug must beat the previous standard of care, a
moving bar), the cautious-regulator problem (post-Vioxx FDA
expectations are higher), the throw-money-at-it problem (large
companies over-spend on programs because they can afford to), and the
basic-research-bias problem (large companies optimise for citations,
not approvals) — has been the standard reference for a decade. None
of those reasons is something computational methods can affect
directly. *Computational methods cheapen the upper bands of the
funnel — virtual screening, ADMET triage, structure-based design —
but most spend is in the lower bands*, and most attrition is in Phase 2,
where the question is whether the biology works in humans. The
unfortunate corollary is that an AI pipeline that produces a Phase 1
candidate in 18 months instead of 5 years saves a few tens of
millions of dollars, against a clinical-trial bill of a billion
dollars or more.

The optimist's case is layered. First, the AI-led era has not yet
reached the stage where its successes or failures will dominate the
approval statistics — INS018_055 will hit Phase 2 readouts in 2025,
the second-generation AI-designed compounds in 2027 and 2028, and so
on. Second, AI may shift the *composition* of the target portfolio
toward harder-to-drug proteins, which has knock-on effects on the
pre-clinical work even when the clinical attrition is unchanged.
Third, *combination therapy* — designing two complementary compounds
together, where one rescues the off-target liability of the other —
is exactly the kind of multi-objective search the ML stack is good
at and the classical stack is not. Fourth, *failure prediction* is
itself a useful target — if an ML model can flag a Phase 2 failure
six months earlier than a human reviewer, the savings are real even
when the approval rate is unchanged.

The pessimist's case is shorter. *Drug attrition is biology, not
chemistry*, and ML has not yet demonstrated a meaningful effect on
the biology. Most of the published AI-discovered compounds belong to
chemotype families and target classes that were already well-served
by classical structure-based design. The 2024 trough of
disillusionment for AI drug discovery — multiple high-profile
program failures, layoffs at Recursion and BenevolentAI, public
recriminations — is consistent with the optimist's view (premature
expectations on long-cycle work) and with the pessimist's view (the
gains were small and concentrated in the upper funnel). The next five
years will resolve which is closer to the truth.


== Hit-to-Lead and Onwards <sec:hit-to-lead>

A confirmed hit has measurable activity, typically in the
single-digit micromolar range. A *lead* has, in addition to that
activity, sub-micromolar potency, 10×–100× selectivity over the most
relevant off-targets, plausible ADMET, and a defensible synthetic
route. The *hit-to-lead* phase is where most of the
medicinal-chemistry work happens: build a *structure-activity
relationship* (SAR) by synthesising and measuring small variations
around the hit, learn which substituents tune which properties, and
move along the analogue series until a small number of compounds clear
the lead bar.

#figure(
  image("../../diagrams/lecture-26/11-hit-to-lead.svg", width: 92%),
  caption: [
    A hit-to-lead analogue series. Six substituent variations on a
    confirmed micromolar hit explore three orders of magnitude of
    potency and an order of magnitude of selectivity. The best lead
    enters the lead-optimisation phase that follows.
  ],
) <fig:h2l>

The ML contribution at this stage is *active learning*. An active
predictor — a GNN trained on the SAR data accumulated so far —
estimates not only the predicted activity for an unmade analogue but
also the predicted uncertainty. The next batch of analogues to
synthesise is then chosen to maximise expected information gain:
spread choices across the activity-uncertainty frontier rather than
just picking the predicted-most-potent. This is *Bayesian
optimisation* applied to chemistry, exactly the procedure used in
adaptive experimental design in any other engineering discipline.
The payoff is order-of-magnitude reductions in the number of
compounds that have to be synthesised to reach a given lead-quality
threshold. Modern pharma groups routinely run loops where ML proposes
the next set, chemists synthesise it, assays measure it, and the model
is updated for the next round.

*Lead optimisation* is the next stage and the most expensive of the
pre-clinical phases. Iterative medicinal chemistry on each lead aims
to improve potency a further 10×–100×, optimise the ADMET profile,
remove specific liabilities (a hERG hit, a CYP3A4 induction signal,
a hepatotoxicity flag), and arrive at a defensible intellectual
property position. *Free-energy perturbation* (FEP+) is the
computational workhorse — it predicts the relative binding affinity
of two closely related compounds to within 1 kcal/mol on favourable
systems and reduces the synthesis burden by triaging analogues
in silico before they are made. *Crystallography* on the lead-target
complex gives structural feedback for the next analogue round.
*Pre-clinical pharmacokinetics* and toxicology in rodents and a
non-human primate confirm — or refute — the predicted ADMET. The
output is the *development candidate*: one or at most three
compounds that earn an IND filing.

The clinical work is, computationally, less interesting. Phase 1 in
healthy volunteers establishes the human safety profile and the
pharmacokinetics. Phase 2 in patients establishes efficacy in the
target indication. Phase 3 confirms efficacy in a much larger
population and confirms the safety profile under realistic conditions
of use. Each phase takes two to four years. The 90 % Phase 2 failure
rate — a number that has been roughly stable for thirty years
despite enormous investment in pre-clinical predictive methods —
remains the central problem in the economics of drug discovery. It
is also the one place where the computational toolkit has had
essentially no measurable impact. The biology of the target in the
human disease is the thing the model does not know.


== Pitfalls, in a List <sec:pitfalls>

A working drug-discovery ML project carries the same kind of
recurring failure modes that the Chapter 16 design-review checklist
catches in genomics. The shape of the failure modes is recognisable
across the field, and the same questions ask before the experiment
and re-ask before reading another team's headline number.

#figure(
  image("../figures/ch26/f4-pitfall-checklist.svg", width: 95%),
  caption: [
    A design-review checklist for a drug-discovery ML project. Five
    questions; five failure modes. Most "AI for drug discovery"
    controversies reduce to one of them.
  ],
) <fig:ch26-checklist>

*Representation choice.* Fingerprints for fast similarity search
and as a strong baseline. Graph neural networks for SAR with held-out
scaffolds, where the inductive bias matches the data. 3-D equivariant
networks for binding-mode questions and conformer-dependent
properties. SMILES-LMs for generation, where the autoregressive
formulation matters. Each representation has a domain where it is
strongest; using the wrong one is the closest analogue in this
chapter to the architecture-mismatch failure in Chapter 16. If your
input is a 3-D pocket and you handed it to a SMILES-LM, the network
has to rediscover the pocket from sequence; it will, but it will
need ten times the data.

*Loss specification.* Predict $log("IC50")$ — equivalently $-log("IC50")$,
written pIC50 — rather than the raw IC50. The dynamic range of IC50
spans five to six orders of magnitude, and an MSE on raw values is
dominated by the few millimolar inactive compounds. Handle censoring
explicitly: an assay measurement reported as "> 10 µM" is a one-sided
constraint, not a single number, and treating it as the latter biases
the predictor toward optimism. Use Huber or quantile losses when the
assay's outlier distribution warrants it. The Chapter 16 lesson on
matched likelihoods — Gaussian for continuous, negative binomial for
counts, ZINB for zero-inflated counts — extends here to Huber for
censored ADMET data and ordinal regression for ranked toxicity
endpoints.

*Split protocol.* The most-published-on pitfall in this section.
*Random splits leak via scaffold homology* — Bemis – Murcko scaffolds
group the molecules in any drug-discovery dataset into a few hundred
classes, and a random 80/20 split places most test compounds within
the same scaffold as a training example. The bare-minimum remedy is a
*scaffold split* (Bemis – Murcko clustering, then split at the cluster
level), used by MoleculeNet and the Therapeutics Data Commons. The
stronger remedy is a *time split*: train on compounds with assay
dates before some cut-off, test on compounds added after. The
strongest remedy is a *target-leave-out* protocol when the question
is cross-target generalisation. The protocol asymmetry between
scaffold and random splits is typically 10–20 percentage points of
AUC; results that look striking under random and collapse under
scaffold split are almost always due to leakage rather than learning.

*Docking-pose plausibility.* A high Vina score on an implausible
pose is a Vina artefact, not a binding prediction. Inspect the top
ten poses by eye before trusting the rank. Check for clashes, for
buried unsatisfied polar groups, for ligand conformations that
require torsion angles no organic chemist would accept. Use
*re-docking RMSD* on a known co-crystal as a system-specific
calibration: if the docker cannot recover the experimental pose
within 2 Å on the training target, its predictions on novel ligands
should be discounted accordingly. The DiffDock and AlphaFold-3
era has not removed this requirement; it has shifted it from
"inspect every pose" to "inspect a strong sample of poses and
back-stop with re-docking RMSD on the held-out co-crystals."

*ADMET as a multi-objective constraint.* The temptation in lead
optimisation is to pick the single most potent analogue from each
round. The trap is that the most potent analogue is often the one
with the worst hERG, the worst solubility, or the worst microsomal
stability — chemical features that drive potency (lipophilicity,
specific reactive groups, rigid scaffolds) drive the off-target
liabilities too. The fix is to track the *Pareto frontier* in the
multi-objective space and pick the next batch from the frontier
rather than from a single axis. *Uncertainty-aware* predictors
(deep ensembles, MC-dropout, Gaussian-process regression on
fingerprints) make the active-learning selection rigorous: pick
points where expected improvement is highest given the multi-axis
uncertainty.

#warn[
  The pitfall list in @fig:ch26-checklist deliberately mirrors the
  Chapter 16 genomics-ML checklist. The five categories — wrong
  representation, wrong loss, leaky split, weak baseline / weak
  validation, multi-objective collapse — generalise across most
  applied-ML failure modes in this book. The vocabulary differs by
  field; the structural shape of the failures does not.
]


== Summary <sec:ch26-summary>

- Drug discovery is constrained multi-objective optimisation in a
  discrete chemical space of roughly $10^60$ drug-like molecules.
  The constraints are potency on the target, selectivity against
  off-targets, and the ADMET cluster of pharmacokinetic and toxicity
  properties. The classical pipeline runs 10–15 years and costs
  \$1–3 billion per approved drug, with about 90 % Phase 2 attrition.
- Molecules become SMILES (storage), molecular graphs (working
  representation), and fingerprints (ML features). Morgan / ECFP
  fingerprints (Rogers and Hahn, 2010) are the modern default for
  similarity-based methods; they are a locality-sensitive hash of
  the molecular graph and share a procedural lineage with the
  Weisfeiler – Lehman test that modern graph neural networks
  generalise.
- The Tanimoto coefficient — Jaccard's 1901 index on bit vectors —
  is the standard fingerprint-similarity metric. Threshold conventions
  are operational rather than universal: $T >= 0.85$ for
  near-duplicates, $T >= 0.6$ for same chemical series,
  $T approx 0.3$ for random drug-like pairs.
- AutoDock Vina (Trott and Olson, 2010) is the open-source virtual
  screening default. It is iterated local search on an empirical
  six-term scoring function, accurate to $plus.minus 2$ kcal/mol —
  good for ranking, weak for absolute affinity. ML rescorers (GNINA),
  MM-GBSA re-scoring, and end-to-end diffusion-based docking
  (DiffDock, RoseTTAFold All-Atom, AlphaFold-3) extend the toolkit.
  Honest evaluation is by prospective enrichment factor, not by
  retrospective DUD-E ROC.
- Lipinski's Rule of Five (1997) and QED (Bickerton, 2012) are the
  druglikeness heuristics that filter libraries before order. Modern
  ADMET prediction is ML-driven (Chemprop, MoleculeNet, Therapeutics
  Data Commons); foundation-model fine-tuning (ChemBERTa, MoLFormer)
  pays off in the low-label regime that defines most internal
  endpoints.
- Graph neural networks (Gilmer et al., 2017; Chemprop / Yang et al.,
  2019) match molecular structure prediction tasks. Generative
  models (REINVENT, JT-VAE, MolDiff, DiffSBDD) invert the inference
  problem: given a desired property profile and a target pocket,
  design molecules that satisfy them. AlphaFold-driven pocket-conditioned
  generation is the current frontier; Insilico Medicine's INS018_055
  is the first compound from an end-to-end ML pipeline to reach
  Phase 2 patient testing.
- Eroom's law (Scannell et al., 2012) names the long-run trend:
  drug-discovery productivity per R&D dollar has been halving every
  nine years for sixty years. ML cheapens the upper bands of the
  attrition funnel — virtual screening, ADMET triage,
  structure-based design — and has not yet measurably affected
  the dominant cost at the bottom, where Phase 2 attrition is
  driven by human biology rather than by chemistry.
- The recurring failure modes are familiar from Chapter 16: wrong
  representation, mis-specified loss, scaffold-leaky split, weak
  baseline, multi-objective collapse to a single axis. The
  design-review checklist in @fig:ch26-checklist asks the five questions
  worth running before any program and before reading any paper's
  headline number.


== Exercises <sec:ch26-exercises>

#strong[1.] #emph[Representation conversion.] Take five drugs of
your choice — aspirin, ibuprofen, atorvastatin, sildenafil, and a
fifth of your choosing. For each, write the SMILES from memory or
from a database, parse with RDKit, compute molecular weight, logP,
TPSA, HBD, HBA, and the rotatable-bond count. Compute the Morgan
fingerprint at radius 2 and 1024 bits. Round-trip each compound
through RDKit's canonical-SMILES output and confirm the canonical
form is stable. Report which descriptors differ from PubChem's
published values and by how much.

#strong[2.] #emph[Tanimoto distribution.] Download the FDA-approved
drugs subset of DrugBank (about 2 000 compounds). Compute pairwise
Morgan fingerprints and Tanimoto similarities. Plot the distribution
of similarities. Identify the highest-similarity pair that is not a
stereoisomer or salt of the same active ingredient — what do the
two compounds have in common? Cluster the library at $T = 0.65$ and
report the five largest clusters; identify the chemical series each
represents.

#strong[3.] #emph[Vina virtual screen.] Pick a kinase target with a
co-crystal structure in the PDB — ABL1 in complex with imatinib
(PDB 1IEP) is a friendly choice. Prepare the receptor with AutoDock
Tools or `meeko`. Dock a small library (100 ZINC compounds plus
imatinib as a positive control). Rank by predicted $Delta G$ and
report the top ten. Compute the Vina pose for imatinib and the
RMSD between Vina's pose and the crystallographic pose. Comment on
whether the docker recovered the experimental pose; if not, identify
which torsion angle is most off.

#strong[4.] #emph[Scaffold split.] Take the MoleculeNet BACE
(beta-secretase 1 inhibition) dataset. Train a Chemprop GNN under
two splits: a random 80/20 and a Bemis – Murcko scaffold split.
Report ROC-AUC under both. Predict, before running, how large the
gap will be; reflect on the difference between your prediction and
the measurement. Compare to the random-forest-on-ECFP baseline on
the same splits.

#strong[5.] #emph[Lipinski and QED on a generated library.] Run
REINVENT (or a similar SMILES-LM generator) for 1 000 sampling steps
with a Tanimoto-similarity reward against a known kinase inhibitor.
Compute Lipinski Rule-of-Five compliance and QED on the generated
set. Compare the distributions to the FDA-approved subset of
DrugBank. Where does the generator over- or under-produce relative
to the drug-like prior? Hypothesise why.

#strong[6.] #emph[ADMET prediction audit.] For 20 compounds from
ChEMBL that have experimental hERG IC50 measurements, run an
open-source hERG predictor (the ADMET-AI service or the Therapeutics
Data Commons predictor). Plot predicted vs measured pIC50.
Compute the Spearman rank correlation. Identify the two compounds
with the largest disagreement and propose, in one sentence each,
why the predictor failed.

#strong[7.] #emph[Active learning sketch.] You have an ML model that
predicts IC50 for a series of analogues against your target, with
uncertainty estimates. Your synthesis budget is 20 compounds in the
next round. Describe — in one page — the policy you would use to
choose those 20, justifying the trade-off between exploitation
(synthesising the highest predicted-potency compounds) and exploration
(synthesising the highest-uncertainty compounds). Compare your policy
to expected-improvement and Thompson-sampling baselines.

#strong[8.] #emph[(Open-ended.)] Pick a 2023–2025 paper claiming a
striking result from an AI-led drug-discovery pipeline. Run the
five-question design-review checklist (@fig:ch26-checklist) against the
paper. Identify the weakest link. Propose one experiment that would
materially change your confidence in the headline result, and predict
whether it would raise or lower the number. Cite the paper.


== Further Reading <sec:ch26-further-reading>

- *Weininger, D.* (1988). "SMILES, a chemical language and information
  system." _Journal of Chemical Information and Computer Sciences_
  28: 31 – 36. The original SMILES paper. Useful for the rule set;
  the Daylight reference manual is the practical companion.
- *Rogers, D., & Hahn, M.* (2010). "Extended-connectivity
  fingerprints." _Journal of Chemical Information and Modeling_ 50:
  742 – 754. The Morgan / ECFP paper. Read for the algorithm and for
  the discussion of why circular signatures generalise.
- *Trott, O., & Olson, A. J.* (2010). "AutoDock Vina: improving the
  speed and accuracy of docking with a new scoring function,
  efficient optimization, and multithreading." _Journal of
  Computational Chemistry_ 31: 455 – 461. The Vina paper. Still the
  most-cited docking reference; the empirical scoring function and
  the iterated-local-search algorithm are both clearly described.
- *Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J.*
  (1997). "Experimental and computational approaches to estimate
  solubility and permeability in drug discovery and development
  settings." _Advanced Drug Delivery Reviews_ 23: 3 – 25. The
  original Rule-of-Five paper. The thirty-year retrospective is worth
  reading; the rules have aged less well than the framing.
- *Bickerton, G. R., Paolini, G. V., Besnard, J., Muresan, S., &
  Hopkins, A. L.* (2012). "Quantifying the chemical beauty of drugs."
  _Nature Chemistry_ 4: 90 – 98. The QED paper. The desirability-function
  formalism is the part to internalise; the specific weights are
  domain-specific.
- *Gilmer, J., Schoenholz, S. S., Riley, P. F., Vinyals, O., & Dahl,
  G. E.* (2017). "Neural message passing for quantum chemistry."
  _Proceedings of the 34th International Conference on Machine
  Learning_, PMLR 70: 1263 – 1272. The MPNN paper. Pair with Yang et
  al. (2019) on Chemprop for the directed variant that became the
  cheminformatics standard.
- *Olivecrona, M., Blaschke, T., Engkvist, O., & Chen, H.* (2017).
  "Molecular de-novo design through deep reinforcement learning."
  _Journal of Cheminformatics_ 9: 48. The original REINVENT paper.
  Read the 2024 Loeffler-et-al. update for the production-grade
  successor.
- *Corso, G., Stärk, H., Jing, B., Barzilay, R., & Jaakkola, T.*
  (2023). "DiffDock: diffusion steps, twists, and turns for
  molecular docking." _International Conference on Learning
  Representations_. The DiffDock paper. Read with the 2024
  RoseTTAFold All-Atom and AlphaFold-3 papers for the current
  state of end-to-end docking.
- *Scannell, J. W., Blanckley, A., Boldon, H., & Warrington, B.*
  (2012). "Diagnosing the decline in pharmaceutical R&D efficiency."
  _Nature Reviews Drug Discovery_ 11: 191 – 200. The Eroom's-law
  paper. The four reasons hold up; the AI-led pipeline is the
  natural fifth chapter that has not yet been written.
- *Wu, Z., Ramsundar, B., Feinberg, E. N., et al.* (2018).
  "MoleculeNet: a benchmark for molecular machine learning."
  _Chemical Science_ 9: 513 – 530. The benchmark paper. The
  Therapeutics Data Commons (Huang et al., 2021) is the modern
  successor; both are worth reading for the split protocols as much
  as for the headline tables.
