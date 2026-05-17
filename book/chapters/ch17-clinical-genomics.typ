#import "../theme/book-theme.typ": *

= Clinical Genomics, Variant Interpretation, and Ethics <ch:clinical-genomics>

#matters[
  Every preceding chapter treated bioinformatics as research. This one
  treats it as regulated software. When a pipeline's output informs a
  diagnosis, a dosing decision, or a surgical plan, the rules change
  underneath you. The accuracy targets are individual rather than
  aggregate; the reruns must be bit-identical years later; the change
  log is the audit trail; and every report carries the signature of a
  physician whose license is at stake. The technical content of this
  chapter — the ACMG/AMP rule table, the evidence ecosystem,
  pharmacogenomics, the regulatory landscape, ancestry bias, indigenous
  data sovereignty — is the technical content of an industry that
  treats your code as a medical device. If you intend to ship genomics
  software into care, the chapter is the price of admission.
]

A clinical genomics laboratory is a research bioinformatics pipeline
with the regulatory shell of a medical device. The sequencing chemistry
is the same; the aligners and variant callers are usually the same; the
reference databases — ClinVar, gnomAD, REVEL, AlphaMissense, SpliceAI —
are the same. What changes is the discipline around them. Every tool
version is pinned. Every reference file is hashed. Every run has a UUID
that traces back through the chain of custody to the tube of blood that
arrived in the receiving room. Every classification has a written
rationale. Every report has a signatory.

The chapter walks the stack the same way a clinical molecular geneticist
walks it. We start at the regulatory boundary and the engineering
discipline it forces. We then move into the ACMG/AMP variant
classification framework — the rule table that has, since 2015, served
as the operating standard for almost every clinical lab in the world.
We unpack the evidence ecosystem the rules consume — population
frequencies, in-silico predictors, splice predictors, prior
classifications, literature — and the working profession of the variant
scientist. We turn to pharmacogenomics, which is the cleanest example
of clinical genomics already operating as a lookup table on a
state-machine, with immunogenomics as the natural extension into HLA
typing and neoantigen prediction. We close on the parts that an
EE-trained mind tends to under-weight: incidental findings and the
right not to know; the FDA Laboratory-Developed Test rule of 2024 and
the EU In-Vitro Diagnostic Regulation; the gaps in the Genetic
Information Nondiscrimination Act; the systematic ancestry bias in
clinical resources; and the Havasupai precedent that began the
modern conversation about indigenous data sovereignty.

The previous sixteen chapters teach you to build the pipeline. This
one teaches you what is required when somebody is going to read its
output and act.


== From Research to Clinic <sec:ch17-research-to-clinic>

The shift from research bioinformatics to clinical bioinformatics is
not a matter of polish. It is a structural change in five dimensions
that constrains every design decision downstream of it.

Accuracy is the first. A research pipeline that reports 99.8 % recall
on a benchmark is excellent. A clinical pipeline that misses
0.2 % of pathogenic variants has missed a patient — and, eventually,
several patients, because the failure mode is silent and repeatable.
Reproducibility is the second. Every classification call must be
re-derivable years later from the same input file, the same tool
versions, the same reference databases, the same hashes. Turnaround
time is the third. Research analyses run for weeks; clinical analyses
run in three to fourteen days for most indications and in hours for
the urgent ones. Regulation is the fourth. A research lab is bound by
the institutional review board and by good scientific practice; a
clinical lab is bound by CLIA in the United States, by ISO 15189
internationally, by IVDR in the European Union, and by national
accreditation bodies on top of that. Sign-off is the fifth. A research
paper has authors; a clinical report has a named physician signatory
whose license to practise is the thing on the line.

#figure(
  image("../../diagrams/lecture-17/01-research-vs-clinical.svg", width: 92%),
  caption: [
    Two pipelines with the same shape, different discipline. The
    research track reuses every clinical stage; the clinical track
    documents every step, pins every version, and routes the report
    through a named physician.
  ],
) <fig:research-vs-clinical>

The lab ecosystem in the United States is structured around three
acronyms. *CLIA*, the Clinical Laboratory Improvement Amendments of
1988, is federal regulation of every clinical laboratory. Any lab that
returns a health-related result on a human sample must be CLIA-certified,
and genomic tests are uniformly classified as high-complexity. *CAP
accreditation* from the College of American Pathologists is a
voluntary-but-expected peer assessment that goes beyond the CLIA
minimums; most respected clinical labs carry it. *NYS CLEP* is the
New York State Clinical Laboratory Evaluation Program, which sets
the strictest jurisdictional bar in the country — passing "NY approval"
on a new assay is a meaningful milestone. Internationally, ISO 15189
is the equivalent medical-laboratory quality standard, and UKAS in the
United Kingdom and DAkkS in Germany are the national accreditation
bodies that audit against it.

Bioinformatics inside a CLIA lab looks unfamiliar to anyone whose only
experience is academic research computing. Every pipeline release has a
written *validation document* — colloquially the "vali-doc" — that
demonstrates accuracy on a known benchmark sample (typically the GIAB
reference HG002, or matched tumour-normal cell lines for somatic
assays). Every release has a *standard operating procedure* with
explicit version tracking, and any change to a tool — even a patch
bump — triggers either re-validation or formal documentation justifying
its absence. Every batch carries internal quality controls (known
positive and negative samples that should call as expected) and the lab
participates in external quality assessment programmes (blind
proficiency tests such as GenQA) that catch systematic drift.

The assay regulatory classification adds a second axis. A *Laboratory-Developed
Test*, or LDT, is a test developed, validated, and performed within a
single clinical laboratory. Historically the FDA exercised "enforcement
discretion" — that is, declined to regulate LDTs — and so the bulk of
clinical genomic testing operated under CLIA alone. The 2024 final rule
ended that discretion; the implications run through the rest of this
chapter and through @sec:ch17-regulatory in particular. An
*In-Vitro Diagnostic device*, or IVD, is an FDA-cleared or FDA-approved
assay sold as a self-contained kit of reagents, instruments, and
software; the regulatory bar is much higher, and the FoundationOne CDx
and Illumina TSO500 panels are canonical examples. *Research Use Only*
reagents — RUO — are labelled explicitly for research; using them in a
clinical workflow is a regulatory violation unless the lab has
performed and documented its own clinical validation, which is the
move by which RUO inputs feed validated LDT outputs.

#figure(
  image("../../diagrams/lecture-17/11-regulatory-landscape.svg", width: 92%),
  caption: [
    The clinical-lab stack as three layers. Physical infrastructure
    below, certified lab in the middle, classified test output on top.
    The 2024 FDA LDT rule reaches into the amber-coloured LDT tier that
    used to operate under enforcement discretion.
  ],
) <fig:regulatory-stack>

#note[
  Most modern clinical genomic tests are LDTs built on a mixed
  infrastructure of FDA-cleared instruments (Illumina sequencers,
  IVD-marked library-preparation kits) and in-house analysis pipelines.
  This is the architecture the 2024 FDA rule was written for. The
  practical effect is that the analysis pipeline — the software a
  bioinformatician writes — sits at the inflection point where
  regulatory scrutiny is changing fastest.
]

Concretely, a clinical software engineer notices a handful of
differences from research practice. Production tooling runs one to
two years behind the research frontier; stability beats novelty.
Containerised execution is the norm, and the specific container image
hash is part of the validation. Every auxiliary file — the reference
genome, ClinVar, gnomAD, HGMD, the predictor models — is pinned by MD5.
Every run has a UUID, and every output file has a provenance chain back
to the raw reads. Pipelines run periodically on the GIAB reference
samples and a regression-test harness compares the output against the
last known-good run; an unexpected diff triggers a release block.
Failures are *fail-closed*: a variant caller that silently outputs an
empty VCF on a low-coverage region must throw a loud warning, not
emit a quietly empty file that downstream consumers treat as "no
findings."

#tip[
  The mental image for a clinical pipeline is closer to FDA-regulated
  medical-device firmware than to a research notebook. The engineering
  is similar in kind to embedded software in a 510(k)-cleared ECG, or
  firmware in an FCC-certified broadcast transmitter — every change
  requires re-validation; every deployment has a version-locked
  release; every failure mode is characterised in a documented
  limitations section. If you have worked in automotive, aviation,
  or medical-device EE, the culture will feel familiar.
]


== ACMG/AMP Variant Classification <sec:ch17-acmg>

Before 2015 clinical labs classified variants independently, and the
same variant could be called "pathogenic" at one lab and "variant of
uncertain significance" at another. The discordance was a problem for
patients moving between providers, for meta-analyses across labs, and
for multidisciplinary tumour boards trying to act on a common set of
calls. In 2015 the American College of Medical Genetics and Genomics
and the Association for Molecular Pathology jointly published a
twenty-eight-page consensus document — Richards, Aziz, Bale, and
colleagues, _Genetics in Medicine_ — that defined a standardised
classification framework. It encoded three things: a five-class output
schema, a set of evidence codes, and a rule table that combined the
evidence codes into a class. The document did not invent a new
classification system; it codified a common language across the
systems that GeneDx, Ambry, ARUP, Partners' Laboratory for Molecular
Medicine, Myriad, Invitae, and other major labs had each been running
internally. The ACMG/AMP guidelines have since become the operating
standard in almost every clinical genomics lab in the world.

The five-class output is the user-facing surface of the framework.
*Pathogenic* (P) means at least 99 % probability that the variant
causes the patient's disease; *likely pathogenic* (LP) is 90 to 99 %;
*variant of uncertain significance* (VUS) covers 10 to 90 %; *likely
benign* (LB) is 0.1 to 10 %; *benign* (B) is below 0.1 %. The
intervals are not approximations — they are the thresholds against
which the rule table is calibrated. Each class also carries a distinct
clinical action. P and LP findings are acted on: a patient with a
pathogenic BRCA1 variant enters enhanced screening, a pathogenic MLH1
carrier enters Lynch-syndrome surveillance, a pathogenic Long-QT
channelopathy variant changes anaesthesia decisions. VUS findings are
reported with a note and explicitly not acted on; reclassifying them as
evidence accumulates is a permanent background task across the field.
LB and B findings are usually not reported at all in a clinical
context.

#figure(
  image("../../diagrams/lecture-17/02-acmg-classes.svg", width: 92%),
  caption: [
    The five classes are five distinct decisions, not a slider with
    VUS as a middle band. The action attached to each class — act,
    wait, ignore — is the operational consequence of the probability
    interval.
  ],
) <fig:acmg-classes>

The VUS problem is structural. The majority of variants observed in
most genes are VUS at first encounter, and for rare genes with little
population-frequency data the VUS rate approaches 50 %. Reclassifying
VUS findings toward P/LP or LB/B as functional studies, family
segregation, and additional case reports accumulate is the central
activity of the variant-curation profession.

=== Evidence codes and their strengths

The ACMG/AMP framework defines twenty-eight evidence codes, organised
by direction (pathogenic versus benign) and by strength (very strong,
strong, moderate, supporting on the pathogenic side; stand-alone,
strong, supporting on the benign side). The pathogenic codes are PVS1
at the top, PS1 through PS4 in the strong tier, PM1 through PM6 in
the moderate tier, and PP1 through PP5 in the supporting tier. The
benign codes mirror them: BA1 stand-alone, BS1 through BS4 strong,
BP1 through BP7 supporting. Two of the original codes — PP5 and BP6,
which counted a "reputable source already classifies this" as
supporting evidence — were deprecated in the 2018 ClinGen update
because they introduced a feedback loop into ClinVar.

The most consequential code is PVS1: a *loss-of-function* variant in a
gene where loss-of-function is a known disease mechanism. Nonsense
mutations, frameshifts, canonical splice-site disruptions, and
single-exon deletions in genes like CFTR or BRCA1 or DMD apply this
code. The caveats are also important: a stop-gain in the very last
exon may escape nonsense-mediated decay, in which case it is not
truly loss-of-function and PVS1 does not apply; the ClinGen Sequence
Variant Interpretation working group has published explicit decision
trees for when PVS1 applies at full strength versus when it should be
downgraded to PVS1_strong or PVS1_moderate. The other strong codes
cover same-amino-acid changes that match a known pathogenic variant
(PS1), confirmed de novo events with parentage verified (PS2),
functional studies showing damaging effect (PS3), and case-control
enrichment (PS4). The moderate codes cover hotspot location (PM1),
absence from gnomAD (PM2), trans-configuration with another pathogenic
variant in a recessive disease (PM3), protein-length changes (PM4),
novel amino-acid changes at a position where a different change is
known pathogenic (PM5), and assumed-de-novo without parentage
confirmation (PM6). The supporting codes cover family segregation
(PP1), constrained-gene location (PP2), in-silico-predictor agreement
(PP3), and phenotype specificity (PP4).

#figure(
  image("../../diagrams/lecture-17/03-evidence-codes.svg", width: 95%),
  caption: [
    The twenty-eight-code reference. Strength gradients indicate which
    codes carry more weight in the rule combinations; the two deprecated
    codes (PP5 and BP6) are crossed out.
  ],
) <fig:evidence-codes>

The benign codes work symmetrically. BA1 is stand-alone: an allele
frequency at or above 5 % in a control population is enough by itself
to classify a variant as benign. BS1 covers a frequency inconsistent
with disease prevalence; BS2 observation in healthy adults at an allele
count inconsistent with disease; BS3 functional studies showing no
damaging effect; BS4 lack of segregation in affected family members.
The supporting benign codes cover missense variants in genes where
truncating variants rather than missense cause disease (BP1),
in-trans observation with another pathogenic variant in a
phenotypically healthy individual (BP2), in-frame indels in repetitive
regions (BP3), in-silico-predictor agreement that the variant is
benign (BP4), an alternative molecular diagnosis already in hand
(BP5), and silent variants at non-splice positions (BP7).

=== Rule combination

Once the evidence codes are assigned the *combinatorial rules* produce
the class. Pathogenic requires either one PVS1 paired with at least one
PS, PM, or PP at specific weighting; or at least two PS; or one PS plus
at least three PM; or one PS plus two PM plus two PP; or one PS plus
at least four PP; or at least three PM plus two PP. Likely pathogenic
relaxes each of those combinations by one strength tier. Benign is
either one BA1 alone or at least two BS codes. Likely benign is one BS
plus one BP, or at least two BP. Anything that does not satisfy a rule
or that has conflicting pathogenic and benign codes is a VUS by
default. The specific combinations form a lookup table that variant
scientists memorise and that automated tools — InterVar, VarSome,
Franklin — encode literally. The hard problem is not the lookup; it is
the upstream judgement that assigns the codes.

#note[
  The ACMG/AMP rule table is an *explicit, auditable rule-based
  classifier*. The choice to keep it rule-based rather than to swap in
  end-to-end machine learning was deliberate. Three reasons. First,
  auditability: every call traces back to specific evidence with a
  written rationale, which a neural network cannot provide.
  Second, regulatory acceptance: FDA and EU IVDR approval processes
  expect explainable outputs, particularly for software classified as
  a medical device under @sec:ch17-regulatory. Third, expert
  agreement: the framework encodes a generation of clinical-genetics
  consensus that a learned model would have to rediscover from much
  smaller labelled sets than the consensus draws on. The matched
  filter from Chapter 16 has a sibling here: the auditable
  finite-state classifier wins in the regulated domain even when a
  data-driven model could in principle achieve higher raw accuracy.
]

The Bayesian reformulation makes the structure clear. Tavtigian and
colleagues (_Genetics in Medicine_, 2018) showed that the 2015 rule
table is well approximated by a multiplicative likelihood-ratio model
in which each evidence code is assigned a likelihood ratio
$"LR"^+$: PVS1 corresponds to a ratio of about 350, a strong code
to about 19, a moderate code to about 4.3, a supporting code to about
2.1, and the benign codes are the mirror downweights. Combining
independent evidence is multiplication of likelihood ratios, and the
rule table encodes which combinations cross the 99 %, 90 %, 10 %, and
0.1 % posterior-probability thresholds that separate the five classes
under a prior of about 1.5 × 10⁻⁴ pathogenic on a randomly
chosen rare variant. The rule-based form exists because clinicians can
audit it; the Bayesian form exists because that is what the rule
table is mathematically equivalent to.

#figure(
  image("../figures/ch17/f1-acmg-bayesian-tavtigian.svg", width: 95%),
  caption: [
    Evidence codes as likelihood ratios under the Tavtigian
    reformulation. Each code's $"LR"^+$ value combines multiplicatively
    with the prior pathogenic odds; the four posterior thresholds at
    0.001, 0.1, 0.9, and 0.99 partition the result into the five
    classes. The rule table is what this calculus looks like once it is
    discretised into auditable categories.
  ],
) <fig:acmg-bayesian>

The EE analogue is *Neyman–Pearson multi-sensor detection*. Independent
sensors with known likelihood ratios combine into a posterior over the
hypothesis of interest, and decision thresholds partition the posterior
into actions. The variant scientist is operating exactly this kind of
detector — the sensors are the evidence codes, the likelihood ratios
are the Tavtigian numbers, and the decisions are P/LP/VUS/LB/B. The
only special feature is that the inputs are not produced by hardware
but by an upstream chain of curation: gnomAD frequencies,
ClinVar lookups, predictor scores, literature reviews. The next
section is about that chain.

=== Gene-specific refinements

The 2015 framework is deliberately gene-agnostic, but disease biology
is not. For high-interest genes the ClinGen consortium has published
*gene-specific refinements* that supersede the generic rules. BRCA1
and BRCA2 have refined PVS1 logic that depends on exon location, and
founder-population frequency thresholds that diverge from the generic
PM2/BA1 cutoffs. TP53 has hotspot definitions that change which
positions qualify for PM1. The mismatch-repair genes — MLH1, MSH2,
MSH6, PMS2 — have specific thresholds and additional evidence codes
calibrated to Lynch syndrome. RASopathy genes have modifications for
de novo handling. The cardiac channelopathy genes SCN5A and KCNQ1
carry modifications that reflect their unusually high VUS burden. The
rule when curating a variant is to apply the gene-specific refinement
if one exists and to fall back to generic ACMG/AMP otherwise.

The Invitae-developed *Sherloc* framework (Nykamp et al., 2017)
extends ACMG/AMP with finer evidence-code subcategories and explicit
quantitative thresholds — points per code, summed to a score — designed
for high-throughput reclassification. Most Invitae and Ambry reports
are derived from Sherloc internally but communicate the result in
ACMG/AMP five-class language. The framework is a refinement, not a
replacement.

=== ClinVar as community truth

*ClinVar*, hosted by the NIH, is the community's aggregate
variant-classification database. As of 2024 it contains roughly 2.5
million unique variants with multiple submissions each from clinical
labs, research groups, and expert panels. Submissions are classified
under ACMG/AMP or Sherloc and ClinVar aggregates them into a
star-rated consensus: a single-star submission is one lab's opinion;
a two-star variant has multiple concordant submissions; three- and
four-star ratings require expert-panel review and, for the highest
tier, a clinical practice guideline that cites the variant.

#figure(
  image("../../diagrams/lecture-17/04-clinvar.svg", width: 92%),
  caption: [
    ClinVar at scale. Submissions grow log-linearly since the
    database's 2013 launch. The plurality of unique variants remain
    VUS — the ongoing reclassification work of the field.
  ],
) <fig:clinvar>

Clinical variant scientists treat ClinVar as a starting point rather
than as ground truth. The first questions on any new variant are
whether anybody has classified it before, whether the existing
classifications are concordant, and — critically — whether the
classification has changed over time. A variant that was likely
pathogenic in 2018 and was reclassified to VUS in 2022 is a different
object from a variant that has been LP since first submission, and the
history matters more than the current label. Single-star submissions
sometimes disagree with four-star expert panels; the expert panel is
usually but not always right. The skill is in calibrating that
"usually."


== The Variant Evidence Ecosystem <sec:ch17-evidence>

A clinical variant classification is built from a set of standardised
inputs that have themselves taken twenty years to mature. The four
that matter most are population frequency, in-silico functional
prediction, splice prediction, and prior classification. Each of them
has a lineage, a current state, and a set of known failure modes; each
of them is consumed by the ACMG/AMP rule table as inputs to particular
evidence codes.

=== Population frequency

The single most-used resource is *gnomAD*, the Genome Aggregation
Database led from the Broad Institute. The 2020 release (Karczewski et
al., _Nature_) tabulated allele frequencies across roughly 125,748
exomes and 15,708 whole genomes from unrelated individuals stratified
by eight continental ancestry groups. The 2024 v4 release expanded the
exome cohort to 730,947 individuals and the genome cohort to 76,215.
gnomAD is the descendant of the Exome Variant Server and the Exome
Aggregation Consortium (ExAC, 2016), themselves descendants of the
1000 Genomes Project. The lineage matters because frequency-based
evidence codes are calibrated to the resource's coverage: PM2 and BA1
thresholds shifted slightly between ExAC and gnomAD v2, and again
between v2 and v4, and any clinical pipeline must pin the version it
queries.

The codes consume gnomAD as follows. BA1 fires when allele frequency
exceeds 5 %: stand-alone benign. BS1 fires when the frequency exceeds
the gene's disease-prevalence-derived threshold: strong benign. PM2
fires when the variant is absent or extremely rare — the precise
threshold is gene-specific in the modern guidelines but is roughly
1 / 200,000 alleles for fully-penetrant autosomal-dominant diseases
and 1 / 20,000 for recessive diseases. The ClinGen panels have
published gene-by-gene thresholds for BRCA1 and BRCA2, the
mismatch-repair genes, the long-QT genes, and several others.

#figure(
  image("../../diagrams/lecture-17/05-gnomad-thresholds.svg", width: 92%),
  caption: [
    The gnomAD frequency axis. PM2 lives at the rare end with
    gene-dependent cutoffs; BA1 lives at the common end at 5 %. The
    grey band in between contributes no frequency-based evidence on
    its own.
  ],
) <fig:gnomad-thresholds>

Ancestry stratification is essential. A variant that is rare in
European cohorts but common in African cohorts is still benign for the
African patient; using European-only frequencies is a known source of
mis-classification. gnomAD reports per-ancestry-group frequencies for
this reason, and the BA1 / BS1 / PM2 evaluations should use the
patient's ancestry-group frequency rather than the global
"popmax" frequency, except in the specific case where popmax exceeds
the BA1 threshold (in which case the variant is benign regardless).

=== In-silico functional impact

For missense variants the rule table allows in-silico functional
predictors to contribute the supporting code PP3 if the predictors
agree on damaging effect, or BP4 if they agree on benign effect. The
codes are deliberately supporting: a single predictor score is never
sufficient by itself.

The lineage runs back to the early 2000s. *SIFT* (Ng and Henikoff,
2001) classified missense variants by sequence conservation alone, on
the assumption that positions strongly conserved across species are
intolerant of change. *PolyPhen* and its successor *PolyPhen-2*
(Adzhubei et al., 2010) added structural features alongside
conservation. Both became staples of clinical pipelines and remain
embedded inside more recent ensemble predictors. *CADD* (Kircher et
al., 2014) combined dozens of conservation, regulatory, and structural
features in a logistic regression trained on the difference between
derived and ancestral alleles, producing a phred-scaled score from
roughly 0 to 40. *REVEL* (Ioannidis et al., 2016) ensembled thirteen
older predictors — SIFT, PolyPhen-2, MutationTaster, GERP, and others —
in a random forest trained on rare missense variants and produced a
0-to-1 score with ACMG/AMP-relevant thresholds of 0.7 or above for PP3
and 0.15 or below for BP4. *AlphaMissense* (Cheng et al., _Science_
2023) was the first predictor trained on the AlphaFold-2 structural
representations from Chapter 15; it predicts likely-pathogenic,
ambiguous, or likely-benign for every possible amino-acid substitution
in the human proteome and currently sits at or near the state of the
art on most benchmarks.

#figure(
  image("../../diagrams/lecture-17/06-predictors.svg", width: 95%),
  caption: [
    Three predictors that contribute to the PP3 / BP4 supporting
    evidence line. Different inputs, different architectures, the same
    role: each is a sensor with a known likelihood ratio that
    independently votes the variant pathogenic or benign.
  ],
) <fig:predictors>

The honest reading is that these predictors are *decision supports*,
not deciders. PP3 is supporting evidence; it cannot produce a
pathogenic classification on its own. A common misuse pattern — seen
in both automated variant-classification tools and hastily reviewed
pipelines — is to treat the predictor score as directly actionable:
"AlphaMissense says pathogenic, report as LP." This is wrong. A single
supporting-tier code cannot produce LP under the ACMG/AMP rule table,
and the predictors are known to over-call pathogenic for rare benign
variants in under-represented populations.

#warn[
  Predictor agreement matters. Using a single predictor is discouraged
  by the framework; the ClinGen specifications generally require two
  or three independent predictors to concur before PP3 or BP4 fires.
  When the predictors disagree, the evidence is not "weakly
  pathogenic" — it is *no evidence at all* on this line. Reading a
  disagreement as a weighted average is the single most common
  software-side mistake in the ACMG/AMP space.
]

=== Splicing

Splice-altering variants are an under-appreciated pathogenicity class.
A variant fifty bases inside an intron can create a new splice site,
skip an exon, and produce a truncated protein — and conventional
predictors trained on protein-level effects miss them entirely.
*SpliceAI* (Jaganathan et al., _Cell_, 2019) is the workhorse. It
takes a ten-kilobase window of genomic sequence centred on the
variant and outputs four per-position probabilities: acceptor gain,
acceptor loss, donor gain, donor loss. A maximum score across the
window above 0.5 is conventionally treated as supporting evidence of a
splice effect, contributing PP3 or, in canonical splice-site
disruptions, the strong-tier PVS1_strong. A score above 0.8 is
treated as strong evidence (PS3_moderate) and routes the variant into
splicing-functional follow-up.

The architecture is a deep convolutional network with a receptive
field that spans both intron and exon, trained on roughly 20,000
canonical splice sites from GENCODE. SpliceAI was the first
deep-learning predictor to reach clinical utility on splice prediction
and remains the standard; subsequent models (Pangolin, SpliceTransformer)
sharpen specific regions but have not displaced it.

=== Literature mining

Literature contributes evidence to several codes: PS3 (functional
studies showing damaging effect), PS4 (case-control enrichment), PP1
(family segregation), PP4 (phenotype specificity). Finding the
relevant literature is increasingly automated. *LitVar*, hosted by
NCBI, links variants to PubMed articles automatically using normalised
variant identifiers and saves a variant scientist hours of hand
searching. *Mastermind*, from Genomenon, is a commercial system with
broader coverage and subscription pricing. ClinVar submissions
themselves frequently include literature citations and serve as a
back-door index into the relevant papers.

=== A working day

A clinical variant-scientist workflow for one variant typically runs
through nine steps. Query gnomAD for ancestry-stratified frequencies;
query ClinVar for prior classifications and submission history; query
functional predictors for REVEL, AlphaMissense, and SpliceAI scores;
search the literature through LitVar or Mastermind; query phenotype
databases (OMIM, Orphanet) for the gene-disease association; assemble
the evidence codes; apply the ACMG/AMP rule table; write a
classification rationale in plain text; enter the classification into
the lab's LIMS and generate the patient-facing report. The total time
per variant is thirty minutes for the simple cases — a clearly
benign common variant, a clearly pathogenic loss-of-function variant
in a well-characterised gene — and four hours for the difficult ones.
A clinical variant scientist handles between ten and thirty variants
per day.

#figure(
  image("../figures/ch17/f2-evidence-dossier.svg", width: 95%),
  caption: [
    A worked evidence dossier for a candidate BRCA1 missense variant.
    The five rows are the five standard resources; each fires zero,
    one, or several evidence codes; the rule table at the bottom
    consumes them and produces the class. This is the data structure
    that lives inside a clinical LIMS for every variant the lab calls.
  ],
) <fig:evidence-dossier>

The profession that lives in this dossier is not one an AI replaces.
ML-based predictors inform individual evidence codes — PP3 and BP4,
in particular — but the synthesis across resources, the judgement
about which gene-specific refinement applies, the framing of the
written rationale, and the legal sign-off all require the human
variant scientist. The role is M.Sc. or Ph.D. level, two to three
years of certification-track supervised practice, and a board
examination from ACMG or its international equivalents. The
ML-augmented future is one in which the variant scientist's
productivity roughly doubles as automated tools pre-fill evidence
codes and literature synthesis; the human still signs the report.


== Pharmacogenomics and Immunogenomics <sec:ch17-pgx>

Pharmacogenomics is the cleanest example of clinical genomics already
operating at scale as a lookup table on a state machine. The patient
is measured (genotype across a handful of cytochrome-P450 and related
loci); the measurement is mapped to a discrete phenotype (poor,
intermediate, normal, rapid, ultra-rapid metaboliser); the phenotype
indexes a fixed dosing recommendation in the CPIC guideline tables.
No ML, no learning, no inference at runtime — because the
evidence-generation machinery has already done the optimisation and
the table is the compiled answer.

The underlying biology is that the cytochrome P450 enzymes — CYP2D6,
CYP2C9, CYP2C19, CYP3A4, and a dozen others — and several other
pharmacokinetic enzymes carry common loss-of-function or
gain-of-function variants at population frequencies that differ
sharply by ancestry. About 20 drug–gene pairs have CPIC Level A
evidence — "use this test result to change dosing" — and another 40
or so sit at Level B with actionable recommendations that are less
universally adopted. Roughly half of FDA-approved drugs carry some
pharmacogenomic label, though only a fraction are routinely
genotyped in clinical practice.

=== Four canonical cases

*CYP2D6 and codeine* is the textbook example. Codeine is a prodrug;
CYP2D6 metabolises it to morphine, which provides analgesia. Poor
metabolisers — about 7 % of European-ancestry populations — get no
pain relief from codeine. Ultra-rapid metabolisers, who are about 3 %
of Europeans and up to 30 % of some North African populations,
metabolise codeine so fast that morphine levels spike and respiratory
depression follows. A small number of paediatric deaths from codeine
in ultra-rapid metabolisers led the FDA to restrict codeine in
children. The pharmacogenomic action is to know the patient's CYP2D6
metaboliser status before prescribing and to choose an alternative
opioid (oxycodone via CYP3A4, hydromorphone) for poor and
ultra-rapid metabolisers.

*Warfarin* (VKORC1 + CYP2C9) is the second canonical case. Warfarin
is the classical anticoagulant, and its dosing window is famously
narrow — too little and the patient clots, too much and the patient
bleeds. The c.-1639G>A variant in VKORC1 sets the enzyme-sensitivity
locus; the \*2 and \*3 alleles in CYP2C9 set the metabolism locus.
Genotype-guided dosing algorithms (Gage et al., 2008; the
International Warfarin Pharmacogenetics Consortium, 2009) combine the
two with age, weight, and INR target into a predicted steady-state
dose, reducing time-to-stable-INR by roughly 15 %. Adoption has been
limited, partly because the direct-acting oral anticoagulants
(rivaroxaban, apixaban) have displaced warfarin for most indications.

*Abacavir* (HLA-B\*57:01) is the third. Abacavir is an HIV
reverse-transcriptase inhibitor; patients carrying the HLA-B\*57:01
allele have approximately a 50 % rate of severe hypersensitivity
reaction on first exposure. Without typing, roughly 5 % of patients
experience a life-threatening reaction. Pre-prescription HLA-B\*57:01
genotyping is now standard of care in every HIV clinic globally, and
has reduced abacavir hypersensitivity incidence by more than 95 %.
The model — single allele, single drug, near-perfect screening
performance — is what the rest of pharmacogenomics aspires to.

*Thiopurines* (TPMT, NUDT15) is the fourth. Thiopurines —
azathioprine, mercaptopurine — are used for acute lymphoblastic
leukaemia, autoimmune disease, and inflammatory bowel disease. TPMT
metabolises them; TPMT deficiency (homozygous \*3A, \*3B, or \*3C)
causes severe myelosuppression at standard doses. The Asian-prevalent
NUDT15 variants have a similar phenotype on a parallel locus.
Pre-treatment TPMT typing in European-ancestry patients and NUDT15
typing in Asian-ancestry patients reduces severe toxicity by an order
of magnitude.

#figure(
  image("../../diagrams/lecture-17/07-pgx-examples.svg", width: 95%),
  caption: [
    Four canonical pharmacogenomic pairs. Drug to gene to variants to
    metaboliser phenotype to clinical action. Once the table is built
    the runtime is one lookup.
  ],
) <fig:pgx-examples>

=== Star alleles and metaboliser phenotypes

Pharmacogenomics uses *star-allele* nomenclature to describe
haplotypes. CYP2D6\*1 is the reference (wild-type) allele; CYP2D6\*4 is
a common loss-of-function haplotype that carries a splicing variant
plus several SNPs; CYP2D6\*10 is a reduced-function haplotype common
in East Asian populations; CYP2D6\*17 is a reduced-function haplotype
common in African populations. Each star allele is a defined
combination of variants — a phased haplotype — and the canonical
definitions are maintained at the *PharmVar* consortium.

Metaboliser phenotypes derive from the patient's two star alleles by
an activity-score sum. Each star allele is assigned an activity score
on a fixed scale (commonly 0, 0.25, 0.5, 1, or 2 depending on the
gene); summing the two yields a score that maps to a phenotype band.
CYP2D6 \*1/\*1 is normal metaboliser (score 2); \*1/\*4 is intermediate
(score 1); \*4/\*4 is poor (score 0); \*1/\*2 ×N where N is a gene
duplication is ultra-rapid (score 3 or above). The star-allele
abstraction hides the underlying variants from the clinician — the
clinician sees only "which stars does this patient carry?" and the
software pipeline does the rest.

#tip[
  When implementing a star-allele caller, two failure modes recur.
  First, novel haplotypes that do not match any defined star allele
  must be flagged rather than silently mapped to the closest known
  allele — the closest match may have a different activity score.
  Second, CYP2D6 in particular has well-known structural variants
  including whole-gene deletions and tandem duplications that
  short-read sequencing cannot reliably resolve; a CNV-aware caller or
  a long-read assay is required for full coverage.
]

=== PharmGKB and CPIC

*PharmGKB* is the community-maintained pharmacogenomics knowledge
base. For each drug–gene pair it lists the known variants and their
effects, the metaboliser phenotype definitions, the evidence levels,
and the links to dosing guidelines. *CPIC*, the Clinical
Pharmacogenetics Implementation Consortium, drafts and publishes the
dosing guidelines themselves in a standardised format: "if patient
genotype is X, give dose Y." The typical clinical workflow is to
order a pharmacogenomic test (a targeted Illumina panel, a microarray,
or whole-genome sequencing in some pre-emptive programs), let
software translate the genotype into star alleles and into a
metaboliser phenotype, present the clinician with a dosing
recommendation pulled from CPIC, and let the clinician adjust the
prescription. Several health systems — Vanderbilt PREDICT, St. Jude
PG4KDS, Mayo Clinic — pre-emptively genotype patients on enrolment so
that the PGx results are available for any future prescribing
decision.

#figure(
  image("../figures/ch17/f3-pgx-decision-flow.svg", width: 95%),
  caption: [
    The pharmacogenomic decision flow as a discrete state machine.
    Patient genotype to star alleles to activity score to metaboliser
    phenotype to CPIC dosing table to clinical action. Each arrow is a
    table lookup. The compiled-table form is what makes this
    deployable; the science was building the tables.
  ],
) <fig:pgx-flow>

The EE framing is *discrete-state control policy lookup*. Measure the
patient state (genotype), map it to a phenotype, look up the action in
a fixed table, apply. This is the simplest possible closed-loop
control: one measurement, a discrete state space, table lookup. The
evidence-generation machinery — clinical trials, mechanistic studies,
expert review — has already done the optimisation, and the CPIC table
is the compiled answer. The EE analogue is a mode-based controller
where each mode has pre-optimised parameters; the runtime's only job
is to identify which mode applies. It is also why pharmacogenomics is
relatively easy to productionise once the tables exist and impossible
for drugs where the tables do not yet exist — the whole game is in
building the lookup table.

Despite the structure, the field is under-utilised. Reimbursement is
patchy; many physicians do not see the genotype when prescribing
because EHR integration is uneven; the actionable list is growing but
each addition requires plumbing. The bottleneck is sociotechnical, not
technical.

=== Immunogenomics: HLA, TCR, BCR

Pharmacogenomics treats HLA-B\*57:01 / abacavir as a single drug–gene
pair. The HLA system is much larger and increasingly clinically
relevant on its own. *Immunogenomics* is the bioinformatics of three
highly polymorphic gene families — HLA, T-cell receptors, and B-cell
receptors / immunoglobulins — that together encode the human adaptive
immune system.

*HLA typing* is the first piece. The Major Histocompatibility Complex
on chromosome 6 contains the HLA genes — the most polymorphic locus in
the human genome, with tens of thousands of catalogued alleles. HLA
class I (A, B, C) presents intracellular peptides to CD8 T cells;
HLA class II (DR, DQ, DP) presents extracellular peptides to CD4 T
cells. HLA typing matters for transplant matching (mismatched grafts
trigger rejection), for drug hypersensitivity beyond
abacavir (HLA-B\*15:02 / carbamazepine, HLA-B\*58:01 / allopurinol),
for autoimmune-disease association (HLA-DRB1\*04 with rheumatoid
arthritis, HLA-DQ2/DQ8 with celiac, HLA-B\*27 with ankylosing
spondylitis), and for cancer immunotherapy (HLA loss-of-heterozygosity
in tumours is a known immune-evasion mechanism). Computational HLA
typing from short-read whole-exome or whole-genome sequencing is
non-trivial because the locus is too divergent for standard alignment;
specialised tools — *OptiType*, *HLA-LA*, *arcasHLA* — graph-align or
assemble reads against the curated IMGT/HLA database.

*T-cell receptor* and *B-cell receptor* repertoires arise from V(D)J
recombination — somatic rearrangement of variable, diversity, and
joining gene segments with random nucleotide insertion at junctions.
Each lymphocyte gets a unique receptor; the repertoire in a healthy
adult holds roughly $10^8$ distinct clones. Bulk repertoire sequencing
amplifies the TCR β chain (or both chains) from sorted T cells and
counts clones; single-cell paired-chain protocols (10× Genomics VDJ)
link the α/β pair to the gene-expression profile of the same cell.
Tools — *MiXCR*, *IgBLAST*, *TRUST4* — reconstruct CDR3 sequences,
V/D/J segment usage, and abundance distributions. Diversity metrics
(Shannon, Simpson, D50), clonality, and the detection of "public"
clones (identical CDR3 across individuals, typically pathogen-specific)
are the standard analytic outputs. Targeted clone tracking — finding
a tumour-specific TCR clone in peripheral blood as a minimal-residual-disease
marker — is now standard practice in several haematology indications.

*Neoantigen prediction* is the bridge from immunogenomics to cancer.
Tumour somatic mutations create novel peptides not present in normal
tissue; if those peptides bind the patient's HLA, they may be
presented to T cells and become targets for personalised immunotherapy.
The pipeline runs tumour-normal sequencing to identify somatic missense
mutations, HLA-types the patient's germline DNA, translates each
mutation into a window of 9-to-11-mer peptides, predicts HLA binding
affinity with *NetMHCpan* or *MHCflurry*, filters for tumour
expression (the peptide must be present on an expressed transcript) and
strong binding (rank below 0.5 % is the convention), and outputs a
ranked candidate list. The pipeline drives personal neoantigen
vaccines — Moderna's mRNA-4157 / Merck's V940 in melanoma, currently
in Phase III — TCR-T cell therapy, and pembrolizumab-response
biomarkers (high neoantigen load correlates with response).

#figure(
  image("../../diagrams/lecture-17/13-immunogenomics.svg", width: 95%),
  caption: [
    Immunogenomics in three layers. Top: HLA peptide presentation
    cascade from somatic mutation through proteasomal cleavage and
    HLA loading to T-cell recognition. Middle: V(D)J recombination of
    the TCR β locus producing the $approx 10^15$-clone theoretical
    repertoire. Bottom: the five-stage neoantigen pipeline cascading
    tumour mutations through HLA binding to clinically actionable
    candidates.
  ],
) <fig:immunogenomics>

The EE framing here is *generative-model inference under a structured
channel*. The V(D)J process is a Markov-like cascade with
conditional probabilities for V → D → J segment usage and
Bernoulli/geometric models for nucleotide insertion. Tools like
*OLGA* and *SONIA* (Quigley, Walczak, and colleagues, 2019–2021) fit
explicit generative models to the repertoire and compute per-clone
generation probabilities — the equivalent of a likelihood under the
V(D)J channel. Public clones are then identified as outliers in the
generation-probability-by-abundance plane: high abundance with
surprisingly low generation probability suggests selection. The
techniques transfer almost directly from speech-coding and
channel-modelling work.


== Secondary Findings and the Right Not to Know <sec:ch17-secondary>

A patient comes in for cardiomyopathy testing. The lab sequences their
exome. In addition to the indication-relevant variant the exome
contains a pathogenic BRCA1 mutation — unrelated to the test
indication but with serious downstream implications. Should the lab
report it? Should the lab have looked? Should the patient have been
told what the lab would look at? These are the *secondary findings*
questions, and they arise naturally whenever broad sequencing is
performed.

The ACMG's position has evolved. In 2013 the original recommendations
proposed an opt-out model: labs were to report a minimum gene list
unless the patient declined. Pushback was immediate, and by 2014 the
recommendation had shifted to opt-in. The 2021 update added explicit
patient-facing consent specifically acknowledging secondary findings,
and the 2024 SF v3.2 list is the operating standard now — about 81
genes selected for inclusion against a single criterion: *actionability*.
A finding is on the SF list because there is a known medical
intervention — surveillance, prophylactic surgery, lifestyle change,
drug avoidance — that changes outcome when the variant is known. The
list is dominated by hereditary cancer syndromes (BRCA1/2,
MLH1/MSH2/MSH6/PMS2/EPCAM for Lynch syndrome, TP53 for Li-Fraumeni,
APC, MUTYH, NF1/2, PTEN, RET, VHL), inherited arrhythmia syndromes
(SCN5A, KCNQ1, KCNH2, RYR2), cardiomyopathies (MYH7, MYBPC3, TNNT2),
familial hypercholesterolemia (LDLR, APOB, PCSK9), malignant
hyperthermia susceptibility (RYR1, CACNA1S), and a handful of
miscellaneous actionable conditions.

#figure(
  image("../../diagrams/lecture-17/08-acmg-sf.svg", width: 92%),
  caption: [
    The ACMG SF v3 list at a glance. Cancer and cardiac genes dominate;
    the inclusion criterion is actionability — a known intervention
    that changes outcome when the variant is known.
  ],
) <fig:acmg-sf>

The right-not-to-know is real and contested. The arguments for
reporting are simple: actionable findings save lives; patients benefit
from information about their own genome; the incremental cost is
near-zero once the sequencing data exist. The arguments against
universal return are also concrete. Psychological distress in the
absence of clear intervention is well documented for some genes
(Huntington's disease testing has decades of experience here).
Downstream costs — surveillance imaging, prophylactic surgery — may
not be reimbursed by insurance. Family members are implicated by the
patient's result without their consent. And not all "actionable"
findings clearly benefit patients; the evidence is genuinely mixed
for some genes, particularly the lower-penetrance ones.

The operational workflow, when an incidental finding is suspected, is
a three-gate decision. First, the variant must be classified as P or
LP under ACMG/AMP — not VUS, however suggestive. Second, the gene
must be on the current ACMG SF list (or otherwise clearly actionable
under documented lab policy). Third, the patient must have consented
to receive secondary findings. All three gates must pass for the
finding to be reported, and the variant is then re-confirmed (ideally
with a different chemistry), routed through a genetic counsellor
before the physician, and documented in a distinct "Secondary
Findings" section of the clinical report.

#figure(
  image("../figures/ch17/f4-sf-three-gate.svg", width: 95%),
  caption: [
    The secondary-findings three-gate decision. Gene on SF list AND
    variant classified P/LP AND patient consented; any failure blocks
    the report. The structure is what makes the decision defensible in
    audit because every path can be traced to specific data and
    documents.
  ],
) <fig:sf-gate>

#warn[
  The failure modes are bad in both directions. A VUS in a
  cancer-predisposition gene reported as actionable can drive an
  unnecessary prophylactic mastectomy. A pathogenic SF-list variant
  missed because it was not queried can cost a patient years of
  preventable surveillance. An SF-list variant reported without
  consent is both an ethical failure and a legal liability. The
  three-gate structure exists because each failure mode has happened.
]

The discussion is also age-stratified. A 65-year-old patient with an
adult-onset movement disorder indication who carries a pathogenic
BRCA1 variant will probably benefit from learning it — they are still
in the surveillance-effective window. An 85-year-old in the same
situation may not. A 15-year-old with a rare seizure disorder
typically does not learn adult-onset BRCA1 findings because the
intervention is not recommended in childhood and the patient cannot
yet make their own decision about their adult-onset risks. The
field's emerging consensus is to defer most adult-onset
secondary-finding disclosures in paediatric patients until the
patient reaches an age where they can consent for themselves.


== The Regulatory Landscape <sec:ch17-regulatory>

The regulatory shell around clinical genomics has been tightening for
a decade and reached an inflection point in 2024. The shell has three
relevant components: the FDA pathway in the United States, the
European IVDR, and the direct-to-consumer carve-out.

The FDA clears or approves sequencing assays through two main
pathways. *510(k) clearance* certifies that the assay is "substantially
equivalent" to a previously cleared predicate device; the process is
faster and less rigorous and is used for moderate-risk devices.
*Premarket Approval* (PMA) is the rigorous pathway for high-risk
Class III devices and is the route for companion-diagnostic
approvals — assays that determine eligibility for a specific
therapy. The major FDA-cleared and approved genomic assays include
FoundationOne CDx (Foundation Medicine / Roche, PMA 2017), a
324-gene tumour-profiling NGS panel approved as a companion
diagnostic for multiple targeted therapies; MSK-IMPACT (Memorial
Sloan Kettering, authorization 2017), a 468-gene tumour panel that
unusually received FDA authorization as an LDT and laid the
groundwork for the 2024 rule change; Illumina TruSight Oncology 500;
and ThermoFisher's Oncomine Dx Target Test for lung cancer. For
germline testing FDA-cleared kits are far fewer; much of the
clinical market is LDTs.

=== The 2024 FDA LDT rule

Historically the FDA claimed regulatory authority over LDTs but
exercised *enforcement discretion* — in effect, not regulating them.
This is how most clinical genomic testing in the United States
operated for two decades. In April 2024 the FDA finalised a rule
ending that discretion. The phased implementation runs four years:
adverse-event reporting requirements in year one, quality-system
regulation compliance in year two, registration and listing
requirements in year three, and premarket review for high-risk LDTs
in year four. Industry response is mixed; large clinical labs (Mayo
Clinic, Quest, LabCorp) are preparing actively, while smaller
academic-affiliated labs are concerned about the documentation
burden. A Congressional review of the rule is ongoing.

#figure(
  image("../../diagrams/lecture-17/09-fda-pathway.svg", width: 92%),
  caption: [
    The regulatory map for a sequencing assay. RUO is research-only;
    LDTs are tightening under the 2024 rule; the IVD branch splits
    between 510(k) substantial-equivalence and PMA full approval.
  ],
) <fig:fda-pathway>

The practical effect for bioinformatics software is that clinical
pipelines will need more explicit regulatory documentation. The
existing validation practices in @sec:ch17-research-to-clinic will
need to formalise further: predetermined change-control plans,
documented design history files, software-of-unknown-provenance
analyses for upstream dependencies, cybersecurity threat modelling
for cloud deployments. The FDA's *Software as a Medical Device*
framework now applies to most clinical pipelines, and the standards
that govern medical-device firmware engineering — IEC 62304 for
software lifecycle, IEC 62366-1 for usability engineering, ISO 14971
for risk management — increasingly govern bioinformatics shipped into
clinical use.

=== EU IVDR

The European Union's *In-Vitro Diagnostic Medical Devices Regulation*
(IVDR) replaced the previous IVD Directive in May 2022, after multiple
delays. It classifies IVDs into four risk classes (A through D), with
genomic tests for severe hereditary disease typically Class C and
tumour-profiling assays Class C or D. The requirements above the old
Directive include notified-body certification for Class C and D
devices, clinical-performance studies demonstrating accuracy on the
intended population, active post-market surveillance, and the CE mark
for market access. Notified-body capacity has been the
bottleneck — many European labs that operated under the previous
Directive have struggled to complete IVDR certification on time, and
the EU has repeatedly extended deadlines. As of 2025 IVDR is a
significant operational pain point for European clinical labs.

=== Direct-to-consumer testing

DTC genomics is regulated differently from clinical testing because
the test result returns to the consumer without a physician
interpreter. *23andMe*, founded 2007, began returning health-related
reports without FDA involvement and received a warning letter in
November 2013 demanding that the health reporting cease. The company
complied and shifted to ancestry-only reports. In April 2017 the FDA
authorised 23andMe to return limited health-related information — a
tightly scoped set of variants rather than general clinical
interpretation. In March 2018 the BRCA scope was expanded but only
for the three Ashkenazi-founder mutations, not general BRCA testing.
The current regulatory equilibrium permits the ancestry-plus-light-health
model with explicit disclaimers and optional genetic counselling, but
the equilibrium remains contested. Should a consumer directly learn
they carry a BRCA variant without a physician interpreter? The
compromise treats the question as a balance between autonomy and the
risk that an isolated DTC result without context drives an
inappropriate clinical decision.


== Ethics, Equity, and Data Sovereignty <sec:ch17-ethics>

The technical content of clinical genomics is one half of the chapter;
the surrounding ethics and policy is the other. Three threads matter
most: the gaps in legal protection against genetic discrimination, the
systematic ancestry bias in clinical resources, and the rights of
communities whose biological samples have driven much of the field's
research.

=== GINA and its limits

The *Genetic Information Nondiscrimination Act* (GINA) was signed
into US federal law in May 2008 after a thirteen-year legislative
history. It passed the Senate 95–0 and the House 414–1 — extremely
rare bipartisan unanimity. It covers two domains. *Employment*:
employers cannot use genetic information in hiring, firing, or
promotion decisions. *Health insurance*: group and individual health
insurers cannot use genetic information for coverage or pricing
decisions. The bipartisan support came after the 2003 Human Genome
Project completion made the discrimination potential concrete enough
to overcome industry resistance, and the carve-outs were left in as
deliberate compromises.

The carve-outs are large. GINA does not cover *life insurance*: a
life insurer can ask about and use genetic test results. It does not
cover *long-term care insurance* or *disability insurance*. It does
not cover *small employers* (fewer than 15 employees). It does not
cover the *US military*, which operates under a separate regulatory
regime. The practical consequence is that a BRCA1-positive
patient may be unable to obtain affordable life insurance, may be
declined long-term-care coverage, and is exposed to disability-insurance
underwriting that uses the result. A handful of states (California,
Florida, Vermont) have stricter state-level protections that fill
parts of the gap, but the federal coverage is narrower than most
patients assume.

#figure(
  image("../../diagrams/lecture-17/10-gina-coverage.svg", width: 92%),
  caption: [
    GINA's reach. Health insurance and large-employer hiring are
    covered. Life, disability, long-term care, small employers, and
    the military sit outside; a handful of states fill the gap
    partially.
  ],
) <fig:gina>

Bills to close the carve-outs have been introduced periodically and
have not passed; the insurance-industry resistance that GINA's
original compromise reflected has not weakened.

=== Ancestry bias in clinical databases

Clinical-genomics resources are dominated by people of European
ancestry, and the bias cascades through every downstream calculation.
ClinVar: roughly 70 % of submissions concern variants seen in European
patients. gnomAD v4: 56 % European overall, with the African-ancestry
subset growing but still smaller. The functional predictors (REVEL,
AlphaMissense, CADD) are trained on labelled sets where European
variants dominate. Polygenic risk scores — Chapter 13 — are
calibrated on overwhelmingly European training cohorts. The GWAS
Catalog tracked diversity directly for over a decade: about 89 % of
its participants were of European ancestry as of mid-2020s estimates
(Sirugo, Williams, and Tishkoff, _Cell_ 2019; updated periodically).

The consequence for non-European patients is concrete and measurable.
*Higher VUS rates*: variants common in African populations are more
often "rare and unknown" in European-dominated ClinVar, and they get
flagged as candidate pathogens that the lab cannot confidently
classify. *Lower predictor accuracy*: missense predictors trained
predominantly on European labelled data underperform on African
variants — measured F1 gaps of roughly 0.05 to 0.10 across recent
benchmarks. *Worse clinical actionability*: fewer firm classifications
translates directly into fewer actionable findings and worse care.

#figure(
  image("../../diagrams/lecture-17/12-ancestry-bias.svg", width: 92%),
  caption: [
    Three layers of bias. ClinVar submissions are more European than
    gnomAD coverage; predictor F1 follows the same ordering.
    Closing the F1 gap is a data problem, not an algorithm tweak.
  ],
) <fig:ancestry-bias>

The EE reading is that ancestry bias is the training-distribution
problem from Chapter 16 in its starkest form. A classifier — the
combination of ACMG/AMP rules plus predictors plus ClinVar lookups —
trained predominantly on European data does not generalise to
African patients. The algorithm is working as designed on the data
it has; the fix is to change the data. Initiatives to close the gap
are now structural infrastructure. *All of Us* (NIH, 2018+) targets
one million-plus US participants with an explicit minority
representation mandate and enrolment at roughly 700,000 as of 2024.
*Our Future Health* (UK, 2020+) targets five million UK participants.
*H3Africa* (2012+) builds African genomic research infrastructure.
*TOPMed* (NHLBI, 2014+) is the multi-ethnic WGS cohort. *PAGE*
(Population Architecture using Genomics and Epidemiology, 2008+) is
the US minority-focused consortium. *Million Veteran Program* (VA,
2011+) recruits a more diverse US cohort than UK Biobank
equivalents. *Biobank Japan*, *Biobank Korea*, and *Taiwan Biobank*
extend the East Asian coverage. The data these cohorts generate flows
back into ClinVar, gnomAD, AlphaMissense training sets, and PRS
portability models. Progress is measurable — gnomAD v4 in 2024 is
materially more diverse than v2 in 2020 — but full equity remains
years away.

#warn[
  Any clinical pipeline you deploy today has systematic accuracy
  differences by ancestry. This is not hypothetical — it is
  measurable and published. When your pipeline signs out a VUS for a
  patient of African ancestry, the VUS rate is higher in absolute
  terms than it would be for a matched European patient. Reporting
  ancestry bias as a known limitation in the test's standard
  documentation is honest engineering. Failing to acknowledge it is
  failing the patient.
]

=== Havasupai and data sovereignty

The *Havasupai Tribe v. Arizona Board of Regents* case (filed 2004,
settled 2010) is the defining precedent for indigenous data
sovereignty in US research. In 1989, Arizona State University
researchers collected blood samples from the Havasupai Tribe — a small
community living in the Grand Canyon region — for a diabetes-genetics
study. Over the following decade the samples were used for additional
research on schizophrenia, inbreeding, and population migration — none
of which the Tribe had consented to, and several of which the Tribe
considered culturally harmful. In 2004 tribal members learned about
the additional research. The Tribe sued; the case settled in 2010 with
ASU returning the blood samples (a significant symbolic act), paying
\$700,000, and building a health clinic.

The case became a defining precedent for *indigenous data
sovereignty*. NIH strengthened tribal-consultation requirements; the
All of Us program includes explicit tribal governance provisions in
its consent framework; the broader research community has internalised
the principle that consent for one study does not authorise
repurposing for another.

The *CARE Principles for Indigenous Data Governance*, articulated by
the Global Indigenous Data Alliance in 2019, codify the position.
*Collective benefit* — data ecosystems facilitate positive outcomes
for Indigenous Peoples. *Authority to control* — Indigenous Peoples'
rights to determine how their data are used. *Responsibility* — those
working with Indigenous data have a duty to share benefits and
acknowledge the relationship. *Ethics* — Indigenous Peoples' rights
and wellbeing are the primary concern. CARE is complementary to the
*FAIR Principles* (Findable, Accessible, Interoperable, Reusable)
from the open-data movement; FAIR addresses technical interoperability
and CARE adds the governance layer. Modern genomic studies involving
Indigenous communities increasingly incorporate both, with
community-level consent, restricted data access, and benefit-sharing
agreements as standard practice.

The broader point generalises beyond Indigenous communities. Genomic
data is uniquely identifying, uniquely informative about kin, and
uniquely difficult to anonymise. The standards developed for
community-level consent in the Havasupai aftermath are progressively
becoming relevant for other vulnerable cohorts — undocumented
populations, refugee groups, communities affected by historic
research abuse. Engineering practices that lock data sharing to the
specific consent scope, that flag repurposing attempts, and that
enable revocation of consent when feasible are the technical
correlate of this policy shift.


== Summary <sec:ch17-summary>

- Clinical bioinformatics is regulated software, not research
  scientific computing. CLIA in the United States, ISO 15189
  internationally, EU IVDR, and FDA Software-as-a-Medical-Device
  frameworks each impose validation, version-pinning,
  change-control, and physician sign-off requirements that shape
  every design decision downstream.
- The ACMG/AMP 2015 framework is the global classification standard.
  Twenty-eight evidence codes combine via a rule table into five
  classes (P, LP, VUS, LB, B). The table is well approximated by a
  Bayesian likelihood-ratio model (Tavtigian 2018); the rule-based
  form exists for auditability.
- The evidence ecosystem is structured. gnomAD supplies frequency;
  REVEL, AlphaMissense, and CADD supply missense-impact prediction;
  SpliceAI supplies splice prediction; ClinVar supplies prior
  classification. Predictors are decision supports for PP3/BP4, never
  deciders.
- Pharmacogenomics already operates as a discrete-state control
  policy. CYP2D6/codeine, VKORC1+CYP2C9/warfarin, HLA-B\*57:01/abacavir,
  TPMT/thiopurines are the canonical cases. PharmGKB and PharmVar
  provide the reference tables; CPIC publishes the dosing guidelines.
- Immunogenomics extends pharmacogenomics into HLA typing, V(D)J
  repertoire analysis, and neoantigen prediction. The neoantigen
  pipeline (tumour-normal sequencing, HLA typing, NetMHCpan binding
  prediction, expression filtering) drives modern personalised
  immunotherapy.
- Secondary findings require an explicit three-gate decision: gene on
  ACMG SF list (v3.2, about 81 genes), variant classified P/LP,
  patient consented. Right-not-to-know is real and is age-stratified
  in practice.
- Regulation tightened sharply in 2024. The FDA LDT rule ends
  enforcement discretion in a phased four-year implementation; EU
  IVDR is fully in force with ongoing notified-body bottlenecks; DTC
  testing operates under a narrowly scoped FDA exception.
- GINA protections are narrower than most assume. Health insurance
  and large-employer hiring are covered; life, disability, long-term
  care, small employers, and the military are not. Patient counsel
  should know the gap.
- Ancestry bias in clinical resources is systematic and measurable.
  ClinVar, gnomAD, and predictor accuracy all skew toward European
  ancestry. The fix is diverse cohorts (All of Us, Our Future Health,
  H3Africa, TOPMed, Million Veteran Program, Biobank Japan/Korea),
  not algorithm tweaks.
- Indigenous data sovereignty is a settled requirement after the
  Havasupai case. CARE Principles complement FAIR; consent is
  community-level, repurposing requires re-consent, and benefit
  sharing is a structural expectation.


== Exercises <sec:ch17-exercises>

#strong[1.] #emph[Bayesian reading of the rule table.] Take Tavtigian's
likelihood-ratio assignments — PVS1 about 350, PS about 19, PM
about 4.3, PP about 2.1, with the benign codes as mirror downweights —
and a prior pathogenic probability of about 1.5 × 10⁻⁴
for a randomly chosen rare variant. Compute the posterior
pathogenic probability for each of the following code combinations and
classify it under the five-class thresholds (0.001, 0.1, 0.9, 0.99):
(a) PVS1 alone; (b) PS3 + PM2; (c) PM1 + PM2 + PP3; (d) BS1 + BP4;
(e) PVS1 + BS1 (conflicting evidence). Verify that your results match
the published rule table.

#strong[2.] #emph[gnomAD threshold check.] Pick five published
pathogenic variants from ClinVar for genes you find interesting (one
from BRCA1, one from CFTR, one from MLH1, one from MYBPC3, one from
LDLR). Query gnomAD v4 for each variant's overall and per-ancestry-group
allele frequencies. For each variant, evaluate whether BA1 fires,
whether BS1 fires under the disease-specific prevalence threshold, and
whether PM2 fires. Report any disagreements with the published
classification and propose an explanation.

#strong[3.] #emph[Predictor agreement.] For 100 missense variants in a
gene of your choice that have known ACMG/AMP classifications in
ClinVar, pull REVEL, AlphaMissense, and CADD scores. Define
"agreement" as all three predictors pointing in the same direction at
their standard thresholds (REVEL above 0.7 or below 0.15,
AlphaMissense pathogenic/benign label, CADD scaled above 25). Tabulate
agreement rate, classification rate among agreed-predictors variants,
and the disagreement breakdown. Compare to the published agreement
rates in the AlphaMissense paper.

#strong[4.] #emph[Star-allele caller.] Implement a minimal CYP2D6
star-allele caller. Input: a sample's genotype at a curated set of
diagnostic positions for \*1, \*2, \*4, \*5 (full-gene deletion), \*10,
\*17, \*2xN (duplication). Output: the diplotype (e.g., \*4/\*10) and the
predicted metaboliser phenotype (PM, IM, NM, RM, UM) using activity
scores from the published CYP2D6 table. Run your caller on a synthetic
dataset of 1,000 simulated diplotypes; compare to PharmCAT or a
similar reference implementation.

#strong[5.] #emph[Secondary-findings audit.] You are tasked with
reviewing 500 historical exome reports from a clinical lab for
secondary-finding compliance. Define a precise audit protocol: which
variants in which genes should have been flagged; which gates (gene
on SF list, classification P/LP, consent on file) should be checked;
what counts as a finding-disclosure error versus a finding-detection
error. Write the protocol as a checklist that a junior auditor could
follow.

#strong[6.] #emph[GINA gap analysis.] A 35-year-old patient has just
received a positive BRCA1 result from a clinical lab. She is shopping
for life insurance, long-term-care coverage, and is interviewing for
an engineering position at a 12-employee startup. For each of the
three contexts, write a one-paragraph summary of the relevant federal
and state-level protections (or absence thereof) that bear on her
situation. Indicate which advice you would route to a genetic
counsellor and which to an attorney.

#strong[7.] #emph[Ancestry-bias quantification.] Pick a published
missense predictor (REVEL or AlphaMissense). Find a published
benchmark with per-ancestry-group performance. Compute the F1 gap
between European and African groups. Then design (in pseudocode) an
ablation: how would you train a predictor with equal weight per
ancestry group, and what dataset would you need? Estimate the size of
the labelled-data shortfall in the under-represented group.

#strong[8.] #emph[(Open-ended.)] Pick one published clinical-genomics
controversy from the past five years — a famous reclassification, a
guideline change, a regulatory action, a high-profile data-sovereignty
case. Identify the technical and policy layers it touches across this
chapter (which evidence codes were involved, which guidelines applied,
which regulatory bodies acted, which ethical principles were at
stake). Write a two-page case study with citations. Propose one
concrete engineering or policy change that would have changed the
outcome.


== Further Reading <sec:ch17-further-reading>

- *Richards, S., Aziz, N., Bale, S., et al.* (2015). "Standards and
  guidelines for the interpretation of sequence variants: a joint
  consensus recommendation of the American College of Medical Genetics
  and Genomics and the Association for Molecular Pathology."
  _Genetics in Medicine_ 17: 405 — 424. The original ACMG/AMP
  framework. Read in full at least once; the twenty-eight evidence
  codes are catalogued there with examples.
- *Tavtigian, S. V., Greenblatt, M. S., Harrison, S. M., et al.*
  (2018). "Modeling the ACMG/AMP variant classification guidelines as
  a Bayesian classification framework." _Genetics in Medicine_ 20:
  1054 — 1060. The likelihood-ratio reformulation of the rule table.
  The mathematics behind the cells in the lookup.
- *Karczewski, K. J., Francioli, L. C., Tiao, G., et al.* (2020).
  "The mutational constraint spectrum quantified from variation in
  141,456 humans." _Nature_ 581: 434 — 443. The gnomAD v2 paper. Pair
  with the v4 data release (2024) for the current resource.
- *Cheng, J., Novati, G., Pan, J., et al.* (2023). "Accurate
  proteome-wide missense variant effect prediction with
  AlphaMissense." _Science_ 381: eadg7492. The current state-of-the-art
  in missense prediction; the AlphaFold-derived structural prior
  carries this paper.
- *Jaganathan, K., Panagiotopoulou, S. K., McRae, J. F., et al.*
  (2019). "Predicting splicing from primary sequence with deep
  learning." _Cell_ 176: 535 — 548. The SpliceAI paper. The
  receptive-field discussion is exemplary EE-style writing about a
  biology problem.
- *Miller, D. T., Lee, K., Chung, W. K., et al.* (2022). "ACMG SF
  v3.1 list for reporting of secondary findings in clinical exome and
  genome sequencing." _Genetics in Medicine_ 24: 1407 — 1414. Read
  alongside the SF v3.2 update for the current list and inclusion
  criteria.
- *Sirugo, G., Williams, S. M., & Tishkoff, S. A.* (2019). "The
  missing diversity in human genetic studies." _Cell_ 177: 26 — 31.
  The canonical reference on ancestry bias. The 89 % European figure
  for the GWAS Catalog originates here.
- *Mello, M. M., & Wolf, L. E.* (2010). "The Havasupai Indian Tribe
  case — lessons for research involving stored biologic samples."
  _New England Journal of Medicine_ 363: 204 — 207. Short, clear
  account of the case and its implications for consent practice.
- *FDA Final Rule on Laboratory-Developed Tests* (2024). 89
  _Federal Register_ 37286. The 2024 LDT rule itself. Read at least
  the preamble; the four-year phased implementation schedule is the
  operative table.
