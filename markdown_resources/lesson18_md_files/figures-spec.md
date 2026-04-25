# Lecture 18 — Figures Specification

> **Scope**: Static diagrams for Lecture 18 (Cancer Genomics: Integrated Capstone).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-18.md`.

---

## 0. Conventions for This Lecture

- Figures are custom SVG; heavy on workflow / architecture diagrams, some quantitative panels (VAF distributions, signatures, CCF plots).
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Cancer-cell-vs-normal-cell colouring: tumour `--warning` red; normal `--accent-soft` cobalt; subclones progressively lighter tints.
- Mutation substitution colours: C>A cobalt, C>G amber, C>T red, T>A grey-dark, T>C teal, T>G violet (COSMIC convention).
- Tier colouring (AMP/ASCO/CAP): Tier I green, Tier II amber-light, Tier III grey, Tier IV muted.
- Biomarker colours: TMB cobalt, MSI amber, HRD teal.
- Typography: Inter for UI labels; JetBrains Mono for tool names, variant HGVS, gene names, percentages, VAFs, activity scores, COSMIC signature IDs (SBS1, SBS2, ...).
- Escape `&`, `<`, `>` as XML entities.

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Hallmarks of cancer | Part 1 | Custom SVG |
| 2 | Driver vs passenger + clonal evolution | Part 1 | Custom SVG |
| 3 | Tumour-normal, FFPE, and ctDNA designs | Part 2 | Custom SVG |
| 4 | Mutect2 workflow | Part 3 | Custom SVG |
| 5 | COSMIC mutational signatures | Part 3 | Custom SVG |
| 6 | Cancer SV detection — short vs long reads | Part 4 | Custom SVG |
| 7 | Gene fusion visualisation | Part 4 | Custom SVG |
| 8 | Subclonal reconstruction and evolutionary tree | Part 5 | Custom SVG |
| 9 | AMP/ASCO/CAP four-tier system | Part 6 | Custom SVG |
| 10 | TMB / MSI / HRD biomarker trio | Part 6 | Custom SVG |
| 11 | Precision oncology therapy decision tree | Part 7 | Custom SVG |
| 12 | Industry ecosystem for cancer genomics | Wrap | Custom SVG |

---

## Figure 1 — Hallmarks of cancer

**File**: `diagrams/lecture-18/01-hallmarks.svg`
**Lecture anchor**: §1.2 Hallmarks of cancer
**ViewBox**: `0 0 1200 720`

### Purpose

The Hanahan–Weinberg "hallmarks of cancer" schematic updated to the 2022 revision. Shows ~14 hallmark capabilities as spokes, annotated with example molecular alterations.

### Content

**Central cell glyph**. Stylised cancer cell with disorganised chromatin.

**10–14 labelled spokes radiating outward**, each a rounded rectangle:

- Sustaining proliferative signalling — e.g. activating RAS / RAF / MEK mutations.
- Evading growth suppressors — e.g. TP53, RB1 loss.
- Resisting cell death — e.g. BCL2 overexpression.
- Enabling replicative immortality — e.g. TERT promoter mutations.
- Inducing angiogenesis — e.g. VEGF overexpression.
- Activating invasion and metastasis — e.g. EMT programs.
- Reprogramming energy metabolism — e.g. HIF1α upregulation (Warburg effect).
- Evading immune destruction — e.g. PD-L1 upregulation.
- Genome instability — e.g. BRCA1/2 loss, MMR loss.
- Tumour-promoting inflammation — e.g. NF-κB activation.
- Phenotypic plasticity (2022) — e.g. EMT / dedifferentiation.
- Senescent cells (2022) — e.g. SASP signalling.
- Non-mutational epigenetic reprogramming (2022) — e.g. chromatin remodelling.
- Polymorphic microbiomes (2022) — e.g. oral / gut microbiome shifts.

**Colour coding**. Original (2000) hallmarks: cobalt. 2011 additions: amber. 2022 additions: teal.

**Top-left small inset**. Timeline: "Hanahan & Weinberg 2000 → 2011 update → 2022 revision".

### Style notes

- Radial layout, spokes evenly spaced.
- Central cell glyph stylised (not photorealistic).
- Each spoke label in Inter with its "example alteration" in JetBrains Mono italic.

---

## Figure 2 — Driver vs passenger + clonal evolution

**File**: `diagrams/lecture-18/02-drivers-clonal.svg`
**Lecture anchor**: §1.3 Driver vs passenger; §1.4 Clonal evolution
**ViewBox**: `0 0 1200 640`

### Purpose

Visualise the driver-vs-passenger distinction and show a clonal-evolution tree of subclones with annotated drivers.

### Content

**Left panel — driver vs passenger** (top half of figure).

A horizontal strip representing a tumour's somatic mutations: ~10,000 tick marks in grey ("passenger"); ~10 tick marks in red ("driver"), labelled with genes (TP53, KRAS G12V, EGFR T790M, etc.).

Below the strip: "Most mutations are passengers (neutral). A handful drive tumour fitness and are therapeutic targets."

**Right panel — clonal evolution tree** (bottom half of figure).

A phylogenetic tree with:

- Root ("MRCA / trunk") with clonal mutations labelled (TP53 inactivation, CDKN2A loss).
- First branch ("clone 1") with CCF 100% for trunk events, plus subclonal drivers (KRAS G12V at CCF 60%).
- Second branch ("clone 2") with the same trunk events and a different subclonal driver (EGFR T790M at CCF 25%).
- Leaves labelled "dominant clone" and "minor clone".

**Sidebar annotation** (narrow vertical box at the right):

- "Trunk: 100% CCF; targetable in all cells."
- "Branches: < 100% CCF; resistance-risk subclones."
- "Single-time-point bulk sequencing: VAF distribution reveals this."

### Style notes

- Driver ticks: red; passenger ticks: grey.
- Tree nodes: coloured by clone.
- Branch labels: gene names in monospace.

---

## Figure 3 — Tumour-normal, FFPE, and ctDNA designs

**File**: `diagrams/lecture-18/03-sequencing-designs.svg`
**Lecture anchor**: §2.1 Tumour-only vs tumour-normal; §2.2 FFPE; §2.3 ctDNA
**ViewBox**: `0 0 1200 560`

### Purpose

Three side-by-side workflow diagrams covering the three canonical cancer sequencing designs.

### Content

**Panel 1 — Tumour-normal paired**. Top: tumour block icon → FASTQ_T. Below: blood draw icon → FASTQ_N. Both → alignment → matched variant caller (Mutect2, Strelka2) → somatic VCF. Annotation: "gold standard; clean germline / somatic separation".

**Panel 2 — FFPE tumour-only**. FFPE block icon (with "10 years old" label) → damaged / fragmented DNA → library prep (with "UMIs + enrichment" tag) → tumour-only somatic caller + aggressive germline filtering (gnomAD frequency) → somatic VCF. Annotation: "cheaper, clinical-friendly; sacrifices ~5–10% specificity".

**Panel 3 — Liquid biopsy / ctDNA**. Blood tube icon → cell-free DNA → UMI-barcoded library → deep sequencing (30,000×) → low-VAF variant calling (sub-0.1% VAF) → minimal residual disease flag or tumour-genotype fingerprint. Annotation: "non-invasive; serial monitoring; needs UMI-based error correction".

**Bottom caption band**. "Design choice driven by: sample availability, budget, clinical question, turnaround time."

### Style notes

- Icons (block / blood tube / FASTQ / VCF) rendered as simple geometric glyphs.
- Each panel with its own subtle colour tint.

---

## Figure 4 — Mutect2 workflow with cancer-specific filters

**File**: `diagrams/lecture-18/04-mutect2.svg`
**Lecture anchor**: §3.2 Mutect2 architecture revisited
**ViewBox**: `0 0 1200 720`

### Purpose

Expand Mutect2's architecture with explicit cancer-specific filter stages.

### Content

**Top — inputs**. Three rounded rectangles: Tumour BAM, Normal BAM, Panel-of-Normals + Reference.

**Middle — local reassembly + joint somatic call**. Rectangle: "Local reassembly → haplotype likelihoods". Outputs "raw somatic VCF".

**Middle — cascade of cancer-specific filters** (vertical stack of 5 rectangles):

1. Panel-of-Normals (PoN) filter — removes recurrent artefacts.
2. Orientation-bias filter — removes FFPE deamination.
3. Contamination filter (CalculateContamination) — accounts for cross-sample mixing.
4. Tumour-in-normal filter — rejects any variant > 5% in normal.
5. Low-VAF artefact filter — subclonal sanity.

Arrows connecting them in sequence.

**Bottom — output**. Final filtered somatic VCF with PASS/FAIL annotations.

**Right-side annotation box**. "Key differences from germline HaplotypeCaller: (a) per-read likelihood over a continuous allele-fraction axis, not just {0, 0.5, 1}; (b) PoN filter; (c) orientation-bias filter; (d) contamination estimation."

### Style notes

- Filter cascade: vertical stack with each filter as a rectangle.
- Inputs: left side; outputs: bottom.
- Arrow flow: top-to-bottom.

---

## Figure 5 — COSMIC mutational signatures

**File**: `diagrams/lecture-18/05-cosmic-signatures.svg`
**Lecture anchor**: §3.4 Mutational signatures
**ViewBox**: `0 0 1200 640`

### Purpose

Canonical per-signature 96-bar spectrum view for a handful of clinically-relevant signatures.

### Content

**Grid of 4–6 signature panels**, each showing a 96-bar histogram with the standard COSMIC layout.

**Panel 1 — SBS1 (clock-like, 5-methylcytosine deamination)**:

- Dominant peaks in C>T at CpG (NCG contexts).
- Colour: red (C>T).
- Label: "SBS1 — spontaneous 5mC deamination; clock-like; present in all tumours".

**Panel 2 — SBS4 (tobacco smoking)**:

- Dominant peaks in C>A at specific contexts.
- Colour: cobalt (C>A).
- Label: "SBS4 — tobacco smoking; C>A at specific flanking contexts".

**Panel 3 — SBS7a (UV light)**:

- Dominant peaks in C>T at CC and TC dinucleotides (UV-dimer signatures).
- Colour: red (C>T).
- Label: "SBS7a — UV exposure; C>T at CC/TC dipyrimidines".

**Panel 4 — SBS3 (HRD / BRCA-deficient)**:

- Flat / diffuse spectrum across all 96 contexts.
- Colour: mixed (uniform).
- Label: "SBS3 — HRD; PARP inhibitor sensitive".

**Panel 5 — SBS13 (APOBEC)**:

- Dominant peaks in C>G at TCW contexts.
- Colour: amber (C>G).
- Label: "SBS13 — APOBEC deamination; common in breast, bladder cancers".

**Panel 6 — SBS6 (MMR deficiency)**:

- Peaks at multiple substitution types; bias toward microsatellite repeat contexts.
- Label: "SBS6 — MMR deficiency; MSI-H tumours; pembrolizumab-responsive".

**Bottom annotation**. "A tumour's 96-bin mutation spectrum ≈ Σ (exposure × signature). NMF decomposition recovers which signatures contributed and how much."

### Style notes

- Each signature panel: standard 96-bar layout with substitution-type colour bands.
- Compact grid to fit all 6 panels.
- Labels below each panel.

---

## Figure 6 — Cancer SV detection — short vs long reads

**File**: `diagrams/lecture-18/06-sv-detection.svg`
**Lecture anchor**: §4.2 SV detection from short vs long reads
**ViewBox**: `0 0 1200 560`

### Purpose

Compare short-read vs long-read resolution of a complex cancer SV (chromothripsis example).

### Content

**Top panel — reference + true SV topology**. A reference chromosome bar with 5 breakpoints indicated. The true rearrangement: "segments A, B, C, D, E shuffled into a scrambled order D-B-E-A-C".

**Middle panel — short-read evidence**. Discordant paired-end reads mapped as arcs connecting distant positions. Split reads rendered as pairs of half-arrows. Annotation: "5 breakpoints called with difficulty; 2 breakpoints missed due to repeats". Per-breakpoint confidence low (amber).

**Bottom panel — long-read evidence**. A handful of 20 kb PacBio HiFi reads that individually span multiple breakpoints. Each breakpoint directly resolved. Annotation: "all 5 breakpoints clearly resolved; structure reconstructed unambiguously". Per-breakpoint confidence high (green).

**Right-side annotation box**. "Short-read SV calling: Manta, GRIDSS — ~60–80% sensitivity on cancer SVs. Long-read: Sniffles2, Severus — ~95%+. Complex events (chromothripsis, ecDNA) favour long reads."

### Style notes

- Arcs for discordant-pair evidence.
- Long reads as horizontal bars with breakpoint tick marks.
- Colour per breakpoint confidence level.

---

## Figure 7 — Gene fusion visualisation

**File**: `diagrams/lecture-18/07-gene-fusion.svg`
**Lecture anchor**: §4.3 Gene fusions from RNA-seq
**ViewBox**: `0 0 1200 560`

### Purpose

Anatomy of a gene fusion from two parent genes to a functional fusion protein, using EML4-ALK as the example.

### Content

**Top — Chromosomal context**. Two different chromosomes drawn side-by-side:

- Chromosome 2 with EML4 gene (exons colour-coded green).
- Chromosome 2 (not different, actually!) — no wait, EML4-ALK is a paracentric inversion on chr 2.
- Actually two genes on chromosome 2p: EML4 and ALK. Inversion joins EML4 exons 1–6 with ALK exons 20–29.

Arrows show the inversion / rearrangement event.

**Middle — fusion transcript**. A spliced mRNA showing EML4 exons 1–6 (green) joined to ALK exons 20–29 (cobalt) with a single junction.

**Bottom — fusion protein**. Domain-level diagram: EML4's TAPE domain + coiled-coil domain + ALK's tyrosine-kinase domain. Annotation: "coiled-coil forces dimerisation → constitutively active ALK kinase → oncogenic signalling".

**Right-side clinical context**. "EML4-ALK in ~5% of lung adenocarcinomas. Treatable with crizotinib, alectinib, lorlatinib. Resistance mutations reported in ALK kinase domain → second-generation inhibitors."

### Style notes

- Exons colour-coded per parent gene.
- Domain boxes labelled with monospace.
- Breakpoint line emphasised.

---

## Figure 8 — Subclonal reconstruction and evolutionary tree

**File**: `diagrams/lecture-18/08-subclonal-tree.svg`
**Lecture anchor**: §5.1 Subclonal architecture; §5.2 Phylogenetic reconstruction
**ViewBox**: `0 0 1200 640`

### Purpose

End-to-end view of subclonal reconstruction: VAF distribution → CCF clustering → phylogenetic tree.

### Content

**Left panel — VAF distribution**. Histogram of ~500 somatic variants. Multiple modes visible: a dominant cluster near 0.30 (clonal in a 60%-pure tumour), a smaller cluster near 0.15, a tail near 0.05. Colour-coded modes.

**Middle panel — PyClone clustering**. Same 500 points now placed in CCF space (recomputed from VAF with purity + ploidy). Three distinct clusters labelled "Clone A (CCF 100%)", "Clone B (CCF 50%)", "Clone C (CCF 20%)".

**Right panel — Phylogenetic tree**. Standard rooted tree:

- Root: "MRCA / trunk" (clonal mutations listed — TP53 loss, CDKN2A loss).
- First split: clone A continues + clone B branches off.
- Second split (inside clone A's lineage): clone C branches off.
- Leaves labelled with driver gene per branch (KRAS G12V in B; EGFR T790M in C).

**Bottom caption band**. "Single-sample bulk reconstruction is possible but CCF resolution is limited; multi-sample (multi-region, longitudinal) sharpens recovery."

### Style notes

- VAF histogram: colour per cluster.
- CCF plot: same colour coding.
- Tree branches colour-matched to clusters.

---

## Figure 9 — AMP/ASCO/CAP four-tier system

**File**: `diagrams/lecture-18/09-amp-asco-cap.svg`
**Lecture anchor**: §6.1 AMP/ASCO/CAP tier system
**ViewBox**: `0 0 1200 560`

### Purpose

Four stacked tiers with canonical examples and clinical-actionability annotations.

### Content

**Four horizontal tiers stacked top-to-bottom**, widest = highest actionability:

- **Tier I — Variants of strong clinical significance** (green).
  - Level A: FDA-approved therapy in this cancer type.
  - Level B: Professional-guidelines-recommended.
  - Example: EGFR L858R in NSCLC → erlotinib.
- **Tier II — Variants of potential clinical significance** (amber-light).
  - Level C: FDA-approved therapy in a different cancer type.
  - Level D: Clinical trials; preclinical strong evidence.
  - Example: BRAF V600E in a non-melanoma tumour → off-label vemurafenib trial.
- **Tier III — Variants of unknown clinical significance** (grey).
  - In a cancer gene but no therapeutic implication.
  - Example: PIK3CA rare missense outside a known hotspot.
- **Tier IV — Benign / passenger variants** (muted).
  - Not clinically relevant.
  - Example: common SNP in gnomAD at > 5%.

**Right-side comparison box**. "ACMG/AMP (germline, L17) classifies by pathogenicity probability. AMP/ASCO/CAP (somatic, L18) classifies by clinical action. Different questions → different frameworks."

### Style notes

- Tier bands stacked vertically with their canonical colours.
- Example variants in JetBrains Mono.

---

## Figure 10 — TMB / MSI / HRD biomarker trio

**File**: `diagrams/lecture-18/10-biomarker-trio.svg`
**Lecture anchor**: §6.3–6.5
**ViewBox**: `0 0 1200 560`

### Purpose

Three side-by-side mini-dashboards, one per biomarker, showing typical distributions, clinical thresholds, and therapy pairings.

### Content

**Panel 1 — TMB**:

- Histogram of TMB per tumour (mutations per Mb), log-scale x-axis 0.1 – 100.
- Vertical threshold line at TMB = 10 mut/Mb.
- Annotation arrow: "High TMB → pembrolizumab (tumour-agnostic)".
- Note: "thresholds are assay-specific".

**Panel 2 — MSI**:

- Bar chart of MSI-H vs MSI-L vs MSS frequency across a cohort.
- MSI-H cluster highlighted.
- Annotation arrow: "MSI-H / dMMR → pembrolizumab tumour-agnostic".
- Small detail: "IHC (MMR protein loss) confirms the sequencing call".

**Panel 3 — HRD**:

- Violin plot of HRD scores across ovarian cancers.
- HRD-positive threshold line annotated (score ≥ 42 for one assay).
- Annotation arrow: "HRD-positive → PARP inhibitor (olaparib, niraparib, rucaparib)".
- Note: "SBS3 signature + LOH + TAI + LST composite".

**Bottom caption band**. "Three distinct biomarker classes; each FDA-linked to specific therapy; each has assay-specific thresholds — validate per lab."

### Style notes

- Panels consistent width, separated by narrow dividers.
- Colour per biomarker (TMB cobalt, MSI amber, HRD teal).

---

## Figure 11 — Precision oncology therapy decision tree

**File**: `diagrams/lecture-18/11-therapy-tree.svg`
**Lecture anchor**: §7.1 Targeted therapy decisions
**ViewBox**: `0 0 1200 640`

### Purpose

Decision tree walking from a somatic variant report to a therapy recommendation, as used in a tumour board.

### Content

**Root node**: "Somatic variant report received at tumour board".

**First decision**: "Tier I variant present?"

- Yes → "FDA-approved targeted drug" (green leaf; examples: EGFR L858R → osimertinib; BRAF V600E + melanoma → dabrafenib + trametinib).
- No → proceed to next decision.

**Second decision**: "Tier II variant in a known actionable gene?"

- Yes → "Clinical trial enrolment" (amber leaf; e.g. basket trials for rare actionable variants).
- No → next.

**Third decision**: "TMB ≥ 10 mut/Mb or MSI-H?"

- Yes → "Immune checkpoint inhibitor (pembrolizumab)" (teal leaf).
- No → next.

**Fourth decision**: "HRD-positive?"

- Yes → "PARP inhibitor (olaparib)" (teal leaf).
- No → "Standard-of-care chemotherapy / exploratory trial" (grey leaf).

**Right-side annotation box**. "In practice: tumour boards weigh all factors jointly (performance status, prior therapies, co-morbidities). Decision trees are simplifications, not algorithms."

### Style notes

- Decision nodes as diamonds; outcomes as rounded rectangles.
- Edges labelled "yes / no".
- Leaves colour-coded by therapy class.

---

## Figure 12 — Industry ecosystem for cancer genomics

**File**: `diagrams/lecture-18/12-industry-ecosystem.svg`
**Lecture anchor**: Wrap-up, industry view
**ViewBox**: `0 0 1200 720`

### Purpose

Map of industry sectors and roles for graduates entering cancer genomics.

### Content

**Four-quadrant layout**:

- **Clinical labs** (top-left; cobalt). Examples: Foundation Medicine, Guardant, Tempus, Natera, Invitae Oncology. Typical roles: clinical bioinformatics pipeline engineer, variant scientist, regulatory lead.
- **Pharma / biotech** (top-right; amber). Examples: Roche, Novartis, AstraZeneca, Pfizer, Moderna, Illumina Clinical Genomics, mid-size biotechs. Typical roles: computational biologist, biomarker scientist, ML engineer for drug discovery, clinical-trial genomics lead.
- **Sequencing / platform** (bottom-left; teal). Examples: Illumina, PacBio, ONT, 10× Genomics, Element Biosciences, Ultima Genomics. Typical roles: bioinformatics R&amp;D, platform pipeline engineer, field applications scientist.
- **Academic cancer centres** (bottom-right; violet). Examples: MSK, Dana-Farber, MD Anderson, Broad Institute, Wellcome Sanger, St Jude, Fred Hutch. Typical roles: translational research scientist, core-facility lead, academic staff.

**Centre of the diagram**: "The cancer-genomics graduate's career paths — often cycling between these sectors over a career".

**Bottom annotation band**. "Hybrid roles (clinical + research + ML) are increasingly the norm. Skills that transfer: pipeline engineering, statistical genetics, ML for genomics, regulatory literacy, clinical-report writing."

### Style notes

- Quadrant boxes with coloured borders.
- Example-company names in Inter.
- Role descriptions in JetBrains Mono italic.
- Central hub connects all four quadrants with thin arrows.
