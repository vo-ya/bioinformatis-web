# Lecture 27 (proposed L11) — Figures Specification

> **Scope**: Static diagrams for Lecture 27 (Mass-Spectrometry Proteomics Primer).
> **Companion files**: `lecture-27.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Spectrum colours: MS1 cobalt, MS2 amber, fragment ions per b/y type (b cobalt, y red).
- Acquisition-mode colours: DDA amber, DIA cobalt, targeted (PRM) green.
- Quantification colours: LFQ grey, SILAC cobalt+amber pair, TMT/iTRAQ multi-colour spectrum.
- Typography: Inter for chrome; JetBrains Mono for peptide sequences, m/z values, protein IDs.

## Figure budget — 12 figures

| # | Title | Part |
|---|---|---|
| 1 | LC-MS/MS workflow | Part 1 |
| 2 | DDA vs DIA acquisition | Part 2 |
| 3 | MS2 spectrum identification + FDR | Part 3 |
| 4 | Quantification methods comparison | Part 4 |
| 5 | Plasma proteome biomarker landscape | Part 6 |
| 6 | Chemoproteomics workflow | Part 6 |
| 7 | Proteomics analysis workflow | Part 7 |
| 8 | Peptide fragmentation: b- and y-ions | Part 1 |
| 9 | Target-decoy FDR estimation | Part 3 |
| 10 | PTM peptide mass shifts | Part 3 |
| 11 | Missing-value patterns in proteomics | Part 5 |
| 12 | Phosphoproteomics enrichment workflow | Part 3 |

---

## Figure 1 — LC-MS/MS workflow

**File**: `diagrams/lecture-27/01-lc-ms-ms.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart with annotated stages:

1. Sample (cell / tissue / plasma).
2. Lysis + protein extraction.
3. Tryptic digestion (cleavage after K and R).
4. Peptide cleanup.
5. LC column separation (~1-3 hr gradient).
6. ESI ionisation.
7. MS1 scan (parent peptide masses).
8. MS2 fragmentation (HCD).
9. Detector → spectra.
10. Database search → identifications.
11. Quantification → protein abundance matrix.

Side annotations: typical scales (~50,000 peptides identified per LC-MS/MS run, $300/sample).

---

## Figure 2 — DDA vs DIA acquisition

**File**: `diagrams/lecture-27/02-dda-vs-dia.svg`
**ViewBox**: `0 0 1200 600`

Two side-by-side LC-MS/MS chromatograms:

- Top: DDA — vertical lines marking only the most-intense MS2 selections; many peptides skipped.
- Bottom: DIA — densely-packed vertical bands marking systematic m/z windows; comprehensive but chimeric.

Side panel: peptide identification overlap across replicates — DDA ~70%, DIA ~85%.

Annotation: "DIA = compressed sensing — under-determined acquisition recovered by spectral library deconvolution."

---

## Figure 3 — MS2 spectrum identification + FDR

**File**: `diagrams/lecture-27/03-ms2-id.svg`
**ViewBox**: `0 0 1200 720`

Top: example MS2 spectrum — m/z (x-axis) × intensity (y-axis) — with annotated b- and y-ion peaks identifying peptide LSDPYHRGSP.

Middle: target-decoy curve — score distribution for target hits (cobalt) and decoy hits (red).

Bottom: FDR vs score threshold — sliding line shows recovery rate at each threshold; 1% FDR cutoff marked.

Side annotation: "modern tools use ML scorers (Percolator, MS2Rescore) trained on target-decoy features."

---

## Figure 4 — Quantification methods comparison

**File**: `diagrams/lecture-27/04-quantification.svg`
**ViewBox**: `0 0 1200 600`

Four side-by-side panels:

- Panel 1 — Label-free (LFQ): MS1 peak integration; ratio across samples.
- Panel 2 — SILAC: heavy/light peptide pairs; ratio = abundance.
- Panel 3 — TMT: 11-plex reporter ions in MS2; ratios calculated per multiplex.
- Panel 4 — Targeted (PRM): selected reaction monitoring; absolute quantification.

Comparison table at bottom: precision, multiplexing, cost, applicability.

---

## Figure 5 — Plasma proteome biomarker landscape

**File**: `diagrams/lecture-27/05-biomarkers.svg`
**ViewBox**: `0 0 1200 720`

A horizontal log-scale axis: protein abundance in plasma (10⁻¹³ to 10⁻³ M).

Annotated proteins at characteristic abundances:

- Albumin: ~10⁻⁴ M (most abundant).
- IgG: ~10⁻⁵ M.
- Cardiac troponin: ~10⁻¹⁰ M.
- BNP: ~10⁻¹¹ M.
- Pancreatic biomarkers: 10⁻¹² M.
- Amyloid-β-42: 10⁻¹³ M (Alzheimer's).

Side panel: detection technology by abundance (Olink at 10⁻⁹+, Simoa at 10⁻¹³+).

Annotation: "the plasma proteome spans 10 orders of magnitude; abundance is the limiting factor for clinical biomarker validation."

---

## Figure 6 — Chemoproteomics workflow

**File**: `diagrams/lecture-27/06-chemoproteomics.svg`
**ViewBox**: `0 0 1200 600`

Top-to-bottom flow:

1. Drug + reactive linker → drug-probe.
2. Apply to live cells.
3. UV-crosslink: drug bound to nearby proteins fixed.
4. Lyse cells → digest.
5. Pull-down crosslinked peptides.
6. LC-MS/MS identifies modified peptides.
7. Map back to source proteins → drug-target list (and off-targets).

Side panel: example outcome — kinase inhibitor X identified ~20 binders, 5 expected (intended targets), 15 off-targets.

---

## Figure 7 — Proteomics analysis workflow

**File**: `diagrams/lecture-27/07-workflow.svg`
**ViewBox**: `0 0 1200 720`

Horizontal banded flowchart:

- Sample acquisition + extraction.
- LC-MS/MS (DDA or DIA).
- Database search (MaxQuant / FragPipe / DIA-NN).
- FDR control.
- Quantification.
- Statistical analysis (limma / DEqMS).
- Pathway enrichment.
- Network analysis.

Each stage with reproducibility-critical step annotations (database version, search engine settings, FDR cutoff, statistical model).

---

## Figure 8 — Peptide fragmentation: b- and y-ions

**File**: `diagrams/lecture-27/08-fragmentation.svg`
**ViewBox**: `0 0 1200 540`

A peptide LSDPYHRGSP shown horizontally with fragmentation cleavage points marked at each peptide bond.

For each cleavage:

- b-ion (N-terminal fragment) labelled with mass.
- y-ion (C-terminal fragment) labelled with mass.

A separate panel shows the predicted MS2 spectrum: peaks at b-ion masses (cobalt), peaks at y-ion masses (red).

Annotation: "the full b/y series uniquely identifies the peptide sequence."

---

## Figure 9 — Target-decoy FDR estimation

**File**: `diagrams/lecture-27/09-target-decoy.svg`
**ViewBox**: `0 0 1200 540`

Two overlapping histograms:

- Target hits (cobalt) — score distribution from real-protein-DB matches.
- Decoy hits (red) — score distribution from reverse-DB matches.

Vertical line at any score threshold defines:

- Target hits accepted = (target hits with score ≥ threshold).
- Decoy hits accepted = (decoy hits with score ≥ threshold).
- FDR = decoy / target.

1% FDR threshold marked.

Annotation: "decoys are by definition false positives; their rate scales with target false-positive rate."

---

## Figure 10 — PTM peptide mass shifts

**File**: `diagrams/lecture-27/10-ptm-shifts.svg`
**ViewBox**: `0 0 1200 540`

Table showing common PTMs and their mass shifts:

- Phosphorylation: +79.97 Da.
- Acetylation: +42.01 Da.
- Methylation: +14.02 Da.
- Mono-ubiquitination: +114.04 Da (GG remnant).
- Methionine oxidation: +15.99 Da.
- Glycosylation: complex (varies by glycan).

Each PTM highlighted on a peptide diagram showing where on the residue it attaches.

Side panel: search-space expansion as variable modifications are added (combinatorial growth).

---

## Figure 11 — Missing-value patterns in proteomics

**File**: `diagrams/lecture-27/11-missing-values.svg`
**ViewBox**: `0 0 1200 600`

Heat map: ~50 proteins (rows) × 30 samples (columns); cells coloured by abundance, white if missing.

Two clear patterns visible:

- Some proteins detected in all samples (top half of matrix).
- Some proteins detected in subset of samples (bottom half) — partly missing-not-at-random (low abundance), partly missing-at-random (stochastic).

Side panel: imputation strategies + comparison.

---

## Figure 12 — Phosphoproteomics enrichment workflow

**File**: `diagrams/lecture-27/12-phospho-enrich.svg`
**ViewBox**: `0 0 1200 600`

Top-to-bottom:

1. Tryptic digest (mostly non-phosphopeptides).
2. TiO2 / IMAC affinity column → enriched phosphopeptides.
3. LC-MS/MS as standard.
4. Variable phospho modification search.
5. Phosphosite localisation (Mascot / Andromeda confidence scoring).
6. Output: phosphosite × sample abundance matrix.

Side panel: typical scale (~20,000 phosphosites identified per experiment).
