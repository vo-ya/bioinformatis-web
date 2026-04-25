# Lecture 15 — Figures Specification

> **Scope**: Static diagrams for Lecture 15 (Protein Structure Prediction, AlphaFold-era).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-15.md`.

---

## 0. Conventions for This Lecture

- Figures are custom SVG; content is architecture-diagram and geometry-heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Architecture diagrams: rounded rectangles for tensor / module blocks with tensor-shape annotations in JetBrains Mono.
- Protein structure depictions: simplified cartoon glyphs — helices as cylinders or ribbons, β-strands as arrows, loops as smooth curves, disordered regions as wiggly lines.
- Colour palette: `--accent` cobalt for MSA/pair representations; amber `#b45309` for structure/coordinates; teal `#0d7377` for confidence metrics; `--fg-muted` grey for legacy/baseline content.
- pLDDT colouring: low (<50) red `#c4342c`, medium (50–70) orange `#b45309`, confident (70–90) cobalt, high (>90) teal-dark.
- Typography: Inter for UI labels; JetBrains Mono for tensor shapes, amino-acid strings, residue indices, numerical values.
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Four levels of protein structure | Part 1 | Custom SVG |
| 2 | CASP GDT_TS progression | Part 2 | Custom SVG |
| 3 | Coevolution and 3D contacts | Part 3 | Custom SVG |
| 4 | DCA contact prediction | Part 3 | Custom SVG |
| 5 | AlphaFold2 architecture overview | Part 4 | Custom SVG |
| 6 | Evoformer block | Part 4 | Custom SVG |
| 7 | Structure module + IPA | Part 4 | Custom SVG |
| 8 | pLDDT and PAE confidence metrics | Part 4 | Custom SVG |
| 9 | ProteinMPNN inverse folding | Part 6 | Custom SVG |
| 10 | RFDiffusion protein design pipeline | Part 6 | Custom SVG |
| 11 | AlphaFold Database scale | Part 7 | Custom SVG |
| 12 | AlphaFold3 diffusion structure module | Part 5 | Custom SVG |

---

## Figure 1 — Four levels of protein structure

**File**: `diagrams/lecture-15/01-structure-hierarchy.svg`
**Lecture anchor**: §1.1 Why structure matters
**ViewBox**: `0 0 1080 520`

### Purpose

Introduce the primary → quaternary hierarchy using a single illustrative protein (haemoglobin is a good choice — multi-subunit and visually canonical).

### Content

**Four horizontal panels stacked vertically**, each labelled with its level:

1. **Primary**. A linear strip of ~50 amino-acid letter-boxes reading "V H L T P E E K S A V T A L W G K V N V D E V G G E A L G …" using standard single-letter codes. Each box coloured weakly by amino-acid class (hydrophobic grey, polar teal, charged+ red, charged- blue).

2. **Secondary**. Same sequence; local motifs highlighted — a ~15-residue stretch rendered as a helical spiral; another stretch as an arrow (β-strand). Labels: "α-helix", "β-strand", "loop". Annotation: "backbone H-bond patterns".

3. **Tertiary**. A 3D cartoon of one haemoglobin α-subunit: helices as cylinders, loops as smooth curves, fold organisation visible. Annotation: "globular 3D fold; side chains packed inside".

4. **Quaternary**. Four subunits (2 α + 2 β) arranged in the assembled haemoglobin tetramer, with a heme group visible in each. Annotation: "functional unit; multi-chain".

### Style notes

- Each panel separated by a thin horizontal rule.
- Chain trace in panel 3: smooth Bezier with cylinders for helices.
- Panel 4: four subunits with their chain labels (α₁, α₂, β₁, β₂).

---

## Figure 2 — CASP GDT_TS progression

**File**: `diagrams/lecture-15/02-casp-progression.svg`
**Lecture anchor**: §2.3 CASP: the community benchmark
**ViewBox**: `0 0 1200 480`

### Purpose

Show the trajectory of structure-prediction quality across CASP editions from 1994 to 2024, highlighting the AlphaFold inflection.

### Content

**Line plot.** X-axis: CASP edition, labelled CASP1 (1994) through CASP16 (2024). Y-axis: median best GDT_TS on Free-Modelling targets, 0–100.

**Curve**:

- CASP1–CASP8: flat trajectory ~15–25.
- CASP9–CASP12 (2010–2016): gentle climb to ~35–40 (Rosetta era).
- CASP13 (2018): jump to ~60 (**AlphaFold1**, annotated).
- CASP14 (2020): jump to ~90 (**AlphaFold2**, annotated with an amber spotlight).
- CASP15–CASP16: plateau around 90–92.

**Era annotations above the curve**:

- "Rosetta fragment assembly era" (spanning CASP7–CASP12).
- "Contact-prediction era" (CASP12–CASP13).
- "AlphaFold2 era" (CASP14+).

**Horizontal reference lines**:

- GDT_TS = 50 ("recognisable fold").
- GDT_TS = 70 ("good prediction").
- GDT_TS = 90 ("experimental quality") — dashed amber.

### Style notes

- Dots at each CASP edition, line connecting them.
- Era bands: thin horizontal stripes in `--bg-subtle` for each era.
- AlphaFold milestones: amber starburst markers.

---

## Figure 3 — Coevolution and 3D contacts

**File**: `diagrams/lecture-15/03-coevolution.svg`
**Lecture anchor**: §3.1 Why coevolving residues are spatially close
**ViewBox**: `0 0 1080 540`

### Purpose

Visualise the mental link between MSA column co-mutations and 3D spatial contacts.

### Content

**Top panel — MSA snippet.** ~10 rows × 40 columns. Use monospace amino-acid letters with a pale tinted background per letter. Two specific columns (e.g. index 12 and 34) highlighted with vertical amber tint strips. Within the tint strips:

- Row 1: L at col 12, V at col 34.
- Row 2: L at col 12, V at col 34.
- Row 3: I at col 12, T at col 34.
- Row 4: I at col 12, T at col 34.
- Row 5: V at col 12, A at col 34.
- Row 6: L at col 12, V at col 34.
- (And so on — illustrating the L↔V and I↔T co-mutation pattern.)

Below the MSA: "Columns 12 & 34 co-mutate: when one changes, the other usually does."

**Middle panel — arrow transition.** A labelled downward arrow: "The MSA co-mutation pattern implies 3D proximity".

**Bottom panel — 3D cartoon.** A folded protein rendered as a backbone trace. Residues 12 and 34 shown as filled spheres touching each other in space (in close contact), elsewhere along the chain. Annotation arrow from each sphere pointing outward: "residue 12", "residue 34". A dashed line between the two spheres labelled "3D contact (~5 Å)".

### Style notes

- MSA: JetBrains Mono amino-acid letters with subtle per-residue-class tints.
- 3D backbone: smooth cobalt curve.
- Contact residues: filled amber spheres with amber dashed line.

---

## Figure 4 — DCA contact prediction

**File**: `diagrams/lecture-15/04-dca-contacts.svg`
**Lecture anchor**: §3.3 Direct Coupling Analysis (DCA)
**ViewBox**: `0 0 1200 440`

### Purpose

Contrast naive mutual information vs inverse-covariance DCA for contact recovery, against a ground-truth experimental contact map.

### Content

**Three side-by-side square matrices**, each ~100×100 residues.

**Left — raw MI matrix**. Heatmap with a diffuse, noisy pattern: diagonal band strong, many off-diagonal bright blobs smeared across wide regions. Label: "Mutual information (MI)". Annotation below: "Indirect correlations confuse the signal".

**Middle — DCA inverse-covariance contact scores**. Heatmap with crisper features: diagonal band visible, several isolated bright spots off-diagonal marking direct contacts, most of the MI noise suppressed. Label: "DCA / plmDCA". Annotation: "Inverse covariance isolates direct contacts".

**Right — native contact map from PDB**. Binary heatmap (black = contact if $d < 8$ Å, white otherwise). Label: "Experimental ground truth".

**Annotation above all three panels**: "~80% of top-ranked DCA contacts are true contacts in the PDB structure".

### Style notes

- All three matrices same size, aligned.
- Colour scale: sequential cobalt for MI and DCA; binary for ground truth.
- Small axis labels: "residue index" on x and y of each panel.

---

## Figure 5 — AlphaFold2 architecture overview

**File**: `diagrams/lecture-15/05-af2-overview.svg`
**Lecture anchor**: §4.2 Overall architecture
**ViewBox**: `0 0 1200 620`

### Purpose

Top-level block diagram of the AlphaFold2 network from input to output.

### Content

**Top — Input band**. Three blocks side by side:

- "Protein sequence (length L)" — cobalt.
- "MSA (up to 5120 × L)" — violet.
- "Templates (optional, K × L)" — amber.

**Middle — Evoformer block**. Large rounded rectangle labelled "Evoformer (48 blocks)". Inside, two smaller boxes:

- "MSA representation (s × L × 256)".
- "Pair representation (L × L × 128)".

Arrows between them indicating bidirectional updates. Annotation at right: "axial attention + triangle updates".

**Middle — Structure module**. Rounded rectangle labelled "Structure module (8 iterations)". Inside: "Single representation (L × 384)", "Invariant Point Attention", "SE(3) frame updates".

**Right — Output band**. Three output blocks:

- "3D coordinates (L × atoms × 3)" — amber.
- "Per-residue pLDDT (L)" — teal.
- "PAE matrix (L × L)" — teal.

**Feedback arrow**. Labelled "recycling (×3–4)", curved from output back to input of Evoformer.

**Bottom-left parameter-count annotation**. "~93M trainable parameters. Training: 11 days on 128 TPUv3 cores."

### Style notes

- Blocks with rounded corners, consistent border thickness.
- Tensor shapes in monospace within each block.
- Arrows with arrow-head markers.
- Recycling feedback arrow dashed.

---

## Figure 6 — Evoformer block

**File**: `diagrams/lecture-15/06-evoformer-block.svg`
**Lecture anchor**: §4.3 The Evoformer
**ViewBox**: `0 0 1200 620`

### Purpose

Show the six sub-operations within one Evoformer block and how MSA / pair representations are bidirectionally updated.

### Content

**Two horizontal bands**.

**Top band — MSA representation flow**. Rectangles left-to-right:

- Input MSA (s × L × 256).
- Row-wise attention with pair bias.
- Column-wise attention.
- Transition (feedforward).
- Output MSA.

Complexity annotations below each: O(s × L²), O(s² × L), O(s × L), etc.

**Bottom band — Pair representation flow**. Rectangles:

- Input pair (L × L × 128).
- Outer-product mean ← MSA (arrow from top band down).
- Triangle multiplicative update.
- Triangle self-attention.
- Transition.
- Output pair.

**Cross-connections**:

- Pair → row-wise MSA attention (as bias).
- MSA → pair (via outer-product-mean).

**Right-side annotation box**. Title "Why this works". Bullets:

- Axial attention = quadratic in one axis at a time (vs full 2D attention).
- Triangle updates = respect the (i,j)+(j,k)→(i,k) geometry.
- Block repeats 48 times; each round refines both representations.

### Style notes

- Two-band layout with cross-connections drawn as diagonal arrows.
- Each sub-operation is a labelled rectangle.
- Annotations in small italic near each operation.

---

## Figure 7 — Structure module + IPA

**File**: `diagrams/lecture-15/07-structure-module.svg`
**Lecture anchor**: §4.5 The structure module + IPA
**ViewBox**: `0 0 1200 520`

### Purpose

Visualise the iterative SE(3)-equivariant refinement performed by the structure module, and the specific geometry of Invariant Point Attention.

### Content

**Top band — iterative refinement timeline**. Three snapshots left-to-right:

- Iteration 0: all residues at origin, random frame orientations (cartoon of chaos).
- Iteration 4: partial convergence — backbone has emerged roughly, some loops out of place.
- Iteration 8: fully-folded cartoon structure.

Arrows between snapshots labelled "IPA update".

**Bottom band — IPA detail**. A cartoon of two residues *i* and *j* with their SE(3) frames drawn as local coordinate axes. Query points emanate from residue *i*'s frame; key points from residue *j*'s frame. A dashed bracket shows "distance computed in a common reference frame". Annotation above: "attention weight = f(distance(query_i, key_j))".

**Bottom-right annotation box**. Title "SE(3) equivariance". Bullets:

- Rotate input structure by R ∈ SO(3) → IPA output rotates by R.
- No data augmentation by rotation needed during training.
- Side chains built at the end from per-residue torsion outputs.

### Style notes

- Top band: three cartoon structures evolving left-to-right.
- Bottom band: residues drawn with small coordinate-axis glyphs.
- Dashed arrows for reference-frame projection.

---

## Figure 8 — pLDDT and PAE confidence metrics

**File**: `diagrams/lecture-15/08-plddt-pae.svg`
**Lecture anchor**: §4.6 Recycling and confidence
**ViewBox**: `0 0 1200 520`

### Purpose

Show how pLDDT and PAE are computed and visualised for a typical multi-domain protein; teach the reader to read both.

### Content

**Left panel — 3D cartoon of a two-domain protein**. Two folded globular domains connected by a flexible linker. Backbone coloured by pLDDT:

- Domain 1 core: teal-dark (very high, >90).
- Domain 1 surface loops: cobalt (high, 70–90).
- Linker: red (<50) with wavy disorder-like styling.
- Domain 2 core: teal-dark.
- Domain 2 surface loops: cobalt.

**Middle panel — PAE matrix**. L × L heatmap. Diagonal blocks of low PAE (cobalt-dark) corresponding to each domain. Off-diagonal blocks (inter-domain) bright amber showing high PAE. Linker residues: horizontal + vertical amber stripes.

**Right panel — interpretation summary text**. Boxed. Bullet list:

- "Individual domains: confidently folded (high pLDDT, low intra-domain PAE)."
- "Linker: disordered; ignore its exact coordinates."
- "Inter-domain orientation: **not** constrained by the prediction. Use each domain independently for function inference."

### Style notes

- pLDDT colour scale legend below the 3D panel.
- PAE colour scale legend below the matrix.
- Interpretation box: `--bg-subtle` background, amber border.

---

## Figure 9 — ProteinMPNN inverse folding

**File**: `diagrams/lecture-15/09-proteinmpnn.svg`
**Lecture anchor**: §6.2 ProteinMPNN
**ViewBox**: `0 0 1200 520`

### Purpose

Show ProteinMPNN's backbone-to-sequence mapping via a GNN with k-nearest-neighbour edges.

### Content

**Left — Backbone input**. 3D cartoon of a small (~50 residue) helix bundle backbone. Residue nodes as spheres at Cα positions; no side chains. Backbone drawn as a smooth cobalt curve.

**Middle — GNN message passing**. Same nodes, edges drawn as k-nearest-neighbour lines (each node has ~10 edges to spatially nearby residues). Small "message-passing" glyphs (tiny arrows along edges). Annotation: "geometric features per edge: distance, angle, orientation".

**Right — Per-residue softmax output**. 20-column bar chart per residue position showing amino-acid probability. Top positions in bar chart (tallest bars) highlighted.

**Bottom band — Sampled sequence** (shown beneath the 3D structure):

- Sequence: "M A V L K E R I Q S T L G N D P M L …"
- Label: "Autoregressive decode: residue i sampled given residues 1..i-1".

**Right-side annotation box**. "~50% of sampled sequences fold correctly to target (AF2 validated). ~10% with Rosetta baseline."

### Style notes

- Backbone in cobalt.
- KNN edges in light grey.
- Per-residue softmax: 20 horizontal bars with AA-letter labels.
- Decoded sequence in monospace at the bottom.

---

## Figure 10 — RFDiffusion protein design pipeline

**File**: `diagrams/lecture-15/10-rfdiffusion.svg`
**Lecture anchor**: §6.3 De novo protein design
**ViewBox**: `0 0 1200 480`

### Purpose

Horizontal pipeline showing the design loop from constraint specification through experimental validation.

### Content

**Five-panel horizontal workflow**.

**Panel 1 — Constraints**. A small schematic of user-specified constraints: "bind Zn²⁺ ion", "coordinate via His-His-His-Glu motif", "overall topology: 4-helix bundle". Annotations as bullet list.

**Panel 2 — RFDiffusion backbone**. A novel backbone cartoon generated by RFDiffusion, fulfilling the constraints. Metal-binding residue positions marked.

**Panel 3 — ProteinMPNN sequence**. The generated amino-acid sequence overlaid as a monospace string below the backbone. Colour-coded by residue class.

**Panel 4 — AlphaFold2 validation**. Side-by-side: target backbone (from panel 2) and AF2-predicted backbone (from the sequence). Both overlaid with RMSD = 1.2 Å noted. Green checkmark.

**Panel 5 — Wet-lab validation**. A cartoon of a protein expressed in *E. coli*, purified, tested for binding. Label: "2 of 12 designs bound Zn²⁺ successfully".

**Bottom caption band**. "10–40% success rate for RFDiffusion + ProteinMPNN designs. Pre-AlphaFold-era success rate: <1%."

### Style notes

- Arrows between panels with small step-number badges.
- Backbone cartoons consistent across panels 2, 4.
- Sequence panel: monospace string with colour-by-class.

---

## Figure 11 — AlphaFold Database scale

**File**: `diagrams/lecture-15/11-af-database.svg`
**Lecture anchor**: §7.1 The AlphaFold Database
**ViewBox**: `0 0 1080 480`

### Purpose

Visualise the scale shift from PDB (experimental) to the AlphaFold Database (predicted).

### Content

**Top — comparison bar chart**. Two bars:

- PDB experimental structures (2021): ~180,000. Small bar.
- AlphaFold DB predicted structures (2024): ~200,000,000. Bar ~1000× longer. Log-scale x-axis.

Bar colouring: PDB amber, AF DB cobalt.

**Middle — AF DB quality distribution**. Horizontal stacked bar showing fraction by pLDDT bucket:

- >90 (very high): ~40%. Teal-dark.
- 70–90 (confident): ~30%. Cobalt.
- 50–70 (low): ~20%. Amber.
- <50 (disordered / orphan): ~10%. Red.

Legend on the right with pLDDT colour conventions.

**Bottom — impact bullets**.

- "Every UniProt protein: predicted structure available, queryable by UniProt ID."
- "Drug targets: historically a handful with experimental structures → now tens of thousands."
- "Functional annotation: structural homology searches (foldseek) enable function transfer across distant homologues."

### Style notes

- Bars drawn with log-scale x-axis; tick marks at 10³, 10⁴, 10⁵, 10⁶, 10⁷, 10⁸.
- Bullet box: `--bg-subtle` background.

---

## Figure 12 — AlphaFold3 diffusion structure module

**File**: `diagrams/lecture-15/12-af3-diffusion.svg`
**Lecture anchor**: §5.1 AlphaFold3 architecture
**ViewBox**: `0 0 1200 520`

### Purpose

Show the diffusion-based coordinate generation used by AlphaFold3, in contrast to AlphaFold2's IPA iterative updates, and highlight the multi-component input handling.

### Content

**Top half — forward diffusion process**. Four snapshots left-to-right:

- $t = 0$: clean folded structure.
- $t = T/4$: slight coordinate noise.
- $t = 3T/4$: heavy noise, recognisable contours fading.
- $t = T$: pure Gaussian noise.

Arrow labelled "add noise, T timesteps".

**Bottom half — reverse (learned denoising) process**. Four snapshots right-to-left (reversing the top):

- $t = T$: noise.
- $t = 3T/4$: starting to cluster.
- $t = T/4$: clearly folding.
- $t = 0$: clean structure.

Each transition labelled "denoise (conditioned on Pairformer)".

**Right-side inset — multi-component input**. A small diagram of: one protein chain (cobalt ribbon) + one double-stranded DNA (rainbow helix) + one small-molecule ligand (ball-and-stick) — all being denoised together in a shared diffusion process.

**Bottom caption band**. "Diffusion handles variable-composition inputs (add or remove a ligand; add a chain) and multi-modal output distributions (multiple conformations) — which the AF2 IPA module cannot."

### Style notes

- Top/bottom split clearly demarcated.
- Noise levels visualised as progressively blurrier coordinate scatters.
- Inset: smaller multi-component cartoon.

### Required XML escape

- Use `&amp;` / `&lt;` / `&gt;` in any embedded text.
- Entity `T` used in timestep annotation is fine as plain letter.
