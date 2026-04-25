# Lecture 15 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 15 (Protein Structure Prediction, AlphaFold-era).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-15.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-15/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- MSA / pair-repr colours: MSA violet, pair cobalt, structure amber.
- Amino-acid letter colouring: hydrophobic grey, polar teal, positive `#c4342c`, negative `#1e3a8a`, cysteine/special `#b45309`.
- pLDDT bands: >90 teal-dark, 70–90 cobalt, 50–70 amber, <50 red.
- Protein cartoons: helices as cylinders, strands as arrows, loops as smooth lines (simplified 2D rendering fine — no GL required).
- Typography: Inter for UI chrome; JetBrains Mono for amino-acid sequences, tensor shapes, residue indices, PAE values, pLDDT scores.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters (MSA depth, noise levels, AF confidence threshold).
- **Step / Run / Reset** — for iterative simulations (Evoformer blocks, IPA iterations, diffusion timesteps).
- **Re-simulate** — for stochastic simulations (sampled sequences, diffusion realisations).
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- MSA Coevolution and Contact Prediction → precision/recall of predicted contacts vs ground-truth; MI vs DCA comparison.
- Evoformer Axial-Attention Visualiser → per-step tensor shapes; trace of how pair representation refines over blocks.
- IPA / Structure Module Walkthrough → per-iteration RMSD to target; equivariance check (rotate input, output rotates the same).
- pLDDT & PAE Confidence Interpreter → per-region trust verdict + recommended downstream use.
- ProteinMPNN Inverse-Folding Demo → diversity of sampled sequences; self-consistency score (predicted fold vs target).
- AlphaFold Database Explorer → filtered set size; pLDDT distribution of the filtered set.
- Protein Structure Level Explorer → cross-level highlighting (click residue → highlight across all levels).

### 1.5 Feasibility gate on user input (required where input is free-form)

- All seven artifacts use preset / slider / dropdown controls only. No free-form sequence input beyond validated amino-acid strings (A–Z alphabet, length ≤ 500).
- Numeric inputs clamp to their valid range.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted for inline math; otherwise none).
- No external data files > 100 KB (presets, example MSAs, tensor traces all hardcoded).
- No WebGL / 3D rendering required — 2D cartoon views are sufficient; any "3D" view is a rotatable 2D projection.

---

## 2. Artifact #1 — MSA Coevolution and Contact Prediction

**File**: `artifacts/lecture-15/01-coevolution-contacts.html`
**Lecture anchor**: §3.3 Direct Coupling Analysis (DCA)
**EE framing reinforced**: MI vs inverse-covariance; graphical-model estimation from categorical samples.

### Teaching purpose

Build a toy protein family under a user-specified contact topology. Simulate an MSA by sampling from a Potts model with those contacts. Compute mutual information and DCA contact scores. Compare against ground truth. Show that DCA recovers contacts cleanly; MI is confounded by indirect correlations.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Contact topology:  [ 1D chain ▾ / 2D sheet / helix / random]│
│ # residues:        [──●── 30 ]                              │
│ MSA depth:         [──●── 500 ]                             │
│ Coupling strength: [──●── 0.8 ]                             │
│ [Re-simulate MSA]                                           │
├─────────────────────────────────────────────────────────────┤
│ Row: three matrices side by side                            │
│   MI (noisy)       DCA (clean)       Ground truth           │
├─────────────────────────────────────────────────────────────┤
│ Precision/recall at top-K contacts:                         │
│   K=20: MI 40% / DCA 90%                                    │
│   K=50: MI 30% / DCA 80%                                    │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ DCA recovers 18 / 20 top contacts correctly.           │   │
│ │ MI recovers 8 / 20 — half are indirect correlations.   │   │
│ │ Reduce MSA depth to 50 → DCA precision drops; MI barely│   │
│ │   changes (already noisy).                             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Contact topology: 1D chain / 2D sheet / α-helix / random (default α-helix).
- Residue count: 10 – 50, default 30.
- MSA depth: 30 – 2000 (log slider), default 500.
- Coupling strength: 0 – 1.5, default 0.8.
- Re-simulate button.

### What they see

- Three contact-score heatmaps side-by-side (MI, DCA, ground truth).
- Precision/recall table at top-K.
- Outcome banner.

### Target aha

DCA is clean; MI is noisy. Reduce MSA depth → DCA's edge shrinks (DCA needs enough samples). Increase coupling strength → both improve but DCA dominates at moderate strengths. Student sees inverse-covariance estimation is the right framework; MI is confounded by transitivity.

### Technical notes

- Pure JS, seeded PRNG.
- Sample MSA from a Potts model using Metropolis-Hastings: Hamiltonian includes user-specified contacts with couplings ± coupling strength.
- MI: standard per-pair from joint + marginals of the 20-letter alphabet (or 4-letter simplified alphabet for speed).
- DCA: implement mean-field approximation (mfDCA). Invert the coupling matrix from frequency statistics.
- Ground truth: exactly the contacts in the user-chosen topology.
- Precision/recall: top-K predicted contacts vs ground-truth set.

### Acceptance criteria

- [ ] Default (α-helix, 30 res, 500 depth, 0.8 strength) → DCA precision > 80% at top-20; MI precision < 50%.
- [ ] Reducing MSA depth to 30 → both drop; DCA still beats MI.
- [ ] Switching topology to random → still holds.
- [ ] Setting coupling to 0 → both near random.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Evoformer Axial-Attention Visualiser

**File**: `artifacts/lecture-15/02-evoformer-attention.html`
**Lecture anchor**: §4.3 The Evoformer
**EE framing reinforced**: axial attention as quadratic-in-one-axis approximation to full 2D attention; triangle updates as geometric regularisation.

### Teaching purpose

Step-by-step walk through one Evoformer block's six sub-operations. Show tensor shapes before and after; show how pair representation sharpens as more Evoformer blocks run.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Input preset: [ toy random ▾ / helix-like / sheet-like ]    │
│ # residues L:   [──●── 20 ]                                 │
│ MSA depth s:    [──●── 16 ]                                 │
│ # Evoformer blocks to run: [──●── 8 / 16 / 32 / 48 ]       │
│ [Step one op]  [Run all]  [Reset]                          │
├─────────────────────────────────────────────────────────────┤
│ Current operation: "row-wise MSA attention (with pair bias)"│
│ Inputs / outputs shown as small heatmaps with shape labels: │
│   MSA repr (s × L × c) before / after                       │
│   Pair repr (L × L × c) before / after                      │
├─────────────────────────────────────────────────────────────┤
│ Pair-repr "contact-score" snapshot over blocks:             │
│   After block 0:  diffuse                                   │
│   After block 8:  partly resolved                           │
│   After block 16: sharper contacts visible                  │
│   Animated progression                                      │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Currently after 8 blocks. Pair repr already shows the  │   │
│ │ dominant contact pattern; more blocks just sharpen it. │   │
│ │ Switch to helix preset → diagonal bands emerge.        │   │
│ │ Switch to sheet preset → cross-strand contacts emerge. │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Input preset: toy random / helix-like / sheet-like.
- Residue count L: 8 – 40, default 20.
- MSA depth s: 4 – 64, default 16.
- # Evoformer blocks: 8 / 16 / 32 / 48.
- Step / Run / Reset buttons.

### What they see

- Current-operation highlight with tensor-shape labels.
- Before / after heatmap snapshots per step.
- Pair-repr evolution across blocks.
- Outcome banner.

### Target aha

With each block, the pair representation's "contact-score slice" sharpens toward the intended fold topology. Axial attention (row-wise then column-wise) is cheaper than full 2D but captures enough. Triangle updates respect geometric consistency. Student sees the *iterative refinement* nature of the Evoformer.

### Technical notes

- Pure JS.
- Simplified arithmetic (reduced channel count, fewer heads) that approximates Evoformer behaviour without full model weights.
- Preset topologies provide ground-truth-like signal in the starting pair repr; operations sharpen it.
- Heatmaps via Canvas or SVG.
- Animations between steps ≤ 400 ms.

### Acceptance criteria

- [ ] Each sub-operation individually steppable.
- [ ] Pair-repr sharpening is visible after 8 blocks.
- [ ] Preset switching shows different contact patterns.
- [ ] Tensor shapes displayed correctly at each step.
- [ ] Opens pre-rendered with pair-repr at block 0.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — IPA / Structure Module Walkthrough

**File**: `artifacts/lecture-15/03-ipa-walkthrough.html`
**Lecture anchor**: §4.5 The structure module + IPA
**EE framing reinforced**: SE(3) equivariance; attention-based geometric refinement.

### Teaching purpose

Visualise the iterative SE(3)-equivariant refinement performed by 8 IPA iterations. Show how residue frames (position + orientation) evolve. Demonstrate equivariance: rotating the input by any R ∈ SO(3) produces outputs rotated by the same R.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Target structure preset: [ 2-helix ▾ / β-hairpin / bundle ] │
│ # residues: [──●── 20 ]                                     │
│ Iteration: [──●── 0 / 1 / … / 8 ]                           │
│ [Step]  [Run all 8]  [Reset]                                │
│ Rotation test: [Rotate input by R] → check output           │
├─────────────────────────────────────────────────────────────┤
│ 2D projection panel (top-down view):                        │
│   Residue positions + frame axes (small glyphs)             │
│   Target structure shown as dashed outline                  │
│   RMSD to target: 3.2 Å at iteration 3, 0.8 Å at iteration 8│
├─────────────────────────────────────────────────────────────┤
│ IPA detail inset:                                           │
│   Residue i's frame; query points drawn                     │
│   Residue j's frame; key points drawn                       │
│   Common-reference-frame projection → attention weight      │
├─────────────────────────────────────────────────────────────┤
│ Equivariance check panel:                                   │
│   Output of IPA on input vs R·input                         │
│   Difference should be R·(original output) → show residual  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 8 IPA iterations: RMSD 10.2 Å → 0.8 Å.                 │   │
│ │ Equivariance check: max residual 2e-6 (numerical).     │   │
│ │ Rotate input → output rotates identically; no R-specific│  │
│ │ parameter training needed.                             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Target structure preset: 2-helix / β-hairpin / small bundle.
- # residues: 10 – 30, default 20.
- Iteration number (slider or play button).
- Rotate button (applies R from user-specified angle; re-checks equivariance).

### What they see

- 2D top-down view of residue positions + frames.
- IPA-detail inset showing one query-key pair.
- Equivariance residual panel.
- Outcome banner with RMSD trajectory + equivariance.

### Target aha

Every iteration reduces RMSD to the target. Rotating the input rotates the output by the same R (residual near numerical zero). Student sees: IPA is attention in 3D, equivariant by construction; the structure module is iterative refinement of geometry.

### Technical notes

- Pure JS.
- Simplified IPA: per-residue frames as (2D position, rotation angle) for visualisation purposes.
- Toy force-law driving frames toward target: gradient of a pairwise distance loss.
- Equivariance check: compute output on input, then on R·input; compare with R applied to the original output; report residual.

### Acceptance criteria

- [ ] Default 2-helix preset converges to RMSD < 1 Å by iteration 8.
- [ ] Equivariance residual < 1e-4 across random rotations.
- [ ] Switching preset re-initialises and converges.
- [ ] Iteration slider correctly renders intermediate states.
- [ ] Opens pre-rendered at iteration 0.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — pLDDT and PAE Confidence Interpreter

**File**: `artifacts/lecture-15/04-plddt-pae.html`
**Lecture anchor**: §4.6 Recycling and confidence
**EE framing reinforced**: pLDDT as calibrated per-residue confidence; PAE as pair-wise aligned error.

### Teaching purpose

Load preset AlphaFold outputs representing common scenarios. Render the 3D structure coloured by pLDDT, the PAE matrix, and a natural-language verdict. Teach the reader to diagnose confidence patterns and decide which parts to trust.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Preset: [ well-folded single domain ▾ / multi-domain /      │
│           disordered region / orphan (shallow MSA) /        │
│           multi-state / membrane protein ]                   │
├─────────────────────────────────────────────────────────────┤
│ 3D cartoon (2D projection) coloured by pLDDT:               │
│   domain 1: teal-dark (high)                                │
│   linker: red (low, disordered)                             │
│   domain 2: teal-dark                                       │
├─────────────────────────────────────────────────────────────┤
│ PAE matrix (L × L heatmap):                                 │
│   two dark diagonal blocks = domains well-defined            │
│   bright off-diagonal = inter-domain orientation uncertain  │
├─────────────────────────────────────────────────────────────┤
│ Verdict card (interpretation):                              │
│   Within domain 1 (residues 1–100): trust coordinates.      │
│   Linker (101–130): disordered, ignore exact coords.         │
│   Within domain 2 (131–250): trust coordinates.             │
│   Inter-domain arrangement: DO NOT use for binding-site      │
│     prediction or docking.                                  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ High per-residue pLDDT ≠ reliable multi-domain orient. │   │
│ │ Always read pLDDT AND PAE; they answer different       │   │
│ │   questions.                                            │   │
│ │ Switch to orphan preset → pLDDT drops below 50 across  │   │
│ │   most residues; prediction unreliable.                │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Preset dropdown: six canonical scenarios.
- Threshold slider: pLDDT cutoff for "trustworthy" region highlighting (default 70).

### What they see

- 3D cartoon (2D projection) coloured by pLDDT.
- PAE matrix with axis ticks.
- Verdict card with per-region interpretation.
- Outcome banner.

### Target aha

pLDDT tells you about local coordinate accuracy; PAE tells you about relative position accuracy. Both matter; they answer different questions. Orphan / disordered / multi-domain proteins produce characteristic confidence-metric patterns — learn to read them.

### Technical notes

- Pure JS.
- Preset data: per-residue pLDDT arrays and L × L PAE matrices hardcoded to reflect realistic patterns.
- 2D projection: pre-computed 2D residue positions for each preset.
- Verdict card text: preset-specific interpretation strings with parameterised residue ranges.

### Acceptance criteria

- [ ] Default (multi-domain) renders 2D cartoon with domain colouring.
- [ ] PAE matrix visibly shows two diagonal blocks.
- [ ] Verdict card is per-preset correct.
- [ ] Orphan preset shows widespread low pLDDT.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — ProteinMPNN Inverse-Folding Demo

**File**: `artifacts/lecture-15/05-proteinmpnn.html`
**Lecture anchor**: §6.2 ProteinMPNN
**EE framing reinforced**: inverse folding as ill-posed inverse with learned structural prior; autoregressive decoding.

### Teaching purpose

Pick a target backbone. Run a simplified ProteinMPNN-style inverse-folding network. Sample multiple sequences. Show sequence diversity and a toy self-consistency score (predicted fold vs target).

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Target backbone: [ 4-helix bundle ▾ / β-barrel / random ]   │
│ # residues: [──●── 50 ]                                     │
│ Sampling temperature: [──●── 0.5 ]                          │
│ # samples: [──●── 10 ]                                      │
│ [Generate sequences]                                        │
├─────────────────────────────────────────────────────────────┤
│ 2D backbone view (top) with residue nodes + k-NN edges      │
├─────────────────────────────────────────────────────────────┤
│ Sampled sequences (up to 10, one per row):                  │
│   MAVLKER... [self-consistency 0.87]                        │
│   MATIKES... [self-consistency 0.82]                        │
│   IAVLKEK... [self-consistency 0.79]                        │
│   ...                                                       │
│ Native-recovery rate: 45% average (positions matching native)│
├─────────────────────────────────────────────────────────────┤
│ Sequence-diversity panel:                                   │
│   per-position amino-acid frequency across the 10 samples   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 10 sampled sequences; mean self-consistency 0.82.       │   │
│ │ Many distinct amino-acid choices at most positions —    │   │
│ │ different sequences, same fold.                         │   │
│ │ Raise temperature to 1.5 → diversity up, consistency    │   │
│ │   down. Lower to 0.1 → nearly deterministic.           │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Target backbone preset.
- # residues: 30 – 80, default 50.
- Sampling temperature: 0.1 – 2.0, default 0.5.
- # samples: 5 – 20, default 10.
- Generate button.

### What they see

- 2D backbone view with k-NN edges.
- List of sampled sequences + per-sample self-consistency.
- Per-position amino-acid frequency chart.
- Outcome banner.

### Target aha

Many distinct sequences can fold to the same backbone (ill-posed inverse). ProteinMPNN's learned structural prior picks the ones likely to fold correctly. Temperature trades diversity for fidelity.

### Technical notes

- Pure JS.
- Simplified ProteinMPNN: per-residue marginal amino-acid probabilities conditioned on local backbone geometry. Autoregressive sampling in a fixed order.
- Self-consistency: approximate by a learned lookup table mapping "local backbone features + amino acid" to "expected to fold correctly?" probability. Or a simpler proxy: hydrophobicity vs solvent exposure matching score.
- Native recovery: if the preset backbone was derived from a known protein, compare to its native sequence.

### Acceptance criteria

- [ ] Default → 10 sampled sequences produced; self-consistency ≥ 0.7 on average.
- [ ] Native-recovery rate ≈ 40–50% at temperature 0.5.
- [ ] Higher temperature → higher diversity, lower self-consistency.
- [ ] Switching backbone preset changes the generated-sequence distribution.
- [ ] Opens pre-rendered with 4-helix bundle default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — AlphaFold Database Explorer

**File**: `artifacts/lecture-15/06-afdb-explorer.html`
**Lecture anchor**: §7.1 The AlphaFold Database
**EE framing reinforced**: pLDDT as filterable confidence metadata; the DB is queryable.

### Teaching purpose

Present a small catalogue (~20–50 entries) of preset AlphaFold predictions spanning the confidence spectrum. Student filters by pLDDT, length, disorder fraction, function category. See which subsets of the database are trustworthy and which aren't.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Filters:                                                    │
│   Min mean pLDDT:  [──●── 70 ]                              │
│   Max disorder %:  [──●── 30 ]                              │
│   Length:          [──●── 50–500 ]                          │
│   Function category: [ any ▾ / enzyme / transporter /       │
│                       disordered / orphan ]                  │
├─────────────────────────────────────────────────────────────┤
│ Results (filtered from 50 preset entries):                  │
│   12 entries match                                          │
│   Grid of thumbnails: each showing 2D cartoon + pLDDT dist  │
├─────────────────────────────────────────────────────────────┤
│ Summary panel:                                              │
│   Mean pLDDT across filtered: 85                            │
│   pLDDT distribution: 40% >90, 30% 70–90, 30% 50–70         │
│   Use cases enabled: drug docking, structural homology       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 12 / 50 entries pass filters (high-quality).           │   │
│ │ Of the excluded: 30% disordered, 20% shallow MSA       │   │
│ │   (orphan), 10% low confidence for other reasons.      │   │
│ │ Relax filters to see the full distribution.             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Min mean pLDDT slider.
- Max disorder fraction slider.
- Length range slider.
- Function-category dropdown.

### What they see

- Grid of thumbnails matching current filter.
- Summary statistics.
- Outcome banner describing filtered subset.

### Target aha

The AlphaFold Database is a searchable catalogue. pLDDT-aware filtering isolates trustworthy subsets. Different protein classes distribute differently across pLDDT — orphans and disordered proteins cluster at low confidence; well-annotated enzymes cluster high.

### Technical notes

- Pure JS.
- ~50 preset entries; each with metadata (UniProt-like ID, length, mean pLDDT, disorder %, function).
- 2D cartoon thumbnails: pre-rendered SVG with pLDDT-colour rainbow.
- Filter logic: straightforward range + category intersection.

### Acceptance criteria

- [ ] Default filter returns a reasonable subset (~10–15 entries).
- [ ] Thumbnails render clearly.
- [ ] Filter changes update the grid live.
- [ ] Summary panel accurately reflects the filtered set.
- [ ] Opens pre-rendered with default filters.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Protein Structure Level Explorer

**File**: `artifacts/lecture-15/07-structure-levels.html`
**Lecture anchor**: §1.2 Primary, secondary, tertiary, quaternary
**EE framing reinforced**: hierarchy of representations; different analyses operate at different levels.

### Teaching purpose

Show a single protein at each structural level with cross-linked highlighting. Click a residue in any level → it highlights across all levels. Reinforce that the levels are the same molecule seen through different abstractions.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Protein: [ haemoglobin α-subunit ▾ / myoglobin / ribonuclease]│
├─────────────────────────────────────────────────────────────┤
│ Primary (sequence strip):                                   │
│   V H L T P E E K S A V T A L W G K V N V D E ...           │
├─────────────────────────────────────────────────────────────┤
│ Secondary (annotated sequence):                             │
│   V H L   | helix 1 | ... | strand 1 | ... | helix 3 | ...  │
│   (sequence letters under their structure annotation)       │
├─────────────────────────────────────────────────────────────┤
│ Tertiary (2D cartoon of folded single chain):               │
│   helices as ribbons, strands as arrows, loops smooth       │
├─────────────────────────────────────────────────────────────┤
│ Quaternary (assembled tetramer with α,α,β,β subunits):      │
│   the selected chain highlighted; others faded              │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Click residue 58 (red highlight propagates through all │   │
│ │ four levels). In haemoglobin α, residue 58 is the      │   │
│ │ proximal histidine binding the heme iron — visible as   │   │
│ │ a sphere in the tertiary view.                          │   │
│ │ Switch to myoglobin → similar architecture, different   │   │
│ │ quaternary (monomer).                                   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Protein-preset dropdown.
- Residue-index input or click-to-select on any level.

### What they see

- Four stacked panels, one per structural level.
- Cross-level highlight on the selected residue.
- Outcome banner with biological context.

### Target aha

A protein isn't one-level; it's a hierarchy of representations. Different analyses operate at different levels: sequence alignments at primary, secondary-structure prediction at secondary, docking at tertiary, assembly interactions at quaternary. Selecting a residue and watching it propagate across levels teaches the hierarchical unity.

### Technical notes

- Pure JS.
- ~3 preset proteins with hardcoded sequence, secondary-structure annotations (helix/strand/loop per residue), 2D residue positions for tertiary, quaternary subunit arrangement.
- Cross-level highlight: selected residue index shared across panels.
- Residue-click event handler on each panel.

### Acceptance criteria

- [ ] Default (haemoglobin α) renders all four levels.
- [ ] Clicking a residue in any panel highlights it in all others.
- [ ] Preset switching updates all four panels consistently.
- [ ] Banner displays biologically correct note for selected residue.
- [ ] Opens pre-rendered with a default residue selected.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Amino-acid letter colouring consistent across artifacts (§1.2 palette).
- pLDDT colour bands identical across artifacts #4, #6.
- Backbone cartoons use consistent 2D conventions (helices as cylinders, strands as arrows, loops smooth).
- Every artifact emits an **outcome banner** per convention §1.4.
- KaTeX (where used) loaded from the same CDN tag.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls beyond KaTeX (where used).
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
