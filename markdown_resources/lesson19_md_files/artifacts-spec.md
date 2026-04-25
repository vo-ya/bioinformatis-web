# Lecture 19 (proposed L3) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 19 (BLAST and Sequence Search Statistics).
> **Companion files**: `lecture-19.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-19/NN-name.html`.
- No build step. Vanilla HTML / CSS / JavaScript.
- Tokens via `../_shared/artifact-theme.css`.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Substitution matrix colours: positive scores green, neutral grey, negative red.
- Algorithm-stage colours: seed cobalt, extension amber, scoring red.
- Typography: Inter for chrome; JetBrains Mono for sequences, k-mers, scores, E-values.
- Default state is instructive: opens with a pre-computed example.
- Every artifact emits an outcome banner with the named aha verdict.
- Every artifact has an "Educational tool — not a replacement for production NCBI BLAST" disclaimer.

## Artifact budget

Seven interactive tools.

| # | Title | Lecture anchor |
|---|---|---|
| 1 | Substitution-matrix score calculator | §2.4 |
| 2 | Seed neighbourhood explorer | §3.1 |
| 3 | E-value calculator | §4.2 |
| 4 | X-drop extension simulator | §3.3 |
| 5 | PSI-BLAST iteration explorer | §5.2 |
| 6 | BLAST-based function transfer | §7.3 |
| 7 | DIAMOND vs BLASTP runtime racer | §6.3 |

---

## Artifact #1 — Substitution-Matrix Score Calculator

**File**: `artifacts/lecture-19/01-substitution-score.html`
**Anchor**: §2.4

### Teaching purpose

Pick two amino acids; the artifact returns the score under each of BLOSUM62 / BLOSUM45 / BLOSUM80 / PAM30 / PAM250 and explains the score's log-odds interpretation.

### UI layout

- Two amino-acid dropdowns (residue 1 × residue 2).
- One matrix dropdown.
- Numeric score readout with sign and verbal interpretation (positive → match-prone, negative → mismatch-prone).
- Heat-map sidebar showing the full matrix's row for the selected residue 1.
- Outcome banner: "score interpretation" (e.g., "BLOSUM62(D,E) = 2 → D and E are conservatively substituted (acidic side chains).").

### Target aha

Substitution matrices are log-odds; positive score means residue pair is observed more than chance, encoding biochemical similarity.

### Acceptance criteria

- All 20 × 20 entries available across 5 matrices.
- BLOSUM62(D,E) → 2; BLOSUM62(W,W) → 11; PAM30(W,W) → 13.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — Seed Neighbourhood Explorer

**File**: `artifacts/lecture-19/02-seed-neighbourhood.html`
**Anchor**: §3.1

### Teaching purpose

Enter a 3-mer (or pick a default). Slide the threshold T. The artifact shows all 3-mers in the neighbourhood (BLOSUM62-score ≥ T with the query 3-mer).

### UI layout

- Input: 3-mer text box (validated to 3 amino acids).
- Slider: threshold T (range 5–20).
- Output: scrollable list of neighbourhood 3-mers, sorted by score.
- Count badge: "47 of 8000 possible 3-mers in neighbourhood".
- Outcome banner: "as T decreases, neighbourhood grows; sensitivity up, speed down".

### Target aha

The threshold T controls a sensitivity-speed trade-off. Lower T → broader seeds → more extensions → higher recall but slower.

### Acceptance criteria

- Default 3-mer "KVL", T = 11 → ~50 neighbours.
- T = 5 → ~500 neighbours.
- T = 20 → mostly self-only (KVL itself).
- Exposes the size scaling visually (neighbourhood size on a log scale).
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — E-value Calculator

**File**: `artifacts/lecture-19/03-evalue-calculator.html`
**Anchor**: §4.2

### Teaching purpose

Slide bit score, query length, and database size. The artifact returns E and P, plus a verdict (significant / borderline / not significant).

### UI layout

- Slider 1: bit score S' (range 20–200).
- Slider 2: query length m (range 50–1000).
- Slider 3: database size n (slider with markers at SwissProt, UniRef90, nr, UniProt full).
- Outputs: E-value, P-value, "significance" verdict.
- Plot: E-value vs S' for the chosen (m, n).
- Outcome banner: "at S' = X with database = SwissProt, E = ... → ..."

### Target aha

E-value scales linearly in database × query, exponentially in bit score. Understanding this scaling tells you when to trust a hit.

### Acceptance criteria

- m=100, n=10⁹, S'=30 → E ≈ 100 (not significant).
- m=100, n=10⁹, S'=50 → E ≈ 10⁻⁴ (significant).
- m=100, n=10⁹, S'=100 → E ≈ 10⁻¹⁹ (gold standard).
- Database-size slider auto-updates E with linear scaling.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — X-Drop Extension Simulator

**File**: `artifacts/lecture-19/04-xdrop-extension.html`
**Anchor**: §3.3

### Teaching purpose

Simulate ungapped extension on a synthetic database sequence with embedded homology. The X-drop trajectory is shown live; user controls X.

### UI layout

- Pre-loaded query and subject sequences with one true homologous region.
- Slider: X-drop threshold (range 5–50).
- Visualisation: sequence alignment with running score at each position; coloured trajectory; X-drop boundary marked.
- Statistics: HSP length and score recovered at this X.
- Outcome banner: "X = 5 cuts extension early → recover only short HSPs. X = 50 over-extends into noise → false long HSPs."

### Target aha

X is a sensitivity-specificity knob. Too small → miss true extensions. Too large → false-positive extensions.

### Acceptance criteria

- Default X = 22; recovers a 50 aa HSP with bit score 60.
- X = 5; recovers a 15 aa HSP with bit score 25.
- X = 50; over-extends into noise; final HSP has score 70 with low identity.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — PSI-BLAST Iteration Explorer

**File**: `artifacts/lecture-19/05-psi-blast.html`
**Anchor**: §5.2

### Teaching purpose

Walk through 5 iterations of PSI-BLAST on a small toy database. Visualize how the PSSM and hit list evolve.

### UI layout

- Iteration counter (current / 5).
- Step-forward button.
- PSSM heat map (residues × positions) updates each iteration.
- Hit list sortable by E-value, with new hits highlighted.
- "Drift detection" indicator: shown if a new hit has low identity to the original query.
- Outcome banner: "after 3 iterations, profile sensitivity is improved; one drift candidate detected at iteration 2".

### Target aha

PSI-BLAST is iterative profile refinement; sensitivity improves at the cost of drift risk; cap at 3 iterations.

### Acceptance criteria

- 5 iterations defined in advance with realistic-looking PSSMs.
- Drift introduces at iteration 4 in default scenario.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — BLAST-Based Function Transfer

**File**: `artifacts/lecture-19/06-function-transfer.html`
**Anchor**: §7.3

### Teaching purpose

Pick a query protein from preset list. Artifact shows simulated BLAST hits, then offers transferred function annotations (GO terms) from the top hits — with confidence scores based on E-value and identity.

### UI layout

- Query protein dropdown (5 presets across kinase, transcription factor, enzyme, structural, signalling).
- Hits table: subject ID, E-value, identity %, source organism, GO terms.
- Transferred annotation panel: aggregated GO terms with confidence.
- Outcome banner: "function transfer at 60% identity → high confidence; below 30% → unreliable for function".

### Target aha

Function transfer is a chain of inferred-from-inferred-from. Each step contributes uncertainty; track it.

### Acceptance criteria

- 5 preset queries each with realistic-looking hit lists.
- Confidence scoring: identity > 50% AND E < 10⁻¹⁰ → high; identity > 30% → medium; below → low.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — DIAMOND vs BLASTP Runtime Racer

**File**: `artifacts/lecture-19/07-diamond-blast-race.html`
**Anchor**: §6.3

### Teaching purpose

Show side-by-side simulated runtime + sensitivity for BLASTP, DIAMOND --sensitive, DIAMOND --fast, MMseqs2, USEARCH on a chosen query length and database size.

### UI layout

- Slider: query length (50–10,000 aa).
- Slider: database size (10⁶–10¹² residues).
- Bar chart: runtime (log-y) for each tool.
- Sensitivity readout per tool.
- Outcome banner: "for query length X against database Y, DIAMOND --sensitive is the right balance".

### Target aha

Choice depends on the workload: small DB → BLASTP fine; large DB → DIAMOND or MMseqs2; high sensitivity needed → BLASTP or MMseqs2.

### Acceptance criteria

- Realistic runtime estimates calibrated against published benchmarks.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Substitution-matrix colour conventions identical across #1, #2, #4.
- Every artifact emits an outcome banner per convention.
- All artifacts explicitly state they are educational and not clinical-grade or research-grade.

## Testing checklist (per artifact)

- [ ] Opens standalone in browser, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner / verdict line visible at end of any user interaction.
