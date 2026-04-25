# Lecture 19 (proposed L3) — Figures Specification

> **Scope**: Static diagrams for Lecture 19 (BLAST and Sequence Search Statistics).
> **Companion files**: `lecture-19.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Substitution-matrix colours: positive scores green `#059669`, neutral grey `#525252`, negative red `#c4342c`.
- Algorithm-stage colours: seed cobalt `#1e40af`, extension amber `#b45309`, scoring red `#c4342c`.
- Speed-vs-sensitivity Pareto: BLAST classical cobalt, DIAMOND amber, MMseqs2 teal, USEARCH violet.
- Typography: Inter for labels; JetBrains Mono for sequences, k-mers, scores, E-values.

## Figure budget

Twelve figures for a ~3h 30min lecture.

| # | Title | Part |
|---|---|---|
| 1 | BLAST search architecture | Part 1 |
| 2 | Seed-extend-score worked example | Part 3 |
| 3 | Karlin-Altschul: bit score → E-value | Part 4 |
| 4 | PSI-BLAST iteration as EM | Part 5 |
| 5 | Search-speed-vs-sensitivity Pareto | Part 6 |
| 6 | Substitution matrix heat map (BLOSUM62) | Part 2 |
| 7 | PAM vs BLOSUM use-case map | Part 2 |
| 8 | Two-hit method: spacing constraint | Part 3 |
| 9 | X-drop extension trajectory | Part 3 |
| 10 | E-value scaling with database size | Part 4 |
| 11 | Reciprocal best hits ortholog identification | Part 7 |
| 12 | BLAST output format dissection | Part 4 |

---

## Figure 1 — BLAST search architecture

**File**: `diagrams/lecture-19/01-blast-pipeline.svg`
**ViewBox**: `0 0 1200 600`

Top-to-bottom flowchart with annotation. Stages:

1. Query sequence (single peptide / nucleotide).
2. Word-list construction: k-mer enumeration + neighbourhood expansion (T threshold).
3. Database scan: hash-table lookups; counts of database hits per query position.
4. Two-hit pairing: spacing within $A$ residues triggers extension.
5. Ungapped extension: X-drop boundaries, HSP scoring.
6. Gapped extension: bounded Smith-Waterman.
7. Score statistics: bit-score and E-value.
8. Significance filter: $E < E_{thresh}$ → reported hits.

Annotations on each stage with computational complexity. Sidebar callout: "BLAST trades exactness for tractability — and recovers statistical interpretability."

---

## Figure 2 — Seed-extend-score worked example

**File**: `diagrams/lecture-19/02-seed-extend.svg`
**ViewBox**: `0 0 1200 540`

Three-panel walkthrough:

- Panel 1: Query peptide "MKVLWAGCERFPN" with one 3-mer "KVL" highlighted; database scan finds match at position 47 of subject "...AGCKVLWAQ...".
- Panel 2: Neighbourhood of "KVL" — three matrix-equivalent 3-mers ("KVI", "RVL", "KVM") also match in subject; visual neighbourhood radius shown as a circle in scoring space.
- Panel 3: Ungapped extension from the seed in both directions; running score curve; X-drop stop boundaries; final HSP highlighted.

---

## Figure 3 — Karlin-Altschul: bit score → E-value

**File**: `diagrams/lecture-19/03-evalue-curve.svg`
**ViewBox**: `0 0 1200 600`

Two-panel figure:

- Panel 1: Gumbel (extreme-value) distribution of random-database bit scores, with quantile thresholds marked at S' = 30, 50, 100. Annotation: "tail decays as $e^{-\lambda S}$".
- Panel 2: Log-y plot of E-value vs bit score for three database sizes ($10^6$, $10^9$, $10^{12}$). Horizontal reference lines at E = 1, 0.001, $10^{-10}$, $10^{-50}$. Annotations on canonical thresholds.

---

## Figure 4 — PSI-BLAST iteration as EM

**File**: `diagrams/lecture-19/04-psi-blast.svg`
**ViewBox**: `0 0 1200 540`

Circular iteration diagram:

- Iteration 0: query → BLAST → initial hits.
- Iteration 1: hits → MSA → PSSM → BLAST.
- Iteration 2-5: more hits → updated PSSM → more hits.
- Convergence label.
- Side panel: PSSM heat map (20 amino acids × 100 positions) showing conserved positions as bright stripes.
- Bottom annotation: "EM-flavoured profile refinement; cap at 3 iterations to limit drift."

---

## Figure 5 — Search-speed-vs-sensitivity Pareto

**File**: `diagrams/lecture-19/05-pareto.svg`
**ViewBox**: `0 0 1200 600`

2D scatter:

- Y-axis: sensitivity (recall at 25% identity), 0–1.
- X-axis: queries-per-second (log scale), $10^0$ to $10^4$.
- Points: BLASTP classical (10⁰ qps, 0.95 sens), DIAMOND --sensitive (10² qps, 0.92), DIAMOND --fast (10³ qps, 0.85), MMseqs2 (10² qps, 0.96), USEARCH (10⁴ qps, 0.65).
- Dashed Pareto frontier line.
- Annotation: "MMseqs2 sits on the Pareto frontier across most query types; pick by what you can tolerate to lose."

---

## Figure 6 — BLOSUM62 substitution matrix heat map

**File**: `diagrams/lecture-19/06-blosum62.svg`
**ViewBox**: `0 0 1000 1000`

20 × 20 colour-coded heat map of BLOSUM62 entries:

- Rows: amino acids in standard order (ARNDCQEGHILKMFPSTWYV).
- Columns: same.
- Colour scale: positive scores green to dark green; zero white; negative scores light red to dark red.
- Numeric value in each cell (font-mono).
- Side legend with property groupings (hydrophobic, polar, charged, etc.).

---

## Figure 7 — PAM vs BLOSUM use-case map

**File**: `diagrams/lecture-19/07-pam-blosum.svg`
**ViewBox**: `0 0 1200 480`

Horizontal axis: percent identity of expected homologs, 100% → 10%.

Multiple horizontal bars showing applicability range of PAM30, PAM70, PAM250, BLOSUM45, BLOSUM62, BLOSUM80, BLOSUM90:

- PAM30 / BLOSUM90: high-identity (close homologs).
- BLOSUM62: moderate (default — typical search).
- PAM250 / BLOSUM45: distant homologs.

Annotations on which scenarios prefer which matrix. Bottom annotation: "BLOSUM62 is the universal default; PAM250 is the legacy classical."

---

## Figure 8 — Two-hit method: spacing constraint

**File**: `diagrams/lecture-19/08-two-hit.svg`
**ViewBox**: `0 0 1200 480`

Sketch a database sequence with various seed hits along its length. Annotated: hits A, B, C, D positioned at residues 100, 145, 200, 350.

- Pair A–B (distance 45): satisfies spacing constraint A ≤ 40 → triggers extension.
- Pair C alone (distance to nearest neighbour > 40): no extension triggered.
- Pair B–C: satisfies spacing → extension.
- Result: extensions on regions with multiple close seeds, reducing false-positive extensions on isolated noise hits.

Annotation: "the two-hit constraint cuts extension cost ~10× while preserving sensitivity for true homologs (which have multiple close seeds)."

---

## Figure 9 — X-drop extension trajectory

**File**: `diagrams/lecture-19/09-x-drop.svg`
**ViewBox**: `0 0 1200 480`

X-axis: residue position relative to seed.
Y-axis: cumulative score during extension.

Three score trajectories:

- Trajectory 1 (true homology): score climbs steadily, reaches 90 at right boundary; X-drop never triggered.
- Trajectory 2 (mid-strength): score climbs to 60 then drifts; X-drop triggers when score drops by X = 22 from peak; extension stops.
- Trajectory 3 (noise): score wobbles around 30; X-drop triggers immediately.

Annotation: "X-drop stops extension when the running score drops X below current peak — the matched-filter analog of leaving the peak region."

---

## Figure 10 — E-value scaling with database size

**File**: `diagrams/lecture-19/10-database-scaling.svg`
**ViewBox**: `0 0 1200 540`

A line plot showing E-value at fixed bit score (S' = 50) as database size grows from $10^6$ to $10^{12}$ residues:

- The E-value increases linearly in $n$ (since $E \propto m \cdot n$).
- Visual examples: SwissProt (~$2 \cdot 10^8$ residues), UniRef90 (~$10^{11}$), nr (~$10^{11}$), full UniProt (~$10^{12}$).
- Reference databases marked with vertical lines.
- Annotation: "Doubling the database doubles the E-value at any fixed bit score."

---

## Figure 11 — Reciprocal best hits for ortholog identification

**File**: `diagrams/lecture-19/11-rbh.svg`
**ViewBox**: `0 0 1200 480`

Two species' protein sets shown as columns of nodes.

- Forward search: gene A (species 1) → BLAST against species 2 → best hit gene B.
- Reverse search: gene B (species 2) → BLAST against species 1 → best hit gene A.
- → A and B are reciprocal best hits = orthologs.

Counter-example: gene A → best hit X in species 2; X → best hit C (≠ A) → not RBH.

Annotation: "RBH catches the ~80% of orthologs that are well-separated; falls down on recent paralogs (use OrthoFinder for those)."

---

## Figure 12 — BLAST output format dissection

**File**: `diagrams/lecture-19/12-blast-output.svg`
**ViewBox**: `0 0 1200 720`

A typical BLAST output snippet, with each line annotated:

```
> tr|A0A024R7T9|A0A024R7T9_HUMAN  Some Protein
Length=247
Score = 187 bits (475),  Expect = 4e-49, Method: Composition-based stats.
Identities = 95/120 (79%),  Positives = 102/120 (85%),  Gaps = 0/120 (0%)

Query  10   MKLVWA...
            MKLVWA
Sbjct  47   MKLVWA...
```

Each piece annotated:

- "187 bits": bit score (matrix-independent).
- "(475)": raw score (matrix-dependent).
- "Expect = 4e-49": E-value.
- "Method: Composition-based stats": local-composition correction in use.
- "Identities = 95/120": exact-match count and proportion.
- "Positives = 102/120": positive-substitution count (D↔E counts even though not identical).
- Query / Sbjct alignment: query position offset, subject position offset, residue alignment.

Side panel: "What to ignore" (low-Identities high-bit hits where the alignment is in a low-complexity region).
