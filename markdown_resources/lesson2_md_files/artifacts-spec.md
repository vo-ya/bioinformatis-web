# Lecture 2 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 2 (Read Alignment).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below only override them when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-02/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified in the per-artifact section.
- The file must render standalone when opened directly in a browser — no server, no CORS gymnastics.
- Artifact is embedded in the lecture page via `<iframe>` loaded lazily when it scrolls near the viewport.

### 1.2 Visual design

- Use the design tokens from `diagram-style-guide.md` §3. Inline them as CSS custom properties at the top of each file.
- DNA bases, when rendered, use the `--base-*` palette exactly. No substitutes.
- Typography: **Inter** for UI chrome, **JetBrains Mono** for sequences and algorithmic state (suffix arrays, BWT strings, CIGAR strings, matrix cells). Same as the diagrams.
- Default state is instructive: the artifact opens showing a meaningful example, not a blank canvas. A student who never touches a control should still see the teaching point.
- Controls live in a strip above or left of the visualization. The visualization gets the bulk of the real estate.
- No animations longer than ~400 ms. No decorative transitions. Motion only when it carries information (e.g., watching an interval contract, a pointer advance).

### 1.3 Interaction model

Artifacts follow a consistent control vocabulary:

- **Input string / query** — editable text field, monospace, validated against the alphabet (`ACGT$` unless noted).
- **Step** — advance one algorithmic step.
- **Play / Pause** — run through steps at ~1 step/sec.
- **Reset** — return to the initial example.
- **Speed** — optional slider, 0.25×–4×. Default 1×.

Any illegal input shows a quiet inline message (`--fg-muted`), not a modal. Never a `window.alert`.

### 1.4 Pedagogical constraint

Every artifact must produce a **specific realization** — the "target aha moment" in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed and should be revised. This is the acceptance test, not a feature count.

### 1.5 Out of scope

- No logins, no accounts, no persistence between sessions.
- No telemetry.
- No analytics.
- No loading of external data files larger than ~50 KB. Reference strings are baked in as string literals.

---

## 2. Artifact #1 — Brute-Force Aligner

**File**: `artifacts/lecture-02/01-brute-force.html`
**Lecture anchor**: §2.1 Brute force and its complexity
**EE framing reinforced**: sliding-window correlation; the O(|R|·|r|) cost is made visceral via a live comparison counter.

### Teaching purpose

Establish the baseline. Every subsequent algorithm in the lecture must beat this one. The student should leave with a concrete number of character comparisons burned into their memory — something like "1.2 million comparisons for a 1 kb reference and a 100 bp query" — so the later improvements feel like real wins, not abstract Big-O deltas.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Reference:  [AAGCTCAGTCGATCGATCGATCAGTCGATCGATCAGTCGA...]  │ <- editable, 200–2000 chars
│ Query:      [ATCGATCAG]                                     │ <- editable, 5–50 chars
│ [Reset] [Step] [Play ▶]   Speed: (slider)                   │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  AAGCTCAGTCGATCGATCGATCAGTCGATCGATCAGTCGA...                │
│       ATCGATCAG                                             │   <- query slides left→right
│       ^  ^                                                  │   <- cell-by-cell compare
│                                                             │
├─────────────────────────────────────────────────────────────┤
│ Position: 12 / 984   Comparisons: 847   Matches found: 3    │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Edit the reference (length cap 2000 characters so the comparison count stays comprehensible).
- Edit the query.
- Step / play / reset.
- (Optional) toggle "early exit on mismatch" — the realistic speedup that still doesn't change the Big-O.

### What they see

- The query slides along the reference one position at a time.
- At each position, each compared character flashes briefly: green if match, red-muted if mismatch.
- A counter tallies **total character comparisons**, not just positions checked. This is the metric that matters.
- Matches are logged in a small results panel with their starting positions.

### Target aha moment

"A 10 kb reference and a 50 bp read takes half a million comparisons, and the human genome is 300,000× bigger than that reference." Optionally, show an estimated wall-clock time extrapolated to 3 Gb — it runs into hours per read.

### Technical notes

- Pure JavaScript, no dependencies.
- Comparison visualization throttled to ≤60 fps; at high speed, batch updates.
- Cap reference at 2000 chars to keep the animation legible; mention the cap in a hint.

### Acceptance criteria

- [ ] Default input demonstrates ≥ 3 match positions.
- [ ] Comparison counter visibly climbs during Play.
- [ ] Stepping one position at a time works correctly at the reference boundary.
- [ ] Legible at 720px artifact width.

---

## 3. Artifact #2 — Suffix Array Explorer

**File**: `artifacts/lecture-02/02-suffix-array.html`
**Lecture anchor**: §2.2 Suffix arrays and binary search
**EE framing reinforced**: sorted lookup table + successive approximation (the suffix array is the LUT; binary search is the SAR ADC).

### Teaching purpose

Make the suffix array concrete — not as a concept but as a tangible data structure the student has seen get built. Then show binary search as the interval-narrowing move that turns O(|R|) into O(log |R|).

### UI layout

Two linked panels.

```
┌────────────────────────────┬──────────────────────────────────┐
│  Build                     │  Query                           │
│  Input: [BANANA$______]    │  Pattern: [ANA___]               │
│  [Build SA]                │  [Search] [Step]                 │
├────────────────────────────┼──────────────────────────────────┤
│   i  SA[i]  suffix         │  Binary search trace:            │
│   0    6    $              │    lo=0  hi=6  mid=3             │
│   1    5    A$             │       compare "ANA" vs "ANA$"    │
│   2    3    ANA$           │       ... interval contracts     │
│   3    1    ANANA$         │                                  │
│   4    0    BANANA$        │  Matching range: [2, 3]          │
│   5    4    NA$            │  → positions in R: {3, 1}        │
│   6    2    NANA$          │                                  │
└────────────────────────────┴──────────────────────────────────┘
```

### Student controls

- Left panel: type an input string up to ~30 characters, press Build SA to see the rotation list sort itself into the suffix array (a brief animated sort is a win).
- Right panel: type a query pattern and run binary search. Step mode advances `lo`, `mid`, `hi` one comparison at a time.

### What they see

- The suffix array builds visibly by sorting the list of suffixes (a short animated insertion or merge sort is fine).
- Binary search highlights three rows (`lo`, `mid`, `hi`). At each step, the compared suffix's prefix is compared to the pattern character-by-character; the `lo`/`hi` pointers jump to the new interval.
- The final answer highlights the contiguous row range whose suffixes start with the query, and maps them back to positions in the original string.

### Target aha moment

"The suffix array sorts all the places the read could start. Binary search on that sorted list is log-time. The read is a *prefix of a range* of suffixes." The last sentence is the key phrase — if the student can explain it back, the artifact worked.

### Technical notes

- Cap input at 30 chars so the rotation list fits on screen.
- Sort visualization is pedagogical, not performance-critical — use a simple sort with a short delay per swap.
- Include the `$` terminator — students should see it and understand why.

### Acceptance criteria

- [ ] Building BANANA$ produces the canonical SA `[6,5,3,1,0,4,2]`.
- [ ] Binary search on "ANA" in BANANA$ returns the range [2,3], positions {3,1}.
- [ ] Invalid pattern (not a substring) ends with an empty range and a muted "not found" note.
- [ ] Pointer highlights use the accent color plus a weight change (not color alone).

---

## 4. Artifact #3 — BWT Builder

**File**: `artifacts/lecture-02/03-bwt-builder.html`
**Lecture anchor**: §3.1 The Burrows-Wheeler Transform
**EE framing reinforced**: reversible transform moving data to a domain where the target operations are cheap — same structural move as FFT/DCT.

### Teaching purpose

Turn BWT from an incantation ("sort rotations, take the last column") into something the student has watched happen. The rotation matrix is the full object; the BWT is the projection onto one column. Once they've seen the rotations line up alphabetically, the LF property (next artifact) stops being magic.

### UI layout

```
┌────────────────────────────────────────────────────────────┐
│ Input string: [ACAACG$___]  [Build BWT]                    │
├────────────────────────────────────────────────────────────┤
│   Rotations (unsorted)     →   Rotations (sorted)          │
│                                                            │
│   ACAACG$                        $ACAACG                   │
│   CAACG$A                        AACG$AC                   │
│   AACG$AC                        ACAACG$                   │
│   ACG$ACA                        ACG$ACA                   │
│   CG$ACAA                        CAACG$A                   │
│   G$ACAAC                        CG$ACAA                   │
│   $ACAACG                        G$ACAAC                   │
│                                     ↑F              ↑L     │
│                                                            │
│   BWT = L column = "GC$AAAC"                               │
└────────────────────────────────────────────────────────────┘
```

### Student controls

- Input string up to ~12 characters, alphabet `ACGT`. The `$` terminator is appended automatically and visually flagged.
- Build BWT button animates the sort.
- Hover a row in either panel to see it highlighted in the other.

### What they see

- The raw rotation matrix on the left.
- The sorted rotation matrix on the right, built by a short animated sort.
- F (first) and L (last) columns flagged with small uppercase kicker labels ("F", "L") above them, in the accent color.
- The BWT string rendered below, each character colored to match its row, so the student sees where each output character came from.

### Target aha moment

"Sorting the rotations puts strings with the same following-context next to each other. The last column then holds characters with similar right-contexts clustered — which is why BWT compresses well and also why we can search it." The compression intuition and the search intuition are the same observation at different scales.

### Technical notes

- Input alphabet restricted to `ACGT` plus the implicit `$`.
- Cap length at 12 so all rotations fit visibly.
- Color each rotation row by a consistent hue per starting-position so hovering makes the F↔L relationship visible.

### Acceptance criteria

- [ ] BWT of "ACAACG$" is "GC$AAAC" (verify with a worked example).
- [ ] Hover on a sorted row highlights the same string in the unsorted panel.
- [ ] F and L columns are visually flagged.
- [ ] No dependencies beyond vanilla JS.

---

## 5. Artifact #4 — FM-Index Backward Search

**File**: `artifacts/lecture-02/04-fm-backward-search.html`
**Lecture anchor**: §3.3 Backward search: the FM index
**EE framing reinforced**: operating in the transformed (BWT) domain — you never reconstruct the original reference; you update interval bounds using rank queries, like successive downconversion in the frequency domain.

### Teaching purpose

Backward search is the single most counterintuitive move in the lecture: we process the query **right-to-left** and maintain an interval `[sp, ep]` in the suffix-array space using only `rank` queries on L. If the student can step through `GATTACA` once and watch the interval shrink, the FM index stops being folklore.

### UI layout

Three regions:

1. **Top**: the reference string (e.g., `GATTACAGATTACAT$`), its BWT, and the F column. A small C-table (cumulative count) sits under F.
2. **Middle**: the query as a row of boxes, processed right-to-left. The current character glows in the accent color.
3. **Bottom**: a live display of `sp`, `ep`, and the matching SA rows highlighted in the F column. A history trace shows every (character, sp, ep) tuple so far.

### Student controls

- Pick a preset reference (3 canonical examples) or enter a short one (≤ 20 chars).
- Pick / type a query.
- Step backward one character at a time, or Play.

### What they see

- The current character highlighted in the query.
- The matrix identity on screen: `sp_new = C[c] + rank(c, sp)` and `ep_new = C[c] + rank(c, ep+1) − 1`, with the values for this step substituted live.
- The F-column interval `[sp, ep]` highlighted as a band; when the query is consumed, the band shows all occurrences (mapped to SA positions and shown on the reference).
- If the interval collapses to empty at any step, the query is declared not found, and the step that killed it is flagged — diagnostically useful.

### Target aha moment

"The reference isn't in memory in any searchable form — only the BWT and the rank arrays are. The interval `[sp, ep]` is the entire search state. Each character of the query shrinks it in O(1)." The student should be able to articulate: *we never look at the reference itself during search.*

### Technical notes

- Reference cap: 20 chars so the SA, BWT, F column, and C-table all fit on one screen.
- Precompute full rank arrays for this artifact (checkpoints are artifact-level overkill; saved for the next section's diagram).
- The (sp, ep) formula displayed should update its values character-by-character, not just once.

### Acceptance criteria

- [ ] Searching `ACA` in `ACAACG$` ends with exactly one occurrence at position 2.
- [ ] Searching a pattern absent from the reference flags the step where the interval collapsed.
- [ ] Each backward step shows the substituted formula with actual numbers.
- [ ] Step direction is visibly right-to-left; no ambiguity.

---

## 6. Artifact #5 — Smith-Waterman Matrix

**File**: `artifacts/lecture-02/05-smith-waterman.html`
**Lecture anchor**: §4.3 Smith-Waterman
**EE framing reinforced**: dynamic programming on a 2D trellis — direct analog of Viterbi decoding; the traceback is the argmax path.

### Teaching purpose

Make the Smith-Waterman matrix fully legible. Every cell is a max over four candidates; every cell's value has a reason. The traceback is a path, not a mystery. This artifact is the single most important one for the lecture's Viterbi framing.

### UI layout

```
┌───────────────────────────────────────────────────────────────┐
│ Reference: [GATTACA_____]   Query: [GCATGCA____]              │
│ Match: +2   Mismatch: -1   Gap: -2                            │
│ [Fill matrix ▶]  [Step-fill]   [Traceback]                    │
├───────────────────────────────────────────────────────────────┤
│        ''   G   A   T   T   A   C   A                         │
│    ''   0   0   0   0   0   0   0   0                         │
│    G    0   2   0   0   0   0   0   0                         │
│    C    0   0   1   0   0   0   2   0                         │
│    A    0   0   2   0   0   2   0   4                         │
│    T    0   ...                                               │
│    ...                                                        │
│                                                               │
│ Score = 7 at (3, 7).  Traceback:                              │
│   GCAT-GCA                                                    │
│   G-ATTACA  → CIGAR = 1=1D2=1I3=                              │
└───────────────────────────────────────────────────────────────┘
```

### Student controls

- Edit reference and query (≤ 15 chars each, to keep the matrix readable).
- Adjust scoring parameters: match reward, mismatch penalty, gap penalty (linear for v1; affine is a stretch goal).
- Three fill modes: full-fill, step-fill (one cell at a time), or paint-fill (Play).
- Traceback button highlights the optimal path from the matrix max.

### What they see

- Matrix cells filled with values; hovering a cell shows the four candidates being maxed over (`diag + score(i,j)`, `up + gap`, `left + gap`, `0`), with the winning candidate bolded.
- The global maximum cell highlighted.
- The traceback path drawn as a contiguous sequence of arrows from the max back to a zero.
- Below the matrix, the alignment rendered on two lines with gaps, and the CIGAR string derived from the path.

### Target aha moment

Two, ideally both:
1. "Every number in this matrix is the best score of a local alignment ending at that cell. The `0` floor is what makes it local rather than global."
2. "This is the same algorithm as Viterbi. The trellis is 2D because we have two sequences. Traceback is argmax-path." If the student has seen Viterbi, the connection should land the moment they watch the traceback.

### Technical notes

- Sequence length cap: 15 chars per axis. A 15×15 matrix is the comfort zone.
- Cell hover card shows the four candidates with arrows pointing in the direction each came from.
- Traceback uses a distinct visual treatment (accent color + thicker stroke) — and also arrow glyphs, so it's not color-alone.

### Acceptance criteria

- [ ] Default example produces a non-trivial alignment with at least one indel.
- [ ] Changing the gap penalty visibly changes the traceback.
- [ ] Hover on any filled cell reveals the four-way max derivation.
- [ ] CIGAR string generation matches the traceback exactly.

---

## 7. Artifact #6 — Seed-and-Extend Visualizer

**File**: `artifacts/lecture-02/06-seed-and-extend.html`
**Lecture anchor**: §5.2 Seed-and-extend in real aligners
**EE framing reinforced**: coarse acquisition + fine tracking. The FM index is the coarse detector; Smith-Waterman on small windows is the fine estimator. Classical hierarchical search.

### Teaching purpose

Tie the whole lecture together. The student should see the FM index from §3 produce seeds, then Smith-Waterman from §4 extend those seeds. If the previous five artifacts are individual instruments, this one is the orchestra.

### UI layout

```
┌──────────────────────────────────────────────────────────────────┐
│ Reference:  [  ~400 bp pre-loaded example  ]                     │
│ Read:       [  50 bp read with one SNV and one 2bp indel  ]      │
│ Seed length k: (slider, 8–20, default 12)                        │
│ [Find seeds] → [Extend] → [Report alignment]                     │
├──────────────────────────────────────────────────────────────────┤
│ Reference, zoomed:                                               │
│  ...ATCGAT|seed1|TTAGC...ATCG|seed2|CCGTA...|seed3|...          │
│           ↑ seed at pos 42 (FM interval size 1)                  │
│           ↑ seed at pos 217 (FM interval size 5, ambiguous)      │
├──────────────────────────────────────────────────────────────────┤
│ Extension windows (Smith-Waterman local on ±30 bp):              │
│   pos 42:  score 94, CIGAR 20=1X12=2D15=  ✓ best                 │
│   pos 217: score 52, no extension beats threshold                 │
├──────────────────────────────────────────────────────────────────┤
│ Final alignment reported: pos 42, MAPQ 60                        │
└──────────────────────────────────────────────────────────────────┘
```

### Student controls

- A small set of preset (reference, read) pairs with named scenarios: clean match, SNV only, indel, repetitive region (two good seeds), chimeric read. Student picks one and steps through.
- Seed length slider. The student sees what happens when k is too small (many ambiguous seeds, slow extension) or too large (no seeds at all in a slightly mismatched region).
- Three-button pipeline: Find seeds → Extend → Report.

### What they see

- **Stage 1 (Find seeds)**: k-mers from the read are queried against the reference's FM index; hits are drawn as vertical markers on a horizontal reference track. The FM interval size for each seed is shown.
- **Stage 2 (Extend)**: for each seed, a small Smith-Waterman matrix is run on a ±30 bp window around the seed position. The student can click a seed to open the local SW matrix (reused from Artifact #5).
- **Stage 3 (Report)**: the winning alignment is reported with CIGAR and MAPQ. Losing extensions are shown greyed out with their scores.

### Target aha moment

"The FM index is fast but doesn't tolerate errors. Smith-Waterman tolerates errors but is slow. Seed-and-extend uses each where it's strong. Every production aligner you've heard of is some variant of this." The student should be able to articulate why a 50 bp read with a single error still aligns — because at least one k-mer in it is error-free, so the seed stage catches it even though a whole-read exact query wouldn't.

### Technical notes

- Reference is 400 bp so the seed plot stays legible. Use a single hardcoded reference excerpt; preset reads are derived from it with controlled edits.
- Clicking a seed to open its SW window reuses the Artifact #5 component (factor it out if possible).
- Seed length slider exposes the core tradeoff — document it with a short in-artifact note.

### Acceptance criteria

- [ ] Default preset finds at least two seeds and produces a valid alignment with a non-trivial CIGAR.
- [ ] Repetitive-region preset produces ≥ 3 candidate seeds and only one wins extension.
- [ ] Seed length k = 8 noticeably clutters the seed plot; k = 18 returns zero seeds on the SNV preset. Both illustrate the tradeoff.
- [ ] The SW window invoked from a seed matches Artifact #5's behavior.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS — extract it into `artifacts/lecture-02/_shared.css` (inlined into each HTML file at build time so the standalone property holds).
- All artifacts share the same DNA base colors from the design tokens. Never substitute.
- All artifacts open with a meaningful default; none is blank on load.
- The Smith-Waterman component from Artifact #5 is reused inside Artifact #6. Implement it as a self-contained JS module so it can be dropped into both.

## 9. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no errors in the console.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720px width; degrades gracefully at 1200px.
- [ ] No reliance on color alone for meaning (§9 of the diagram style guide applies here too).
- [ ] No `alert()`, no console spam, no external calls.
