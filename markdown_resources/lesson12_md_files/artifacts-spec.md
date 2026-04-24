# Lecture 12 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 12 (Population Genetics Fundamentals).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-12.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-12/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Allele-frequency plots: allele A cobalt `#1e3a8a`, allele a crimson `#c4342c`.
- Ancestry / ADMIXTURE components: categorical palette cobalt / crimson / amber / teal / violet / slate.
- Coalescent trees: lineages charcoal `#0a0a0a`; coalescence events highlighted with amber `#b45309`.
- Selection statistics: under-neutrality grey `#525252`; extreme tail (top 1%) amber `#b45309`.
- Typography: Inter for UI chrome; JetBrains Mono for numerical values, statistics, genotype strings, parameter labels.
- Default state is instructive: opens with a meaningful example pre-rendered.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters with sensible ranges.
- **Step / Run / Reset** — for iterative simulations (Wright-Fisher generations, coalescent merges).
- **Re-simulate** — for stochastic simulations (draws a new random realisation from the same parameters).
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Hardy-Weinberg Explorer → chi-squared HWE test result for the current observed-vs-expected counts; verdict "consistent with HWE" vs "deviates significantly".
- Wright-Fisher Drift + Selection Simulator → fixation probability and median time-to-fixation over the simulated ensemble; regime label (drift-dominated / moderate / selection-dominated) based on $2N_e s$.
- LD Decay Explorer → fitted decay constant $2ct$; half-max recombination distance.
- Coalescent Tree Simulator → observed TMRCA, total tree length, and comparison to analytical expectation $4N_e(1 - 1/n)$.
- PSMC Demographic Reconstructor → recovered $N_e(t)$ curve overlaid on ground truth; RMSE in log-$N_e$; regions where inference is well-constrained vs under-constrained.
- ADMIXTURE Component Explorer → cross-validation error per K; recovered-vs-true ancestry match quality (correlation / Frobenius distance).
- Selection Scan — iHS and Tajima's D → per-window statistic value; flagged windows at the 1% tail; true/false-positive breakdown relative to planted selection site.

### 1.5 Feasibility gate on user input (required where input is free-form)

- All seven artifacts here use slider / dropdown / numeric-input controls only (no free-form text or uploads).
- Numeric inputs clamp to their valid range and display the clamped value.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted for in-page math; otherwise none).
- No external data files > 100 KB (per-artifact hardcoded demographic / population presets).

---

## 2. Artifact #1 — Hardy-Weinberg Explorer

**File**: `artifacts/lecture-12/01-hardy-weinberg.html`
**Lecture anchor**: §1.2 Hardy-Weinberg equilibrium
**EE framing reinforced**: HWE as stationary distribution of a memoryless Markov chain; chi-squared test as detection of memory effects.

### Teaching purpose

Slide the allele frequency $p$ and watch the expected Hardy-Weinberg genotype frequencies $(p^2, 2pq, q^2)$ update live. Overlay an "observed" dataset that the student can perturb toward excess homozygosity or excess heterozygosity. Compute a chi-squared HWE test; see when and why the test rejects.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Allele frequency p:   [──●── 0.50 ]                         │
│ Cohort size:          [──●── 500 diploids ]                 │
│ Deviation knob:       [──●── 0.0  ]  (− = excess het;       │
│                                      + = excess homo)       │
│ [Re-sample]                                                  │
├─────────────────────────────────────────────────────────────┤
│ Genotype-frequency plot:                                     │
│   Expected (HWE) curves for AA, Aa, aa across p              │
│   Current-p bars overlaid; observed counts as dots           │
├─────────────────────────────────────────────────────────────┤
│ Observed vs expected table:                                  │
│   AA: exp=125 obs=148                                       │
│   Aa: exp=250 obs=204                                       │
│   aa: exp=125 obs=148                                       │
│   χ² = 12.97   df=1   p = 3.2e-4                             │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ p = 0.50; cohort = 500                                 │   │
│ │ Deviation: +0.20 (excess homozygosity)                 │   │
│ │ χ² = 12.97 (p = 3.2e-4) → REJECTS HWE at α = 0.05      │   │
│ │ Interpretation: population structure or inbreeding     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Allele frequency slider $p$: 0.01 – 0.99 (default 0.50).
- Cohort-size slider $N$: 20 – 5000 (default 500).
- Deviation slider: −0.3 to +0.3 (default 0). Negative excess heterozygotes; positive excess homozygotes. Interpreted as inbreeding coefficient $F$.
- Re-sample button (re-draws observed counts from multinomial with deviation-adjusted expected frequencies).

### What they see

- Three shared-axis curves of $p^2, 2pq, q^2$ over $p \in [0,1]$ with current-$p$ vertical line.
- Bars for the three genotype counts: expected (translucent) and observed (solid) side by side.
- Observed-vs-expected table with chi-squared computation.
- Outcome banner with verdict.

### Target aha

With deviation = 0, observed matches expected (χ² small, fail to reject). Crank deviation up (excess homo) → chi-squared grows, p-value drops — student sees that the **deviation**, not the allele frequency itself, is what the test detects. Shrink the cohort to 20 diploids → even with deviation = 0, the test sometimes rejects by chance, and sometimes fails to reject even with strong deviation — student sees the role of sample size in detection power.

### Technical notes

- Pure JS, seeded mulberry32.
- Expected counts: $N \cdot (p^2, 2pq, q^2)$.
- Deviation applied as an inbreeding coefficient $F$: genotype frequencies $(p^2 + Fpq, 2pq(1-F), q^2 + Fpq)$.
- Observed counts drawn from multinomial with those probabilities.
- Chi-squared statistic: $\sum (O - E)^2 / E$ with df = 1 (one allele-frequency parameter estimated).
- KaTeX used for the inline formulas; defer-loaded from CDN.

### Acceptance criteria

- [ ] Default (p=0.5, N=500, F=0) shows χ² small, fail to reject.
- [ ] Setting F=+0.2 → χ² large, reject at α=0.05.
- [ ] Setting N=20 and F=0 → rejection fraction on repeated Re-samples is roughly 5%.
- [ ] Setting p=0.99 → the heterozygote bar is visibly tiny; AA bar dominates.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Wright-Fisher Drift + Selection Simulator

**File**: `artifacts/lecture-12/02-wright-fisher.html`
**Lecture anchor**: §2.3 Selection
**EE framing reinforced**: WF as discrete Markov chain with bias ($s$) + diffusion ($1/N$); the critical scaling is $2Ns$, not $s$ or $N$ separately.

### Teaching purpose

Simulate Wright-Fisher + selection for many replicate populations. Plot allele-frequency trajectories; compute fixation probability and time to fixation; show that $2Ns$ (not $s$ or $N$ alone) governs whether a locus behaves as drift-dominated, moderate, or selection-dominated.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Population size N:        [──●── 1000 ]                     │
│ Selection coefficient s:  [──●── 0.01 ]  (negative allowed) │
│ Dominance h:              [──●── 0.5  ]                     │
│ Initial frequency p0:     [──●── 0.05 ]                     │
│ Generations:              [──●── 2000 ]                     │
│ Replicates:               [──●── 100  ]                     │
│ [Run]  [Reset]                                              │
├─────────────────────────────────────────────────────────────┤
│ Trajectory plot:                                            │
│   100 coloured traces of p vs t; mean trajectory bolded     │
│   0 and 1 absorbing lines shaded                            │
├─────────────────────────────────────────────────────────────┤
│ Histograms:                                                 │
│   left: fixation-time distribution (for traces that fixed)  │
│   right: final-frequency histogram                          │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 2Ns = 20  →  SELECTION-DOMINATED regime                │   │
│ │ P(fix) = 0.78  (vs neutral expectation 0.05)           │   │
│ │ Median fixation time: 420 generations                   │   │
│ │ Try N=100, s=0.01 → 2Ns=2 (moderate); or s=0 (neutral) │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- $N$: 20 – 10000 (log slider), default 1000.
- $s$: −0.05 – +0.1, default +0.01.
- $h$: 0 (recessive) – 1 (dominant), default 0.5 (additive).
- $p_0$: 0.001 – 0.5, default 0.05.
- Generations: 100 – 5000, default 2000.
- Replicates: 10 – 500, default 100.
- Run / Reset buttons.

### What they see

- Trajectory plot with one line per replicate; mean trajectory emphasised.
- Fixation-time histogram (among replicates that fixed).
- Final-frequency histogram (binned over replicates).
- Outcome banner reporting the regime via $2Ns$, fixation probability, and median fixation time.

### Target aha

Start with $N=1000, s=0.01$ → $2Ns=20$, strong sweep, most replicates fix. Reduce to $N=50, s=0.01$ → $2Ns=1$, drift dominates, many loci are lost despite positive selection. Increase to $N=50, s=0.2$ → $2Ns=20$ again, selection wins. Student sees that **the product $2Ns$, not either parameter alone, sets the regime.**

### Technical notes

- Pure JS, seeded PRNG for reproducibility.
- WF + selection update per generation:
  - Compute $\bar{w} = p^2 w_{AA} + 2pq w_{Aa} + q^2 w_{aa}$.
  - Deterministic $p' = p(p w_{AA} + q w_{Aa}) / \bar{w}$.
  - Stochastic: $X_{t+1} \sim \text{Binomial}(2N, p')$; $p_{t+1} = X_{t+1}/(2N)$.
- $w_{AA}=1+s, w_{Aa}=1+hs, w_{aa}=1$ (selection favouring $A$).
- Trajectories truncated at fixation (absorbing).
- Histograms use ≤ 25 bins each.

### Acceptance criteria

- [ ] Default (N=1000, s=0.01, p0=0.05) → P(fix) > 0.5 over 100 replicates; banner shows "selection-dominated".
- [ ] Setting s=0 → P(fix) ≈ p0 (neutral expectation).
- [ ] Setting N=50, s=0.01 → P(fix) near neutral; banner shows "drift-dominated" or "moderate".
- [ ] Negative s with p0=0.5 → trajectories drift down on average.
- [ ] Opens pre-rendered with the default replicate ensemble.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — LD Decay Explorer

**File**: `artifacts/lecture-12/03-ld-decay.html`
**Lecture anchor**: §3.3 LD decay and haplotype blocks
**EE framing reinforced**: LD as autocorrelation; decay governed by $2ct$ (recombination rate × time).

### Teaching purpose

Simulate a short chromosome region where a starting LD pattern (e.g. founder haplotype) decays under recombination across generations. Plot mean $r^2$ vs genomic distance; fit decay curve; extract decay constant. Show how recombination rate and time since founding interact.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Recombination rate c:   [──●── 1.0 cM/Mb ]                  │
│ Generations since founding: [──●── 100 ]                    │
│ Population size N:       [──●── 10000 ]                     │
│ Region length:           [──●── 1 Mb ]                      │
│ SNPs in region:          [──●── 100 ]                       │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ LD decay curve:                                             │
│   x: pairwise distance (log)                                │
│   y: mean r² per distance bin                               │
│   Overlay: fitted exp(-2ct · d) curve                       │
├─────────────────────────────────────────────────────────────┤
│ LD heatmap:                                                 │
│   100×100 pairwise r² matrix, triangular block structure    │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Fitted decay constant: 2ct = 2.0 × 10⁻⁶ / bp           │   │
│ │ Half-max distance: ~350 kb                             │   │
│ │ Increase t to 1000 → r² drops to noise beyond ~35 kb   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Recombination rate $c$: 0.1 – 5 cM/Mb, default 1 cM/Mb.
- Generations $t$: 10 – 10000 (log slider), default 100.
- Population size $N$: 500 – 100000 (log slider), default 10000.
- Region length: 100 kb – 5 Mb, default 1 Mb.
- SNPs in region: 50 – 200, default 100.
- Re-simulate.

### What they see

- LD decay curve with fitted exponential overlay.
- LD heatmap (triangular block-diagonal) for the simulated region.
- Outcome banner with fitted decay constant and half-max distance.

### Target aha

Start default: decay curve with characteristic half-max near 300–500 kb. Crank generations to 10000 → LD now drops to noise within a few kb. Drop generations to 10 and crank $c$ down → LD persists across the whole region. Student sees the decay constant is the product $2ct$ — either axis (time or recombination rate) moves it.

### Technical notes

- Pure JS, seeded PRNG.
- Simplified model: start with a single founder haplotype repeated $N$ times; each generation, each chromosome picks a parent and undergoes recombination at Poisson rate $c$ per unit distance; SNPs are introduced at a fixed mutation rate so that allele frequencies are non-trivial.
- Compute pairwise $r^2$ across all SNP pairs after $t$ generations.
- Bin pairwise distances (log-spaced); plot mean $r^2$ per bin.
- Fit $r^2(d) = \exp(-2ct \cdot d)$ by least squares on log space.

### Acceptance criteria

- [ ] Default (c=1, t=100) shows plausible decay with half-max ~300 kb.
- [ ] Increasing $t$ 10× shrinks half-max ~10×.
- [ ] Heatmap shows block structure consistent with the decay curve.
- [ ] Fitted decay constant changes with $c \cdot t$.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Coalescent Tree Simulator

**File**: `artifacts/lecture-12/04-coalescent-tree.html`
**Lecture anchor**: §4.2 The standard coalescent
**EE framing reinforced**: coalescent as backward-time Markov chain; exponential inter-coalescence waiting times scaled by $\binom{k}{2}/(2N)$.

### Teaching purpose

Simulate Kingman's standard coalescent for user-specified $n$ and $N_e$. Draw the resulting tree with branch lengths proportional to coalescent waiting times. Show that TMRCA ≈ $4N_e$ regardless of $n$, because the **last** coalescence dominates the total time.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Sample size n:    [──●── 10  ]                              │
│ Effective pop N:  [──●── 10000 ]                            │
│ Replicates:       [──●── 200 ]                              │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Tree panel:                                                 │
│   One realised coalescent tree, present at bottom            │
│   Branch lengths ∝ waiting times T_k                        │
│   Each coalescence labelled with its T_k value               │
├─────────────────────────────────────────────────────────────┤
│ Waiting-time decomposition:                                 │
│   T_n:   180 gens        bar: █                             │
│   T_{n-1}: 240 gens      bar: ██                            │
│   …                                                         │
│   T_2: 14500 gens        bar: ████████████████████          │
├─────────────────────────────────────────────────────────────┤
│ Ensemble histogram:                                         │
│   TMRCA distribution over 200 replicates                    │
│   expected 4N(1 − 1/n) = 36000 marked                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Observed mean TMRCA: 35820 generations                 │   │
│ │ Expected 4N(1 − 1/n) = 36000  (match within 1%)        │   │
│ │ Last coalescence (T₂) ~= 40% of total tree length      │   │
│ │ Double n to 20 → TMRCA barely changes                  │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Sample size $n$: 2 – 50, default 10.
- Effective population size $N_e$: 100 – 100000 (log slider), default 10000.
- Replicates: 20 – 1000, default 200.
- Re-simulate.

### What they see

- A single realised coalescent tree drawn as an SVG cladogram (present at bottom, past at top).
- Waiting-time decomposition chart: each $T_k$ as a horizontal bar labelled with its value.
- Ensemble histogram of TMRCA with analytical expectation marked.
- Outcome banner.

### Target aha

Simulate with $n=10$ → TMRCA around $4N_e$. Increase $n$ to 50 → TMRCA barely changes. Zoom into the waiting-time bars: $T_2$ (last coalescence) is enormous; $T_{50}$ (first) is tiny. Student sees: the **last** pair of lineages takes most of the time, and that's why TMRCA is nearly flat in $n$.

### Technical notes

- Pure JS, seeded PRNG.
- Standard Kingman algorithm: at each step with $k$ active lineages, draw $T_k \sim \text{Exp}(\binom{k}{2}/(2N))$; pick a random pair to coalesce.
- Tree layout: recursive assignment of x-positions to leaves; y = coalescent time.
- Ensemble: compute TMRCA distribution over the replicate set.
- Compare observed mean to analytical $4N(1 - 1/n)$.

### Acceptance criteria

- [ ] Default (n=10, N=10000, 200 reps) gives mean TMRCA within 10% of $4N(1 - 1/n)$.
- [ ] Doubling $n$ changes TMRCA by less than ~10% (vs doubling $N$ which doubles it).
- [ ] Waiting-time bars show monotonically increasing lengths (T_n smallest, T_2 largest).
- [ ] Tree redraws on Re-simulate.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — PSMC Demographic Reconstructor

**File**: `artifacts/lecture-12/05-psmc.html`
**Lecture anchor**: §5.1 PSMC — from a single genome
**EE framing reinforced**: PSMC as HMM whose hidden state is coalescent time; inference of $N_e(t)$ as state-occupancy inversion.

### Teaching purpose

Simulate a single diploid genome under a user-defined demographic history $N_e(t)$ (step function). Run a simplified PSMC-style HMM inference and recover $N_e(t)$. Compare recovered to ground truth. Show PSMC's strength (captures ancient bottlenecks from **one** diploid) and weakness (poor resolution in very recent and very ancient past).

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Demographic history (N_e vs time):                          │
│   Step editor with draggable breakpoints                    │
│   Presets: [ out-of-Africa ▾ / bottleneck / constant /      │
│             recent expansion / custom ]                      │
├─────────────────────────────────────────────────────────────┤
│ Simulation parameters:                                      │
│   Mutation rate μ:   [──●── 1.2e-8 /bp/gen ]                │
│   Recombination r:   [──●── 1.0e-8 /bp/gen ]                │
│   Genome size L:     [──●── 100 Mb ]                        │
│   Window size:       [──●── 100 kb ]                        │
│   Time discretisation: [──●── 40 bins ]                     │
│ [Simulate genome]  [Run PSMC-like inference]                │
├─────────────────────────────────────────────────────────────┤
│ N_e(t) plot:                                                │
│   Ground truth as solid step function                        │
│   Recovered as overlaid broken-line with uncertainty band   │
│   Log-log axes                                              │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Recovered N_e(t) matches truth within 20% for t ∈      │
│ │   [30 kya, 500 kya]; biased in [0, 20 kya] and beyond  │
│ │   1 Mya.                                               │   │
│ │ Log-Ne RMSE over reliable range: 0.18                  │   │
│ │ Bottleneck at 100 kya correctly recovered              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Demographic-history preset dropdown + step-function editor (drag breakpoint to change $N_e$ or time).
- Mutation rate $\mu$: 0.5e-8 – 2e-8 per bp per generation, default 1.2e-8.
- Recombination rate $r$: 0.5e-8 – 2e-8, default 1e-8.
- Genome size $L$: 10 Mb – 1 Gb (log), default 100 Mb.
- Window size: 10 kb – 1 Mb (log), default 100 kb.
- Time discretisation bins: 10 – 64, default 40.
- Simulate / Run buttons.

### What they see

- Step editor for the demographic history (ground truth).
- $N_e(t)$ plot with ground truth as solid line and recovered curve overlaid with an uncertainty band.
- Outcome banner with RMSE and reliable-range annotation.

### Target aha

Start with the out-of-Africa preset: simulate + run → recovered curve tracks the bottleneck around 100 kya. Change the bottleneck location to 10 kya → **recovered curve misses it** (too recent for PSMC). Move it to 5 Mya → also missed (too ancient — no coalescences old enough). Student sees: PSMC has a reliable time window determined by the coalescent-time distribution of a single diploid genome.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate a single diploid genome: along a 100 Mb sequence, draw a piecewise-constant coalescent-time trajectory (approximating the SMC') using the user's $N_e(t)$; emit per-window heterozygosity as Poisson($2\mu \cdot$ coalescent time $\cdot$ window-size).
- Inference: discretise time into bins; fit a Gaussian-emission HMM where the hidden state is the coalescent time bin; transition matrix approximates the SMC' correlations between adjacent windows.
- Baum-Welch for a few EM iterations to recover per-bin emission rates → invert to $N_e(t)$.
- Uncertainty band: parametric-bootstrap estimate from final log-likelihood curvature.

### Acceptance criteria

- [ ] Default preset (out-of-Africa) recovers the ancient bottleneck with log-$N_e$ RMSE ≤ 0.3 in the $[30$ kya, $500$ kya$]$ range.
- [ ] Very recent bottleneck (< 20 kya) systematically biased.
- [ ] Very ancient bottleneck (> 2 Mya) systematically biased or unrecovered.
- [ ] Step editor produces a valid monotonic-time piecewise-constant $N_e(t)$.
- [ ] Opens pre-rendered with a default simulation + inference already computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — ADMIXTURE Component Explorer

**File**: `artifacts/lecture-12/06-admixture.html`
**Lecture anchor**: §5.3 Admixture inference
**EE framing reinforced**: admixture inference as mixture-model decomposition; K as model-order selection (not biological truth).

### Teaching purpose

Simulate a multi-population cohort where the true number of ancestral populations is user-controlled ($K_{\text{true}} \in \{2, 3, 4\}$), with some admixed individuals. Run a simplified ADMIXTURE-style EM for user-specified $K_{\text{fit}}$. Show that fitting too large a $K$ just splits a true cluster in half — the likelihood keeps improving but the biology doesn't.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ True number of populations K_true:  [──●── 3  ]             │
│ Individuals per population:        [──●── 40 ]              │
│ Admixed fraction:                   [──●── 10% ]            │
│ SNPs:                               [──●── 500 ]            │
│ Fitted K:                           [──●── 3  ]             │
│ [Re-simulate]  [Run ADMIXTURE (EM)]                         │
├─────────────────────────────────────────────────────────────┤
│ Stacked-bar ancestry plot:                                  │
│   one column per individual, K stacked component colours    │
│   individuals sorted by majority component                  │
├─────────────────────────────────────────────────────────────┤
│ Recovered vs true mapping:                                  │
│   Hungarian-matching to align recovered ↔ true components   │
│   Match table with per-component correlation                │
├─────────────────────────────────────────────────────────────┤
│ K-sweep panel:                                              │
│   Plot CV error (from 5-fold CV) vs K_fit ∈ {2..7}          │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ K_true=3, K_fit=3: 95% component-correlation match     │   │
│ │ Over-fit K_fit=5: one true cluster split into two      │   │
│ │   (recovered components 3 and 4 both inherit it)       │   │
│ │ CV error minimum at K=3 (matches truth)                │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- $K_{\text{true}}$: 2 – 4, default 3.
- Individuals per population: 10 – 100, default 40.
- Admixed fraction: 0 – 50%, default 10%.
- SNPs: 100 – 2000 (log), default 500.
- $K_{\text{fit}}$: 2 – 7, default 3.
- Re-simulate / Run buttons.

### What they see

- Stacked-bar ancestry plot (per-individual mixture).
- Recovered-vs-true component mapping table (Hungarian-matched).
- CV error sweep across $K_{\text{fit}}$.
- Outcome banner.

### Target aha

Simulate $K_{\text{true}}=3$, fit $K=3$ → clean recovery. Fit $K=5$ → a true cluster gets split into two recovered components. CV error vs $K$: clear minimum at $K=3$, noisy afterwards. Student sees: **K is a hyperparameter**, and choosing it is a model-selection problem — no biological "true K" is being recovered.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate allele frequencies per population (draw from Dirichlet-like distribution for differentiation); generate genotypes; admixed individuals have mixture-of-population ancestry.
- EM for ADMIXTURE-style objective: Q (per-individual ancestry) × P (per-component allele frequencies).
- Hungarian matching (simple approximation) to align recovered components with truth.
- 5-fold CV: mask 20% of genotypes, refit, predict masked.

### Acceptance criteria

- [ ] Default (K_true=3, K_fit=3) recovers matching components with > 0.9 correlation.
- [ ] Over-fitting K_fit=5 → at least one true cluster is split (evident in bar plot).
- [ ] CV minimum near $K_{\text{true}}$.
- [ ] Re-simulate yields fresh data.
- [ ] Opens pre-rendered with default fit.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Selection Scan — iHS and Tajima's D

**File**: `artifacts/lecture-12/07-selection-scan.html`
**Lecture anchor**: §6.2 iHS and EHH; §6.1 Tajima's D
**EE framing reinforced**: selection scans as multi-hypothesis detection with controlled FAR; sensitivity tied to null-distribution estimation.

### Teaching purpose

Simulate a coalescent region under neutrality (null) and a separate region under recent positive selection (signal). Compute Tajima's D and a simplified iHS in sliding windows. Plot both as Manhattan-style scans; mark the planted sweep; report precision/recall at user-controlled significance thresholds.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulation:                                                 │
│   Region length:   [──●── 2 Mb ]                            │
│   N_e:             [──●── 10000 ]                           │
│   Sample size:     [──●── 50 ]                              │
│   Selection s:     [──●── 0.05 ] (set to 0 for pure null)   │
│   Sweep frequency: [──●── 0.7 ] (current allele freq)       │
│   Sweep position:  [──●── 1 Mb ] (centre of the region)     │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Scan plot (Manhattan):                                      │
│   Top:    Tajima's D vs window position                     │
│   Bottom: iHS vs window position                            │
│   Threshold line (user-adjustable percentile)               │
│   Planted sweep site marked with vertical bar               │
├─────────────────────────────────────────────────────────────┤
│ Detection summary:                                          │
│   Flagged windows: 12                                       │
│   True positives:  1  (the planted sweep)                   │
│   False positives: 11                                       │
│   Precision = 8%,  Recall = 100% (for the one sweep)        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Planted sweep detected at genome-wide 0.01 threshold   │   │
│ │ 11 false positives (expected for 200 windows × 1%)     │   │
│ │ Drop threshold to 0.001 → 2 FPs; still flags sweep     │   │
│ │ Set s=0 → no sweep to detect; FPs still happen          │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Region length: 500 kb – 5 Mb, default 2 Mb.
- $N_e$: 1000 – 50000 (log), default 10000.
- Sample size: 20 – 200, default 50.
- Selection coefficient $s$: 0 – 0.1, default 0.05.
- Sweep current frequency: 0.1 – 0.95, default 0.7.
- Sweep position within the region: default centre.
- Significance threshold (percentile): default 99th.
- Re-simulate.

### What they see

- Two-panel Manhattan plot (Tajima's D top, iHS bottom) with threshold line and the planted sweep marked.
- Detection summary box with TP/FP counts and precision/recall.
- Outcome banner.

### Target aha

With $s=0.05$, planted sweep → both statistics peak near the sweep site, but so do several unrelated windows (false positives). Drop $s$ to 0 → **the sweep is gone, but the FPs are still there** (by construction, under null). Student sees: a single-statistic significant window isn't a claim about selection — it's a candidate that must be validated against properly-calibrated neutral null, ideally with multiple statistics.

### Technical notes

- Pure JS, seeded PRNG.
- Approximate coalescent simulation for the neutral background (ms-style, simplified to 1–2 recombination breakpoints per Mb per generation).
- Selection simulated as a stochastic sweep: place the beneficial allele, run WF forward with selection to reach the target frequency, embed in the surrounding neutral tree.
- Tajima's D per window: $\pi$, $S$, Tajima formula from §6.1.
- Simplified iHS: compute extended-haplotype-homozygosity decay from the sweep allele; standardise.
- Windowing: non-overlapping ~10 kb windows.
- Threshold line at user-selected percentile.

### Acceptance criteria

- [ ] Default (s=0.05) → planted sweep is flagged by at least one of the two statistics.
- [ ] Pure neutral (s=0) → roughly 1% of windows exceed the 99th-percentile threshold (sanity).
- [ ] Tightening the threshold reduces FPs as expected.
- [ ] Re-simulate yields fresh scan.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Allele / ancestry / selection colours are consistent across artifacts (palette defined in §1.2 above).
- Coalescent / WF / LD artifacts use a common "replicates ensemble" pattern (draw many, summarise).
- Every artifact emits an **outcome banner** per convention §1.4.
- Parameter scales (log sliders for $N_e$, mutation rate, recombination rate) are consistent across artifacts.

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
