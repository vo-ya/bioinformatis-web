"""Build exercise.ipynb for L12 — Population Genetics Fundamentals.

Run from this directory:
    python3 build_notebook.py

Emits exercise.ipynb. Re-running overwrites the file.
"""

from __future__ import annotations

import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def md(text: str):
    """Markdown cell."""
    return new_markdown_cell(text)


def code(source: str):
    """Plain (visible) code cell."""
    return new_code_cell(source)


def hidden(source: str):
    """Code cell whose source is collapsed by default in Colab / Jupyter."""
    cell = new_code_cell(source)
    cell.metadata = {
        "jupyter": {"source_hidden": True},
        "cellView": "form",
    }
    return cell


# ---------------------------------------------------------------------------
# Cell sources
# ---------------------------------------------------------------------------


TITLE_MD = """# L12 — Population Genetics Fundamentals

In this exercise you simulate the **Wright-Fisher** Markov chain, the workhorse
forward-time model of allele frequencies in a finite population. You will run
100 replicate trajectories of neutral drift, recover the effective population
size `Nₑ` from the observed variance, then layer a weak selection coefficient
on top and watch the deterministic trend emerge from stochastic noise.
"""


AHA_MD = """> **Aha.** Drift is a **1-D random walk** whose variance grows like
> `p(1-p)·t/Nₑ`. The √-time spreading is the same square-root dispersion that
> shows up in Brownian motion, shot noise, and Wiener processes. Weak
> selection adds a **deterministic drift term** to that stochastic noise — the
> signal-to-noise ratio scales as **Nₑ · s**: bigger populations or stronger
> selection give the trend room to dominate the noise.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Deterministic for the whole notebook.
MASTER_SEED = 42
rng = np.random.default_rng(MASTER_SEED)

print(f"NumPy {np.__version__}, master seed = {MASTER_SEED}")
"""


STEP1_MD = """## Step 1 (8 min) — One Wright-Fisher generation: binomial resampling

The Wright-Fisher model is the simplest forward-time discrete-generation model
in population genetics. With `N` diploid individuals there are `2N` allele
copies in the population. Each generation the next-generation pool is built
by sampling `2N` chromosomes **with replacement** from the current pool;
equivalently, the next-generation count of allele `A` is

$$X_{t+1} \\sim \\mathrm{Binomial}(2N,\\; p_t)$$

where `p_t = X_t / 2N` is the current allele frequency. That is the whole
model: a Markov chain on `{0, 1, …, 2N}` with absorbing barriers at 0 and 2N
(fixation / loss). Implement one generation and verify it gives the right
mean and variance over many resamples.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — One Wright-Fisher generation.
# ----------------------------------------------------------------------

N_DIPLOID = 1000     # population size
TWO_N     = 2 * N_DIPLOID
P0        = 0.5      # starting allele frequency


def wf_step(p: float, two_n: int, rng_local) -> float:
    """Sample one generation of Wright-Fisher drift; return next frequency."""
    # TODO: draw a Binomial(two_n, p) and convert back to a frequency.
    raise NotImplementedError


# Sanity check: simulate 5000 single-generation steps from p0=0.5 and verify
# E[Delta p] ~= 0 and Var[Delta p] ~= p(1-p) / (2N).
# TODO: run the check and print the empirical vs theoretical numbers.
'''

STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

N_DIPLOID = 1000
TWO_N     = 2 * N_DIPLOID
P0        = 0.5


def wf_step(p: float, two_n: int, rng_local) -> float:
    x_next = rng_local.binomial(two_n, p)
    return x_next / two_n


# Sanity check.
rng_check = np.random.default_rng(7)
n_trials = 5000
deltas = np.empty(n_trials)
for i in range(n_trials):
    p_next = wf_step(P0, TWO_N, rng_check)
    deltas[i] = p_next - P0

mean_emp = deltas.mean()
var_emp  = deltas.var(ddof=1)
var_theory = P0 * (1 - P0) / TWO_N  # Binomial variance / (2N)^2

print(f"E[Delta p]   empirical = {mean_emp:+.6f}   theoretical = 0.000000")
print(f"Var[Delta p] empirical = {var_emp:.6e}   theoretical = {var_theory:.6e}")
print(f"  ratio empirical/theoretical = {var_emp/var_theory:.3f}  (should be ~1)")
'''


STEP2_MD = """## Step 2 (15 min) — 100 replicates of neutral drift; plot 5 trajectories

The single-generation kernel above is now iterated for `T = 200` generations
across `R = 100` independent replicate populations, all starting from the
same `p_0 = 0.5`. Each replicate is one trajectory of the Markov chain.
Plot 5 of them on the same axes and overlay the theoretical drift envelope
`p_0 ± √(p_0(1-p_0)·t/(2N))`: this is the **one-σ band** of a Brownian-motion
approximation to drift, accurate while no replicate gets too close to the
fixation / loss boundaries.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Replicate Wright-Fisher trajectories.
# ----------------------------------------------------------------------

T_GENS    = 200      # generations per replicate
R_REPS    = 100      # number of replicate populations
P0        = 0.5


def simulate_drift(p0: float, two_n: int, T: int, R: int,
                   seed: int = 0) -> np.ndarray:
    """Return an (R, T+1) array of allele-frequency trajectories.

    Row r, column t = allele frequency in replicate r at generation t.
    Column 0 = p0 for every replicate.
    """
    # TODO: allocate an (R, T+1) array; fill column 0 with p0; iterate.
    raise NotImplementedError


# trajectories = simulate_drift(P0, TWO_N, T_GENS, R_REPS, seed=1)
# Plot 5 trajectories + the +/- 1 sigma drift envelope.
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

T_GENS = 200
R_REPS = 100
P0     = 0.5


def simulate_drift(p0: float, two_n: int, T: int, R: int,
                   seed: int = 0) -> np.ndarray:
    rng_local = np.random.default_rng(seed)
    traj = np.empty((R, T + 1), dtype=np.float64)
    traj[:, 0] = p0
    # Vectorised: at each generation, draw R independent Binomial(2N, p_r).
    for t in range(T):
        traj[:, t + 1] = rng_local.binomial(two_n, traj[:, t]) / two_n
    return traj


trajectories = simulate_drift(P0, TWO_N, T_GENS, R_REPS, seed=1)
print(f"Trajectories shape: {trajectories.shape}")

# At t = T_GENS, count fixations (p=1), losses (p=0), and segregating reps.
final = trajectories[:, -1]
n_fixed = int(np.sum(final >= 1.0 - 1e-12))
n_lost  = int(np.sum(final <= 1e-12))
n_seg   = R_REPS - n_fixed - n_lost
print(f"After {T_GENS} generations: {n_fixed} fixed, {n_lost} lost, {n_seg} still segregating")

# Plot 5 trajectories + +/- 1 sigma drift envelope.
ts = np.arange(T_GENS + 1)
sigma_t = np.sqrt(P0 * (1 - P0) * ts / TWO_N)
upper   = np.clip(P0 + sigma_t, 0, 1)
lower   = np.clip(P0 - sigma_t, 0, 1)

fig, ax = plt.subplots(figsize=(9, 4.5))
for r in range(5):
    ax.plot(ts, trajectories[r], lw=1.3, alpha=0.8, label=f"replicate {r}")
ax.fill_between(ts, lower, upper, color="grey", alpha=0.2,
                label=r"$p_0 \\pm \\sqrt{p_0(1-p_0)\\,t/(2N)}$")
ax.axhline(P0, color="black", lw=0.5, alpha=0.6)
ax.set_xlabel("generation")
ax.set_ylabel("allele frequency p")
ax.set_title(f"Wright-Fisher neutral drift  (N={N_DIPLOID}, T={T_GENS}, R={R_REPS})")
ax.set_ylim(0, 1)
ax.legend(loc="upper right", fontsize=9)
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (12 min) — Recover Nₑ from the variance of Δp

The neutral Wright-Fisher kernel has a closed-form per-generation variance:

$$\\mathrm{Var}[\\,p_{t+1} \\mid p_t\\,] = \\frac{p_t(1-p_t)}{2N_e}$$

If you observe many replicates evolving in parallel, you can **invert** that
relationship to recover `Nₑ` from the empirical variance in `Δp = p_{t+1} − p_t`:

$$\\hat{N_e} = \\frac{p_t(1-p_t)}{2\\,\\widehat{\\mathrm{Var}}[\\Delta p]}$$

Compute the estimate at every generation `t` (using the across-replicate
variance), then scatter the estimated `Nₑ` against the true `N = 1000`.
The estimator is noisy at any single `t`, but the time-averaged estimate
should land in a tight band around the truth.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Variance-based effective-population-size estimator.
# ----------------------------------------------------------------------

def estimate_ne(trajectories: np.ndarray) -> np.ndarray:
    """For each generation t in [0, T-1], estimate Nₑ from across-replicate
    variance of (p_{t+1} - p_t). Returns a length-T array of Nₑ estimates."""
    # TODO:
    # 1. Compute deltas = traj[:, 1:] - traj[:, :-1]; shape (R, T).
    # 2. p_t = traj[:, :-1]; mean across replicates -> p_bar_t (length T).
    # 3. var_t = deltas.var(axis=0, ddof=1).
    # 4. ne_t = p_bar_t * (1 - p_bar_t) / (2 * var_t), masked where var_t==0.
    raise NotImplementedError


# ne_estimates = estimate_ne(trajectories)
# Scatter (generation, estimate) and overlay the true N as a horizontal line.
# Also print mean(ne_estimates) over t where it's well-defined.
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

def estimate_ne(trajectories: np.ndarray) -> np.ndarray:
    deltas = np.diff(trajectories, axis=1)               # (R, T)
    p_t    = trajectories[:, :-1]                        # (R, T)
    p_bar  = p_t.mean(axis=0)                            # (T,)
    var_t  = deltas.var(axis=0, ddof=1)                  # (T,)
    # Avoid div-by-zero where every replicate is fixed/lost (var_t -> 0 and
    # p_bar(1-p_bar) -> 0 too; ratio is undefined).
    with np.errstate(divide="ignore", invalid="ignore"):
        ne_t = np.where(var_t > 0,
                        p_bar * (1 - p_bar) / (2.0 * var_t),
                        np.nan)
    return ne_t


ne_estimates = estimate_ne(trajectories)
mask = np.isfinite(ne_estimates) & (ne_estimates > 0)
ne_clean = ne_estimates[mask]

# A robust point estimate: median over generations where the variance is well-defined.
ne_median = float(np.median(ne_clean))
ne_mean   = float(np.mean(ne_clean))

print(f"True N           = {N_DIPLOID}")
print(f"Estimated Nₑ     median  = {ne_median:7.1f}")
print(f"Estimated Nₑ     mean    = {ne_mean:7.1f}")
print(f"  (R = {R_REPS} replicates, T = {T_GENS} generations; "
      f"{mask.sum()} of {len(ne_estimates)} generations had finite var)")

# Scatter true vs estimated.
fig, ax = plt.subplots(figsize=(9, 4))
ts = np.arange(len(ne_estimates))[mask]
ax.scatter(ts, ne_clean, s=10, alpha=0.4, label=r"per-generation $\\hat{N_e}$")
ax.axhline(N_DIPLOID, color="red", lw=2, label=f"true N = {N_DIPLOID}")
ax.axhline(ne_median, color="black", lw=1.5, ls="--",
           label=fr"median $\\hat{{N_e}}$ = {ne_median:.0f}")
ax.set_xlabel("generation t")
ax.set_ylabel(r"$\\hat{N_e}$  =  $p_t(1-p_t) / (2\\,\\widehat{Var}[\\Delta p])$")
ax.set_title("Variance-based effective-population-size estimator")
ax.set_yscale("log")
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Add a selection coefficient s = 0.01

Now make allele `A` slightly fitter. Under viability selection with relative
fitnesses `w_{AA} = 1 + 2s`, `w_{Aa} = 1 + s`, `w_{aa} = 1` (genic / additive
selection, the standard textbook setup), the expected frequency after
selection in an infinite population is

$$p' = \\frac{p\\,(1 + s + sp)}{1 + 2sp}$$

That is the standard one-generation deterministic update for additive
viability selection: numerator is `p · w_A` (allele `A`'s marginal fitness),
denominator is the mean fitness `w̄ = 1 + 2sp`. In a finite population we
apply this deterministic update first, then binomial-sample as in Step 1.
Overlay 5 drift-only trajectories and 5 drift+selection trajectories. For
`Nₑ · s = 10` (here `N = 1000, s = 0.01`) selection should visibly bias the
mean trajectory upward over a few hundred generations; for `Nₑ · s ≪ 1` it
would be invisible behind drift noise.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Wright-Fisher with additive viability selection.
# ----------------------------------------------------------------------

SELECTION_S = 0.01


def wf_step_selection(p: float, two_n: int, s: float, rng_local) -> float:
    """One generation of WF with additive viability selection."""
    # TODO:
    # 1. Deterministic selection update: p_sel = p(1 + s + s*p) / (1 + 2*s*p).
    # 2. Stochastic binomial sampling at p_sel.
    raise NotImplementedError


def simulate_drift_selection(p0: float, two_n: int, s: float,
                             T: int, R: int, seed: int = 0) -> np.ndarray:
    """Return (R, T+1) trajectories under drift + additive selection."""
    # TODO: same loop as Step 2, but call wf_step_selection.
    raise NotImplementedError


# Compare neutral vs s=0.01.
# traj_neutral = simulate_drift(P0, TWO_N, T_GENS, R_REPS, seed=3)
# traj_sel     = simulate_drift_selection(P0, TWO_N, SELECTION_S, T_GENS, R_REPS, seed=3)
# Plot 5 of each on the same axes.
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

SELECTION_S = 0.01


def wf_step_selection(p: float, two_n: int, s: float, rng_local) -> float:
    # Deterministic additive viability selection.
    # Fitnesses w_aa = 1, w_Aa = 1 + s, w_AA = 1 + 2s;
    # marginal fitness of A: w_A = p(1+2s) + (1-p)(1+s) = 1 + s + s*p;
    # mean fitness: w_bar = 1 + 2*s*p.
    p_sel = p * (1.0 + s + s * p) / (1.0 + 2.0 * s * p)
    p_sel = float(np.clip(p_sel, 0.0, 1.0))
    return rng_local.binomial(two_n, p_sel) / two_n


def simulate_drift_selection(p0: float, two_n: int, s: float,
                             T: int, R: int, seed: int = 0) -> np.ndarray:
    rng_local = np.random.default_rng(seed)
    traj = np.empty((R, T + 1), dtype=np.float64)
    traj[:, 0] = p0
    for t in range(T):
        # Vectorised: deterministic selection update across replicates, then
        # an R-vector binomial draw.
        p_t = traj[:, t]
        p_sel = p_t * (1.0 + s + s * p_t) / (1.0 + 2.0 * s * p_t)
        p_sel = np.clip(p_sel, 0.0, 1.0)
        traj[:, t + 1] = rng_local.binomial(two_n, p_sel) / two_n
    return traj


traj_neutral = simulate_drift(P0, TWO_N, T_GENS, R_REPS, seed=3)
traj_sel     = simulate_drift_selection(P0, TWO_N, SELECTION_S, T_GENS, R_REPS, seed=3)

# Mean trajectory across replicates is the cleanest summary.
mean_neutral = traj_neutral.mean(axis=0)
mean_sel     = traj_sel.mean(axis=0)
ts = np.arange(T_GENS + 1)

# Deterministic infinite-population trajectory under additive selection
# (closed-form solution of the recursion ignoring drift).
det = np.empty_like(ts, dtype=float)
det[0] = P0
for t in range(T_GENS):
    pt = det[t]
    det[t + 1] = pt * (1.0 + SELECTION_S + SELECTION_S * pt) / (1.0 + 2.0 * SELECTION_S * pt)

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)

# Left: 5 neutral vs 5 with selection, individual replicates.
ax = axes[0]
for r in range(5):
    ax.plot(ts, traj_neutral[r], color="steelblue", alpha=0.5, lw=1)
    ax.plot(ts, traj_sel[r],     color="crimson",   alpha=0.5, lw=1)
# Dummy lines just for legend entries.
ax.plot([], [], color="steelblue", lw=1.5, label="neutral (s=0)")
ax.plot([], [], color="crimson",   lw=1.5, label=f"selection (s={SELECTION_S})")
ax.axhline(P0, color="black", lw=0.5, alpha=0.6)
ax.set_xlabel("generation")
ax.set_ylabel("allele frequency p")
ax.set_title("Five replicates each: drift only vs drift + selection")
ax.set_ylim(0, 1)
ax.legend(loc="upper right", fontsize=9)

# Right: mean across all 100 replicates, with the deterministic curve.
ax = axes[1]
ax.plot(ts, mean_neutral, color="steelblue", lw=2, label="mean p, neutral")
ax.plot(ts, mean_sel,     color="crimson",   lw=2, label=f"mean p, s={SELECTION_S}")
ax.plot(ts, det, "k--", lw=1.5, label="deterministic limit (no drift)")
ax.axhline(P0, color="black", lw=0.5, alpha=0.6)
ax.set_xlabel("generation")
ax.set_title(f"Mean of R = {R_REPS} replicates  (Nₑ·s = {N_DIPLOID*SELECTION_S:.0f})")
ax.set_ylim(0, 1)
ax.legend(loc="lower right", fontsize=9)

plt.tight_layout()
plt.show()

# Print fixation counts: how often did A reach fixation under each regime?
fix_neutral = int(np.sum(traj_neutral[:, -1] >= 1.0 - 1e-12))
fix_sel     = int(np.sum(traj_sel[:, -1]     >= 1.0 - 1e-12))
print(f"Fixations of A after {T_GENS} generations:")
print(f"  neutral          : {fix_neutral} / {R_REPS}")
print(f"  s = {SELECTION_S} (Nₑs = {N_DIPLOID*SELECTION_S:.0f}) : {fix_sel} / {R_REPS}")
'''


STEP5_MD = """## Step 5 (10 min) — Summary table: per-generation mean and variance

Finally, assemble the canonical population-genetics summary table for the
neutral run: `(generation, mean p, var p)` across replicates, plus a
mini-table comparing neutral vs selection at a handful of generations. The
variance row visualises the **√t drift law**: the variance grows linearly in
`t` (so `σ(p)` grows like √t) while the mean stays put.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Per-generation neutral summary + neutral-vs-selection comparison.
# ----------------------------------------------------------------------

# TODO:
# 1. Build a DataFrame of (generation, mean_p, var_p) for the neutral run, at
#    a sparse set of generations e.g. [0, 25, 50, 100, 150, 200].
# 2. Add columns mean_p_selection and var_p_selection from the selection run.
# 3. Print with reasonable formatting.
# 4. Plot variance vs generation (neutral): empirical points + theoretical
#    line  Var[p_t] = p_0(1-p_0) * (1 - (1 - 1/(2N))^t).
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.

sample_t = [0, 25, 50, 100, 150, 200]

rows = []
for t in sample_t:
    rows.append({
        "generation":         t,
        "mean_p_neutral":     float(traj_neutral[:, t].mean()),
        "var_p_neutral":      float(traj_neutral[:, t].var(ddof=1)),
        "mean_p_selection":   float(traj_sel[:, t].mean()),
        "var_p_selection":    float(traj_sel[:, t].var(ddof=1)),
    })
summary = pd.DataFrame(rows)

# Theoretical variance under neutral drift.
ts_all   = np.arange(T_GENS + 1)
var_theory = P0 * (1 - P0) * (1 - (1 - 1.0 / TWO_N) ** ts_all)

# Pretty print.
print("Per-generation summary  (R = 100 replicates, N = 1000)")
print(summary.to_string(index=False, float_format=lambda v: f"{v:8.4f}"))

# Variance trajectory plot.
var_neutral_all = traj_neutral.var(axis=0, ddof=1)
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(ts_all, var_neutral_all, "o-", ms=2, lw=1.2, alpha=0.7,
        label="empirical Var[p_t]  (R = 100)")
ax.plot(ts_all, var_theory, "r--", lw=2,
        label=r"theory: $p_0(1-p_0)\\,[1 - (1 - 1/(2N))^t]$")
ax.set_xlabel("generation t")
ax.set_ylabel("Var[p_t]")
ax.set_title("Variance of allele frequency under neutral drift")
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout()
plt.show()

# Quick numerical agreement check.
ratio = var_neutral_all[1:] / np.maximum(var_theory[1:], 1e-12)
print(f"Empirical / theoretical variance ratio: "
      f"median = {np.median(ratio):.3f}, min = {ratio.min():.3f}, max = {ratio.max():.3f}")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline. If
you ran the reference solutions above they should all pass; if you wrote
your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. wf_step is unbiased: mean of (p_next - p) ~= 0 over many resamples.
rng_check = np.random.default_rng(101)
deltas_check = np.array([wf_step(0.4, TWO_N, rng_check) - 0.4 for _ in range(2000)])
assert abs(deltas_check.mean()) < 5e-3, (
    f"wf_step appears biased: mean Delta p = {deltas_check.mean():+.4f}"
)

# 2. wf_step variance matches the Binomial-derived theoretical value.
var_emp_check    = deltas_check.var(ddof=1)
var_theory_check = 0.4 * 0.6 / TWO_N
assert 0.5 * var_theory_check < var_emp_check < 2.0 * var_theory_check, (
    f"wf_step variance off: empirical {var_emp_check:.2e} vs theory {var_theory_check:.2e}"
)

# 3. simulate_drift returns the right shape and starts at p0.
traj_check = simulate_drift(0.5, TWO_N, 50, 30, seed=11)
assert traj_check.shape == (30, 51), f"shape mismatch: {traj_check.shape}"
assert np.allclose(traj_check[:, 0], 0.5), "first column should equal p0"

# 4. Effective-population-size estimator recovers N within +/-30 %.
traj_for_ne = simulate_drift(0.5, TWO_N, 200, 200, seed=12)
ne_check    = estimate_ne(traj_for_ne)
ne_check_clean = ne_check[np.isfinite(ne_check) & (ne_check > 0)]
ne_med = float(np.median(ne_check_clean))
assert 0.7 * N_DIPLOID <= ne_med <= 1.3 * N_DIPLOID, (
    f"Nₑ estimate off: median {ne_med:.0f} vs true {N_DIPLOID}"
)
print(f"Nₑ recovery: median estimate {ne_med:.0f} vs true {N_DIPLOID}  (within ±30 %)")

# 5. Selection visibly shifts the mean trajectory at Nₑ·s = 10.
# Independent seeds for the two runs so they are not correlated.
traj_n = simulate_drift          (0.5, TWO_N, T_GENS, 200, seed=13)
traj_s = simulate_drift_selection(0.5, TWO_N, SELECTION_S, T_GENS, 200, seed=29)
mean_n = traj_n[:, -1].mean()
mean_s = traj_s[:, -1].mean()
assert mean_s > mean_n + 0.10, (
    f"selection did not shift mean as expected: "
    f"selection mean p_T = {mean_s:.3f}, neutral mean p_T = {mean_n:.3f}"
)
print(f"Selection lifts mean p by "
      f"{mean_s - mean_n:+.3f} at generation {T_GENS} "
      f"(Nₑs = {N_DIPLOID*SELECTION_S:.0f}).")

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — random walks, drift envelopes, SNR ∝ Nₑ · s

You implemented the canonical population-genetics forward-time model and saw
three EE-flavoured truths:

1. **Drift is a 1-D random walk.** The per-generation increment `Δp` has
   mean 0 and variance `p(1-p)/(2N)`. After `t` generations the variance has
   accumulated to `p_0(1-p_0)·[1 - (1 - 1/(2N))^t] ≈ p_0(1-p_0)·t/(2N)` for
   small `t`. That is **√t spreading**, the discrete-time cousin of Brownian
   motion / a Wiener process. The `±√(p_0(1-p_0)·t/(2N))` envelope in Step 2
   is the one-σ band of that random walk before either absorbing boundary
   makes its presence felt.

2. **Nₑ is an inverse-variance estimator.** Step 3 inverted the kernel:
   `N̂ₑ = p(1-p)/(2·Var[Δp])`. That is exactly a moment-method estimator —
   same idea as inferring the noise variance of an amplifier from the
   spread of its zero-input output. With `R = 100` parallel replicates the
   per-generation estimate is noisy, but the median across many generations
   is well-concentrated around the true `N`.

3. **Selection vs drift is an SNR problem.** Step 4 added a deterministic
   trend of size `≈ s·p(1-p)` per generation on top of stochastic noise of
   size `√(p(1-p)/(2N))`. The signal-to-noise ratio per generation is
   roughly `s · √(2N · p(1-p))`; over `t` generations the trend grows like
   `s·t·p(1-p)` while the noise grows like `√(t·p(1-p)/(2N))`, so signal /
   noise ∝ **`√(Nₑ · s² · t)`**. The combination `Nₑ · s` is the diffusion-
   approximation **scaled selection coefficient**, the population-genetics
   analogue of an SNR in dB: when `Nₑ · s ≫ 1` selection dominates; when
   `Nₑ · s ≪ 1` the trajectory is effectively neutral and selection is
   invisible behind drift noise.

Everything coalescent-, PSMC-, and selection-scan-flavoured in the rest of
Lecture 12 is built on top of the kernel you just simulated: it's the same
process run backward (coalescent), the same process with state-dependent
rates (selection scans), or the same process averaged over time (Nₑ
inference from real data).
"""


# ---------------------------------------------------------------------------
# Assemble + write
# ---------------------------------------------------------------------------


def build():
    nb = new_notebook()
    nb.cells = [
        md(TITLE_MD),
        md(AHA_MD),
        code(PREAMBLE),
        code(IMPORTS),

        md(STEP1_MD),
        code(STEP1_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP5_SOLUTION),

        md(SELFCHECK_MD),
        code(SELFCHECK),

        md(EE_MD),
    ]
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
        "colab": {"provenance": [], "toc_visible": True},
    }
    return nb


if __name__ == "__main__":
    nb = build()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "exercise.ipynb")
    with open(out_path, "w", encoding="utf-8") as f:
        nbformat.write(nb, f)
    nb2 = nbformat.read(out_path, as_version=4)
    n_md = sum(1 for c in nb2.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb2.cells if c.cell_type == "code")
    n_hidden = sum(
        1 for c in nb2.cells
        if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden")
    )
    print(f"Wrote {out_path}")
    print(f"  cells: {len(nb2.cells)} total  ({n_md} md, {n_code} code, {n_hidden} hidden-solution)")
