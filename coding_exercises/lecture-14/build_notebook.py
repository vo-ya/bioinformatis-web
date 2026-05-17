"""Build exercise.ipynb for L14 — Data Engineering, File Formats, and Reproducibility.

Run from this directory:
    python3 build_notebook.py

Emits exercise.ipynb. Re-running overwrites the file.

The exercise is recast from "full Snakemake pipeline" (too slow on Colab) into
two pieces that together capture the L14 aha:

  (a) GIAB-style VCF benchmarking — load a synthetic "predicted" VCF + a
      "truth" VCF, match by (CHROM, POS, REF, ALT), compute TP/FP/FN, sweep
      QUAL to draw a precision-recall / ROC curve.

  (b) Workflow-DAG resumption logic — implement a content-addressable cache,
      run a 4-task linear DAG, manually invalidate one task's output, show
      which downstream tasks must re-run, and observe wall-clock savings.
"""

from __future__ import annotations

import os
# Make the shared Colab-form helper importable from the parent dir.
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from apply_colab_form import apply_colab_form  # noqa: E402

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


TITLE_MD = """# L14 — Data Engineering, File Formats, and Reproducibility

In this exercise you build the two engineering primitives that sit under
every modern bioinformatics analysis:

1. **GIAB-style VCF benchmarking** — load a "predicted" and a "truth" VCF,
   match variants by `(CHROM, POS, REF, ALT)`, tabulate TP/FP/FN, sweep the
   QUAL threshold to draw a precision-recall curve. This is exactly what
   `hap.py` computes when a paper says "97 % SNV recall on GIAB HG002."
2. **Workflow-DAG resumption logic** — implement a tiny content-addressable
   cache (input hash → output blob), wire four tasks into a linear DAG,
   manually invalidate one task, and observe that everything downstream
   re-runs while everything upstream is skipped. This is what Snakemake,
   Nextflow, and `make` do under the hood.

The notebook is fully self-contained: two VCFs are generated in-notebook
with NumPy from deterministic seeds, and the workflow tasks are toy
functions that simulate compute by sleeping. An optional extension cell
attempts to fetch a real GIAB chr20 slice; if the network call fails the
synthetic path independently passes self-check.
"""


AHA_MD = """> **Aha.** Workflow managers are **content-addressable caches plus a
> topological scheduler**. The cache key is `hash(inputs, code, tool
> version)`. If the key hits, skip the task. If not, run it. Reproducibility
> falls out of declarative dependency graphs: the runtime, not the
> scientist, decides what to re-run after a failure.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 matplotlib==3.8.4 -q
"""


IMPORTS = """import io
import os
import json
import time
import hashlib
import pickle
import tempfile
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
np.random.seed(42)

# Workspace for the cache demo — wiped on each run.
CACHE_DIR = os.path.join(tempfile.gettempdir(), "l14_cache")
os.makedirs(CACHE_DIR, exist_ok=True)
for f in os.listdir(CACHE_DIR):
    os.remove(os.path.join(CACHE_DIR, f))

print(f"Workspace ready. Cache dir: {CACHE_DIR}")
"""


STEP1_MD = """## Step 1 (8 min) — Generate two VCFs and inspect the schema

Real GIAB benchmarks are too big for a 5-minute Colab budget, so we
synthesise two VCFs that have the same schema as a real call set:

- **`truth.vcf`** — 300 variants, the "ground truth".
- **`pred.vcf`** — 500 variants. Built so that ~70 % of the truth set
  overlaps (true positives), with the remainder split across false negatives
  (in truth, missed by pred) and false positives (in pred, not in truth).
  Each pred variant also carries a `QUAL` score (the caller's confidence).

The schema follows VCF v4.2:

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    101  .   A    G    35.0  PASS    DP=42
```

This is the same data model that `bcftools`, `hap.py`, and every variant
caller consumes. The point of this step is just to confirm you can parse
the header + record lines into a DataFrame.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Synthesise truth + pred VCFs and parse them into DataFrames.
# ----------------------------------------------------------------------

N_TRUTH = 300
N_PRED  = 500
RECALL_TARGET = 0.70   # fraction of truth variants we want pred to recover


def synth_vcfs(n_truth: int, n_pred: int, recall_target: float, seed: int = 42) -> tuple[str, str]:
    """Return (truth_vcf_text, pred_vcf_text). Both deterministic given seed."""
    # TODO:
    # 1. Pick a chromosome panel (e.g. ["chr20", "chr21", "chr22"]).
    # 2. Sample n_truth (chrom, pos, ref, alt) tuples with no duplicates.
    # 3. Take floor(recall_target * n_truth) as "shared" truth variants;
    #    sample the rest of pred as fresh tuples (the false positives).
    # 4. Assign each pred variant a QUAL drawn from a mixture:
    #      shared truth variants -> Normal(40, 8)
    #      false positives       -> Normal(20, 8)
    #    This gives QUAL real discriminative power.
    # 5. Emit header lines + tab-delimited record lines in VCF v4.2 form.
    raise NotImplementedError


def parse_vcf(text: str) -> pd.DataFrame:
    """Parse VCF text -> DataFrame with columns [CHROM, POS, REF, ALT, QUAL]."""
    # TODO: skip lines starting with '##'; the first '#CHROM' line is the column
    # header. Each subsequent line is a tab-split record. Coerce POS -> int,
    # QUAL -> float (treat '.' as NaN).
    raise NotImplementedError


# truth_text, pred_text = synth_vcfs(N_TRUTH, N_PRED, RECALL_TARGET)
# truth = parse_vcf(truth_text)
# pred  = parse_vcf(pred_text)
'''


SOLN_HEADER = """*Click ▶ to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution — Step 1.

N_TRUTH = 300
N_PRED  = 500
RECALL_TARGET = 0.70
CHROMS = ["chr20", "chr21", "chr22"]

VCF_HEADER = """##fileformat=VCFv4.2
##source=L14-synth
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO"""


def _sample_unique_sites(n: int, rng: np.random.Generator) -> list[tuple[str, int, str, str]]:
    """Sample n unique (chrom, pos, ref, alt) tuples."""
    seen = set()
    out = []
    bases = ["A", "C", "G", "T"]
    while len(out) < n:
        chrom = CHROMS[int(rng.integers(0, len(CHROMS)))]
        pos = int(rng.integers(1, 50_000_000))
        ref = bases[int(rng.integers(0, 4))]
        # ALT != REF
        alt = bases[int(rng.integers(0, 4))]
        while alt == ref:
            alt = bases[int(rng.integers(0, 4))]
        key = (chrom, pos, ref, alt)
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def synth_vcfs(n_truth: int, n_pred: int, recall_target: float, seed: int = 42):
    rng = np.random.default_rng(seed)

    # Truth set.
    truth_sites = _sample_unique_sites(n_truth, rng)

    # Shared TP slice + fresh FP slice.
    n_shared = int(np.floor(recall_target * n_truth))
    rng.shuffle(truth_sites)  # in-place; shuffles tuples list
    shared = list(truth_sites[:n_shared])

    # Distinct FP sites (avoid clashing with truth).
    truth_keys = set(truth_sites)
    fp_sites = []
    bases = ["A", "C", "G", "T"]
    while len(fp_sites) + n_shared < n_pred:
        chrom = CHROMS[int(rng.integers(0, len(CHROMS)))]
        pos = int(rng.integers(1, 50_000_000))
        ref = bases[int(rng.integers(0, 4))]
        alt = bases[int(rng.integers(0, 4))]
        while alt == ref:
            alt = bases[int(rng.integers(0, 4))]
        key = (chrom, pos, ref, alt)
        if key in truth_keys:
            continue
        truth_keys.add(key)
        fp_sites.append(key)

    # QUAL mixture: TPs are "good", FPs are "weak".
    tp_quals = rng.normal(40.0, 8.0, size=n_shared).clip(0, 99)
    fp_quals = rng.normal(20.0, 8.0, size=len(fp_sites)).clip(0, 99)

    pred_records = [(s, q) for s, q in zip(shared, tp_quals)]
    pred_records += [(s, q) for s, q in zip(fp_sites, fp_quals)]
    rng.shuffle(pred_records)  # shuffle the line order

    def fmt_truth(site):
        chrom, pos, ref, alt = site
        # Truth VCFs typically have QUAL = '.', FILTER = PASS.
        return f"{chrom}\\t{pos}\\t.\\t{ref}\\t{alt}\\t.\\tPASS\\tDP=50"

    def fmt_pred(site, qual):
        chrom, pos, ref, alt = site
        return f"{chrom}\\t{pos}\\t.\\t{ref}\\t{alt}\\t{qual:.1f}\\tPASS\\tDP=42"

    # Sort truth by (chrom, pos) for prettiness.
    truth_sorted = sorted(truth_sites, key=lambda s: (CHROMS.index(s[0]), s[1]))
    truth_text = VCF_HEADER + "\\n" + "\\n".join(fmt_truth(s) for s in truth_sorted)

    pred_sorted = sorted(pred_records, key=lambda r: (CHROMS.index(r[0][0]), r[0][1]))
    pred_text = VCF_HEADER + "\\n" + "\\n".join(fmt_pred(s, q) for s, q in pred_sorted)

    return truth_text, pred_text


def parse_vcf(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            continue
        parts = line.split("\\t")
        # Standard VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO ...
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        try:
            qual = float(parts[5])
        except ValueError:
            qual = np.nan
        rows.append((chrom, pos, ref, alt, qual))
    return pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT", "QUAL"])


truth_text, pred_text = synth_vcfs(N_TRUTH, N_PRED, RECALL_TARGET)
truth = parse_vcf(truth_text)
pred = parse_vcf(pred_text)

print("Truth VCF (first 4 lines):")
print("\\n".join(truth_text.splitlines()[:4]))
print()
print(f"truth rows: {len(truth)}    pred rows: {len(pred)}")
print()
print("Truth head:")
print(truth.head().to_string(index=False))
print()
print("Pred head (QUAL is the caller's confidence):")
print(pred.head().to_string(index=False))
'''


STEP2_MD = """## Step 2 (12 min) — TP / FP / FN matching by `(CHROM, POS, REF, ALT)`

The simplest variant-benchmark matching rule: two records are the same
variant iff their `(CHROM, POS, REF, ALT)` tuples are equal. (Real
benchmarks also handle indel left-alignment and complex-variant
decomposition; we deliberately skip that — it is `bcftools norm` territory,
not the L14 aha.)

Implement matching as a set intersection:

- **TP** = pred ∩ truth — records present in both.
- **FP** = pred − truth — predicted variants with no truth support.
- **FN** = truth − pred — truth variants the caller missed.

Then print the confusion-matrix-style summary.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — TP / FP / FN matching and confusion matrix.
# ----------------------------------------------------------------------


def variant_key(row) -> tuple[str, int, str, str]:
    """A 4-tuple is unique enough for our synthetic call set."""
    return (row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])


def match_vcfs(truth: pd.DataFrame, pred: pd.DataFrame) -> dict:
    """Return {'tp': DataFrame, 'fp': DataFrame, 'fn': DataFrame}."""
    # TODO:
    # 1. Build set of variant_key for truth and for pred.
    # 2. Compute intersection (TP), pred - truth (FP), truth - pred (FN).
    # 3. Return DataFrames sliced by membership.
    raise NotImplementedError


# m = match_vcfs(truth, pred)
# print(f"TP={len(m['tp'])}  FP={len(m['fp'])}  FN={len(m['fn'])}")
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.


def variant_key(row) -> tuple[str, int, str, str]:
    return (row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])


def match_vcfs(truth: pd.DataFrame, pred: pd.DataFrame) -> dict:
    truth_keys = set(map(tuple, truth[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None)))
    pred_keys  = set(map(tuple, pred[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None)))

    tp_keys = truth_keys & pred_keys
    fp_keys = pred_keys - truth_keys
    fn_keys = truth_keys - pred_keys

    def slice_by_keys(df, keys):
        if not keys:
            return df.iloc[0:0].copy()
        df_keys = list(map(tuple, df[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None)))
        mask = [k in keys for k in df_keys]
        return df[mask].reset_index(drop=True)

    return {
        "tp": slice_by_keys(pred, tp_keys),    # carry pred's QUAL into the TP slice
        "fp": slice_by_keys(pred, fp_keys),
        "fn": slice_by_keys(truth, fn_keys),
    }


m = match_vcfs(truth, pred)
tp, fp, fn = len(m["tp"]), len(m["fp"]), len(m["fn"])
print(f"Confusion-matrix counts (no QUAL filter):")
print(f"  TP (in both)             = {tp}")
print(f"  FP (pred only)           = {fp}")
print(f"  FN (truth only, missed)  = {fn}")
print()
print("Sanity: TP + FN should equal the truth set size.")
print(f"  TP + FN = {tp + fn}   (expected {len(truth)})")
print(f"  TP + FP = {tp + fp}   (expected {len(pred)})")
print()
print("TP head (QUAL came from pred):")
print(m["tp"].head().to_string(index=False))
'''


STEP3_MD = """## Step 3 (15 min) — Precision, recall, F1, and the QUAL-swept ROC

With the confusion matrix in hand:

$$\\text{precision} = \\frac{TP}{TP + FP}, \\qquad
\\text{recall} = \\frac{TP}{TP + FN}, \\qquad
F_1 = \\frac{2 \\cdot \\text{precision} \\cdot \\text{recall}}{\\text{precision} + \\text{recall}}$$

Then **sweep the QUAL threshold**. A caller with a well-calibrated QUAL
should let you trade precision for recall: at high QUAL you keep only the
confident calls (high precision, low recall); at low QUAL you keep
everything (high recall, low precision). The trace of (recall, precision)
as QUAL drops from 99 to 0 is the **precision-recall curve** — the EE
analogue is an ROC curve on a detection problem.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — precision, recall, F1; sweep QUAL threshold for a PR curve.
# ----------------------------------------------------------------------


def precision_recall_f1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    """Return (precision, recall, F1). Define 0/0 = 0 for clarity."""
    # TODO
    raise NotImplementedError


def sweep_qual(truth: pd.DataFrame, pred: pd.DataFrame,
               thresholds: np.ndarray) -> pd.DataFrame:
    """For each threshold q, filter pred by QUAL >= q and recompute (P, R, F1)."""
    # TODO: at each q,
    #   1. filter pred to QUAL >= q,
    #   2. re-match against truth,
    #   3. record (q, tp, fp, fn, precision, recall, F1).
    # Return a DataFrame sorted by threshold ascending.
    raise NotImplementedError


# qs = np.linspace(0, 60, 31)
# curve = sweep_qual(truth, pred, qs)
# Plot precision + recall vs QUAL on one axis, F1 on a twin axis.
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.


def precision_recall_f1(tp: int, fp: int, fn: int):
    p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0.0
    return p, r, f1


def sweep_qual(truth: pd.DataFrame, pred: pd.DataFrame, thresholds: np.ndarray) -> pd.DataFrame:
    rows = []
    truth_keys = set(map(tuple, truth[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None)))
    for q in thresholds:
        sub = pred[pred["QUAL"] >= q]
        pred_keys = set(map(tuple, sub[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None)))
        tp = len(truth_keys & pred_keys)
        fp = len(pred_keys - truth_keys)
        fn = len(truth_keys - pred_keys)
        p, r, f1 = precision_recall_f1(tp, fp, fn)
        rows.append({
            "QUAL_min": float(q),
            "TP": tp, "FP": fp, "FN": fn,
            "precision": p, "recall": r, "F1": f1,
        })
    return pd.DataFrame(rows).sort_values("QUAL_min").reset_index(drop=True)


qs = np.linspace(0, 60, 31)
curve = sweep_qual(truth, pred, qs)

# Pick the QUAL that maximises F1.
best_idx = curve["F1"].idxmax()
best = curve.loc[best_idx]
print(f"Best F1 = {best['F1']:.3f} at QUAL >= {best['QUAL_min']:.1f}")
print(f"  precision = {best['precision']:.3f}    recall = {best['recall']:.3f}")
print(f"  TP={int(best['TP'])} FP={int(best['FP'])} FN={int(best['FN'])}")

# --- Plot 1: precision / recall / F1 vs QUAL threshold --------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4))

ax[0].plot(curve["QUAL_min"], curve["precision"], label="precision", color="#1f77b4", lw=2)
ax[0].plot(curve["QUAL_min"], curve["recall"],    label="recall",    color="#d62728", lw=2)
ax[0].plot(curve["QUAL_min"], curve["F1"],        label="F1",        color="#2ca02c", lw=2, ls="--")
ax[0].axvline(best["QUAL_min"], color="grey", ls=":", lw=1)
ax[0].set_xlabel("QUAL threshold")
ax[0].set_ylabel("metric")
ax[0].set_ylim(0, 1.05)
ax[0].set_title("Sweep QUAL: precision vs recall trade-off")
ax[0].legend(loc="lower left")

# --- Plot 2: precision-recall curve ---------------------------------------
ax[1].plot(curve["recall"], curve["precision"], "o-", color="#9467bd", lw=2)
ax[1].plot([best["recall"]], [best["precision"]], "*", color="#2ca02c",
           markersize=16, label=f"best F1 @ QUAL>={best['QUAL_min']:.0f}")
ax[1].set_xlabel("recall")
ax[1].set_ylabel("precision")
ax[1].set_xlim(0, 1.05)
ax[1].set_ylim(0, 1.05)
ax[1].set_title("Precision-recall curve")
ax[1].legend(loc="lower left")

plt.tight_layout()
plt.show()

# Save these for the self-check.
pr_curve = curve
best_metrics = best
'''


STEP4_MD = """## Step 4 (15 min) — Workflow-DAG resumption with a content-addressable cache

Workflow managers (Snakemake, Nextflow, `make`) are two ideas glued
together:

1. **A dependency DAG** — tasks declare their inputs and outputs; the
   runtime topologically sorts them.
2. **A content-addressable cache** — each task's output is stored under a
   key derived from `hash(inputs, code, tool version)`. If the key exists,
   skip the task; otherwise run it and write the result to the cache.

Implement a tiny version. We have **four tasks** wired in a chain:

```
load_truth → load_pred → match → metrics
```

Each task `sleep`s briefly to simulate compute. Run the whole DAG once and
measure wall-clock. Then **invalidate** `load_pred` (e.g. the caller was
re-run with new parameters). Re-run the DAG. Observe that:

- `load_truth` is a **cache hit** — its inputs did not change.
- `load_pred`, `match`, and `metrics` are **cache misses** — anything
  downstream of an invalidated task is invalid too.

This is the same logic behind Snakemake's `-R` flag and Nextflow's
`-resume`.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Content-addressable cache + 4-task DAG with resume logic.
# ----------------------------------------------------------------------


@dataclass
class Task:
    name: str
    func: Callable                  # python callable, returns the artifact
    inputs: dict                    # named upstream task results / static params
    code_version: str = "1.0"       # bump to invalidate without touching inputs


def content_hash(payload: Any) -> str:
    """Stable SHA-256 hex digest of a JSON-serialisable payload."""
    # TODO: json.dumps with sort_keys=True; encode utf-8; sha256().hexdigest().
    raise NotImplementedError


def cache_key(task: Task, upstream_keys: dict[str, str]) -> str:
    """Combine task name + code version + upstream cache keys into one key."""
    # TODO: hash the dict {"name": task.name, "code": task.code_version,
    # "upstream": upstream_keys, "static_inputs": <non-task inputs>}.
    raise NotImplementedError


def run_dag(tasks: list[Task], cache_dir: str = CACHE_DIR, log: list = None) -> dict:
    """Run tasks in declared order with cache lookup. Return {name: result}."""
    # TODO:
    # For each task in order:
    #   - resolve upstream inputs from the results dict
    #   - compute cache_key from this task + its upstream tasks' cache keys
    #     (NOT their output hashes — that's how real workflow managers do it,
    #     so a code-version bump invalidates downstream even if the output
    #     bytes happen to match)
    #   - if cache_dir/<key>.pkl exists -> load it (cache hit; log "HIT")
    #   - else -> run task.func, pickle to cache (cache miss; log "MISS")
    # Return (results_dict, log_list).
    raise NotImplementedError


# Define the 4 tasks (load_truth, load_pred, match, metrics) and wire them up.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


@dataclass
class Task:
    name: str
    func: Callable
    inputs: dict
    code_version: str = "1.0"


def content_hash(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _summarise_for_hash(obj) -> Any:
    """JSON-friendly summary used in cache keys (DataFrames -> shape + checksum)."""
    if isinstance(obj, pd.DataFrame):
        # Cheap stable fingerprint without serialising the whole frame.
        return {
            "shape": list(obj.shape),
            "cols": list(map(str, obj.columns)),
            "checksum": int(pd.util.hash_pandas_object(obj, index=True).sum() % (10**16)),
        }
    if isinstance(obj, dict):
        return {k: _summarise_for_hash(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_summarise_for_hash(v) for v in obj]
    return obj


def cache_key(task: Task, upstream_keys: dict[str, str]) -> str:
    # Separate task-result inputs (referenced by name) from raw static inputs.
    static = {k: _summarise_for_hash(v) for k, v in task.inputs.items()
              if not (isinstance(v, str) and v.startswith("@"))}
    payload = {
        "name":     task.name,
        "code":     task.code_version,
        "upstream": dict(sorted(upstream_keys.items())),
        "static":   static,
    }
    return content_hash(payload)


def run_dag(tasks: list[Task], cache_dir: str = CACHE_DIR, log: list = None):
    if log is None:
        log = []
    results: dict[str, Any] = {}
    # Real workflow managers (Snakemake, Nextflow) key downstream tasks off
    # the upstream task's *cache key* (a hash of its declared inputs + code
    # version), not off the upstream task's output content. This is what makes
    # "bump tool version, rerun everything downstream" work even when the
    # outputs happen to be byte-identical.
    upstream_keys: dict[str, str] = {}

    for task in tasks:
        ups = {}
        resolved_inputs = {}
        for k, v in task.inputs.items():
            if isinstance(v, str) and v.startswith("@"):
                ref = v[1:]
                resolved_inputs[k] = results[ref]
                ups[ref] = upstream_keys[ref]
            else:
                resolved_inputs[k] = v

        key = cache_key(task, ups)
        path = os.path.join(cache_dir, f"{task.name}_{key[:16]}.pkl")

        t0 = time.time()
        if os.path.exists(path):
            with open(path, "rb") as f:
                value = pickle.load(f)
            status = "HIT "
        else:
            value = task.func(**resolved_inputs)
            with open(path, "wb") as f:
                pickle.dump(value, f)
            status = "MISS"
        dt = time.time() - t0

        upstream_keys[task.name] = key
        results[task.name] = value
        log.append({"task": task.name, "status": status, "seconds": dt, "key": key[:12]})

    return results, log


# ---- Task implementations -------------------------------------------------

def task_load_truth(seed: int = 42, n: int = N_TRUTH, recall: float = RECALL_TARGET) -> pd.DataFrame:
    time.sleep(0.5)  # simulate "loading from disk"
    text, _ = synth_vcfs(n, N_PRED, recall, seed=seed)
    return parse_vcf(text)


def task_load_pred(seed: int = 42, n_truth: int = N_TRUTH, n_pred: int = N_PRED,
                   recall: float = RECALL_TARGET) -> pd.DataFrame:
    time.sleep(0.5)
    _, text = synth_vcfs(n_truth, n_pred, recall, seed=seed)
    return parse_vcf(text)


def task_match(truth: pd.DataFrame, pred: pd.DataFrame) -> dict:
    time.sleep(0.5)
    return match_vcfs(truth, pred)


def task_metrics(match: dict) -> dict:
    time.sleep(0.5)
    tp = len(match["tp"]); fp = len(match["fp"]); fn = len(match["fn"])
    p, r, f1 = precision_recall_f1(tp, fp, fn)
    return {"TP": tp, "FP": fp, "FN": fn, "precision": p, "recall": r, "F1": f1}


tasks = [
    Task("load_truth", task_load_truth, {"seed": 42, "n": N_TRUTH, "recall": RECALL_TARGET}),
    Task("load_pred",  task_load_pred,  {"seed": 42, "n_truth": N_TRUTH, "n_pred": N_PRED, "recall": RECALL_TARGET}),
    Task("match",      task_match,      {"truth": "@load_truth", "pred": "@load_pred"}),
    Task("metrics",    task_metrics,    {"match": "@match"}),
]

# --- Run 1: cold cache ----------------------------------------------------
print("Run 1 — cold cache")
print("-" * 50)
results, log1 = run_dag(tasks)
total1 = sum(e["seconds"] for e in log1)
for e in log1:
    print(f"  {e['task']:<12}  {e['status']}  {e['seconds']*1000:6.0f} ms  key={e['key']}")
print(f"  TOTAL: {total1*1000:.0f} ms")
print()
print("metrics output:")
print(json.dumps(results["metrics"], indent=2))

# --- Run 2: warm cache, nothing changed -----------------------------------
print()
print("Run 2 — warm cache, no changes")
print("-" * 50)
_, log2 = run_dag(tasks)
total2 = sum(e["seconds"] for e in log2)
for e in log2:
    print(f"  {e['task']:<12}  {e['status']}  {e['seconds']*1000:6.0f} ms  key={e['key']}")
print(f"  TOTAL: {total2*1000:.0f} ms   (speedup {total1/max(total2,1e-9):.0f}x)")

# --- Run 3: invalidate load_pred (e.g. caller re-run with new params) ------
print()
print("Run 3 — invalidate load_pred (bump code_version)")
print("-" * 50)
tasks[1] = Task("load_pred", task_load_pred,
                {"seed": 42, "n_truth": N_TRUTH, "n_pred": N_PRED, "recall": RECALL_TARGET},
                code_version="2.0")  # <-- this changes the cache key
_, log3 = run_dag(tasks)
total3 = sum(e["seconds"] for e in log3)
for e in log3:
    print(f"  {e['task']:<12}  {e['status']}  {e['seconds']*1000:6.0f} ms  key={e['key']}")
print(f"  TOTAL: {total3*1000:.0f} ms")
print()
print("Note: load_truth is still a HIT (its inputs are unchanged),")
print("but load_pred, match, and metrics all MISS — invalidation propagates.")

# Save logs for the self-check.
log_run1, log_run2, log_run3 = log1, log2, log3
'''


STEP5_MD = """## Step 5 (10 min) — When does resume actually matter?

A markdown summary cell of when the cache-and-resume machinery pays off, and
an **optional** extension cell that tries to fetch a real GIAB v4.2.1 chr20
slice and re-runs the matching logic against it. The synthetic path above
must independently pass self-check; the GIAB cell prints a friendly notice
if the public URL is unreachable.

**When resume matters most:**

- **Long-running tasks.** A 6-hour alignment that fails at hour 5 should
  resume at hour 5 — not hour 0. The cache pays for itself the first time
  this happens.
- **Large intermediates.** A 200 GB BAM that took 90 minutes to align is
  the kind of artifact you do *not* want to regenerate because a downstream
  bug surfaces.
- **Iterative analysis.** When you change *only the last step* (filter
  threshold, summary plot), everything upstream is unchanged and skipping
  it gives an order-of-magnitude turnaround speedup.
- **Provenance.** Cache keys are a free audit trail: the key changes if
  and only if something that could affect the result changed.

**When it doesn't help:**

- One-off scripts that run in seconds. The cache bookkeeping isn't free.
- Non-deterministic tasks (random seeds left unfixed). The same inputs
  hash the same, but the outputs differ — the cache will silently return
  stale results. The fix: pin seeds and treat them as inputs.
- Tasks that depend on hidden state (system time, network responses).
  The cache key only sees what you declare. Hidden inputs → silent staleness.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Optional: fetch a real GIAB chr20 slice and re-run matching.
# ----------------------------------------------------------------------

# OPTIONAL extension. This cell tries to download a small public GIAB chr20
# slice. If the fetch fails (no network, URL moved, etc.), it prints a
# friendly notice and bails — the synthetic path above is enough to pass the
# self-check below.

GIAB_TRUTH_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
    "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)

def try_fetch_giab(url: str = GIAB_TRUTH_URL, max_bytes: int = 2_000_000) -> str | None:
    """Best-effort: HEAD-and-skip if the host is unreachable; else stream
    the first max_bytes for a quick smoke test only."""
    # TODO (optional): urllib + gzip; return decoded text or None on failure.
    return None


text = try_fetch_giab()
if text is None:
    print("Skipping the GIAB extension: synthetic VCFs already exercise the same logic.")
else:
    print(f"Fetched {len(text)} bytes of GIAB truth VCF (chr20 slice).")
'''


STEP5_SOLUTION = '''# Reference solution — Step 5 (optional GIAB fetch).
import gzip

GIAB_TRUTH_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
    "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)


def try_fetch_giab(url: str = GIAB_TRUTH_URL, max_bytes: int = 2_000_000) -> str | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L14-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            buf = resp.read(max_bytes)
        # File is gzipped; decompress what we have.
        try:
            text = gzip.decompress(buf).decode("utf-8", errors="replace")
        except Exception:
            # Truncated gzip stream — fall back to decoding raw bytes.
            text = buf.decode("utf-8", errors="replace")
        return text
    except Exception as exc:
        print(f"  GIAB fetch failed: {exc}")
        return None


text = try_fetch_giab()
if text is None:
    print("Skipping the GIAB extension: synthetic VCFs already exercise the same logic.")
else:
    # Parse only the chr20 records we successfully decompressed.
    head_lines = [ln for ln in text.splitlines() if ln.startswith("##")][:8]
    print("GIAB header sample:")
    for ln in head_lines:
        print(f"  {ln}")
    # Try to parse — handle the case where the trailing record was cut off.
    safe_text = "\\n".join(text.splitlines()[:-1])
    try:
        giab = parse_vcf(safe_text)
        print()
        print(f"Parsed {len(giab):,} GIAB records (chr20-prefix subset, may be truncated)")
        print(giab.head().to_string(index=False))
    except Exception as exc:
        print(f"  Parse failed (truncated stream): {exc}")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces. If you ran the reference
solutions above they all pass; if you wrote your own and one fails, revisit
the matching step.
"""


SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. VCF generation + parse round-trips to the right shape.
assert len(truth) == N_TRUTH, f"truth length {len(truth)} != {N_TRUTH}"
assert len(pred)  == N_PRED,  f"pred length {len(pred)} != {N_PRED}"
assert set(truth.columns) == {"CHROM", "POS", "REF", "ALT", "QUAL"}

# 2. TP + FN = truth count, TP + FP = pred count (no QUAL filter).
m = match_vcfs(truth, pred)
tp, fp, fn = len(m["tp"]), len(m["fp"]), len(m["fn"])
assert tp + fn == N_TRUTH, f"TP+FN={tp+fn} should equal truth size {N_TRUTH}"
assert tp + fp == N_PRED,  f"TP+FP={tp+fp} should equal pred size {N_PRED}"

# 3. Recall at QUAL=0 should land near the design target of 0.70.
p0, r0, _ = precision_recall_f1(tp, fp, fn)
assert 0.60 <= r0 <= 0.80, f"recall@QUAL=0 = {r0:.2f} outside [0.60, 0.80]"

# 4. Sweeping QUAL must improve precision monotonically (more or less) and
#    drop recall monotonically. Check the endpoints.
low  = pr_curve.iloc[0]
high = pr_curve.iloc[-1]
assert high["precision"] >= low["precision"] - 0.02, \\
    f"high-QUAL precision {high['precision']:.2f} should be >= low-QUAL {low['precision']:.2f}"
assert high["recall"] <= low["recall"] + 0.02, \\
    f"high-QUAL recall {high['recall']:.2f} should be <= low-QUAL {low['recall']:.2f}"

# 5. F1 maximiser should be away from both endpoints (i.e. the trade-off is real).
best_q = best_metrics["QUAL_min"]
assert 5.0 < best_q < 55.0, f"best F1 at QUAL={best_q} is suspicious for this design"
assert best_metrics["F1"] > 0.70, f"best F1 = {best_metrics['F1']:.2f} should exceed 0.70"

# 6. DAG resume: run 2 should be all HITs; run 3 should HIT load_truth and
#    MISS the other three (load_pred invalidated, propagation downstream).
hits_run2 = [e for e in log_run2 if e["status"] == "HIT "]
assert len(hits_run2) == 4, f"Run 2 had {len(hits_run2)} hits, expected 4 (full warm cache)"

statuses_run3 = {e["task"]: e["status"].strip() for e in log_run3}
assert statuses_run3["load_truth"] == "HIT",  "load_truth should still be cached in run 3"
assert statuses_run3["load_pred"]  == "MISS", "load_pred was invalidated; should MISS"
assert statuses_run3["match"]      == "MISS", "match is downstream of load_pred; should MISS"
assert statuses_run3["metrics"]    == "MISS", "metrics is downstream of match; should MISS"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — caches, hashes, and dataflow scheduling

You built the same primitive that lives under Snakemake, Nextflow, GNU
`make`, Apache Spark, and TensorFlow's `tf.function`:

1. **Content-addressable storage.** The cache key is `hash(inputs, code,
   tool version)`. It is **content-addressable** because the address is
   derived from the *content*, not from the file path. Git uses the same
   trick on commits, IPFS on every blob, BLAKE3 in modern build systems.
2. **Dependency DAG.** Tasks declare their inputs; the runtime builds a
   DAG and runs it in topological order. The DAG also tells the runtime
   which tasks can run in **parallel** (sibling subtrees), and which
   downstream tasks become **invalid** when an upstream task changes.
3. **Detection-theoretic benchmarking.** Precision / recall / F1 are the
   same statistics a radar engineer uses for `P_d` and `P_FA`. The QUAL
   sweep is a sliding **detection threshold** on a soft score; the
   precision-recall curve is the operating-characteristic curve of the
   variant caller.

So when a workflow paper says "nf-core's resume cut our runtime from 4 hours
to 4 minutes after a config tweak," it isn't magic. It is hash + DAG +
disk. The same machinery the EE world has been calling **memoisation** for
fifty years.
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
        md(SOLN_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(SOLN_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(SOLN_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(SOLN_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(SOLN_HEADER),
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
        apply_colab_form(nb)
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
