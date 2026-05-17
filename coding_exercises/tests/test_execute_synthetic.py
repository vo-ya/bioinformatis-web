"""Execute the synthetic-only notebooks end-to-end.

Targets only the 17 lectures that generate all their data in-cell with NumPy:
L01, L03–L06, L09–L14, L18, L22–L25, L27. The other 10 lectures fetch external
data (Scanpy, HuggingFace, UniProt, …) and are covered by ``test_fetches.py``
plus the on-demand full execution suite.

Per-notebook procedure:
  1. Load JSON, drop pip-install cells (would shell out, slow + brittle).
  2. Execute every remaining cell with ``allow_errors=True`` so a stubbed
     ``NotImplementedError`` in a TODO cell doesn't abort the kernel.
  3. Walk cell outputs and fail on any error *other* than NotImplementedError
     in a visible (non-hidden-solution) cell. NotImplementedError in a TODO
     stub is the *designed* failure mode — students hit it when they reach
     the cell they need to fill in.

This mirrors a student's run: they execute all cells, hit NotImplementedError
on the TODO they need to write, then expand the hidden solution which
redefines the function and reruns the downstream visualization.

Usage::

    python coding_exercises/tests/test_execute_synthetic.py
    python coding_exercises/tests/test_execute_synthetic.py -l 22  # one lecture
    python -m pytest coding_exercises/tests/test_execute_synthetic.py -v

Dependencies:
    pip install nbclient nbformat ipykernel
    pip install numpy scipy pandas matplotlib networkx statsmodels scikit-learn biopython
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass

from _helpers import (  # type: ignore[import-not-found]
    SYNTHETIC_LECTURES,
    discover_notebooks,
    filter_cells_for_execution,
    is_hidden_solution,
    is_synthetic,
    lecture_number,
    load_notebook,
)

PER_CELL_TIMEOUT_SEC = 180  # generous — most synthetic cells finish in <5 s


def _import_nbclient():
    try:
        import nbformat
        from nbclient import NotebookClient
        from nbclient.exceptions import CellExecutionError
        return nbformat, NotebookClient, CellExecutionError
    except ImportError as e:
        raise RuntimeError(
            "nbclient + nbformat + ipykernel required to run notebooks. "
            "Install with: pip install nbclient nbformat ipykernel"
        ) from e


@dataclass
class ExecutionResult:
    path: str
    ok: bool
    elapsed_s: float
    n_cells_executed: int
    n_cells_skipped: int
    message: str = ""


def execute_notebook(path: str) -> ExecutionResult:
    nbformat, NotebookClient, _CellExecutionError = _import_nbclient()

    started = time.monotonic()
    nb_dict = load_notebook(path)
    filtered_cells, counts = filter_cells_for_execution(nb_dict)
    nb_dict["cells"] = filtered_cells
    nb = nbformat.from_dict(nb_dict)

    # Resolve the working dir to the notebook's directory so any relative-path
    # file IO inside the notebook resolves correctly.
    resources = {"metadata": {"path": os.path.dirname(path)}}
    # allow_errors=True keeps the kernel running past a cell error so we can
    # inspect *all* cells' outputs and distinguish designed NotImplementedError
    # stubs from real bugs.
    client = NotebookClient(nb, timeout=PER_CELL_TIMEOUT_SEC,
                            allow_errors=True, resources=resources)
    try:
        client.execute()
    except Exception as e:  # noqa: BLE001 — kernel crash, timeout, etc.
        elapsed = time.monotonic() - started
        return ExecutionResult(
            path=path, ok=False, elapsed_s=elapsed,
            n_cells_executed=len(filtered_cells),
            n_cells_skipped=counts.get("pip-skipped", 0),
            message=f"{type(e).__name__}: {e}",
        )

    # Walk cell outputs. A cell error is acceptable iff:
    #   (a) it's NotImplementedError, AND
    #   (b) the cell is a visible TODO cell (not a hidden-solution cell).
    # Anything else is a real bug.
    real_errors: list[str] = []
    for idx, cell in enumerate(nb.cells):
        if cell.get("cell_type") != "code":
            continue
        for out in cell.get("outputs", []):
            if out.get("output_type") != "error":
                continue
            ename = out.get("ename", "?")
            evalue = out.get("evalue", "")
            allowed = (ename == "NotImplementedError"
                       and not is_hidden_solution(cell))
            if not allowed:
                real_errors.append(f"cell[{idx}] {ename}: {evalue.splitlines()[0] if evalue else ''}".strip())

    elapsed = time.monotonic() - started
    return ExecutionResult(
        path=path, ok=not real_errors, elapsed_s=elapsed,
        n_cells_executed=len(filtered_cells),
        n_cells_skipped=counts.get("pip-skipped", 0),
        message="; ".join(real_errors),
    )


def synthetic_notebooks() -> list[str]:
    return [p for p in discover_notebooks() if is_synthetic(p)]


# ─── pytest hook ─────────────────────────────────────────────────────────────

def pytest_generate_tests(metafunc):  # noqa: D401
    if "synthetic_path" in metafunc.fixturenames:
        paths = synthetic_notebooks()
        ids = [f"L{lecture_number(p):02d}" for p in paths]
        metafunc.parametrize("synthetic_path", paths, ids=ids)


def test_execute(synthetic_path: str):
    r = execute_notebook(synthetic_path)
    assert r.ok, f"{os.path.basename(os.path.dirname(synthetic_path))}: {r.message}"


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lecture", "-l", action="append", type=int,
                    help="restrict to lecture number(s); may repeat (e.g. -l 1 -l 3)")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="show full traceback for failures")
    args = ap.parse_args(argv)

    targets = synthetic_notebooks()
    if args.lecture:
        wanted = set(args.lecture)
        for n in wanted:
            if n not in SYNTHETIC_LECTURES:
                print(f"  warn: L{n:02d} is not in the synthetic set; "
                      f"synthetic = {SYNTHETIC_LECTURES}", file=sys.stderr)
        targets = [p for p in targets if lecture_number(p) in wanted]

    if not targets:
        print("No notebooks to execute.", file=sys.stderr)
        return 2

    print(f"Executing {len(targets)} synthetic notebook(s) "
          f"(per-cell timeout {PER_CELL_TIMEOUT_SEC}s)\n")
    n_ok = 0
    total_elapsed = 0.0
    results: list[ExecutionResult] = []
    for path in targets:
        r = execute_notebook(path)
        results.append(r)
        total_elapsed += r.elapsed_s
        tag = "OK " if r.ok else "FAIL"
        n_ok += int(r.ok)
        line = (f"  [{tag}] L{lecture_number(path):02d}  "
                f"{r.elapsed_s:6.1f}s  cells: "
                f"{r.n_cells_executed} run / {r.n_cells_skipped} skipped")
        print(line)
        if not r.ok:
            print(f"        {r.message}")

    n_fail = len(results) - n_ok
    print(f"\n{n_ok}/{len(results)} ok, {n_fail} failed  "
          f"(total elapsed {total_elapsed:.1f}s)")
    if n_fail:
        print("\nFailures:")
        for r in results:
            if not r.ok:
                print(f"  L{lecture_number(r.path):02d}  {r.message}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
