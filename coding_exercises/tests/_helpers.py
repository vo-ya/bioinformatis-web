"""Shared helpers for the notebook test suites.

Conventions inherited from the build scripts:
- Every notebook is `coding_exercises/lecture-NN/exercise.ipynb`.
- Each step has the cell pattern: step-md → visible TODO code → "Click ▶" md →
  hidden-solution code (``metadata.jupyter.source_hidden = True``).
- A TODO stub is identified by either ``raise NotImplementedError`` or a
  literal ``# TODO`` marker in the visible code, combined with a
  hidden-solution code cell within the next 2 cells.
"""
from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass
from typing import Iterable

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
NOTEBOOKS_GLOB = os.path.join(REPO_ROOT, "coding_exercises", "lecture-*", "exercise.ipynb")


def discover_notebooks() -> list[str]:
    """All 27 notebook paths, sorted by lecture number."""
    return sorted(glob.glob(NOTEBOOKS_GLOB))


def lecture_number(path: str) -> int:
    """Extract the integer lecture number from a notebook path."""
    parent = os.path.basename(os.path.dirname(path))  # e.g. 'lecture-15'
    return int(parent.split("-", 1)[1])


def load_notebook(path: str) -> dict:
    """Parse a .ipynb as JSON. Avoids the nbformat dep for the static suite."""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def cell_source(cell: dict) -> str:
    src = cell.get("source", "")
    return "".join(src) if isinstance(src, list) else src


def is_hidden_solution(cell: dict) -> bool:
    if cell.get("cell_type") != "code":
        return False
    return cell.get("metadata", {}).get("jupyter", {}).get("source_hidden") is True


def is_pip_install_cell(cell: dict) -> bool:
    """A cell whose every non-blank, non-comment line is a ``!pip`` / ``%pip``.

    Handles backslash continuations: a ``!pip install foo \\`` line plus its
    indented continuation is treated as a single pip statement.
    """
    if cell.get("cell_type") != "code":
        return False
    in_magic_cont = False
    saw_pip = False
    for raw in cell_source(cell).splitlines():
        ln = raw.strip()
        if not ln:
            in_magic_cont = False  # blank ends a continuation
            continue
        if in_magic_cont:
            in_magic_cont = raw.rstrip().endswith("\\")
            continue
        if ln.startswith("#"):
            continue
        if ln.startswith(("!pip ", "%pip ")):
            saw_pip = True
            in_magic_cont = raw.rstrip().endswith("\\")
            continue
        return False
    return saw_pip


def neutralize_pip_lines(src: str) -> str:
    """Comment out ``!pip``/``%pip`` lines so a cell doesn't shell out under test.

    Used for mixed cells that have a pip-install line plus real Python.
    """
    out = []
    for line in src.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("!pip ", "%pip ")):
            out.append("# (pip skipped in test) " + line)
        else:
            out.append(line)
    return "\n".join(out)


@dataclass
class CellTriage:
    """Per-cell verdict after the test filter pass."""
    keep: bool
    reason: str  # 'normal', 'hidden-solution', 'todo-stub-skipped', 'pip-skipped'


def triage_cell(cells: list[dict], i: int) -> CellTriage:
    """Decide whether a cell should be executed.

    Drop only pure pip-install cells (they shell out, slow the test, and risk
    version skew). Keep everything else — including TODO stubs.

    Why keep the TODO stubs? They typically define module-level constants
    (``N_READS = 2000`` etc.) along with function stubs that raise
    ``NotImplementedError``. The stubs only fail when *called*; subsequent
    hidden-solution cells redefine the same names, so by the time anything
    actually invokes them, the real implementation is in scope. Skipping the
    TODO cells loses the constants and the notebook fails to execute.
    """
    c = cells[i]
    if is_pip_install_cell(c):
        return CellTriage(False, "pip-skipped")
    return CellTriage(True, "hidden-solution" if is_hidden_solution(c) else "normal")


def filter_cells_for_execution(nb: dict) -> tuple[list[dict], dict[str, int]]:
    """Return (cells_to_execute, counts) with TODO stubs + pip-installs stripped.

    Mutates each kept cell's source to neutralize any inline ``!pip`` lines.
    """
    cells = nb.get("cells", [])
    out: list[dict] = []
    counts = {"normal": 0, "hidden-solution": 0, "todo-stub-skipped": 0, "pip-skipped": 0}
    for i, c in enumerate(cells):
        verdict = triage_cell(cells, i)
        counts[verdict.reason] = counts.get(verdict.reason, 0) + 1
        if not verdict.keep:
            continue
        if c.get("cell_type") == "code":
            new_src = neutralize_pip_lines(cell_source(c))
            c = dict(c)  # shallow copy so we don't mutate the original notebook
            c["source"] = new_src
        out.append(c)
    return out, counts


# Lectures that fetch zero external data (purely synthetic in-cell generation).
# These are the only notebooks we attempt to execute end-to-end in the synthetic
# test suite.
SYNTHETIC_LECTURES = (
    1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 18, 22, 23, 24, 25, 27,
)


def is_synthetic(path: str) -> bool:
    return lecture_number(path) in SYNTHETIC_LECTURES
