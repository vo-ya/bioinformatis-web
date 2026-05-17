"""Static structure + syntax check for every notebook.

This suite makes no network calls and imports none of the notebooks' runtime
dependencies. It just parses each ``.ipynb`` and checks:

  1. Valid JSON + minimal notebook schema (nbformat 4, has cells, etc.).
  2. Every code cell ``compile()``s — catches syntax errors after a build-script
     edit before students see them.
  3. Cell structure looks right: at least one ``!pip install`` preamble, at
     least one self-check assert cell, every visible TODO has a corresponding
     hidden-solution cell within 2 cells.
  4. Each hidden-solution cell has both ``source_hidden`` and a non-empty body.

Usage::

    python coding_exercises/tests/test_structure.py
    python -m pytest coding_exercises/tests/test_structure.py -v

Runs in <1 second across all 27 notebooks; safe to wire into CI.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field

from _helpers import (  # type: ignore[import-not-found]
    cell_source,
    discover_notebooks,
    is_hidden_solution,
    is_pip_install_cell,
    lecture_number,
    load_notebook,
)


def _sanitize_magics(src: str) -> str:
    """Comment out IPython magics so ``compile()`` accepts the cell.

    Handles backslash line-continuations: a ``!pip install foo \`` line plus its
    continuation line(s) all get commented out as one logical magic call.
    """
    out: list[str] = []
    in_magic = False
    for line in src.splitlines():
        stripped = line.lstrip()
        if not in_magic and stripped.startswith(("!", "%")):
            out.append("# " + line)
            in_magic = line.rstrip().endswith("\\")
        elif in_magic:
            out.append("# " + line)
            in_magic = line.rstrip().endswith("\\")
        else:
            out.append(line)
    return "\n".join(out)


@dataclass
class StructureReport:
    path: str
    ok: bool
    n_cells: int
    n_code: int
    n_solutions: int
    n_todos_paired: int
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _has_self_check(cells: list[dict]) -> bool:
    """Self-check is any code cell with an ``assert`` plus the green-check string."""
    for c in cells:
        if c.get("cell_type") != "code":
            continue
        src = cell_source(c)
        if "assert" in src and ("Self-check" in src or "self-check" in src or "✅" in src):
            return True
    return False


def _todo_has_solution(cells: list[dict], i: int) -> bool:
    """Is there a hidden-solution code cell within 2 cells of ``i``?"""
    for j in range(i + 1, min(i + 3, len(cells))):
        if is_hidden_solution(cells[j]):
            return True
    return False


def check_notebook(path: str) -> StructureReport:
    rep = StructureReport(path=path, ok=True, n_cells=0, n_code=0,
                          n_solutions=0, n_todos_paired=0)
    try:
        nb = load_notebook(path)
    except Exception as e:  # noqa: BLE001
        rep.ok = False
        rep.errors.append(f"JSON parse failed: {e}")
        return rep

    if nb.get("nbformat") != 4:
        rep.errors.append(f"nbformat is {nb.get('nbformat')!r}, expected 4")
        rep.ok = False

    cells = nb.get("cells")
    if not isinstance(cells, list) or not cells:
        rep.errors.append("notebook has no cells")
        rep.ok = False
        return rep

    rep.n_cells = len(cells)
    rep.n_code = sum(1 for c in cells if c.get("cell_type") == "code")
    rep.n_solutions = sum(1 for c in cells if is_hidden_solution(c))

    # 1. compile() every code cell
    for idx, c in enumerate(cells):
        if c.get("cell_type") != "code":
            continue
        src = cell_source(c)
        sanitized = _sanitize_magics(src)
        try:
            compile(sanitized, f"{path}::cell[{idx}]", "exec")
        except SyntaxError as e:
            rep.ok = False
            rep.errors.append(f"cell[{idx}] SyntaxError: {e.msg} at line {e.lineno}")

    # 2. TODO → hidden-solution pairing
    todo_cells_found = 0
    unpaired_todos: list[int] = []
    for i, c in enumerate(cells):
        if c.get("cell_type") != "code" or is_hidden_solution(c):
            continue
        src = cell_source(c)
        if "raise NotImplementedError" in src or "# TODO" in src:
            todo_cells_found += 1
            if _todo_has_solution(cells, i):
                rep.n_todos_paired += 1
            else:
                unpaired_todos.append(i)
    if unpaired_todos:
        rep.errors.append(f"{len(unpaired_todos)} TODO cell(s) without hidden "
                          f"solution within 2 cells: indices {unpaired_todos}")
        rep.ok = False

    # 3. Soft expectations — warnings, not failures.
    if not any(is_pip_install_cell(c) for c in cells):
        rep.warnings.append("no !pip install preamble cell found")
    if not _has_self_check(cells):
        rep.warnings.append("no self-check assert cell found")
    if rep.n_solutions < 3:
        rep.warnings.append(f"only {rep.n_solutions} hidden-solution cells "
                            f"(expected ≥3 — usually 5 per the spec)")

    # 4. Each hidden-solution cell has a non-empty body.
    for i, c in enumerate(cells):
        if not is_hidden_solution(c):
            continue
        if not cell_source(c).strip():
            rep.ok = False
            rep.errors.append(f"cell[{i}] is a hidden-solution cell with empty body")

    return rep


def run_all() -> list[StructureReport]:
    return [check_notebook(p) for p in discover_notebooks()]


# ─── pytest hook ─────────────────────────────────────────────────────────────

def pytest_generate_tests(metafunc):  # noqa: D401
    if "notebook_path" in metafunc.fixturenames:
        paths = discover_notebooks()
        ids = [f"L{lecture_number(p):02d}" for p in paths]
        metafunc.parametrize("notebook_path", paths, ids=ids)


def test_structure(notebook_path: str):
    rep = check_notebook(notebook_path)
    assert rep.ok, "; ".join(rep.errors) or "(no error message)"


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="print per-notebook detail even on success")
    ap.add_argument("--warnings-as-errors", action="store_true",
                    help="treat warnings as failures")
    args = ap.parse_args(argv)

    reports = run_all()
    n_ok = 0
    n_fail = 0
    n_warn = 0
    print(f"Checking {len(reports)} notebook(s)\n")
    for rep in reports:
        lec = lecture_number(rep.path)
        tag = "OK " if rep.ok and not (args.warnings_as_errors and rep.warnings) else "FAIL"
        if tag == "OK ":
            n_ok += 1
        else:
            n_fail += 1
        line = (f"  [{tag}] L{lec:02d}  cells={rep.n_cells:>3}  "
                f"code={rep.n_code:>2}  solutions={rep.n_solutions}  "
                f"TODOs-paired={rep.n_todos_paired}")
        if rep.warnings:
            n_warn += len(rep.warnings)
        print(line)
        if rep.errors:
            for e in rep.errors:
                print(f"        ERROR: {e}")
        if rep.warnings and (args.verbose or args.warnings_as_errors):
            for w in rep.warnings:
                print(f"        warn:  {w}")

    print(f"\n{n_ok}/{len(reports)} ok, {n_fail} failed, {n_warn} warnings")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
