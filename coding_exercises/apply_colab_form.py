"""Apply Colab's collapsible-form pattern to every hidden-solution cell.

Why this exists
---------------
Our hidden-solution cells originally only had ``metadata.jupyter.source_hidden
= true``, which is the **JupyterLab** convention. Colab ignores that key and
renders the cell as plain visible code — so students would see the reference
solution immediately on opening the notebook.

Colab's actual hide mechanism is:
  1. ``metadata.cellView = "form"``
  2. The source's first non-blank line must be ``#@title <label>``.

When both are present, Colab renders the cell as a collapsed widget showing
only the title; clicking the title reveals the body. Reference:
https://colab.research.google.com/notebooks/forms.ipynb

This module
-----------
``apply_colab_form(nb)`` walks an in-memory notebook and patches every
hidden-solution cell to have both pieces. Idempotent — safe to re-run.

Running this file directly walks every ``coding_exercises/lecture-*/exercise.ipynb``
and rewrites them in place. Build scripts can also import the helper and call
it just before ``nbformat.write`` so the property survives regeneration.
"""
from __future__ import annotations

import glob
import os
import sys

TITLE_LINE = "#@title 🔓 Reference solution — click to reveal\n"


def _is_hidden_solution(cell) -> bool:
    if cell.get("cell_type") != "code":
        return False
    return cell.get("metadata", {}).get("jupyter", {}).get("source_hidden") is True


def _source_to_str(source) -> str:
    return source if isinstance(source, str) else "".join(source)


def apply_colab_form(nb) -> int:
    """Patch every hidden-solution cell in ``nb`` for Colab compatibility.

    Returns the number of cells modified. Accepts either an ``nbformat``
    notebook object (with attribute access) or a plain dict.
    """
    n_patched = 0
    cells = nb["cells"] if isinstance(nb, dict) else nb.cells
    for cell in cells:
        if not _is_hidden_solution(cell):
            continue
        # Ensure cellView: form.
        md = cell["metadata"] if isinstance(cell, dict) else cell.metadata
        if md.get("cellView") != "form":
            md["cellView"] = "form"
        # Ensure source starts with #@title.
        src = _source_to_str(cell["source"] if isinstance(cell, dict) else cell.source)
        if not src.lstrip().startswith("#@title"):
            new_src = TITLE_LINE + src
            if isinstance(cell, dict):
                cell["source"] = new_src
            else:
                cell.source = new_src
        n_patched += 1
    return n_patched


def _main(argv: list[str]) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    paths = sorted(glob.glob(os.path.join(here, "lecture-*", "exercise.ipynb")))
    if not paths:
        print("No notebooks found under coding_exercises/lecture-*/", file=sys.stderr)
        return 2

    import json
    total = 0
    for path in paths:
        with open(path, encoding="utf-8") as f:
            nb = json.load(f)
        n = apply_colab_form(nb)
        total += n
        with open(path, "w", encoding="utf-8") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
            f.write("\n")
        rel = os.path.relpath(path, here)
        print(f"  {rel}: patched {n} cell(s)")
    print(f"\nTotal cells patched: {total} across {len(paths)} notebook(s)")
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
