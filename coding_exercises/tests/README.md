# coding_exercises/tests

Three suites — fast to slow. All three pass under `pytest` if you'd rather run
them that way.

| File                          | Probes / cases | Runtime | What it catches                                                       |
|-------------------------------|----------------|---------|-----------------------------------------------------------------------|
| `test_structure.py`           | 27             | <1 s    | Cell-structure regressions, syntax errors, missing TODO ↔ solution pair |
| `test_fetches.py`             | 20             | ~15 s   | External-endpoint drift (URL moved, Pfam family gone, model retired)  |
| `test_execute_synthetic.py`   | 17             | ~70 s   | "The solution actually runs" — for the 17 pure-synthetic lectures     |

Run them all:

```bash
python -m pytest coding_exercises/tests/ -q
```

The 10 network-bound lectures (L02, L07–L08, L15–L17, L19–L21, L26) are
covered by `test_fetches.py` at the endpoint level; we don't execute them
end-to-end because their dependencies are heavy (scanpy, transformers, rdkit)
and the CI cost is not worth it. To run those manually, install the deps and
open the notebook in Jupyter.

---

## `test_structure.py` — static, no deps

For every notebook in `coding_exercises/lecture-*/exercise.ipynb`:

- valid JSON + nbformat 4
- every code cell `compile()`s (IPython `!` / `%` magics are stubbed before
  compilation, including backslash line-continuations)
- every TODO cell (one that contains `raise NotImplementedError` or `# TODO`)
  is followed by a hidden-solution cell within 2 cells
- every hidden-solution cell has a non-empty body
- soft warnings: missing `!pip install` preamble, missing self-check assert

```bash
python coding_exercises/tests/test_structure.py
python coding_exercises/tests/test_structure.py -v   # show warnings too
python -m pytest coding_exercises/tests/test_structure.py -v
```

## `test_fetches.py` — public-endpoint reachability

Probes every public URL the notebooks fetch from at runtime. A failing probe
is informational, not catastrophic: every notebook ships a deterministic
synthetic fallback. The point is to catch endpoint drift before students do.

```bash
python coding_exercises/tests/test_fetches.py
python coding_exercises/tests/test_fetches.py --verbose
python coding_exercises/tests/test_fetches.py -l L15 -l L20   # one lecture
python coding_exercises/tests/test_fetches.py --json
python -m pytest coding_exercises/tests/test_fetches.py -v
```

Coverage:

| Lecture | Endpoints                                                  |
|---------|------------------------------------------------------------|
| L02     | NCBI EFetch (E. coli K-12 genome)                          |
| L07     | Scanpy `pbmc3k` CDN                                        |
| L08     | Scanpy `paul15` + `pbmc3k` CDNs                            |
| L15     | InterPro Pfam seed, UniProt, AlphaFold-DB PAE              |
| L16     | HuggingFace (DNABERT-2 + ESM2), Ensembl, UniProt           |
| L17     | NCBI ClinVar esummary, gnomAD GraphQL                      |
| L19     | UniProt search                                             |
| L20     | UniProt globin family, InterPro Pfam seed                  |
| L21     | Ensembl region, InterPro Pfam HMM, UniProt kinase          |
| L26     | ChEMBL approved-drug search                                |

Adding a probe: append a `Probe(...)` row to `CATALOG` in `test_fetches.py`.
- `method="HEAD"` — just check reachability (fastest)
- `method="GET"` + `expect_substring=...` — check the body contains a marker
- `method="JSON"` + `expect_json_key=...` — parse JSON, assert a top-level key

A 405 on HEAD counts as reachable; some endpoints reject HEAD even when GET works.

## `test_execute_synthetic.py` — execute the synthetic lectures

End-to-end execution of the 17 lectures that need no network: **L01, L03–L06,
L09–L14, L18, L22–L25, L27**. Each notebook is executed in a fresh kernel via
`nbclient`.

```bash
python coding_exercises/tests/test_execute_synthetic.py
python coding_exercises/tests/test_execute_synthetic.py -l 22       # one lecture
python -m pytest coding_exercises/tests/test_execute_synthetic.py -v
```

Strategy:

- `!pip install` cells are stripped (the test env has the deps; reinstalling
  is slow + risks version skew).
- All other cells run, including the TODO stubs. A `NotImplementedError`
  raised inside a visible TODO cell is the *designed* failure mode and is
  silently allowed; any other exception in any cell fails the test.
- A hidden-solution cell raising any exception (including
  `NotImplementedError`) is a real bug and fails the test.

Dependencies: `nbclient`, `nbformat`, `ipykernel`, plus the science stack used
by the synthetic notebooks (`numpy scipy pandas matplotlib networkx
statsmodels scikit-learn biopython`).

```bash
pip install nbclient nbformat ipykernel
pip install numpy scipy pandas matplotlib networkx statsmodels scikit-learn biopython
```
