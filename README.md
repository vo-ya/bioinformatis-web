# bioinformatics-web

Lecture material and companion website for a bioinformatics course taught at a School of Electrical Engineering. Each lecture is a long scrollable HTML page with inline interactive artifacts and static figures — it's the same document the lecturer projects in class and the students read afterwards.

**Live site:** <https://vo-ya.github.io/bioinformatis-web/>

## Run it locally

No build step. Any static file server works:

```
cd bioinformatics-web
python3 -m http.server 8000
```

Open <http://localhost:8000/>. Alternatives: `npx serve .` if you have Node, or the VS Code "Live Server" extension. Opening files directly from disk (`file://`) mostly works but iframe auto-resize and some browser security policies are more reliable over `http://`.

## Repository layout

```
.
├── index.html                            # Landing page
├── lectures/
│   ├── lecture-01.html ... lecture-27.html   # Long-form lecture pages (27 shipped)
├── artifacts/
│   ├── _shared/                          # Shared design tokens + postMessage resizer
│   ├── lecture-01/ ... lecture-27/       # 6–8 interactive HTML artifacts per lecture
├── diagrams/
│   ├── lecture-01/ ... lecture-27/       # Static SVG figures (10–13 per lecture)
├── coding_exercises/                     # One 60-min Colab notebook per lecture
│   ├── lecture-01/ ... lecture-27/       # exercise.ipynb + build_notebook.py
│   ├── apply_colab_form.py               # Patcher: wires hidden solutions to Colab's #@title form pattern
│   └── tests/                            # Structural + endpoint + execution test suites
├── assets/
│   ├── styles.css                        # Site-wide tokens + homepage
│   ├── lecture.css                       # Lecture-page layout, callouts, figures
│   └── logo/                             # Mark, wordmark, favicon variants
├── scripts/
│   └── phase-a-check.py                  # Automated Phase A exit-check runner
├── markdown_resources/                   # Specs + source-of-truth content
│   ├── website-spec.md  homepage-spec.md  logo-spec.md         # Course-wide specs
│   ├── lecture-style-guide.md  diagram-style-guide.md
│   ├── generate_new_lesson_flow.md                              # Phase A/B/C recipe
│   └── lessonN_md_files/                                        # Per-lecture sources (lesson1 ... lesson27)
│       ├── artifacts-spec.md  figures-spec.md  lecture-NN.md
└── proof_html_resources/                 # Standalone-viewable design references
    ├── fastq-anatomy-proof.html
    └── logo-proof-v2.html
```

Each lecture embeds its artifacts via `<iframe>` from `artifacts/lecture-NN/` and its figures via `<img>` from `diagrams/lecture-NN/`. All paths are relative — the site runs under any subpath.

## Specs & authoring guides

These are the sources of truth when adding content. Read the relevant one before editing code.

| File | Use for |
|---|---|
| `markdown_resources/website-spec.md` | Site architecture, folder layout, design tokens, embedding contracts |
| `markdown_resources/homepage-spec.md` | Landing page layout and lecture-card list |
| `markdown_resources/logo-spec.md` | Course mark and wordmark |
| `markdown_resources/lecture-style-guide.md` | Voice, structure, callout types, formatting for lecture prose |
| `markdown_resources/diagram-style-guide.md` | House style for SVG diagrams (colour, stroke, typography, arrows) |
| `markdown_resources/generate_new_lesson_flow.md` | End-to-end recipe for generating a new lecture (phases, steps, validations) |
| `scripts/phase-a-check.py` | Automated runner for the Phase A exit check — run before every Phase B |
| `markdown_resources/lessonN_md_files/` | Per-lecture sources (lesson1 through lesson27) |
| `proof_html_resources/fastq-anatomy-proof.html` | Canonical SVG reference — copy arrow markers and typography patterns from here |
| `proof_html_resources/logo-proof-v2.html` | Logo visual reference; course mark is locked to Variant A3 |

## Design system in one glance

- **Light mode only.** No dark-mode scaffolding, no `prefers-color-scheme` queries. See `website-spec.md` §5.7 for why.
- **Three fonts.** Source Serif 4 for headings, Inter for body/UI, JetBrains Mono for any biological sequence or code.
- **One accent.** Deep cobalt (`--accent: #1e3a8a`). Used sparingly.
- **Base colours are locked.** `--base-a`/`--base-t`/`--base-g`/`--base-c`/`--base-u`/`--base-n` — students build visual memory across artifacts.
- **Zero build step.** Plain HTML, CSS, vanilla JS. Chart.js and KaTeX via CDN where needed.

## Coding exercises (Colab)

Every lecture has a companion 60-minute Jupyter notebook under
`coding_exercises/lecture-NN/exercise.ipynb`. The lecture page exposes them
via the "Open in Colab" badge in its meta strip — students click it, Colab
spins up a free CPU session pointing at the notebook on `main`, and they
hit File → Save a copy in Drive if they want to edit and keep their work.

Each notebook has the same 5-step structure: a markdown step prompt, a
visible `# TODO` scaffold, and a hidden reference solution (rendered in
Colab as a collapsible "🔓 Reference solution — click to reveal" form).
Data is either synthetic-in-cell (NumPy seeded, 17 lectures) or fetched
from public endpoints with a deterministic synthetic fallback when the
network call fails (10 lectures, covering NCBI / UniProt / Ensembl /
InterPro–Pfam / AlphaFold-DB / HuggingFace / Scanpy / ClinVar / gnomAD /
ChEMBL).

Re-generation is from `build_notebook.py` in each lecture folder — never
hand-edit the `.ipynb`. Three test suites live in `coding_exercises/tests/`:

- `test_structure.py`           — static cell + syntax check on all 27 (<1 s)
- `test_fetches.py`             — probes every public endpoint (~15 s)
- `test_execute_synthetic.py`   — actually runs the 17 synthetic notebooks end-to-end (~70 s)

All three pass under `python -m pytest coding_exercises/tests/`.

## Current status

**Shipped: all 27 lectures · 186 interactive tools · 320 figures · 27 Colab coding exercises · ~96 contact hours.**

| # | Lecture | Time | Figures | Tools |
|---|---|---|---|---|
| 01 | Foundations: From Cells to Sequences to FASTQ | 3h 35m | 12 | 6 |
| 02 | Read Alignment: From Brute Force to FM Index and Back | 3h 50m | 10 | 6 |
| 03 | DNA Sequence Assembly: From Reads to Genomes | 3h 25m | 11 | 6 |
| 04 | Variant Calling: From Aligned Reads to Called Differences | 3h 30m | 11 | 6 |
| 05 | Bulk RNA-seq: From Reads to Transcript Abundances | 3h 50m | 12 | 6 |
| 06 | Differential Expression and Count Statistics | 3h 30m | 11 | 6 |
| 07 | Single-Cell RNA-seq Fundamentals | 3h 30m | 11 | 6 |
| 08 | Advanced Single-Cell: Trajectories, Integration, Multi-Modal | 3h 30m | 10 | 6 |
| 09 | ChIP-seq, ATAC-seq, and Peak Calling | 3h 30m | 12 | 7 |
| 10 | Methylation, Hi-C, and 3D Genome Organisation | 3h 30m | 12 | 8 |
| 11 | Long Reads and the Pangenome | 3h 30m | 12 | 8 |
| 12 | Population Genetics Fundamentals | 3h 30m | 12 | 7 |
| 13 | GWAS and Statistical Genetics | 3h 30m | 12 | 7 |
| 14 | Data Engineering, File Formats, and Reproducibility | 3h 30m | 12 | 7 |
| 15 | Protein Structure Prediction in the AlphaFold Era | 3h 30m | 12 | 7 |
| 16 | ML in Genomics: Architectures, Pitfalls, Frontiers | 3h 55m | 13 | 8 |
| 17 | Clinical Genomics, Variant Interpretation, and Ethics | 4h 00m | 13 | 8 |
| 18 | Cancer Genomics: Integrated Capstone | 3h 30m | 12 | 7 |
| 19 | BLAST and Sequence Search Statistics | 3h 30m | 12 | 7 |
| 20 | Multiple Sequence Alignment, Phylogenetics, and Comparative Genomics | 3h 55m | 13 | 8 |
| 21 | HMMs, Profile HMMs, and Gene Finding | 3h 45m | 12 | 7 |
| 22 | Network Biology and Pathway Analysis | 3h 15m | 12 | 7 |
| 23 | Metagenomics and the Microbiome | 3h 30m | 12 | 7 |
| 24 | CRISPR Functional Screens and DepMap | 3h 20m | 12 | 7 |
| 25 | Causal Inference and Mendelian Randomisation | 3h 10m | 12 | 7 |
| 26 | Drug Discovery and Chemoinformatics | 3h 25m | 12 | 7 |
| 27 | Mass-Spectrometry Proteomics + Metabolomics | 3h 25m | 13 | 7 |

## Adding a new lecture

The full recipe is in `markdown_resources/generate_new_lesson_flow.md`. Short version:

1. **Phase A — specs.** Author `lecture-NN.md`, `figures-spec.md`, `artifacts-spec.md` under `markdown_resources/lessonN_md_files/`. Run `python3 scripts/phase-a-check.py markdown_resources/lessonN_md_files/` — it must exit 0 before Phase B starts.
2. **Phase B — build.** B2 lecture HTML · B3 SVG figures · B4 interactive artifacts (parallelise with background agents for 5+ artifacts) · B5 homepage card + stats strip. One commit per step.
3. **Phase C — quality.** C6 smoke test (local server, asset-link walk, JS `node --check`, SVG XML parse) · C7 review pass (within-part redundancy, "every section earns its seat", callouts, responsive, equations) · C8 push.

The `audit-lecture` skill (`.claude/skills/audit-lecture`) automates the structural checks: SVG XML validity, JS syntax, exactly-one resize.js per artifact, theme CSS reference, disclaimer + outcome banner per artifact, figure references resolving, lecture meta vs disk count parity.

## License

Apache 2.0 (course code and tooling). See `LICENSE` for the full text.

The course mark and wordmark carry their own identity restrictions — see `assets/logo/LICENSE.md`.
