---
name: audit-lecture
description: Static audit of a lecture's figures (SVGs) and artifacts (HTML interactives) against the conventions established for this course. Runs read-only checks; reports findings as a punch list. Pass a lecture number (e.g. `/audit-lecture 13`), a range (`5-18`), or `all` to scan everything.
---

# Audit a lecture's figures and artifacts

Static, read-only audit. Does not modify files. Reports issues you should investigate and fix manually.

## Conventions enforced

### Figures (`diagrams/lecture-N/*.svg`)

- File parses as valid XML.
- Each SVG referenced from `lectures/lecture-N.html` exists on disk.
- Each `<img>` referencing the lecture's diagrams folder has a non-empty `alt` attribute.
- Each figure embedded as `<figure>...<figcaption>` rather than a bare `<img>`.
- File count in `diagrams/lecture-N/` matches the count claimed in the homepage card and the lecture's TOC / wrap-up.
- No mismatched-tag XML errors (we have hit `</title>` after `<text>` before).

### Artifacts (`artifacts/lecture-N/*.html`)

- Exactly **one** `<script src="../_shared/resize.js" defer></script>` near `</body>` (C6 smoke gate).
- The `<style>` / `<link>` block references `../_shared/artifact-theme.css`.
- Inline `<script>` passes `node --check` (catches the `time = '12–18 mo',` comma-operator bug we hit, plus typos and stray punctuation).
- An educational disclaimer is present (regex search for "Educational" / "Not for clinical use" / "Not a benchmark" / similar).
- An outcome banner / verdict pattern is present (regex search for `id="banner"` / `class="verdict"` / `class="banner"` / `class="outcome"`).
- The artifact is referenced from `lectures/lecture-N.html` as `<iframe ... src="../artifacts/lecture-N/<file>.html">`.

### Cross-references

- Lecture meta (e.g. `12 figures · 7 interactive tools`) matches actual file counts.
- Homepage card meta matches the lecture's own meta.

### Known bug patterns to grep for

These are concrete bugs we have found in past reviews. Surfacing them automatically:

- **Comma-operator assignment** in JS: `^\s*\w+ = ['"][^'"]*['"],\s*$` followed by another assignment on the next line. (Fixed in L18 artifact 6.)
- **Internal cap below slider max**: search `Math.min(\s*\d+\s*,\s*N\b)` patterns where the literal is far below realistic slider ranges. (Fixed in L13 artifact 7.)
- **Hardcoded numbers in lecture prose that should match the artifact** — flag pairs like `\$\d+` and `~\d+ GB` for manual cross-check against the corresponding artifact data. (Fixed in L1 cost table; L2 dense-rank arithmetic.)
- **Empty-input → NaN** patterns: search `Math.floor(.*length / K)` or `\.repeat\(.*length` without a length-guard above. (Fixed in L1 nanopore.)
- **Missing references**: a citation in prose (e.g. "FINEMAP", "CAVIAR") not present in the recommended-reading list. (Fixed in L13.)

The bug-pattern checks are best-effort heuristics; treat each match as "look here", not "this is a bug".

## Procedure

1. Resolve which lectures to audit:
   - `$1 == "all"` → 01..18
   - `$1` matches `\d+-\d+` → that range
   - `$1` is a single number → just that lecture
   - Default if no arg → ask user for a number.

   **Naming convention**: directories and files use zero-padded two-digit numbering everywhere (`lectures/lecture-05.html`, `diagrams/lecture-09/`, `artifacts/lecture-12/`). Always pass numbers through `printf "%02d"` before composing paths.

2. For each target lecture `N`:
   a. **SVG XML validity**: parse every `.svg` in `diagrams/lecture-N/`.
   b. **Reference completeness**: list every `<img src="../diagrams/lecture-N/..."` in `lectures/lecture-N.html`; confirm the referenced file exists; confirm `alt` is non-empty.
   c. **Reverse**: list every `.svg` on disk; confirm it's referenced in the lecture.
   d. **Artifact resize.js count**: each `.html` in `artifacts/lecture-N/` has exactly one `<script src="../_shared/resize.js"` reference.
   e. **Artifact JS syntax**: extract inline JS from `<script>...</script>` blocks, run `node --check`.
   f. **Artifact theme CSS**: each `.html` references `../_shared/artifact-theme.css`.
   g. **Disclaimer + banner**: each `.html` matches an educational-disclaimer regex AND an outcome-banner regex.
   h. **Iframe linkage**: each `.html` referenced from the lecture; each iframe in the lecture points at an existing artifact.
   i. **Bug-pattern grep**: run the regex patterns from the "Known bug patterns" section.
   j. **Counts cross-check**: lecture's stated `≈ 3h Xmin · NN figures · NN interactive tools` matches actual file counts; homepage card matches lecture meta.

3. Report per-lecture as a concise punch list. Group by **HIGH** (broken: XML invalid, JS fails to parse, broken reference, missing artifact) / **MEDIUM** (heuristic match suggesting a likely bug) / **LOW** (cosmetic: missing alt, missing disclaimer string, count mismatch by 1) / **OK** (everything passes).

4. Do NOT auto-fix. Leave fixes to a follow-up edit pass — a flagged item may not be a real bug, and the regexes are heuristic.

## Reference bash for each check

```bash
# Resolve repo root
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

# For each lecture N in the requested range:
for N in <range>; do
  LEC="lectures/lecture-${N}.html"
  DIAG="diagrams/lecture-${N}"
  ART="artifacts/lecture-${N}"

  # a. SVG XML validity
  for f in "$DIAG"/*.svg; do
    [[ -f "$f" ]] || continue
    python3 -c "import xml.etree.ElementTree as ET; ET.parse('$f')" 2>&1
  done

  # b. Image references
  grep -oE '\.\./diagrams/lecture-'"$N"'/[^"]+\.svg' "$LEC" | sort -u

  # c. Reverse: list SVGs on disk
  ls "$DIAG"/*.svg 2>/dev/null

  # d. resize.js count per artifact
  for f in "$ART"/*.html; do
    [[ -f "$f" ]] || continue
    n=$(grep -c 'resize\.js' "$f")
    [[ "$n" == "1" ]] || echo "FAIL $f resize.js count = $n"
  done

  # e. Artifact JS syntax
  for f in "$ART"/*.html; do
    [[ -f "$f" ]] || continue
    sed -n '/<script>/,/<\/script>/p' "$f" \
      | grep -v -E '^<script|^</script' \
      > "/tmp/audit_${N}_$(basename $f .html).js"
    node --check "/tmp/audit_${N}_$(basename $f .html).js"
  done

  # f. Theme CSS reference
  for f in "$ART"/*.html; do
    [[ -f "$f" ]] || continue
    grep -q 'artifact-theme\.css' "$f" || echo "FAIL $f missing theme css"
  done

  # g. Disclaimer + banner
  for f in "$ART"/*.html; do
    [[ -f "$f" ]] || continue
    grep -qiE 'educational|not for clinical|not a benchmark|illustrative only' "$f" || echo "MEDIUM $f no disclaimer"
    grep -qE 'class="banner|class="verdict|id="banner|class="outcome' "$f" || echo "MEDIUM $f no outcome banner"
  done

  # h. Iframe linkage
  grep -oE '\.\./artifacts/lecture-'"$N"'/[^"]+\.html' "$LEC" | sort -u

  # i. Bug-pattern grep — comma-operator assignment in JS
  for f in "$ART"/*.html; do
    [[ -f "$f" ]] || continue
    awk '/^\s*[A-Za-z_]\w* = ['\''\"][^'\''\"]*['\''\"],\s*$/ {print FILENAME":"NR":"$0}' "$f"
  done

  # j. Counts cross-check
  nfig=$(ls "$DIAG"/*.svg 2>/dev/null | wc -l)
  nart=$(ls "$ART"/*.html 2>/dev/null | wc -l)
  meta=$(grep -oE 'class="lecture-card-meta"[^<]*<' index.html | grep -A0 "lecture-${N}\.html" || true)
  echo "L${N}: ${nfig} figures, ${nart} artifacts on disk"
done
```

## Output format

For each lecture, emit a block like:

```
═══ Lecture 13 ═══
✓ 12 SVGs parse as valid XML
✓ 12 figures referenced; all have alt + figcaption
✓ 8 artifacts: resize.js × 1, JS parses, theme css present, disclaimer present
⚠ MEDIUM artifacts/lecture-13/06-burden-skat.html: comma-operator pattern at line 487
ℹ LOW   stats strip says "8 interactive tools" but disk has 7 (counted off-by-one)
```

End with a top-of-file summary: `N lectures audited · X HIGH · Y MEDIUM · Z LOW`.

## When NOT to use

- Do not use this skill to verify pedagogical correctness — only structural / mechanical checks.
- Do not use it to verify actual visual rendering — only XML / count / regex level.
- Do not auto-fix; report only.
