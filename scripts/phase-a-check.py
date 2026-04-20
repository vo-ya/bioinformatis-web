#!/usr/bin/env python3
"""
Phase A exit check for a new-lecture spec bundle.

Usage: python3 scripts/phase-a-check.py markdown_resources/lessonN_md_files/

Runs the Phase A exit check described in markdown_resources/generate_new_lesson_flow.md.
Exits 0 if every check passes; exits 1 if any check fails (with a human-readable report).

Checks applied:
  1. Artifact count matches between lecture.md and artifacts-spec.md.
  2. Figure count matches between lecture.md and figures-spec.md.
  3. Artifact and figure numbering is contiguous 1..N in both files.
  4. Every embed-marker filename is well-formed (NN-kebab.{html,svg} under the
     expected lecture-NN/ directory) AND matches the spec-file's `File:` field.
  5. Every `Lecture anchor: §N.M` in the spec files points at a real subsection.
  6. All callout labels use one of the approved set; distribution falls within
     the per-lecture target ranges.
  7. Top-matter block is present and names the correct output HTML file.
  8. Top-matter Duration equals the sum of Part kickers + Wrap-up kicker.
"""

import re
import sys
from pathlib import Path

ALLOWED_CALLOUTS = {
    'Intuition box', 'EE framing', 'Historical pointer',
    'Discussion prompt', 'Warning box',
}
TARGETS = {
    'Intuition box': (3, 4),
    'EE framing': (5, 7),
    'Historical pointer': (1, 3),
    'Discussion prompt': (1, 3),
    'Warning box': (2, 5),
}
FILENAME_SHAPE = re.compile(
    r'^(artifacts|diagrams)/lecture-(\d{2})/(\d{2})-[a-z0-9]+(?:-[a-z0-9]+)*\.(html|svg)$'
)

def parse_duration(text):
    """Parse '3h 30min', '≈3h 30min', '30 min', etc. → integer minutes."""
    text = text.replace('≈', '').replace('~', '').strip()
    m = re.match(r'(?:(\d+)\s*h)?\s*(\d+)?\s*min', text)
    if not m: return None
    h = int(m.group(1) or 0)
    mm = int(m.group(2) or 0)
    return h * 60 + mm

def main(spec_dir):
    spec_dir = Path(spec_dir)
    if not spec_dir.is_dir():
        print(f"error: {spec_dir} is not a directory", file=sys.stderr)
        return 2

    m = re.match(r'lesson(\d+)_md_files', spec_dir.name)
    if not m:
        print(f"error: directory name must be lessonN_md_files/, got {spec_dir.name}", file=sys.stderr)
        return 2
    N = int(m.group(1))
    NN = f"{N:02d}"

    lec_path = spec_dir / f"lecture-{NN}.md"
    fig_path = spec_dir / "figures-spec.md"
    art_path = spec_dir / "artifacts-spec.md"
    for p in (lec_path, fig_path, art_path):
        if not p.is_file():
            print(f"error: missing {p}", file=sys.stderr)
            return 2

    L = lec_path.read_text()
    F = fig_path.read_text()
    A = art_path.read_text()

    checks = []  # list of (name, ok, detail)

    # ── 1. artifact count ──────────────────────────────────────────
    lec_arts = re.findall(r'\*\*EMBED — Artifact #(\d+):', L)
    spec_arts = re.findall(r'^## \d+\. Artifact #(\d+)', A, re.M)
    ok = len(lec_arts) == len(spec_arts) and len(lec_arts) > 0
    checks.append((
        '1. artifact count',
        ok,
        f"lecture markers={len(lec_arts)} spec sections={len(spec_arts)}"
    ))

    # ── 2. figure count ────────────────────────────────────────────
    lec_figs = re.findall(r'\*\*FIGURE — Figure #(\d+):', L)
    spec_figs = re.findall(r'^## Figure (\d+) —', F, re.M)
    ok = len(lec_figs) == len(spec_figs) and len(lec_figs) > 0
    checks.append((
        '2. figure count',
        ok,
        f"lecture markers={len(lec_figs)} spec sections={len(spec_figs)}"
    ))

    # ── 3. contiguous numbering ────────────────────────────────────
    def is_contig(nums):
        s = sorted(int(x) for x in nums)
        return s == list(range(1, len(s) + 1))
    a_ok = is_contig(lec_arts) and is_contig(spec_arts)
    f_ok = is_contig(lec_figs) and is_contig(spec_figs)
    checks.append((
        '3. contiguous numbering',
        a_ok and f_ok,
        f"artifacts {'ok' if a_ok else 'gap'} · figures {'ok' if f_ok else 'gap'}"
    ))

    # ── 4. filename alignment + well-formedness ────────────────────
    lec_art_files = dict(re.findall(
        r'\*\*EMBED — Artifact #(\d+): [^*]+\*\* → `([^`]+)`', L))
    spec_art_files = {m.group(1): m.group(2) for m in re.finditer(
        r'## \d+\. Artifact #(\d+)[^\n]+\n\n\*\*File\*\*: `([^`]+)`', A)}
    lec_fig_files = dict(re.findall(
        r'\*\*FIGURE — Figure #(\d+): [^*]+\*\* → `([^`]+)`', L))
    spec_fig_files = {m.group(1): m.group(2) for m in re.finditer(
        r'## Figure (\d+) —[^\n]+\n\n\*\*File\*\*: `([^`]+)`', F)}

    issues_4 = []
    for label, d_lec, d_spec, expected_dir in [
        ('Artifact', lec_art_files, spec_art_files, f'artifacts/lecture-{NN}'),
        ('Figure',   lec_fig_files, spec_fig_files, f'diagrams/lecture-{NN}'),
    ]:
        for k in sorted(set(d_lec) | set(d_spec), key=int):
            l = d_lec.get(k); s = d_spec.get(k)
            if l is None:
                issues_4.append(f"{label} #{k}: missing from lecture.md")
                continue
            if s is None:
                issues_4.append(f"{label} #{k}: missing from spec.md")
                continue
            if l != s:
                issues_4.append(f"{label} #{k}: lecture={l!r} spec={s!r}")
                continue
            shape = FILENAME_SHAPE.match(l)
            if not shape:
                issues_4.append(f"{label} #{k}: {l!r} does not match NN-kebab.{{html,svg}} under {expected_dir}/ (stray character or non-kebab name?)")
                continue
            if shape.group(1) + '/lecture-' + shape.group(2) != expected_dir:
                issues_4.append(f"{label} #{k}: {l!r} should be under {expected_dir}/")
                continue
            if int(shape.group(3)) != int(k):
                issues_4.append(f"{label} #{k}: filename NN prefix {shape.group(3)!r} does not match entry number {k}")
    ok = not issues_4
    checks.append((
        '4. filename alignment + shape',
        ok,
        ('all paths well-formed and match' if ok else '\n     ' + '\n     '.join(issues_4))
    ))

    # ── 5. anchor alignment ────────────────────────────────────────
    subs = set(re.findall(r'^### (\d+\.\d+) ', L, re.M))
    anchor_a = re.findall(r'\*\*Lecture anchor\*\*: §(\d+\.\d+)', A)
    anchor_f = re.findall(r'\*\*Lecture anchor\*\*: §(\d+\.\d+)', F)
    missing = [a for a in anchor_a + anchor_f if a not in subs]
    checks.append((
        '5. anchor alignment',
        not missing,
        f"all anchors resolve" if not missing else f"missing subsections: {missing}"
    ))

    # ── 6. callout labels + distribution ───────────────────────────
    labels = re.findall(r'^> \*\*([^*]+?)\*\*', L, re.M)
    # Strip subtitle suffixes: "EE framing — something" / "EE framing (dB)"
    def canon(l):
        l = l.strip()
        for k in ALLOWED_CALLOUTS:
            if l == k or l.startswith(k + ' —') or l.startswith(k + ' ('):
                return k
        return l
    label_canons = [canon(l) for l in labels if l.strip() not in {'Duration', 'Audience', 'File'}]
    bad = [l for l in label_canons if l not in ALLOWED_CALLOUTS]
    totals = {k: label_canons.count(k) for k in ALLOWED_CALLOUTS}
    out_of_range = {k: v for k, v in totals.items() if not (TARGETS[k][0] <= v <= TARGETS[k][1])}
    ok = not bad and not out_of_range
    detail_parts = [f"{k}: {v}" for k, v in totals.items()]
    if bad:
        detail_parts.append(f"invented labels: {set(bad)}")
    if out_of_range:
        detail_parts.append(
            f"out of target range: " +
            ", ".join(f"{k}={v} (target {TARGETS[k][0]}–{TARGETS[k][1]})"
                      for k, v in out_of_range.items())
        )
    checks.append(('6. callouts', ok, ' · '.join(detail_parts)))

    # ── 7. top-matter ──────────────────────────────────────────────
    tm = re.match(
        rf'^# Lecture {N}[^\n]*\n\n> \*\*Duration\*\*:\s*([^\n]+)\n'
        rf'> \*\*Audience\*\*:[^\n]+\n'
        rf'> \*\*File\*\*: [^`]*`lectures/lecture-{NN}\.html`', L)
    if not tm:
        checks.append(('7. top-matter', False,
            f"missing or malformed; expected blockquote with Duration/Audience/File → lectures/lecture-{NN}.html"))
        top_dur_min = None
    else:
        top_dur_min = parse_duration(tm.group(1))
        checks.append(('7. top-matter', True,
            f"Duration={tm.group(1).strip()} · File→lectures/lecture-{NN}.html"))

    # ── 8. duration sum consistency ───────────────────────────────
    part_kickers = re.findall(r'^## Part \d+ —[^\n(]+\(([^)]+)\)', L, re.M)
    wrap_kicker = re.search(r'^## Wrap-up\s*\(([^)]+)\)', L, re.M)
    part_mins = [parse_duration(p) for p in part_kickers]
    wrap_min = parse_duration(wrap_kicker.group(1)) if wrap_kicker else None
    if None in part_mins or wrap_min is None:
        checks.append(('8. duration sum', False,
            f"could not parse all Part/Wrap-up kickers (parts={part_mins} wrap={wrap_min})"))
    elif top_dur_min is None:
        checks.append(('8. duration sum', False, "top-matter Duration unparseable"))
    else:
        summed = sum(part_mins) + wrap_min
        ok = summed == top_dur_min
        checks.append(('8. duration sum', ok,
            f"top-matter={top_dur_min}min · parts+wrap={summed}min "
            f"(parts {'+'.join(str(x) for x in part_mins)} +{wrap_min} wrap)"))

    # ── report ─────────────────────────────────────────────────────
    print(f"Phase A exit check for {spec_dir}")
    print("─" * 60)
    all_pass = True
    for name, ok, detail in checks:
        marker = 'PASS' if ok else 'FAIL'
        if not ok: all_pass = False
        print(f"  {marker}  {name}")
        if detail:
            for line in detail.split('\n'):
                print(f"         {line}")
    print("─" * 60)
    print(f"{'ALL 8 CHECKS PASS' if all_pass else 'SOME FAILED — fix the spec files and re-run'}")
    return 0 if all_pass else 1

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
