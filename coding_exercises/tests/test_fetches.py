"""Probe every external data source the coding-exercise notebooks depend on.

Each notebook either generates synthetic data in-cell or fetches from a public
endpoint at runtime. This suite probes those endpoints so we know in advance
which notebooks will hit their fallback path on any given day.

A failing probe is informational, not catastrophic: every notebook ships a
deterministic synthetic fallback, so the exercise still runs offline. The
purpose of this suite is to catch endpoint drift (renamed URLs, removed Pfam
families, dead model IDs) before students do.

Usage
-----
    python coding_exercises/tests/test_fetches.py            # summary report
    python coding_exercises/tests/test_fetches.py --verbose  # full detail
    python coding_exercises/tests/test_fetches.py --json     # machine readable
    python -m pytest coding_exercises/tests/test_fetches.py  # pytest mode
"""
from __future__ import annotations

import argparse
import gzip
import json
import socket
import ssl
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Callable

DEFAULT_TIMEOUT = 20  # seconds — Ensembl + InterPro can be slow under load
USER_AGENT = "bioinformatics-course-test-fetches/1.0 (vojislav@nephronomics.com)"


def _build_ssl_context() -> ssl.SSLContext:
    """Use certifi's CA bundle if present.

    python.org's macOS installer ships without a wired-up system trust store,
    so plain ``urllib`` calls die with CERTIFICATE_VERIFY_FAILED on every
    HTTPS endpoint. certifi sidesteps this — and is already a transitive dep
    of nearly every notebook (via requests / huggingface-hub).
    """
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


_SSL_CTX = _build_ssl_context()


@dataclass
class Probe:
    lecture: str
    label: str
    url: str
    method: str = "HEAD"  # HEAD | GET | JSON
    expect_status: tuple[int, ...] = (200, 301, 302, 303, 307, 308, 405)
    expect_substring: str | None = None  # only for GET
    expect_json_key: str | None = None   # only for JSON
    headers: dict[str, str] = field(default_factory=dict)
    # 405 is allowed because some endpoints reject HEAD; we treat that as
    # "reachable" since GET from the notebook would still work.


# ─── Catalog ─────────────────────────────────────────────────────────────────
# One row per public endpoint touched by a notebook. The notebook number
# indicates which lecture's notebook would call this. Synthetic-only lectures
# (L01, L03-L06, L09-L14, L18, L22-L25, L27) are intentionally absent.

CATALOG: list[Probe] = [
    # L02 — NCBI EFetch (E. coli K-12 genome)
    Probe("L02", "NCBI EFetch ping",
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"),
    Probe("L02", "NCBI EFetch E. coli genome",
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id=U00096.3&rettype=fasta&retmode=text",
          method="GET", expect_substring=">"),

    # L07 — Scanpy pbmc3k
    Probe("L07", "Scanpy pbmc3k CDN",
          "https://falexwolf.de/data/pbmc3k_raw.h5ad"),

    # L08 — Scanpy paul15 + pbmc3k (shares L07 probe)
    Probe("L08", "Scanpy paul15 CDN",
          "https://falexwolf.de/data/paul15.h5"),

    # L15 — Pfam PF00042 seed + UniProt P69905 + AlphaFold-DB PAE
    Probe("L15", "InterPro Pfam PF00042 seed alignment",
          "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/"
          "?annotation=alignment:seed",
          method="GET", expect_substring="STOCKHOLM"),
    Probe("L15", "UniProt P69905 JSON",
          "https://rest.uniprot.org/uniprotkb/P69905.json",
          method="JSON", expect_json_key="primaryAccession"),
    # AlphaFold-DB rotates the version suffix on each release; v4 was retired
    # in 2026. The notebook fetches v6 first and falls back through v5/v4.
    Probe("L15", "AlphaFold-DB PAE JSON (v6)",
          "https://alphafold.ebi.ac.uk/files/"
          "AF-P69905-F1-predicted_aligned_error_v6.json",
          method="JSON"),

    # L16 — HuggingFace models + Ensembl + UniProt
    Probe("L16", "HF DNABERT-2 model metadata",
          "https://huggingface.co/api/models/zhihan1996/DNABERT-2-117M",
          method="JSON", expect_json_key="modelId"),
    Probe("L16", "HF ESM2 model metadata",
          "https://huggingface.co/api/models/facebook/esm2_t6_8M_UR50D",
          method="JSON", expect_json_key="modelId"),
    Probe("L16", "Ensembl REST ping",
          "https://rest.ensembl.org/info/ping?content-type=application/json",
          method="JSON"),
    Probe("L16", "UniProt P69905 (shared w/ L15)",
          "https://rest.uniprot.org/uniprotkb/P69905.json",
          method="JSON", expect_json_key="primaryAccession"),

    # L17 — ClinVar + gnomAD
    Probe("L17", "NCBI ClinVar esummary",
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
          "?db=clinvar&id=17661&retmode=json",
          method="JSON"),
    Probe("L17", "gnomAD GraphQL endpoint reachable",
          "https://gnomad.broadinstitute.org/api"),

    # L19 — UniProt BLOSUM-style query
    Probe("L19", "UniProt search for insulin",
          "https://rest.uniprot.org/uniprotkb/search"
          "?query=accession:P01308&format=json&size=1",
          method="JSON", expect_json_key="results"),

    # L20 — UniProt globin set + InterPro Pfam seed
    Probe("L20", "UniProt globin family search",
          "https://rest.uniprot.org/uniprotkb/search"
          "?query=family:globin+AND+reviewed:true&format=json&size=8",
          method="JSON", expect_json_key="results"),
    Probe("L20", "InterPro Pfam PF00042 (shared w/ L15)",
          "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/"
          "?annotation=alignment:seed",
          method="GET", expect_substring="STOCKHOLM"),

    # L21 — Ensembl BRCA1 region + InterPro Pfam HMM
    Probe("L21", "Ensembl human BRCA1 region",
          "https://rest.ensembl.org/sequence/region/human/"
          "17:43044295..43046294?content-type=application/json",
          method="JSON", expect_json_key="seq"),
    Probe("L21", "InterPro Pfam PF00042 HMM",
          "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00042/"
          "?annotation=hmm",
          method="GET", expect_substring="HMMER"),
    Probe("L21", "UniProt kinase search",
          "https://rest.uniprot.org/uniprotkb/search"
          "?query=family:kinase+AND+reviewed:true&format=json&size=5",
          method="JSON", expect_json_key="results"),

    # L26 — ChEMBL approved-drug query
    Probe("L26", "ChEMBL approved-drug molecule search",
          "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
          "?molecule_properties__mw_freebase__range=150,500"
          "&max_phase=4&limit=5",
          method="JSON", expect_json_key="molecules"),
]


# ─── Probe runner ────────────────────────────────────────────────────────────

@dataclass
class Result:
    probe: Probe
    ok: bool
    status: int | None
    elapsed_ms: float
    message: str  # empty if ok, otherwise reason


def _http(probe: Probe, timeout: float) -> Result:
    started = time.monotonic()
    headers = {"User-Agent": USER_AGENT, **probe.headers}
    method = "GET" if probe.method in ("GET", "JSON") else "HEAD"
    req = urllib.request.Request(probe.url, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as r:
            status = r.status
            body = r.read() if probe.method in ("GET", "JSON") else b""
            # InterPro serves Pfam alignments with Content-Encoding: gzip
            # (transport-level) and HMMs with Content-Type: application/gzip
            # (payload-level). urllib doesn't auto-decode either; do it here.
            ctype = (r.headers.get("Content-Type") or "").lower()
            cenc = (r.headers.get("Content-Encoding") or "").lower()
            if body and (cenc == "gzip" or "gzip" in ctype or body[:2] == b"\x1f\x8b"):
                try:
                    body = gzip.decompress(body)
                except OSError:
                    pass  # body wasn't actually gzip; leave as-is
    except urllib.error.HTTPError as e:
        elapsed = (time.monotonic() - started) * 1000
        ok = e.code in probe.expect_status
        msg = "" if ok else f"HTTP {e.code} {e.reason}"
        return Result(probe, ok, e.code, elapsed, msg)
    except (urllib.error.URLError, socket.timeout, TimeoutError) as e:
        elapsed = (time.monotonic() - started) * 1000
        return Result(probe, False, None, elapsed, f"network: {e}")
    except Exception as e:  # noqa: BLE001 — surface anything unexpected
        elapsed = (time.monotonic() - started) * 1000
        return Result(probe, False, None, elapsed, f"unexpected: {type(e).__name__}: {e}")

    elapsed = (time.monotonic() - started) * 1000
    if status not in probe.expect_status:
        return Result(probe, False, status, elapsed,
                      f"status {status} not in {probe.expect_status}")

    if probe.method == "GET" and probe.expect_substring:
        text = body.decode("utf-8", errors="replace")
        if probe.expect_substring not in text:
            return Result(probe, False, status, elapsed,
                          f"body missing substring {probe.expect_substring!r}")
    if probe.method == "JSON":
        try:
            payload = json.loads(body.decode("utf-8", errors="replace"))
        except json.JSONDecodeError as e:
            return Result(probe, False, status, elapsed, f"invalid JSON: {e}")
        if probe.expect_json_key and probe.expect_json_key not in payload:
            return Result(probe, False, status, elapsed,
                          f"JSON missing key {probe.expect_json_key!r}")
    return Result(probe, True, status, elapsed, "")


def run_all(catalog: list[Probe] = CATALOG, timeout: float = DEFAULT_TIMEOUT,
            on_result: Callable[[Result], None] | None = None) -> list[Result]:
    results: list[Result] = []
    for probe in catalog:
        r = _http(probe, timeout)
        results.append(r)
        if on_result is not None:
            on_result(r)
    return results


# ─── pytest hooks ────────────────────────────────────────────────────────────

def _pytest_id(p: Probe) -> str:
    return f"{p.lecture}::{p.label}"


def pytest_generate_tests(metafunc):  # noqa: D401
    if "probe" in metafunc.fixturenames:
        metafunc.parametrize("probe", CATALOG, ids=[_pytest_id(p) for p in CATALOG])


def test_endpoint(probe: Probe):
    """Each catalog entry becomes one pytest case."""
    r = _http(probe, DEFAULT_TIMEOUT)
    assert r.ok, f"{probe.lecture} {probe.label}: {r.message} ({probe.url})"


# ─── CLI driver ──────────────────────────────────────────────────────────────

def _print_human(r: Result, verbose: bool) -> None:
    glyph = "OK " if r.ok else "FAIL"
    line = f"  [{glyph}] {r.probe.lecture:<4} {r.probe.label:<48} {r.elapsed_ms:6.0f} ms"
    if not r.ok or verbose:
        line += f"  status={r.status}"
        if r.message:
            line += f"  ({r.message})"
        if verbose:
            line += f"\n        {r.probe.url}"
    print(line)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="show URL + status for every probe (not just failures)")
    ap.add_argument("--json", action="store_true",
                    help="emit machine-readable JSON instead of a human report")
    ap.add_argument("--lecture", "-l", action="append",
                    help="restrict to one lecture (e.g. -l L15 -l L20). May repeat.")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    args = ap.parse_args(argv)

    catalog = CATALOG
    if args.lecture:
        wanted = {x.upper() for x in args.lecture}
        catalog = [p for p in CATALOG if p.lecture in wanted]
        if not catalog:
            print(f"No probes match {args.lecture}", file=sys.stderr)
            return 2

    if args.json:
        results = run_all(catalog, timeout=args.timeout)
        out = [{
            "lecture": r.probe.lecture,
            "label": r.probe.label,
            "url": r.probe.url,
            "method": r.probe.method,
            "ok": r.ok,
            "status": r.status,
            "elapsed_ms": round(r.elapsed_ms, 1),
            "message": r.message,
        } for r in results]
        print(json.dumps(out, indent=2))
        return 0 if all(r.ok for r in results) else 1

    print(f"Probing {len(catalog)} endpoint(s) (timeout={args.timeout}s)\n")
    results = run_all(catalog, timeout=args.timeout,
                      on_result=lambda r: _print_human(r, args.verbose))
    n_ok = sum(1 for r in results if r.ok)
    n_fail = len(results) - n_ok
    print(f"\n{n_ok}/{len(results)} reachable, {n_fail} failed")
    if n_fail:
        print("\nFailed probes:")
        for r in results:
            if not r.ok:
                print(f"  - {r.probe.lecture} {r.probe.label}: {r.message}")
                print(f"      {r.probe.url}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
