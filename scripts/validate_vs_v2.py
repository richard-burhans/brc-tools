#!/usr/bin/env python3
"""Validate Galaxy-produced pggb outputs against v2 native run.

Compares execution/outputs/ to Pv4-pangenome/v2/pggb_out/.
- GFA line counts (S/L/P) — sanity check.
- Canonical GFA md5 (sort S/L/P lines) — bonus stability check.
- odgi stats on .og files (built on-the-fly from GFA if missing) —
  expects node/edge counts within +-0.5%.

Writes report to execution/validation_report.md.
"""

from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GALAXY_OUT = REPO / "execution" / "outputs"
V2_OUT = Path("/home/anton/git/Pv4-pangenome/v2/pggb_out")
REPORT = REPO / "execution" / "validation_report.md"

ODGI_IMG = "quay.io/biocontainers/odgi:0.9.4--h077b44d_0"


def open_gfa(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def gfa_line_counts(path: Path) -> dict:
    counts = Counter()
    with open_gfa(path) as fh:
        for line in fh:
            counts[line[:1]] += 1
    return dict(counts)


def canonical_gfa_md5(path: Path) -> str:
    h_lines, s_lines, l_lines, p_lines = [], [], [], []
    with open_gfa(path) as fh:
        for line in fh:
            tag = line[:1]
            if tag == "H":
                h_lines.append(line.rstrip())
            elif tag == "S":
                s_lines.append(line.rstrip())
            elif tag == "L":
                l_lines.append(line.rstrip())
            elif tag == "P":
                p_lines.append(line.rstrip())
    blob = "\n".join(h_lines + sorted(s_lines) + sorted(l_lines) + sorted(p_lines))
    return hashlib.md5(blob.encode()).hexdigest()


def gfa_to_og_via_docker(gfa: Path, work: Path) -> Path:
    """Build .og from a (gzipped) GFA via odgi biocontainer."""
    gfa_local = work / "graph.gfa"
    if gfa.suffix == ".gz":
        with gzip.open(gfa, "rb") as src, gfa_local.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy(gfa, gfa_local)
    og_local = work / "graph.og"
    subprocess.check_call([
        "docker", "run", "--rm",
        "-v", f"{work}:/work", "-w", "/work",
        "-u", f"{__import__('os').getuid()}:{__import__('os').getgid()}",
        ODGI_IMG,
        "odgi", "build", "-g", gfa_local.name, "-o", og_local.name,
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return og_local


def odgi_stats(og: Path) -> dict:
    out = subprocess.check_output([
        "docker", "run", "--rm",
        "-v", f"{og.parent.resolve()}:/work:ro", "-w", "/work",
        ODGI_IMG,
        "odgi", "stats", "-i", og.name, "-S",
    ], text=True)
    lines = out.strip().splitlines()
    header, vals = lines[0], lines[1]
    # ⛔ strict=True, deliberately. This parses `odgi stats -S`, and a header row and a value row
    # of differing width is a malformed stat block, not something to silently truncate to the
    # shorter of the two -- which would yield a partial, plausible-looking stats dict.
    return dict(zip(header.split("\t"), vals.split("\t"), strict=True))


def main() -> int:
    galaxy_gfa = next(GALAXY_OUT.glob("*.smooth.final.gfa.gz"), None) \
                 or next(GALAXY_OUT.glob("*.smooth.fix.gfa.gz"), None) \
                 or next(GALAXY_OUT.glob("*.gfa.gz"), None)
    v2_gfa = next(V2_OUT.glob("*.smooth.fix.gfa.gz"), None)

    lines = ["# Galaxy vs v2-native pggb validation\n"]
    lines.append(f"- galaxy_gfa: `{galaxy_gfa}`")
    lines.append(f"- v2_gfa:     `{v2_gfa}`")
    if galaxy_gfa is None or v2_gfa is None:
        lines.append("\n**Missing inputs; cannot validate yet.**")
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text("\n".join(lines))
        print("\n".join(lines))
        return 1

    g_counts = gfa_line_counts(galaxy_gfa)
    v_counts = gfa_line_counts(v2_gfa)
    lines.append("\n## GFA line counts")
    lines.append("| tag | galaxy | v2 | diff %% |")
    lines.append("|----:|-------:|---:|-------:|")
    overall_ok = True
    for tag in sorted(set(g_counts) | set(v_counts)):
        gc = g_counts.get(tag, 0)
        vc = v_counts.get(tag, 0)
        pct = 100 * (gc - vc) / vc if vc else (float("inf") if gc else 0)
        ok = abs(pct) < 0.5
        lines.append(f"| {tag} | {gc} | {vc} | {pct:+.3f} |")
        if tag in ("S", "L", "P") and not ok:
            overall_ok = False
    lines.append(f"\n**S/L/P counts within +-0.5%%: {overall_ok}**")

    gmd5 = canonical_gfa_md5(galaxy_gfa)
    vmd5 = canonical_gfa_md5(v2_gfa)
    lines.append("\n## Canonical GFA md5")
    lines.append(f"- galaxy: `{gmd5}`")
    lines.append(f"- v2:     `{vmd5}`")
    lines.append(f"- match:  **{gmd5 == vmd5}**")

    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        g_work = tdir / "galaxy"
        g_work.mkdir()
        v_work = tdir / "v2"
        v_work.mkdir()
        try:
            g_og = gfa_to_og_via_docker(galaxy_gfa, g_work)
            v_og = gfa_to_og_via_docker(v2_gfa, v_work)
            gs = odgi_stats(g_og)
            vs = odgi_stats(v_og)
            lines.append("\n## odgi stats (built from GFA)")
            lines.append("| metric | galaxy | v2 | diff %% |")
            lines.append("|--------|-------:|---:|-------:|")
            stats_ok = True
            for k in gs:
                try:
                    gv = float(gs[k])
                    vv = float(vs[k])
                    pct = 100 * (gv - vv) / vv if vv else 0
                    if abs(pct) >= 0.5:
                        stats_ok = False
                    lines.append(f"| {k} | {gv:.0f} | {vv:.0f} | {pct:+.3f} |")
                except (ValueError, KeyError):
                    lines.append(f"| {k} | {gs.get(k)} | {vs.get(k)} | - |")
            lines.append(f"\n**odgi stats within +-0.5%%: {stats_ok}**")
        # ValueError belongs here beside CalledProcessError: `odgi_stats` zips strict, so a
        # malformed stat block raises from inside this block. Without it the exception escapes
        # past REPORT.write_text below and the run produces NO report at all -- discarding the
        # GFA counts and md5 sections already computed above, which do not depend on odgi.
        except (subprocess.CalledProcessError, ValueError) as e:
            lines.append(f"\n**odgi build/stats failed: {e}**")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
