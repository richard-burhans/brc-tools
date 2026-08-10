"""Core logic for converting LexicMap TSV to UCSC BED + bedGraph tracks.

Coordinate conventions:
  - LexicMap coordinates are 1-based inclusive.
  - BED is 0-based half-open.
  - The density track counts the number of HSPs overlapping each window.
  - HSPs on the minus strand (sstr == '-') are encoded as strand '-' in BED;
    query coordinates are always reported on the forward strand.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass, field
from typing import Iterator


@dataclass
class HspRow:
    """One row from a LexicMap search TSV."""
    query: str
    qlen: int
    qstart: int
    qend: int
    sgenome: str
    sseqid: str
    hsp: str
    bitscore: float
    sstr: str


@dataclass
class ChromInfo:
    """Per-chromosome metadata accumulated during streaming."""
    name: str
    qlen: int
    intervals: list[tuple[int, int]] = field(default_factory=list)
    n_hsps: int = 0


def read_lexicmap(path: str) -> Iterator[dict]:
    """Read a LexicMap TSV file, yielding one dict per data row.

    Skips duplicate header lines that LexicMap sometimes emits.
    """
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("query") == "query":
                continue
            yield row


def _parse_row(row: dict) -> HspRow:
    """Convert a raw dict row into an HspRow with typed fields."""
    return HspRow(
        query=row["query"],
        qlen=int(row["qlen"]),
        qstart=int(row["qstart"]),
        qend=int(row["qend"]),
        sgenome=row["sgenome"],
        sseqid=row["sseqid"],
        hsp=row["hsp"],
        bitscore=float(row["bitscore"]),
        sstr=row.get("sstr", "+"),
    )


def _bed_fields(hsp: HspRow, query_name: str) -> list:
    """Convert an HspRow to a BED6 row [chrom, start, end, name, score, strand]."""
    start = min(hsp.qstart, hsp.qend) - 1
    end = max(hsp.qstart, hsp.qend)
    name = f"{hsp.sgenome}|{hsp.sseqid}|hsp{hsp.hsp}"
    score = min(int(hsp.bitscore), 1000)
    strand = "-" if hsp.sstr == "-" else "+"
    return [query_name, start, end, name, score, strand]


def compute_density_from_intervals(
    intervals: list[tuple[int, int]],
    query_name: str,
    qlen: int,
    window: int,
) -> list[tuple]:
    """Compute per-window HSP overlap count, returning bedGraph rows.

    Args:
        intervals: list of (start, end) tuples, 1-based inclusive.
        query_name: chromosome name for the bedGraph chrom column.
        qlen: length of the query sequence in bp.
        window: window size in bp.

    Returns:
        List of (chrom, start, end, value) tuples (0-based half-open).
    """
    if not intervals:
        return []

    n_windows = (qlen + window - 1) // window

    events: list[tuple[int, int]] = []
    for s, e in intervals:
        events.append((s, 1))
        events.append((e + 1, -1))
    events.sort()

    result: list[tuple] = []
    ei = 0
    active = 0
    n_events = len(events)

    for wi in range(n_windows):
        w_start = wi * window + 1  # 1-based
        w_end = min((wi + 1) * window, qlen)

        # Process all events up to w_start
        while ei < n_events and events[ei][0] <= w_start:
            active += events[ei][1]
            ei += 1

        # Track max active count within this window
        max_active = active
        while ei < n_events and events[ei][0] <= w_end:
            active += events[ei][1]
            if active > max_active:
                max_active = active
            ei += 1

        if max_active > 0:
            result.append((query_name, w_start - 1, w_end, max_active))

    return result


def write_track_stanza(
    prefix: str,
    query_names: list[str],
    window: int,
    n_hsps: int,
    min_bitscore: float,
    bed_exists: bool,
    bedgraph_exists: bool,
) -> str:
    """Generate a UCSC custom track stanza string."""
    qdesc = ", ".join(query_names) if query_names else "unknown"

    filter_desc = ""
    if min_bitscore > 0:
        filter_desc = f", bitscore>={min_bitscore}"

    lines: list[str] = []

    if bedgraph_exists:
        bg_name = f"{prefix}.bedGraph"
        lines.append(
            f'track type=bedGraph name="{prefix}_density" '
            f'description="Lexicmap HSP density (window={window}bp){filter_desc} on {qdesc}" '
            f"visibility=full color=200,0,0 altColor=100,0,0 priority=1"
        )
        lines.append(f"# Upload {bg_name} as a file, or host it and provide a URL below:")
        lines.append(f"# bigDataUrl=https://your-server/{bg_name}")
        lines.append("")

    if bed_exists:
        bed_name = f"{prefix}.bed"
        lines.append(
            f'track name="{prefix}_hsps" '
            f'description="Lexicmap HSP hits ({n_hsps} HSPs{filter_desc}) on {qdesc}" '
            f"visibility=pack itemRgb=On priority=2"
        )
        lines.append(f"# Upload {bed_name} as a file, or host it and provide a URL below:")
        lines.append(f"# bigDataUrl=https://your-server/{bed_name}")
        lines.append("")

    return "\n".join(lines)


def convert(
    input_paths: list[str],
    output_prefix: str,
    window: int = 1000,
    min_bitscore: float = 0.0,
    query_override: str | None = None,
) -> dict:
    """Convert one or more LexicMap TSV files to UCSC track files.

    Streams each input file once, writing BED lines while collecting
    per-chromosome intervals for density computation. Produces:
      - {prefix}.bed
      - {prefix}.bedGraph
      - {prefix}.track.txt

    Args:
        input_paths: list of LexicMap TSV file paths (one or more).
        output_prefix: prefix for output files.
        window: window size in bp for density track.
        min_bitscore: minimum bitscore to include an HSP.
        query_override: override chromosome name in BED (single-file only).

    Returns:
        Dict with summary stats: n_hsps, n_filtered, chromosomes, bed_path, etc.
    """
    bed_path = f"{output_prefix}.bed"
    bg_path = f"{output_prefix}.bedGraph"
    track_path = f"{output_prefix}.track.txt"

    chromosomes: dict[str, ChromInfo] = {}
    total_hsps = 0
    total_filtered = 0

    with open(bed_path, "w", newline="") as bed_f:
        bed_writer = csv.writer(bed_f, delimiter="\t")

        for input_path in input_paths:
            for raw_row in read_lexicmap(input_path):
                hsp = _parse_row(raw_row)

                if hsp.bitscore < min_bitscore:
                    total_filtered += 1
                    continue

                qname = query_override or hsp.query

                if qname not in chromosomes:
                    chromosomes[qname] = ChromInfo(name=qname, qlen=hsp.qlen)
                chrom = chromosomes[qname]

                bed_writer.writerow(_bed_fields(hsp, qname))

                lo = min(hsp.qstart, hsp.qend)
                hi = max(hsp.qstart, hsp.qend)
                chrom.intervals.append((lo, hi))
                chrom.n_hsps += 1
                total_hsps += 1

    if total_hsps == 0:
        sys.exit(f"No HSPs pass the filter ({total_filtered} filtered out).")

    # Write bedGraph with per-chromosome density
    total_windows = 0
    with open(bg_path, "w", newline="") as bg_f:
        bg_writer = csv.writer(bg_f, delimiter="\t")
        for qname in sorted(chromosomes.keys()):
            chrom = chromosomes[qname]
            density = compute_density_from_intervals(
                chrom.intervals, qname, chrom.qlen, window
            )
            for row in density:
                bg_writer.writerow(row)
            total_windows += len(density)

    # Write track stanza
    query_names = sorted(chromosomes.keys())
    stanza = write_track_stanza(
        output_prefix, query_names, window, total_hsps, min_bitscore,
        bed_exists=True, bedgraph_exists=True,
    )
    with open(track_path, "w") as f:
        f.write(stanza + "\n")

    return {
        "n_hsps": total_hsps,
        "n_filtered": total_filtered,
        "n_chromosomes": len(chromosomes),
        "n_windows": total_windows,
        "bed_path": bed_path,
        "bedgraph_path": bg_path,
        "track_path": track_path,
        "chromosomes": query_names,
    }
