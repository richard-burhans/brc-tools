"""CLI entry point for lexicmap-to-ucsc."""

from __future__ import annotations

import argparse
import sys

from . import convert


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lexicmap-to-ucsc",
        description=(
            "Convert LexicMap search results (TSV) to UCSC Genome Browser tracks. "
            "Accepts one or more input TSV files (e.g. one per chromosome) and "
            "produces a single BED, bedGraph, and track stanza file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "inputs",
        nargs="+",
        help="one or more LexicMap output TSV files",
    )
    p.add_argument(
        "prefix",
        help="output file prefix (e.g. pf8_hits → pf8_hits.bed, .bedGraph, .track.txt)",
    )
    p.add_argument(
        "--window",
        type=int,
        default=1000,
        help="window size in bp for density track (default: 1000)",
    )
    p.add_argument(
        "--min-bitscore",
        type=float,
        default=0.0,
        help="minimum bitscore to include an HSP (default: 0, no filter)",
    )
    p.add_argument(
        "--query",
        default=None,
        help="override chromosome name in BED (single input only; ignored for multiple)",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    query_override = args.query if len(args.inputs) == 1 else None

    result = convert(
        input_paths=args.inputs,
        output_prefix=args.prefix,
        window=args.window,
        min_bitscore=args.min_bitscore,
        query_override=query_override,
    )

    print(f"Wrote BED: {result['bed_path']} ({result['n_hsps']} HSPs)", file=sys.stderr)
    if result["n_filtered"] > 0:
        print(f"  ({result['n_filtered']} filtered by bitscore)", file=sys.stderr)
    print(
        f"Wrote bedGraph: {result['bedgraph_path']} "
        f"({result['n_windows']} windows with coverage, "
        f"{result['n_chromosomes']} chromosome(s))",
        file=sys.stderr,
    )
    print(f"Wrote track stanza: {result['track_path']}", file=sys.stderr)
    print(f"\nDone! To view in UCSC:", file=sys.stderr)
    print(f"  1. Go to the UCSC Genome Browser for your assembly", file=sys.stderr)
    print(f"  2. click 'add custom tracks' / 'My Data → Custom Tracks'", file=sys.stderr)
    print(f"  3. Upload {result['bed_path']} and {result['bedgraph_path']}", file=sys.stderr)


if __name__ == "__main__":
    main()
