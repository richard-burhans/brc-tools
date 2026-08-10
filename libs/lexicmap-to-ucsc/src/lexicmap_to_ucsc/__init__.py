"""Convert LexicMap search results (TSV) to UCSC Genome Browser tracks.

Public API:
    convert() — main function, handles one or more input TSV files
    compute_density_from_intervals() — per-window HSP density
    read_lexicmap() — streaming TSV reader
    HspRow, ChromInfo — data classes
"""

from .core import (
    convert,
    compute_density_from_intervals,
    read_lexicmap,
    write_track_stanza,
    HspRow,
    ChromInfo,
)

__all__ = [
    "convert",
    "compute_density_from_intervals",
    "read_lexicmap",
    "write_track_stanza",
    "HspRow",
    "ChromInfo",
]
