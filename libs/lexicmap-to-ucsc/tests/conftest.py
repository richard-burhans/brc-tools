"""Test fixtures: small LexicMap TSV files for testing."""

import os
import tempfile
import textwrap


HEADER = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore"


def make_tsv(rows: list[str]) -> str:
    """Create a temporary TSV file with the given data rows and return its path."""
    fd, path = tempfile.mkstemp(suffix=".tsv")
    os.close(fd)
    with open(path, "w") as f:
        f.write(HEADER + "\n")
        for row in rows:
            f.write(row + "\n")
    return path


def cleanup(path: str) -> None:
    if os.path.exists(path):
        os.remove(path)


def make_chr1_tsv() -> str:
    """A small TSV with 3 HSPs on chr1 (1000 bp)."""
    rows = [
        # query, qlen, hits, sgenome, sseqid, qcovGnm, cls, hsp, qcovHSP, alenHSP, pident, gaps, qstart, qend, sstart, send, sstr, slen, evalue, bitscore
        "chr1\t1000\t5\tSRR001\tSEQ_A\t100\t1\t1\t100\t200\t98.5\t0\t100\t300\t50\t250\t+\t500\t1e-50\t500.0",
        "chr1\t1000\t5\tSRR002\tSEQ_B\t100\t1\t2\t100\t150\t95.0\t2\t500\t650\t10\t160\t-\t300\t1e-30\t300.0",
        "chr1\t1000\t5\tSRR003\tSEQ_C\t100\t1\t3\t100\t100\t90.0\t0\t800\t900\t1\t100\t+\t200\t1e-20\t150.0",
    ]
    return make_tsv(rows)


def make_chr2_tsv() -> str:
    """A small TSV with 2 HSPs on chr2 (2000 bp)."""
    rows = [
        "chr2\t2000\t3\tSRR004\tSEQ_D\t100\t1\t1\t100\t250\t99.0\t0\t100\t350\t1\t250\t+\t400\t1e-60\t800.0",
        "chr2\t2000\t3\tSRR005\tSEQ_E\t100\t1\t2\t100\t180\t92.0\t1\t1500\t1680\t20\t200\t-\t300\t1e-25\t250.0",
    ]
    return make_tsv(rows)


def make_low_bitscore_tsv() -> str:
    """A TSV with mixed bitscores for testing min_bitscore filter."""
    rows = [
        "chrX\t500\t2\tSRR001\tSEQ_A\t100\t1\t1\t100\t200\t98.5\t0\t100\t300\t50\t250\t+\t500\t1e-50\t500.0",
        "chrX\t500\t2\tSRR002\tSEQ_B\t100\t1\t2\t100\t150\t95.0\t2\t100\t250\t10\t160\t+\t300\t1e-10\t50.0",
        "chrX\t500\t2\tSRR003\tSEQ_C\t100\t1\t3\t100\t100\t90.0\t0\t400\t500\t1\t100\t+\t200\t1e-5\t25.0",
    ]
    return make_tsv(rows)
