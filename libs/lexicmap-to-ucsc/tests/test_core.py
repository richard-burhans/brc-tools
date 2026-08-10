"""Tests for lexicmap_to_ucsc.core."""

import os
import tempfile

import pytest

from lexicmap_to_ucsc import (
    convert,
    compute_density_from_intervals,
    read_lexicmap,
    HspRow,
    ChromInfo,
)
from lexicmap_to_ucsc.core import _parse_row, _bed_fields

from conftest import make_chr1_tsv, make_chr2_tsv, make_low_bitscore_tsv, cleanup


class TestReadLexicmap:
    def test_reads_rows(self):
        path = make_chr1_tsv()
        try:
            rows = list(read_lexicmap(path))
            assert len(rows) == 3
            assert rows[0]["query"] == "chr1"
            assert rows[0]["bitscore"] == "500.0"
        finally:
            cleanup(path)

    def test_skips_duplicate_header(self, tmp_path):
        path = str(tmp_path / "test.tsv")
        from conftest import HEADER
        with open(path, "w") as f:
            f.write(HEADER + "\n")
            f.write("chr1\t1000\t1\tSRR001\tSEQ_A\t100\t1\t1\t100\t200\t98.5\t0\t100\t300\t50\t250\t+\t500\t1e-50\t500.0\n")
            # Duplicate header
            f.write("query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tevalue\tbitscore\n")
            f.write("chr1\t1000\t1\tSRR002\tSEQ_B\t100\t1\t2\t100\t150\t95.0\t2\t500\t650\t10\t160\t-\t300\t1e-30\t300.0\n")
        rows = list(read_lexicmap(path))
        assert len(rows) == 2


class TestParseRow:
    def test_parses_correctly(self):
        raw = {
            "query": "chr1", "qlen": "1000", "qstart": "100", "qend": "300",
            "sgenome": "SRR001", "sseqid": "SEQ_A", "hsp": "1",
            "bitscore": "500.0", "sstr": "+",
        }
        hsp = _parse_row(raw)
        assert hsp.query == "chr1"
        assert hsp.qlen == 1000
        assert hsp.qstart == 100
        assert hsp.qend == 300
        assert hsp.bitscore == 500.0
        assert hsp.sstr == "+"


class TestBedFields:
    def test_forward_strand(self):
        hsp = HspRow("chr1", 1000, 100, 300, "SRR001", "SEQ_A", "1", 500.0, "+")
        fields = _bed_fields(hsp, "chr1")
        assert fields == ["chr1", 99, 300, "SRR001|SEQ_A|hsp1", 500, "+"]

    def test_reverse_strand(self):
        hsp = HspRow("chr1", 1000, 500, 650, "SRR002", "SEQ_B", "2", 300.0, "-")
        fields = _bed_fields(hsp, "chr1")
        assert fields[5] == "-"

    def test_score_capped_at_1000(self):
        hsp = HspRow("chr1", 1000, 100, 300, "SRR001", "SEQ_A", "1", 5000.0, "+")
        fields = _bed_fields(hsp, "chr1")
        assert fields[4] == 1000

    def test_swapped_coords(self):
        hsp = HspRow("chr1", 1000, 300, 100, "SRR001", "SEQ_A", "1", 500.0, "+")
        fields = _bed_fields(hsp, "chr1")
        assert fields[1] == 99
        assert fields[2] == 300


class TestComputeDensity:
    def test_empty_intervals(self):
        result = compute_density_from_intervals([], "chr1", 1000, 200)
        assert result == []

    def test_single_interval(self):
        result = compute_density_from_intervals([(100, 300)], "chr1", 1000, 200)
        # Window 0: 1-200, interval 100-300 overlaps → count 1
        # Window 1: 201-400, interval 100-300 overlaps → count 1
        assert len(result) == 2
        assert result[0] == ("chr1", 0, 200, 1)
        assert result[1] == ("chr1", 200, 400, 1)

    def test_overlapping_intervals(self):
        result = compute_density_from_intervals(
            [(100, 300), (150, 250)], "chr1", 1000, 200
        )
        # Both overlap window 0 (1-200) and window 1 (201-400)
        assert result[0] == ("chr1", 0, 200, 2)
        assert result[1] == ("chr1", 200, 400, 2)

    def test_non_overlapping_intervals(self):
        result = compute_density_from_intervals(
            [(100, 200), (800, 900)], "chr1", 1000, 200
        )
        # Window 0 (1-200): interval 100-200 overlaps
        # Window 3 (601-800): interval 800-900 starts at 800, overlaps
        # Window 4 (801-1000): interval 800-900 overlaps
        assert len(result) == 3
        assert result[0] == ("chr1", 0, 200, 1)
        assert result[1] == ("chr1", 600, 800, 1)
        assert result[2] == ("chr1", 800, 1000, 1)

    def test_window_at_end(self):
        # qlen=1000, window=300 → 4 windows: 1-300, 301-600, 601-900, 901-1000
        result = compute_density_from_intervals([(950, 1000)], "chr1", 1000, 300)
        assert len(result) == 1
        assert result[0] == ("chr1", 900, 1000, 1)


class TestConvert:
    def test_single_file(self, tmp_path):
        tsv = make_chr1_tsv()
        prefix = str(tmp_path / "out")
        try:
            result = convert([tsv], prefix, window=200)
            assert result["n_hsps"] == 3
            assert result["n_chromosomes"] == 1
            assert result["chromosomes"] == ["chr1"]
            assert os.path.exists(f"{prefix}.bed")
            assert os.path.exists(f"{prefix}.bedGraph")
            assert os.path.exists(f"{prefix}.track.txt")

            # Check BED content
            with open(f"{prefix}.bed") as f:
                lines = f.read().strip().split("\n")
            assert len(lines) == 3
            fields = lines[0].split("\t")
            assert fields[0] == "chr1"
            assert fields[1] == "99"  # 100-1 = 99 (0-based)
            assert fields[2] == "300"
            assert fields[5] == "+"  # forward strand

            # Check reverse strand in BED
            fields = lines[1].split("\t")
            assert fields[5] == "-"
        finally:
            cleanup(tsv)

    def test_multiple_files(self, tmp_path):
        tsv1 = make_chr1_tsv()
        tsv2 = make_chr2_tsv()
        prefix = str(tmp_path / "multi")
        try:
            result = convert([tsv1, tsv2], prefix, window=200)
            assert result["n_hsps"] == 5  # 3 + 2
            assert result["n_chromosomes"] == 2
            assert result["chromosomes"] == ["chr1", "chr2"]

            # BED should have 5 lines
            with open(f"{prefix}.bed") as f:
                lines = f.read().strip().split("\n")
            assert len(lines) == 5

            # bedGraph should have entries for both chromosomes
            with open(f"{prefix}.bedGraph") as f:
                bg_lines = f.read().strip().split("\n")
            chroms_in_bg = set(l.split("\t")[0] for l in bg_lines)
            assert "chr1" in chroms_in_bg
            assert "chr2" in chroms_in_bg
        finally:
            cleanup(tsv1)
            cleanup(tsv2)

    def test_min_bitscore_filter(self, tmp_path):
        tsv = make_low_bitscore_tsv()
        prefix = str(tmp_path / "filtered")
        try:
            result = convert([tsv], prefix, window=200, min_bitscore=100.0)
            assert result["n_hsps"] == 1  # only the 500.0 bitscore hit
            assert result["n_filtered"] == 2

            with open(f"{prefix}.bed") as f:
                lines = f.read().strip().split("\n")
            assert len(lines) == 1
        finally:
            cleanup(tsv)

    def test_query_override(self, tmp_path):
        tsv = make_chr1_tsv()
        prefix = str(tmp_path / "override")
        try:
            result = convert([tsv], prefix, window=200, query_override="custom_chr")
            assert result["chromosomes"] == ["custom_chr"]

            with open(f"{prefix}.bed") as f:
                first_line = f.readline()
            assert first_line.startswith("custom_chr\t")
        finally:
            cleanup(tsv)

    def test_no_hsps_exits(self, tmp_path):
        from conftest import HEADER
        path = str(tmp_path / "empty.tsv")
        with open(path, "w") as f:
            f.write(HEADER + "\n")
        with pytest.raises(SystemExit):
            convert([path], str(tmp_path / "empty_out"), window=200)

    def test_track_stanza_content(self, tmp_path):
        tsv = make_chr1_tsv()
        prefix = str(tmp_path / "stanza")
        try:
            convert([tsv], prefix, window=200)
            with open(f"{prefix}.track.txt") as f:
                content = f.read()
            assert "type=bedGraph" in content
            assert "track name=" in content
            assert "chr1" in content
        finally:
            cleanup(tsv)
