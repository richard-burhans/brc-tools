"""Tests for the CLI interface."""

import subprocess
import sys
import os

from conftest import make_chr1_tsv, make_chr2_tsv, cleanup


class TestCLI:
    def test_single_file(self, tmp_path):
        tsv = make_chr1_tsv()
        prefix = str(tmp_path / "cli_out")
        try:
            result = subprocess.run(
                [sys.executable, "-m", "lexicmap_to_ucsc", tsv, prefix, "--window", "200"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0
            assert os.path.exists(f"{prefix}.bed")
            assert os.path.exists(f"{prefix}.bedGraph")
            assert os.path.exists(f"{prefix}.track.txt")
            assert "3 HSPs" in result.stderr
        finally:
            cleanup(tsv)

    def test_multiple_files(self, tmp_path):
        tsv1 = make_chr1_tsv()
        tsv2 = make_chr2_tsv()
        prefix = str(tmp_path / "cli_multi")
        try:
            result = subprocess.run(
                [sys.executable, "-m", "lexicmap_to_ucsc", tsv1, tsv2, prefix, "--window", "200"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0
            assert "5 HSPs" in result.stderr
            assert "2 chromosome(s)" in result.stderr
        finally:
            cleanup(tsv1)
            cleanup(tsv2)

    def test_min_bitscore(self, tmp_path):
        from conftest import make_low_bitscore_tsv
        tsv = make_low_bitscore_tsv()
        prefix = str(tmp_path / "cli_filter")
        try:
            result = subprocess.run(
                [sys.executable, "-m", "lexicmap_to_ucsc", tsv, prefix,
                 "--window", "200", "--min-bitscore", "100"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0
            assert "1 HSPs" in result.stderr
            assert "filtered" in result.stderr
        finally:
            cleanup(tsv)
