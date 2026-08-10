# Libraries

This directory contains reusable Python libraries for genomic data processing, each published independently to PyPI and Bioconda.

## `genome-io`

Lightweight library for parsing and processing genomic data formats (BED, GFF3, FASTA, MAF, orthology tables, etc.).

- **PyPI**: `genome-io`
- **Bioconda**: `genome-io`
- **Status**: Stable (v0.1.0+)
- **Dependencies**: Minimal (pyfaidx optional)

See `genome-io/README.md` for details.

## `lexicmap-to-ucsc`

Convert LexicMap search results (TSV) to UCSC Genome Browser tracks (BED + bedGraph).
Accepts multiple input files (e.g. one per chromosome) and produces a single set of
output files. Pure Python, no dependencies beyond stdlib.

- **PyPI**: `lexicmap-to-ucsc`
- **Bioconda**: `lexicmap-to-ucsc`
- **Status**: Alpha (v0.1.0)
- **Dependencies**: None (Python stdlib only)

See `lexicmap-to-ucsc/README.md` for details.

## Future: `pangenome-helpers`

Orchestration logic for pangenome workflows (manifest loading, orthogroup filtering, CDS grouping, etc.). Will depend on `genome-io`.

- **PyPI**: `pangenome-helpers`
- **Bioconda**: `pangenome-helpers`
- **Status**: Planned
- **Dependencies**: `genome-io>=0.1.0`

---

## Development

Each library has its own `pyproject.toml`, tests, and documentation. To work on a specific library:

```bash
cd genome-io
pip install -e ".[test]"
pytest
```

To run all tests:

```bash
for lib in genome-io pangenome-helpers; do
  [ -d "$lib" ] && (cd "$lib" && pytest) || true
done
```
