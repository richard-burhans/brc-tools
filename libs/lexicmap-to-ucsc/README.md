# lexicmap-to-ucsc

Convert [LexicMap](https://github.com/shenwei356/LexicMap) search results (TSV) to
[UCSC Genome Browser](https://genome.ucsc.edu/) custom tracks (BED + bedGraph).

Accepts one or more input TSV files (e.g. one per chromosome) and produces a single
set of output files:

- **`.bed`** — BED6 with one interval per HSP (detailed hit view)
- **`.bedGraph`** — Hit density per genomic window (count of overlapping HSPs)
- **`.track.txt`** — UCSC custom track stanza for easy uploading

## Install

```bash
pip install lexicmap-to-ucsc
```

Or from source:

```bash
cd libs/lexicmap-to-ucsc
pip install -e ".[test]"
```

## Usage

### Single input file

```bash
lexicmap-to-ucsc lexicmap_output.tsv my_hits --window 200 --min-bitscore 200
```

### Multiple input files (e.g. one per chromosome)

```bash
lexicmap-to-ucsc chr1.tsv chr2.tsv chr3.tsv pf8_hits --window 200
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `inputs` (positional) | — | One or more LexicMap output TSV files |
| `prefix` (positional) | — | Output file prefix |
| `--window` | 1000 | Window size in bp for density track |
| `--min-bitscore` | 0 | Minimum bitscore to include an HSP |
| `--query` | None | Override chromosome name in BED (single input only) |

## Python API

```python
from lexicmap_to_ucsc import convert

result = convert(
    input_paths=["chr1.tsv", "chr2.tsv"],
    output_prefix="pf8_hits",
    window=200,
    min_bitscore=200,
)
print(f"{result['n_hsps']} HSPs across {result['n_chromosomes']} chromosome(s)")
```

## Coordinate conventions

- LexicMap coordinates are 1-based inclusive; BED is 0-based half-open.
- The density track counts the number of HSPs overlapping each window.
- HSPs on the minus strand (`sstr == '-'`) are encoded as strand `-` in BED;
  query coordinates are always reported on the forward strand.

## Run tests

```bash
cd libs/lexicmap-to-ucsc
pytest
```
