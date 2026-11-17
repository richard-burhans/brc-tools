# Test data provenance

Source: Cannabis sativa cs10, chromosome NC_044378.1 (71,238,074 bp), as NCBI distributes
it, with its NCBI soft-masking (lower case) kept. Local copy:
`data/kmer_poly_work/cs10_NC_044378.1.fa`.

Each record was cut with `samtools faidx` (1-based, inclusive coordinates). The only edit is
the header, renamed to `cs10_<start>_<end>`. The 50-column line width comes from samtools.

| File | Element id | Record | NC_044378.1 region | Lower-case bases |
|------|------------|--------|--------------------|------------------|
| `cs10_A.fa` | `cs10_A` | `cs10_12000001_12020000` | 12,000,001-12,020,000 | 25,440 of 40,000 (both records) |
| `cs10_A.fa` | `cs10_A` | `cs10_28000001_28020000` | 28,000,001-28,020,000 | |
| `cs10_B.fa` | `cs10_B` | `cs10_57000001_57020000` | 57,000,001-57,020,000 | 15,724 of 20,000 |

Neither file has N bases.

The regions were chosen by running the four detectors with the wrapper commands and the
pinned versions (blast 2.17.0, tantan 51, FasTAN 0.8, FastGA 1.5.20260729):

- `cs10_A`: all four detectors fire (dustmasker 59, WindowMasker 156, tantan 180,
  FasTAN 15 intervals).
- `cs10_B`: FasTAN finds no tandem array, so its BED is empty. This tests the
  empty-detector path through the union and the percent-masked table.

SHA-1:

```
651a21d6ba51eec5ee4487bd4bd501c21e344673  cs10_A.fa
c73c12db98edcc92b5f0e5cc09ac3de81edb2345  cs10_B.fa
```
