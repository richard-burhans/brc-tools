# Validation — reproduction of the reference outputs

Every number below was produced by running the CLIs in `scripts/` on real data and
diffing against the reference outputs. Nothing here is estimated.

Interpreter: `/mnt/ssd/pv4_full/.driver-venv/bin/python` (pysam 0.24.0).

## 1. `part_haplotype_caller.py` vs `/home/anton/pyeuk-bench/lofreq_arm/calls_long.tsv`

Input: the marker-aligned filtered BAMs under `/home/anton/a1-cyclospora/results/<SPEC>/`.
Comparison: the caller-A rows of `calls_long.tsv`, matched on `(part, seq)`, with every
field (`span`, `reads`, `freq`, `name`, `match`, `nearest`, `edit`) required to agree.

Run on **all 203 specimens** of `specimens_203.txt`, not a subset:

```
specimens compared                 : 203
specimens identical (row set + all fields) : 203   (100%)
reference caller-A rows            : 4854
rows emitted by this CLI           : 4854
extra rows / missing rows          : 0 / 0
field disagreements                : 0
NOVEL rows reproduced (name = <PART>_NOV_<md5[:6]>) : 391
```

Per-specimen detail (specimen | reference rows | my rows | verdict):

| C_IA013_18 | 20 | 20 | identical |
| C_IA018_18 | 27 | 27 | identical |
| C_IA029_18 | 12 | 12 | identical |
| C_IA031_18 | 21 | 21 | identical |
| C_IA034_18 | 22 | 22 | identical |
| C_IA039_18 | 25 | 25 | identical |
| C_IA040_18 | 25 | 25 | identical |
| C_IA045_18 | 22 | 22 | identical |
| C_IA048_18 | 23 | 23 | identical |
| C_IA052_18 | 28 | 28 | identical |
| C_IA058_18 | 28 | 28 | identical |
| C_IA059_18 | 19 | 19 | identical |
| C_IA062_18 | 25 | 25 | identical |
| C_IA064_18 | 28 | 28 | identical |
| C_IA068_18 | 28 | 28 | identical |
| C_IA069_18 | 28 | 28 | identical |
| C_IA075_18 | 27 | 27 | identical |
| C_IA076_18 | 18 | 18 | identical |
| C_IA077_18 | 27 | 27 | identical |
| C_IA078_18 | 25 | 25 | identical |
| C_IA082_18 | 24 | 24 | identical |
| C_IA083_18 | 28 | 28 | identical |
| C_IA085_18 | 20 | 20 | identical |
| C_IA086_18 | 20 | 20 | identical |
| C_IA087_18 | 26 | 26 | identical |
| C_IA091_18 | 28 | 28 | identical |
| C_IA092_18 | 27 | 27 | identical |
| C_IA093_18 | 21 | 21 | identical |
| C_IA094_18 | 26 | 26 | identical |
| C_IA102_18 | 29 | 29 | identical |
| C_IA104_18 | 28 | 28 | identical |
| C_IA107_18 | 26 | 26 | identical |
| C_IA110_18 | 31 | 31 | identical |
| C_IA111_18 | 28 | 28 | identical |
| C_IA156_18 | 12 | 12 | identical |
| C_IL004_18 | 18 | 18 | identical |
| C_IL007_18 | 13 | 13 | identical |
| C_IL012_18 | 23 | 23 | identical |
| C_IL013_18 | 27 | 27 | identical |
| C_IL018_18 | 28 | 28 | identical |
| C_IL019_18 | 26 | 26 | identical |
| C_IL020_18 | 23 | 23 | identical |
| C_IL022_18 | 28 | 28 | identical |
| C_IL023_18 | 26 | 26 | identical |
| C_IL024_18 | 26 | 26 | identical |
| C_IL025_18 | 24 | 24 | identical |
| C_IL031_18 | 26 | 26 | identical |
| C_IL032_18 | 25 | 25 | identical |
| C_IL033_18 | 28 | 28 | identical |
| C_IL034_18 | 28 | 28 | identical |
| C_IL036_18 | 19 | 19 | identical |
| C_IL037_18 | 27 | 27 | identical |
| C_IL039_18 | 25 | 25 | identical |
| C_IL040_18 | 24 | 24 | identical |
| C_IL042_18 | 28 | 28 | identical |
| C_IL044_18 | 25 | 25 | identical |
| C_IL046_18 | 27 | 27 | identical |
| C_IL047_18 | 28 | 28 | identical |
| C_IL048_18 | 24 | 24 | identical |
| C_IL049_18 | 17 | 17 | identical |
| C_IL050_18 | 28 | 28 | identical |
| C_IL054_18 | 21 | 21 | identical |
| C_IL055_18 | 26 | 26 | identical |
| C_IL056_18 | 28 | 28 | identical |
| C_IL057_18 | 14 | 14 | identical |
| C_IL058_18 | 28 | 28 | identical |
| C_IL059_18 | 18 | 18 | identical |
| C_IL060_18 | 28 | 28 | identical |
| C_IL062_18 | 26 | 26 | identical |
| C_IL065_18 | 28 | 28 | identical |
| C_IL066_18 | 24 | 24 | identical |
| C_IL067_18 | 19 | 19 | identical |
| C_IL068_18 | 18 | 18 | identical |
| C_IL071_18 | 25 | 25 | identical |
| C_IL072_18 | 28 | 28 | identical |
| C_IL074_18 | 18 | 18 | identical |
| C_IL076_18 | 28 | 28 | identical |
| C_IL077_18 | 28 | 28 | identical |
| C_IL079_18 | 28 | 28 | identical |
| C_IL080_18 | 23 | 23 | identical |
| C_IL081_18 | 22 | 22 | identical |
| C_IL082_18 | 14 | 14 | identical |
| C_IL085_18 | 18 | 18 | identical |
| C_IL086_18 | 25 | 25 | identical |
| C_IL087_18 | 23 | 23 | identical |
| C_IL097_18 | 18 | 18 | identical |
| C_IL099_18 | 21 | 21 | identical |
| C_IL113_18 | 15 | 15 | identical |
| C_IL119_18 | 9 | 9 | identical |
| C_IL120_18 | 22 | 22 | identical |
| C_IL126_18 | 18 | 18 | identical |
| C_MO003_18 | 28 | 28 | identical |
| C_OH003_18 | 26 | 26 | identical |
| C_OH006_18 | 12 | 12 | identical |
| C_OH014_18 | 23 | 23 | identical |
| C_WI008_18 | 19 | 19 | identical |
| C_WI011_18 | 22 | 22 | identical |
| C_WI013_18 | 22 | 22 | identical |
| C_WI020_18 | 22 | 22 | identical |
| C_WI024_18 | 20 | 20 | identical |
| C_WI026_18 | 28 | 28 | identical |
| C_WI031_18 | 23 | 23 | identical |
| C_WI033_18 | 22 | 22 | identical |
| C_WI036_18 | 19 | 19 | identical |
| C_WI037_18 | 18 | 18 | identical |
| C_WI038_18 | 23 | 23 | identical |
| C_WI039_18 | 23 | 23 | identical |
| C_WI040_18 | 22 | 22 | identical |
| C_WI041_18 | 20 | 20 | identical |
| C_WI043_18 | 28 | 28 | identical |
| C_WI044_18 | 28 | 28 | identical |
| C_WI045_18 | 22 | 22 | identical |
| C_WI046_18 | 24 | 24 | identical |
| C_WI047_18 | 22 | 22 | identical |
| C_WI048_18 | 25 | 25 | identical |
| C_WI049_18 | 28 | 28 | identical |
| C_WI050_18 | 28 | 28 | identical |
| C_WI051_18 | 28 | 28 | identical |
| C_WI052_18 | 28 | 28 | identical |
| C_WI053_18 | 27 | 27 | identical |
| C_WI054_18 | 28 | 28 | identical |
| C_WI055_18 | 28 | 28 | identical |
| C_WI056_18 | 28 | 28 | identical |
| C_WI057_18 | 28 | 28 | identical |
| C_WI058_18 | 23 | 23 | identical |
| C_WI059_18 | 24 | 24 | identical |
| C_WI060_18 | 27 | 27 | identical |
| C_WI062_18 | 26 | 26 | identical |
| C_WI063_18 | 27 | 27 | identical |
| C_WI065_18 | 27 | 27 | identical |
| C_WI067_18 | 27 | 27 | identical |
| C_WI068_18 | 26 | 26 | identical |
| C_WI069_18 | 27 | 27 | identical |
| C_WI070_18 | 28 | 28 | identical |
| C_WI072_18 | 21 | 21 | identical |
| C_WI073_18 | 15 | 15 | identical |
| C_WI074_18 | 26 | 26 | identical |
| C_WI075_18 | 24 | 24 | identical |
| C_WI076_18 | 26 | 26 | identical |
| C_WI078_18 | 27 | 27 | identical |
| C_WI079_18 | 29 | 29 | identical |
| C_WI080_18 | 15 | 15 | identical |
| C_WI081_18 | 25 | 25 | identical |
| C_WI082_18 | 14 | 14 | identical |
| C_WI083_18 | 23 | 23 | identical |
| C_WI084_18 | 28 | 28 | identical |
| C_WI085_18 | 28 | 28 | identical |
| C_WI087_18 | 29 | 29 | identical |
| C_WI090_18 | 27 | 27 | identical |
| C_WI092_18 | 28 | 28 | identical |
| C_WI093_18 | 29 | 29 | identical |
| C_WI094_18 | 25 | 25 | identical |
| C_WI095_18 | 25 | 25 | identical |
| C_WI096_18 | 28 | 28 | identical |
| C_WI097_18 | 27 | 27 | identical |
| C_WI098_18 | 26 | 26 | identical |
| C_WI099_18 | 29 | 29 | identical |
| C_WI100_18 | 26 | 26 | identical |
| C_WI101_18 | 27 | 27 | identical |
| C_WI102_18 | 29 | 29 | identical |
| C_WI103_18 | 26 | 26 | identical |
| C_WI107_18 | 27 | 27 | identical |
| C_WI109_18 | 28 | 28 | identical |
| C_WI111_18 | 20 | 20 | identical |
| C_WI114_18 | 27 | 27 | identical |
| C_WI118_18 | 19 | 19 | identical |
| C_WI121_18 | 18 | 18 | identical |
| C_WI122_18 | 25 | 25 | identical |
| C_WI123_18 | 20 | 20 | identical |
| C_WI125_18 | 24 | 24 | identical |
| C_WI128_18 | 21 | 21 | identical |
| C_WI129_18 | 29 | 29 | identical |
| C_WI130_18 | 17 | 17 | identical |
| C_WI131_18 | 24 | 24 | identical |
| C_WI133_18 | 28 | 28 | identical |
| C_WI134_18 | 15 | 15 | identical |
| C_WI135_18 | 23 | 23 | identical |
| C_WI138_18 | 21 | 21 | identical |
| C_WI139_18 | 19 | 19 | identical |
| C_WI145_18 | 17 | 17 | identical |
| C_WI150_18 | 22 | 22 | identical |
| C_WI159_18 | 28 | 28 | identical |
| C_WI172_18 | 26 | 26 | identical |
| C_WI173_18 | 23 | 23 | identical |
| C_WI176_18 | 16 | 16 | identical |
| C_WI177_18 | 26 | 26 | identical |
| C_WI178_18 | 28 | 28 | identical |
| C_WI187_18 | 28 | 28 | identical |
| C_WI194_18 | 21 | 21 | identical |
| C_WI197_18 | 28 | 28 | identical |
| C_WI198_18 | 29 | 29 | identical |
| C_WI200_18 | 11 | 11 | identical |
| C_WI202_18 | 10 | 10 | identical |
| C_WI208_18 | 28 | 28 | identical |
| C_WI209_18 | 16 | 16 | identical |
| S_MN002_18 | 27 | 27 | identical |
| S_MN003_18 | 28 | 28 | identical |
| S_MN013_18 | 19 | 19 | identical |
| S_MN015_18 | 28 | 28 | identical |
| S_MN020_18 | 28 | 28 | identical |
| S_MN023_18 | 29 | 29 | identical |
| S_MN024_18 | 28 | 28 | identical |
| S_MN026_18 | 15 | 15 | identical |
## 2. `junction_caller.py` vs `/home/anton/a1-cyclospora/results/<SPEC>/junction.tsv`

### Which input this needs — stated plainly

`a1_junction.py` was run by `a1_specimen.sh` on the **trimmed FASTQ pair** (`t1/t2.fastq.gz`),
which that script then deleted. It cannot be run from the BAMs in `results/`: those are
`bwa mem` alignments against `markers.fa`, which has 7 contigs
(Nu_360i2, Nu_378, Nu_CDS1-4, Mt_MSR) and **no junction contig**, filtered to MAPQ ≥ 20
proper pairs with the unmapped fraction discarded (`idxstats` `*  0  0  0`). Junction reads
are therefore absent from them. `a1_map.sh` says so in a comment, and it is measurable:

| input for C_IA031_18 | reads scanned | left_anchored | flank_reads | call |
|---|---|---|---|---|
| trimmed FASTQ pair | 384,530 | 1,111 | 1,294 | Mt_Cmt199.A_Junction_Hap_17, 8 fragments |
| uBAM of the same reads (`samtools import`) | 384,530 | 1,111 | 1,294 | identical, byte for byte |
| `C_IA031_18.filtered.bam` | 52,048 | 0 | **0** | `NO_CALL / marker_absent` |

The junction caller needs reads: the trimmed FASTQ, an unaligned BAM, or the unmapped
fraction. The marker BAM is not a valid input, and it fails *silently* — it looks like a
clean negative for every specimen.

### How the comparison was made

The trimmed FASTQ no longer exist on disk, so they were regenerated exactly:
raw FASTQ re-fetched from the ENA URLs in
`/home/anton/a1-cyclospora/manifest/gold_benchmark.tsv`, then
`fastp 0.23.4 --detect_adapter_for_pe --length_required 50` from
`/home/anton/a1-cyclospora/sif/fastp.sif` — the same binary and the same flags
`a1_map.sh` used. That the regeneration was faithful is itself confirmed by the result:
`reads_scanned`, `left_anchored`, `left_anchor_only`, `right_anchor_only`,
`spanning_reads`, `discordant_fragments`, `spanning_fragments` and `flank_reads` all
match the archived `junction_diag.json` exactly, for all 25 specimens where that file
was regenerated alongside the TSV.

### Result: 26 specimens, 26 byte-identical `junction.tsv`

Deliberately chosen to cover every code path: three length classes, the one NOVEL array
in the archive, and all four `flag` values including each `NO_CALL` sub-case.

| specimen | reference call | repeats | fragments | total spanning | flank reads | `diff` vs reference |
|---|---|---|---|---|---|---|
| C_IA013_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 10 | 10 | 86 | identical |
| C_IA018_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 5 | 5 | 179 | identical |
| C_IA029_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 1483 | 1494 | 31622 | identical |
| C_IA031_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 8 | 8 | 1294 | identical |
| C_IA045_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 4 | 4 | 2388 | identical |
| C_IA058_18 | NO_CALL | NA | 0 | 0 | 0 | identical |
| C_IA064_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 330 | 340 | 24262 | identical |
| C_IA075_18 | NO_CALL | NA | 0 | 3 | 1417 | identical |
| C_IA082_18 | NO_CALL | NA | 0 | 0 | 163 | identical |
| C_IA083_18 | NO_CALL | NA | 0 | 1 | 495 | identical |
| C_IL004_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 50 | 51 | 4103 | identical |
| C_IL012_18 | NO_CALL | NA | 0 | 0 | 4 | identical |
| C_IL019_18 | NO_CALL | NA | 0 | 2 | 196 | identical |
| C_IL020_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 3 | 4 | 1532 | identical |
| C_IL025_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 86 | 88 | 7971 | identical |
| C_IL031_18 | NO_CALL | NA | 0 | 0 | 260 | identical |
| C_IL034_18 | NO_CALL | NA | 0 | 3 | 7625 | identical |
| C_IL036_18 | NO_CALL | NA | 0 | 2 | 6972 | identical |
| C_IL040_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 176 | 179 | 8711 | identical |
| C_IL048_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 295 | 304 | 11864 | identical |
| C_WI130_18 | NOVEL | 6 | 9 | 9 | 577 | identical |
| C_WI208_18 | Mt_Cmt214.A_Junction_Hap_20 | 7 | 4 | 4 | 1222 | identical |
| C_WI177_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 25 | 25 | 890 | identical |
| C_WI200_18 | Mt_Cmt199.A_Junction_Hap_17 | 6 | 218 | 232 | 12209 | identical |
| S_MN002_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 21 | 21 | 572 | identical |
| S_MN013_18 | Mt_Cmt169.A_Junction_Hap_8 | 4 | 13 | 14 | 262 | identical |

```
identical : 26 / 26   (diff -q, whole file, including the header)
differs   :  0 / 26
```

Coverage of the edge cases in that set:

| case | specimen(s) |
|---|---|
| `Mt_Cmt169.A_Junction_Hap_8` (4 repeats) | C_IA013_18, C_IA018_18, C_IA029_18, C_WI177_18, S_MN002_18, S_MN013_18 |
| `Mt_Cmt199.A_Junction_Hap_17` (6 repeats) | 9 specimens |
| `Mt_Cmt214.A_Junction_Hap_20` (7 repeats) | C_WI208_18 |
| `NOVEL` array, 6 repeats, 2 mismatches from Hap_17 | C_WI130_18 |
| flag `marker_absent` (no flank reads at all) | C_IA058_18 |
| flag `marker_present_but_no_spanning_fragment` | C_IA082_18, C_IL012_18, C_IL031_18 |
| flag `below_threshold` (1-3 spanning, under `--min-reads 3` after fragment collapse) | C_IA075_18, C_IA083_18, C_IL019_18, C_IL034_18, C_IL036_18 |
| call sitting exactly on the floor (3 fragments) | C_IL020_18 |
| very high depth (1,483 fragments / 31,622 flank reads) | C_IA029_18 |
| discordant-mate fragments discarded | C_WI208_18 (1 of 5 in the test-data extract) |

The `--bam` path was verified separately on C_IA031_18: a uBAM built with
`samtools import` from the same trimmed FASTQ produces a `junction.tsv` byte-identical
to both the FASTQ run and the archived reference.

## 3. `build_hds_sheet.py` vs `/home/anton/pyeuk-bench/lofreq_arm/sheets/`

Built from the 203 per-specimen TSVs of section 1 and diffed against the sheets the
benchmark's own `build_sheets.py` wrote:

| sheet | dimensions | `diff` |
|---|---|---|
| `sheet_A_153_withnovel.txt` | 153 x 59 | identical |
| `sheet_A_153_knownonly.txt` (`--drop-novel`) | 153 x 31 | identical |
| `sheet_A_203_withnovel.txt` | 203 x 63 | identical |

The three-state encoding was checked directly on a 26-specimen sheet carrying both PART
and junction calls: `X` present, empty absent, entirely empty locus block = not called.
All four junction columns fall into the single `Mt_Junction` locus block, so a specimen
with no junction call has that whole block empty (C_IA058_18, C_IA075_18, C_IA082_18,
C_IA083_18, C_IL012_18, C_IL019_18, C_IL031_18, C_IL034_18, C_IL036_18), while
C_IA029_18 and C_WI200_18 show empty `Nu_360i2_PART_*` blocks from amplicon dropout with
their junction block populated.

## 4. What is NOT validated here

* Caller B (the LoFreq-informed arm) is not ported, by design.
* The junction validation used 26 specimens, not all 203, because it required
  re-downloading and re-trimming the raw reads (~150 MB/specimen). The 26 were chosen to
  exhaust the code paths rather than to be a random sample. Extending to all 203 is a
  matter of running the same two scripts over the full manifest.
* Downstream stages (PyEuk distance and clustering, the k=2 / ARI >= 0.94 target) are not
  exercised by these three CLIs and are not claimed here.
