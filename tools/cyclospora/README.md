# Cyclospora cayetanensis genotyping — standalone callers

Three argparse CLIs that turn Illumina amplicon reads of *Cyclospora cayetanensis* into
a CDC-format haplotype data sheet (HDS). They are the hardened, path-free versions of
the one-off benchmark implementations, written so a Galaxy wrapper can call them
directly.

| script | stage | in | out |
|---|---|---|---|
| `scripts/part_haplotype_caller.py` | 3a — per-PART haplotypes | marker-aligned BAM | long TSV of calls |
| `scripts/junction_caller.py` | 3b — Mt-Junction | **reads** (FASTQ or uBAM) | junction TSV |
| `scripts/build_hds_sheet.py` | 6 — wide sheet | the TSVs above | CDC HDS sheet |

Only `part_haplotype_caller.py` needs a third-party library (`pysam`); the junction
caller needs it only on its `--bam` path.

## Provenance

* `part_haplotype_caller.py` — `/home/anton/pyeuk-bench/lofreq_arm/extract_reads.py` +
  `call_haplotypes.py`, **caller A** (the direct read-collapse path). Caller B, the
  LoFreq-informed arm, is intentionally not ported: it was the benchmark's comparison
  arm and needs a per-specimen VCF.
* `junction_caller.py` — `/home/anton/a1-cyclospora/scripts/a1_junction.py`.
* `build_hds_sheet.py` — `/home/anton/pyeuk-bench/lofreq_arm/build_sheets.py`,
  generalised to read per-specimen files and extended to carry the junction locus.

## Which input goes where — the one thing that is easy to get wrong

`part_haplotype_caller.py` takes the **marker-aligned, filtered BAM** (bwa mem vs
`markers.fa`, MAPQ ≥ 20, proper pairs).

`junction_caller.py` takes **reads**, not that BAM. The junction locus is *not a contig
in `markers.fa`*, so junction reads never reach the marker BAM. Feeding it the filtered
BAM is not an error that announces itself — it silently returns
`NO_CALL … marker_absent` with `flank_reads=0` for every specimen. Measured on
`C_IA031_18`:

| input | reads scanned | flank_reads | call |
|---|---|---|---|
| trimmed FASTQ pair | 384,530 | 1,294 | `Mt_Cmt199.A_Junction_Hap_17`, 8 fragments |
| uBAM of the same reads | 384,530 | 1,294 | identical |
| `C_IA031_18.filtered.bam` | 52,048 | **0** | `NO_CALL / marker_absent` |

So: trimmed FASTQ, an unaligned BAM, or the unmapped fraction of a mapping run — any of
those. Never the marker BAM.

## Gates — pre-registered, do not tune

From `/home/anton/pyeuk-bench/lofreq_arm/THRESHOLD_PREREGISTRATION.md`:

* PART locus is **called** at `--min-span 50` fully-spanning reads. Below that the
  caller emits *no rows at all* for that PART; downstream that becomes an entirely empty
  locus block, which is how the HDS format encodes "not called".
* A PART haplotype is kept at `--min-freq 0.05` **and** `--min-reads 10`.
* The junction uses `--min-reads 3` spanning *fragments* and `--min-freq 0.05`. The low
  floor is not a relaxed gate: Tn5 cuts inside the AT-rich repeat array, so typical
  spanning depth is 2–15 fragments even at ~900k reads while thousands of reads carry
  one flank.

## The three-state sheet encoding

```
X      haplotype present
""     haplotype absent, at a locus that WAS called
""     ... for every column of one locus  ->  locus NOT CALLED (amplicon dropout)
```

There is no separate missing code; the empty block *is* the code, and PyEuk's dropout
handling keys off it. Its one irreducible ambiguity, inherited from CDC: a locus that
was called but whose haplotypes all failed the frequency/read gates also yields an empty
block. `build_hds_sheet.py --summary` reports empty blocks per specimen so this can be
audited, and `part_haplotype_caller.py --summary` reports the actual spanning depth per
PART, which resolves it.

Junction columns (`Mt_Cmt<len>.<x>_Junction_Hap_<n>`) all belong to the **single** locus
`Mt_Junction`. The `Cmt<len>` in a name is a length class, not a locus.

## Two caveats the junction caller must keep honouring

Both are verified against `MAPPING_JUNCTION_WITH_PRIMERS_FEB_2020.fasta` and both
silently corrupt results if dropped:

1. `Mt_Cmt127.A_Junction_Hap_2` is **135 bp on disk, not 127** — a truncated right flank
   plus 11 nt of Nextera adapter. Never key off the `Cmt<length>` in a name. The caller
   keys off repeat count and array sequence only, and holds the right anchor at ≤ 39 nt
   so Hap_2 still matches.
2. CDC's 2022 nomenclature contains classes absent from the Feb-2020 20-reference set
   (`Cmt139` = 2 repeats, `Cmt229` = 8). Any repeat count is therefore permitted, and an
   unmatched array is reported as `NOVEL` with its class — never dropped, never forced
   onto the nearest reference.

The anchor length is bounded on both sides: ≤ 39 for caveat 1, and ≥ 14 because at 12 nt
the left anchor also occurs *inside* `Mt_Cmt154.D`'s repeat array, which would corrupt
the phase. `--anchor` is range-checked.

## Usage

```bash
# 3a  per-PART haplotypes
part_haplotype_caller.py \
    --bam SPEC.filtered.bam --markers markers.fa --parts parts.bed \
    --known-haplotypes haplotypes78.fa \
    --min-span 50 --min-freq 0.05 --min-reads 10 \
    --specimen SPEC --out SPEC.parts.tsv [--summary SPEC.parts.summary.tsv]

# 3b  Mt-Junction  (READS, not the marker BAM)
junction_caller.py \
    --junction-ref junction.fa --specimen SPEC --out SPEC.junction.tsv \
    --min-reads 3 --min-freq 0.05 \
    --fastq trim_R1.fq.gz trim_R2.fq.gz          # or: --bam unaligned.bam
    [--emit-specimen-column] [--diag SPEC.junction.json]

# 6  wide CDC HDS sheet
build_hds_sheet.py --calls calls/ --out sheet.txt \
    [--specimens specimens_153.txt] [--drop-novel] [--summary sheet_summary.tsv]
```

`--emit-specimen-column` prepends a `specimen` column to the junction TSV. It is **off**
by default so the output stays byte-identical to the reference implementation; turn it
**on** for Galaxy, where the file reaching `build_hds_sheet.py` is called
`dataset_<n>.dat` and the id cannot be recovered from the path.

## Validation

`selftest.sh` runs all three on `test-data/` and diffs against the checked-in expected
outputs. The real proof is in `VALIDATION.md`: byte-identical reproduction of the
reference outputs on 203 specimens (PART) and 26 specimens (junction).
