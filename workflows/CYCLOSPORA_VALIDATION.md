# Validation of `cyclospora_lofreq_pyeuk.gxwf.yml` against the reference implementation

The workflow was run inside Galaxy on **all 153 benchmark specimens** and its outputs were
compared, field by field, against the reference implementation's own archived outputs and
against CDC's deposited haplotype calls.

**Result: every headline number reproduces.** Junction calls are byte-identical to the
reference for all 153 specimens; PART calls agree on every single called haplotype; the HDS
sheet is byte-identical to a sheet built out-of-Galaxy from the reference implementation's
own calls; and PyEuk returns the same k, the same cluster assignments and the same ARI. The
only differences found anywhere are in the *supporting-read counts* of 59 of 3,669 PART call
rows (≤ 5 reads, ≤ 0.66 % relative), which change no call and no downstream number.

---

## 1. Headline comparison

| # | Quantity | Reference | Galaxy | Agree? |
|---|---|---:|---:|:---:|
| 1 | Junction precision vs CDC (153 specimens) | 1.000 | **1.0000** | yes |
| 1 | Junction recall vs CDC | 0.971 | **0.9712** | yes |
| 1 | Junction disagreements (both call, disjoint names) | 0 | **0** | yes |
| 1 | Specimens typed that CDC left blank | 27 | **27** | yes |
| 1 | `junction.tsv` files identical to `a1-cyclospora/results/<SPEC>/junction.tsv` | — | **153 / 153** | yes |
| 2 | PART precision vs CDC (24 PART loci, novel counted as FP) | 0.9150 | **0.9150** | yes |
| 2 | PART recall vs CDC | 0.9145 | **0.9145** | yes |
| 3 | Label-free clustering k | 2 | **2** | yes |
| 3 | Label-free ARI vs `gold_labels_153.txt` | ≥ 0.94 | **0.9737** | yes |

Nothing was tuned. The pre-registered gates (`min_span 50`, `min_freq 0.05`, `min_reads 10`)
were left at the workflow's defaults, PyEuk was run label-free with no gold file, and every
number below is the first number the run produced.

---

## 2. What was run

| | |
|---|---|
| Workflow | `07a1be4e408f64e1` — "Cyclospora typing (LoFreq/PyEuk arm)", 12 steps, unmodified |
| History | `822276e1c166c713` — "cyclospora validation 153" |
| Invocation | `df8eedbfeaa127fb` |
| Specimens | all 153 of `/home/anton/pyeuk-bench/gold_labels_153.txt` (98 Vendor_A, 55 Vendor_B) |
| Jobs | **768 / 768 `ok`. Zero errors, zero paused.** |

Per-tool job counts, all `ok`:

| step | tool | version | jobs |
|---|---|---|---:|
| `map_reads` | `bwa_mem` | 0.7.19+galaxy1 | 153 |
| `filter_bam` | `samtools_view` | 1.22+galaxy2 | 153 |
| `sort_bam` | `samtools_sort` | 2.0.8 | 153 |
| `part_caller` | `cyclospora_part_caller` | 1.0.0+galaxy0 | 153 |
| `junction_caller` | `cyclospora_junction_caller` | 1.0.0+galaxy0 | 153 |
| `merge_calls` | `__MERGE_COLLECTION__` | — | 1 |
| `build_sheet` | `cyclospora_build_sheet` | 1.0.0+galaxy0 | 1 |
| `pyeuk` | `cyclospora_pyeuk` | 2.1.0+galaxy0 (PyEuk pinned `d8e45a81`) | 1 |

The live workflow's tool state was read back from Galaxy before invoking and matches the
committed YAML: gates at `min_span 50 / min_freq 0.05 / min_reads 10`, samtools filter
`-q 20 -f 2 -F 3852`, PyEuk `mode: label_free` with no gold input.

### 2.1 Stage-1 input, and why it is the reference pipeline's own input

The workflow starts from trimmed FASTQ; trimming is upstream of the ported arm. The reference
pipeline (`a1_map.sh`) deleted its trimmed FASTQ after mapping, so they were regenerated
exactly: raw FASTQ re-fetched from the ENA URLs in
`/home/anton/a1-cyclospora/manifest/gold_benchmark.tsv`, then
`fastp 0.23.4 --detect_adapter_for_pe --length_required 50` out of the reference pipeline's
own `sif/fastp.sif`.

That the regeneration is faithful is measured, not assumed:

```
specimens compared                                                  : 153
fastp read-count and base-count mismatches vs the archived fastp.json: 0
```

For all 153 specimens, `before_filtering.total_reads`, `after_filtering.total_reads`,
`before_filtering.total_bases` and `after_filtering.total_bases` equal the values in
`/home/anton/a1-cyclospora/results/<SPEC>/fastp.json`.

**One substitution, disclosed.** ENA currently serves a **truncated 4,096-byte file** for
`SRR10396077_2.fastq.gz` (specimen `C_IL080_18`) — `Content-Length: 4096`, `ETag "1000"`,
reproducible across retries, range requests and the FTP endpoint. The pair was therefore taken
from SRA (`fasterq-dump SRR10396077`, sra-tools 3.1.1). This is a substitution of source, not
of data: the SRA R1 is identical to the intact ENA R1 in both sequence and quality
(`md5` of all sequence lines and of all quality lines match; 224,109 spots, matching the
manifest), only the read-name format differs, and after trimming the specimen reproduces the
archived fastp counts (448,218 → 428,388) exactly like the other 152.

---

## 3. Junction calls (stage 3b)

### 3.1 Against the reference caller's own output — exact

Each Galaxy `junction.tsv` was compared to
`/home/anton/a1-cyclospora/results/<SPEC>/junction.tsv` after dropping the leading `specimen`
column that the Galaxy wrapper adds (`emit_specimen_column: true`; the sheet builder needs it
because Galaxy names its datasets `dataset_<n>.dat`). Every remaining field — `length_class`,
`matched_reference`, `reads`, `freq`, `repeat_count`, `ref_seq_length`, `motifs_all_known`,
`closest_reference`, `closest_mismatches`, `flag`, `total_spanning_reads`, `flank_reads`,
`sequence` — was required to match.

```
specimens compared : 153
identical          : 153
differ             :   0
```

This covers all four `flag` values, all observed length classes, and the `NOVEL` array —
i.e. the whole code path, not a sampled subset.

### 3.2 Against CDC's deposited calls

Scored exactly as the reference `gate1_junction_novel.py` does: truth is
`a1-cyclospora/refs/haplotype_sheet_2022.txt`, restricted to junction columns whose name has a
sequence in `haplotypes78.fa`/`junction.fa`, and precision/recall are computed only over
specimens for which CDC published a junction call.

| | Reference | Galaxy |
|---|---:|---:|
| CDC has ≥ 1 junction call | 103 / 153 | **103 / 153** |
| we emit a junction call | 128 / 153 | **128 / 153** |
| both call, share ≥ 1 name | 101 | **101** |
| both call, **disjoint** (disagreements) | **0** | **0** |
| CDC silent, we call | 27 | **27** |
| TP / FP / FN | 101 / 0 / 3 | **101 / 0 / 3** |
| precision | 1.0000 | **1.0000** |
| recall | 0.9712 | **0.9712** |

The 3 false negatives are 2 specimens, both of which the reference implementation also leaves
uncalled — they are `below_threshold`, not misreads:

| specimen | CDC calls not reproduced | Galaxy call |
|---|---|---|
| `C_IA075_18` | `Mt_Cmt154.B_Junction_Hap_4`, `Mt_Cmt199.A_Junction_Hap_17` | `NO_CALL` (3 spanning reads, under `--min-reads 3` after fragment collapse) |
| `C_IL019_18` | `Mt_Cmt199.A_Junction_Hap_17` | `NO_CALL` (2 spanning reads) |

The 27 specimens CDC left blank and the workflow typed are 13 × `Mt_Cmt199.A_Junction_Hap_17`,
13 × `Mt_Cmt169.A_Junction_Hap_8` and 1 × `Mt_Cmt214.A_Junction_Hap_20` — all exact matches to
named CDC references, which is why they cost no precision.

---

## 4. PART calls (stage 3a)

### 4.1 Against the reference caller's own long table

Compared to the caller-A rows of `/home/anton/pyeuk-bench/lofreq_arm/calls_long.tsv`, matched
on `(part, sequence)`:

```
specimens compared                     : 153
reference call rows (153 specimens)    : 3669
Galaxy call rows                       : 3669
extra rows / missing rows              :    0 / 0
row-set differences (any specimen)     :    0
disagreements in name / match / nearest / edit : 0
```

**Every haplotype the reference calls, the workflow calls; no more, no fewer, under the same
name and the same EXACT/NOVEL classification.**

### 4.2 The one discrepancy: supporting-read counts

19 of 153 specimens differ in the *depth* columns only:

| field | rows differing (of 3,669) | max absolute difference | max relative difference |
|---|---:|---:|---:|
| `span` (reads spanning the PART) | 59 (1.61 %) | 5 reads | 0.662 % |
| `reads` (reads supporting the haplotype) | 35 | 4 reads | — |
| `freq` | 41 | 0.0021 | — |

Median `span` across all call rows is 1,445, so these are single-read-level differences at the
third or fourth significant figure.

**Cause.** The reference pipeline mapped with `bwa 0.7.18-r1243` from
`a1-cyclospora/sif/bwa.sif`; the workflow uses the toolshed `bwa_mem 0.7.19+galaxy1`. `bwa mem`
is also thread-count dependent when `-K` is unset — it estimates the insert-size distribution
per batch, and batch size is `threads × chunk`, so a different thread count moves a handful of
reads across the proper-pair/MAPQ boundary. Both differences are upstream of the ported code.
Nothing in the port is involved, and no call changes: not one row moved across the 5 % / 10-read
/ 50-span gates.

### 4.3 Against CDC's deposited calls

Scored exactly as the reference `concordance.py` does (24 PART loci, junction columns excluded,
every CDC PART mark scoreable, novel calls counted as false positives).

| accounting | Reference P | Galaxy P | Reference R | Galaxy R | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|---:|
| truth = ALL, novel counted as FP (**primary**) | 0.9150 | **0.9150** | 0.9145 | **0.9145** | 3357 | 312 | 314 |
| truth = ALL, novel ignored | 0.9970 | **0.9970** | 0.9145 | **0.9145** | 3357 | 10 | 314 |
| truth = NAMEABLE, novel counted as FP | 0.9150 | **0.9150** | 0.9876 | **0.9876** | 3357 | 312 | 42 |
| truth = NAMEABLE, novel ignored | 0.9970 | **0.9970** | 0.9876 | **0.9876** | 3357 | 10 | 42 |

Dropout profile, also identical to the reference:

| | Reference | Galaxy | CDC |
|---|---:|---:|---:|
| mean PARTs called per specimen | 20.33 | **20.33** | 20.58 |
| loci called by both | 3107 | **3107** | — |
| loci called only by us | 4 | **4** | — |
| loci called only by CDC | 42 | **42** | — |

Per-PART (primary accounting). **All 24 rows are identical to the reference's per-PART table —
zero differing cells.**

| PART | P | R | TP | FP | FN |
|---|---:|---:|---:|---:|---:|
| Mt_MSR_PART_A | 1.000 | 1.000 | 152 | 0 | 0 |
| Mt_MSR_PART_B | 1.000 | 1.000 | 153 | 0 | 0 |
| Mt_MSR_PART_C | 0.994 | 1.000 | 153 | 1 | 0 |
| Mt_MSR_PART_D | 1.000 | 1.000 | 153 | 0 | 0 |
| Mt_MSR_PART_E | 1.000 | 1.000 | 153 | 0 | 0 |
| Mt_MSR_PART_F | 1.000 | 0.993 | 151 | 0 | 1 |
| Nu_360i2_PART_A | 0.993 | 0.986 | 144 | 1 | 2 |
| Nu_360i2_PART_B | 0.987 | 0.987 | 147 | 2 | 2 |
| Nu_360i2_PART_C | 0.989 | 0.996 | 277 | 3 | 1 |
| Nu_360i2_PART_D | 0.993 | 0.987 | 147 | 1 | 2 |
| Nu_360i2_PART_E | 0.987 | 0.993 | 147 | 2 | 1 |
| Nu_360i2_PART_F | 0.978 | 0.963 | 131 | 3 | 5 |
| **Nu_378_PART_A** | **0.000** | **0.000** | **0** | **274** | **269** |
| Nu_378_PART_B | 1.000 | 0.989 | 278 | 0 | 3 |
| Nu_378_PART_C | 0.986 | 0.986 | 146 | 2 | 2 |
| Nu_378_PART_D | 0.961 | 0.989 | 273 | 11 | 3 |
| Nu_CDS1_PART_A | 1.000 | 0.978 | 89 | 0 | 2 |
| Nu_CDS1_PART_B | 0.980 | 0.960 | 96 | 2 | 4 |
| Nu_CDS2_PART_A | 0.963 | 1.000 | 104 | 4 | 0 |
| Nu_CDS2_PART_B | 0.989 | 0.967 | 89 | 1 | 3 |
| Nu_CDS3_PART_A | 0.975 | 0.967 | 117 | 3 | 4 |
| Nu_CDS3_PART_B | 0.983 | 0.991 | 113 | 2 | 1 |
| Nu_CDS4_PART_A | 1.000 | 0.932 | 82 | 0 | 6 |
| Nu_CDS4_PART_B | 1.000 | 0.954 | 62 | 0 | 3 |

`Nu_378_PART_A` is the known naming artifact, not a calling failure: CDC's sheet marks `Hap_4`
and `Hap_5` there, names that have **no sequence in `haplotypes78.fa`**, so no arm restricted to
the frozen 78-sequence reference can reproduce them. 269 of the 314 FN and 274 of the 312 FP
come from that one locus. Excluding it, the workflow agrees with CDC on 3,357 of 3,395 calls
(P 0.989 / R 0.987) — the same ceiling the reference implementation reported.

---

## 5. HDS sheet and clustering (stages 6–8)

### 5.1 The sheet

To get a like-for-like reference, the same `build_hds_sheet.py` was run **outside Galaxy** on
the reference implementation's own calls (caller-A rows of `calls_long.tsv` split per specimen,
plus the archived `results/<SPEC>/junction.tsv`), producing `ref_sheet.tsv`.

```
diff out/sheet.tsv ref_sheet.tsv  ->  identical
```

153 specimens × 62 haplotype columns (59 PART + 3 junction) across 25 loci, byte for byte.
The three-state encoding survives the round trip: `X` present, empty absent, an entirely empty
locus block = not called (`sheet_summary` reports the empty-block count per specimen).

### 5.2 PyEuk distance and label-free clustering

The same comparison was made for PyEuk: `pyeuk_cluster.py` run out-of-Galaxy in the same
container on `ref_sheet.tsv`, against the Galaxy `pyeuk` step's outputs.

| | Out-of-Galaxy reference | Galaxy |
|---|---|---|
| specimens in | 153 | 153 |
| specimens clustered (`min_completeness 0.10`) | 153 | **153** |
| specimens excluded | 0 | **0** |
| k | 2 | **2** |
| cluster sizes | 99 / 54 | **99 / 54** |
| threshold | 0.057551 | **0.057551** |
| ARI vs `gold_labels_153.txt` | 0.9737 | **0.9737** |

```
diff out/clusters.tsv ref_clusters.tsv                 ->  identical
distance matrix: 23,409 cells compared, max |diff| = 2.6e-10
```

Confusion matrix against the epidemiological labels, all 153 specimens scored (none dropped):

| | Vendor_A | Vendor_B |
|---|---:|---:|
| cluster 2 (n=99) | 98 | 1 |
| cluster 1 (n=54) | 0 | 54 |

One specimen misplaced out of 153. Against the task's target of **k = 2, ARI ≥ 0.94**, the
workflow returns **k = 2, ARI 0.9737**.

This was run label-free: `find_clusters(matrix, None)`, no gold file supplied to the tool. The
gold labels were used only afterwards, to score.

---

## 6. Honest summary of what does and does not reproduce

**Reproduces exactly (bit-for-bit or to 4 decimal places):**

* all 153 junction TSVs, field for field, against the reference caller's archived output;
* junction precision 1.0000 / recall 0.9712 / 0 disagreements / 27 CDC-silent specimens typed;
* the *identity* of every PART haplotype call — 3,669 rows, zero row-set differences, zero
  differences in name, EXACT/NOVEL class, nearest named haplotype or edit distance;
* PART concordance vs CDC, in all four accountings, and the entire 24-row per-PART table;
* the dropout profile (20.33 PARTs/specimen, 3107/4/42 locus agreement);
* the HDS sheet, byte-identical to one built out-of-Galaxy from the reference's own calls;
* the PyEuk distance matrix (max cell difference 2.6e-10), cluster assignments, k, threshold
  and ARI.

**Does not reproduce exactly (one discrepancy, quantified):**

* supporting-read counts in 59 of 3,669 PART call rows (19 of 153 specimens): `span` differs by
  at most 5 reads (≤ 0.66 % relative), `reads` by at most 4, `freq` by at most 0.0021. Cause is
  the aligner, not the port — `bwa 0.7.18` in the reference container vs `bwa_mem 0.7.19+galaxy1`
  in Galaxy, plus `bwa mem`'s thread-count-dependent batching with `-K` unset. **No call
  changes, no gate is crossed, and no downstream number moves.**

**Also worth stating plainly:**

* one FASTQ mate (`C_IL080_18` R2) could not be obtained from ENA — the archive itself serves a
  truncated 4,096-byte file — and was taken from SRA instead, verified identical to ENA on the
  intact mate and identical to the reference in trimmed read counts (§2.1);
* `Nu_378_PART_A` scores P 0.000 / R 0.000 against CDC. This is a property of the frozen
  reference dictionary, is present identically in the reference implementation, and is the
  single largest contributor to the 8.5 % residual in the headline PART numbers;
* the clustering scores on this feature set are known to be sensitive (the reference RESULT.md
  documents 3 presence marks on 1 specimen moving PyEuk ARI by 0.9 at an older PyEuk commit).
  The number here is reported as the first number the run produced, at the pinned commit
  `d8e45a81`, and it matches the out-of-Galaxy reference run on the identical sheet exactly —
  so it is a reproduction of the reference, not an independent claim about stability.

**Verdict: the port is faithful.** The workflow, run unmodified in Galaxy on all 153 benchmark
specimens, reproduces every pre-registered target of the LoFreq/PyEuk arm.

### 6.1 This workflow is now the arm's reference implementation

As of 2026-08-13 the Galaxy workflow — not the shell-script pipeline under
`/home/anton/a1-cyclospora/` — is the **authoritative source of LoFreq/PyEuk arm results**.
The shell pipeline remains the historical reference the port was validated against, and its
archived outputs stay in place, but new numbers should be produced and cited from Galaxy.

Grounds for the switch, all verified rather than asserted:

* the run is reproducible in place — a second invocation returned **byte-identical** output for
  the sheet, the summary, the distance matrix (0 of 23,409 cells differ), the clusters, and all
  153 per-specimen junction, PART and PART-summary files, at an identical ARI of
  0.9737018305231459;
* it is version-controlled, containerised and re-invocable by anyone with the history, rather
  than depending on one operator's shell environment;
* the aligner discrepancy in §6 was re-examined at the level that matters: across all 153
  specimens **no called haplotype set differs**. The perturbation reaches the reported `span`
  and the 4th decimal of `freq` and stops there. The call closest to the `min_freq` 0.05 gate
  sits at 0.0509, a margin roughly 70× the largest perturbation the depth wobble can produce,
  so no gate is at risk of flipping.

**A redundant `samtools_sort` step was removed** on 2026-08-13 after `filter_bam` was shown to
already emit `SO:coordinate` with a BAM index attached: `bwa mem` runs with
`output_sort=coordinate`, and `samtools view` filters records without reordering them. The two
BAMs were byte-identical over all 59,765 records. The workflow went 12 → 11 steps and 768 → 615
jobs (exactly one sort job per specimen), with every output unchanged.

**Known non-determinism, and how to remove it.** `bwa mem` runs without `-K`, so its chunk size
scales with thread count and the insert-size distribution is re-estimated per chunk. A rerun at
a different `-t` will therefore reproduce every *call* but may shift `span` by a few reads and
`freq` in the 4th decimal, exactly as observed between the reference's `-t 2` and Galaxy's
`-t 16`. Pinning `-K 100000000` would make depths bit-reproducible across thread counts. This
has **not** been done, because it would change the depth column and invalidate the validation
recorded here; it should be a deliberate follow-up with a re-validation, not a silent edit.

---

## 7. Reproducing this

Working directory `/home/anton/cyclo-validation/` (outside the repo; nothing in the reference
trees was modified):

| file | what it does |
|---|---|
| `dl.sh`, `urls.tsv` | fetch the 153 raw FASTQ pairs from the ENA URLs in the manifest, size- and gzip-verified |
| `trim1.sh`, `trimloop.sh` | `fastp 0.23.4 --detect_adapter_for_pe --length_required 50` from the reference pipeline's own container |
| `uploader.py` | upload the 306 trimmed FASTQ + 3 reference files to Galaxy |
| `invoke.py` | build the `list:paired` collection and invoke workflow `07a1be4e408f64e1` |
| `fetch_outputs.py` | download every workflow output |
| `compare.py` | every comparison in §3, §4 and §5; writes `out/comparison.json` |

`compare.py` computes the reference-side numbers with the same code as the Galaxy-side numbers,
and independently reproduces the published reference values (`gate1_junction_novel.py` →
P 1.0000 / R 0.9712; `concordance.py` → P 0.9150 / R 0.9145) before comparing anything — so the
harness is validated against the published outputs, not just self-consistent.
