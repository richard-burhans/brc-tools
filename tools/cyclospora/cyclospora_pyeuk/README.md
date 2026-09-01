# cyclospora_pyeuk

Stage 7–8 of the *Cyclospora cayetanensis* amplicon typing arm: a CDC haplotype
data sheet goes in, a specimen × specimen distance matrix and a table of
outbreak cluster assignments come out.

The wrapper makes exactly two PyEuk calls and does no arithmetic of its own:

```python
matrix = PyEukDistanceEngine(epsilon=0.3072,
                             min_completeness=0.10).compute_revised_wibs_matrix(sheet_df)
clusters, k, threshold = CyclosporaClusterFinder(stringency=95.0,
                                                 robust=True).find_clusters(matrix, None)
```

`find_clusters(matrix, None)` — the second positional argument is the
gold-standard file — is the **label-free** path and the default. `pyeuk_cluster.py`
is a marshaller: read the TSV, validate it, call those two, write two TSVs,
print diagnostics.

## Inputs and outputs

| | name | format | notes |
|---|---|---|---|
| in | `sheet` | tabular | CDC HDS wide sheet. Column 1 `Seq_ID`, one column per haplotype name, `X` = present, empty = absent, an entirely empty locus block = **not called** |
| in | `calibration.gold` | tabular | optional, **off by default**. `Seq_ID` + `Cluster_alias`, with header |
| out | `distance_matrix` | tabular | square symmetric wIBS matrix, zero diagonal, `Seq_ID` header row and index column |
| out | `clusters` | tabular | `Seq_ID`, `Assigned_cluster` (1-indexed; `-1` = excluded by the completeness filter) |

Advanced: `min_completeness` (0.10), `k_min` (2), `k_max` (50),
`report_excluded` (on). Supervised only: `stringency` (95.0), `robust` (on),
`default_threshold` (0.05).

`epsilon` is deliberately **not** exposed. It parameterises PyEuk's *Bayesian*
distance; `compute_revised_wibs_matrix` never reads it, so a knob for it would
do nothing. It is fixed at the engine default, 0.3072.

## Reference data: this tool needs none, on purpose

The brief asks each tool in this arm to decide deliberately between bundling the
reference files (`markers.fa`, `parts.bed`, `haplotypes78.fa`, `junction.fa`)
with `<required_files>` and taking them as workflow inputs.

**For `cyclospora_pyeuk` the answer is neither, and that is not a dodge.** PyEuk
operates on haplotype *names* only. Sequences were consumed upstream: the PART
caller and the Mt-Junction caller resolved reads to `haplotypes78.fa` names (or
minted `<PART>_NOV_<md5[:6]>`), and the sheet builder froze those names into
column headers. By the time the sheet exists, the reference FASTA and BED can no
longer change any number this tool produces. Bundling them would be dead weight
that implies a dependency the tool does not have, and exposing them as inputs
would invite a user to supply a reference set inconsistent with the one the
sheet was named against, with no way for this tool to detect it.

The one file this tool *does* need at runtime is `pyeuk_cluster.py`, and that is
bundled with `<required_files>` in the house style.

The recommendation for the arm as a whole is bundling: the four reference files
total a few hundred kilobytes, they are static, and a workflow that lets a user
swap `haplotypes78.fa` silently changes what every haplotype name *means* across
specimens and across runs. The naming reference is part of the method, not a
parameter of it.

## Container, not conda

One Apptainer image for the whole arm:
`/data/containers/cyclospora-typing-1.0.0.sif`, built from
`../container/cyclospora-typing.def`. PyEuk is not on bioconda, and more to the
point two or more `<requirement>` entries in a Galaxy tool resolve to a merged
`mulled-v1-<hash>` conda environment that may never have been built; the job
then dies with exit 1 and a completely empty stderr.

Galaxy runs the image contained, with a clean environment and with `tmp` left
unmounted, so `/tmp` is not writable on a compute node and `$TMPDIR` arrives
empty. The `@SCRATCH@` macro therefore starts every command block with

```
mkdir -p wrk && export TMPDIR="$(pwd)/wrk" && export HOME=... NUMBA_CACHE_DIR=... MPLCONFIGDIR=...
```

`HOME` and `NUMBA_CACHE_DIR` matter here specifically: PyEuk's wIBS kernel is
Numba-JIT compiled, and Numba writes to a cache directory under `$HOME` if it is
not told otherwise.

This also means the destination has to enable Singularity. The default TPV rule
does not, and without `singularity_enabled: true` Galaxy silently ignores the
`<container>` tag and runs the job against whatever `python` the compute node
happens to have. The rule that covers the arm is in
`/srv/galaxy/config/tpv_rules_local.yml`:

```yaml
  (.*cyclospora.*):
    cores: 1
    mem: 8
    params:
      singularity_enabled: true
      singularity_volumes: "$job_directory:rw,$tool_directory:ro,$job_directory/outputs:rw,$working_directory:rw,/data:rw"
```

## Version pin

PyEuk is pinned to commit `d8e45a81820f28b6b3abf326e62870c346cd2890`. The
package's own `__version__` there is `2.1.0` (the brief called the release
2.1.7; the commit is the one intended). It is the last commit before
`f6e301d feat(distance_engine): Implement SNP-Weighted KING-wIBS`, i.e. the last
**binary** KING-wIBS engine.

Do not move the pin to v3.x. The SNP-weighted engine was benchmarked on the same
CDC sheets and regresses: AUC drops, 4 of 6 sheets collapse to k = 1, and the
null controls fall from 12/12 to 10/12.

## Label-free is the default, and why

Supplying gold-standard labels does not merely *score* the clustering — it
changes it. In the supervised path PyEuk calibrates the maximum tolerated
within-cluster distance from the labelled pairs (median + 3 × 1.4826 × MAD) and
then accepts the smallest k at which 95 % of within-cluster pairs sit below it.
A clustering fitted to labels and then scored against the same labels is not
evidence of anything.

The label-free path never sees a label: it reads the merge heights of the Ward
dendrogram, takes the largest relative gap above a 0.22 noise floor whose
resulting partition has no cluster smaller than max(5, 10 % of n), and cuts
there. That is the number reported below.

## Measured behaviour

On the CDC reference sheet `/home/anton/pyeuk-bench/haplotype_sheet_153.txt`,
label-free:

```
[PyEuk] Filtered dataset: 144 / 153 specimens passed completeness criteria.
[PyEuk-wIBS] Executing Numba JIT C-kernel on 144 specimens across 25 locus windows...
[ClusterFinder] Dendrogram Merge Height Gap Knee Detection: Optimal k = 2
    (Height Gap = 0.00347, Rel Gap = 0.3014, Min Cluster Size = 52 >= 14, Threshold = 0.00979)
k = 2, sizes {1: 52, 2: 92}, ARI = 0.9721 against the Vendor_A / Vendor_B labels
```

The ARI is computed *after* the fact, from `gold_labels_153.txt`, and is not
available to the clusterer.

Run through Galaxy on 2026-08-13 (job 35266, Slurm 34353, **node02**) the tool
reproduced that exactly: 144 / 153 passed the filter, k = 2, sizes {1: 52,
2: 92}, height gap 0.00347, rel gap 0.3014, threshold 0.00979; the 144 cluster
labels agree with the reference row for row and in the same order; post hoc
against the gold labels, ARI 0.9721, sensitivity 0.9905, specificity 0.9811.
The nine excluded specimens (`C_IA156_18`, `C_IL007_18`, `C_IL057_18`,
`C_IL082_18`, `C_IL119_18`, `C_OH006_18`, `C_WI134_18`, `C_WI200_18`,
`S_MN026_18`) appear as cluster `-1`.

## Test data

`test-data/sheet_40.tsv` is a **row subset** of `haplotype_sheet_153.txt`: 19
Vendor_A specimens, 19 Vendor_B specimens and the two lowest-completeness
specimens in the cohort (`C_IL119_18`, `C_OH006_18`), keeping the full
165-column header. Keeping every column is what makes the subset a faithful
test: it preserves the completeness denominator, so those two specimens still
fail the 0.10 filter and the `-1` reporting path is exercised at default
settings. Columns that are empty for all 40 specimens contribute nothing —
their whole locus window is uncalled, so PyEuk's pairwise-complete mask drops
them from every pair.

`expected_matrix_40.tsv` and `expected_clusters_40.tsv` are that subset's output
at default settings, kept for eyeballing; the `<tests>` block asserts on
contents rather than diffing whole files, because `numpy.linalg.eigh` in the PSD
projection is BLAS-dependent and the last couple of digits of the matrix are not
guaranteed to be identical across CPU architectures.

## What has and has not been verified

Verified on 2026-08-13, all four jobs on **node02** (a compute node), Slurm
34353 / 34363 / 34364 / 34365, `singularity_enabled: true` on the destination:

- the 153-specimen reference sheet reproduces the reference result exactly
  (see above);
- all three interface paths run through Galaxy on `test-data/sheet_40.tsv` and
  reproduce the same runs done by hand inside the container: label-free
  (k = 2, sizes 20/18, two specimens at `-1`), supervised with
  `gold_40.tsv` (AUC 0.7490, threshold 0.00260, k = 2, same partition), and
  `report_excluded = false` (39 lines, no `-1` rows). The Cheetah conditional
  is doing its job: `--gold` appears on the command line only in the supervised
  case;
- the container really is in use — `python3 -c "import cyclospora_pyeuk"` on
  node02's system Python fails with `ModuleNotFoundError`, so a job that
  imports it successfully cannot have run outside the image.

**Not** verified: the `<tests>` block has never been executed by planemo, which
is not installed on this host. Its expected values were taken from real runs of
the identical command inside the identical container and cross-checked against
the Galaxy runs above, but the test harness itself is unexercised.

## Floating-point reproducibility

The matrix is written with `float_format="%.10g"`. Ten significant digits is far
more than the analysis uses and keeps the file readable, but do not expect
byte-identical matrices across CPU generations: the PSD projection runs
`numpy.linalg.eigh` on the Gram matrix, and OpenBLAS picks kernels per CPU.

Measured, same container and same input sheet, head node (Core Ultra 7 265)
versus compute node02 (Xeon W-2255):

```
max |head - node02| over all 144 x 144 entries = 4.3e-11   (matrix values span 0 .. 6.2e-3)
cluster assignments                             identical, all 144 specimens, same order
```

That is ~7e-9 relative, larger than the ~1e-16 you would expect from pure
summation reordering, and the reason is structural rather than alarming: the
projection clips near-zero eigenvalues at exactly 0, which is discontinuous, and
then takes a square root, which amplifies a perturbation δ by 1/(2·d) ≈ 80 at
these distances. **Compare matrices with `atol=1e-9`, not byte-for-byte.**
Cluster assignments are integers and *are* expected to match exactly; if they
ever do not, the partition is sitting on a knife edge and the finding, not the
arithmetic, is what needs attention.
