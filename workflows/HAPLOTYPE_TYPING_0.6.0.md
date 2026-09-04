# The Galaxy workflow that runs PyEuk 0.6.0

State as of 2026-09-01, on `anton-slim`. Companion to
`tools/cyclospora/DEVELOPMENT-TIMELINE.md`, which records why the versions below are
what they are.

There is no PyEuk 0.6.9. The current release is **0.6.0**, tag `v0.6.0`, commit
`c3aefae`.

## What it does

Types a cohort of amplicon specimens from reads, end to end, with no curated haplotype
catalogue. The only inputs are paired FASTQ per specimen and the panel FASTA. Everything
else is derived from the cohort's own data.

Paired FASTQ → bwa-mem2 against the panel → MAPQ ≥ 20 and proper pairs → windows sized
from the cohort's own fragment-length distribution → per-specimen haplotypes read off
reads that span each window end to end → one sheet → PyEuk distance and clustering.

Three properties are deliberate and worth knowing before changing anything:

`define_windows` is a **reduction over the whole BAM collection**, not a per-specimen
step. Windows must be identical for every specimen, or the sheet's columns do not
correspond across rows. The collection enters whole; the resulting BED fans back out.

A haplotype is read off a read that spans the window **end to end**. It is an observation
on one molecule, not an inference across molecules. Per-site variant calling records
mutations independently of the molecule carrying them, so it cannot tell one strain
carrying N mutations from a mixture where a second strain contributes them.

`min_freq` is a **haplotype** frequency, not a per-site allele frequency. Thresholding
per-site frequency can accept some sites of a minor haplotype and reject others,
assembling a genotype no molecule carries. Thresholding haplotype frequency cannot.

## The authoritative file

`pyeuk: galaxy/haplotype_window_typing.gxwf.yml`, on the `galaxy-workflows` branch.

The workflow moved into the library repository on 2026-09-01. PyEuk 0.5.0 had already
moved the amplicon front end into the library, so every computational step except read
alignment is a `pyeuk` subcommand; the workflow that wires them belongs beside them. It is
not copied into this repository — one source, no drift. Three sibling variants sit beside
it:

| file | when to use it |
|---|---|
| `haplotype_window_typing.gxwf.yml` | windows derived from the cohort's own reads. **The default.** |
| `haplotype_window_typing_bed.gxwf.yml` | windows supplied as a BED, for cohorts whose source publication defines its own intervals |
| `haplotype_window_typing_frombam.gxwf.yml` | starts from aligned BAMs instead of reads |
| `haplotype_window_typing_frombam_bed.gxwf.yml` | both of the above |

Two rules, both learned the hard way:

`haplotype_window_typing.ga.SUPERSEDED` **must not be imported.** It still carries
`spanning_target: 0.7`, which forces the deprecated read-length percentile rule, and
`min_maf: 0.05`, which destroys called-locus information under the deployed weighting.

Do not round-trip the YAML through `yaml.dump`. It silently rewrites `build_sheet` to
`in: None, state: None`. Edit it as text.

## Tool versions, as deployed

Verified against the running server on 2026-09-01.

| step | tool id | version |
|---|---|---|
| 1 | `bwa_mem2_idx` | `2.3+galaxy0` |
| 2 | `bwa_mem2` | `2.3+galaxy0` |
| 3 | `samtools_view` | `1.22+galaxy2` |
| 4 | `haplotype_define_windows` | `0.3.0+galaxy3` |
| 5 | `haplotype_window_caller` | `0.3.0+galaxy3` |
| 6 | `haplotype_window_sheet` | `0.3.0+galaxy3` |
| 7 | `haplotype_pyeuk` | **`0.6.0+galaxy0`** |

Only step 7 depends on PyEuk. It resolves to the Singularity image
`/data/containers/haplotype-pyeuk-0.5.sif`, which carries PyEuk 0.6.0 pinned to `c3aefae`.
Steps 4-6 depend on `pysam` alone and have not moved since they were written.

The gxformat2 files pin no `tool_version` for the local tools. They pick up whatever the
toolbox currently offers, which is why importing today yields 0.6.0 without editing
anything.

The toolbox currently offers **one** version of `haplotype_pyeuk`. A stored workflow
pinning `0.5.0+galaxy0` will resolve to 0.6.0 if re-run. That is safe here: 0.6.0 is a
rename with no numeric change, verified by byte-identical output from both images on the
same 153-specimen sheet.

## Parameters you set per cohort

Exposed as workflow inputs. Defaults shown are what the file ships.

| input | default | what it controls |
|---|---|---|
| `reads` | — | `list:paired` collection. **The element identifier is the specimen id.** It becomes column 1 of the calls and the row id of the sheet, so it must match `[A-Za-z0-9_.-]+` |
| `panel` | — | Amplicon panel FASTA. Used twice on purpose: as the bwa reference, and as the sequence haplotype names are described against. A different dataset here silently renames every haplotype |
| `min_span` | 30 | Minimum total fully-spanning reads before a window is called at all. The first gate, and on *Cyclospora* the binding one |
| `min_freq` | 0.05 | Minimum haplotype frequency |
| `min_reads` | 10 | Minimum supporting reads per haplotype |
| `min_completeness` | 0.1 | Minimum fraction of loci a specimen must have called to enter the distance matrix |
| `cut` | `count` | `count` or `distance` |
| `linkage_threshold` | unset | Dissimilarity for `cut=distance`. Leave unset to calibrate from the data |
| `project_psd` | `false` | Project the distance matrix onto the PSD cone |

Pinned in the tool state, not exposed. Change these by editing the YAML:

| step | setting | value |
|---|---|---|
| `define_windows` | `sample` | 20 BAMs scanned |
| | `min_spanning` | 0.30 |
| | `window_min` | 40 |
| | `window_max` | 0 (unbounded; 100 on panels where linkage/error trade-off matters) |
| | `width_step` | 10 |
| `call_haplotypes` | `denoise_edits` | 1 |
| | `denoise_ratio` | 8.0 |
| `build_sheet` | `min_freq`, `min_reads`, `min_maf` | 0, 0, **0** |
| `pyeuk` | `k_min`, `k_max` | 2, 50 |
| | `report_excluded` | true |

`build_sheet`'s `min_maf` is **0, i.e. off**, and must stay that way. It was 0.05. Under
the corrected heterozygosity weighting the filter is redundant, and removing it took
*Cyclospora* ARI from -0.0057 to 0.9734.

## Before running on a new cohort

Five things do not transfer, in rough order of how much they bite.

**Sweep `min_span`.** 30 is a *Cyclospora* result, not a constant. On that cohort 24-42
is a flat plateau at ARI 0.9737 with one specimen misassigned; between 42 and 44 it falls
to 0.7536 with ten misassigned. The old default of 50 sat on the wrong side of that
cliff. A sweep needs no re-mapping — one permissive caller pass serves the whole grid.
Note this gate is harder on **mixed** specimens: k co-infecting genotypes split the
spanning reads k ways, so a mixture needs k times the depth of a clonal sample to clear
the same value.

**Check `window_max` against your amplicon length.** Window length trades linkage against
sequencing error. On a *Cryptosporidium* titration with known proportions: at 250 bp a
pure control matched its own string in 62% of reads and both 75:25 mixtures were reported
as a single haplotype, the minor component having fragmented below the frequency gate. At
100 bp the control reached 79% and every mixture resolved.

**Lower `min_freq` if you expect low-frequency components.** 0.05 came from a titration
whose minor components were 25-75%, so nothing below the gate was ever exercised. On a
cohort with 1% components, 0.05 makes them undetectable *by construction* — every
observed frequency then piles up just above the gate. Set 0.005 or lower and let
`min_reads` suppress noise instead.

**Pick the cut mode from the expected structure, not from the score.** `count` splits
into k groups and suits a closed investigation where every specimen belongs somewhere.
`distance` cuts at a fixed dissimilarity and returns unrelated specimens as singletons,
which suits surveillance. A cluster count cannot represent a structure that is mostly
singletons: on the CDC cohort, whose published truth is 93 groups with 79 singletons, the
count rule rejects every k that could reproduce it and returns 1.

**Raise `min_completeness` if there is a low-completeness tail.** Those specimens
otherwise produce distances pinned at the engine's no-shared-data ceiling. On CDC
AmpliSeq five specimens sit at 0.46-0.48 against a median of 0.945, and they generate six
pairs at distance 1.0 where the real maximum is 0.2888.

`project_psd` stays `false` with `cut=count`. Set it `true` together with
`cut=distance`, where the regularised geometry helps: CDC scores ARI 0.7952 with the
projection on and 0.6805 with it off.

## Outputs

| output | what it is |
|---|---|
| `filtered_bam` | Per-specimen BAM, MAPQ ≥ 20 and proper pairs, coordinate sorted |
| `windows` | The windows this cohort was typed on. **Part of the result, not an intermediate** — the haplotype names only mean anything relative to these intervals |
| `calls` | Per-specimen window haplotype calls, with read counts and frequencies |
| `sheet` | Rows = specimens, columns = observed haplotypes, `X` = present |
| `haplotype_map` | Column name → window, interval, and the content-derived haplotype string |
| `calls_long` | Long-format calls with read frequency. **This is where mixtures live**; the binary sheet cannot carry them |
| `distance_matrix` | PyEuk wIBS distances |
| `clusters` | Cluster assignment. Specimens failing completeness appear as cluster `-1` |

An entirely empty locus block in the sheet means **not called**. "Amplified and
reference-identical" is the ordinary haplotype `=`, so the two are distinguishable.

## Importing and running it

The instance is `https://anton-slim.tailf1b947.ts.net` externally, `http://localhost:8080`
from the head node.

The documented API route is `POST /api/workflows {"from_path": ...}`. **It does not work
right now.** The key in `/mnt/ssd/pv4_full/configs/env.sh` still holds the one created
2026-06-11; two newer keys were made on 25 and 26 August, and creating a key retires the
previous one. Every authenticated call returns
`{"err_msg":"Provided API key has expired.","err_code":401001}`. Copy the current key from
User → Preferences → Manage API Key into that file before relying on the API.

Until then, import through the UI: Workflows → Import → Upload file → the `.gxwf.yml`.
That is how the current copy was imported — stored workflow **73**, workflow **123**,
whose `haplotype_pyeuk` step resolved to `0.6.0+galaxy0`.

Verify what actually landed rather than trusting the toolbox:

```bash
# what the server believes the tool is
curl -s http://localhost:8080/api/tools/haplotype_pyeuk | python3 -m json.tool

# what version an imported workflow's steps resolved to
sudo -u galaxy psql -d galaxy -c \
  "SELECT order_index, tool_id, tool_version FROM workflow_step
   WHERE workflow_id = (SELECT latest_workflow_id FROM stored_workflow
                        WHERE deleted=false ORDER BY update_time DESC LIMIT 1)
     AND tool_id IS NOT NULL ORDER BY order_index;"
```

Two failure modes on this instance, both previously observed. The toolbox cache will
report a stale tool version as though a deployment had not happened. And a Galaxy job can
report `state=ok` having produced nothing — check output bytes, not job state.

## Related

- `tools/cyclospora/DEVELOPMENT-TIMELINE.md` — why each version exists, and what was measured
- `~/galaxy-ops-log.md` — every change ever made to shared instance state
- `pyeuk: galaxy/README.md` — the workflow files themselves, and what each parameter does
- `BRC-research: cyclospora/variant-encoding/galaxy/` — the copies the frozen runs used
- `~/git/pyeuk/CHANGELOG.md` — what changed in each library release
