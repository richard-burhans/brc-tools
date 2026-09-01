# Development timeline: PyEuk 0.3.0 to 0.6.0, and the Galaxy tools that consume it

Written 2026-09-01, at the end of the development run that produced PyEuk 0.6.0.

## What this document is

This records **why** each version exists and **what was measured** to force it.

It does not record what changed. `pyeuk/CHANGELOG.md` already does that, release by
release, and it is complete for 0.3.0 through 0.6.0. Read that first. This file is the
half that is not written down anywhere else: the reasoning, the evidence, the rejected
options, and how the library releases reached a running Galaxy server.

Every claim below names a source. Where a source is a commit, the repository is given.
Where a source is a file outside git, that is stated, because those are the ones a reset
will lose.

## Where things live

The work spans four repositories and one directory that is not in git at all.

| repository | branch | what it holds |
|---|---|---|
| `spond/pyeuk` (local `~/git/pyeuk`) | `master` | the library. Tags `v0.3.0` to `v0.6.0`. `CHANGELOG.md`. `origin` is spond, `fork` is nekrut |
| `BRC-research` | `main` | `cyclospora/variant-encoding/` — Galaxy tool wrappers, container def files, gxformat2 workflows, cohort results, reviews |
| `brc-tools` | `cyclospora-lofreq-pyeuk` | the earlier LoFreq/PyEuk arm, the ToolShed-shaped tools, and this file |
| `PyEuk-paper` | `main` | the manuscript, its figures, and `HANDOVER.md` |
| `galaxy-ops` | `main` | `galaxy-ops-log.md`, one entry per change to shared instance state |

Not in git, and therefore at risk:

| path | what it holds |
|---|---|
| `/data/containers/haplotype-pyeuk-0.{1..5}.sif` | the built images. Def files for 0.3, 0.4 and 0.5 are tracked in `BRC-research`; 0.1 and 0.2 are not |
| `/data/galaxy_local/tools/` | the tool wrappers Galaxy actually reads. These are **copies**. Editing a repo changes nothing |
| `/data/cutadapt/RESULT.md` | the cutadapt evaluation write-up |
| `/data/plasmodium/*.tsv` | haplotype sheets used for regression checks |
| Galaxy PostgreSQL | the imported workflows and every frozen run |

## Version lineage

| date | PyEuk | container | `haplotype_pyeuk` tool | window tools |
|---|---|---|---|---|
| 2026-08-14 | v0.3.0 | 0.1, 0.2 (`master`) | — | — |
| 2026-08-25 | v0.4.0 | 0.3 (tag `v0.4.0`) | 0.4.0 | 0.3.0+galaxy0..3 |
| 2026-08-27 | v0.5.0 | 0.4 (commit `3dfb526`) | 0.5.0+galaxy0 | 0.3.0+galaxy3 |
| 2026-09-01 | v0.6.0 | 0.5 (commit `c3aefae`) | 0.6.0+galaxy0 | 0.3.0+galaxy3 |

The container numbering and the PyEuk version numbering are not the same sequence and
never were. Image 0.3 carries PyEuk 0.4.0. Image 0.5 carries PyEuk 0.6.0. Every image
was checked directly by importing the package inside it.

The window tools have not moved since they were written. They depend on `pysam` only,
not on PyEuk, which is why the PyEuk line moved four times underneath them without
touching them.

## Release by release

### v0.3.0 (2026-08-14) — the library becomes multi-pathogen

Extended beyond *Cyclospora* to *Cryptosporidium* and *Giardia*. Added `--ploidy`,
de novo haplotype learning, long-read ingestion, and the Gram PSD projection. Packaged
with `pyproject.toml` for the first time, which is what made a container pin possible.

Containers 0.1 and 0.2 both carry this version. Their def files were never tracked.

### v0.4.0 (2026-08-25) — three defects that made the distance wrong

This release exists because the distance was measurably wrong in three independent ways.
The rationale is preserved in the header of
`BRC-research: cyclospora/variant-encoding/galaxy/containers/haplotype-pyeuk-0.3.def`.

**Unobserved allele columns inflated the denominator** (`37f59a3`, spond/pyeuk#8).
Without the fix, distances are compressed roughly fivefold.

**The weighting was inverted for presence/absence input** (`1f61e15`, spond/pyeuk#11).
The KING form `1/sqrt(p(1-p))` is a standardisation for centred genotype dosages. On
presence/absence columns it is minimised at `p = 0.5` and diverges as `p` approaches 0,
so it weights rare alleles hardest while the discriminating signal sits in balanced
columns. The default became heterozygosity, `2p(1-p)`. On the 153-specimen *Cyclospora*
sheet this moved ARI from -0.0057 to the 0.95-1.00 range.

**The minimum-cluster-size guard made `k > 1` unreachable for small cohorts**
(`13b808d`). The guard was `max(5, ...)`, so no cohort of `n <= 9` could ever return more
than one group, and `n = 10` accepted only an exact 5/5 split. It became
`max(2, min(5, int(0.10 * n)))`. A 10-specimen, three-species cohort had been collapsing
to `k = 1` on that integer comparison alone.

**`k_min` was ignored in label-free mode** (`a27b403`, spond/pyeuk#12). The 0.2200
relative-gap floor also became overridable. Without this, `k = 3` was unreachable on the
*Cryptosporidium* three-species cohort: its relative gap of 0.0417 fell below the floor
and was rejected.

### v0.5.0 (2026-08-27) — the front end moves into the library, and a second way to cut

Two things happened.

The amplicon front end became library subcommands rather than only Galaxy scripts:
`define-windows`, `call-haplotypes`, `build-sheet`, behind the `[amplicon]` extra so the
core still installs without `pysam`.

A second cut mode was added, and this one was forced by a cohort the first mode cannot
represent. The rationale is in the header of `haplotype-pyeuk-0.4.def`. A cluster count
cannot describe a cohort whose true structure is mostly singletons. On the CDC
*P. vivax* AmpliSeq cohort the count rule returns `k = 1` and ARI 0.0000 against the
published assignment. The distance cut returns ARI 0.7952. Count mode is for a closed
investigation where every specimen belongs somewhere. Distance mode is for surveillance.

### v0.6.0 (2026-08-27, deployed 2026-09-01) — a rename, and nothing else

The import package became `pyeuk`. `cyclospora_pyeuk` survives as a deprecation shim.

**This release changes no number.** The diff from the 0.5.0 pin (`3dfb526`) to `v0.6.0`
is three commits. Changed lines outside the rename:

| file | changed lines | what they are |
|---|---:|---|
| `clustering.py` | 0 | — |
| `micro_assembly.py` | 0 | — |
| `naming.py` | 0 | — |
| `ont_processor.py` | 0 | — |
| `distance_engine.py` | 2 | import path |
| `haplotype_sheet.py` | 2 | import path |
| `cli.py` | 14 | import paths |
| `amplicon/define_windows.py` | 8 | text of the "pysam is missing" error |
| `amplicon/window_haplotypes.py` | 8 | same |

Verified on real data rather than argued from the diff. Both images were run against
`/data/plasmodium/r_cyclo_sheet.tsv`, 153 specimens, `min_completeness 0.1`, `k 2..50`,
`project_psd false`, `cut count`. The distance matrix is byte-identical at 178,021 bytes.
The cluster assignment is byte-identical at 2,013 bytes. Both give `k = 2`, threshold
2.2512, sizes 56 and 97, nothing excluded.

The rename forced one wrapper change. The shim calls
`warnings.warn(..., DeprecationWarning, stacklevel=2)`, and Python's default filter shows
that for a `__main__` caller. The line lands on stderr, and `haplotype_pyeuk` runs under
`detect_errors="aggressive"`. So `pyeuk_cluster.py` now imports `pyeuk` directly, in both
the deployed copy and the ToolShed copy.

## Decisions taken on the Galaxy side

These never reached the library. They are workflow and wrapper decisions, and each is
recorded in a `BRC-research` commit.

**Window placement is searched, not assumed** (`4eab54a`, `df039e6`). The first
implementation tiled from the first covered base and only chose a width. That is wrong
for an amplicon panel, because the two strands do not start at the same base. Placement
is now searched together with width, and candidates are scored by expected spanned bases
rather than by width alone.

**`min_span` default 50 to 30, and exposed** (`735e06a`). The value changes results, so
it stopped being a hidden constant.

**`min_maf` default 0** (`3995241`, and earlier `5260d80`, `610a0e3`). Minor-allele
filtering was the single largest defect found in adversarial review. On the *Cyclospora*
cohort it took ARI from -0.0057 to 0.9734 once removed. Under the corrected
heterozygosity weighting the filter is redundant (`d3fcd7c`).

**Denoising default of 1 edit** (`e322185`), set from the *Cryptosporidium* mixture
titration rather than chosen.

**Reads are streamed per window** (`4a7a3d5`). Peak RSS went from 815 MB to 30 MB. The
old behaviour was killing Galaxy jobs.

**Haplotypes are named by reference coordinate, not window offset** (`7d139db`).

**A specimen with no aligned reads must not fail the job** (`ea0e0c0`).

**bwa-mem2, with the index built once** (`9cca974`).

**The PSD projection is off by default** (`4714a7e`). The projection makes the distances
embed in Euclidean space, which Ward assumes, but it compresses the closest pairs — and
tight pairs are what a relatedness measure exists to resolve. Off, *Cyclospora* went from
3 errors in 153 to 1, and PvAmpSeq from 5 in 277 to 2. Turn it on together with a
distance cut, where the regularised geometry helps instead.

## Tried and rejected

**cutadapt as a pre-alignment step.** Tested 2026-08-28 on *Teladorsagia* beta-tubulin
(3 specimens) and PvAmpSeq (4 specimens). Write-up at `/data/cutadapt/RESULT.md`, which
is **not in git**. *Teladorsagia* calls were identical, 3 of 3 cells. PvAmpSeq matched on
47 of 48 cells. Cutadapt discarded 11.7% of PvAmpSeq pairs as untrimmed. It does remove a
real 1 bp deletion artefact at *Teladorsagia* position 194, present in 100% of raw reads,
but that position lies outside every called window. Cutadapt does not detect primers; the
sequences must be supplied by hand, which would weaken the claim that the pipeline needs
nothing beyond reads and a reference. Verdict: keep it out, and fix the aligner edge
effect in `define_windows` instead.

**`spanning_target 0.7`.** The percentile rule it forces is deprecated. It survives only
in `haplotype_window_typing.ga.SUPERSEDED`, which `galaxy/WORKFLOWS.md` says must not be
imported.

**A variant-frequency gate at 2%** (`df503d7`). It makes things worse, not better.

**Count-mode clustering for surveillance cohorts.** See v0.5.0 above. ARI 0.0000.

## Things to carry forward

**The headline ARI is in-sample.** The label-free ARI of 0.9737 is an in-sample number.
`sweep/honest_sweep.py` puts the out-of-sample estimate at **0.86-0.90**. Both
`RESTART.md` and `review/adversarial-benchmark.md` in `BRC-research` state this
correctly. The out-of-sample figure is the honest headline.

**Galaxy tool wrappers are copies.** The repo is not what runs. Deploying means: copy
into `/data/galaxy_local/tools/`, bump the version, restart the service, then verify
through the API. The toolbox cache will otherwise report the old version as if nothing
happened.

**Verify bytes, not job state.** A Galaxy job can report `state=ok` having produced
nothing.

**Gate sweeps need no remap.** One permissive caller pass serves a whole parameter grid.

**The API key in `/mnt/ssd/pv4_full/configs/env.sh` is stale.** It holds the key created
2026-06-11. Two newer keys were created on 25 and 26 August, and creating a key retires
the previous one, so that file has been dead since 25 August. Every call returns
`{"err_msg":"Provided API key has expired.","err_code":401001}`. Reads can be done
against PostgreSQL directly; the fix is to copy the current key from
User to Preferences to Manage API Key.

**The toolbox now offers one version of `haplotype_pyeuk`.** The 0.6.0 bump overwrote the
file rather than adding a second one, so stored workflows pinning `0.5.0+galaxy0` — the
frozen paper runs among them — resolve to 0.6.0 if re-run. The byte-identical check above
is what makes that safe. The old wrapper is backed up beside the live one as
`haplotype_pyeuk.xml.bak-2026-09-01-1444`, and `haplotype-pyeuk-0.4.sif` is untouched, so
an exactly-pinned 0.5.0 can be restored as a second tool file if it is ever wanted.

**Bioconda `pyeuk` is not published yet.** PR #68487 is awaiting a maintainer, and the
IWC submission is queued behind it. The ToolShed copies already pin
`<requirement type="package" version="0.6.0">pyeuk</requirement>`, so they will work the
day it lands and not before.

## How to check any line in this document

```bash
# the library, release by release
git -C ~/git/pyeuk log --decorate --date=short --pretty="%h %ad %d %s" | head -20
cat ~/git/pyeuk/CHANGELOG.md

# the reasoning behind each container pin, in the def file headers
head -30 ~/git/BRC-research/cyclospora/variant-encoding/galaxy/containers/haplotype-pyeuk-0.{3,4,5}.def

# what each image actually contains
singularity exec /data/containers/haplotype-pyeuk-0.5.sif \
  python3 -c "import pyeuk; print(pyeuk.__version__)"
singularity inspect /data/containers/haplotype-pyeuk-0.5.sif | grep PyEuk_Ref

# the Galaxy-side decisions, newest first
git -C ~/git/BRC-research log --date=short --pretty="%h %ad %s" -- cyclospora/variant-encoding

# what the running server believes right now
curl -s http://localhost:8080/api/tools/haplotype_pyeuk | python3 -m json.tool

# every change ever made to shared instance state
cat ~/galaxy-ops-log.md
```
