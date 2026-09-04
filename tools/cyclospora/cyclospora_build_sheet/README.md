# cyclospora_build_sheet — CDC HDS wide sheet from per-specimen calls

Stage 6 of the *Cyclospora cayetanensis* LoFreq/PyEuk typing arm. Takes a **list
collection** of per-specimen caller TSVs and writes the single wide haplotype data
sheet (HDS) that PyEuk consumes.

| | |
|---|---|
| tool id | `cyclospora_build_sheet` |
| version | `1.0.0+galaxy0` |
| container | `/data/containers/cyclospora-typing-1.0.0.sif` |
| script | `build_hds_sheet.py` (bundled via `required_files`) |

## Inputs / outputs

| direction | name | type | required |
|---|---|---|---|
| in | `calls` | `data_collection`, `collection_type="list"`, `format="tabular"` | yes |
| in | `specimens` | `data`, `format="txt,tabular"` | no |
| in | `drop_novel` | `boolean` → `--drop-novel` | no, default off |
| in | `emit_summary` | `boolean` | no, default off |
| out | `sheet` | `tabular` | always |
| out | `summary` | `tabular` | only when `emit_summary` |

Collection elements may be PART-caller output, junction-caller output, or both for
the same specimens; the type of each file is decided from its **header line**, never
from its name.

## Reference data: none, and that is a decision

`markers.fa`, `parts.bed`, `haplotypes78.fa` and `junction.fa` are neither bundled
in this wrapper nor exposed as workflow inputs. **This tool needs none of them.**

Naming a haplotype — the only step in the arm that consults a reference — happens
upstream, in `part_haplotype_caller.py` (which matches against `haplotypes78.fa`) and
`junction_caller.py` (against `junction.fa`). By the time a call reaches this tool it
already carries its final name, whether that is `Nu_360i2_PART_A_Hap_1`, a novel
`Nu_360i2_PART_A_NOV_d850eb`, or a junction `Mt_Cmt199.A_Junction_Hap_17`. Building
the sheet is then pure aggregation over strings: group columns into loci by name,
mark presence, sort.

Handing this step a reference input would add a way for the sheet to disagree with
the calls it was built from — a second, later opinion about naming, with nothing to
reconcile it against — and would buy nothing. So the only file bundled here is
`build_hds_sheet.py`.

(For the record, where reference data *is* needed — the two callers — bundling via
`required_files` is the wrong call: `haplotypes78.fa` and `junction.fa` are CDC
nomenclature releases that change independently of the wrapper, and pinning them
inside a tool version would silently freeze the nomenclature. Those belong as
workflow inputs. That is a different tool's README to write.)

## Why the container, for a stdlib-only script

`build_hds_sheet.py` imports nothing outside the Python standard library, so on the
face of it this tool needs no dependency at all. It still declares the shared
Apptainer image, for two reasons:

1. The arm's mandate is one image for all four tools. Multiple conda
   `<requirement>` entries resolve to a merged `mulled-v1-<hash>` environment that
   may never have been built, and the job then exits 1 with **completely empty
   stderr** — a failure mode that is close to undebuggable across parallel jobs.
2. The sheet and the distance matrix PyEuk computes from it then come out of the
   same interpreter (Python 3.11.15), so the string handling that produces column
   names is the same code that later reads them.

The cost is a 223 MiB image on a step that would run under any Python 3. That is
cheaper than a second dependency mechanism in the same pipeline.

## Where the specimen id comes from

Galaxy names every dataset `dataset_<n>.dat`, so filenames carry no identity. In
order of preference:

1. **the file's own `specimen` column.** PART-caller output always has one. The
   junction caller writes one only when run with `--emit-specimen-column`
   (**Emit specimen column** in its wrapper) — turn that on for Galaxy. The column
   wins over everything else.
2. **the element identifier**, which the wrapper stages onto disk as
   `<identifier>.tsv` before invoking the script, so the script's filename fallback
   still works. Characters outside `[A-Za-z0-9_.-]` are folded to `_`; if two
   identifiers collide after folding, `ln -s` fails and the job stops, rather than
   one file being silently overwritten.

A single file carrying two different specimen ids is rejected, not merged.

## The three-state encoding, and the one thing not to break

```
X       haplotype present
""      haplotype absent, at a locus that WAS called
""      ... for EVERY column of one locus  ->  locus NOT CALLED (amplicon dropout)
```

The entirely empty locus block *is* the missing-data code; PyEuk's dropout handling
keys off it. Never write a placeholder into an empty block.

Its one irreducible ambiguity, inherited from CDC: a locus that was called but whose
haplotypes all failed the frequency/read gates also yields an empty block.
`emit_summary` reports empty blocks per specimen so this is auditable; the PART
caller's own summary output resolves which case it was, since it carries the actual
spanning depth.

Junction columns (`Mt_Cmt<len>.<x>_Junction_Hap_<n>`) all belong to the **single**
locus `Mt_Junction`. The `Cmt<len>` in a name is a length class, not a locus.

## Verification

Run through Galaxy on real data, not on the test-data fixtures:

- 153 per-specimen PART TSVs — produced by `../scripts/part_haplotype_caller.py` from
  the filtered BAMs under `/home/anton/a1-cyclospora/results/`, at the pre-registered
  gates (`--min-span 50 --min-freq 0.05 --min-reads 10`) — staged as one `list`
  collection whose element identifiers are the specimen ids, plus
  `lofreq_arm/specimens_153.txt` as the specimen list.
- Output compared byte for byte against the reference sheet.

| run | reference | result |
|---|---|---|
| `drop_novel` off | `lofreq_arm/sheets/sheet_A_153_withnovel.txt` | `cmp` clean, md5 `56d1f4be…`; 153 rows × 59 columns / 24 loci |
| `drop_novel` on | `lofreq_arm/sheets/sheet_A_153_knownonly.txt` | `cmp` clean; 153 rows × 31 columns / 23 loci |

Both jobs ran under `singularity -s exec --contain --cleanenv --ipc --pid --no-mount tmp`
on `/data/containers/cyclospora-typing-1.0.0.sif`, dispatched by Slurm to **node02** —
a compute node, not the head node. `__instrument_core_container` in the job directory
records `{"container_id": "/data/containers/cyclospora-typing-1.0.0.sif",
"container_type": "singularity"}`.

The "one file, one specimen" guard was checked live as well: a TSV carrying two
different ids in its `specimen` column fails the job with exit 1 and

```
calls/MIXED.tsv mixes 2 specimens; one specimen per file, or use a long table split beforehand
```

rather than silently merging two specimens into one row.

## Keeping the bundled script in step

`build_hds_sheet.py` here is a verbatim copy of `../scripts/build_hds_sheet.py`;
`required_files` cannot reach outside the tool directory, so the file has to live in
both places. They must not drift:

```bash
cmp tools/cyclospora/scripts/build_hds_sheet.py \
    tools/cyclospora/cyclospora_build_sheet/build_hds_sheet.py
```

## Deployment

```bash
sudo mkdir -p /data/galaxy_local/tools/cyclospora_build_sheet
sudo cp -r tools/cyclospora/cyclospora_build_sheet/. \
           /data/galaxy_local/tools/cyclospora_build_sheet/
sudo chown -R galaxy:galaxy /data/galaxy_local/tools/cyclospora_build_sheet
# register in /data/galaxy_local/local_tool_conf.xml, then, serialised:
flock /tmp/gxreload.lock -c 'curl -s -X PUT -H "x-api-key: $GALAXY_API_KEY" \
    "$GALAXY_URL/api/configuration/toolbox"'
```

The deployed copy is a build artifact; this directory is the source of truth.

TPV needs the job to run the image, which the default rule does not enable. The rule
covering the whole arm lives in `/srv/galaxy/config/tpv_rules_local.yml`:

```yaml
  (.*cyclospora.*):
    cores: 1
    mem: 8
    params:
      singularity_enabled: true
      singularity_volumes: "$job_directory:rw,$tool_directory:ro,$job_directory/outputs:rw,$working_directory:rw,/data:rw"
```
