# cyclospora_part_caller

Stage 3a of the *Cyclospora cayetanensis* LoFreq/PyEuk typing arm: one
marker-aligned BAM in, one long-format table of per-PART haplotype calls out.
Map it over a BAM collection to type a run.

The wrapper is a thin shell around `../scripts/part_haplotype_caller.py`, which is
the ported and byte-validated caller **A** of the benchmark implementation
(`extract_reads.py` + `call_haplotypes.py`). `part_haplotype_caller.py` in this
directory is a **symlink** to that shared script so the two cannot drift; the
deployment step dereferences it (`cp -rL`), so the copy under
`/data/galaxy_local/tools/cyclospora_part_caller/` is a real file.

| | |
|---|---|
| tool id | `cyclospora_part_caller` |
| version | `1.0.0+galaxy0` |
| container | `/data/containers/cyclospora-typing-1.0.0.sif` (pysam 0.23.3) |
| inputs | BAM (+ optional specimen id, reference set, gates, summary switch) |
| outputs | `calls` (tabular, 10 col); `summary` (tabular, 5 col, optional) |

## Reference data: bundled, not workflow inputs

`markers.fa`, `parts.bed` and `haplotypes78.fa` are committed under `refs/` and
shipped with the tool via `<required_files>`. A *From history* branch exists, but
it is not the default.

The three files total **13.7 KB** and are a fixed panel definition, not data:

```
808a1929c7aff841b000b4a74a72316aebe34185f4a54c31a37a9e669e52dcc8  refs/markers.fa
4448421ddcca223829383686700d123688ef663a7e27321ebcc97c27c19cc204  refs/parts.bed
46f8bf91a504939a761009e5a78153bc6c4ab9f390b40a030e838265f7e9d344  refs/haplotypes78.fa
```

(copied verbatim from `/home/anton/a1-cyclospora/refs/`, which is what every
measured number in this port was produced against).

Bundling wins here for one reason that outweighs flexibility: **the three files
are not independent, and a mismatched trio fails silently rather than loudly.**

- `parts.bed` names its intervals, and those names are what the haplotype records
  in `haplotypes78.fa` are keyed by (`<PART>_Hap_<k>`). Supply a BED whose column
  4 uses any other convention and every call comes back `NOVEL` — a plausible,
  well-formed, completely wrong answer.
- `parts.bed`'s coordinates only mean anything relative to the exact `markers.fa`
  the BAM was aligned against. A one-base difference shifts every PART and
  changes every called sequence.
- The failure surfaces as *degraded concordance*, not as an error, and the
  concordance target (precision 0.9150 / recall 0.9145 against CDC's own calls)
  is the only thing that would catch it.

The caller does defend the part it can: `--markers` is not used to build any
called base, only to assert that every BED contig exists in the reference and
every interval lies inside it. That turns "wrong BED entirely" into a hard error.
It cannot catch "right contigs, wrong convention in column 4", which is the case
bundling actually removes.

The flexible half is kept because the panel is not eternal: CDC's 2022
nomenclature already contains haplotype classes absent from the 2020 reference
set, and a workflow that needs to type a different marker panel should not need a
new tool. Choosing *From history* exposes all three as data params, and they must
then be supplied as a matched set.

## Specimen id

Column 1 of the output is the specimen id, and it becomes the row id of the CDC
HDS sheet downstream. Leave the **Specimen id** field empty and the tool uses
`$input.element_identifier` (with a trailing `.bam` stripped), which is the right
answer when the tool is mapped over a named BAM collection. Set it explicitly for
a single dataset, whose element identifier is the Galaxy dataset name and usually
carries `.filtered.bam` or worse.

The id is charset-checked (`[A-Za-z0-9_.-]+`) in the command block before anything
else runs, so a surprising element identifier fails immediately instead of
producing a sheet row nobody can join on.

## Gates

Exposed, defaulted to the pre-registered values, and documented as pre-registered:

| gate | default | effect |
|---|---|---|
| `--min-span` | 50 | fewer spanning reads ⇒ the PART emits **no rows at all** |
| `--min-freq` | 0.05 | fraction of the PART's spanning reads |
| `--min-reads` | 10 | both floors must hold |

They are exposed so the tool is honest about what it applies, not so they can be
searched over. Every published figure for this arm holds at these values only.

## Why an empty PART matters

A PART below the span floor produces zero rows. Downstream that becomes an
entirely empty locus block in the HDS sheet, which is how the CDC format encodes
*not called* (amplicon dropout) as distinct from *called, haplotype absent*.
PyEuk's dropout handling keys off exactly that. Emitting a placeholder row would
destroy the distinction irrecoverably, so the tool does not.

The encoding has one inherited ambiguity: a PART that *was* called but whose
haplotypes all failed the frequency/read floors also yields an empty block. The
optional **per-PART depth summary** resolves it — it reports the real spanning
depth and the `called` flag for every PART in the BED, zero-coverage ones
included.

## Deployment

The repo is the source of truth; `/data/galaxy_local/tools/cyclospora_part_caller/`
is a build artifact.

```bash
sudo rsync -a --delete --copy-links --exclude test-data --exclude README.md \
    tools/cyclospora/cyclospora_part_caller/ \
    /data/galaxy_local/tools/cyclospora_part_caller/
sudo chown -R galaxy:galaxy /data/galaxy_local/tools/cyclospora_part_caller
# new tool ⇒ toolbox reload, serialised against the other agents on this host
flock /tmp/gxreload.lock -c \
  'curl -s -X PUT -H "x-api-key: $GALAXY_API_KEY" "$GALAXY_URL/api/configuration/toolbox"'
```

`--copy-links` is required: `part_haplotype_caller.py` is a symlink in the repo.

Registered in `/data/galaxy_local/local_tool_conf.xml` under the
`cyclospora` section.

### Job destination

The tool needs the Apptainer image, and on this instance a destination only runs
containers when its params say so. The TPV rule that makes that true for every
tool of this arm lives in `/srv/galaxy/config/tpv_rules_local.yml`:

```yaml
  (.*cyclospora.*):
    cores: 1
    mem: 4
    params:
      singularity_enabled: true
      singularity_volumes: "$job_directory:rw,$tool_directory:ro,$job_directory/outputs:rw,$working_directory:rw,/data:rw"
```

`/data:rw` is not optional: Galaxy's datasets (`/data/datasets`), job directories
(`/data/jobs`) and this tool directory all live under `/data`, and the default
volume list does not include it. Without it the job cannot see its own input BAM.

## Verifying by hand

`../selftest.sh` runs the underlying CLI on `../test-data/` and diffs against
checked-in expected output. The stronger evidence is `../VALIDATION.md`: the
caller reproduces the reference `calls_long.tsv` caller-A rows on **all 203
specimens**, 4854/4854 rows, zero field disagreements, including all 391 `NOVEL`
names.
