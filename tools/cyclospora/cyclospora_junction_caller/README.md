# cyclospora_junction_caller — Galaxy wrapper

Stage 3b of the Cyclospora cayetanensis LoFreq/PyEuk typing arm: the Mt-Junction
(8th CDC marker) repeat-length caller. One read set in, one junction TSV out.

| | |
|---|---|
| tool id | `cyclospora_junction_caller` |
| version | `1.0.0+galaxy0` |
| engine | `junction_caller.py` (symlink to `../scripts/junction_caller.py`) |
| container | `/data/containers/cyclospora-typing-1.0.0.sif` |

`junction_caller.py` here is a **symlink** to `../scripts/junction_caller.py`. The
script in `tools/cyclospora/scripts/` is the single source of truth, exercised by
`tools/cyclospora/selftest.sh`; the wrapper must not fork it. Deploy with `cp -rL`
so the deployed tool directory gets a real file.

## Inputs

The `Reads to scan` conditional has three arms:

| arm | params | staged as |
|---|---|---|
| `fastq_pair` | `fastq1`, `fastq2` (`fastqsanger`, `fastqsanger.gz`) | `reads_R1.<ext>`, `reads_R2.<ext>` |
| `fastq_collection` | one `paired` collection | same |
| `bam` | `bam_input` (`unsorted.bam`, `qname_sorted.bam`, `bam`) | `reads.bam` |

Datasets are staged under names that carry a real extension because the caller
picks gzip vs plain text off the suffix; a Galaxy `dataset_N.dat` would be read
as uncompressed text and yield zero reads.

Other parameters: `specimen_name` (blank means "use the element identifier"),
the `reference` conditional, `min_reads` (3), `min_freq` (0.05),
`emit_specimen_column` (on), `emit_diag` (off).

`emit_specimen_column` defaults **on** here although the standalone CLI defaults
it off. In Galaxy the file that reaches the sheet builder is called
`dataset_N.dat` and the specimen id cannot be recovered from the path, so the id
has to travel inside the file. Turn it off only to reproduce the reference
implementation byte for byte.

## Outputs

* `junction` — tabular, always written. Columns: `specimen` (optional),
  `length_class`, `matched_reference`, `reads`, `freq`, `repeat_count`,
  `ref_seq_length`, `motifs_all_known`, `closest_reference`,
  `closest_mismatches`, `flag`, `total_spanning_reads`, `flank_reads`,
  `sequence`.
* `diag` — JSON read counters, only when `emit_diag` is on.

`reads` counts spanning **fragments**, not reads: both mates cover the same short
insert, so arrays are collapsed per fragment.

## Reference data: bundled, with an escape hatch

`refs/junction.fa` (4 KB, CDC `MAPPING_JUNCTION_WITH_PRIMERS_FEB_2020`, 20
records) ships **inside the tool** via `<required_files>` and is the default.

Why bundled rather than a workflow input:

1. **It is not a parameter, it is part of the algorithm.** The caller derives
   the LEFT/RIGHT constant flanks and both 18 nt anchors *from this file* — the
   shortest record is taken as the zero-repeat class. A different file silently
   changes what "spanning" means. That is not a knob a workflow should expose by
   default.
2. **The failure mode of a wrong file is invisible.** Anchors that do not occur
   in the reads produce `NO_CALL / marker_absent` for every specimen: a clean,
   plausible-looking negative. There is no crash to notice.
3. **It is static and tiny.** One CDC release since Feb 2020, 4 KB. Making every
   workflow carry a reference dataset, and every user pick the right one, buys
   nothing.
4. **Reproducibility is version-pinned.** The reference travels with the tool
   version, so `1.0.0+galaxy0` means one exact reference set.

The escape hatch is the `reference` conditional: choose *From my history* to
supply a newer CDC release without waiting for a tool version bump. Nothing about
the algorithm is hard-coded to the bundled contents; only the default is.

The two caveats the reference file forces on the caller (`Hap_2` is 135 bp on
disk despite its `Cmt127` name, and CDC's 2022 nomenclature has repeat classes
absent from this 20-record set) are documented in `../README.md` and in the
script's module docstring.

## Do not feed it the marker BAM

The junction locus is not a contig in `markers.fa`, so junction reads never reach
the marker-aligned, MAPQ≥20 BAM. That BAM produces `flank_reads=0` and
`NO_CALL / marker_absent` for every specimen and looks exactly like a real
negative. The second `<test>` pins that behaviour: wrong input must stay visibly
wrong, never an empty file.

## Verification

Beyond the two `<tests>`, the deployed tool was run once by hand on real
specimen `C_IA031_18` (SRR10395990, adapter-trimmed with the pipeline's own
`fastp 0.23.4 --detect_adapter_for_pe --length_required 50`). The Galaxy output,
with the `specimen` column removed, is **byte-identical** to the archived
reference `/home/anton/a1-cyclospora/results/C_IA031_18/junction.tsv`, and the
diagnostics JSON is identical to the archived `junction_diag.json`
(`reads_scanned=384530`, `flank_reads=1294`, `spanning_fragments=8`).

## Deployment

```bash
sudo cp -rL tools/cyclospora/cyclospora_junction_caller \
            /data/galaxy_local/tools/cyclospora_junction_caller
sudo chown -R galaxy:galaxy /data/galaxy_local/tools/cyclospora_junction_caller
# register in /data/galaxy_local/local_tool_conf.xml, then, under the shared lock:
flock /tmp/gxreload.lock -c 'curl -s -X PUT -H "x-api-key: $GALAXY_API_KEY" \
    "$GALAXY_URL/api/configuration/toolbox"'
```

### The container only runs if the destination enables it

`<container type="singularity">` is ignored unless the job destination sets
`singularity_enabled: true`. On this instance that comes from the
`(.*cyclospora.*)` rule in `/srv/galaxy/config/tpv_rules_local.yml`, and TPV
loads its config **at Galaxy startup**: `watch_job_rules` defaults to `false`,
so a rule added to that file after startup is not picked up and TPV keeps
serving the config it read when the process began.

This fails quietly. A non-containerized job runs `python3` from Galaxy's own
venv, which happens to have `pysam`, so both the FASTQ and the BAM path succeed
and produce correct numbers — with an unpinned interpreter. Check which one ran
instead of assuming:

```bash
curl -s -H "x-api-key: $GALAXY_API_KEY" "$GALAXY_URL/api/jobs/<job_id>?full=true" \
  | python3 -c "import json,sys; print(json.load(sys.stdin)['command_version'])"
```

`version_command` prints the interpreter version: **3.11.15** is the container,
anything else is not. Fixing it needs a Galaxy restart, not a toolbox reload.

Bump `@VERSION_SUFFIX@` on every interface change and verify with
`/api/tools/cyclospora_junction_caller?io_details=true`; a reload reports success
and keeps serving the old interface when the version has not moved.
