# A pullable TOGA2 image

## Why

`hillerlab/TOGA2` publishes **no** Docker image, **no** bioconda package and **no** biocontainer —
only `supply/containers/apptainer.def`. Three consequences follow, and all three are why this
directory exists:

1. `tools/toga2/toga2.xml` must name a **filesystem path** (`/data/containers/toga2-2.0.9.sif`)
   rather than an image. Every deployment builds its own SIF and puts it on shared storage.
2. Nothing verifies that the SIF a job actually ran matches the def the wrapper was written
   against. A path is not a digest.
3. `tools/toga2/toga2.xml` exists and is registered (`galaxy_config_local_tool_conf.xml:52`), but
   **cannot run anywhere the SIF has not been built and staged first**. That is a deployment
   barrier, not a code one.

A published image removes those three — it is a **deployment** fix and nothing more.

⚠ **What this does NOT do: make TOGA2 a User-Defined Tool.** An earlier draft of this file claimed
it would, and that was wrong. The container was never the only obstacle. TOGA2 drives Nextflow,
which spawns its own job hierarchy, wants many cores and writable state, and runs for hours, while
a UDT is a single sandboxed `shell_command` that gets `GALAXY_SLOTS=1` and ~3.8 GB on
usegalaxy.org. That is not a matter of degree: TOGA2's preprocessing step requests **6 GB**
(`toga_main.py:1183`, not settable from the CLI) and dies outright — `CRITICAL: All preprocessing
step jobs died`, `req: 6 GB; avail: 3.7 GB`. It can be forced under the ceiling with
`--cesar_memory_limit 3`, but that **discards every projection above the limit before alignment**,
which is precisely the hard CESAR2 re-alignments the WF-C.4 rescue pass exists to recover — silent
at test scale, lossy at real scale. The TOGA2 rescue pass is absent from
`project_annotations_udt.gxwf.yml` for reasons a published image does not address.

## What this image is, and what it is not

⛔ **It is NOT a transcription of upstream's `%post`.** An earlier version of this file said it was,
in order, with deviations marked. That stopped being true when the Dockerfile was reworked: it is a
**re-derivation from upstream's own declared environment**, which the def itself does not use.

| | upstream's def | this image |
|---|---|---|
| base | `continuumio/anaconda3` | `condaforge/miniforge3:24.11.3-2` |
| python / nextflow / rust | whatever the base carries; Nextflow from `get.nextflow.io` at build time | `mamba env create` from upstream's **`conda.yaml`** |
| python packages | `pip install -r requirements.txt` (via `make`) | `uv pip install --strict` from the same file, then `make` |
| TOGA2 revision | default branch | `--branch ${TOGA2_TAG}` |
| `NXF_HOME` | `$TMPDIR/.nextflow`, which collapses to an unwritable `/.nextflow` when `TMPDIR` is empty | set to a populated path; the wrapper overrides it to the job dir |
| `%runscript` | appends SLURM entries to `/etc/passwd` with **typographic** quotes (U+201C/U+201D), so the lines are malformed | not carried over |
| `ENTRYPOINT` | n/a | explicitly reset to `[]` |

⚠ `Acquire::ForceIPv4` is **not** upstream's — it is in the *local* `toga2-2.0.9.def` because one
machine has no IPv6 egress. It is deliberately not here: baking it into a published image ships
somebody's network as a dependency. (A Docker build on the default bridge has no global IPv6
address, so `AI_ADDRCONFIG` never returns AAAA and the line would be a no-op anyway.)

⚠ `ENTRYPOINT []` is a **reset, not an omission**. `condaforge/miniforge3` sets
`ENTRYPOINT ["tini","--"]`, so a previous version of this file asserting "NOT ADDED · ENTRYPOINT"
described a property the image did not have — `docker inspect` reported `["tini","--"]`. `SHELL` is
reset for the same reason: the `conda run …` used for build steps is otherwise baked into
`.Config.Shell` and inherited by any downstream `FROM`.

## The Nextflow parser bug: fixed upstream in v2.0.9, and pinning was never the fix

Through **v2.0.8**, TOGA2's own `nextflow_configs/execute_joblist.nf` put its argument guard and
channel construction at the **top level**, which Nextflow's strict parser (25.10 and later) rejects:
`Statements cannot be mixed with script declarations`. **v2.0.9** moved both inside `workflow {}`.

| script | Nextflow 24.10.5 | Nextflow 26.04.6 |
|---|---|---|
| v2.0.6 / v2.0.7a / v2.0.8 | 3 of 3 ✔ | compilation failed |
| **v2.0.9** | 3 of 3 ✔ | 3 of 3 ✔ |

So the fix is the **tag**, not a pin. `conda.yaml` asks for `nextflow>=25.10.0` — a floor whose
satisfying set includes the versions that broke v2.0.8, which is why "honour the constraint" was not
by itself sufficient and why an intermediate revision of this image shipped a broken Nextflow. An
even earlier revision pinned `NXF_VER=24.10.5`, which *also* violated that floor. Both are wrong
answers to the same question. If a parser failure recurs, **check whether upstream fixed the script
before pinning anything.**

Upstream issue [#42](https://github.com/hillerlab/TOGA2/issues/42) (filed 2026-08-05) describes the
v2.0.8 failure; v2.0.9 was tagged 2026-08-18 and fixes it. The issue is still open.

⚠ v2.0.9 also **moved the environment files**: `conda.yaml` → `assets/venv/`, `requirements.txt` →
`assets/requirements/`. The Dockerfile discovers both rather than hardcoding either, because
`build.bash` builds every tag.

## Verification

Two gates, and only the second one is worth much.

**Upstream's `check_dependencies.py`** (`essentials`, `python`, `third_party`, plus a non-empty
`missing_third_party.txt` treated as an error, since upstream's own checker warns and exits 0).
⚠ This is a **weak** instrument: it verifies that three third-party programs are *present*, compares
only *minimum* versions, and touches neither Nextflow, `toga2.py`, CESAR2.0, the Rust binaries, the
Cython modules nor the UCSC binaries. It was measured to **pass on an image with `nextflow`,
`toga2.py` and `cesar` deleted.** Note also that the `SpliceAiChecker` in that file is defined and
never used — `third_party` asserts nothing about SpliceAI. `make`'s own `check` target is weaker
still: it runs *before* anything is compiled and omits `check_python`/`check_third_party`.

**A complete TOGA2 run, at build time (~60 s).** A 6 kb reference, a 6 kb query with a 200 bp
deletion, one two-exon in-frame BED12 transcript, a one-line isoform table and a one-line chain,
run through the wrapper's own command line, asserting all four declared outputs are non-empty. This
exercises Nextflow's parser, `toga2.py`, CESAR2.0, the compiled modules and `faToTwoBit`. The same
fixtures are the Galaxy `<test>` in `tools/toga2/toga2.xml` (`tools/toga2/test-data/`).

A third assertion checks that `toga2.py --version` emits `TOGA2, version …` and not a
`SyntaxWarning`. That is not cosmetic: the wrapper's `<version_command>` pipes it through
`head -1`, and because `uv pip install` does not byte-compile while `pip` does, `ete3`'s pre-3.12
escape sequences made **Galaxy record the warning text as the tool version**. Hence the
`compileall` step — under a read-only SIF the `.pyc` can never be written at runtime.

## Build and publish

```
./build.bash --tags v2.0.9                                    # build only
QUAY_ROBOT_USER=… QUAY_ROBOT_PASSWORD=… ./build.bash --push   # build every tag and push
```

`build.bash` logs in through a **disposable `DOCKER_CONFIG`**, not `~/.docker/config.json`: a plain
`docker login` leaves the robot token on disk as base64 with no credential helper, and `docker
logout` is not a reliable cleanup (it prints success, exits 0, and on a read-only config dir leaves
the file byte-identical with the secret still in it). It also disables `xtrace` across the
credential region, because `--password-stdin` protects the process table but does nothing about
`set -x` — which an exported `SHELLOPTS` can switch on from outside.

The digest is taken from `docker push`'s own output, not `docker inspect`: a never-pushed image
already has a `RepoDigests` entry, and `{{index .RepoDigests 0}}` can name a repository the script
never pushed to.

Then point `@TOGA_IMAGE@` in `tools/toga2/macros.xml` at the published reference. ⚠ Prefer a
**digest** (`@sha256:…`) over a tag: the whole point is that a job's container can be identified
after the fact, and a tag can move.

## ⚠ Before trusting it

The wrapper's CLI mapping and output schema were verified against **v2.0.8**. The pin has moved to
**v2.0.9**. The build-time smoke test confirms the four-file output *schema* still holds, but it is
a liveness-and-shape test on inputs too small for anything to be rejected for memory — it cannot
detect projection loss at real scale. Re-run the real WF-C.4 comparison before treating v2.0.9
output as interchangeable with v2.0.8.

## Known, and deliberately not fixed here

- The image is large (~5 GB compressed): the build stage is never copied out of, so the shipped
  image carries `src/rust/target`, `rustlib`, the JDK, `build-essential` and `uv`. Only the clone's
  `.git` is removed.
- `make` reaches unpinned upstream inputs: six UCSC binaries over plain **HTTP** with no checksum,
  PRANK cloned at its default branch HEAD, IqTree2 from a GitHub release. Fixing those belongs
  upstream.
- `train_models` **retrains** the XGBoost orthology models at build time. `RANDOM_STATE` is seeded,
  but the models are a function of the resolved xgboost/numpy/sklearn — they are not shipped
  artifacts.
- Base image and `uv` are pinned by **tag**, not digest.
