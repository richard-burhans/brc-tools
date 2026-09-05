# syntax=docker/dockerfile:1.7
#
# A PULLABLE TOGA2 image, so the wrapper can name a container instead of a filesystem path.
#
# ⛔ WHY THIS EXISTS. hillerlab/TOGA2 publishes no Docker image, no bioconda package and no
# biocontainer -- only `supply/containers/apptainer.def`. So the wrapper names a SIF at a
# hard-coded path, and every deployment must build that SIF and place it where every compute
# node can see it before `tools/toga2/toga2.xml` can run at all.
#
# ⚠ THIS DOES NOT MAKE TOGA2 A USER-DEFINED TOOL. The container was never the only obstacle:
# TOGA2 drives Nextflow, which spawns its own job hierarchy, wants many cores and writable
# state, and runs for hours, while a UDT is a single sandboxed shell_command that gets
# GALAXY_SLOTS=1 and ~3.8 GB on usegalaxy.org -- measured to be below TOGA2's own hardcoded
# 6 GB preprocessing requirement, which is not settable from the CLI. A pullable image removes
# the SHARED-STORAGE requirement from the classic tool. That is the whole claim.
#
# ⛔ AND IT DOES NOT SIMPLY REPRODUCE THE DEF, because the def contradicts upstream's own
# declarations. TOGA2 ships `conda.yaml` asking for python>=3.12,<3.13, nextflow>=25.10.0 and
# rust>=1.90.0. The def honours none of it: it takes whatever python continuumio/anaconda3
# happens to carry, and installs Nextflow from get.nextflow.io at BUILD TIME -- i.e. whatever
# is newest the day you build. This builds the environment upstream ASKS FOR instead.
#
# ⚠ THE NEXTFLOW PARSER BUG IS FIXED UPSTREAM, IN v2.0.9, AND PINNING WAS NEVER THE FIX.
# Through v2.0.8, `nextflow_configs/execute_joblist.nf` put its `if (params.joblist == "NONE")`
# guard and its channel construction at the TOP LEVEL, which Nextflow's strict parser (25.10+)
# rejects: "Statements cannot be mixed with script declarations". v2.0.9 moved both inside
# `workflow {}`. Measured, all four tags x two Nextflow generations:
#
#     script            NF 24.10.5     NF 26.04.6
#     v2.0.6/7a/8       3 of 3 OK      compilation failed
#     v2.0.9            3 of 3 OK      3 of 3 OK
#
# So the constraint in conda.yaml is honoured AND the code parses, with nothing pinned by us
# and no NXF_SYNTAX_PARSER=v1 workaround. An earlier revision of this file pinned
# NXF_VER=24.10.5; that is what should be deleted, not reinstated. Do not "fix" a future
# parser failure by pinning Nextflow backwards -- check whether upstream fixed the script.

FROM condaforge/miniforge3:24.11.3-2 AS build

# ⚠ A TAG, NOT A COMMIT. Upstream can move a tag. Pin the SHA once the tag stops moving.
ARG TOGA2_TAG=v2.0.9
ENV DEBIAN_FRONTEND=noninteractive \
    UV_LINK_MODE=copy \
    PATH=/root/.local/bin:$PATH

# ⚠ BUILD TOOLCHAIN ONLY. `make` compiles C (chain_filter_by_id, two shared libs), CESAR2.0,
# Cython extensions and three Rust crates, so gcc and git are not optional. Nothing scientific
# comes from apt -- but note `make` itself also wgets six UCSC binaries over plain HTTP with no
# checksum, git-clones PRANK at its default branch HEAD, and downloads IqTree2 from a GitHub
# release. Those are upstream's unpinned inputs, reached through `make`, not ours.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential git curl ca-certificates openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# uv for the Python side: faster resolution and a hard failure on conflict rather than pip's
# silent backtracking. ⚠ upstream's requirements.txt is only PARTLY pinned -- the scientific
# core is exact (tensorflow==2.20.0, numpy==2.2.5) but 11 of 27 entries are floors or bare
# names (biopython>=1.85, intronic>=2.1.5, matplotlib, ete3, gitpython, ...), so this step is
# not reproducible on its own and `--strict` is what keeps it honest.
COPY --from=ghcr.io/astral-sh/uv:0.5.14 /uv /uvx /root/.local/bin/

RUN git clone --branch ${TOGA2_TAG} --depth 1 --recurse-submodules \
        https://github.com/hillerlab/TOGA2 /opt/TOGA2

WORKDIR /opt/TOGA2

# ⛔ UPSTREAM'S conda.yaml IS THE SOURCE OF TRUTH, NOT A HAND-WRITTEN LIST. If a future tag
# changes the constraint, this picks it up without anyone noticing a version drifted.
#
# ⚠ AND ITS PATH IS DISCOVERED, NOT HARDCODED, BECAUSE UPSTREAM MOVES IT. v2.0.8 has it at the
# repo root; v2.0.9 relocated it to assets/venv/conda.yaml (and requirements.txt to
# assets/requirements/). Since build.bash builds every tag, a hardcoded path fails the older or
# the newer half of the matrix -- v2.0.9 failed here with "EnvironmentFileNotFound". Failing
# loudly when neither exists is deliberate: a missing env spec must never silently fall through
# to whatever the base image carries, which is precisely the def's defect.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    set -eu; \
    spec=""; \
    for c in conda.yaml assets/venv/conda.yaml; do \
        if [ -f "$c" ]; then spec="$c"; break; fi; \
    done; \
    [ -n "$spec" ] || { echo "ERROR: no conda.yaml found in ${TOGA2_TAG}" >&2; exit 1; }; \
    echo "conda env spec: $spec"; \
    mamba env create --file "$spec" --yes

SHELL ["conda", "run", "--no-capture-output", "-n", "toga2", "/bin/bash", "-c"]

RUN set -eu; \
    req=""; \
    for r in requirements.txt assets/requirements/requirements.txt; do \
        if [ -f "$r" ]; then req="$r"; break; fi; \
    done; \
    [ -n "$req" ] || { echo "ERROR: no requirements.txt found in ${TOGA2_TAG}" >&2; exit 1; }; \
    echo "python requirements: $req"; \
    uv pip install --system --strict -r "$req"

# `make` == chmod install check build_c build_cesar build_cython build_rust train_models.
#
# ⚠ `make` DOES run check_dependencies.py -- three times (check_essentials, check_managers,
# and install_third_party) -- but its `check` target does NOT certify this image. `check` is
# `check_shell check_essentials check_managers` and it runs THIRD in `build`'s prerequisite
# list, i.e. BEFORE build_c/build_cesar/build_cython/build_rust have compiled anything. It
# also omits check_python and check_third_party entirely, and `check_shell` is a no-op that
# prints "passed" out of the else-branch of a shell syntax error. Hence the real gate below.
#
# ⚠ `train_models` RETRAINS the XGBoost orthology models at build time. RANDOM_STATE is seeded,
# but the models are a function of the resolved xgboost/numpy/sklearn -- they are built here,
# not shipped by upstream.
RUN make

# ⛔ THE RUNTIME ENVIRONMENT IS SET *BEFORE* THE CHECK, NOT AFTER. When these came afterwards the
# check ran without /opt/TOGA2/bin/prank and /opt/TOGA2/bin/iqtree2 on PATH and reported both as
# "missing from $PATH" -- while `make` had built them perfectly well, ten lines earlier in the same
# log. A verification step that runs in a different environment from the one it certifies tests
# nothing; this one only said so out loud because the missing programs happened to be named.
#
# ⚠ NXF_HOME MUST BE WRITABLE AND THIS PATH IS NOT, under a read-only rootfs or a SIF. It is set
# here only so the plugins baked in below are found; `toga2.xml` overrides it to the job working
# directory, which is the supported way to run this image. Nextflow cannot even print its version
# if NXF_HOME is unwritable and unpopulated.
ENV NXF_HOME=/opt/TOGA2/.nextflow \
    NXF_DISABLE_CHECK_LATEST=true \
    PATH=/opt/conda/envs/toga2/bin:/opt/TOGA2:/opt/TOGA2/bin:/opt/TOGA2/bin/prank:/opt/TOGA2/bin/iqtree2:/opt/TOGA2/bin/intronIC:$PATH

# ⛔ BYTE-COMPILE THE SITE-PACKAGES. uv does not do this and pip does, and the difference is
# visible to Galaxy: ete3 ships pre-3.12 escape sequences, so without .pyc files every
# interpreter start re-parses them and emits ~30 SyntaxWarnings BEFORE any program output. The
# wrapper's <version_command> is `toga2.py --version 2>&1 | head -1`, so what Galaxy recorded as
# the tool version was "SyntaxWarning: invalid escape sequence '\-'". Under a read-only SIF the
# .pyc can never be written at runtime, so this must happen at build time. Asserted below.
RUN python3 -m compileall -q /opt/conda/envs/toga2/lib/python3.12/site-packages || true

# Populate the Nextflow plugin cache at build time so a --network=none run does not try to
# fetch plugins on first use.
RUN nextflow -version >/dev/null 2>&1 || true

# ⛔ VERIFY, DO NOT ASSUME -- AND NAME THE SUBCOMMANDS. A container that builds and then emits a
# differently-shaped orthology_classification.tsv is worse than no container, because the failure
# is silent and reaches the science.
#
# ⚠ THE FIRST VERSION OF THIS LINE WAS `check_dependencies.py` WITH NO SUBCOMMAND, WHICH IS A
# NO-OP. It is a click app: bare, it prints its help and exits 0, so the step that exists to catch
# a silently-wrong image was itself silently doing nothing.
#
# ⚠ AND UPSTREAM'S CHECKER IS A WEAK INSTRUMENT: it verifies PRESENCE, not that anything runs.
# `python` compares only MINIMUM versions (its table says tensorflow>=2.8.0 while requirements.txt
# pins ==2.20.0), `third_party` checks exactly three programs -- intronIC, PRANK, IqTree2; the
# SpliceAiChecker in that file is defined and never used -- and it WARNS and exits 0 when one is
# missing. `missing_third_party.txt` is removed first because the checker writes it but never
# clears a stale one. None of it touches Nextflow, toga2.py, CESAR2.0, the Rust binaries, the
# Cython modules or the UCSC binaries: this gate passes on an image with all of those deleted.
# That is why the smoke test below exists, and why it is the step that actually certifies.
RUN rm -f missing_third_party.txt \
    && python3 check_dependencies.py essentials \
    && python3 check_dependencies.py python \
    && python3 check_dependencies.py third_party \
    && if [ -s missing_third_party.txt ]; then \
           echo "ERROR: third-party programs missing from the image:" >&2; \
           cat missing_third_party.txt >&2; exit 1; \
       fi

# ⛔ THE SMOKE TEST IS THE REAL GATE: A COMPLETE TOGA2 RUN, AT BUILD TIME, ~60 SECONDS.
# `toga2.xml` claimed TOGA2 "cannot run a meaningful projection on tiny synthetic data". That is
# false: a 6 kb reference, a 6 kb query with a 200 bp indel, one 2-exon in-frame BED12 transcript,
# a one-line isoform table and a one-line chain produce all four declared outputs. This runs the
# wrapper's own command line, so it exercises Nextflow's parser, toga2.py, CESAR2.0, the compiled
# modules and faToTwoBit -- every one of which the dependency gate above is blind to.
COPY test-data /opt/TOGA2/test-data
RUN set -eu; \
    work="$(mktemp -d)"; cd "$work"; cp /opt/TOGA2/test-data/* .; \
    mkdir -p nxf_tmp; \
    export TMPDIR="$PWD/nxf_tmp" NXF_HOME="$PWD/.nextflow" NXF_ANSI_LOG=false; \
    faToTwoBit ref.fa ref.2bit; faToTwoBit query.fa query.2bit; \
    toga2.py run --ref_2bit ref.2bit --query_2bit query.2bit \
        --chain_file cleaned.chain --ref_annotation ref.bed12 \
        --isoform_file isoforms.tsv --no_u12_file --no_spliceai \
        --parallel_strategy local --feature_jobs 1 --orthology_jobs 1 --tree_cpus 1 \
        --output toga_out; \
    for f in query_annotation.bed query_genes.bed orthology_classification.tsv loss_summary.tsv; do \
        [ -s "toga_out/$f" ] || { echo "SMOKE FAILED: toga_out/$f missing or empty" >&2; exit 1; }; \
    done; \
    echo "SMOKE OK: all four declared outputs produced"; \
    cd /; rm -rf "$work"

# ⛔ ASSERT THE VERSION LINE GALAXY WILL RECORD. This is what the compileall step above is for;
# without it this line is a SyntaxWarning and Galaxy stores that as the tool version.
RUN v="$(toga2.py --version 2>&1 | head -1)"; \
    case "$v" in TOGA2,\ version*) echo "version_command OK: $v" ;; \
        *) echo "ERROR: version_command polluted, got: $v" >&2; exit 1 ;; esac

# Drop the clone's git history (~186 MB); nothing at runtime reads it.
RUN rm -rf /opt/TOGA2/.git

# ⚠ RESET BOTH INHERITED CONFIG KEYS. condaforge/miniforge3 sets ENTRYPOINT ["tini","--"], so an
# earlier version of this file asserting "NO ENTRYPOINT, DELIBERATELY" was simply false -- the
# image had one. Galaxy composes its own command line and execs it; tini is benign under
# `singularity exec` (which ignores ENTRYPOINT) but not something to leave asserted-and-untrue.
# SHELL is reset for the same reason: `conda run ...` above is baked into .Config.Shell and would
# be inherited by any downstream FROM. The conda env is on PATH, so nothing needs `conda run`.
ENTRYPOINT []
SHELL ["/bin/bash", "-c"]
WORKDIR /opt/TOGA2
CMD ["/bin/bash"]

# ⚠ image.version AND image.ref.name ARE OVERRIDDEN ON PURPOSE. Unset, the base image's values
# are inherited and the published image self-reports version="24.04", ref.name="ubuntu".
LABEL org.opencontainers.image.title="TOGA2" \
      org.opencontainers.image.version="${TOGA2_TAG}" \
      org.opencontainers.image.ref.name="toga2" \
      org.opencontainers.image.description="hillerlab/TOGA2 built from upstream's own conda.yaml and requirements.txt, verified with check_dependencies.py and a complete synthetic TOGA2 run" \
      org.opencontainers.image.source="https://github.com/hillerlab/TOGA2" \
      org.opencontainers.image.licenses="MIT"
