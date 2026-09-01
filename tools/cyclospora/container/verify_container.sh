#!/usr/bin/env bash
# Verify /data/containers/cyclospora-typing-1.0.0.sif.
#
# Runs the image with the same flags Galaxy uses (--contain --cleanenv
# --no-mount tmp), so /tmp is not writable and $TMPDIR is empty, exactly as on a
# compute node.
#
# Usage: verify_container.sh [sif] [workdir]
set -euo pipefail

SIF="${1:-/data/containers/cyclospora-typing-1.0.0.sif}"
WORK="${2:-$(mktemp -d)}"
BENCH="/home/anton/pyeuk-bench"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$WORK/out"
cp "$HERE/verify_pyeuk.py" "$WORK/"

apptainer exec --contain --cleanenv --no-mount tmp \
    --bind "$WORK":/work --bind "$BENCH":/bench:ro --pwd /work "$SIF" \
    bash -c '
        mkdir -p wrk && export TMPDIR="$(pwd)/wrk" &&
        python -c "import pysam, numpy, pandas, scipy, sklearn, numba; print(\"deps ok\")" &&
        python -c "from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine; print(\"pyeuk ok\")" &&
        python -c "from cyclospora_pyeuk.clustering import CyclosporaClusterFinder; print(\"cluster ok\")" &&
        python /work/verify_pyeuk.py /bench/haplotype_sheet_153.txt /bench/gold_labels_153.txt /work/out
    '

echo "work dir: $WORK"
