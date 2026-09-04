#!/bin/bash
# Run the three callers on test-data/ and diff against the checked-in expected outputs.
#   PYTHON=/path/to/python-with-pysam ./selftest.sh
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3}"
T="$HERE/test-data"
C="$HERE/scripts"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT
rc=0

check() {  # check <label> <got> <want>
  if diff -q "$2" "$3" > /dev/null; then echo "ok    $1"; else echo "FAIL  $1"; diff "$3" "$2" | head -20; rc=1; fi
}

"$PY" "$C/part_haplotype_caller.py" \
    --bam "$T/parts_reads.bam" --markers "$T/markers.fa" --parts "$T/parts.bed" \
    --known-haplotypes "$T/haplotypes.fa" --specimen TEST --out "$W/parts.tsv" 2>/dev/null \
  && check part_haplotype_caller "$W/parts.tsv" "$T/expected_parts.tsv" \
  || { echo "FAIL  part_haplotype_caller (exit $?)"; rc=1; }

"$PY" "$C/junction_caller.py" \
    --specimen TEST --junction-ref "$T/junction.fa" --emit-specimen-column \
    --fastq "$T/junction_R1.fastq.gz" "$T/junction_R2.fastq.gz" --out "$W/junction.tsv" 2>/dev/null \
  && check junction_caller "$W/junction.tsv" "$T/expected_junction.tsv" \
  || { echo "FAIL  junction_caller (exit $?)"; rc=1; }

"$PY" "$C/build_hds_sheet.py" \
    --calls "$W/parts.tsv" "$W/junction.tsv" --out "$W/sheet.txt" 2>/dev/null \
  && check build_hds_sheet "$W/sheet.txt" "$T/expected_sheet.txt" \
  || { echo "FAIL  build_hds_sheet (exit $?)"; rc=1; }

# The marker-aligned BAM is the WRONG input for the junction caller and must fail loudly
# in the output rather than quietly look like a real negative.
"$PY" "$C/junction_caller.py" --specimen TEST --junction-ref "$T/junction.fa" \
    --bam "$T/parts_reads.bam" --out "$W/wrong.tsv" 2>/dev/null
if grep -q "marker_absent" "$W/wrong.tsv" && [ "$(cut -f12 "$W/wrong.tsv" | tail -1)" = "0" ]; then
  echo "ok    junction_caller on a marker BAM yields flank_reads=0 / marker_absent"
else
  echo "FAIL  junction_caller marker-BAM behaviour changed"; rc=1
fi

exit $rc
