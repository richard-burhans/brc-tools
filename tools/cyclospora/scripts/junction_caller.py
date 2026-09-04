#!/usr/bin/env python3
"""Cyclospora Mt-Junction (8th CDC marker) repeat-length caller.

Ported from /home/anton/a1-cyclospora/scripts/a1_junction.py. The algorithm, the
constants (ARRAY_START=64, MOTIF_LEN=15, ANCHOR=18) and the output columns are
unchanged; what changed is the interface (--out is a FILE, not a directory) and the BAM
path, which now goes through pysam instead of shelling out to `samtools fastq`.

INPUT: THIS CALLER NEEDS READS, NOT AN AMPLICON ALIGNMENT
--------------------------------------------------------
The junction locus is NOT a contig in markers.fa, so junction reads are absent from the
marker-aligned BAM -- a filtered marker BAM yields flank_reads=0 and NO_CALL for every
specimen. Feed this caller either the trimmed FASTQ pair, or a BAM that still holds the
unaligned reads (an unaligned/uBAM, or the unmapped fraction of a mapping run). Reads
are scanned as raw strings on both strands; alignment information is never used.

Structure of the reference set (MAPPING_JUNCTION_WITH_PRIMERS_FEB_2020, 20 records):

    [ LEFT 64 nt constant ][ array = n x 15 nt ][ RIGHT 45 nt constant ]

with canonical amplicon length 109 + 15*n.

TWO CAVEATS, BOTH VERIFIED AGAINST THE REFERENCE FILE, BOTH LOAD-BEARING
  * Mt_Cmt127.A_Junction_Hap_2 is 135 bp ON DISK, not 127. Its repeat count is 1
    (canonical length 124); the surplus is a truncated RIGHT flank plus 11 nt of Nextera
    adapter. NEVER key off the Cmt<length> in the name. This caller keys off repeat
    count and array sequence only, and the right anchor is kept at <= 39 nt so Hap_2's
    truncated right flank still matches.
  * CDC's 2022 nomenclature contains classes (Cmt139 = 2 repeats, Cmt229 = 8) that are
    absent from the Feb-2020 20-reference set. Any repeat count is therefore permitted,
    and an array with no exact reference is reported as NOVEL with its class -- never
    dropped, never forced onto the nearest reference.

Method
------
A read is SPANNING if it contains the 18 nt left anchor (the last 18 nt of LEFT, which
ends exactly at the array start) and, after it, the 18 nt right anchor (the first 18 nt
of RIGHT). Reads matching only the reverse complement of the left anchor are reverse
complemented first. The intervening string is the repeat array; its length must be a
multiple of 15 and it must be N-free. Anchor length must stay <= 39 (Hap_2) and >= 14
(at 12 nt the left anchor also occurs inside Mt_Cmt154.D's array, which would corrupt
the phase). Because both mates of a fragment cover the same short insert, arrays are
collapsed per fragment; mates that disagree are discarded as discordant, so `reads` in
the output is really spanning FRAGMENTS.

Why --min-reads defaults to 3: the array is extremely AT-rich, which is what Tn5
prefers, so Nextera XT tagmentation cuts INSIDE the array and most fragments stop
there. Typical spanning depth is 2-15 fragments even at ~900k reads, while thousands of
reads carry one flank. That is a property of the library prep, not of this code, and it
is why CDC has no junction call for 50 of the 153 benchmark specimens. flank_reads is
reported on every row so "marker absent" can be distinguished from "marker present but
never spanned".

All variants passing the gates are emitted: mixed infections are real, and multiple
length classes in one specimen are expected output, not an error.

Usage
-----
  junction_caller.py --specimen SPEC --junction-ref junction.fa --out SPEC.junction.tsv \\
      --fastq trim_R1.fq.gz trim_R2.fq.gz
  junction_caller.py --specimen SPEC --junction-ref junction.fa --out SPEC.junction.tsv \\
      --bam unmapped.bam
"""
import argparse
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

MOTIF_LEN = 15
ARRAY_START = 64          # verified: the repeat array begins at position 64 in the refs
ANCHOR = 18               # see module docstring: must be <= 39 and >= 14

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

COLUMNS = ["length_class", "matched_reference", "reads", "freq", "repeat_count",
           "ref_seq_length", "motifs_all_known", "closest_reference",
           "closest_mismatches", "flag", "total_spanning_reads", "flank_reads",
           "sequence"]


def revcomp(s):
    return s.translate(_COMP)[::-1]


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else None


class JunctionRefs:
    """junction.fa -> LEFT / RIGHT constants, anchors, and an array -> reference map."""

    def __init__(self, path, array_start=ARRAY_START, anchor=ANCHOR):
        self.seqs = read_fasta(path)
        if not self.seqs:
            raise SystemExit("no sequences in %s" % path)
        self.array_start = array_start
        # The shortest record is the zero-repeat class, so it is exactly LEFT+RIGHT.
        shortest = min(self.seqs.values(), key=len)
        self.LEFT = shortest[:array_start]
        self.RIGHT = shortest[array_start:]
        self.LA = self.LEFT[-anchor:]
        self.RA = self.RIGHT[:anchor]
        self.LA_rc = revcomp(self.LA)
        self.RA_rc = revcomp(self.RA)

        self.array2refs = defaultdict(list)
        self.ref_array = {}
        self.warnings = []
        for name, seq in self.seqs.items():
            if not seq.startswith(self.LEFT):
                self.warnings.append("%s: does not start with the common LEFT flank" % name)
                continue
            idx = -1
            j = seq.find(self.RA, array_start)
            while j >= 0:                       # first IN-PHASE right anchor wins
                if (j - array_start) % MOTIF_LEN == 0:
                    idx = j
                    break
                j = seq.find(self.RA, j + 1)
            if idx < 0:
                self.warnings.append("%s: right anchor not found in phase" % name)
                continue
            arr = seq[array_start:idx]
            self.ref_array[name] = arr
            self.array2refs[arr].append(name)
        self.motifs = set()
        for arr in self.ref_array.values():
            for i in range(0, len(arr), MOTIF_LEN):
                self.motifs.add(arr[i:i + MOTIF_LEN])

    def canonical_length(self, arr):
        return len(self.LEFT) + len(arr) + len(self.RIGHT)

    def match(self, arr):
        """(exact reference names, closest same-length reference, its mismatches)."""
        exact = self.array2refs.get(arr, [])
        best, bestd = "", None
        for name, ra in self.ref_array.items():
            d = hamming(arr, ra)
            if d is None:
                continue
            if bestd is None or d < bestd:
                best, bestd = name, d
        return exact, best, bestd


def iter_fastq(paths):
    for path in paths:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            i, name = 0, None
            for line in fh:
                m = i & 3
                if m == 0:
                    name = line[1:].split()[0]
                elif m == 1:
                    yield name, line.strip().upper()
                i += 1


def iter_bam(path):
    """Every primary record's sequence. Works on unaligned/unsorted/unindexed BAM.

    Mirrors `samtools fastq -n`, which likewise drops secondary and supplementary
    records. Strand does not matter: scan() searches both orientations anyway.
    """
    import pysam
    save = pysam.set_verbosity(0)               # silence "no index" chatter on uBAM
    try:
        af = pysam.AlignmentFile(path, "rb", check_sq=False)
    finally:
        pysam.set_verbosity(save)
    with af:
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            seq = rec.query_sequence
            if seq:
                yield rec.query_name, seq.upper()


def strip_mate(name):
    return name[:-2] if name.endswith("/1") or name.endswith("/2") else name


def scan(refs, reads):
    """(Counter of array -> spanning fragments, diagnostics Counter)."""
    LA, RA, LA_rc = refs.LA, refs.RA, refs.LA_rc
    nla = len(LA)
    frag = defaultdict(set)
    stats = Counter()
    for name, seq in reads:
        stats["reads_scanned"] += 1
        i = seq.find(LA)
        if i < 0:
            if seq.find(LA_rc) >= 0:
                seq = revcomp(seq)
                i = seq.find(LA)
            else:
                if refs.RA in seq or refs.RA_rc in seq:
                    stats["right_anchor_only"] += 1
                continue
        stats["left_anchored"] += 1
        j = seq.find(RA, i + nla)
        if j < 0:
            stats["left_anchor_only"] += 1
            continue
        arr = seq[i + nla:j]
        if len(arr) % MOTIF_LEN:
            stats["out_of_phase"] += 1
            continue
        if "N" in arr:
            stats["array_with_N"] += 1
            continue
        stats["spanning_reads"] += 1
        frag[strip_mate(name)].add(arr)

    counts = Counter()
    for arrays in frag.values():
        if len(arrays) == 1:
            counts[next(iter(arrays))] += 1
        else:
            stats["discordant_fragments"] += 1
    stats["spanning_fragments"] = sum(counts.values())
    return counts, stats


def call(refs, counts, min_reads, min_freq, flank_reads=0):
    total = sum(counts.values())
    rows = []
    ordered = counts.most_common()
    passing = [(a, c) for a, c in ordered
               if c >= min_reads and total and c / total >= min_freq]
    for arr, c in passing:
        exact, closest, dist = refs.match(arr)
        n = len(arr) // MOTIF_LEN
        known = all(arr[i:i + MOTIF_LEN] in refs.motifs
                    for i in range(0, len(arr), MOTIF_LEN))
        # Relationship to the dominant variant: PCR stutter / sequencing error triage.
        # Advisory only -- it never suppresses a row.
        flag = "ok"
        if passing and arr != passing[0][0]:
            top = passing[0][0]
            if abs(len(arr) - len(top)) == MOTIF_LEN:
                flag = "stutter_candidate_of:%drep" % (len(top) // MOTIF_LEN)
            elif len(arr) == len(top) and hamming(arr, top) == 1:
                flag = "err_candidate_of:top"
        name = ",".join(exact) if exact else "NOVEL"
        ref_len = ";".join(str(len(refs.seqs[e])) for e in exact) if exact else "NA"
        rows.append({
            "length_class": refs.canonical_length(arr),
            "matched_reference": name,
            "reads": c,
            "freq": round(c / total, 4) if total else 0.0,
            "repeat_count": n,
            "ref_seq_length": ref_len,
            "motifs_all_known": "yes" if known else "no",
            "closest_reference": closest if not exact else name.split(",")[0],
            "closest_mismatches": (0 if exact else (dist if dist is not None else "NA")),
            "flag": flag,
            "total_spanning_reads": total,
            "flank_reads": flank_reads,
            "sequence": refs.LEFT + arr + refs.RIGHT,
        })
    if not rows:
        rows.append({
            "length_class": "NA", "matched_reference": "NO_CALL", "reads": 0,
            "freq": 0.0, "repeat_count": "NA", "ref_seq_length": "NA",
            "motifs_all_known": "NA", "closest_reference": "NA",
            "closest_mismatches": "NA",
            "flag": ("below_threshold" if total else
                     ("marker_present_but_no_spanning_fragment" if flank_reads
                      else "marker_absent")),
            "total_spanning_reads": total, "flank_reads": flank_reads, "sequence": "NA",
        })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Cyclospora Mt-Junction repeat-length caller (anchor-and-count)")
    ap.add_argument("--specimen", required=True)
    ap.add_argument("--out", required=True, help="output TSV (always written)")
    ap.add_argument("--fastq", nargs="+", default=[],
                    help="read FASTQ(.gz); pass R1 and R2")
    ap.add_argument("--bam", help="BAM of READS instead of FASTQ (unaligned/uBAM or the "
                                  "unmapped fraction -- NOT the marker-aligned BAM, "
                                  "which contains no junction reads)")
    ap.add_argument("--junction-ref", required=True, help="junction reference FASTA")
    ap.add_argument("--min-reads", type=int, default=3,
                    help="minimum spanning FRAGMENTS for a call; spanning depth here is "
                         "inherently 2-15 because Tn5 fragments the AT-rich array")
    ap.add_argument("--min-freq", type=float, default=0.05)
    ap.add_argument("--array-start", type=int, default=ARRAY_START)
    ap.add_argument("--anchor", type=int, default=ANCHOR)
    ap.add_argument("--diag", help="also write this diagnostics JSON")
    ap.add_argument("--emit-specimen-column", action="store_true",
                    help="prepend a 'specimen' column. OFF by default so the output is "
                         "byte-identical to the reference implementation; turn it ON "
                         "whenever the file will be aggregated by build_hds_sheet.py "
                         "under a name that does not carry the specimen id (Galaxy "
                         "datasets are named dataset_<n>.dat)")
    args = ap.parse_args(argv)

    if not args.fastq and not args.bam:
        ap.error("give --fastq or --bam")
    if args.fastq and args.bam:
        ap.error("give --fastq or --bam, not both")
    if not 14 <= args.anchor <= 39:
        ap.error("--anchor must be between 14 and 39 (see the module docstring)")
    for p in list(args.fastq) + ([args.bam] if args.bam else []) + [args.junction_ref]:
        if not os.path.exists(p):
            raise SystemExit("no such file: %s" % p)

    refs = JunctionRefs(args.junction_ref, args.array_start, args.anchor)
    for w in refs.warnings:
        sys.stderr.write("[warn] junction reference: %s\n" % w)
    sys.stderr.write("[info] %d/%d references parsed; LEFT=%d RIGHT=%d anchor=%d motifs=%d\n"
                     % (len(refs.ref_array), len(refs.seqs), len(refs.LEFT),
                        len(refs.RIGHT), args.anchor, len(refs.motifs)))

    reads = iter_bam(args.bam) if args.bam else iter_fastq(args.fastq)
    counts, stats = scan(refs, reads)
    flank_reads = stats["left_anchored"] + stats["right_anchor_only"]
    rows = call(refs, counts, args.min_reads, args.min_freq, flank_reads)

    outdir = os.path.dirname(os.path.abspath(args.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    prefix = ["specimen"] if args.emit_specimen_column else []
    with open(args.out, "w") as fh:
        fh.write("\t".join(prefix + COLUMNS) + "\n")
        for r in rows:
            cells = ([args.specimen] if prefix else []) + [str(r[c]) for c in COLUMNS]
            fh.write("\t".join(cells) + "\n")

    diag = {"specimen_id": args.specimen, "min_reads": args.min_reads,
            "min_freq": args.min_freq, **dict(stats),
            "flank_reads": flank_reads, "distinct_arrays": len(counts),
            "top_arrays": [{"repeats": len(a) // MOTIF_LEN, "reads": c}
                           for a, c in counts.most_common(10)]}
    if args.diag:
        with open(args.diag, "w") as fh:
            json.dump(diag, fh, indent=2)
    sys.stderr.write(json.dumps(diag) + "\n")
    sys.stderr.write("[info] wrote %s (%d row(s))\n" % (args.out, len(rows)))


if __name__ == "__main__":
    main()
