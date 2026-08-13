#!/usr/bin/env python3
"""Per-PART haplotype caller for Cyclospora cayetanensis amplicon data.

Ported from the two-stage benchmark implementation
  /home/anton/pyeuk-bench/lofreq_arm/extract_reads.py   (per-read, per-refpos extraction)
  /home/anton/pyeuk-bench/lofreq_arm/call_haplotypes.py (caller "A", direct read-collapse)
into a single argparse CLI with no hardcoded paths, so a Galaxy wrapper can call it.

Only caller A is ported. Caller B (the LoFreq-informed variant that overwrites
non-variant positions with the reference base) was the benchmark's comparison arm and
needs a per-specimen VCF; it is deliberately out of scope here.

--markers is therefore NOT used to build any called sequence -- caller A reads the
haplotype straight off the alignment. It is used to validate that every BED interval
exists in, and lies inside, the reference the BAM was mapped against, which turns a
silently-empty result (wrong BED, or a BED for a different marker set) into an error.

Method
------
There is NO assembly. Each amplicon in the marker FASTA is cut into ~100 bp PART
intervals by the BED; reads are long enough to span a PART end to end, so the
haplotype of a PART is simply the aligned subsequence of a spanning read, and
identical subsequences are collapsed and counted.

A read counts as SPANNING a PART [start, end) only if
    read.reference_start <= start  and  read.reference_end >= end
i.e. the alignment covers every reference position of the interval. Secondary,
supplementary and unmapped records are skipped. (The upstream BAM is expected to be
already filtered to MAPQ >= 20 proper pairs, as produced by the mapping stage; this
script does not re-apply those filters, exactly like the reference implementation.)

Per spanning read the interval is rendered as one string per reference position:
  * an aligned base (CIGAR M/=/X)              -> that base
  * a deleted / skipped position (CIGAR D/N)   -> the empty string
  * an insertion (CIGAR I)                     -> appended, lowercased, to the string
                                                  of the preceding reference position
  * soft/hard clips and pads                   -> consume query/nothing, contribute nothing
The haplotype sequence is the concatenation of those per-position strings, uppercased,
so indels inside a PART change the called sequence length. A read whose alignment
leaves any position of the interval unset is discarded.

Gates (PRE-REGISTERED in
/home/anton/pyeuk-bench/lofreq_arm/THRESHOLD_PREREGISTRATION.md -- do not tune)
  --min-span  50    a PART is CALLED only if >= 50 reads fully span it. Below that the
                    PART emits no rows at all, which downstream is read as "not called"
                    (an entirely empty locus block in the HDS sheet), NOT as "absent".
  --min-freq  0.05  a haplotype is kept at >= 5% of the spanning reads, AND
  --min-reads 10    at >= 10 supporting reads. Both must hold.
Frequency denominator is the spanning-read count of the PART, not the sum of the
haplotypes that survived the gates.

Naming
------
Exact string match against the known-haplotype FASTA (haplotypes78.fa). Records there
are named <PART>_Hap_<k>, so the PART a record belongs to is its name with the trailing
"_Hap_<k>" removed; records without "_PART_" in the name are ignored (the junction
references live in a different file). No match -> <PART>_NOV_<md5(seq)[:6]>, and the
nearest known haplotype of the same PART plus its Levenshtein distance are reported for
triage only -- the distance never affects the call.

Usage
-----
  part_haplotype_caller.py --bam SPEC.filtered.bam --markers markers.fa \\
      --parts parts.bed --known-haplotypes haplotypes78.fa \\
      --specimen SPEC --out SPEC.parts.tsv

Output: one TSV row per kept haplotype, columns
  specimen part span reads freq name match nearest edit seq
"""
import argparse
import hashlib
import os
import sys

import pysam

COLUMNS = ["specimen", "part", "span", "reads", "freq", "name", "match",
           "nearest", "edit", "seq"]


def read_fasta(path):
    """FASTA -> {full header line (minus '>'): uppercased sequence}."""
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                out[name] = []
            elif name is not None:
                out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}


def load_parts(path):
    """BED -> [(part_name, chrom, start, end)] in file order."""
    parts = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            parts.append((f[3], f[0], int(f[1]), int(f[2])))
    if not parts:
        raise SystemExit("no 4-column intervals in %s" % path)
    return parts


def load_known(path):
    """haplotypes78.fa -> ({part: {seq: name}}, {part: [(name, seq)]})."""
    by_seq, by_part = {}, {}
    for name, seq in read_fasta(path).items():
        if "_PART_" not in name:
            continue
        part = name.rsplit("_Hap_", 1)[0]
        by_seq.setdefault(part, {})[seq] = name
        by_part.setdefault(part, []).append((name, seq))
    return by_seq, by_part


def levenshtein(a, b):
    """Edit distance; Hamming shortcut when the lengths already agree."""
    if a == b:
        return 0
    if len(a) == len(b):
        return sum(1 for x, y in zip(a, b) if x != y)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def read_part_repr(read, start, end):
    """Per-reference-position strings for [start, end), or None if not fully spanned."""
    if read.reference_start > start or read.reference_end < end:
        return None
    seq = read.query_sequence
    if seq is None or read.cigartuples is None:
        return None
    n = end - start
    cells = [None] * n
    qpos = 0
    rpos = read.reference_start
    last_idx = None
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):                    # M / = / X
            for k in range(ln):
                i = rpos + k - start
                if 0 <= i < n:
                    cells[i] = seq[qpos + k]
                    last_idx = i
            qpos += ln
            rpos += ln
        elif op in (2, 3):                     # D / N
            for k in range(ln):
                i = rpos + k - start
                if 0 <= i < n:
                    cells[i] = ""
                    last_idx = i
            rpos += ln
        elif op == 1:                          # I -- attach to previous refpos
            if last_idx is not None and 0 <= last_idx < n:
                cells[last_idx] = cells[last_idx] + seq[qpos:qpos + ln].lower()
            qpos += ln
        elif op == 4:                          # S
            qpos += ln
        elif op in (5, 6):                     # H / P
            pass
        else:
            return None
    if any(c is None for c in cells):
        return None
    return cells


def collect(bam_path, parts):
    """{part: (span, {haplotype_sequence: count})} over all spanning reads."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for name, chrom, start, end in parts:
            obs, span = {}, 0
            try:
                it = bam.fetch(chrom, start, end)
            except ValueError:
                # contig absent from this BAM's header -> zero coverage, not an error
                out[name] = (0, {})
                continue
            for read in it:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                cells = read_part_repr(read, start, end)
                if cells is None:
                    continue
                span += 1
                s = "".join(cells).upper()
                obs[s] = obs.get(s, 0) + 1
            out[name] = (span, obs)
    return out


def call(observed, parts, known_by_seq, known_by_part,
         min_span, min_freq, min_reads, specimen):
    rows = []
    for part, _chrom, _start, _end in parts:
        span, obs = observed[part]
        if span < min_span:
            continue                            # locus NOT CALLED
        kn = known_by_seq.get(part, {})
        for seq, cnt in sorted(obs.items(), key=lambda x: (-x[1], x[0])):
            freq = cnt / span
            if freq < min_freq or cnt < min_reads:
                continue
            if seq in kn:
                name, match, nearest, edit = kn[seq], "EXACT", kn[seq], 0
            else:
                best, bd = "", 10 ** 9
                for nm, ks in known_by_part.get(part, []):
                    d = levenshtein(seq, ks)
                    if d < bd:
                        bd, best = d, nm
                nearest, edit = best, (bd if best else -1)
                name = "%s_NOV_%s" % (part, hashlib.md5(seq.encode()).hexdigest()[:6])
                match = "NOVEL"
            rows.append([specimen, part, span, cnt, "%.4f" % freq,
                         name, match, nearest, edit, seq])
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Cyclospora per-PART amplicon haplotype caller (direct read-collapse)")
    ap.add_argument("--bam", required=True,
                    help="coordinate-sorted, indexed BAM aligned to --markers "
                         "(expected pre-filtered to MAPQ>=20 proper pairs)")
    ap.add_argument("--markers", required=True, help="marker/amplicon FASTA")
    ap.add_argument("--parts", required=True,
                    help="BED4 of PART intervals; column 4 is the PART name")
    ap.add_argument("--known-haplotypes", required=True,
                    help="FASTA of named reference haplotypes (haplotypes78.fa)")
    ap.add_argument("--min-span", type=int, default=50,
                    help="spanning reads required to CALL a PART (pre-registered: 50)")
    ap.add_argument("--min-freq", type=float, default=0.05,
                    help="haplotype frequency floor (pre-registered: 0.05)")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="haplotype read floor (pre-registered: 10)")
    ap.add_argument("--specimen", required=True, help="specimen id written to column 1")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--summary", help="optional TSV of per-PART span / called status")
    args = ap.parse_args(argv)

    for p in (args.bam, args.markers, args.parts, args.known_haplotypes):
        if not os.path.exists(p):
            raise SystemExit("no such file: %s" % p)

    parts = load_parts(args.parts)
    markers = read_fasta(args.markers)
    missing = sorted({c for _n, c, _s, _e in parts} - set(markers))
    if missing:
        raise SystemExit("BED references contigs absent from --markers: %s"
                         % ", ".join(missing))
    for name, chrom, start, end in parts:
        if end > len(markers[chrom]) or start < 0 or start >= end:
            raise SystemExit("interval %s (%s:%d-%d) is outside %s (%d bp)"
                             % (name, chrom, start, end, chrom, len(markers[chrom])))

    known_by_seq, known_by_part = load_known(args.known_haplotypes)
    if not known_by_part:
        sys.stderr.write("[warn] no _PART_ records in %s; every call will be NOVEL\n"
                         % args.known_haplotypes)

    observed = collect(args.bam, parts)
    rows = call(observed, parts, known_by_seq, known_by_part,
                args.min_span, args.min_freq, args.min_reads, args.specimen)

    with open(args.out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("specimen\tpart\tspan\tcalled\tn_haplotypes\n")
            kept = {}
            for r in rows:
                kept[r[1]] = kept.get(r[1], 0) + 1
            for part, _c, _s, _e in parts:
                span = observed[part][0]
                fh.write("%s\t%s\t%d\t%s\t%d\n"
                         % (args.specimen, part, span,
                            "yes" if span >= args.min_span else "no", kept.get(part, 0)))

    called = sum(1 for p, _c, _s, _e in parts if observed[p][0] >= args.min_span)
    sys.stderr.write("[info] %s: %d/%d PARTs called, %d haplotype rows -> %s\n"
                     % (args.specimen, called, len(parts), len(rows), args.out))


if __name__ == "__main__":
    main()
