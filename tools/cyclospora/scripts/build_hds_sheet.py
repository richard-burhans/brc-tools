#!/usr/bin/env python3
"""Build a CDC HDS wide haplotype sheet from per-specimen caller output.

Ported from /home/anton/pyeuk-bench/lofreq_arm/build_sheets.py, generalised to consume
the per-specimen TSVs written by part_haplotype_caller.py and junction_caller.py instead
of one pre-aggregated long table, and extended to carry the junction locus (the
benchmark script could not, because its arm had no junction caller).

FORMAT -- three states, exactly as in /home/anton/pyeuk-bench/haplotype_sheet_153.txt
  header : Seq_ID <TAB> haplotype column names, sorted
  "X"    : haplotype PRESENT in that specimen
  ""     : haplotype ABSENT, the locus having been called
  an ENTIRELY EMPTY LOCUS BLOCK (every column of one locus empty for that specimen)
         : locus NOT CALLED -- amplicon dropout
There is no separate missing-data code; the empty block IS the code, and PyEuk's
dropout handling keys off it. That is why the callers must emit nothing at all for a
locus that failed its coverage gate, and why this script must never write a placeholder
into an empty block.

The one irreducible ambiguity of the format, inherited from CDC: a locus that WAS
called but had no haplotype clear the frequency/read gates also produces an empty
block, and is therefore indistinguishable from dropout. --summary reports, per
specimen, how many locus blocks came out empty so this can be audited.

Column -> locus mapping
  *_Junction_*        -> the single locus "Mt_Junction" (all length classes are one
                         locus; the Cmt<n> in the name is a length class, NOT a locus)
  *_PART_<X>_...      -> "<marker>_PART_<X>"
  anything else       -> the name with a trailing _Hap_*/_NOV_* removed

Inputs
  --calls accepts any mix of files and directories. Directories are searched (one level
  deep, plus <dir>/<specimen>/junction.tsv) for *.tsv. File type is detected from the
  header line, never from the filename:
    part-caller TSV : header starts with "specimen\tpart\t..."   (specimen in column 1)
    junction TSV    : header contains "length_class". The specimen id comes from a
                      leading "specimen" column when junction_caller.py was run with
                      --emit-specimen-column; otherwise it is inferred from the filename
                      (<specimen>.junction.tsv) or, for a file literally named
                      junction.tsv, from the parent directory name. Prefer the explicit
                      column -- Galaxy names its datasets dataset_<n>.dat.
  Rows whose matched_reference is NO_CALL contribute no column, by construction.

Usage
  build_hds_sheet.py --calls results/ --out sheet.txt
  build_hds_sheet.py --calls a.parts.tsv a.junction.tsv b.parts.tsv --out sheet.txt \\
      --specimens specimens_153.txt --summary sheet_summary.tsv
"""
import argparse
import csv
import hashlib
import os
import re
import sys

PART_RE = re.compile(r"^(.*_PART_[A-Za-z0-9]+)")


def locus_of(colname):
    """Locus a haplotype column belongs to (see module docstring)."""
    if "_Junction" in colname:
        return "Mt_Junction"
    m = PART_RE.match(colname)
    if m:
        return m.group(1)
    for sep in ("_Hap_", "_NOV_"):
        if sep in colname:
            return colname.rsplit(sep, 1)[0]
    return colname


def specimen_from_path(path):
    """Specimen id for a file that carries no specimen column."""
    base = os.path.basename(path)
    if base == "junction.tsv":
        return os.path.basename(os.path.dirname(os.path.abspath(path)))
    for suffix in (".junction.tsv", "_junction.tsv", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def gather_inputs(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            for entry in sorted(os.listdir(p)):
                full = os.path.join(p, entry)
                if os.path.isfile(full) and entry.endswith(".tsv"):
                    out.append(full)
                elif os.path.isdir(full):
                    for entry2 in sorted(os.listdir(full)):
                        if entry2.endswith(".tsv"):
                            out.append(os.path.join(full, entry2))
        elif os.path.isfile(p):
            out.append(p)
        else:
            raise SystemExit("no such file or directory: %s" % p)
    if not out:
        raise SystemExit("--calls matched no .tsv files")
    return out


def parse_file(path, include_novel):
    """(specimen, {haplotype names}) or None if the file is not caller output."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return None
    header = rows[0]
    names = set()

    if header[:2] == ["specimen", "part"]:
        idx = {h: i for i, h in enumerate(header)}
        specimens = set()
        for r in rows[1:]:
            if len(r) < len(header):
                continue
            specimens.add(r[idx["specimen"]])
            if r[idx["match"]] == "NOVEL" and not include_novel:
                continue
            names.add(r[idx["name"]])
        if len(specimens) > 1:
            raise SystemExit("%s mixes %d specimens; one specimen per file, or use a "
                             "long table split beforehand" % (path, len(specimens)))
        specimen = specimens.pop() if specimens else specimen_from_path(path)
        return specimen, names

    if header and header[0] in ("length_class", "specimen") and "length_class" in header:
        idx = {h: i for i, h in enumerate(header)}
        # junction_caller.py --emit-specimen-column puts the id in the file; otherwise
        # it has to be inferred from the path, which is why that flag exists.
        specimens = set()
        for r in rows[1:]:
            if len(r) <= idx["matched_reference"]:
                continue
            if "specimen" in idx:
                specimens.add(r[idx["specimen"]])
            ref = r[idx["matched_reference"]]
            if ref in ("NO_CALL", ""):
                continue
            if ref == "NOVEL":
                if not include_novel:
                    continue
                # A novel array has no CDC name. Label it by its length class, keeping
                # "_Junction" so it lands in the junction block, and disambiguate by
                # hashing the called sequence -- two distinct novel arrays of the SAME
                # length class must not collapse into one column.
                seq = r[idx["sequence"]] if "sequence" in idx else ""
                h = hashlib.md5(seq.encode()).hexdigest()[:6]
                names.add("Mt_Cmt%s.NOV_Junction_Hap_NOV_%s"
                          % (r[idx["length_class"]], h))
                continue
            # exact matches may be a comma-separated tie: every tied name is present
            for nm in ref.split(","):
                if nm:
                    names.add(nm)
        if len(specimens) > 1:
            raise SystemExit("%s mixes %d specimens; one specimen per file"
                             % (path, len(specimens)))
        return (specimens.pop() if specimens else specimen_from_path(path)), names

    return None


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build a CDC HDS wide haplotype sheet from per-specimen caller TSVs")
    ap.add_argument("--calls", nargs="+", required=True,
                    help="caller TSVs and/or directories containing them")
    ap.add_argument("--out", required=True, help="output sheet (TSV)")
    ap.add_argument("--specimens",
                    help="file of specimen ids, one per line: fixes the row order and "
                         "forces rows for specimens with no calls at all")
    ap.add_argument("--drop-novel", action="store_true",
                    help="keep only haplotypes that matched a named reference")
    ap.add_argument("--summary", help="optional per-specimen audit TSV")
    args = ap.parse_args(argv)

    include_novel = not args.drop_novel
    per_spec = {}
    order = []
    n_files = 0
    for path in gather_inputs(args.calls):
        parsed = parse_file(path, include_novel)
        if parsed is None:
            sys.stderr.write("[warn] not caller output, skipped: %s\n" % path)
            continue
        n_files += 1
        specimen, names = parsed
        if specimen not in per_spec:
            per_spec[specimen] = set()
            order.append(specimen)
        per_spec[specimen] |= names

    if args.specimens:
        specimens = [l.strip() for l in open(args.specimens) if l.strip()]
        unknown = [s for s in per_spec if s not in set(specimens)]
        if unknown:
            sys.stderr.write("[warn] %d specimen(s) in --calls but not in --specimens, "
                             "dropped: %s\n"
                             % (len(unknown), ", ".join(sorted(unknown)[:5])))
    else:
        specimens = sorted(order)

    cols = sorted({c for s in specimens for c in per_spec.get(s, ())})
    loci = {}
    for c in cols:
        loci.setdefault(locus_of(c), []).append(c)

    with open(args.out, "w") as fh:
        fh.write("Seq_ID\t" + "\t".join(cols) + "\n")
        for s in specimens:
            present = per_spec.get(s, set())
            fh.write(s + "\t" + "\t".join("X" if c in present else "" for c in cols) + "\n")

    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("specimen\tn_haplotypes\tn_loci_called\tn_loci_empty\tempty_loci\n")
            for s in specimens:
                present = per_spec.get(s, set())
                empty = sorted(L for L, cs in loci.items()
                               if not any(c in present for c in cs))
                fh.write("%s\t%d\t%d\t%d\t%s\n"
                         % (s, len(present), len(loci) - len(empty), len(empty),
                            ",".join(empty) if empty else "-"))

    sys.stderr.write("[info] %d file(s) -> %d specimens x %d haplotype columns "
                     "across %d loci -> %s\n"
                     % (n_files, len(specimens), len(cols), len(loci), args.out))


if __name__ == "__main__":
    main()
