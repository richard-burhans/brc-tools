#!/usr/bin/env python3
"""Sketch every assembly in a collection, then compare -- labelling by element identifier.

The identifiers cannot be read from the collection inside a job (see the module docstring of
scripts/build_inventory_udts.py); they are supplied as a file, in collection order.
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--rendered", required=True, help="file holding the rendered collection input")
ap.add_argument("--ids", required=True, help="element identifiers, one per line, in collection order")
ap.add_argument("--ksize", default="31")
ap.add_argument("--scaled", default="1000")
a = ap.parse_args()

# ⛔ THE KEYS ARE QUOTED. Galaxy renders a collection input as JSON -- {"path": "/..."} -- not as a
# JavaScript object literal, so a pattern written for a bare path key matches NOTHING and this
# script reported "0 assemblies" against a perfectly good render: the heredoc held six well-formed
# File records and every one was missed, which is why WF-A could never sketch anything. Parse it as
# JSON, and fall back to the pattern only if that fails, so a future change in render shape
# degrades loudly instead of silently finding zero.
_raw = pathlib.Path(a.rendered).read_text()
_json_ok = True
try:
    _recs = json.loads(_raw[_raw.index("["):_raw.rindex("]") + 1])
    paths = [r["path"] for r in _recs if isinstance(r, dict) and r.get("path")]
except Exception:                                   # noqa: BLE001 -- see below
    # ⚠ THE FALLBACK IS A DIAGNOSTIC, NOT A SECOND PARSER. Once it scraped even ONE path the
    # render fault below stopped firing and the count check took over -- so a render shaped
    # {"class":"File","path":...} instead of a list, which is exactly the "future change in render
    # shape" this fallback anticipates, was reported as "1 assemblies but 2 identifiers. The
    # identifier file must come from the SAME collection" and sent the operator to inspect the one
    # input that was correct. Remember that the JSON parse failed, and say so.
    _json_ok = False
    paths = re.findall(r'"?path"?\s*:\s*"?([^",}\s]+)', _raw)
# ⚠ utf-8-SIG, NOT utf-8, and the same identifier checks the two sibling helpers already make.
# A BOM survives `.strip()` and lands in the first name, so the matrix comes out labelled
# `﻿cs10` -- which then joins against nothing in the self-pair or relabel files built from the
# SAME file. Measured. A tab is rejected there and was accepted here; the CSV header took it.
ids = [x.strip() for x in
       pathlib.Path(a.ids).read_text(encoding="utf-8-sig").splitlines() if x.strip()]
_bad = [i for i in ids if "	" in i]
if _bad:
    sys.exit(f"identifier(s) contain a tab, which lands in the similarity matrix header and breaks "
             f"every column contract downstream: {_bad[:3]}")

# ⛔ THE ASSERTION IS THE POINT. Pairing by position is correct only while the two lists describe
# the same collection; if they do not, every row of the matrix is mislabelled and nothing
# downstream can tell. Refuse instead.
# ⛔ DIAGNOSE THE RENDER BEFORE BLAMING THE IDENTIFIERS. This comparison used to run first, so an
# unparseable or empty rendered block reported "0 assemblies but N identifiers -- the identifier file
# must come from the SAME collection", pointing the operator at the one input that was correct. A
# render this cannot read is a DIFFERENT fault and says so.
if not _json_ok or not paths:
    sys.exit(f"could not parse the rendered collection input as JSON; the fallback pattern "
             f"scraped {len(paths)} path(s). That is a "
             "RENDER problem, not an identifier problem -- the identifier file is not implicated. "
             "The block should be a JSON array of File records; check what the tool actually "
             "received before changing anything about the identifiers.")
# ⛔ DUPLICATE IDENTIFIERS SILENTLY DESTROY A ROW. Each signature is staged at `stage/{name}.sig`,
# so two elements sharing a name overwrite one file -- and the count assertion below still passes,
# because the COUNTS match. The result is a fully populated matrix in which one genome does not
# appear and another appears twice, which is exactly the "runs and mislabels its output" failure
# this script's assertions exist to prevent.
if len(set(ids)) != len(ids):
    _dupes = sorted({i for i in ids if ids.count(i) > 1})
    sys.exit(f"duplicate element identifier(s) {_dupes[:3]}: a repeat would put two genomes in one "
             f"row of the matrix. Counts alone cannot detect this.")
# ⚠ AND THESE NAMES ARE STILL REFUSED, THOUGH THE REASON HAS CHANGED. While signatures were staged
# at `stage/{name}.sig` this was a data-loss guard: `gA` and `./gA` are distinct strings naming the
# SAME file, so one sketch overwrote the other, the matrix came out `./gA,./gA` with an
# off-diagonal 1.0 that is a self-comparison, and a genome was absent -- exit 0. `../escaped` wrote
# the signature outside the job directory. Staging by INDEX (below) removed both hazards outright.
# What is left is narrower and worth keeping anyway: such a name still travels into `--name` and
# becomes a COLUMN LABEL in similarity.csv, where a leading dash or an embedded path reads as a
# malformed strain and joins against nothing downstream. Refusing early says so; it is no longer
# the difference between a correct matrix and a corrupt one.
_bad = [i for i in ids if "/" in i or i in (".", "..") or i.startswith("-")]
if _bad:
    sys.exit(f"element identifier(s) {_bad[:3]} contain a path separator, are a directory alias, or "
             f"begin with a dash. Galaxy allows them; this tool will not stage them, because such a "
             f"name can silently overwrite another element's signature or write outside the job.")
if len(paths) != len(ids):
    sys.exit(f"refusing to guess: {len(paths)} assemblies but {len(ids)} identifiers. "
             f"The identifier file must come from the SAME collection, via the IUC "
             f"collection_element_identifiers tool.")
subprocess.run(["mkdir", "-p", "stage"], check=True)
for _i, (path, name) in enumerate(zip(paths, ids, strict=True)):
    sig = f"stage/{_i:04d}.sig"
    subprocess.run(["sourmash", "sketch", "dna", "-p", f"k={a.ksize},scaled={a.scaled}",
                    "--name", name, "-o", sig, path], check=True)
    # ⛔ AN EMPTY SKETCH IS INDISTINGUISHABLE FROM AN UNRELATED GENOME, AND sourmash SAYS NOTHING.
    # At --scaled 1000 a small enough sequence contributes NO hashes: measured, a 500 bp assembly
    # and 180 kb of N both sketched to 0 mins, and the similarity.csv that came out was
    # BYTE-IDENTICAL to the run where that element was a genuine 200 kb unrelated genome -- row
    # `0.0, 1.0, 0.0`, exit 0. So a panel member that failed upstream, or any small-genome member
    # (an organelle, a plasmid, an apicoplast), reads as "shares nothing with anyone" and this
    # matrix feeds WF-I's fold order. The count assertion above cannot see it; only the hashes can.
    _n = sum(len(s.get("mins") or []) for rec in json.loads(pathlib.Path(sig).read_text())
             for s in (rec.get("signatures") or []))
    # ⚠ AND THE BOUNDARY IS NOT ZERO. The harm this describes -- a row that reads "shares nothing
    # with anyone" -- is a property of a SMALL sketch, not of an empty one. Measured at
    # scaled=1000, k=31: a 360 nt viroid gives 1 hash, a 2 kb plasmid 3, a 5 kb virus 7, and each
    # of their similarity rows is BYTE-IDENTICAL to a genuinely unrelated 200 kb genome's. Worse,
    # a 30 kb truncation of a 3 Mb assembly -- 19 hashes, 100% contained in it -- scored 0.0064
    # against its own parent, indistinguishable from 0.0 in a heatmap and in WF-I's fold order.
    # Refusing at 0 and waving through 1 draws the line in the one place it does not belong.
    _floor = max(20, int(a.scaled) // 50)
    if _n == 0:
        sys.exit(f"{name} sketched to ZERO hashes at scaled={a.scaled}, so it can only appear in "
                 f"the matrix as a genome sharing nothing with anyone -- which is not a failure "
                 f"the run would otherwise report. Lower --scaled, or drop the element. ⚠ It is "
                 f"N-masking that empties a sketch, NOT soft-masking: a 200 kb genome in "
                 f"all-lowercase gives 195 hashes, measured.")
    if _n < _floor:
        print(f"⚠ {name} sketched only {_n} hash(es) at scaled={a.scaled} (floor {_floor}). Its "
              f"row will look like an unrelated genome's whatever it actually is -- a 5 kb genome "
              f"at 7 hashes is byte-identical in the matrix to something that shares nothing. "
              f"Treat its similarities as unmeasured, or lower --scaled for the whole panel.",
              file=sys.stderr)
    print(f"sketched {name} <- {path}  ({_n} hashes)", file=sys.stderr)

sigs = [f"stage/{i:04d}.sig" for i in range(len(ids))]
subprocess.run(["sourmash", "compare", "--ksize", a.ksize, "-o", "cmp", "--csv", "similarity.csv",
                *sigs], check=True)
subprocess.run(["sourmash", "plot", "--labels", "cmp"], check=True)
