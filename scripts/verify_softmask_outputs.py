#!/usr/bin/env python3
"""Verify a completed softmask run by recomputing its mask two independent ways.

⛔ "34 JOBS OK" IS NOT EVIDENCE THE OUTPUT IS RIGHT, AND THIS PIPELINE HAS ALREADY PROVED IT. Every
step of the committed workflow exited 0 while `maskfasta` was being fed the ORIGINAL assembly
instead of the uppercased one, so the published `softmasked_fasta` carried NCBI's pre-existing mask
unioned with this workflow's, indistinguishably -- and the invocation summary was entirely green.
An invocation summary reports scheduling, not correctness.

THE CHECK. For each strain, take two numbers that are derived from different files by different
tools and must agree if the pipeline is sound:

    A. bases covered by the merged union BED   (from the intervals, via bedtools)
    B. lowercase bases in the soft-masked FASTA (from the sequence, via bedtools maskfasta)

A == B exactly, or the mask that was APPLIED is not the mask that was COMPUTED. Under the old
wiring B exceeded A by however much the assembly arrived carrying -- for cs10 that is tens of
percent, not a rounding difference.

⚠ It also asserts the uppercased input really is uppercase, because if it were not, A == B could
hold while both were wrong.

    python3 scripts/verify_softmask_outputs.py --history <history_id>
"""
from __future__ import annotations

import argparse
import os
import sys

from bioblend.galaxy import GalaxyInstance

#: Collection names produced by workflows/softmask/softmask_udt.gxwf.yml.
UPPER = "Uppercased FASTA"
MERGED_PREFIX = "Merged "
MASKED_PREFIX = "bedtools MaskFastaBed"


def connect() -> GalaxyInstance:
    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    return GalaxyInstance(url=url.rstrip("/"), key=key)


def collections(gi: GalaxyInstance, history: str) -> list[dict]:
    return [i for i in gi.histories.show_history(history, contents=True)
            if i.get("history_content_type") == "dataset_collection"]


def elements(gi: GalaxyInstance, collection_id: str) -> dict[str, str]:
    d = gi.dataset_collections.show_dataset_collection(collection_id)
    return {e["element_identifier"]: e["object"]["id"] for e in d["elements"]}


def text(gi: GalaxyInstance, dataset_id: str) -> str:
    b = gi.datasets.download_dataset(dataset_id, use_default_filename=False)
    return b.decode("utf-8", "replace") if isinstance(b, bytes) else str(b)


def fasta_stats(t: str) -> tuple[int, int]:
    """(residues, lowercase residues)."""
    res = low = 0
    for line in t.splitlines():
        if line.startswith(">"):
            continue
        res += len(line)
        low += sum(1 for c in line if "a" <= c <= "z")
    return res, low


def bed_span(t: str) -> int:
    """Bases covered. Assumes MERGED input -- overlapping intervals would double-count."""
    total = 0
    for line in t.splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            total += int(f[2]) - int(f[1])
    return total


def pick(cols: list[dict], predicate) -> dict | None:
    hits = [c for c in cols if predicate(c["name"])]
    return max(hits, key=lambda c: c["hid"]) if hits else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--history", required=True, help="history id of a completed softmask run")
    args = ap.parse_args()

    gi = connect()
    cols = collections(gi, args.history)
    upper = pick(cols, lambda n: n == UPPER)
    merged = pick(cols, lambda n: n.startswith(MERGED_PREFIX))
    masked = pick(cols, lambda n: n.startswith(MASKED_PREFIX))
    missing = [lbl for lbl, c in (("uppercased", upper), ("merged union", merged),
                                  ("masked FASTA", masked)) if c is None]
    if missing:
        sys.exit(f"history {args.history} has no {', '.join(missing)} collection -- "
                 f"is this a completed softmask run?")

    up_e, mg_e, mk_e = (elements(gi, c["id"]) for c in (upper, merged, masked))
    strains = sorted(set(up_e) & set(mg_e) & set(mk_e))
    if not strains:
        sys.exit("no strain appears in all three collections; element identifiers do not line up.")
    for label, e in (("uppercased", up_e), ("merged", mg_e), ("masked", mk_e)):
        if set(e) != set(strains):
            print(f"  ⚠ {label} collection covers {sorted(set(e) ^ set(strains))} not in the "
                  f"intersection; those strains are NOT checked below.")

    print("  union BED coverage vs masked-FASTA lowercase -- derived independently\n")
    failures = 0
    for s in strains:
        covered = bed_span(text(gi, mg_e[s]))
        res, low = fasta_stats(text(gi, mk_e[s]))
        ures, ulow = fasta_stats(text(gi, up_e[s]))
        ok_upper = (ulow == 0 and ures == res)
        ok_mask = (covered == low)
        failures += (not ok_upper) + (not ok_mask)
        print(f"    {s}")
        print(f"      uppercased input : {ures:,} nt, {ulow} lowercase"
              f"        {'ok' if ok_upper else '⛔ NOT UPPERCASE / LENGTH CHANGED'}")
        print(f"      merged union BED : {covered:,} nt ({covered/res:.2%})")
        print(f"      masked FASTA     : {low:,} lowercase ({low/res:.2%})")
        print(f"      delta            : {low - covered:+,}"
              f"        {'MATCH' if ok_mask else '⛔ applied mask != computed mask'}")
    print()
    if failures:
        print(f"  ⛔ {failures} check(s) FAILED across {len(strains)} strain(s)")
        return 1
    print(f"  ✅ {len(strains)} strain(s): the mask applied is exactly the mask computed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
