#!/usr/bin/env python3
"""Verify a softmask invocation by recomputing its mask two independent ways.

⛔ "42 JOBS OK" IS NOT EVIDENCE THE OUTPUT IS RIGHT, AND THIS PIPELINE HAS ALREADY PROVED IT. Every
step of the committed workflow exited 0 while `maskfasta` was being fed the ORIGINAL assembly
instead of the uppercased one, so the published `softmasked_fasta` carried NCBI's pre-existing mask
unioned with this workflow's, indistinguishably -- and the invocation summary was entirely green.
An invocation summary reports scheduling, not correctness.

THE CHECK. For each strain, two numbers derived from different files by different tools, which must
agree if the pipeline is sound:

    A. bases covered by the merged union BED   (from the intervals, via bedtools merge)
    B. lowercase bases in the soft-masked FASTA (from the sequence, via bedtools maskfasta)

A == B exactly, or the mask APPLIED is not the mask COMPUTED. Under the old wiring B exceeded A by
whatever the assembly arrived carrying -- for cs10 that is tens of percent, not a rounding error.

⚠ It also asserts the uppercased input really is uppercase. Without that, A == B could hold with
both wrong.

⛔ RESOLVE BY INVOCATION, NOT BY SCANNING A HISTORY. The previous version searched a history for
collections whose NAMES matched `"Merged "` / `"bedtools MaskFastaBed"` and took the highest hid of
each, independently. Two consequences, both silent:

  * Those are ToolShed-derived DEFAULT names that nothing in this repo generates or asserts. One
    tool version bump and the script exits "is this a completed softmask run?" -- a false negative
    that reads like the run failed.
  * `build_up_softmask.py` runs all seven tiers into ONE history, and tiers 6 and 7 both produce a
    merged BED and a masked FASTA. Nothing tied the three collections to a single invocation, so a
    partially-failed tier 7 could have the union read from one tier and the FASTA from another --
    and the script would still print "the mask applied is exactly the mask computed".

An invocation exposes `output_collections` keyed by the workflow's OWN declared output names, so
`--invocation` resolves each collection exactly and cannot cross runs.

    python3 scripts/verify_softmask_outputs.py --invocation <invocation_id>
"""
from __future__ import annotations

import argparse
import os
import sys

from bioblend.galaxy import GalaxyInstance

#: Invocation states meaning Galaxy will create no further jobs.
#:
#: ⛔ `completed` IS THE COMMON ONE AND WAS MISSING HERE TOO. This file warned "invocation state is
#: 'completed', not 'scheduled' -- results below may be from an incomplete run" about a run that had
#: finished perfectly, which is the same mistake `run_softmask_udt.py` made in a place where it
#: inverted the verdict rather than just printing a false caution. Measured on this account:
#: `completed` on 22 of 25 invocations.
TERMINAL_INVOCATION_STATES = ("completed", "scheduled", "cancelled", "failed")

#: Declared output names in workflows/softmask/softmask_udt.gxwf.yml.
UPPER_OUT = "uppercased_fasta"
UNION_OUT = "mask_union"
MASKED_OUT = "softmasked_fasta"


def connect() -> GalaxyInstance:
    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    return GalaxyInstance(url=url.rstrip("/"), key=key)


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
    """Bases covered. ⚠ Assumes MERGED input; overlapping intervals would double-count."""
    total = 0
    for line in t.splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            total += int(f[2]) - int(f[1])
    return total


def resolve(gi: GalaxyInstance, invocation_id: str) -> dict[str, str]:
    """Map the three declared output names to collection ids, all from ONE invocation."""
    inv = gi.invocations.show_invocation(invocation_id)
    if inv.get("state") not in TERMINAL_INVOCATION_STATES:
        print(f"  ⚠ invocation state is {inv.get('state')!r}, which is not terminal -- "
              f"results below may be from an incomplete run.")
    cols = inv.get("output_collections") or {}
    missing = [n for n in (UPPER_OUT, UNION_OUT, MASKED_OUT) if n not in cols]
    if missing:
        sys.exit(f"invocation {invocation_id} does not declare {', '.join(missing)}. "
                 f"It exposes: {sorted(cols)}.\n"
                 f"  A workflow predating the mask_union/uppercased_fasta outputs cannot be "
                 f"verified this way -- re-run the current softmask_udt.gxwf.yml.")
    return {n: cols[n]["id"] for n in (UPPER_OUT, UNION_OUT, MASKED_OUT)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--invocation", required=True,
                    help="invocation id of a softmask run (NOT a history id -- a history can hold "
                         "several runs, and mixing their collections is exactly the bug this "
                         "argument exists to prevent)")
    args = ap.parse_args()

    gi = connect()
    ids = resolve(gi, args.invocation)
    up_e = elements(gi, ids[UPPER_OUT])
    mg_e = elements(gi, ids[UNION_OUT])
    mk_e = elements(gi, ids[MASKED_OUT])

    strains = sorted(set(up_e) & set(mg_e) & set(mk_e))
    if not strains:
        sys.exit("no strain appears in all three collections; element identifiers do not line up.")
    for label, e in ((UPPER_OUT, up_e), (UNION_OUT, mg_e), (MASKED_OUT, mk_e)):
        extra = sorted(set(e) - set(strains))
        if extra:
            print(f"  ⚠ {label} also holds {extra}, absent from another collection -- NOT checked.")

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
        print(f"      merged union BED : {covered:,} nt ({covered / res:.2%})")
        print(f"      masked FASTA     : {low:,} lowercase ({low / res:.2%})")
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
