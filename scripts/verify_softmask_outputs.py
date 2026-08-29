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

#: Invocation states in which Galaxy may still create jobs. Everything else means scheduling is
#: finished, whatever the jobs are doing.
#:
#: ⛔ THE VALUES COME FROM GALAXY'S OWN ENUM, not from watching behaviour. lib/galaxy/schema/
#: invocation.py::InvocationState documents each one, and two of these were missing when this set
#: was written from observation alone:
#:     new                       "Brand new workflow invocation"
#:     ready                     "Workflow ready for another iteration of scheduling."
#:     requires_materialization  "an otherwise NEW or READY workflow that requires inputs to be
#:                                materialized (undeferred)"
#:     cancelling                "invocation scheduler will cancel job in next iteration."
#:
#: ⚠ AND THE SAME FILE SETTLES WHY `completed` CANNOT BE USED AS THE WAIT CONDITION. It defines
#: `scheduled` as "Workflow has been scheduled" and `completed` as "All jobs have reached terminal
#: states" -- so `completed` is the state one WANTS, and it is nevertheless unreliable: measured
#: over 60 invocations on this account, 50 `completed` and 9 `scheduled`, interleaved across the
#: whole timeline, with structurally identical runs landing differently and one sitting `scheduled`
#: for 6.8 days with all ten jobs `ok` and `update_time` frozen at creation. Galaxy records the
#: transition in a separate `workflow_invocation_completion` row (model/__init__.py), so an
#: invocation whose completion hook never fires stays `scheduled` forever. Wait on the JOBS.
SCHEDULING_IN_PROGRESS = ("new", "ready", "requires_materialization", "cancelling")

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


def bed_intervals(t: str) -> set[tuple[str, int, int]]:
    """Parse a BED into a set of (chrom, start, end).

    ⛔ POSITIONS, NOT A TOTAL. This used to return only the summed width, and the check compared
    that against the masked FASTA's lowercase COUNT. Equal totals do not mean the same bases: an
    adversarial pass showed the old form certifying four wrong pipelines as correct -- intervals at
    0-10 against lowercase at 50-60, an off-by-one, a mask applied to the wrong chromosome, and an
    empty BED against an entirely unmasked FASTA. All four printed "delta +0 MATCH".

    ⚠ Track/browser/comment lines are REFUSED rather than skipped. Silently ignoring them would
    under-count, which is the same class of error one level down.
    """
    out: set[tuple[str, int, int]] = set()
    for n, line in enumerate(t.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith(("#", "track", "browser")):
            sys.exit(f"BED line {n} is a {line.split()[0]!r} line; this verifier expects a plain "
                     f"merged BED and will not silently skip it: {line[:60]!r}")
        f = line.split("\t")
        if len(f) < 3:
            sys.exit(f"BED line {n} has {len(f)} fields, need at least 3: {line[:60]!r}")
        try:
            out.add((f[0], int(f[1]), int(f[2])))
        except ValueError:
            sys.exit(f"BED line {n} has non-integer coordinates: {line[:60]!r}")
    return out


def lowercase_runs(t: str) -> set[tuple[str, int, int]]:
    """Maximal lowercase runs of a FASTA, as (chrom, start, end) in BED half-open coordinates.

    ⚠ The chrom is the FIRST WHITESPACE TOKEN of the header, because that is what bedtools writes
    into the BED it produced. `bedtools maskfasta` also TRUNCATES the header to exactly that token,
    so the masked FASTA and its own intervals agree by construction -- but the uppercased input does
    not have truncated headers, and comparing the two naively would mismatch on description text.
    """
    runs: set[tuple[str, int, int]] = set()
    chrom, pos, start = None, 0, None
    for line in t.splitlines():
        if line.startswith(">"):
            if chrom is not None and start is not None:
                runs.add((chrom, start, pos))
            chrom, pos, start = line[1:].split()[0] if len(line) > 1 else "", 0, None
            continue
        for ch in line:
            if "a" <= ch <= "z":
                if start is None:
                    start = pos
            elif start is not None:
                runs.add((chrom, start, pos))
                start = None
            pos += 1
    if chrom is not None and start is not None:
        runs.add((chrom, start, pos))
    return runs


def span(intervals: set[tuple[str, int, int]]) -> int:
    return sum(e - s for _, s, e in intervals)


def resolve(gi: GalaxyInstance, invocation_id: str) -> dict[str, str]:
    """Map the three declared output names to collection ids, all from ONE invocation."""
    inv = gi.invocations.show_invocation(invocation_id)
    state = inv.get("state")
    # ⛔ TERMINAL IS NOT THE SAME AS SOUND. The tuple means "Galaxy will create no further jobs",
    # which is the right wait condition elsewhere -- but `failed` and `cancelled` are terminal too,
    # and certifying a cancelled run's partly-populated collections would be worse than useless.
    if state in ("failed", "cancelled"):
        sys.exit(f"invocation {invocation_id} is in state {state!r}; refusing to verify it.")
    if state in SCHEDULING_IN_PROGRESS:
        print(f"  ⚠ invocation state is {state!r}: Galaxy may still be creating jobs, "
              f"so the collections below may be incomplete.")
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
        union = bed_intervals(text(gi, mg_e[s]))
        masked_text = text(gi, mk_e[s])
        applied = lowercase_runs(masked_text)
        res, low = fasta_stats(masked_text)
        ures, ulow = fasta_stats(text(gi, up_e[s]))
        if res == 0:
            print(f"    {s}\n      ⛔ masked FASTA is EMPTY")
            failures += 1
            continue

        ok_upper = (ulow == 0 and ures == res)
        # ⛔ AN EMPTY MASK IS A FAILURE, NOT A MATCH. Nothing masked against nothing computed is
        # 0 == 0, and the old total-based check printed ✅ for it -- certifying a run in which the
        # mask was never applied at all.
        ok_nonempty = bool(union)
        only_computed = union - applied
        only_applied = applied - union
        ok_mask = ok_nonempty and not only_computed and not only_applied
        failures += (not ok_upper) + (not ok_mask)

        print(f"    {s}")
        print(f"      uppercased input : {ures:,} nt, {ulow} lowercase"
              f"        {'ok' if ok_upper else '⛔ NOT UPPERCASE / LENGTH CHANGED'}")
        print(f"      merged union BED : {span(union):,} nt in {len(union):,} intervals "
              f"({span(union) / res:.2%})")
        print(f"      masked FASTA     : {low:,} lowercase in {len(applied):,} runs "
              f"({low / res:.2%})")
        if not ok_nonempty:
            print("      ⛔ the union is EMPTY -- nothing was masked, which is not a pass")
        elif ok_mask:
            print("      positions        : IDENTICAL interval-for-interval")
        else:
            print(f"      ⛔ computed-but-not-applied: {len(only_computed):,} intervals, "
                  f"{span(only_computed):,} nt")
            print(f"      ⛔ applied-but-not-computed: {len(only_applied):,} intervals, "
                  f"{span(only_applied):,} nt")
            for iv in sorted(only_computed)[:3]:
                print(f"          only in BED  : {iv}")
            for iv in sorted(only_applied)[:3]:
                print(f"          only in FASTA: {iv}")
    print()
    if failures:
        print(f"  ⛔ {failures} check(s) FAILED across {len(strains)} strain(s)")
        return 1
    print(f"  ✅ {len(strains)} strain(s): the mask applied is exactly the mask computed, "
          f"interval for interval")
    return 0


if __name__ == "__main__":
    sys.exit(main())
