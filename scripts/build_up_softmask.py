#!/usr/bin/env python3
"""Invoke the softmask workflow at increasing scope, stopping at the first tier that fails.

WHY TIERS AND NOT ONE RUN. `check_softmask_stages.py` proves each TOOL works on its own bytes. What
it cannot prove is the WIRING -- and every failure this workflow has had was wiring, not tools:
a conditional whose `hist` case demanded a dataset nothing supplied, four `state:` blocks that
imported cleanly and refused at invoke, an output connected from the wrong upstream step. A single
whole-workflow run reports the first of those and hides the rest behind it; running progressively
larger prefixes says exactly how far the graph is sound.

Each tier is the real workflow file with later steps and their outputs removed, so a tier that
passes is evidence about the committed file rather than about a hand-built copy.

    python3 scripts/build_up_softmask.py --fasta chr1.fa
    python3 scripts/build_up_softmask.py --fasta chr1.fa --from-tier 4
"""
from __future__ import annotations

import argparse
import copy
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

ROOT = pathlib.Path(__file__).resolve().parent.parent
WORKFLOW = ROOT / "workflows/softmask/softmask_udt.gxwf.yml"
UDT_DIR = ROOT / "udt"
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

TIER_CEILING = 7200          # 2 h per tier; a real chromosome is slow and single-threaded
UDTS = ("fasta_uppercase", "dustmasker_bed3", "windowmasker_bed3", "tantan_bed3",
        "lc_classify", "samtools_faidx", "fastan_gdb", "fastan_scan", "fastan_bed",
        "masking_row", "masking_header")

#: Cumulative tiers. Each names the steps present; outputs are pruned to whatever survives.
TIERS = [
    ("1 uppercase only", ["uppercase"]),
    ("2 + reference index", ["uppercase", "faidx_ref"]),
    ("3 + one masker + classify", ["uppercase", "faidx_ref", "dustmasker_bed3",
                                   "dustmasker_classify"]),
    ("4 + all four maskers", ["uppercase", "faidx_ref", "dustmasker_bed3", "dustmasker_classify",
                               "windowmasker_bed3", "windowmasker_classify",
                               "tantan_bed3", "tantan_classify",
                               "fastan_gdb", "fastan_scan", "fastan_bed"]),
    ("5 + union/sort/merge", ["uppercase", "faidx_ref", "dustmasker_bed3", "dustmasker_classify",
                              "windowmasker_bed3", "windowmasker_classify",
                              "tantan_bed3", "tantan_classify",
                              "fastan_gdb", "fastan_scan", "fastan_bed",
                              "union_cat", "sort_bed", "merge_bed"]),
    ("6 + maskfasta + faidx", ["uppercase", "faidx_ref", "dustmasker_bed3", "dustmasker_classify",
                               "windowmasker_bed3", "windowmasker_classify",
                               "tantan_bed3", "tantan_classify",
                               "fastan_gdb", "fastan_scan", "fastan_bed",
                               "union_cat", "sort_bed", "merge_bed", "maskfasta", "faidx"]),
    ("7 everything (+ genomecov, table, report)", None),   # None = keep every step
]


def connect() -> GalaxyInstance:
    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    return GalaxyInstance(url=url.rstrip("/"), key=key)


def register_all(gi: GalaxyInstance) -> dict[str, str]:
    out = {}
    for name in UDTS:
        doc = yaml.safe_load((UDT_DIR / f"{name}.gxtool.yml").read_text(encoding="utf-8"))
        created = gi.make_post_request(f"{gi.url}/unprivileged_tools",
                                       payload={"representation": doc}, params={"key": gi.key})
        out[doc["id"]] = created["uuid"]
    print(f"  registered {len(out)} UDT(s)")
    return out


def await_dataset(gi: GalaxyInstance, dataset_id: str, label: str, tries: int = 1200) -> None:
    """Block until a dataset is ready, and DIE if it errored or never settled.

    ⛔ `state in ("ok", "error")` IS NOT A READINESS TEST. It was used as one here, so an upload
    that failed -- bad format detection, truncated transfer -- became a collection element and the
    whole workflow ran on it. Exhausting the retries fell through just as silently.
    """
    for _ in range(tries):
        state = gi.datasets.show_dataset(dataset_id)["state"]
        if state == "ok":
            return
        if state in ("error", "discarded", "failed_metadata"):
            info = gi.datasets.show_dataset(dataset_id).get("misc_info") or ""
            sys.exit(f"{label} is in state {state!r} and cannot be used: {info[:200]}")
        time.sleep(5)
    sys.exit(f"{label} never became ready after {tries * 5}s (last state {state!r}).")


def prune(doc: dict, keep: list[str] | None) -> dict:
    """Keep only `keep` steps, and only outputs whose source step survives."""
    d = copy.deepcopy(doc)
    if keep is None:
        return d
    d["steps"] = {k: v for k, v in d["steps"].items() if k in keep}
    d["outputs"] = {k: v for k, v in (d.get("outputs") or {}).items()
                    if str(v.get("outputSource", "")).split("/")[0] in keep}
    return d


def render(gi: GalaxyInstance, doc: dict, uuids: dict[str, str]) -> str:
    portable = gi.workflows.import_workflow_dict(doc)
    native = gi.workflows.export_workflow_dict(portable["id"])
    gi.workflows.delete_workflow(portable["id"])
    for step in native["steps"].values():
        tid = step.get("tool_id")
        if not tid:
            continue
        if tid in uuids:
            step["tool_uuid"] = uuids[tid]
        else:
            full = gi.tools.show_tool(tid)["id"]
            if full != tid:
                step["tool_id"] = step["content_id"] = full
    return gi.workflows.import_workflow_dict(native)["id"]


def run_tier(gi: GalaxyInstance, wf_id: str, hdca: str, history: str) -> tuple[bool, str]:
    """Invoke one tier and wait it out. True only if every job it created succeeded.

    ⛔ A TIMEOUT IS A FAILURE AND SO IS AN EMPTY JOB SET. The previous version polled to a deadline,
    then computed `bad = {states in (error, paused, deleted)}` and returned `not bad` -- so a tier
    that hit the deadline with `{"running": 3}` returned TRUE and the script advanced, and a tier
    that produced no jobs at all (`states == {}`) also returned TRUE. This is the one script whose
    whole purpose is to say exactly how far the graph is sound, and it would have said "all of it"
    for a run that never finished.

    ⚠ It waits on the INVOCATION leaving `new`/`ready` too, not just on the jobs created so far:
    with a mapped-over collection Galaxy schedules jobs incrementally, so an early `{"ok": 1}` says
    nothing about the rest of the tier.
    """
    handles = gi.workflows.show_workflow(wf_id)["inputs"]
    inputs = {sid: {"src": "hdca", "id": hdca} for sid in handles}
    try:
        inv = gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history,
                                           allow_tool_state_corrections=True)
    except Exception as e:  # noqa: BLE001 - a tier reports whatever refused it, not one class
        return False, f"INVOKE REFUSED: {str(e)[:260]}"

    deadline = time.monotonic() + TIER_CEILING
    timed_out = True
    while time.monotonic() < deadline:
        detail = gi.invocations.show_invocation(inv["id"])
        states = gi.invocations.get_invocation_summary(inv["id"]).get("states", {})
        pending = {k: v for k, v in states.items() if k in ("new", "queued", "running", "paused")}
        if detail.get("state") not in SCHEDULING_IN_PROGRESS and states and not pending:
            timed_out = False
            break
        time.sleep(15)

    states = gi.invocations.get_invocation_summary(inv["id"]).get("states", {})
    url = f"{gi.base_url}/workflows/invocations/{inv['id']}"
    if timed_out:
        return False, f"TIMED OUT after {TIER_CEILING}s, jobs {states or '{}'}  {url}"
    if not states:
        return False, f"NO JOBS created by this tier  {url}"
    bad = {k: v for k, v in states.items() if k != "ok"}
    if bad:
        return False, f"FAILED jobs {states}  {url}"
    return True, f"jobs {states}  {url}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    ap.add_argument("--from-tier", type=int, default=1)
    args = ap.parse_args()

    gi = connect()
    uuids = register_all(gi)
    base = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))

    h = gi.histories.create_history(name=f"softmask build-up: {args.fasta.name}")
    print(f"  history {gi.base_url}/histories/view?id={h['id']}")
    up = gi.tools.paste_content(args.fasta.read_text(encoding="utf-8"), h["id"], file_type="fasta")
    ds = up["outputs"][0]["id"]
    await_dataset(gi, ds, f"upload {args.fasta.name!r}")
    hdca = gi.histories.create_dataset_collection(h["id"], {
        "collection_type": "list", "name": "assemblies",
        "element_identifiers": [{"name": args.fasta.stem, "src": "hda", "id": ds}]})["id"]

    for n, (label, keep) in enumerate(TIERS, start=1):
        if n < args.from_tier:
            continue
        doc = prune(base, keep)
        print(f"\nTIER {label}  ({len(doc['steps'])} steps)")
        try:
            wf_id = render(gi, doc, uuids)
        except Exception as e:  # noqa: BLE001 - same: surface the refusal, do not classify it
            print(f"  ⛔ RENDER FAILED: {str(e)[:240]}")
            return 1
        ok, detail = run_tier(gi, wf_id, hdca, h["id"])
        print(f"  {'✅' if ok else '⛔'} {detail}")
        if not ok:
            print(f"\n  Stopped at tier {n}. Everything below it is sound; this tier's wiring is not.")
            return 1
    print("\n✅ every tier ran clean, including the whole workflow")
    return 0


if __name__ == "__main__":
    sys.exit(main())
