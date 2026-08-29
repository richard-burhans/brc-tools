#!/usr/bin/env python3
"""Register the softmask UDTs, render an instance-resolved workflow, and invoke it.

⛔ THE COMMITTED WORKFLOW CANNOT BE INVOKED AS WRITTEN, AND IMPORTING IT SAYS OTHERWISE. Measured on
usegalaxy.org 2026-08-28: a gxformat2 file whose steps name UDTs imports with every node present and
no errors, and then `POST .../invocations` refuses with

    the following required tools are not installed: brc-dustmasker-bed3 (version None), ...

for TWO independent reasons:

  1. A UDT IS ADDRESSED BY UUID. Its plain id is not in the toolbox that the invocation's
     "required tools" gate consults -- `GET /api/tools` cannot see it at all, while
     `GET /api/unprivileged_tools` lists it happily.
  2. AN UNVERSIONED TOOLSHED ID DOES NOT RESOLVE EITHER. `bedtools_sortbed`, `bedtools_mergebed`
     and `bedtools_maskfastabed` are all short forms; the tool API resolves them, workflow
     invocation does not. `cat1` is fine because it is a built-in with no toolshed path.

⚠ THE ROUND TRIP IS THE MECHANISM, NOT A DETOUR. gxformat2 has no step-level `tool_uuid` -- posting
one returns HTTP 500 -- but Galaxy's NATIVE workflow format does. So: import the portable file
(which resolves the connections and preserves every `state:` block), export the native
representation, resolve the identities there, and re-import THAT. Hand-authoring the native dict
instead does NOT work; it 500s. The technique is rung-intel's, in `galaxy/run_on_galaxy.py`.

The committed workflow stays portable: UUIDs are per-account facts this repository must not hold,
and pinning `+galaxyN` revisions is what the short-id convention exists to avoid.

    python3 scripts/run_softmask_udt.py --register-only
    python3 scripts/run_softmask_udt.py --fasta a.fa --fasta b.fa
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance
from softmask_lib import ROOT, SCHEDULING_IN_PROGRESS, WORKFLOW, await_dataset, connect, register_all

#: Only the tools this workflow uses. `env_probe` lives in udt/ too and is a diagnostic.
POLL_SECONDS = 20
POLL_CEILING = 21600          # 6 h: a 101 Mb chromosome, single-threaded, two-pass windowmasker


def render_and_import(gi: GalaxyInstance, uuids: dict[str, str], work: pathlib.Path) -> str:
    """Portable gxformat2 -> native -> identities resolved -> re-imported. Returns workflow id."""
    portable = gi.workflows.import_workflow_dict(yaml.safe_load(WORKFLOW.read_text(encoding="utf-8")))
    native = gi.workflows.export_workflow_dict(portable["id"])
    gi.workflows.delete_workflow(portable["id"])  # a scaffold, not an artifact

    notes, resolved = [], 0
    for step in native["steps"].values():
        tool_id = step.get("tool_id")
        if not tool_id:
            continue
        if tool_id in uuids:
            step["tool_uuid"] = uuids[tool_id]
            notes.append(f"    {tool_id:28} -> uuid {uuids[tool_id]}")
        else:
            full = gi.tools.show_tool(tool_id)["id"]
            if full != tool_id:
                step["tool_id"] = step["content_id"] = full
                notes.append(f"    {tool_id:28} -> {full}")
            else:
                notes.append(f"    {tool_id:28} (built-in, unchanged)")
        resolved += 1
    print("\n".join(notes))

    native["name"] = "WF-B softmask (UDT edition) — resolved for this instance"
    work.mkdir(parents=True, exist_ok=True)
    (work / "workflow_resolved.ga").write_text(json.dumps(native, indent=2) + "\n", encoding="utf-8")
    print(f"  rendered {resolved} step(s) -> {work / 'workflow_resolved.ga'}")

    imported = gi.workflows.import_workflow_dict(native)
    print(f"  imported {imported['id']}")
    return imported["id"]


#: Suffixes stripped to derive an element identifier, longest first within each group.
COMPRESSION_EXT = (".gz", ".bz2", ".xz", ".zst")
SEQUENCE_EXT = (".fasta", ".fna", ".fas", ".fa", ".seq")


def collection_name(path: pathlib.Path) -> str:
    """The strain key for a collection element: the filename with its sequence suffixes removed.

    ⛔ NOT `pathlib.suffixes`, WHICH SPLITS ON EVERY DOT. `"".join(p.suffixes)` turns
    `GCA_900626175.2_cs10_genomic.fna` into `GCA_900626175` -- it treats `.2_cs10_genomic` as a
    suffix. Two versions of one strain (`PvP01.v1.fa.gz`, `PvP01.v2.fa.gz`) then collapse to the
    SAME identifier, and a list collection with duplicate element identifiers is a silent data
    error: the element identifier is the strain key every downstream track and
    `verify_softmask_outputs.py` joins on.
    """
    name = path.name
    for ext in COMPRESSION_EXT:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    for ext in SEQUENCE_EXT:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name or path.name


def upload_collection(gi: GalaxyInstance, history_id: str, fastas: list[pathlib.Path]) -> str:
    """Upload each FASTA and gather them into a list collection keyed by filename stem.

    ⚠ upload_file, NOT paste_content. paste_content reads the whole file into memory and posts it
    as a form field, which is fine for a test chunk and wrong for a real assembly -- cs10 chromosome
    1 alone is 102 MB and a whole genome is ~900 MB. It also cannot carry gzip: the bytes would be
    decoded as text.

    ⚠ The datatype is chosen from the SUFFIX, and it matters. Declaring a .gz file as `fasta` makes
    Galaxy hand the compressed bytes to the first tool as though they were sequence; declaring it
    `fasta.gz` lets Galaxy decompress or pass through as each tool needs.
    """
    ids = []
    for p in fastas:
        ftype = "fasta.gz" if p.suffix == ".gz" else "fasta"
        up = gi.tools.upload_file(str(p), history_id, file_type=ftype)
        ids.append((collection_name(p), up["outputs"][0]["id"]))
    for nm, ds in ids:
        await_dataset(gi, ds, f"upload {nm!r}")   # a 100 MB upload takes minutes, not seconds
    desc = {"collection_type": "list", "name": "assemblies",
            "element_identifiers": [{"name": n, "src": "hda", "id": i} for n, i in ids]}
    hdca = gi.histories.create_dataset_collection(history_id, desc)
    print(f"  uploaded {len(ids)} FASTA(s) into collection {hdca['id']}")
    return hdca["id"]


def report_upgrade_messages(gi: GalaxyInstance, wf_id: str, inputs: dict, history_id: str) -> dict | None:
    """Print the parameters `allow_tool_state_corrections` is about to silence.

    Returns the invocation if the unflagged attempt SUCCEEDED -- there was nothing to silence, and
    the attempt is therefore already a real run, so use it rather than invoking twice. Otherwise
    None, having printed what the flagged run below will accept.

    ⛔ THE FLAG DOES NOT CORRECT ANYTHING. Galaxy computes the upgrade messages either way; the flag
    only decides whether they are raised or logged. lib/galaxy/workflow/modules.py::
    populate_module_and_state:

        if step.upgrade_messages:
            if allow_tool_state_corrections:
                log.debug('Workflow step "%i" had upgrade messages: %s', ...)
            else:
                raise exceptions.MessageException("Workflow step has upgrade messages", ...)

    `log.debug` goes to Galaxy's server log, not into any response we can read -- so passing the flag
    and walking away means a parameter could start defaulting differently and nothing here would ever
    say so. That is the silent-success shape this project keeps getting bitten by.

    The refusal carries the list in `err_data`, and is raised inside build_workflow_run_configs
    BEFORE any invocation row or job exists, so an unflagged attempt is a free read of what the flag
    would hide.

    ⚠ IT IS A FLOOR, NOT A CENSUS. The loop above raises on the FIRST step carrying messages, so one
    attempt reports one step. That is why fixing these by hand went sort_bed, then genomecov, then
    merge_bed -- three runs, three refusals, each revealing only the next one.

    ⚠ AND `exc.body` IS A JSON STRING, NOT A DICT. Reading it with .get() silently yields nothing and
    the preflight then reports "clean" for a workflow with messages -- measured, not assumed.
    """
    try:
        inv = gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history_id)
    except Exception as exc:  # noqa: BLE001 - any refusal shape is worth reading
        body = getattr(exc, "body", None)
        if isinstance(body, str):
            try:
                body = json.loads(body)
            except ValueError:
                body = None
        data = body.get("err_data") if isinstance(body, dict) else None
        if not data:
            print(f"  upgrade-message preflight: nothing readable ({str(exc)[:140]})")
            return None
        n = sum(len(v) for v in data.values())
        print(f"  upgrade-message preflight: {n} default(s) the run below will accept, "
              f"on the first offending step only")
        for order_index in sorted(data, key=int):
            for param, message in sorted(data[order_index].items()):
                print(f"    step {order_index}: {param}: {message}")
        return None
    # ⚠ A SUCCEEDING PREFLIGHT HAS ALREADY INVOKED. The refusal path costs nothing because Galaxy
    # raises before writing an invocation row; the success path scheduled a real run. Hand it back.
    print("  upgrade-message preflight: none -- the flag is inert for this workflow")
    return inv


def await_invocation(gi: GalaxyInstance, invocation_id: str) -> int:
    """Wait for every job to reach a terminal state, then report EVERY step. 0 iff all jobs are ok.

    ⛔ A TIMEOUT IS A FAILURE, NOT A PASS. The previous version polled `while now < deadline` and
    then fell straight into the reporting block, where an invocation step's `state` of `scheduled`
    counted as OK -- so a run that exhausted POLL_CEILING with jobs still executing exited 0. The
    one thing this function exists to report is whether the run finished, and it reported success
    for a run that had not.

    ⛔ AND THE PASS/FAIL DECISION READS JOB STATES, NOT STEP STATES. An invocation step's `state` is
    its SCHEDULING state (`new`/`ready`/`scheduled`); a step whose job died still reads `scheduled`.
    The step list is printed for orientation only. The verdict comes from the job-state summary.

    ⚠ The wait cannot stop at "all current jobs are terminal" either: with a mapped-over collection
    Galaxy creates jobs incrementally, so `{"ok": 1}` is reachable while most of the graph has not
    been scheduled yet. It waits for the INVOCATION to leave `new`/`ready` as well.
    """
    deadline = time.monotonic() + POLL_CEILING
    timed_out = True
    while time.monotonic() < deadline:
        detail = gi.invocations.show_invocation(invocation_id)
        states = gi.invocations.get_invocation_summary(invocation_id).get("states", {})
        pending = {k: v for k, v in states.items() if k in ("new", "queued", "running", "paused")}
        scheduling_done = detail.get("state") not in SCHEDULING_IN_PROGRESS
        if scheduling_done and states and not pending:
            timed_out = False
            break
        print(f"    ... invocation={detail.get('state')} jobs={states or '{}'}")
        time.sleep(POLL_SECONDS)

    detail = gi.invocations.show_invocation(invocation_id)
    states = gi.invocations.get_invocation_summary(invocation_id).get("states", {})
    for step in detail.get("steps", []):
        label = step.get("workflow_step_label") or f"step {step.get('order_index')}"
        print(f"    {label:24} {step.get('state') or '-'}")
    print(f"  invocation state: {detail.get('state')}  |  jobs: {states or '{}'}")

    if timed_out:
        print(f"  ⛔ TIMED OUT after {POLL_CEILING}s with jobs still pending -- NOT a pass.")
        return 1
    bad = {k: v for k, v in states.items() if k != "ok"}
    if bad:
        print(f"  ⛔ jobs did not all succeed: {bad}")
        return 1
    if not states:
        print("  ⛔ the invocation produced NO jobs at all -- nothing ran.")
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=pathlib.Path, action="append", default=[],
                    help="an assembly FASTA; repeat for more. Becomes the input collection.")
    ap.add_argument("--register-only", action="store_true")
    ap.add_argument("--work", type=pathlib.Path, default=ROOT / "build/softmask_udt")
    args = ap.parse_args()

    gi = connect()
    print("Registering UDTs")
    uuids = register_all(gi)
    if args.register_only:
        return 0
    if not args.fasta:
        sys.exit("Need at least one --fasta (or use --register-only).")

    print("Rendering the instance-resolved workflow")
    wf_id = render_and_import(gi, uuids, args.work)

    history = gi.histories.create_history(name="WF-B softmask (UDT edition)")
    print(f"  history {gi.base_url}/histories/view?id={history['id']}")
    hdca = upload_collection(gi, history["id"], args.fasta)

    handles = gi.workflows.show_workflow(wf_id)["inputs"]
    inputs = {sid: {"src": "hdca", "id": hdca} for sid in handles}
    # ⚠ allow_tool_state_corrections IS REQUIRED, and the reason is not obvious.
    # Galaxy treats every parameter a `state:` block leaves unset as an "upgrade message" and
    # REFUSES the whole invocation over it -- while the import succeeds silently, so the workflow
    # looks correct in the editor and cannot run. gxformat2 `state:` blocks are conventionally
    # partial (the committed softmask.gxwf.yml sets three keys on a tool with a dozen), so without
    # this flag a portable workflow is unrunnable unless every default is transcribed by hand. That
    # was three consecutive refusals here -- sort_bed, then genomecov, then merge_bed -- each
    # naming a different set of untouched booleans.
    #
    # ⛔ IT IS NOT A BLANKET "IGNORE PROBLEMS" SWITCH, but it does mean tool defaults, not this
    # file, decide any parameter the workflow does not name. Parameters that MATTER to the analysis
    # are still set explicitly in the workflow (report_select, max, soft) precisely so that an
    # upstream default change cannot move them silently.
    inv = report_upgrade_messages(gi, wf_id, inputs, history["id"])
    if inv is None:
        inv = gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history["id"],
                                           allow_tool_state_corrections=True)
    base = gi.base_url
    print(f"  INVOKED {inv['id']} -> {base}/workflows/invocations/{inv['id']}")
    return await_invocation(gi, inv["id"])


if __name__ == "__main__":
    sys.exit(main())
