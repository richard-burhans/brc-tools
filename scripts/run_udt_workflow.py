#!/usr/bin/env python3
"""Register the UDTs a workflow needs, render an instance-resolved copy, and invoke it.

This is `run_softmask_udt.py` with the softmask knowledge removed: the workflow is an argument, the
UDT set is DERIVED from the workflow rather than listed, and the inputs are named on the command
line. `run_softmask_udt.py` keeps the softmask-specific parts (uploading assemblies, `--hdca`) and
imports the two generic pieces from here, so there is one copy of each.

    python3 scripts/run_udt_workflow.py \
        --workflow workflows/align_chain_project/align_chain_udt.gxwf.yml \
        --input masked_fastas=hdca:<id> --input sizes=hdca:<id> \
        --upload self_pairs=self_pairs.txt --upload relabel_map=relabel.tsv \
        --history-name "WF-C UDT edition"

⛔ THE COMMITTED WORKFLOW CANNOT BE INVOKED AS WRITTEN, AND IMPORTING IT SAYS OTHERWISE. A gxformat2
file whose steps name UDTs imports with every node present and no errors, then
`POST .../invocations` refuses with "the following required tools are not installed". A UDT is
addressed by UUID and its plain id is not in the toolbox the invocation gate consults. The fix is
the round trip: import the portable file, export Galaxy's NATIVE representation (which does carry
`tool_uuid`), resolve identities there, re-import that. Hand-authoring the native dict 500s.

⚠ THE UDT SET IS DERIVED, NOT LISTED, and that is the difference that makes this generic. Each
`udt/*.gxtool.yml` declares its own `id`; a workflow step naming one of those ids needs it
registered. A hardcoded tuple (`softmask_lib.UDTS`) is correct for exactly one workflow and is
silently wrong for the next -- registering too few fails at invoke with a list of missing tools,
and registering too many leaves unused tools on the account.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance
from softmask_lib import (
    ROOT,
    SCHEDULING_IN_PROGRESS,
    UDT_DIR,
    await_dataset,
    connect,
    invoke,
    register_one,
)

POLL_SECONDS = 20
POLL_CEILING = 21600


def udt_ids() -> dict[str, str]:
    """`{tool id: stem}` for every UDT document in udt/, read from the documents themselves."""
    out = {}
    for path in sorted(UDT_DIR.glob("*.gxtool.yml")):
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        out[doc["id"]] = path.name.removesuffix(".gxtool.yml")
    return out


def needed_udts(workflow: pathlib.Path) -> list[str]:
    """The UDT stems this workflow's steps name, in first-use order (each once)."""
    doc = yaml.safe_load(workflow.read_text(encoding="utf-8"))
    known, seen = udt_ids(), []
    for step in (doc.get("steps") or {}).values():
        stem = known.get(step.get("tool_id"))
        if stem and stem not in seen:
            seen.append(stem)
    return seen


def _fill(tool_inputs, state, connections, prefix=""):
    """Insert each unset parameter's own default into `state`; return the paths filled."""
    filled = []
    for i in tool_inputs:
        name, kind = i["name"], i.get("type")
        path = f"{prefix}{name}"
        if kind == "section":
            filled += _fill(i.get("inputs", []), state.setdefault(name, {}), connections, f"{path}|")
        elif kind == "conditional":
            test = i.get("test_param") or {}
            sub = state.setdefault(name, {})
            if test.get("name") and test["name"] not in sub:
                sub[test["name"]] = test.get("value")
                filled.append(f"{path}|{test['name']}={test.get('value')!r}")
            chosen = sub.get(test.get("name"))
            for case in i.get("cases", []):
                if case.get("value") == chosen:
                    filled += _fill(case.get("inputs", []), sub, connections, f"{path}|")
        elif kind == "repeat":
            continue                      # a repeat has no default instance to invent
        elif kind in ("data", "data_collection"):
            # ⛔ ONLY UNCONNECTED ONES. A connected input carries its edge in input_connections and
            # its state entry is Galaxy's business; writing null over it would cut the wire.
            if name not in state and path not in connections:
                state[name] = None
                filled.append(f"{path}=null (optional, unconnected)")
        elif name not in state and not name.startswith("__"):
            state[name] = i.get("value")
            filled.append(f"{path}={i.get('value')!r}")
    return filled


def fill_step_defaults(gi: GalaxyInstance, native: dict) -> None:
    """Write every unset parameter's own default into each step's tool_state, and SAY SO.

    ⛔ WHY THIS EXISTS RATHER THAN A LONGER WORKFLOW FILE. Galaxy answers an invocation whose step
    leaves a parameter unset with "No value found for X. Using default: Y", and this project's
    invoke() refuses on upgrade messages rather than silencing them with
    allow_tool_state_corrections. Naming every parameter of every step in the workflow would satisfy
    that -- and for KegAlign, whose seed, step, xdrop and scoring set DECIDE WHAT THE ALIGNMENT
    FINDS, the workflow does exactly that. But the nine UCSC chain tools carry dozens of mechanical
    knobs that no reader of the file wants to see, and burying five scientific choices in eighty
    lines of inherited defaults makes the file worse, not more reproducible.

    ⚠ SO THE DEFAULTS ARE FILLED HERE AND RECORDED THERE. Every value comes from the tool's own API
    on the instance being run, each fill is printed, and the resolved `.ga` written beside the run
    is the artifact that says what actually executed. Nothing is silenced: a parameter this cannot
    resolve still reaches invoke() as a refusal.
    """
    for step in native["steps"].values():
        tool_id = step.get("tool_id")
        if not tool_id or not step.get("tool_state"):
            continue
        try:
            spec = gi.tools.show_tool(tool_id, io_details=True)
        except Exception as exc:                       # a UDT is not in /api/tools at all
            print(f"    {tool_id}: defaults not fetched ({str(exc)[:40]})")
            continue
        state = json.loads(step["tool_state"])
        connections = set(step.get("input_connections") or {})
        filled = _fill(spec.get("inputs", []), state, connections)
        if filled:
            step["tool_state"] = json.dumps(state)
            label = step.get("label") or tool_id.split("/")[-2] if "/" in tool_id else tool_id
            print(f"    defaults filled for {label}: {', '.join(filled)}")


def render_and_import(gi: GalaxyInstance, workflow: pathlib.Path, uuids: dict[str, str],
                      work: pathlib.Path, name: str) -> str:
    """Portable gxformat2 -> native -> identities resolved -> re-imported. Returns workflow id."""
    portable = gi.workflows.import_workflow_dict(yaml.safe_load(workflow.read_text(encoding="utf-8")))
    native = gi.workflows.export_workflow_dict(portable["id"])
    gi.workflows.delete_workflow(portable["id"])  # a scaffold, not an artifact

    notes, resolved = [], 0
    for step in native["steps"].values():
        tool_id = step.get("tool_id")
        if not tool_id:
            continue
        if tool_id in uuids:
            step["tool_uuid"] = uuids[tool_id]
            notes.append(f"    {tool_id:44} -> uuid {uuids[tool_id]}")
        else:
            full = gi.tools.show_tool(tool_id)["id"]
            if full != tool_id:
                step["tool_id"] = step["content_id"] = full
                notes.append(f"    {tool_id:44} -> {full}")
            else:
                notes.append(f"    {tool_id:44} (resolves as written)")
        resolved += 1
    print("\n".join(notes))

    fill_step_defaults(gi, native)

    native["name"] = name
    work.mkdir(parents=True, exist_ok=True)
    out = work / f"{workflow.stem}_resolved.ga"
    out.write_text(json.dumps(native, indent=2) + "\n", encoding="utf-8")
    print(f"  rendered {resolved} step(s) -> {out}")

    imported = gi.workflows.import_workflow_dict(native)
    print(f"  imported {imported['id']}")
    return imported["id"]


def await_invocation(gi: GalaxyInstance, invocation_id: str) -> int:
    """Wait for every job to reach a terminal state, then report EVERY step. 0 iff all jobs are ok.

    ⛔ A TIMEOUT IS A FAILURE, NOT A PASS, and ⛔ THE VERDICT READS JOB STATES, NOT STEP STATES: an
    invocation step's `state` is its SCHEDULING state, so a step whose job died still reads
    `scheduled`. Both rules were learned the hard way here; see the history of this function in
    run_softmask_udt.py.
    """
    deadline = time.monotonic() + POLL_CEILING
    timed_out = True
    while time.monotonic() < deadline:
        detail = gi.invocations.show_invocation(invocation_id)
        states = gi.invocations.get_invocation_summary(invocation_id).get("states", {})
        pending = {k: v for k, v in states.items() if k in ("new", "queued", "running", "paused")}
        if detail.get("state") not in SCHEDULING_IN_PROGRESS and states and not pending:
            timed_out = False
            break
        print(f"    ... invocation={detail.get('state')} jobs={states or '{}'}", flush=True)
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


def parse_input(spec: str) -> tuple[str, dict]:
    """`label=hdca:ID` / `label=hda:ID` -> (label, {"src": .., "id": ..})."""
    label, _, ref = spec.partition("=")
    src, _, ident = ref.partition(":")
    if src not in ("hdca", "hda") or not ident:
        sys.exit(f"--input {spec!r}: expected label=hdca:ID or label=hda:ID")
    return label, {"src": src, "id": ident}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workflow", type=pathlib.Path, required=True)
    ap.add_argument("--input", action="append", default=[], metavar="LABEL=SRC:ID",
                    help="bind a workflow input to an existing collection or dataset")
    ap.add_argument("--upload", action="append", default=[], metavar="LABEL=PATH",
                    help="upload a file into the run history and bind it to that input")
    ap.add_argument("--history-name")
    ap.add_argument("--history", metavar="ID", help="run in an EXISTING history instead of a new one")
    ap.add_argument("--register-only", action="store_true")
    ap.add_argument("--work", type=pathlib.Path, default=ROOT / "build/udt_runs")
    args = ap.parse_args()

    gi = connect()
    stems = needed_udts(args.workflow)
    print(f"Registering {len(stems)} UDT(s) named by {args.workflow.name}")
    uuids = {}
    for stem in stems:
        tool_id, version, uuid = register_one(gi, stem)
        uuids[tool_id] = uuid
        print(f"  registered {tool_id:24} v{version} -> {uuid}")
    if args.register_only:
        return 0

    name = args.history_name or f"{args.workflow.stem} (UDT edition)"
    print("Rendering the instance-resolved workflow")
    wf_id = render_and_import(gi, args.workflow, uuids, args.work, f"{name} — resolved")

    history = ({"id": args.history} if args.history
               else gi.histories.create_history(name=name))
    print(f"  history {gi.base_url}/histories/view?id={history['id']}")

    inputs = dict(parse_input(s) for s in args.input)
    for spec in args.upload:
        label, _, path = spec.partition("=")
        up = gi.tools.upload_file(path, history["id"])
        ds = up["outputs"][0]["id"]
        await_dataset(gi, ds, f"upload {label!r}")
        inputs[label] = {"src": "hda", "id": ds}
        print(f"  uploaded {path} -> {label} ({ds})")

    # ⚠ INPUTS ARE BOUND BY LABEL, THEN TRANSLATED TO STEP IDS. show_workflow()["inputs"] is keyed
    # by step id with the label inside, so binding by label here fails loudly on a typo instead of
    # silently leaving an input unfilled -- which Galaxy reports much later, as a scheduling error.
    handles = gi.workflows.show_workflow(wf_id)["inputs"]
    by_label = {v.get("label") or v.get("uuid"): sid for sid, v in handles.items()}
    missing = sorted(set(by_label) - set(inputs))
    unknown = sorted(set(inputs) - set(by_label))
    if unknown:
        sys.exit(f"no such workflow input: {unknown} (this workflow takes {sorted(by_label)})")
    if missing:
        sys.exit(f"unbound workflow input(s): {missing}")
    bound = {by_label[label]: ref for label, ref in inputs.items()}

    inv = invoke(gi, wf_id, bound, history["id"])
    print(f"  INVOKED {inv['id']} -> {gi.base_url}/workflows/invocations/{inv['id']}")
    return await_invocation(gi, inv["id"])


if __name__ == "__main__":
    sys.exit(main())
