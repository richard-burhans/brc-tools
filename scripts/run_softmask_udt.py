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
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

ROOT = pathlib.Path(__file__).resolve().parent.parent
WORKFLOW = ROOT / "workflows/softmask/softmask_udt.gxwf.yml"
UDT_DIR = ROOT / "udt"
#: Only the tools this workflow uses. `env_probe` lives in udt/ too and is a diagnostic.
UDTS = ("fasta_uppercase", "dustmasker_bed3", "windowmasker_bed3", "tantan_bed3",
        "lc_classify", "samtools_faidx")
POLL_SECONDS = 20
POLL_CEILING = 3600


def connect() -> GalaxyInstance:
    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    return GalaxyInstance(url=url.rstrip("/"), key=key)


def register(gi: GalaxyInstance) -> dict[str, str]:
    """Create each UDT, returning {tool_id: uuid}.

    ⚠ THERE IS NO UPDATE. Every create makes a NEW tool, so re-running this leaves the previous
    definition beside the new one. That is Galaxy's behaviour, not a bug here, but it means the
    uuid returned by THIS call is the only one guaranteed to match the YAML on disk -- never reuse
    a uuid recorded by an earlier run after editing a definition.
    """
    key = gi.key
    mapping = {}
    for name in UDTS:
        path = UDT_DIR / f"{name}.gxtool.yml"
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        created = gi.make_post_request(f"{gi.url}/unprivileged_tools",
                                       payload={"representation": doc}, params={"key": key})
        mapping[doc["id"]] = created["uuid"]
        print(f"  registered {doc['id']:24} v{doc['version']} -> {created['uuid']}")
    return mapping


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


def upload_collection(gi: GalaxyInstance, history_id: str, fastas: list[pathlib.Path]) -> str:
    ids = []
    for p in fastas:
        up = gi.tools.paste_content(p.read_text(encoding="utf-8"), history_id, file_type="fasta")
        ids.append((p.stem, up["outputs"][0]["id"]))
    for _, ds in ids:
        for _ in range(80):
            if gi.datasets.show_dataset(ds)["state"] in ("ok", "error"):
                break
            time.sleep(3)
    desc = {"collection_type": "list", "name": "assemblies",
            "element_identifiers": [{"name": n, "src": "hda", "id": i} for n, i in ids]}
    hdca = gi.histories.create_dataset_collection(history_id, desc)
    print(f"  uploaded {len(ids)} FASTA(s) into collection {hdca['id']}")
    return hdca["id"]


def await_invocation(gi: GalaxyInstance, invocation_id: str) -> int:
    """Report EVERY step's state, never a single roll-up.

    A workflow that runs nine of twelve steps and fails three is not "failed" in any useful sense;
    which steps failed is the whole content of the result.
    """
    deadline = time.monotonic() + POLL_CEILING
    while time.monotonic() < deadline:
        summary = gi.invocations.get_invocation_summary(invocation_id)
        states = summary.get("states", {})
        if states and not (states.get("new") or states.get("queued") or states.get("running")):
            break
        print(f"    ... {states}")
        time.sleep(POLL_SECONDS)
    detail = gi.invocations.show_invocation(invocation_id)
    failed = 0
    for step in detail.get("steps", []):
        js = step.get("state") or "-"
        label = step.get("workflow_step_label") or f"step {step.get('order_index')}"
        flag = "" if js in ("scheduled", "ok", "-") else "  <-- NOT OK"
        if flag:
            failed += 1
        print(f"    {label:24} {js}{flag}")
    print(f"  invocation state: {detail.get('state')}  |  jobs: {gi.invocations.get_invocation_summary(invocation_id).get('states')}")
    return 1 if failed else 0


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
    uuids = register(gi)
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
    inv = gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history["id"])
    base = gi.base_url
    print(f"  INVOKED {inv['id']} -> {base}/workflows/invocations/{inv['id']}")
    return await_invocation(gi, inv["id"])


if __name__ == "__main__":
    sys.exit(main())
