#!/usr/bin/env python3
"""Invoke the pggb-pangenome-build workflow on the prepared Pv history.

Usage:
    GALAXY_URL=http://localhost:8080 GALAXY_API_KEY=... \\
    python scripts/run_workflow.py [--workflow workflows/pggb-pangenome-build/*.ga]

Reads history_id + collection_id from execution/history.json
(produced by upload_history.py), invokes the workflow, waits, downloads
the final GFA + .og to execution/outputs/, dumps invocation summary
to execution/invocation.json.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from bioblend.galaxy import GalaxyInstance

REPO = Path(__file__).resolve().parents[1]
HIST_JSON = REPO / "execution" / "history.json"
INV_JSON = REPO / "execution" / "invocation.json"
OUT_DIR = REPO / "execution" / "outputs"


WF_PARAMS = {
    "ref_strain": "GCA_900093555.2",
    "segment_length": 5000,
    "mapping_id": 90,
    "n_haplotypes": 8,
    "min_match": 23,
    "pansn_delim": "#",
}


def wait_for_invocation(gi, inv_id, timeout=8 * 3600, poll=30):
    start = time.time()
    while time.time() - start < timeout:
        inv = gi.invocations.show_invocation(inv_id)
        state = inv["state"]
        # ⛔ `completed` MUST BE HERE, and `scheduled` MUST STAY. Galaxy settles in either one,
        # unpredictably: 50 `completed` and 9 `scheduled` across 60 invocations on a real account,
        # interleaved through the whole timeline, one sitting `scheduled` for 6.8 days with every
        # job ok. lib/galaxy/schema/invocation.py defines `scheduled` as "workflow has been
        # scheduled" and `completed` as "all jobs have reached terminal states", and the transition
        # between them is written from a separate workflow_invocation_completion row -- so an
        # invocation whose completion hook never fires stays `scheduled` forever. Testing for only
        # one of the two means a wholly successful run never returns: the loop falls through to
        # `raise TimeoutError` after the full 8 hours. The job states below are the real test.
        if state in ("ok", "scheduled", "completed"):
            jobs = gi.invocations.get_invocation_summary(inv_id).get("states", {})
            if jobs.get("running", 0) == 0 and jobs.get("new", 0) == 0 and jobs.get("queued", 0) == 0:
                return inv
        if state in ("failed", "cancelled"):
            return inv
        print(f"  state={state} jobs={gi.invocations.get_invocation_summary(inv_id).get('states', {})}",
              flush=True)
        time.sleep(poll)
    raise TimeoutError("Workflow did not finish within timeout")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow", default=str(REPO / "workflows/pggb-pangenome-build/pggb-pangenome-build.ga"))
    parser.add_argument("--history-name", default="pv-pangenome-2026-05-25 (run)")
    args = parser.parse_args()

    url = os.environ.get("GALAXY_URL", "http://localhost:8080")
    key = os.environ.get("GALAXY_API_KEY")
    if not key:
        print("ERROR: set GALAXY_API_KEY", file=sys.stderr)
        return 1
    gi = GalaxyInstance(url, key=key)

    if not HIST_JSON.exists():
        print(f"ERROR: {HIST_JSON} missing; run upload_history.py first", file=sys.stderr)
        return 1
    hist = json.loads(HIST_JSON.read_text())
    history_id = hist["history_id"]
    collection_id = hist["collection_id"]

    wf_path = Path(args.workflow)
    if not wf_path.exists():
        print(f"ERROR: workflow {wf_path} missing", file=sys.stderr)
        return 1

    wf = gi.workflows.import_workflow_from_local_path(str(wf_path))
    wf_id = wf["id"]
    print(f"workflow_id={wf_id}")

    inv = gi.workflows.invoke_workflow(
        wf_id, history_id=history_id,
        inputs={"0": {"src": "hdca", "id": collection_id}},
        params=WF_PARAMS,
    )
    inv_id = inv["id"]
    print(f"invocation_id={inv_id}")

    final = wait_for_invocation(gi, inv_id)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    contents = gi.histories.show_history(history_id, contents=True, details="all")
    interesting = [
        d for d in contents
        if d["state"] == "ok"
        and any(d["name"].endswith(ext) for ext in (".gfa.gz", ".og", ".og.lay", ".vcf"))
    ]
    for d in interesting:
        target = OUT_DIR / d["name"]
        gi.datasets.download_dataset(d["id"], file_path=str(target), use_default_filename=False)
        print(f"  downloaded {target.name}")

    INV_JSON.write_text(json.dumps({
        "url": url,
        "workflow_id": wf_id,
        "invocation_id": inv_id,
        "history_id": history_id,
        "final_state": final["state"],
        "downloaded": [d["name"] for d in interesting],
    }, indent=2))
    print(f"saved {INV_JSON}")
    return 0 if final["state"] == "ok" else 2


if __name__ == "__main__":
    sys.exit(main())
