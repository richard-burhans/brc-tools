#!/usr/bin/env python3
"""Register a User-Defined Tool on a Galaxy server and smoke-test it.

WHY A SCRIPT AND NOT A SNIPPET. Porting a wrapper is a LOOP -- convert, register, run, read the
error, fix, repeat -- and each turn of it needs the same four API calls. Retyping them invites the
small differences that make one attempt not comparable with the last.

⚠ THREE API DETAILS THAT ARE NOT IN THE OBVIOUS PLACE, each of which cost a failed attempt:

  * A UDT is INVOKED BY `tool_uuid`, not by `tool_id`. Posting the id to /api/tools returns
    "Tool not found", which reads like the registration failed when it did not.
  * There is NO UPDATE. Every `create` makes a new tool, so an unchanged `version` leaves the old
    definition sitting beside the new one and `--bump` exists to keep them distinguishable.
  * A tool with NO data inputs fails before reaching a node -- no runner, no stderr, error in about
    a second. Give the smoke test a real input.

    python3 scripts/udt_deploy.py udt/dustmasker.yml --smoke-fasta test.fa
    python3 scripts/udt_deploy.py udt/dustmasker.yml --register-only
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance


def creds(instance: str) -> tuple[str, str]:
    if instance == "main":
        return os.environ["GALAXY_URL"], os.environ["GALAXY_API_KEY"]
    return os.environ["GALAXY_URL_2"], os.environ["GALAXY_API_KEY_2"]


def wait(gi: GalaxyInstance, job_id: str, timeout: int = 1800) -> dict:
    t0 = time.time()
    while time.time() - t0 < timeout:
        j = gi.jobs.show_job(job_id, full_details=True)
        if j.get("state") in ("ok", "error", "deleted", "paused"):
            return j
        time.sleep(10)
    return gi.jobs.show_job(job_id, full_details=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("definition", type=pathlib.Path)
    ap.add_argument("--instance", default="main")
    ap.add_argument("--register-only", action="store_true")
    ap.add_argument("--smoke-fasta", type=pathlib.Path,
                    help="a small FASTA to run the tool on. ⚠ Without an input the job dies before "
                         "reaching a node, with no stderr to read.")
    ap.add_argument("--bump", help="override the version, since create never updates in place")
    ap.add_argument("--input-name", default="input", help="the tool's data input parameter")
    args = ap.parse_args()

    doc = yaml.safe_load(args.definition.read_text())
    if args.bump:
        doc["version"] = args.bump
    url, key = creds(args.instance)
    gi = GalaxyInstance(url=url, key=key)

    r = gi.make_post_request(url.rstrip("/") + "/api/unprivileged_tools",
                             payload={"representation": doc}, params={"key": key})
    uuid = r["uuid"]
    print(f"registered {doc['id']} v{doc['version']} -> {uuid}")
    if args.register_only or not args.smoke_fasta:
        return 0

    h = gi.histories.create_history(name=f"UDT smoke: {doc['id']} {doc['version']}")
    up = gi.tools.paste_content(args.smoke_fasta.read_text(), h["id"], file_type="fasta")
    ds = up["outputs"][0]["id"]
    for _ in range(60):
        if gi.datasets.show_dataset(ds)["state"] == "ok":
            break
        time.sleep(3)

    j = gi.make_post_request(
        url.rstrip("/") + "/api/tools",
        payload={"history_id": h["id"], "tool_uuid": uuid,          # uuid, NOT tool_id
                 "inputs": {args.input_name: {"src": "hda", "id": ds}}}, params={"key": key})
    job = wait(gi, j["jobs"][0]["id"])
    print(f"  job {job['state']}  exit={job.get('exit_code')}")
    print(f"  history: {url.rstrip('/')}/histories/view?id={h['id']}")
    if job["state"] != "ok":
        # ⚠ Print BOTH streams: a container that fails to pull leaves tool_stderr empty and puts
        # nothing useful in stdout either, which is itself the diagnosis.
        print("  stderr:", (job.get("stderr") or "(empty)")[:600])
        print("  stdout:", (job.get("stdout") or "(empty)")[:300])
        return 1
    out = next(iter(job["outputs"].values()))
    d = gi.datasets.show_dataset(out["id"])
    print(f"  output: {d.get('extension')} {d.get('file_size')} bytes")
    body = gi.datasets.download_dataset(d["id"], use_default_filename=False)
    text = body.decode("utf-8", "replace") if isinstance(body, bytes) else str(body)
    print("  first lines:")
    for line in text.splitlines()[:5]:
        print("   ", line[:100])
    if not text.strip():
        print("  ⛔ EMPTY OUTPUT — the job succeeded and produced nothing, which is a silent failure")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
