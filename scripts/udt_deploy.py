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
import json
import os
import pathlib
import sys
import time

import requests
import yaml
from bioblend.galaxy import GalaxyInstance

#: --instance name -> the suffix its $GALAXY_URL/$GALAXY_API_KEY pair carries. "main" is
#: usegalaxy.org and takes no suffix; the sibling stage_wfa_test.py uses the same "_3" for laila.
INSTANCES = {"main": "", "vgp": "_2", "laila": "_3"}


def creds(instance: str) -> tuple[str, str]:
    """⚠ THE SERVER IS NOT INTERCHANGEABLE, AND SLOTS ARE THE REASON TO SAY SO HERE. A UDT gets
    GALAXY_SLOTS=1 on usegalaxy.org (measured, job 79420986) but whatever `local_slots` the local
    job_conf declares on laila -- 16 today. So a wrapper reading $GALAXY_SLOTS runs single-threaded
    on the server we ship to and multi-threaded on the one we test on. Correctness transfers between
    them; runtimes and thread-related behaviour do not.
    """
    try:
        suffix = INSTANCES[instance]
    except KeyError:
        raise SystemExit(
            f"unknown --instance {instance!r}; choose from {', '.join(INSTANCES)}"
        ) from None
    try:
        return os.environ[f"GALAXY_URL{suffix}"], os.environ[f"GALAXY_API_KEY{suffix}"]
    except KeyError as exc:
        raise SystemExit(f"--instance {instance} needs {exc.args[0]} in the environment") from None


def fetch_one(url: str, key: str, history_id: str, path: pathlib.Path) -> str:
    """Upload one file and return its dataset id.

    ⛔ NOT `gi.tools.paste_content`, AND NOT /api/tools. paste_content posts to the classic
    `upload1` tool, which a 26.1 server does not carry -- laila's panel has no `upload1` at all, and
    the failure surfaces as `{"err_msg": "Tool not found.", "err_code": 400014}`, which reads like
    the UDT registration failed when the UDT registered perfectly well. Posting `__DATA_FETCH__` to
    /api/tools is refused too ("must use alternative endpoint"). /api/tools/fetch is the endpoint
    that works, and it is what usegalaxy.org uses as well, so this is not a laila special case.
    """
    with path.open("rb") as fh:
        r = requests.post(
            url.rstrip("/") + "/api/tools/fetch",
            headers={"x-api-key": key},
            data={"history_id": history_id,
                  "targets": json.dumps([{"destination": {"type": "hdas"},
                                          "elements": [{"src": "files", "name": path.name,
                                                        "ext": "fasta"}]}])},
            files={"files_0|file_data": (path.name, fh)},
            timeout=600,
        )
    r.raise_for_status()
    return r.json()["outputs"][0]["id"]


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
    ap.add_argument("--instance", default="main", choices=sorted(INSTANCES))
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
    ds = fetch_one(url, key, h["id"], args.smoke_fasta)
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
