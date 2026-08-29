#!/usr/bin/env python3
"""Upload the 8 P. vivax FASTAs in data/raw/ to a Galaxy history as a list collection.

Usage:
    GALAXY_URL=http://localhost:8080 \\
    GALAXY_API_KEY=... \\
    python scripts/upload_history.py [--history-name pv-pangenome-2026-05-25]

Emits the history_id and collection_id to stdout (and to execution/history.json).
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
RAW = REPO / "data" / "raw"
OUT = REPO / "execution" / "history.json"


def wait_for_dataset(gi: GalaxyInstance, hist_id: str, ds_id: str, timeout: int = 600) -> str:
    start = time.time()
    while time.time() - start < timeout:
        ds = gi.histories.show_dataset(hist_id, ds_id)
        state = ds.get("state")
        if state == "ok":
            return state
        if state in ("error", "failed_metadata", "discarded"):
            raise RuntimeError(f"Dataset {ds_id} ended in state {state}: {ds.get('misc_info')}")
        time.sleep(5)
    raise TimeoutError(f"Dataset {ds_id} did not become ok in {timeout}s")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--history-name", default="pv-pangenome-2026-05-25")
    parser.add_argument("--collection-name", default="pv_assemblies")
    args = parser.parse_args()

    url = os.environ.get("GALAXY_URL", "http://localhost:8080")
    key = os.environ.get("GALAXY_API_KEY")
    if not key:
        print("ERROR: set GALAXY_API_KEY", file=sys.stderr)
        return 1

    gi = GalaxyInstance(url, key=key)

    fastas = sorted(RAW.glob("*.fa.gz"))
    if not fastas:
        print(f"ERROR: no FASTAs in {RAW}; run scripts/fetch_pv_assemblies.sh first", file=sys.stderr)
        return 1

    hist = gi.histories.create_history(name=args.history_name)
    print(f"history_id={hist['id']}")

    ds_ids = []
    for fa in fastas:
        strain = fa.name.replace(".fa.gz", "")
        upload = gi.tools.upload_file(
            str(fa), history_id=hist["id"], file_type="fasta.gz", file_name=strain,
        )
        ds_id = upload["outputs"][0]["id"]
        wait_for_dataset(gi, hist["id"], ds_id)
        ds_ids.append((strain, ds_id))
        print(f"  uploaded {strain} -> {ds_id}")

    coll_payload = {
        "collection_type": "list",
        "element_identifiers": [
            {"src": "hda", "id": ds_id, "name": strain} for strain, ds_id in ds_ids
        ],
        "name": args.collection_name,
    }
    coll = gi.histories.create_dataset_collection(hist["id"], coll_payload)
    print(f"collection_id={coll['id']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "url": url,
        "history_id": hist["id"],
        "collection_id": coll["id"],
        "datasets": ds_ids,
    }, indent=2))
    print(f"saved {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
