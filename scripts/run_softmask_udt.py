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
import pathlib
import sys

from bioblend.galaxy import GalaxyInstance
from run_udt_workflow import await_invocation, render_and_import
from softmask_lib import ROOT, WORKFLOW, await_dataset, connect, invoke, register_all

#: Only the tools this workflow uses. `env_probe` lives in udt/ too and is a diagnostic.
POLL_SECONDS = 20
POLL_CEILING = 21600          # 6 h: a 101 Mb chromosome, single-threaded, two-pass windowmasker


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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=pathlib.Path, action="append", default=[],
                    help="an assembly FASTA; repeat for more. Becomes the input collection.")
    ap.add_argument("--hdca", metavar="ID",
                    help="run on an EXISTING collection instead of uploading. The workflow "
                         "uppercases its input, so a collection that is already uppercased (or "
                         "already soft-masked) gives the same result as the raw assemblies -- the "
                         "mask is always recomputed de novo. Saves re-uploading large genomes.")
    ap.add_argument("--history-name", default="WF-B softmask (UDT edition)")
    ap.add_argument("--register-only", action="store_true")
    ap.add_argument("--work", type=pathlib.Path, default=ROOT / "build/softmask_udt")
    args = ap.parse_args()

    gi = connect()
    print("Registering UDTs")
    uuids = register_all(gi)
    if args.register_only:
        return 0
    if bool(args.fasta) == bool(args.hdca):
        sys.exit("Give exactly one of --fasta (upload) or --hdca (reuse), or --register-only.")

    print("Rendering the instance-resolved workflow")
    wf_id = render_and_import(gi, WORKFLOW, uuids, args.work,
                              "WF-B softmask (UDT edition) — resolved for this instance")

    history = gi.histories.create_history(name=args.history_name)
    print(f"  history {gi.base_url}/histories/view?id={history['id']}")
    if args.hdca:
        # Galaxy copies an input collection into the target history itself, so a run can reuse a
        # collection that lives elsewhere without duplicating gigabytes on the way in.
        hdca = args.hdca
        n = gi.dataset_collections.show_dataset_collection(hdca).get("element_count")
        print(f"  reusing collection {hdca} ({n} element(s))")
    else:
        hdca = upload_collection(gi, history["id"], args.fasta)

    handles = gi.workflows.show_workflow(wf_id)["inputs"]
    inputs = {sid: {"src": "hdca", "id": hdca} for sid in handles}
    # ⛔ NO allow_tool_state_corrections. Every parameter every step can take is named in the
    # workflow, so there is nothing for it to silence -- and it never fixed anything anyway: it only
    # swaps Galaxy's refusal for a log.debug on the server that no response exposes. A refusal here
    # is real news. See softmask_lib.invoke.
    inv = invoke(gi, wf_id, inputs, history["id"])
    base = gi.base_url
    print(f"  INVOKED {inv['id']} -> {base}/workflows/invocations/{inv['id']}")
    return await_invocation(gi, inv["id"])


if __name__ == "__main__":
    sys.exit(main())
