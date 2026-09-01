"""Shared pieces of the softmask UDT scripts.

WHY THIS EXISTS. Four scripts drive the same workflow -- run_softmask_udt (the whole thing),
build_up_softmask (tier by tier), check_softmask_stages (one stage at a time) and
verify_softmask_outputs (the result). They had grown byte-identical copies of `connect`,
`await_dataset` and the invocation-state set, and near-identical copies of the UDT registration
and the FASTA counter.

⛔ THE DUPLICATED COPIES WERE NOT HARMLESS. Both hard-won corrections in this area had to be
applied three times by hand: `await_dataset` learning that `state in ("ok", "error")` is not a
readiness test, and SCHEDULING_IN_PROGRESS gaining `requires_materialization` and `cancelling`
from Galaxy's own enum. A fix that must be repeated N times is a fix that will eventually be
applied N-1 times.

Only things that were ALREADY the same live here. `main`, the renderers and the per-script
assertions differ for real reasons and stay where they are.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

ROOT = pathlib.Path(__file__).resolve().parent.parent
UDT_DIR = ROOT / "udt"
WORKFLOW = ROOT / "workflows/softmask/softmask_udt.gxwf.yml"

#: The UDTs the softmask workflow needs, in dependency order.
UDTS = ("fasta_uppercase", "dustmasker_bed3", "windowmasker_bed3", "tantan_bed3",
        "lc_classify", "samtools_faidx", "fastan_gdb", "fastan_scan", "fastan_bed",
        "masking_row", "masking_header")

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


def connect() -> GalaxyInstance:
    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    return GalaxyInstance(url=url.rstrip("/"), key=key)


def register_one(gi: GalaxyInstance, name: str) -> tuple[str, str, str]:
    """Create ONE UDT from udt/<name>.gxtool.yml, returning (tool_id, version, uuid)."""
    doc = yaml.safe_load((UDT_DIR / f"{name}.gxtool.yml").read_text(encoding="utf-8"))
    created = gi.make_post_request(f"{gi.url}/unprivileged_tools",
                                   payload={"representation": doc}, params={"key": gi.key})
    return doc["id"], str(doc["version"]), created["uuid"]


def register_all(gi: GalaxyInstance, names: tuple[str, ...] = UDTS, verbose: bool = True) -> dict[str, str]:
    """Create each UDT, returning {tool_id: uuid}.

    ⚠ THERE IS NO UPDATE. Every create makes a NEW tool, so re-running this leaves the previous
    definition beside the new one. That is Galaxy's behaviour, not a bug here, but it means the
    uuid returned by THIS call is the only one guaranteed to match the YAML on disk -- never reuse
    a uuid recorded by an earlier run after editing a definition.
    """
    mapping = {}
    for name in names:
        tool_id, version, uuid = register_one(gi, name)
        mapping[tool_id] = uuid
        if verbose:
            print(f"  registered {tool_id:24} v{version} -> {uuid}")
    if not verbose:
        print(f"  registered {len(mapping)} UDT(s)")
    return mapping


def await_dataset(gi: GalaxyInstance, dataset_id: str, label: str, tries: int = 1200) -> None:
    """Block until a dataset is ready, and DIE if it errored or never settled.

    ⛔ `state in ("ok", "error")` IS NOT A READINESS TEST. It was used as one here, so an upload
    that failed -- bad format detection, truncated transfer -- became a collection element and the
    whole workflow ran on it. Exhausting the retries fell through just as silently.
    """
    state = "unknown"
    for _ in range(tries):
        state = gi.datasets.show_dataset(dataset_id)["state"]
        if state == "ok":
            return
        if state in ("error", "discarded", "failed_metadata"):
            info = gi.datasets.show_dataset(dataset_id).get("misc_info") or ""
            sys.exit(f"{label} is in state {state!r} and cannot be used: {info[:200]}")
        time.sleep(5)
    sys.exit(f"{label} never became ready after {tries * 5}s (last state {state!r}).")


def fasta_stats(text: str) -> tuple[int, int, int]:
    """(sequences, residues, lowercase residues).

    Callers that want only two of the three unpack and discard; this was two functions differing
    solely in whether they counted headers.
    """
    seqs = res = low = 0
    for line in text.splitlines():
        if line.startswith(">"):
            seqs += 1
            continue
        res += len(line)
        low += sum(1 for c in line if "a" <= c <= "z")
    return seqs, res, low


class UpgradeMessagesRefused(RuntimeError):
    """Galaxy refused an invocation because some step's `state:` leaves a parameter unset."""

    def __init__(self, data: dict) -> None:
        self.data = data
        lines = [f"step {i}: {k}: {v}"
                 for i in sorted(data, key=int) for k, v in sorted(data[i].items())]
        super().__init__("Galaxy refused the invocation over unset parameters:\n  "
                         + "\n  ".join(lines))


def invoke(gi: GalaxyInstance, wf_id: str, inputs: dict, history_id: str) -> dict:
    """Invoke a workflow WITHOUT allow_tool_state_corrections, reporting what it would have hidden.

    ⛔ THE FLAG WAS NEVER A FIX. `workflow/modules.py::populate_module_and_state` either raises on a
    step's upgrade messages or, with the flag, calls `log.debug` -- to Galaxy's server log, which no
    response exposes. Passing it does not settle which value a parameter takes; it only removes the
    one place that would have told us the question was open. Every parameter this workflow's steps
    can take is now named in softmask_udt.gxwf.yml, so there is nothing to silence, and a refusal
    here is real news rather than noise to be switched off.

    ⚠ THE REFUSAL IS THE ONLY RELIABLE AUDIT. Comparing a step's `state:` against the tool's
    parameter list misses two shapes, both of which really occurred here: a parameter nested inside
    a repeat's conditional, and an OPTIONAL `data` input, which Galaxy counts as unset exactly like
    a required one -- leaving it unconnected is not the same as naming it null. Galaxy also raises
    on the FIRST offending step only, so one refusal is a floor, not a census.
    """
    try:
        return gi.workflows.invoke_workflow(wf_id, inputs=inputs, history_id=history_id)
    except Exception as exc:
        body = getattr(exc, "body", None)
        if isinstance(body, str):
            try:
                body = json.loads(body)
            except ValueError:
                body = None
        data = body.get("err_data") if isinstance(body, dict) else None
        if data:
            raise UpgradeMessagesRefused(data) from exc
        raise
