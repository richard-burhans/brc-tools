#!/usr/bin/env python3
"""Static port/parameter check for a `.gxwf.yml` workflow, before anyone invokes it.

    python3 scripts/check_workflow_ports.py workflows/softmask/softmask_udt.gxwf.yml

⛔ AN INVOCATION IS AN EXPENSIVE WAY TO FIND A TYPO. A workflow naming a step output that does not
exist, or setting a parameter the tool does not have, fails after Galaxy has scheduled it -- or
worse, does not fail: `allow_tool_state_corrections` downgrades the refusal to a `log.debug` on the
SERVER, which no API response exposes, and the job then runs with a default you did not choose.
Everything here is cheap and runs before any of that.

TWO HALVES, AND THE SECOND ONE WITHHOLDS RATHER THAN PASSES.
  * LOCAL -- the `brc-*` UDTs are files in this repo, so every `in:` port, every `state:` key and
    every `from_work_dir` is checkable offline and exactly.
  * REMOTE -- ToolShed and built-in tools are checked against a LIVE Galaxy, because the installed
    parameter tree is the only authority on what a wrapper actually accepts. Needs GALAXY_URL and
    GALAXY_API_KEY. ⚠ Without them this SKIPS and says so; it never reports a clean run over
    checks that did not execute.

Exit status is 1 if any problem was found, 0 otherwise. A SKIP is not a pass and prints as its own
line, so `--strict` is available for CI, where an unrun check should be a failure.
"""
import argparse
import json
import os
import pathlib
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

import yaml

# ⚠ scripts/ IS NOT A PACKAGE AND THIS SCRIPT IS RUN BY PATH, so the sibling module is importable
# only after its directory is on sys.path -- `python3 scripts/x.py` puts scripts/ there, but
# `python3 -m` and a symlinked entry point do not, and this has to work under all three.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from check_udt_definitions import lint as udt_lint      # noqa: I001 -- must follow the path insert

ROOT = pathlib.Path(__file__).resolve().parent.parent

#: ⛔ OUTSIDE THE REPOSITORY, ON PURPOSE. This used to be `ROOT / ".tool_schema_cache"`, which put
#: several hundred JSON files inside the working tree and left them out of commits only for as long
#: as somebody remembered a .gitignore line. That lasted one branch: a later branch cut from main --
#: which did not carry the ignore -- swept the whole cache into a `git add -A` and a one-file change
#: arrived as a hundred-file diff. A cache that lives outside the tree cannot be committed by
#: accident, and needs no cooperation from .gitignore to stay that way.
CACHE = pathlib.Path(
    os.environ.get("XDG_CACHE_HOME") or (pathlib.Path.home() / ".cache")) / "brc-tools" / "tool-schemas"


# ---------------------------------------------------------------- remote schema ----------------
def tool_schema(tid, url, key, use_cache=True):
    """The tool's parameter tree as the running Galaxy reports it."""
    c = CACHE / (tid.replace("/", "__").replace(":", "_") + ".json")
    if use_cache and c.exists():
        return json.loads(c.read_text())
    req = urllib.request.Request(
        f"{url}/api/tools/{urllib.parse.quote(tid, safe='')}?io_details=true",
        headers={"x-api-key": key})
    d = json.load(urllib.request.urlopen(req, timeout=120))
    CACHE.mkdir(parents=True, exist_ok=True)   # ⚠ parents: ~/.cache/brc-tools need not exist
    c.write_text(json.dumps(d))
    return d


def flatten(inputs, prefix=""):
    """Galaxy's nested input tree -> {addressable_name: param}.

    ⚠ A REPEAT IS REPORTED ONCE AND USED MANY TIMES. Galaxy describes `queries` with a single
    `queries_0|input2`, but a workflow may legitimately address `queries_1|…`, `queries_2|…` for
    however many instances it connects. Recording the repeat's PREFIX (`queries_*`) is what lets
    the caller accept those; matching the literal index would flag every multi-input `cat1` in the
    repository as broken, which is exactly what an earlier version of this did.
    """
    out = {}
    for p in inputs or []:
        name = f"{prefix}{p.get('name', '')}"
        out[name] = p
        t = p.get("type")
        if t == "conditional":
            tp = p.get("test_param") or {}
            out[f"{name}|{tp.get('name', '')}"] = tp
            for case in p.get("cases") or []:
                out.update(flatten(case.get("inputs"), f"{name}|"))
        elif t == "repeat":
            out[f"{name}_*"] = p                       # the wildcard the docstring describes
            for k, v in flatten(p.get("inputs"), f"{name}_*|").items():
                out[k] = v
        elif t == "section":
            out.update(flatten(p.get("inputs"), f"{name}|"))
    return out


def source_of(ref):
    """The `step/output` a connection names, in either gxformat2 spelling, else None.

    ⛔ `in:` HAS TWO FORMS AND THE DICT ONE IS CANONICAL. A connection is written either as the
    shorthand `port: step/output` or as `port: {source: step/output}`. An earlier version of this
    script tested `isinstance(ref, str)` and SKIPPED everything else, so every dict-form connection
    went unvalidated AND uncounted for reachability -- which made it report four correctly-wired
    steps in msa/selection/ucsc_hub as UNREACHABLE while silently checking none of their links.
    A checker that cannot see half the edges is worse than no checker, because it reports on them.
    """
    if isinstance(ref, str):
        return ref
    if isinstance(ref, dict):
        s = ref.get("source")
        return s if isinstance(s, str) else None
    return None


def sources_of(ref):
    """EVERY `step/output` a connection names -- a port may take more than one.

    ⛔ `in:` HAS THREE SPELLINGS AND THIS SCRIPT KNEW TWO. Besides `port: step/out` and
    `port: {source: step/out}`, a port that accepts multiple datasets takes a LIST of either form:

        image_content_input:
          - sourmash_panel/heatmap
          - sourmash_panel/dendrogram

    `source_of` returns None for a list, and every caller treated None as "not a connection", so
    both of those edges were invisible -- unvalidated, uncounted for reachability, and missing from
    the graph the collection checks below reason over. It showed up as WF-A's `assemblies` and
    `proteomes` being reported as never meeting when they meet at exactly this port.
    """
    if isinstance(ref, list):
        return [s for r in ref for s in sources_of(r)]
    one = source_of(ref)
    return [one] if one else []


def addressable(key, params):
    """Does `key` name something in `params`, allowing any repeat index?"""
    if key in params:
        return True
    # queries_1|input2 -> queries_*|input2
    import re
    wild = re.sub(r"_\d+(\|)", r"_*\1", key)
    if wild in params or wild.split("|")[0] in params:
        return True
    return key.split("|")[0] in params



# ---------------------------------------------------------------- self-test --------------------
#: (fixture, label, find, replace, expected finding). Offline cases only -- the remote ones need a
#: Galaxy, and a self-test that silently needs credentials is a self-test that silently does not run.
DEFECTS = [
    ("workflows/softmask/softmask_udt.gxwf.yml", "UDT output typo",
     "fasta: uppercase/output", "fasta: uppercase/nope", "BAD-OUTPUT"),
    ("workflows/softmask/softmask_udt.gxwf.yml", "step does not exist",
     "input: uppercase/output", "input: ghost/output", "BROKEN-LINK"),
    ("workflows/softmask/softmask_udt.gxwf.yml", "workflow output names a ghost step",
     "outputSource: faidx/sizes", "outputSource: nosuch/sizes", "BROKEN-OUTPUT"),
    ("workflows/softmask/softmask_udt.gxwf.yml", "UDT parameter typo",
     "      input: assemblies", "      inputt: assemblies", "UNKNOWN-PARAM"),
    ("workflows/inventory/inventory_udt.gxwf.yml", "handoff tag without the name: prefix",
     "{tags: name:wfc_sizes}", "{tags: wfc_sizes}", "PLAIN-TAG"),
    ("workflows/inventory/inventory_udt.gxwf.yml", "a rename string restyled away from its target",
     "'WF-C2 input: anchor_bed12s'", "'WF-C2 input: anchor BED12'", "PARITY-DRIFT"),
    # ⚠ THE `want` NAMES THE INPUT, NOT JUST THE CODE. Both checks below already fire on the
    # unmodified WF-A port (for `anchor_gene_gff3s`, correctly), so a fixture that only greps for
    # the code would pass without the injection doing anything at all -- the classic vacuous test.
    # Naming the input the injection is supposed to newly implicate makes the case real.
    ("workflows/inventory/inventory_udt.gxwf.yml", "a collection nothing receives whole",
     "      assemblies: assemblies\n", "      assemblies: proteomes\n",
     "EMPTY-COLLECTION `assemblies`"),
    # ⚠ AND THE INJECTION HAS TO CUT THE EDGE, NOT RENAME THE PORT. Renaming
    # `image_content_input` left both sources exactly where they were -- the graph does not care
    # what a port is called -- so the fixture changed nothing and the case reported "not detected"
    # while the check was working correctly.
    ("workflows/inventory/inventory_udt.gxwf.yml", "two collections that stop meeting",
     "      image_content_input:\n        - sourmash_panel/heatmap\n"
     "        - sourmash_panel/dendrogram\n", "      image_content_input: []\n",
     "UNRECONCILED `assemblies` and `proteomes`"),
    ("workflows/msa/msa.gxwf.yml", "dict-form connection to a ghost step",
     "source: group_cds_by_og/pep", "source: ghoststep/pep", "BROKEN-LINK"),
]


def self_test() -> int:
    """Prove each defect class is actually detected, and that a clean workflow stays clean."""
    import subprocess
    import tempfile

    me = [sys.executable, str(pathlib.Path(__file__).resolve())]
    env = dict(os.environ, GALAXY_URL="", GALAXY_API_KEY="")     # offline, deterministically
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        # ⛔ THE FIXTURE MUST KEEP ITS NAME AND ITS SIBLING. Written to a generic `case.gxwf.yml`,
        # every fixture silently skipped the whole parity half -- that block runs only for a file
        # ending `_udt.gxwf.yml` and only when the classic sits beside it, so a PARITY-DRIFT case
        # reported "not detected" while the check was in fact never reached. A self-test that can
        # skip the check it is testing is the failure mode this file exists to prevent.
        def stage(rel, text):
            tmp = pathlib.Path(td) / pathlib.Path(rel).name
            tmp.write_text(text)
            sib = (ROOT / rel).with_name(pathlib.Path(rel).name.replace("_udt.gxwf.yml",
                                                                       ".gxwf.yml"))
            if sib.exists() and sib.name != pathlib.Path(rel).name:
                (pathlib.Path(td) / sib.name).write_text(sib.read_text())
            return tmp

        for rel, label, old, new, want in DEFECTS:
            src = (ROOT / rel).read_text()
            if src.count(old) < 1:
                print(f"  BROKEN-FIXTURE  {label}: {old!r} no longer occurs in {rel}")
                failed += 1
                continue
            tmp = stage(rel, src.replace(old, new, 1))
            r = subprocess.run([*me, str(tmp)], capture_output=True, text=True, env=env)
            hit = want in r.stdout
            print(f"  {'pass' if hit else 'FAIL'}  {label:<38} expect {want}")
            failed += not hit
        clean = "workflows/softmask/softmask_udt.gxwf.yml"
        tmp = stage(clean, (ROOT / clean).read_text())
        r = subprocess.run([*me, str(tmp)], capture_output=True, text=True, env=env)
        no_problems = "PROBLEM" not in r.stdout
        print(f"  {'pass' if no_problems else 'FAIL'}  {'unmodified stays clean':<38} expect none")
        failed += not no_problems
    total = len(DEFECTS) + 1
    print(f"\n  {total - failed}/{total} self-test cases behaved correctly")
    return 1 if failed else 0


# ---------------------------------------------------------------- the check --------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workflow", type=pathlib.Path, nargs="?")
    ap.add_argument("--no-cache", action="store_true", help="re-fetch every tool schema")
    ap.add_argument("--self-test", action="store_true",
                    help="inject known defects into this repo's own workflows and confirm each is "
                         "caught; needs no Galaxy. A green checker nobody has seen go red is not "
                         "evidence of anything.")
    ap.add_argument("--strict", action="store_true",
                    help="treat a SKIPPED remote half as a failure (use in CI)")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if args.workflow is None:
        ap.error("a workflow path is required (or --self-test)")

    wf = yaml.safe_load(args.workflow.read_text())
    steps, wf_in, wf_out = wf["steps"], wf.get("inputs") or {}, wf.get("outputs") or {}

    udt = {}
    for p in sorted((ROOT / "udt").glob("*.gxtool.yml")):
        d = yaml.safe_load(p.read_text())
        udt[d["id"]] = {"file": p.name,
                        "inputs": {i["name"]: i for i in d.get("inputs") or []},
                        "outputs": {o["name"]: o for o in d.get("outputs") or []},
                        "cmd": d.get("shell_command", "")}

    problems, notes, skips, unchecked = [], [], [], []

    def bad(kind, where, msg):
        problems.append((kind, where, msg))

    # -- 1. every `in:` reference resolves --------------------------------------------------
    for name, s in steps.items():
        for port, raw in (s.get("in") or {}).items():
          for ref in sources_of(raw):
            if "/" not in ref:
                if ref not in wf_in:
                    bad("BROKEN-LINK", f"{name}.{port}", f"{ref!r} is not a workflow input")
                continue
            src, out = ref.split("/", 1)
            if src not in steps:
                bad("BROKEN-LINK", f"{name}.{port}", f"names step {src!r}, which does not exist")
            elif steps[src].get("tool_id") in udt \
                    and out not in udt[steps[src]["tool_id"]]["outputs"]:
                bad("BAD-OUTPUT", f"{name}.{port}",
                    f"{src} has no output {out!r}; it declares "
                    f"{sorted(udt[steps[src]['tool_id']]['outputs'])}")

    # -- 2. workflow outputs resolve --------------------------------------------------------
    for name, ref in wf_out.items():
        ref = ref.get("outputSource", ref) if isinstance(ref, dict) else ref
        src, out = ref.split("/", 1)
        if src not in steps:
            bad("BROKEN-OUTPUT", f"outputs.{name}", f"names step {src!r}, which does not exist")
        elif steps[src].get("tool_id") in udt \
                and out not in udt[steps[src]["tool_id"]]["outputs"]:
            bad("BAD-OUTPUT", f"outputs.{name}",
                f"{src} has no output {out!r}; it declares "
                f"{sorted(udt[steps[src]['tool_id']]['outputs'])}")

    # -- 3. UDT steps: declared inputs all supplied, no unknown keys ------------------------
    for name, s in steps.items():
        tid = s.get("tool_id", "")
        if tid not in udt:
            continue
        given = set(s.get("in") or {}) | set(s.get("state") or {})
        declared = set(udt[tid]["inputs"])
        for m in sorted(declared - given):
            bad("UNSUPPLIED", f"{name}.{m}", f"{tid} declares input {m!r}; step supplies nothing")
        for e in sorted(given - declared):
            bad("UNKNOWN-PARAM", f"{name}.{e}",
                f"{tid} declares no input {e!r}; it has {sorted(declared)}")

    # -- 4. every UDT this workflow uses, linted as a TOOL -----------------------------------
    #    ⛔ TWO COPIES OF A CHECK DRIFT, AND THE WEAKER ONE IS THE ONE THAT RUNS. This script used
    #    to carry its own phantom-output and nested-`$(` tests. `scripts/check_udt_definitions.py`
    #    then grew the same two, better -- and an adversarial pass found a case the copy HERE let
    #    through (`from_work_dir: anchor.isoforms.tsv`, whose stem `anchor` appears in the command,
    #    so the test passed while the file was written only conditionally and often not at all).
    #    There is now one implementation, and the workflow half calls it for the tools a run will
    #    actually use. Everything it reports -- schema, templating, silent success -- lands here.
    for tid, info in sorted(udt.items()):
        if tid not in {s.get("tool_id") for s in steps.values()}:
            continue
        tp, tn = udt_lint(ROOT / "udt" / info["file"])
        for code, wh, msg in tp:
            bad(code, wh, msg)
        notes += [f"{code:<9} {wh} — {msg}" for code, wh, msg in tn]

    # -- 4a. a handoff tag must carry the `name:` prefix -------------------------------------
    #    ⛔ ONLY `name:` TAGS PROPAGATE. `AUTO_PROPAGATED_TAGS = ["name"]`
    #    (lib/galaxy/model/__init__.py) is the only tag set Galaxy copies onto a job's outputs
    #    (tools/actions/__init__.py). A plain tag therefore sticks to the dataset it was set on and
    #    is GONE from everything derived from it downstream -- and it does not render as a
    #    #hashtag either. Every tag in this repository exists to let the NEXT workflow find its
    #    inputs by tag, so a missing prefix silently turns a handoff into a manual search. WF-A's
    #    UDT edition restored all five of the classic's tag actions without the prefix and read as
    #    full parity to a checker that compared only the action KEYS.
    for name, s in steps.items():
        for act, spec in (s.get("post_job_actions") or {}).items():
            if (spec or {}).get("action_type") != "TagDatasetAction":
                continue
            for tag in str((spec.get("action_arguments") or {}).get("tags") or "").split(","):
                if tag.strip() and not tag.strip().startswith(("name:", "#")):
                    # ⚠ A NOTE, BECAUSE A PLAIN TAG IS A LEGITIMATE TERMINAL LABEL and half the
                    # tags in this repository are deliberately that -- `SoftmaskedFasta` marks the
                    # dataset it is on and is not meant to travel. It is a PROBLEM only when the
                    # workflow being ported used `name:` and this one does not, which is what
                    # PARITY-DRIFT below is for. Say what the prefix does; do not fail the run.
                    notes.append(f"PLAIN-TAG {name}.{act}: tag {tag.strip()!r} has no `name:` "
                                 f"prefix, so Galaxy propagates it to NOTHING derived from this "
                                 f"dataset (AUTO_PROPAGATED_TAGS = ['name']) and it renders as an "
                                 f"ordinary tag, not a #hashtag. Deliberate?")

    # -- 4b. a connection into a repeat needs that repeat instance DECLARED ------------------
    #    ⛔ THIS CHECK EXISTS BECAUSE THE ABSENCE OF ONE SHIPPED A BROKEN WORKFLOW PAST IT. WF-A's
    #    UDT edition connected `results_0|software_cond|input` while declaring no `state: results:`
    #    block at all. Every key it DID declare was valid, so this script reported a clean pass --
    #    and a repeat has no default instance to fall back on (run_udt_workflow's `_fill` skips
    #    repeats explicitly), so MultiQC would have run with its `software` test-param unset and
    #    published an empty report. Validating what is PRESENT says nothing about what is MISSING.
    for name, s in steps.items():
        declared = s.get("state") or {}
        for key in (s.get("in") or {}):
            m = re.match(r"([A-Za-z_][A-Za-z0-9_]*)_(\d+)\|", key)
            if not m:
                continue
            repeat, idx = m.group(1), int(m.group(2))
            have = declared.get(repeat)
            if not isinstance(have, list) or len(have) <= idx:
                bad("REPEAT-UNDECLARED", f"{name}.{key}",
                    f"connects to instance {idx} of the `{repeat}` repeat, but `state.{repeat}` "
                    f"declares {0 if not isinstance(have, list) else len(have)} instance(s). "
                    f"A repeat has no default instance, so the connected value lands on an "
                    f"otherwise-unset entry and the step runs with its own parameters missing.")

    # -- 4c. a *_udt port must not silently lose anything its classic sibling has -------------
    #    ⛔ A PORT THAT QUIETLY DROPS THINGS IS THE FAILURE THIS REPOSITORY KEEPS HITTING. The WF-A
    #    UDT edition lost all five of the classic's post_job_actions -- the renames and tags that
    #    ARE the handoff to WF-C and WF-C2 -- and nothing noticed, because every port it did declare
    #    was valid. An operator then picks collections out of a history by position.
    #
    #    ⚠ A DELIBERATE OMISSION IS FINE AND MUST BE SAID SO. Anything named in the port's own
    #    `doc:` is accepted: this edition drops `signatures` on purpose and explains why. The rule
    #    is not "lose nothing", it is "lose nothing SILENTLY".
    if args.workflow.name.endswith("_udt.gxwf.yml"):
        classic = args.workflow.with_name(args.workflow.name.replace("_udt.gxwf.yml",
                                                                    ".gxwf.yml"))
        if not classic.exists():
            notes.append(f"NO-CLASSIC no sibling {classic.name}; parity not checked")
        else:
            cw = yaml.safe_load(classic.read_text())
            doc = str(wf.get("doc") or "")
            for kind, mine, theirs in (
                    ("input", set(wf_in), set(cw.get("inputs") or {})),
                    ("output", set(wf_out), set(cw.get("outputs") or {}))):
                for lost in sorted(theirs - mine):
                    # ⚠ BACKTICKED, NOT A BARE SUBSTRING. `lost in doc` passed for `signatures`
                    # purely because the word appears in unrelated prose, so deleting the actual
                    # explanation still read as "explained". An omission is documented when the
                    # NAME is named.
                    if f"`{lost}`" in doc:
                        notes.append(f"DROPPED   {kind} `{lost}` — absent here, explained in doc:")
                    else:
                        bad("PARITY-LOSS", f"{kind} {lost}",
                            f"`{classic.name}` declares it and this port does not, and the port's "
                            f"`doc:` never mentions it. Say why, or restore it.")
            for label, cs in (cw.get("steps") or {}).items():
                theirs = set(cs.get("post_job_actions") or {})
                if not theirs:
                    continue
                mine = set((steps.get(label) or {}).get("post_job_actions") or {})
                if label not in steps:
                    continue
                for lost in sorted(theirs - mine):
                    act = (cs["post_job_actions"][lost].get("action_arguments") or {})
                    what = act.get("newname") or act.get("tags") or lost
                    bad("PARITY-LOSS", f"{label}.{lost}",
                        f"the classic sets this post_job_action ({what!r}) and the port does not; "
                        f"downstream staging identifies these datasets by name/tag.")
                # ⛔ AND THE ARGUMENTS, NOT JUST THE KEYS. Comparing action NAMES said the WF-A
                # port had every one of the classic's five renames and tags. It did -- with two
                # of the five renamed strings rewritten ('WF-C2 input: anchor BED12' for
                # 'WF-C2 input: anchor_bed12s') and all five tags stripped of `name:`. The keys
                # are the bookkeeping; the arguments are the thing that has to match.
                for kept in sorted(theirs & mine):
                    was = (cs["post_job_actions"][kept] or {}).get("action_arguments") or {}
                    now = ((steps[label].get("post_job_actions") or {})[kept]
                           or {}).get("action_arguments") or {}
                    for k in sorted(set(was) & set(now)):
                        if was[k] != now[k] and f"`{now[k]}`" not in doc:
                            bad("PARITY-DRIFT", f"{label}.{kept}.{k}",
                                f"the classic sets {was[k]!r} and this port sets {now[k]!r}. These "
                                f"strings ARE the handoff -- the next workflow stages by this "
                                f"exact name or tag. Match it, or say in `doc:` why it differs.")

    # -- 5. dead ends and unreachable steps -------------------------------------------------
    consumed = {x for s in steps.values() for r in (s.get("in") or {}).values()
                for x in sources_of(r)}
    consumed |= {(r.get("outputSource", r) if isinstance(r, dict) else r)
                 for r in wf_out.values()}
    for name, s in steps.items():
        if s.get("tool_id") in udt:
            for o in udt[s["tool_id"]]["outputs"]:
                if f"{name}/{o}" not in consumed:
                    notes.append(f"UNUSED    {name}/{o} — declared, consumed by nothing")

    reachable, frontier = set(), {n for n, s in steps.items()
                                  if any("/" not in x
                                         for r in (s.get("in") or {}).values()
                                         for x in sources_of(r))}
    while frontier:
        n = frontier.pop()
        if n in reachable:
            continue
        reachable.add(n)
        frontier |= {m for m, s in steps.items()
                     if any(x.split("/")[0] == n
                            for r in (s.get("in") or {}).values() for x in sources_of(r))}
    for n in steps:
        if n in reachable:
            continue
        # ⚠ A STEP WITH NO DATASET CONNECTION IS A ROOT, NOT AN ORPHAN. `ucsc_hub`'s `hub_check`
        # declares `hub_url: {default: hub.txt}` -- a parameter, not an edge -- so nothing can
        # reach it and nothing needs to. Calling that UNREACHABLE made a correct workflow fail.
        # Only a step that HAS connections and still cannot be reached is actually stranded.
        if not any(sources_of(r) for r in (steps[n].get("in") or {}).values()):
            notes.append(f"ROOT      {n} — takes no dataset connection, runs on parameters alone")
        else:
            bad("UNREACHABLE", n, "has connections, but none trace back to a workflow input")

    # -- 6. version pinning -----------------------------------------------------------------
    #    ⚠ A BARE ToolShed ID IS NOT PINNED. It resolves to whatever that Galaxy has installed
    #    NEWEST at import time. That is a live hazard for any step carrying a hand-written
    #    `state:` block, because a different wrapper revision can rename a parameter -- and the
    #    state blocks in this repository exist precisely to avoid allow_tool_state_corrections.
    unpinned = [(n, s["tool_id"]) for n, s in steps.items()
                if s.get("tool_id") not in udt and "/" not in s.get("tool_id", "")
                and s.get("tool_id") != "cat1"]          # cat1 is built-in; it has no other form

    # -- 7. remote parameters AND remote OUTPUT NAMES ---------------------------------------
    #    ⚠ THE OUTPUT HALF IS NOT OPTIONAL. An earlier version checked output names only for the
    #    UDTs, because those are files in this repo -- so `maskfasta/absent`, naming a ToolShed
    #    tool, passed clean. The live schema lists a tool's outputs as surely as its inputs, and a
    #    reference to one that does not exist is the single likeliest hand-editing error.
    url, key = os.environ.get("GALAXY_URL", "").rstrip("/"), os.environ.get("GALAXY_API_KEY", "")
    remote = [(n, s) for n, s in steps.items() if s.get("tool_id") not in udt]
    remote_params = {}
    if not (url and key):
        skips.append(f"REMOTE HALF SKIPPED — {len(remote)} non-UDT step(s) unchecked, including "
                     f"every reference to one of their outputs. "
                     f"Set GALAXY_URL and GALAXY_API_KEY to check them. NOT a pass.")
    else:
        remote_outs = {}
        for name, s in remote:
            try:
                sch = tool_schema(s["tool_id"], url, key, use_cache=not args.no_cache)
            except Exception:  # noqa: BLE001  -- reported by the parameter pass just below
                continue
            remote_outs[name] = {o.get("name") for o in sch.get("outputs") or []}

        refs = [(f"{n}.{port}", x) for n, s in steps.items()
                for port, r in (s.get("in") or {}).items() for x in sources_of(r)]
        refs += [(f"outputs.{n}", (r.get("outputSource", r) if isinstance(r, dict) else r))
                 for n, r in wf_out.items()]
        for where, ref in refs:
            if "/" not in ref:
                continue
            src, out = ref.split("/", 1)
            if src in remote_outs and out not in remote_outs[src]:
                bad("BAD-OUTPUT", where,
                    f"{src} ({steps[src]['tool_id'].split('/')[-2] if '/' in steps[src]['tool_id'] else steps[src]['tool_id']}) "
                    f"has no output {out!r}; it declares {sorted(remote_outs[src])}")

        for name, s in remote:
            tid = s["tool_id"]
            try:
                sch = tool_schema(tid, url, key, use_cache=not args.no_cache)
            # ⛔ A TOOL THIS GALAXY DOES NOT HAVE IS AN UNCHECKED STEP, NOT A BROKEN ONE, and
            # conflating the two made this script report every non-UDT workflow in the repo as
            # defective: `workflows/softmask/softmask.gxwf.yml` names `dustmasker`, `tantan`,
            # `fastan` and `masking_table`, which are installed on the BRC instance and 404 on
            # usegalaxy.org -- which is the entire reason the UDT edition exists. Saying "6
            # PROBLEMS" there is a statement about the SERVER I asked, not about the workflow.
            # It withholds, and names the server, so the reader can point it at the right one.
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    unchecked.append(f"{name} ({tid})")
                else:
                    bad("UNRESOLVED", name, f"{tid} -> HTTP {e.code}")
                continue
            # ⚠ Anything else -- a timeout, an auth failure, malformed JSON -- is still REPORTED,
            # never swallowed, and never allowed to abort the remaining checks.
            except Exception as e:  # noqa: BLE001
                bad("UNRESOLVED", name, f"{tid} -> {type(e).__name__}: {e}")
                continue
            params = flatten(sch.get("inputs"))
            remote_params[name] = params
            for key_ in list(s.get("in") or {}) + list(s.get("state") or {}):
                if not addressable(key_, params):
                    bad("NO-SUCH-PARAM", f"{name}.{key_}",
                        f"{tid.split('/')[-2] if '/' in tid else tid} v{sch.get('version')} "
                        f"has no such parameter")
            for n2, t2 in unpinned:
                if n2 == name:
                    vs = sch.get("versions") or []
                    notes.append(f"UNPINNED  {name}: tool_id {t2!r} carries no version; "
                                 f"{url} has {len(vs)} installed and resolves it to "
                                 f"v{sch.get('version')}")

    # -- 8. an EMPTY collection input schedules, runs nothing, and succeeds -------------------
    #    ⛔ ZERO JOBS IS NOT AN ERROR TO GALAXY. A step that maps over a collection with no
    #    elements produces no jobs at all; its output collection is empty, everything downstream of
    #    it is empty, and the invocation reaches `scheduled` with every job `ok` because there were
    #    none to fail. Measured on WF-A: an empty `anchor_gene_gff3s` gives empty `anchor_bed12s`
    #    and `anchor_isoforms` and a green run. The ONLY thing that can notice is a step that
    #    receives the collection WHOLE -- a `data_collection` port, or a remote parameter that
    #    accepts one -- because only that step is entered even when there is nothing in it.
    #
    #    ⚠ THIS IS A NOTE. Mapping over a collection is the normal shape of every workflow here,
    #    and "an operator supplied an empty input" is not a defect in the workflow. What the note
    #    buys is knowing WHICH inputs have nothing standing between them and a silent no-op.
    coll_in = {n for n, spec in wf_in.items()
               if isinstance(spec, dict) and spec.get("type") == "collection"}

    def consumers(seed):
        """(step, port) pairs taking `seed` DIRECTLY -- not through a derived collection."""
        return [(n, port) for n, s in steps.items()
                for port, r in (s.get("in") or {}).items() if seed in sources_of(r)]

    for seed in sorted(coll_in):
        whole, unknown = [], []
        for n, port in consumers(seed):
            tid = steps[n].get("tool_id", "")
            if tid in udt:
                spec = udt[tid]["inputs"].get(port.split("|")[0]) or {}
                if spec.get("type") == "data_collection":
                    whole.append(f"{n}.{port}")
            elif n in remote_params:
                param = remote_params[n].get(port) or {}
                if param.get("type") == "data_collection" or param.get("multiple"):
                    whole.append(f"{n}.{port}")
            else:
                unknown.append(f"{n}.{port}")           # remote half skipped: cannot say
        if whole:
            continue
        tail = (f" ({len(unknown)} consumer(s) are non-UDT steps this run did not check: "
                f"{', '.join(sorted(unknown))})" if unknown else "")
        notes.append(f"EMPTY-COLLECTION `{seed}` is consumed element-wise ONLY, so an empty one "
                     f"schedules zero jobs, publishes empty collections and returns a SUCCESSFUL "
                     f"invocation. Nothing in this workflow sees it whole.{tail}")

    # -- 9. two collections nothing ever compares ---------------------------------------------
    #    ⛔ A CORRECT-LOOKING ANCHOR FOR THE WRONG GENOME PASSES EVERYTHING. WF-A takes
    #    `assemblies` and `anchor_gene_gff3s` and hands both to WF-C/TOGA2, and `anchor_prep` never
    #    sees an assembly: an anchor GFF3 for a different strain produces a fully populated BED12,
    #    twelve green outputs and exit 0. The two facts that make it undetectable are STRUCTURAL
    #    and visible from here -- the two inputs' downstream cones never meet at a step, so there
    #    is no point in the graph at which anything COULD compare them.
    #
    #    ⚠ AND A JOIN IS NOT A COMPARISON, so this is a floor and not a guarantee. WF-A's
    #    `proteomes` and `assemblies` DO meet, at multiqc -- which merges a BUSCO summary keyed by
    #    one collection's identifiers with plots keyed by the other's and compares nothing. This
    #    check can prove that nothing is able to compare them; it cannot prove that something does.
    # ⚠ A FIXED POINT, NOT ONE PASS. `steps` is a mapping and its order is the order somebody
    # typed it, so a single sweep propagates a seed only as far as the next step that happens to
    # come later in the file. The first version of this reported WF-A's `assemblies` and
    # `proteomes` as never meeting -- they meet at multiqc, four steps downstream of both -- which
    # is the check inventing the finding it was written to look for. Iterate until nothing moves.
    feeds = {n: set() for n in steps}
    while True:
        before = {n: set(v) for n, v in feeds.items()}
        for n, s in steps.items():
            for r in (s.get("in") or {}).values():
                for ref in sources_of(r):
                    src = ref.split("/")[0] if "/" in ref else ref
                    feeds[n] |= ({src} & coll_in) | feeds.get(src, set())
        if feeds == before:
            break
    joined = set()
    for seeds in feeds.values():
        for a in seeds:
            for b in seeds:
                if a != b:
                    joined.add(tuple(sorted((a, b))))
    for a, b in sorted({tuple(sorted((x, y))) for x in coll_in for y in coll_in if x != y}
                       - joined):
        notes.append(f"UNRECONCILED `{a}` and `{b}` are both collection inputs and their "
                     f"downstream steps never meet, so NOTHING in this workflow is in a position "
                     f"to compare them -- not their element identifiers, not the sequence names "
                     f"inside them. A file for the wrong genome, or a panel and an anchor set that "
                     f"disagree, runs green.")

    # -- report -----------------------------------------------------------------------------
    print(f"  {args.workflow.name}: {len(steps)} steps, {len(wf_in)} input(s), "
          f"{len(wf_out)} output(s)")
    print(f"  checkable offline: {sum(1 for s in steps.values() if s.get('tool_id') in udt)}"
          f"  ·  needs a live Galaxy: {len(remote)}\n")
    for s in skips:
        print(f"  ⚠ SKIP  {s}")
    if unchecked:
        print(f"  ⚠ SKIP  {len(unchecked)} step(s) NOT CHECKED — {url} has no such tool "
              f"(installed elsewhere?): {', '.join(unchecked)}")
    if problems:
        print(f"  ⛔ {len(problems)} PROBLEM(S)")
        for k, w, m in problems:
            print(f"    {k:<15} {w:<34} {m}")
    elif not (skips or unchecked):
        print("  ✅ every port resolves and every parameter names a real tool input")
    for n in sorted(set(notes)):
        print(f"    {n}")
    return 1 if problems or ((skips or unchecked) and args.strict) else 0


if __name__ == "__main__":
    sys.exit(main())
