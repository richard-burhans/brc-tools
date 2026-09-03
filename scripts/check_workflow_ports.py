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
import sys
import urllib.error
import urllib.parse
import urllib.request

import yaml

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
        tmp = pathlib.Path(td) / "case.gxwf.yml"
        for rel, label, old, new, want in DEFECTS:
            src = (ROOT / rel).read_text()
            if src.count(old) < 1:
                print(f"  BROKEN-FIXTURE  {label}: {old!r} no longer occurs in {rel}")
                failed += 1
                continue
            tmp.write_text(src.replace(old, new, 1))
            r = subprocess.run([*me, str(tmp)], capture_output=True, text=True, env=env)
            hit = want in r.stdout
            print(f"  {'pass' if hit else 'FAIL'}  {label:<38} expect {want}")
            failed += not hit
        clean = ROOT / "workflows/softmask/softmask_udt.gxwf.yml"
        tmp.write_text(clean.read_text())
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
            ref = source_of(raw)
            if ref is None:
                continue
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

    # -- 4. a UDT output must name a file its command writes --------------------------------
    #    ⚠ SUBSTRING, AND DELIBERATELY LOOSE. `FasTAN -oscan` writes `scan.1ano` without the
    #    literal string ever appearing, so an exact match reports a phantom on a correct tool.
    #    Checking the STEM catches the real defect (a renamed output file) without that.
    for _tid, t in udt.items():
        for oname, o in t["outputs"].items():
            fwd = o.get("from_work_dir")
            if fwd and fwd.split(".")[0] not in t["cmd"]:
                bad("PHANTOM-OUTPUT", f"{t['file']}:{oname}",
                    f"from_work_dir {fwd!r}: no sign of {fwd.split('.')[0]!r} in shell_command")

    # -- 5. dead ends and unreachable steps -------------------------------------------------
    consumed = {source_of(r) for s in steps.values() for r in (s.get("in") or {}).values()
                if source_of(r)}
    consumed |= {(r.get("outputSource", r) if isinstance(r, dict) else r)
                 for r in wf_out.values()}
    for name, s in steps.items():
        if s.get("tool_id") in udt:
            for o in udt[s["tool_id"]]["outputs"]:
                if f"{name}/{o}" not in consumed:
                    notes.append(f"UNUSED    {name}/{o} — declared, consumed by nothing")

    reachable, frontier = set(), {n for n, s in steps.items()
                                  if any((source_of(r) or "/") .find("/") < 0
                                         for r in (s.get("in") or {}).values())}
    while frontier:
        n = frontier.pop()
        if n in reachable:
            continue
        reachable.add(n)
        frontier |= {m for m, s in steps.items()
                     if any((source_of(r) or "").split("/")[0] == n
                            for r in (s.get("in") or {}).values())}
    for n in steps:
        if n in reachable:
            continue
        # ⚠ A STEP WITH NO DATASET CONNECTION IS A ROOT, NOT AN ORPHAN. `ucsc_hub`'s `hub_check`
        # declares `hub_url: {default: hub.txt}` -- a parameter, not an edge -- so nothing can
        # reach it and nothing needs to. Calling that UNREACHABLE made a correct workflow fail.
        # Only a step that HAS connections and still cannot be reached is actually stranded.
        if not any(source_of(r) for r in (steps[n].get("in") or {}).values()):
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

        refs = [(f"{n}.{port}", source_of(r)) for n, s in steps.items()
                for port, r in (s.get("in") or {}).items() if source_of(r)]
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
