#!/usr/bin/env python3
"""Static lint for a Galaxy User-Defined Tool: the `udt-authoring` skill, mechanised.

    python3 scripts/check_udt_definitions.py                    # every udt/*.gxtool.yml
    python3 scripts/check_udt_definitions.py udt/anchor_prep.gxtool.yml
    python3 scripts/check_udt_definitions.py --self-test        # prove each check goes red

⛔ A CHECKLIST A HUMAN RUNS IS A CHECKLIST A HUMAN SKIPS. The `udt-authoring` skill ships a
27-row "common mistakes" table and a 13-item pre-submit checklist, and this repository has shipped
five of those mistakes anyway -- three of them in tools that were then invoked on a live Galaxy and
had NEVER been capable of running. The rules were written down the whole time. What was missing was
something that reads them on every file, every time, without being asked.

WHAT THIS IS NOT. It is not `galaxy-tool-util`. The skill's own `scripts/validate.py` runs the
server's real `UserToolSource` model and `lint_user_tool_source`, and where that is installed it is
the better authority for the SCHEMA half below -- run both. This one needs nothing but PyYAML, so
it runs in any checkout and in CI, and it goes further than the schema in the two directions that
actually cost this project time:

  * TEMPLATING -- what `$(...)` means, which the schema does not model at all. Galaxy owns that
    sequence and interpolates it ANYWHERE in a `shell_command`, heredocs included. Every finding
    here was measured on a real job, not inferred.
  * SILENT SUCCESS -- a tool that exits 0 having produced nothing, or the wrong thing. No schema
    can see these; they were found by attacking the tools with adversarial inputs, and each one is
    named with the case that produced it.

⚠ CLEAN HERE DOES NOT MEAN THE TOOL RUNS. Nothing offline can tell you the container image exists
-- a wrong build tag is the single most common real UDT failure ("manifest unknown", at run time,
after a clean create) -- nor that the command works, nor that the server's role and job-config
gates will let it execute. This narrows the search; it does not end it.

Exit status is 1 if any PROBLEM was found. Notes never fail the run.
"""
import argparse
import itertools
import json
import os
import pathlib
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- the schema, as data ----------
#: `UserToolSource` is `extra="forbid"`, so an unknown field is a parse error rather than a warning.
TOP_LEVEL = {"class", "id", "version", "name", "description", "container", "shell_command",
             "inputs", "outputs", "requirements", "configfiles", "citations", "license", "help",
             "tests", "profile", "edam_operations", "edam_topics", "xrefs"}

#: XML/Cheetah habits that leak in. The VALUE is what to do instead -- a checker that only says
#: "rejected" leaves the author to go looking for the replacement.
XML_FIELDS = {
    "command": "the field is `shell_command`",
    "truevalue": "use `value: false` + a ternary in shell_command: $(inputs.x ? '--x' : '')",
    "falsevalue": "use `value: false` + a ternary in shell_command",
    "argument": "not a UDT field; name the flag in shell_command yourself",
    "is_dynamic": "not a UDT field",
    "parameter_type": "the field is `type`",
    "hidden": "hidden parameters do not exist in a UDT",
}
LEAF_TYPES = {"boolean", "integer", "float", "text", "select", "color", "data", "data_collection"}
GROUP_TYPES = {"conditional", "repeat", "section"}
#: Present in XML, rejected by the UDT parser -- each of these is a parse error, not a warning.
XML_TYPES = {"hidden", "drill_down", "data_column", "genomebuild", "group_tag", "baseurl",
             "rules", "directory"}
SCALAR_TYPES = {"boolean", "integer", "float", "text", "select", "color"}
#: What a `data` input actually renders as -- a CWL-style File record. Anything else after the dot
#: evaluates to `undefined`, which interpolates as the literal string and the command then runs on
#: a path that does not exist. `$(inputs.reads.pth)` is the whole failure, and testing only that
#: SOMETHING followed the dot let it through.
FILE_FIELDS = {"path", "basename", "nameroot", "nameext", "class", "location", "size", "format",
               "checksum", "secondaryFiles", "contents", "dirname"}
#: A `pattern` discovery entry has no defaults in the model: all nine or it fails validation.
DISCOVER_KEYS = {"discover_via", "pattern", "directory", "format", "visible", "recurse",
                 "match_relative_path", "assign_primary_output", "sort_key", "sort_comp"}
HELP_FORMATS = {"markdown", "restructuredtext", "plain_text"}
#: ⚠ THESE ARE TAGGED BY BARE VERSION AND THE `--<build>` SMELL DOES NOT APPLY. A biocontainers
#: image for a TOOL carries a build suffix, so a tag without one is usually a guess -- but the
#: base-language images are not built that way, and `quay.io/biocontainers/python:3.12` is what
#: seven tools in this repository have actually run on. A note that fires on a verified image
#: seven times is how a reader learns to skip the notes.
PLAIN_TAG_OK = {"python", "perl", "r-base", "busybox", "bash", "gawk"}
#: Flags whose ARGUMENT is a thread count. Recording $GALAXY_SLOTS while passing a literal here is
#: the mistake the skill calls out by name: "recording the value isn't using it".
THREAD_FLAGS = r"(?:--threads|--cpus|--num[-_]threads|--nproc|--jobs|-@|-p|-t|-j)"


#: Where a container check's answers are remembered, OUTSIDE the working tree -- the same reason
#: check_workflow_ports.py keeps its tool-schema cache there: a cache inside the repo gets
#: committed by accident exactly once, and then needs a .gitignore line forever.
CACHE = pathlib.Path(
    os.environ.get("XDG_CACHE_HOME") or (pathlib.Path.home() / ".cache")) / "brc-tools" / "containers.json"


def container_exists(image, timeout=20):
    """Does this image actually exist? (True/False/None for "could not ask").

    ⛔ THE ONE FAILURE EVERY OTHER CHECK HERE IS BLIND TO. Nothing in the schema, in
    galaxy-tool-util's lint, or in Galaxy's own create step resolves a container image: a UDT
    naming a tag that does not exist REGISTERS cleanly and then dies at job start with
    `manifest unknown`. The udt-authoring skill calls guessing the `--<hash>_<build>` suffix the
    single most common real UDT failure -- and this repository shipped exactly that, in
    `quay.io/biocontainers/tantan:51--h4ac6f70_0`, which meant WF-B's UDT edition had never been
    able to run. This file NAMED that hazard in a comment and did not check it. One HTTP request
    closes it.

    ⚠ AND A FAILED REQUEST IS `None`, NOT `False`. A proxy refusal or a rate limit is a statement
    about this network, not about the image, and reporting it as a missing container would be the
    same laundering of a self-inflicted fetch failure into a fact about the world that has bitten
    this project before.
    """
    if not image.startswith("quay.io/biocontainers/"):
        return None
    repo_tag = image.split("quay.io/biocontainers/", 1)[1]
    if ":" not in repo_tag:
        return None
    name, tag = repo_tag.rsplit(":", 1)
    cache = {}
    if CACHE.exists():
        try:
            cache = json.loads(CACHE.read_text())
        except ValueError:
            cache = {}
    if image in cache:
        return cache[image]
    # depot.galaxyproject.org mirrors every biocontainer and answers a HEAD cheaply; quay's tag
    # API is the second opinion, because depot lags a very new build by a few hours.
    answer = None
    try:
        req = urllib.request.Request(
            f"https://depot.galaxyproject.org/singularity/{urllib.parse.quote(name)}%3A{urllib.parse.quote(tag)}",
            method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            answer = r.status == 200
    except urllib.error.HTTPError as e:
        answer = False if e.code == 404 else None
    except Exception:                                   # noqa: BLE001 -- network, not the image
        answer = None
    if answer is False:
        try:
            url = (f"https://quay.io/api/v1/repository/biocontainers/{urllib.parse.quote(name)}"
                   f"/tag/?onlyActiveTags=true&specificTag={urllib.parse.quote(tag)}")
            with urllib.request.urlopen(url, timeout=timeout) as r:
                answer = bool(json.load(r).get("tags"))
        except Exception:                               # noqa: BLE001
            pass
    if answer is not None:
        cache[image] = answer
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    return answer


def walk_inputs(nodes, prefix=""):
    """Every input node, structural groups included, as (addressed_name, node).

    ⚠ A GROUP IS NOT A LEAF AND ITS CHILDREN ARE NOT TOP-LEVEL. `conditional` nests under
    `whens`, `repeat`/`section` under `parameters`, and a `$(inputs.X)` reference resolves against
    the TOP-LEVEL name only -- a conditional's test parameter is reached as `inputs.cond.…`, never
    as `inputs.test_param`. Flattening the tree into one namespace would accept a reference the
    server rejects, which is the wrong direction for a checker to be wrong in.
    """
    for n in nodes or []:
        if not isinstance(n, dict):
            continue
        name = f"{prefix}{n.get('name', '?')}"
        yield name, n
        t = n.get("type")
        if t == "conditional":
            tp = n.get("test_parameter")
            if isinstance(tp, dict):
                yield f"{name}.{tp.get('name', '?')}", tp
            for w in n.get("whens") or []:
                yield from walk_inputs((w or {}).get("parameters"), f"{name}.")
        elif t in {"repeat", "section"}:
            yield from walk_inputs(n.get("parameters"), f"{name}.")


def expressions(cmd):
    """Every unescaped `$(...)` in the command, as (line, expr, start, end).

    ⚠ `\\$(...)` IS THE DOCUMENTED ESCAPE AND IS NOT AN EXPRESSION. Galaxy passes an escaped form
    to the shell untouched, so `\\$(date)` is correct code -- the skill lists exactly that as the
    FIX for an unescaped literal. A scanner that flags it punishes the idiom it should be teaching.

    ⚠ AND IT MUST COUNT DEPTH, NOT MATCH A REGEX. `\\$\\((?:[^()]*\\$\\()` requires the two openers to
    have no parenthesis between them, so `$(tr -d '()' < '$(inputs.x.path)')` -- a real construct
    from a real tool that never ran -- slips straight through it.
    """
    i = 0
    while (i := cmd.find("$(", i)) != -1:
        if i and cmd[i - 1] == "\\":
            i += 2
            continue
        depth, j = 0, i + 1
        while j < len(cmd):
            if cmd[j] == "(":
                depth += 1
            elif cmd[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        yield cmd[:i].count("\n") + 1, cmd[i + 2:j], i, j
        i = j + 1


def statements(cmd):
    """The command's top-level statements, as (line, text) -- heredoc bodies and quotes skipped.

    ⛔ EVERY UDT IN THIS REPOSITORY IS MOSTLY NOT SHELL. They write an awk or Python program into
    the work dir with a heredoc and then run it, so a naive line split reads a few hundred lines of
    another language as shell statements and every check built on it becomes noise. This tracks
    heredoc terminators, single/double quotes (an awk program is one multi-line quoted string) and
    `\\`-continuations, and yields only what the SHELL would treat as a statement of its own.
    """
    out, i, n = [], 0, len(cmd)
    quote = heredoc = pending = None
    start, line = 0, 1
    while i < n:
        ch = cmd[i]
        if heredoc is not None:                      # inside a heredoc body: consume whole lines
            eol = cmd.find("\n", i)
            eol = n if eol == -1 else eol
            if cmd[i:eol].strip() == heredoc:
                heredoc = None
            i, start, line = eol + 1, eol + 1, line + 1
            continue
        if quote:
            if ch == "\\" and quote == '"' and i + 1 < n:
                i += 2
                continue
            if ch == quote:
                quote = None
            line += ch == "\n"
            i += 1
            continue
        if ch in "'\"":
            quote, i = ch, i + 1
            continue
        if ch == "\\" and i + 1 < n:                 # a continuation joins the two lines
            line += cmd[i + 1] == "\n"
            i += 2
            continue
        if cmd.startswith("<<", i):
            m = re.match(r"<<-?\s*(['\"]?)([A-Za-z_][A-Za-z0-9_]*)\1", cmd[i:])
            if m:
                pending = m.group(2)
                i += m.end()
                continue
        if ch == "\n":
            out.append((line, cmd[start:i]))
            i, start, line = i + 1, i + 1, line + 1
            heredoc, pending = pending, None
            continue
        i += 1
    if start < n:
        out.append((line, cmd[start:]))
    return out


# ---------------------------------------------------------------- the checks -------------------
def authoritative(path, data):
    """The SERVER's own verdict, when `galaxy-tool-util` is installed. (findings, ran?).

    ⛔ THIS IS THE REAL AUTHORITY FOR THE SCHEMA HALF AND THIS SCRIPT IS NOT. `UserToolSource` and
    `lint_user_tool_source` are the exact code a Galaxy runs on create, so where they are available
    they settle every question about structure -- and they will keep doing so through schema
    changes this file has never heard of. Everything re-implemented above exists to run in a
    checkout that does NOT have them, and to cover the two halves they do not model at all.

    ⚠ MEASURED, so that neither is mistaken for the other: over this script's own 16 injected
    defects, 8 are caught by BOTH, 8 by this script ALONE (templating, container pinning, phantom
    and asymmetric outputs) and 0 by galaxy-tool-util alone. Run both; neither is a superset.
    """
    try:
        from galaxy.tool_util.lint import lint_user_tool_source
        from galaxy.tool_util_models import UserToolSource, format_validation_errors
        from pydantic import ValidationError
    except ImportError:
        return [], False
    try:
        tool = UserToolSource.model_validate(data)
    except ValidationError as exc:
        return [("SCHEMA", path.name, b) for b in format_validation_errors(exc)], True
    return [("LINT", path.name, b) for b in lint_user_tool_source(tool)], True


def lint(path, containers=True):                    # this IS a checklist; it is long on purpose
    """Return (problems, notes), each a list of (CODE, where, message)."""
    problems, notes = [], []
    where = path.name

    def bad(code, msg, at=""):
        problems.append((code, f"{where}{at}", msg))

    def note(code, msg, at=""):
        notes.append((code, f"{where}{at}", msg))

    try:
        d = yaml.safe_load(path.read_text())
    except yaml.YAMLError as e:
        bad("UNPARSEABLE", f"not valid YAML: {e}")
        return problems, notes
    if not isinstance(d, dict):
        bad("UNPARSEABLE", "top level is not a mapping")
        return problems, notes

    # -- 1. the schema half: what `UserToolSource` itself would reject ----------------------
    if d.get("class") != "GalaxyUserTool":
        bad("CLASS", f"`class` is {d.get('class')!r}; a UDT is exactly `GalaxyUserTool`")
    for k in sorted(set(d) - TOP_LEVEL):
        bad("XMLISM", f"unknown top-level field `{k}` -- {XML_FIELDS.get(k, 'the schema is '
                                                                        'extra=forbid')}")
    tid = d.get("id")
    if tid is None:
        note("NO-ID", "no `id`; recommended so the tool is addressable by name")
    elif not re.fullmatch(r"[a-z][a-z0-9_-]{2,254}", str(tid)):
        bad("BAD-ID", f"id {tid!r} fails `^[a-z][a-z0-9_-]*$` / length 3-255 "
                      f"(surfaces as string_pattern_mismatch)")
    name = d.get("name") or ""
    if len(str(name).strip()) < 5:
        bad("SHORT-NAME", f"`name` is {name!r}; the model requires 5+ non-blank characters")
    container = d.get("container")
    if not isinstance(container, str) or not container.strip():
        bad("BLANK-CONTAINER", f"`container` must be a non-blank STRING, got {type(container).__name__}"
                               f"{' (a dict is the XML habit)' if isinstance(container, dict) else ''}")
    else:
        tag = container.rsplit(":", 1)[-1] if ":" in container.rsplit("/", 1)[-1] else ""
        if not tag:
            bad("CONTAINER-UNPINNED", f"{container!r} names no tag, so it resolves to `latest` and "
                                      f"the tool is not reproducible")
        elif tag == "latest":
            bad("CONTAINER-UNPINNED", f"{container!r} is pinned to `latest`; a rebuild silently "
                                      f"changes what ran")
        elif containers and (exists := container_exists(container)) is False:
            bad("NO-SUCH-CONTAINER",
                f"{container!r} DOES NOT EXIST -- neither depot.galaxyproject.org nor quay has "
                f"that tag. The tool will REGISTER cleanly, because creating a UDT never resolves "
                f"its image, and then every job dies at start with `manifest unknown`. Never guess "
                f"the `--<hash>_<build>` suffix; look it up.")
        elif containers and exists is None:
            note("CONTAINER-UNCHECKED",
                 f"{container!r}: could not confirm the image exists (not a biocontainer, or the "
                 f"registry was unreachable). Nothing offline can resolve an image, so this is the "
                 f"one failure a clean run here does not rule out.")
        elif ("biocontainers/" in container and "--" not in tag
              and container.rsplit("/", 1)[-1].split(":")[0] not in PLAIN_TAG_OK):
            # ⚠ NOT A PROBLEM -- A SMELL. A biocontainers tag almost always carries a
            # `--<hash>_<build>` suffix, and the skill's #1 runtime failure is an author GUESSING
            # that suffix (or omitting it). Nothing offline can confirm the tag exists; saying
            # which tags look hand-made is the most this can honestly do.
            note("CONTAINER-TAG", f"{container!r}: a biocontainers tag usually carries a "
                                  f"`--<hash>_<build>` suffix. Verify it exists -- a wrong tag "
                                  f"passes every check here and dies with `manifest unknown`.")
    if "version" in d and not str(d["version"]).strip():
        bad("BLANK-STRING", "`version` is blank (dynamic_tool.blank_string)")
    cmd = d.get("shell_command")
    if not isinstance(cmd, str) or not cmd.strip():
        bad("NO-COMMAND", "`shell_command` is missing or empty")
        cmd = ""
    desc = str(d.get("description") or "").strip()
    if not desc:
        note("NO-DESCRIPTION", "no `description`; it is the one line shown in the tool menu")
    elif re.fullmatch(r"(?i)run(s)? (the )?tool\.?|tool|todo.*", desc):
        note("VAGUE-DESCRIPTION", f"description {desc!r} does not say what the tool does")

    # -- 2. help: schema-optional, convention-required -------------------------------------
    help_ = d.get("help")
    if help_ is None:
        bad("HELP-MISSING", "no `help` block. Galaxy renders it under the parameter form, so a UDT "
                            "without one is a black box to everyone including its author later.")
    elif isinstance(help_, str):
        bad("HELP-SHAPE", "`help` is a bare string; the wire shape is an object "
                          "{format, content} and a string is rejected at create")
    elif isinstance(help_, dict):
        fmt, content = help_.get("format"), str(help_.get("content") or "")
        if fmt not in HELP_FORMATS:
            bad("HELP-SHAPE", f"help.format is {fmt!r}; expected one of {sorted(HELP_FORMATS)}")
        if len(content.strip()) < 40 or re.search(r"(?i)\btodo\b|\bfixme\b|<placeholder>", content):
            bad("HELP-STUB", "help.content is a stub or a TODO placeholder -- no better than no help")

    # -- 3. inputs ------------------------------------------------------------------------
    declared = {}
    for iname, node in walk_inputs(d.get("inputs")):
        declared[iname] = node
        for k in sorted(set(node) & set(XML_FIELDS)):
            bad("XMLISM", f"input `{iname}` carries `{k}` -- {XML_FIELDS[k]}")
        t = node.get("type")
        if t in XML_TYPES:
            bad("BAD-INPUT-TYPE", f"input `{iname}` has type {t!r}, which exists in XML tools and "
                                  f"is rejected by the UDT parser")
        elif t not in LEAF_TYPES | GROUP_TYPES:
            bad("BAD-INPUT-TYPE", f"input `{iname}` has type {t!r}; expected one of "
                                  f"{sorted(LEAF_TYPES | GROUP_TYPES)}")
        if t == "select" and not (node.get("options") or []):
            bad("SELECT-NO-OPTIONS", f"select `{iname}` declares no `options` (a UDT select is "
                                     f"static -- there is no dynamic option source)")
        for v in node.get("validators") or []:
            if not isinstance(v, dict) or v.get("type") != "regex":
                continue
            expr = str(v.get("expression") or "")
            if "regex" in v:
                bad("VALIDATOR-FIELD", f"input `{iname}`: a regex validator's pattern goes in "
                                       f"`expression`; `regex` is the XML spelling and the model "
                                       f"rejects it as an extra field")
            # ⛔ `re.match` IS ANCHORED AT THE START ONLY, and Galaxy's own
            # RegexParameterValidatorModel docstring says so: "To enforce a match of the complete
            # value use `$` at the end of the expression." Without it only the FIRST CHARACTER is
            # constrained. Measured against that model, `[A-Za-z0-9_.-]+` accepted `Pk ANKA`,
            # `a$(id)`, `../../etc/passwd` and `Pk_ANKA'; touch pwned.txt; echo '` -- and because a
            # scalar is interpolated into shell_command UNQUOTED and not shlex-quoted, that last
            # one cut a python invocation short of its --output-dir and left the job at EXIT 0
            # with all four of its from_work_dir outputs missing. A validator whose message
            # promises a character set and enforces one character is worse than none.
            elif expr and not expr.rstrip().endswith(("$", r"\Z", r"\z")):
                bad("UNANCHORED-VALIDATOR",
                    f"input `{iname}`: regex {expr!r} is applied with `re.match`, which anchors at "
                    f"the START only -- everything after the first matching run is unconstrained. "
                    f"End it with `$`.")
        if node.get("optional") and "value" not in node and t in SCALAR_TYPES:
            note("NO-DEFAULT", f"optional `{iname}` has no `value`; a user must fill it in to run "
                               f"the tool at all")
        label = str(node.get("label") or "")
        if not label:
            note("NO-LABEL", f"input `{iname}` has no `label`; the form shows the bare name")
        elif label.strip().lower().replace("_", " ") == iname.lower().replace("_", " "):
            note("LABEL-ECHOES-NAME", f"input `{iname}`'s label just repeats its name")

    top_level_inputs = {n for n in declared if "." not in n}

    # -- 4. outputs ------------------------------------------------------------------------
    fwd_targets = []
    for o in d.get("outputs") or []:
        if not isinstance(o, dict):
            bad("BAD-OUTPUT", f"output entry {o!r} is not a mapping")
            continue
        oname, t = o.get("name", "?"), o.get("type")
        if t not in {"data", "collection"}:
            bad("BAD-OUTPUT-TYPE", f"output `{oname}` has type {t!r}; a UDT supports only `data` "
                                   f"and `collection` (the scalar output types in the underlying "
                                   f"model are rejected by the YAML parser)")
        fwd, disc = o.get("from_work_dir"), o.get("discover_datasets")
        if t == "collection" and not disc:
            bad("OUTPUT-UNCLAIMED", f"collection `{oname}` needs `discover_datasets` "
                                    f"(dynamic_tool.output_unclaimed)")
        elif not fwd and not disc:
            bad("OUTPUT-UNCLAIMED", f"output `{oname}` claims no file: give it `from_work_dir` or "
                                    f"`discover_datasets` (dynamic_tool.output_unclaimed). A UDT "
                                    f"does not hand your command an output path to write to.")
        if fwd:
            fwd_targets.append((oname, str(fwd)))
        for entry in disc or []:
            if isinstance(entry, dict) and entry.get("discover_via", "pattern") == "pattern":
                missing = sorted(DISCOVER_KEYS - set(entry))
                if missing:
                    bad("DISCOVER-INCOMPLETE",
                        f"output `{oname}`: a pattern discovery entry has NO defaults in the "
                        f"model, so all of them are required; missing {missing}")
        if t == "data" and not o.get("format") and not o.get("format_source"):
            note("NO-FORMAT", f"output `{oname}` declares neither `format` nor `format_source`; "
                              f"Galaxy will guess the datatype")

    # -- 5. templating: what `$(...)` means -------------------------------------------------
    #    ⛔ THIS HALF IS WHY THE SCRIPT EXISTS. Three tools in this repository were created,
    #    accepted, invoked, and had never been able to run: they wrapped a SHELL command
    #    substitution around a Galaxy expression. Galaxy ate the outer one, the command line came
    #    out empty, and the job died with `exit=None` before the container started. Measured on
    #    usegalaxy.org (26.1) and a local 25.0. The schema cannot see any of this.
    lib_fns = set()
    for r in d.get("requirements") or []:
        for src in (r or {}).get("expression_lib") or []:
            lib_fns |= set(re.findall(r"function\s+([A-Za-z_$][\w$]*)", str(src)))
    referenced = set()
    for line, expr, start, end in expressions(cmd):
        at = f":{line}"
        if "$(" in expr:
            bad("NESTED-SUBSTITUTION",
                f"`$(` nested inside `$(`: $({expr[:60]}…). Galaxy interpolates the OUTER one as "
                f"an expression and the job fails with an empty command line before anything runs. "
                f"Pass the path and let the inlined script read the file, or escape a literal "
                f"shell substitution as \\$(...).", at)
            continue
        m = re.match(r"\s*inputs\.([A-Za-z_][A-Za-z0-9_]*)", expr)
        if not m:
            fn = re.match(r"\s*([A-Za-z_$][\w$]*)\s*\(", expr)
            if fn and fn.group(1) in lib_fns:
                continue
            bad("NOT-A-GALAXY-EXPR",
                f"`$({expr[:50]})` is not an `inputs.<name>` reference and no `javascript` "
                f"requirement declares it. Galaxy evaluates `$(...)` as ECMAScript -- a shell "
                f"command substitution written here does not survive. Escape it as \\$(...).", at)
            continue
        ref = m.group(1)
        referenced.add(ref)
        if ref not in top_level_inputs:
            bad("UNDECLARED-INPUT-REF",
                f"`$(inputs.{ref})` names no declared input (dynamic_tool.undeclared_input_ref); "
                f"this tool declares {sorted(top_level_inputs)}", at)
            continue
        node = declared[ref]
        t, rest = node.get("type"), expr[m.end():]
        # ⚠ THE WHOLE EXPRESSION, NOT JUST ITS FIRST WORD. `$(inputs.family_list ? "--family-list
        # '" + inputs.family_list.path + "'" : "")` is the documented way to make an optional data
        # input into a flag, and reading only the head of it reported four DATA-NO-PATH problems
        # against correct code. A bare `inputs.X` with NOTHING after it is the mistake; an
        # expression that goes on to do something is the author's business.
        tail = rest.strip()
        bare = not tail.startswith(".")
        if t == "data" and not node.get("multiple"):
            if bare and not tail:
                bad("DATA-NO-PATH", f"`$(inputs.{ref})` on a `data` input renders the whole File "
                                    f"OBJECT, not a path. Write `$(inputs.{ref}.path)`.", at)
            elif tail.startswith("."):
                # ⛔ AND THE FIELD HAS TO EXIST. `.startswith(".")` alone accepted
                # `$(inputs.input.pth)`, which evaluates to undefined and hands the command a path
                # that is not there -- the exact failure DATA-NO-PATH was written for, one typo to
                # the right of where it was looking.
                field = re.match(r"\.([A-Za-z_][A-Za-z0-9_]*)", tail)
                if field and field.group(1) not in FILE_FIELDS:
                    bad("NO-SUCH-FILE-FIELD",
                        f"`$(inputs.{ref}.{field.group(1)})`: a `data` input renders a File record, "
                        f"which has no `{field.group(1)}`. It evaluates to undefined and the "
                        f"command runs on a path that does not exist. Fields: "
                        f"{sorted(FILE_FIELDS)}.", at)
        elif t == "data" and node.get("multiple"):
            if tail.startswith(".path"):
                bad("MULTIPLE-NO-MAP", f"`{ref}` is `multiple: true`, so it is an ARRAY of File "
                                       f"objects and `.path` on the array is undefined. Write "
                                       f"`$(inputs.{ref}.map(d => d.path).join(' '))`.", at)
        elif t == "data_collection" and tail.startswith(".path"):
            # ⛔ A COLLECTION IS AN ARRAY TOO. `.path` on it is undefined exactly as it is for a
            # `multiple: true` data input, and the collection branch only ever looked at the
            # no-dot case, so this spelling passed.
            #
            # ⚠ UNLESS THE TOOL SAYS IT MEANS IT, which is the same escape the workflow parity
            # check offers: `collection_probe` renders this spelling ON PURPOSE, because measuring
            # what it yields is the entire tool -- and what it measured is why `sourmash_panel`
            # parses a JSON array instead. Refusing a tool for demonstrating the hazard it exists
            # to document is the guard doing harm, so naming the construct in `help` downgrades it.
            # ⚠ NOT `declared` -- that name is the input table this loop is reading, and shadowing
            # it turned the next iteration into `TypeError: 'bool' object is not subscriptable`.
            said_so = f"inputs.{ref}.path" in str((help_ or {}).get("content") or "")
            (note if said_so else bad)(
                "MULTIPLE-NO-MAP",
                f"`{ref}` is a collection -- an ARRAY of File objects -- so `.path` on it is "
                f"undefined. Write `$(inputs.{ref}.map(d => d.path).join(' '))`."
                + (" (named in `help`, so read as deliberate)" if said_so else ""), at)
        elif t == "data_collection" and not tail:
            # ⚠ DELIBERATE IN THIS REPOSITORY, SO A NOTE. `$(inputs.xs)` on a collection renders a
            # JSON array of File records, and `sourmash_panel` captures exactly that into a heredoc
            # on purpose -- a job container never sees the collection itself. It is still worth
            # saying out loud, because the same spelling on a SINGLE data input is a bug.
            note("COLLECTION-RENDERED", f"`$(inputs.{ref})` renders the collection as a JSON array "
                                        f"of File records. Intended? The idiom for paths is "
                                        f"`$(inputs.{ref}.map(d => d.path).join(' '))`.", at)
        elif t in SCALAR_TYPES and tail.startswith(".path"):
            bad("PATH-ON-SCALAR", f"`{ref}` is a `{t}` -- a value, not a File -- so `.path` is "
                                  f"undefined. Write `$(inputs.{ref})`.", at)
        if t == "text" and not tail:
            # ⛔ SCALARS ARE NOT shlex.quote()d. `$(...)` in a shell_command is substituted as RAW
            # TEXT (unlike a base_command argument), so a text value carrying a space, a quote or a
            # `$` breaks the command -- and injects into it when the value is not yours.
            before = cmd[start - 1] if start else ""
            after = cmd[end + 1] if end + 1 < len(cmd) else ""
            if not (before == after and before in "'\""):
                bad("UNQUOTED-SCALAR",
                    f"free-text `$(inputs.{ref})` is interpolated UNQUOTED. A space or quote in the "
                    f"value breaks the command line (and injects into it). Quote it "
                    f"('$(inputs.{ref})') or pass it through a configfile.", at)
    for cf in d.get("configfiles") or []:
        for line, expr, _s, _e in expressions(str((cf or {}).get("content") or "")):
            m = re.match(r"\s*inputs\.([A-Za-z_][A-Za-z0-9_]*)", expr)
            if m:
                referenced.add(m.group(1))
                if m.group(1) not in top_level_inputs:
                    bad("UNDECLARED-INPUT-REF",
                        f"configfile {cf.get('name')!r} line {line}: `$(inputs.{m.group(1)})` names "
                        f"no declared input -- a configfile counts for this exactly as the command "
                        f"does")
    for unused in sorted(top_level_inputs - referenced):
        note("UNUSED-INPUT", f"input `{unused}` is declared and never referenced; it shows on the "
                             f"form and reaches nothing")

    # -- 6. Cheetah habits and shell-variable spellings --------------------------------------
    for pat, msg in (
            (r"^\s*#(for|if|end|set|while)\b", "Cheetah directive -- UDTs interpolate sandboxed "
                                               "ECMAScript in `$(...)`, there are no #directives"),
            (r"\$\{on_string\}|\$\{tool\.name\}|\$__tool_directory__",
             "a Cheetah macro; not available in a UDT")):
        for line, text in statements(cmd):
            if re.search(pat, text):
                bad("CHEETAH", f"{text.strip()[:60]!r}: {msg}", f":{line}")
    for m in re.finditer(r"(?<!\\)\$\{([A-Za-z_][A-Za-z0-9_]*)\}", cmd):
        # ⛔ MEASURED, NOT INFERRED. `${NAME}` in a shell_command is FATAL on this Galaxy while the
        # bare `$NAME` form works -- the braced spelling collides with the templating pass and the
        # job never builds. `$GALAXY_SLOTS` is right; `${GALAXY_SLOTS}` is a job that does not run.
        bad("BRACED-VAR", f"`${{{m.group(1)}}}`: the braced form does not survive templating. "
                          f"Write `${m.group(1)}`.", f":{cmd[:m.start()].count(chr(10)) + 1}")
    for m in re.finditer(r"(?<![\\$\w])\$([a-z][a-z0-9_]*)", cmd):
        if m.group(1) in top_level_inputs:
            bad("CHEETAH", f"`${m.group(1)}` is Cheetah's spelling and expands as a SHELL variable "
                           f"here -- almost always empty. Write `$(inputs.{m.group(1)})`.",
                f":{cmd[:m.start()].count(chr(10)) + 1}")

    # -- 7. $GALAXY_SLOTS: request it, and actually pass it ----------------------------------
    cores = any((r or {}).get("type") == "resource" and
                ("cores_min" in (r or {}) or "cores_max" in (r or {}))
                for r in d.get("requirements") or [])
    if "$GALAXY_SLOTS" in cmd and not cores:
        note("SLOTS-UNREQUESTED", "the command uses $GALAXY_SLOTS but no `resource` requirement "
                                  "asks for cores, so the job gets whatever the destination "
                                  "defaults to -- usually 1")
    hard = [(line, m.group(0)) for line, text in statements(cmd)
            for m in re.finditer(THREAD_FLAGS + r"[= ]\s*(\d+)", text) if m.group(1) != "1"]
    if hard and "$GALAXY_SLOTS" in cmd:
        bad("SLOTS-IGNORED",
            f"{hard[0][1]!r} hardcodes a thread count while the command also mentions "
            f"$GALAXY_SLOTS. Recording the allocation is not using it -- the flag's VALUE must be "
            f"$GALAXY_SLOTS or the job runs at the literal regardless of what it was given.",
            f":{hard[0][0]}")

    # -- 8. silent success: exit 0, nothing produced -----------------------------------------
    lines = [(ln, s.strip()) for ln, s in statements(cmd)
             if s.strip() and not s.strip().startswith("#")]
    if not re.search(r"^\s*set\s+-\w*e", cmd, re.MULTILINE):
        CLOSERS = ("fi", "done", "esac", "}", "else", "elif", "then", "do", ")", ";;")
        OPENERS = ("&&", "||", "|", ";", "{", "(", "then", "do", "else", "in")
        for (ln_a, a), (_ln_b, b) in itertools.pairwise(lines):
            if a.endswith(OPENERS) or b.startswith(CLOSERS) or a.startswith(CLOSERS):
                continue
            # ⚠ NARROWED, BECAUSE THE FIRST VERSION FLAGGED 16 OF 19 TOOLS AND WOULD HAVE TAUGHT
            # ITS READERS TO SKIP IT. Every UDT here opens by writing an awk or Python program
            # with `cat > prog <<'EOF'`, and a heredoc write that fails takes the very next
            # statement down with it (the interpreter cannot find its script); the same goes for
            # `set`, `cd`, an assignment, and an `echo` that is only talking to the log. What is
            # left is a statement that can fail on its own and be stepped over.
            if re.match(r"(set|export|cd|[A-Za-z_]\w*=)\b", a) or "<<" in a:
                continue
            if re.match(r"(echo|printf)\b", a) and not any(f in a for _o, f in fwd_targets):
                continue
            # ⚠ `|| true` IS THE AUTHOR SAYING SO. The probe tools run a dozen diagnostics that are
            # EXPECTED to fail on some hosts and swallow their status deliberately; flagging those
            # would be flagging the intent.
            if a.endswith("|| true"):
                continue
            # ⚠ A STATEMENT THAT IS ENTIRELY A GALAXY EXPRESSION IS NOT A STATEMENT YET. The
            # optional-input idiom renders either a fragment that already ends in `&&` or the empty
            # string, and the shell sees neither an unchained command nor a blank one -- a newline
            # after `&&` is a continuation. Judging it as written text flagged correct code.
            if a.startswith("$(") and a.endswith(")"):
                continue
            # ⚠ A NOTE, NOT A PROBLEM -- SAY WHAT WAS MEASURED. Attacking the one instance of this
            # in the repository (`samtools faidx seq.fa` / `cut -f1,2 seq.fa.fai`) could not make
            # it fire: samtools writes its index at the end, so a failure leaves no `.fai` and
            # `cut` fails too. The hazard is real and the correctness is an accident of the
            # wrapped tool's write order, which is not a thing to depend on -- but calling it a
            # PROBLEM would be claiming a break nobody has produced.
            note("UNCHAINED-COMMAND",
                 f"{a[:44]!r} is followed by {b[:34]!r} with no `&&` and no `set -e`. The job's "
                 f"exit code is the LAST statement's, so if the first fails the second still runs "
                 f"and the tool can exit 0 having produced nothing -- or a truncated output "
                 f"nothing downstream can tell from a real one. Chain them.", f":{ln_a}")
            break
    #    ⚠ AN AWK VARIABLE IS THE FILENAME. `awk -v out_bed=anchor.bed12` means the literal target
    #    never appears beside the redirect that writes it, so both checks below have to resolve
    #    `-v NAME=VALUE` bindings before they can say anything true about a target.
    aliases = dict(re.findall(r"-v\s+([A-Za-z_]\w*)=([^\s'\"]+)", cmd))
    for oname, fwd in fwd_targets:
        spellings = [fwd, *[k for k, v in aliases.items() if v == fwd]]
            # ⚠ SUBSTRING, AND DELIBERATELY LOOSE. `FasTAN -oscan` writes `scan.1ano` without the
            # literal name ever appearing, so an exact match reports a phantom on a correct tool.
        # ⚠ A from_work_dir MAY NAME A SUBDIRECTORY, and the file is then usually written by an
        # inlined script that only knows the basename -- `--output-dir outdir` plus a script that
        # opens `triage.tsv`. Testing the full path alone reported four phantoms against a tool
        # whose outputs are all produced correctly.
        base = fwd.rsplit("/", 1)[-1]
        if not any(s in cmd for s in [*spellings, base]):
            # ⚠ THE STEM IS A HINT, NOT A MATCH, so it is a NOTE. `FasTAN -oscan` writes
            # `scan.1ano` without the name ever appearing, which is why a stem fallback exists at
            # all -- but the fallback is wide: renaming a target to `chrom.tsv` still matches the
            # stem `chrom` of the `chrom.sizes` sitting in the command, and the check went quiet on
            # an output nothing writes. Split the two answers instead of merging them.
            if fwd.split(".")[0] in cmd or base.split(".")[0] in cmd:
                note("PHANTOM-OUTPUT", f"output `{oname}`: from_work_dir {fwd!r} never appears in "
                                       f"shell_command; only its STEM does, which may be another "
                                       f"file's name. Confirm something writes exactly {base!r}.")
            else:
                bad("PHANTOM-OUTPUT", f"output `{oname}`: from_work_dir {fwd!r} -- neither it, nor "
                                      f"{base!r}, nor their stem appears anywhere in shell_command")
    created = {}
    for oname, fwd in fwd_targets:
        spellings = "|".join([re.escape(fwd)]
                             + [re.escape(k) for k, v in aliases.items() if v == fwd])
        pat = (r"(?:touch|:\s*>|printf\s+(?:\"\"|'')\s*>|echo\s+-n\s+>)\s*[\"']?"
               rf"(?:{spellings})\b")
        created[oname] = bool(re.search(pat, cmd))
    server, ran = authoritative(path, d)
    problems += server
    lint.ran_authoritative = ran

    if any(created.values()) and not all(created.values()):
        # ⛔ MEASURED ON A REAL INPUT. `anchor_prep` opens its BED with `printf "" > out_bed` and a
        # comment explaining why -- awk opens a redirect on FIRST WRITE, so zero rows means no file
        # -- and then writes its isoforms table through `| "sort > out_iso"` with no such guard. A
        # GFF3 whose features are typed `transcript` instead of `mRNA` (AUGUSTUS, StringTie, EVM,
        # and gffread's own output) therefore produced a correct BED12, NO isoforms file at all,
        # and exit 0. An author who pre-creates one output and not the other has not decided the
        # second one is safe; they have forgotten it.
        bad("OUTPUT-CREATE-ASYMMETRY",
            f"{sorted(k for k, v in created.items() if v)} are pre-created before the command can "
            f"fail to write them, and {sorted(k for k, v in created.items() if not v)} are not. "
            f"An unwritten from_work_dir file is a green job with a missing output.")
    return problems, notes


# ---------------------------------------------------------------- self-test --------------------
#: (label, find, replace, expected code). Every case is injected into a REAL tool in this repo, so
#: a check that stops matching the way these files are actually written fails here rather than
#: going quietly green forever.
DEFECTS = [
    ("udt/samtools_faidx.gxtool.yml", "command field is the XML one",
     "shell_command: |", "command: |", "XMLISM"),
    ("udt/samtools_faidx.gxtool.yml", "id with uppercase",
     "id: brc-samtools-faidx", "id: BRC-Samtools", "BAD-ID"),
    ("udt/samtools_faidx.gxtool.yml", "container as a dict",
     "container: quay.io/biocontainers/samtools:1.24--h9dcdb79_1",
     "container:\n  type: docker\n  image: samtools", "BLANK-CONTAINER"),
    ("udt/samtools_faidx.gxtool.yml", "a container tag that does not exist",
     "container: quay.io/biocontainers/samtools:1.24--h9dcdb79_1",
     "container: quay.io/biocontainers/samtools:1.24--hdeadbeef_9", "NO-SUCH-CONTAINER"),
    ("udt/samtools_faidx.gxtool.yml", "container pinned to latest",
     "container: quay.io/biocontainers/samtools:1.24--h9dcdb79_1",
     "container: quay.io/biocontainers/samtools:latest", "CONTAINER-UNPINNED"),
    ("udt/samtools_faidx.gxtool.yml", "output claims no file",
     "    from_work_dir: chrom.sizes\n", "", "OUTPUT-UNCLAIMED"),
    ("udt/samtools_faidx.gxtool.yml", "from_work_dir names a file nothing writes",
     "from_work_dir: chrom.sizes", "from_work_dir: nosuchthing.tsv", "PHANTOM-OUTPUT"),
    ("udt/samtools_faidx.gxtool.yml", "reference to an undeclared input",
     "$(inputs.input.path)", "$(inputs.inputt.path)", "UNDECLARED-INPUT-REF"),
    ("udt/samtools_faidx.gxtool.yml", "data input used without .path",
     "'$(inputs.input.path)'", "'$(inputs.input)'", "DATA-NO-PATH"),
    ("udt/samtools_faidx.gxtool.yml", "shell substitution wrapped round an expression",
     "cut -f1,2 seq.fa.fai > chrom.sizes",
     "cut -f1,2 $(basename '$(inputs.input.path)') > chrom.sizes", "NESTED-SUBSTITUTION"),
    ("udt/samtools_faidx.gxtool.yml", "unescaped literal shell substitution",
     "cut -f1,2 seq.fa.fai > chrom.sizes",
     "echo $(date) && cut -f1,2 seq.fa.fai > chrom.sizes", "NOT-A-GALAXY-EXPR"),
    ("udt/samtools_faidx.gxtool.yml", "braced shell variable",
     "cut -f1,2 seq.fa.fai", "cut -f1,2 ${PWD}/seq.fa.fai", "BRACED-VAR"),
    ("udt/samtools_faidx.gxtool.yml", "help sent as a bare string",
     "help:\n  format: markdown\n  content: |", "help: |", "HELP-SHAPE"),
    ("udt/samtools_faidx.gxtool.yml", "XML boolean fields on an input",
     "  - name: input\n    type: data",
     "  - name: input\n    type: data\n    truevalue: --x", "XMLISM"),
    ("udt/samtools_faidx.gxtool.yml", "an XML-only input type",
     "  - name: input\n    type: data", "  - name: input\n    type: data_column", "BAD-INPUT-TYPE"),
    # ⚠ THE INJECTION ADDS A PRE-CREATE; IT DOES NOT REMOVE ONE. Written the other way round --
    # deleting anchor_prep's `printf "" > out_bed` -- the fixture left NO output pre-created, which
    # is not an asymmetry and correctly did not fire. (anchor_prep itself trips this check for
    # real; that is a finding, not a fixture.)
    ("udt/samtools_faidx.gxtool.yml", "one output pre-created, the other not",
     "  samtools faidx seq.fa", "  touch seq.fa.fai\n  samtools faidx seq.fa",
     "OUTPUT-CREATE-ASYMMETRY"),
    ("udt/phase_c2_triage.gxtool.yml", "a regex validator anchored at one end only",
     'expression: "[A-Za-z0-9_.-]+$"', 'expression: "[A-Za-z0-9_.-]+"', "UNANCHORED-VALIDATOR"),
    ("udt/sourmash_panel.gxtool.yml", "a scalar addressed as a File",
     "--ksize '$(inputs.ksize)'", "--ksize '$(inputs.ksize.path)'", "PATH-ON-SCALAR"),
]


def self_test() -> int:
    """Prove each check goes RED on a real file, and that the repo's own tools stay clean."""
    import tempfile
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td) / "case.gxtool.yml"
        for rel, label, old, new, want in DEFECTS:
            src = (ROOT / rel).read_text()
            if src.count(old) < 1:
                print(f"  BROKEN-FIXTURE  {label}: {old[:40]!r} no longer occurs in {rel}")
                failed += 1
                continue
            tmp.write_text(src.replace(old, new, 1))
            codes = {c for c, _w, _m in lint(tmp)[0]}
            hit = want in codes
            print(f"  {'pass' if hit else 'FAIL'}  {label:<44} expect {want}"
                  f"{'' if hit else '  got ' + str(sorted(codes))}")
            failed += not hit
    print(f"\n  {len(DEFECTS) - failed}/{len(DEFECTS)} injected defects were detected")
    return 1 if failed else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tools", type=pathlib.Path, nargs="*",
                    help="UDT files to check (default: every udt/*.gxtool.yml)")
    ap.add_argument("--self-test", action="store_true",
                    help="inject each defect class into this repo's own tools and confirm it is "
                         "caught. A checker nobody has seen go red is not evidence of anything.")
    ap.add_argument("--notes", action="store_true", help="also print convention notes")
    ap.add_argument("--offline", action="store_true",
                    help="skip the container-existence probe -- the one check here that needs a "
                         "network, and the one failure nothing offline can rule out")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    paths = args.tools or sorted((ROOT / "udt").glob("*.gxtool.yml"))
    problems, notes = [], []
    for p in paths:
        pr, nt = lint(p, containers=not args.offline)
        problems += pr
        notes += nt
    print(f"  {len(paths)} UDT definition(s) checked against the udt-authoring skill")
    if args.offline:
        print("  ⚠ CONTAINER PROBE SKIPPED (--offline) — a UDT naming an image that does not "
              "exist passes every other check here and dies at job start. NOT a pass.")
    if getattr(lint, "ran_authoritative", False):
        print("  + galaxy-tool-util is installed, so the SERVER's own schema validation and lint "
              "ran too\n")
    else:
        # ⚠ A SKIP IS NOT A PASS AND PRINTS AS ITS OWN LINE. Everything above is a
        # re-implementation; the authority for the schema half is the server's model, and this run
        # did not have it. `pip install galaxy-tool-util` and the line changes.
        print("  ⚠ galaxy-tool-util NOT installed — the SCHEMA half here is a re-implementation "
              "and the server's own validator did not run. `pip install galaxy-tool-util`.\n")
    if problems:
        print(f"  ⛔ {len(problems)} PROBLEM(S)")
        for k, w, m in problems:
            print(f"    {k:<24} {w:<38} {m}")
    else:
        print("  ✅ no schema, templating or silent-success problem found")
    if notes:
        print(f"\n  {len(notes)} note(s){'' if args.notes else ' (--notes to show)'}")
        if args.notes:
            for k, w, m in notes:
                print(f"    {k:<24} {w:<38} {m}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
