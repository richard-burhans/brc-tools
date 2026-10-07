#!/usr/bin/env python3
"""Convert a classic Galaxy tool XML in this repo into a User-Defined Tool (UDT) document.

⛔ THIS SCRIPT IS FORK-SIDE AND IS NOT PART OF ANY UPSTREAM PULL REQUEST. It needs
`galaxy-tool-source` and `galaxy-tool-util`, and brc-tools has no dependency manifest at all --
its one generator CI job installs `pyyaml` inline. What upstream receives is the GENERATED
document plus a provenance stamp that `scripts/check_udt_provenance.py` verifies with the standard
library alone. Upstream never has to install or trust this.

Run it without installing anything:

    uv run --with galaxy-tool-source==0.3.7 --with galaxy-tool-util==26.1.1 \
        python scripts/udt_convert.py tools/chainStitchId/chainStitchId.xml --out udt/

⚠ WHY A LIBRARY AND NOT MORE REGEXES. `scripts/xml_to_udt.py` reads the raw XML with ElementTree
and finds `<requirements>` in 14 of the 50 wrappers here. It is not that the other 36 lack them --
they declare them through `<expand macro="requirements"/>`, which that reader cannot follow. Using
`galaxy_tool_source.macros.expanded_detection_root` the count is 48 of 50 with zero errors, and the
`@TOOL_VERSION@` tokens come back already substituted. The same applies to the command: a regex for
surviving `$name` is an approximation of a lexer, and `galaxy_tool_source.cheetah_refs` is the
lexer.

⛔ IT STILL REFUSES RATHER THAN GUESSES, and the refusals are now ACCURATE rather than incidental.
A tool that runs and is wrong costs more than one that does not exist -- so the shapes below are
rejected by name, with the measurement or the reason attached:

  * more than one conda package -- a UDT gets exactly ONE container. The mulled image name IS
    derivable (`galaxy.tool_util.deps.mulled.util.v2_image_name`, one call, contradicting the note
    generated into `udt/lc_classify.gxtool.yml`); what is missing is a PUBLISHED image. This prints
    the computed name so it can be registered with BioContainers (`planemo container_register`).
  * any Cheetah directive in the command -- `#if`, `#for`, `#set`. In a `shell_command` a leading
    `#` is a comment, so a survivor is silently dropped logic.
  * `$__tool_directory__` -- the helper script beside the wrapper does not exist in a container.
  * `element_identifier` on a whole collection -- a job receives paths; identifiers do not reach it
    by any route.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import sys
import urllib.request

from galaxy_tool_source.binding import load_tool
from galaxy_tool_source.cheetah_refs import tool_cheetah_references
from galaxy_tool_source.macros import expanded_detection_root

CONVERTER_VERSION = "0.3.0"
DEPOT = "https://depot.galaxyproject.org/singularity/"

#: References a converted command may still contain. Everything else is a refusal, because the
#: point of the reference model is to enumerate what is there rather than to hope.
PORTABLE_REFS = {"tool.name", "on_string", "GALAXY_SLOTS"}

#: ⛔ `${NAME}` IS FATAL IN A UDT shell_command AND `$NAME` WORKS -- measured, not inferred. So a
#: wrapper's idiomatic `\${GALAXY_SLOTS:-N}` cannot travel as written: the braces have to go, and
#: the `:-N` fallback with them. ⚠ That is only safe because GALAXY_SLOTS IS exported to a UDT job
#: (measured =1), so dropping the default cannot leave an empty argument behind.
SLOTS_FORMS = re.compile(r"\$\{GALAXY_SLOTS(?::-\s*\d+)?\}")


class Refusal(Exception):
    """Raised with the reason a wrapper cannot be converted. The reason is the product."""


def _repo_relative(path: pathlib.Path) -> str:
    """`path` relative to the git repository holding it, so a stamp never records a home directory."""
    p = path.resolve()
    for parent in p.parents:
        if (parent / ".git").exists():
            return str(p.relative_to(parent))
    return p.name


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def requirements(doc) -> list[tuple[str, str]]:
    """`[(package, version)]` with macros expanded and @TOKENS@ substituted."""
    root = expanded_detection_root(doc)
    return [((e.text or "").strip(), e.get("version") or "")
            for e in root.iter() if e.tag == "requirement" and e.get("type") == "package"]


def depot_images(cache: pathlib.Path | None) -> list[str]:
    if cache and cache.exists():
        return cache.read_text().splitlines()
    with urllib.request.urlopen(DEPOT, timeout=120) as fh:      # noqa: S310 - fixed https host
        body = fh.read().decode("utf-8", "replace")
    # ⚠ BOTH SPELLINGS. The depot index writes older images as `pkg%3Aver` and newer ones as
    # `pkg:ver` (measured 2026-09-27: ucsc-axtsort 332-469 encoded, 482 plain). Matching only `%3A`
    # made every recent build read as "no published biocontainer".
    names = [m.group(1) for m in re.finditer(r'>([A-Za-z0-9_.\-]+(?:%3A|:)[^<]+)</a>', body)]
    names = [n.replace("%3A", ":") for n in names]
    if cache:
        cache.write_text("\n".join(names))
    return names


def resolve_container(pkg: str, version: str, images: list[str]) -> str:
    """The biocontainer for one package, VERIFIED against the depot rather than assembled.

    ⚠ The build suffix (`--h0b57e2e_0`) is not derivable from the requirement, which is exactly why
    this looks it up instead of formatting a string. Newest build wins; ties are broken by the
    listing's own order, which is lexical and therefore stable across runs.
    """
    # ⚠ TWO TAG SHAPES, AND MISSING THE SECOND IS A FALSE REFUSAL. Most biocontainers carry a build
    # suffix (`ucsc-chainstitchid:482--h0b57e2e_0`), but some are published under the bare
    # version -- `python:3.12` is, and it is the image four of this repo's UDTs already use. A
    # lookup that only matched `pkg:version--*` reported "no published biocontainer for
    # python=3.12" for tools whose container is sitting in the depot, which is the same class of
    # wrong answer that made xml_to_udt.py's refusals untrustworthy.
    suffixed = sorted(n for n in images if n.startswith(f"{pkg}:{version}--"))
    if suffixed:
        return f"quay.io/biocontainers/{suffixed[-1]}"
    if f"{pkg}:{version}" in images:
        return f"quay.io/biocontainers/{pkg}:{version}"
    raise Refusal(f"no published biocontainer for {pkg}={version} (looked for {pkg}:{version} and "
                  f"{pkg}:{version}--* in the depot)")


def command_text(doc) -> str:
    node = doc.root.find("command")
    if node is None:
        raise Refusal("the wrapper has no <command>")
    return "".join(node.itertext())


def check_translatable(doc, cmd: str) -> list[str]:
    """Refuse every shape this converter cannot carry. Returns the data-reference names it may."""
    if re.search(r"^\s*#(if|else|elif|end|for|set|silent|import|def)\b", cmd, re.M):
        d = re.search(r"^\s*#(\w+)", cmd, re.M).group(1)
        raise Refusal(f"the command uses the Cheetah directive `#{d}`; in a shell_command a leading "
                      f"`#` is a COMMENT, so translating it away would silently drop the logic")
    names = []
    for ref in tool_cheetah_references(doc.root):
        raw = ref.name.strip("${}")
        if raw in PORTABLE_REFS:
            continue
        if raw == "__tool_directory__":
            # ▶ Allowed here and resolved later by `tool_dir_helpers`, which refuses any use that
            # is not `$__tool_directory__/<file>` with that file present.
            continue
        if raw.startswith("__"):
            raise Refusal(f"the command reads `{ref.name}`; a container has no directory beside the "
                          f"wrapper, so the helper it names would not exist at run time")
        if "element_identifier" in raw:
            raise Refusal(f"the command reads `{ref.name}`; a job receives PATHS and identifiers do "
                          f"not reach it by any route -- restructure the tool to run per element")
        names.append(raw.split(".")[0])
    return names


#: `$__tool_directory__/<file>` -- the only form of that reference a UDT can carry.
TOOL_DIR_HELPER = re.compile(r"\$__tool_directory__/([\w.\-]+)")


def helper_delim(fname: str) -> str:
    """The heredoc delimiter for an inlined helper. One definition, used by the check and the
    emitter, so they can never disagree about what would terminate the body early."""
    return "TOOLDIR_" + re.sub(r"[^A-Z0-9]", "_", fname.upper())


def tool_dir_helpers(cmd: str, xml: pathlib.Path) -> list[tuple[str, str]]:
    """`[(filename, text)]` for each helper the command runs from beside the wrapper.

    ▶ A UDT HAS NO DIRECTORY, SO THE HELPER TRAVELS INSIDE IT. The wrapper's idiom is
    `python '$__tool_directory__/x.py'`; a container has no x.py beside it, which is why this was
    a flat refusal. But the converter already carries a <configfile> as a quoted heredoc, and a
    helper script is the same problem: text that must exist as a file before the command runs. So
    it is read from beside the XML and emitted the same way. The repo's hand-written UDTs for
    these tools do exactly this by hand.

    ⛔ ONLY THE `/<file>` FORM. A bare `$__tool_directory__`, or one joined to a glob or a
    subdirectory, is still refused: there is nothing to enumerate, and inlining "whatever is in
    that directory" would quietly change what the tool runs.
    """
    out = []
    for fname in dict.fromkeys(TOOL_DIR_HELPER.findall(cmd)):
        helper = xml.parent / fname
        if not helper.is_file():
            raise Refusal(f"the command runs `$__tool_directory__/{fname}` and no such file sits "
                          f"beside the wrapper, so there is nothing to carry into the UDT")
        text = helper.read_text()
        # ⛔ The heredoc is quoted, so the body passes through literally -- but a line equal to the
        # delimiter would end it early and silently truncate the script.
        delim = helper_delim(fname)
        if any(ln.strip() == delim for ln in text.splitlines()):
            raise Refusal(f"{fname} contains a line equal to its heredoc delimiter {delim}")
        out.append((fname, text))
    return out


def params(doc) -> list[dict]:
    out = []
    section = doc.root.find("inputs")
    if section is None:
        raise Refusal("the wrapper has no <inputs>")
    # ⚠ SCOPED TO <inputs>, because <tests> carries <param> elements too -- and a test param has no
    # `type`, so an unscoped scan refuses every wrapper that HAS tests, for a reason that is false.
    for p in section.iter("param"):
        name = p.get("name") or (p.get("argument") or "").lstrip("-").replace("-", "_")
        kind = p.get("type")
        if kind in ("integer", "float"):
            # ▶ A scalar with a default carries straight across: a UDT input of the same type, read
            # as `$(inputs.<name>)`, its default as `value:` (Galaxy's YamlIntegerParameter; `default:` is
            # rejected as an extra field). Only the two numeric types -- a select or boolean has semantics
            # (options, truevalue/falsevalue) this converter would have to guess at.
            value = p.get("value")
            if value is None:
                raise Refusal(f"parameter `{name}` is {kind} with no default; a UDT input needs one")
            out.append({"name": name, "type": kind, "default": value,
                        "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
            continue
        if kind == "text":
            # ▶ A text input carries across as a UDT text, read as `$(inputs.<name>)`. Its validators do
            # not: the wrapper's script must check the value itself, as it would for any argv.
            out.append({"name": name, "type": "text", "default": p.get("value", ""),
                        "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
            continue
        if kind == "select":
            # ▶ A single select carries across as a UDT select, read as `$(inputs.<name>)`: its value
            # is the option's value, as Cheetah renders it. A multiple select would be a list.
            if p.get("multiple") == "true":
                raise Refusal(f"parameter `{name}` is a multiple select; port it by hand")
            opts = [(o.get("value", ""), "".join(o.itertext()).strip()) for o in p.iter("option")]
            if not opts:
                raise Refusal(f"parameter `{name}` is a select with no static options")
            default = next((v for v, _ in ((o.get("value", ""), o) for o in p.iter("option"))
                            if _.get("selected") == "true"), opts[0][0])
            out.append({"name": name, "type": "select", "options": opts, "default": default,
                        "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
            continue
        if kind == "boolean":
            # ▶ A UDT CARRIES A BOOLEAN, BUT NOT GALAXY'S SEMANTICS FOR ONE. Galaxy renders `$flag`
            # as its truevalue/falsevalue; a UDT's `$(inputs.flag)` is evaluated as ECMAScript and
            # yields the bare lowercase `true`/`false` (measured on usegalaxy.org 26.1 -- see
            # fasta_uppercase.gxtool.yml, the first UDT in the repo to carry one). So the flag has
            # to travel in a TERNARY, which `translate` builds; this only records the two values.
            tv, fv = p.get("truevalue", ""), p.get("falsevalue", "")
            # ⛔ THE TERNARY IS BUILT BY STRING CONCATENATION, so a value holding a quote, a
            # backslash or a `$` would either break the expression or smuggle a second
            # substitution into it. Refuse rather than escape: the escaping rules differ between
            # Galaxy's evaluator and the shell, and getting them subtly wrong is silent.
            for val, which in ((tv, "truevalue"), (fv, "falsevalue")):
                if any(c in val for c in "'\"$`"):
                    raise Refusal(f"parameter `{name}` has {which}={val!r}, which holds a quote, "
                                  f"backslash, `$` or backtick; the UDT ternary is built by "
                                  f"concatenation, so port this one by hand")
            out.append({"name": name, "type": "boolean", "truevalue": tv, "falsevalue": fv,
                        "default": p.get("checked", "false") == "true",
                        "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
            continue
        if kind != "data":
            raise Refusal(f"parameter `{name}` is type={kind!r}; this converter carries data, "
                          f"text, select, integer and float inputs only, so port it by hand")
        out.append({"name": name, "type": "data", "format": p.get("format", "data"),
                    "multiple": p.get("multiple") == "true",
                    "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
    return out


def configfiles(doc) -> list[tuple[str, str]]:
    """`[(name, text)]` for each <configfile>, which a UDT carries as a heredoc.

    ⛔ WHY THIS EXISTS. Galaxy replaces every newline in <command> with a space before running it
    (lib/galaxy/tools/evaluation.py, "Remove newlines from command line"), so a wrapper cannot keep
    an inline script in its command; the idiom is a <configfile>. A UDT's shell_command keeps its
    newlines, so the same text travels as `cat > <name>.cfg <<'DELIM'`.
    ⛔ ONLY A CONFIGFILE WITH NO CHEETAH IN IT. A configfile is itself a Cheetah template; a `$name`
    or `#if` inside one would need the same translation as the command, and a heredoc quoted with
    'DELIM' would pass it through as a literal. So any `$` or directive is a refusal -- the wrapper
    should hand its values to the script as argv instead.
    """
    node = doc.root.find("configfiles")
    if node is None:
        return []
    out = []
    for cf in node.iter("configfile"):
        name, text = cf.get("name"), "".join(cf.itertext())
        # ▶ `\$` and `\#` are Cheetah's escapes for a literal `$` and `#` -- how a POSIX sh script,
        # which cannot avoid `$`, lives in a configfile. They render as `$` and `#`, so they carry
        # across; an UNESCAPED `$` is a Cheetah reference and still a refusal.
        if re.search(r"(?<!\\)\$", text):
            raise Refusal(f"configfile `{name}` contains an unescaped `$` (a Cheetah reference); pass "
                          f"values to it as argv, and write the shell's own `$` as `\\$`")
        if re.search(r"^\s*#(if|else|elif|end|for|set|silent|import|def)\b", text, re.M):
            raise Refusal(f"configfile `{name}` uses a Cheetah directive")
        if re.search(r"^\s*##", text, re.M):
            raise Refusal(f"configfile `{name}` has a `##` line, which Cheetah drops as a comment; "
                          f"escape it as `\\#\\#` if it is meant to reach the script")
        text = text.replace("\\$", "$").replace("\\#", "#")
        # ⛔ `$(` IS GALAXY'S UDT EXPRESSION SYNTAX. A shell command substitution written `$( ... )` in
        # the script is evaluated as JavaScript when the job is built, and the job fails with no command
        # line and no stderr (measured on laila 26.1, 2026-09-28: gt-multiz-roast). Backticks are safe.
        if "$(" in text:
            raise Refusal(f"configfile `{name}` uses `$( ... )`, which a UDT evaluates as a JavaScript "
                          f"expression; write the shell command substitution with backticks")
        delim = f"CONFIGFILE_{name.upper()}"
        if re.search(rf"^{delim}$", text, re.M):
            raise Refusal(f"configfile `{name}` contains its own heredoc delimiter {delim}")
        out.append((name, text))
    return out


def outputs(doc) -> list[dict]:
    out = []
    section = doc.root.find("outputs")
    if section is None:
        raise Refusal("the wrapper has no <outputs>")
    for d in section:
        # ⛔ A COMMENT IS A CHILD TOO. ElementTree yields comment nodes from an iteration over
        # children, and their `.tag` is a FUNCTION, not a string -- so a wrapper that documents
        # its outputs with an XML comment was refused for having an "output element
        # <cyfunction Comment ...>", which names nothing a reader can act on.
        if not isinstance(d.tag, str):
            continue
        if d.tag == "data":
            out.append({"name": d.get("name"), "format": d.get("format", "data"),
                        "label": d.get("label", "")})
        elif d.tag == "collection":
            # ▶ A list discovered by one filename pattern carries across as a UDT collection output
            # (measured working on .org 26.1, 2026-09-08). Anything richer is refused, never dropped:
            # skipping an output the converter cannot read would register a tool with one fewer output.
            disc = list(d.iter("discover_datasets"))
            if d.get("type") != "list" or len(disc) != 1 or not disc[0].get("pattern"):
                raise Refusal(f"output collection `{d.get('name')}` is not a list with one "
                              f"discover_datasets pattern; port it by hand")
            dd = disc[0]
            out.append({"name": d.get("name"), "collection": True, "label": d.get("label", ""),
                        "pattern": dd.get("pattern"), "directory": dd.get("directory"),
                        "format": dd.get("format") or dd.get("ext") or "data",
                        "sort_by": dd.get("sort_by", "filename")})
        else:
            raise Refusal(f"output element <{d.tag}> is not carried; port it by hand")
    return out


def translate(cmd: str, ins: list[dict], outs: list[dict],
              cfgs: list[tuple[str, str]] = (),
              helpers: list[tuple[str, str]] = ()) -> tuple[str, dict[str, str]]:
    """Rewrite a `&&`-joined command into a shell_command, and name each output's work-dir file.

    ⚠ ONLY `&&`-JOINED LINES. Galaxy runs `<command>` through a shell that stops at the first
    failure because the wrapper chains with `&&`; a `shell_command` is an ordinary script, so the
    chain becomes `set -e` plus one statement per line. Any other joiner (`;`, `|` at a line end,
    a bare newline between statements) changes what a failure does, so it is refused.
    """
    # Galaxy drops Cheetah `##` comment lines before running a command; so does this, or they reach
    # the UDT's script as shell comments that describe a template nobody renders.
    body = "\n".join(ln for ln in cmd.strip().splitlines() if not ln.lstrip().startswith("##")).strip()
    if re.search(r";\s*$", body, re.M):
        raise Refusal("the command joins statements with `;`, which does not stop at a failure the "
                      "way the `&&` chain does; port it by hand")
    stmts = [s.strip() for s in re.split(r"&&\s*\n?", body) if s.strip()]
    workfiles = {o["name"]: f"{o['name']}.dat" for o in outs if not o.get("collection")}
    lines = []
    for s in stmts:
        for name, _text in cfgs:
            s = s.replace(f"'${name}'", f"{name}.cfg").replace(f"${name}", f"{name}.cfg")
        for i in ins:
            if i["type"] == "data" and i.get("multiple"):
                # ▶ Cheetah renders a multiple data input as its paths joined by commas; a UDT gets a
                # list, and only a JavaScript expression walks it. Measured on laila 26.1 2026-09-28:
                # `.map(...)` works, an indexed `[0].path` does not, and element_identifier is undefined.
                # ⛔ TWO SHAPES, AND ONE OF THEM RENDERED AS '' WITH NO ERROR. A direct run passes a list of
                # files; a workflow mapping over a nested list passes a list holding ONE object keyed by
                # element identifier ({"tupChi1": {File}, ...}), so `.map(e => e.path)` gave [undefined] and
                # roast ran on nothing (laila, growler_multiz, 2026-09-29). Both are flattened here.
                ref = (f"$(inputs.{i['name']}.reduce(function(a, x){{ return a.concat(x && x.path ? [x] : "
                       f"Object.keys(x).map(function(k){{ return x[k]; }})); }}, [])"
                       f".map(function(e){{ return e.path; }}).join(','))")
            elif i["type"] == "data":
                ref = f"$(inputs.{i['name']}.path)"
            elif i["type"] == "boolean":
                # ⛔ A QUOTED BOOLEAN WOULD PASS AN EMPTY ARGUMENT, NOT NO ARGUMENT. `'$flag'`
                # becomes `''` when the flag is off, and a program receiving one empty argv entry
                # is not the same as receiving none -- several of these tools would read it as an
                # empty positional. Galaxy's own renderer has the same shape, so this is a
                # pre-existing hazard in the wrapper; refuse rather than carry it across silently.
                if f"'${i['name']}'" in s:
                    raise Refusal(f"the command quotes `${i['name']}`, a boolean; when it is off "
                                  f"that renders as an empty argument rather than none. Unquote it "
                                  f"in the wrapper (a flag needs no quoting) and convert again")
                # ⚠ DOUBLE quotes inside the ternary, so it survives being embedded anywhere the
                # command already uses single quotes.
                ref = (f'$(inputs.{i["name"]} ? "{i["truevalue"]}" : "{i["falsevalue"]}")')
            else:
                ref = f"$(inputs.{i['name']})"
            s = s.replace(f"'${i['name']}'", f"'{ref}'")
            s = s.replace(f"${i['name']}", ref)
        for o in outs:
            if o.get("collection"):
                continue
            s = s.replace(f"'${o['name']}'", workfiles[o["name"]])
            s = s.replace(f"${o['name']}", workfiles[o["name"]])
        lines.append(s)
    joined = "\n".join(lines)
    # ▶ The braced GALAXY_SLOTS forms collapse to the bare one BEFORE the surviving-reference
    # check, which is what lets a wrapper keep writing the Galaxy idiom while the UDT gets the
    # only form that works there.
    joined = SLOTS_FORMS.sub("$GALAXY_SLOTS", joined)
    # ▶ The helper now sits in the working directory, so its reference becomes a plain relative
    # path. Done after the `&&` split, so the rewrite never reaches into a heredoc body.
    joined = TOOL_DIR_HELPER.sub(r"\1", joined)
    # The heredocs go FIRST and after translation, so neither the `&&` split nor the reference
    # rewrite ever touches the script text.
    heredocs = [f"cat > {name}.cfg <<'CONFIGFILE_{name.upper()}'\n{text.strip(chr(10))}\n"
                f"CONFIGFILE_{name.upper()}" for name, text in cfgs]
    # ⚠ Helpers go in the SAME place and for the same reason as a configfile: a file the command
    # needs to exist before it runs. Emitted BEFORE them, so a helper can never depend on one.
    heredocs = [f"cat > {fname} <<'{helper_delim(fname)}'\n{text.rstrip(chr(10))}\n"
                f"{helper_delim(fname)}" for fname, text in helpers] + heredocs
    left = re.findall(r"\$\{?[A-Za-z_][\w.]*\}?", joined)
    left = [x for x in left if not x.startswith("$(") and x.strip("${}") not in PORTABLE_REFS]
    if left:
        raise Refusal(f"these references survived translation and would reach the shell as literals: "
                      f"{sorted(set(left))}")
    return "\n".join([*heredocs, joined]), workfiles


def version_of(doc) -> str:
    """The wrapper's own `version`, @TOKENS@ expanded -- e.g. `1.04.52+galaxy1`.

    ⚠ A UDT cannot be updated, so re-registering a CHANGED tool at the SAME version leaves two
    active registrations that nothing can tell apart. Carrying the wrapper's version makes a
    changed wrapper a new version, and it is PEP 440, which usegalaxy.org requires of a UDT.
    """
    root = expanded_detection_root(doc)
    return root.get("version") or "0.1.0"


def clean_label(label: str, fallback: str) -> str:
    """A UDT label is literal text, so the wrapper's Cheetah template must not survive into it.

    `${tool.name} on ${on_string}: stitched chains` is Galaxy templating that a UDT does not
    evaluate -- it would ship to the tool form verbatim -- and the unquoted `: ` would also make
    the YAML a mapping. Keep the human half after the last colon.
    """
    text = re.sub(r"\$\{[^}]*\}", "", label)
    text = text.split(":")[-1].strip(" :")
    return fallback if text.lower() in ("", "on", "of", "from") else text


def yaml_block(text: str, indent: str) -> str:
    return "\n".join(indent + ln if ln.strip() else "" for ln in text.splitlines())


def convert(xml: pathlib.Path, images: list[str], id_prefix: str = "brc-",
            wrapper_version: bool = False, base_image: str = "") -> tuple[str, str]:
    doc = load_tool(str(xml))
    reqs = requirements(doc)
    if len(reqs) != 1:
        if len(reqs) > 1:
            from galaxy.tool_util.deps.mulled.util import build_target, v2_image_name
            name = v2_image_name([build_target(p, v) for p, v in reqs])
            published = [n for n in images if n.startswith(name.split(":")[0] + ":")]
            raise Refusal(
                f"the wrapper needs {len(reqs)} packages ({', '.join(p for p, _ in reqs)}) and a UDT "
                f"gets ONE container. The mulled name IS derivable: {name.split(':')[0]} -- "
                + ("it is published, pin it by hand" if published else
                   "it is NOT published; register it (planemo container_register) or split the tool"))
        if not reqs and base_image:
            # ▶ A WRAPPER WITH NO REQUIREMENT STILL NEEDS AN IMAGE, AND ONLY THE OPERATOR KNOWS
            # WHICH. These are the pure-script tools: they declare no package because they need
            # none, so there is nothing to derive a container from. The repo's hand-written UDTs
            # for exactly these tools all pin `quay.io/biocontainers/python:3.12`, but that is a
            # CHOICE about which interpreter the script is run under, not a fact about the
            # wrapper. ⛔ So it is named on the command line and never defaulted: guessing an
            # interpreter version for someone else's script is how a tool runs under 3.12 and
            # breaks under 3.13 with nothing recording the decision.
            container = base_image
            pkg, version = "", ""
        else:
            raise Refusal("the wrapper declares no conda requirement, so there is no image to "
                          "resolve; pass --base-image to name one (the repo's hand-written UDTs "
                          "for script-only tools use quay.io/biocontainers/python:3.12)")
    if reqs:
        pkg, version = reqs[0]
        container = resolve_container(pkg, version, images)

    cmd = command_text(doc)
    check_translatable(doc, cmd)
    ins, outs, cfgs = params(doc), outputs(doc), configfiles(doc)
    shell, workfiles = translate(cmd, ins, outs, cfgs, tool_dir_helpers(cmd, xml))

    tool_id = doc.root.get("id")
    udt_id = id_prefix + re.sub(r"(?<!^)(?=[A-Z])", "-", tool_id).replace("_", "-").lower()
    desc = (doc.root.findtext("description") or "").strip()
    helpnode = doc.root.find("help")
    helptext = ("".join(helpnode.itertext()).strip() if helpnode is not None else "")

    macros = xml.parent / "macros.xml"
    prov = {
        "source": _repo_relative(xml),
        "tool_id": tool_id,
        "command_sha256": _sha(cmd),
        **({f"configfile_{n}_sha256": _sha(t) for n, t in cfgs}),
        "requirements": f"{pkg}={version}",
        "macros_sha256": _sha(macros.read_text()) if macros.exists() else "(no macros.xml)",
        "container": container,
        "converter": f"scripts/udt_convert.py {CONVERTER_VERSION}",
    }
    # ⚠ A UDT OUTPUT HAS ONE FORMAT. A wrapper's <change_format> cannot be carried, and dropping it
    # silently would mislabel some runs' output, so it is recorded here and printed.
    outs_node = doc.root.find("outputs")
    changed = ([d.get("name") for d in outs_node if d.find("change_format") is not None]
               if outs_node is not None else [])
    if changed:
        prov["change_format_dropped"] = ", ".join(changed)
        print(f"WARNING {xml}: <change_format> on {changed} is not carried; the UDT's output always has "
              f"the declared format", file=sys.stderr)
    head = ["# ⛔ GENERATED by scripts/udt_convert.py (fork-side) -- do not hand-edit.",
            "# Edit the wrapper under tools/ and regenerate. The stamp below is what",
            "# scripts/check_udt_provenance.py verifies, using the standard library only.",
            "#",
            "# provenance:"]
    head += [f"#   {k}: {v}" for k, v in prov.items()]

    doc_lines = head + [
        "class: GalaxyUserTool",
        f"id: {udt_id}",
        f'version: "{version_of(doc) if wrapper_version else "0.1.0"}"',
        f"name: {tool_id} ({'BRC ' if id_prefix == 'brc-' else ''}UDT)",
        f"description: {json.dumps(desc)}" if desc else "",
        f"container: {container}",
        "shell_command: |",
        yaml_block("set -e\n" + shell, "  "),
        "inputs:",
    ]
    for i in ins:
        if i["type"] == "data":
            doc_lines += [f"  - name: {i['name']}", "    type: data", f"    format: {i['format']}"]
            if i.get("multiple"):
                doc_lines.append("    multiple: true")
        elif i["type"] == "text":
            doc_lines += [f"  - name: {i['name']}", "    type: text", f"    value: {json.dumps(i['default'])}"]
        elif i["type"] == "select":
            # ⚠ The default is the option marked `selected: true`; a select has no `value:` key and
            # /api/unprivileged_tools refuses one as extra_forbidden (laila 26.1, 2026-09-28).
            doc_lines += [f"  - name: {i['name']}", "    type: select", "    options:"]
            for value, label in i["options"]:
                doc_lines += [f"      - label: {json.dumps(label)}", f"        value: {json.dumps(value)}"]
                if value == i["default"]:
                    doc_lines.append("        selected: true")
        elif i["type"] == "boolean":
            # ⚠ BEFORE the zero-default branch below: `float(False)` is 0.0, so an unchecked
            # boolean would otherwise be declared `optional: true` with no value and render as
            # `null`, and `null ? a : b` takes the FALSE arm by accident rather than by design.
            doc_lines += [f"  - name: {i['name']}", "    type: boolean",
                          f"    value: {'true' if i['default'] else 'false'}"]
        elif float(i["default"]) == 0:
            # ⛔ A ZERO DEFAULT FAILS GALAXY 26.1's UDT LINT ("TestsCaseValidation: Serious problem
            # parsing tool source"), for 0, "0" and 0.0 alike, measured on laila 2026-09-27; 1 passes.
            # So a zero default is declared OPTIONAL with no value, and an unset optional renders as
            # the literal `null` -- the wrapper's command must read `null` as 0 (measured the same day).
            doc_lines += [f"  - name: {i['name']}", f"    type: {i['type']}", "    optional: true"]
        else:
            doc_lines += [f"  - name: {i['name']}", f"    type: {i['type']}", f"    value: {i['default']}"]
        if i["label"]:
            doc_lines.append(f"    label: {json.dumps(i['label'])}")
        if i["help"]:
            doc_lines.append(f"    help: {json.dumps(i['help'])}")
    doc_lines.append("outputs:")
    for o in outs:
        if o.get("collection"):
            doc_lines += [f"  - name: {o['name']}", "    type: collection", "    collection_type: list",
                          "    discover_datasets:", "      - discover_via: pattern",
                          f"        pattern: {json.dumps(o['pattern'])}",
                          f"        format: {o['format']}", "        visible: false",
                          f"        sort_key: {json.dumps(o['sort_by'])}"]
            if o["directory"]:
                doc_lines.append(f"        directory: {json.dumps(o['directory'])}")
            if o["label"]:
                doc_lines.append(f"    label: {json.dumps(clean_label(o['label'], o['name']))}")
            continue
        doc_lines += [f"  - name: {o['name']}", f"    type: data", f"    format: {o['format']}",
                      f"    from_work_dir: {workfiles[o['name']]}"]
        if o["label"]:
            doc_lines.append(f"    label: {json.dumps(clean_label(o['label'], o['name']))}")
    if helptext:
        doc_lines += ["help:", "  format: markdown", "  content: |",
                      yaml_block(helptext, "    ")]
    return udt_id, "\n".join(ln for ln in doc_lines if ln is not None) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xml", nargs="+", type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, help="write <id>.gxtool.yml here (default: stdout)")
    ap.add_argument("--id-prefix", default="brc-",
                    help="UDT id prefix; `brc-` for this repo's wrappers, another for a sibling repo's")
    ap.add_argument("--wrapper-version", action="store_true",
                    help="version the UDT as the wrapper is versioned, instead of 0.1.0")
    ap.add_argument("--base-image", default="",
                    help="container for a wrapper that declares NO conda requirement (the "
                         "script-only tools). Never defaulted: which interpreter the script runs "
                         "under is a decision, not a fact about the wrapper. The repo's "
                         "hand-written UDTs for these use quay.io/biocontainers/python:3.12")
    ap.add_argument("--depot-cache", type=pathlib.Path,
                    help="file holding a cached depot listing (written if absent)")
    a = ap.parse_args()

    images = depot_images(a.depot_cache)
    rc = 0
    for xml in a.xml:
        try:
            udt_id, text = convert(xml, images, a.id_prefix, a.wrapper_version, a.base_image)
        except Refusal as e:
            print(f"REFUSING {xml}: {e}", file=sys.stderr)
            rc = 1
            continue
        if a.out:
            # named for the UDT id, not the directory, so `udt/` reads as one set:
            # brc-chain-stitch-id -> chain_stitch_id.gxtool.yml
            path = a.out / (udt_id.removeprefix(a.id_prefix).replace("-", "_") + ".gxtool.yml")
            path.write_text(text)
            print(f"wrote {path} ({udt_id})")
        else:
            print(text)
    return rc


if __name__ == "__main__":
    sys.exit(main())
