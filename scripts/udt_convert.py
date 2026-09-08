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
import pathlib
import re
import sys
import urllib.request

from galaxy_tool_source.binding import load_tool
from galaxy_tool_source.cheetah_refs import tool_cheetah_references
from galaxy_tool_source.macros import expanded_detection_root

CONVERTER_VERSION = "0.1.0"
DEPOT = "https://depot.galaxyproject.org/singularity/"

#: References a converted command may still contain. Everything else is a refusal, because the
#: point of the reference model is to enumerate what is there rather than to hope.
PORTABLE_REFS = {"tool.name", "on_string"}


class Refusal(Exception):
    """Raised with the reason a wrapper cannot be converted. The reason is the product."""


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def requirements(doc) -> list[tuple[str, str]]:
    """`[(package, version)]` with macros expanded and @TOKENS@ substituted."""
    root = expanded_detection_root(doc)
    return [((e.text or "").strip(), e.get("version") or "")
            for e in root.iter() if e.tag == "requirement" and e.get("type") == "package"]


#: The depot listing cannot plausibly be smaller than this. It held 124,750 images on 2026-09-08,
#: so the floor is two orders of magnitude below the truth and exists only to tell a real listing
#: from a broken one.
MIN_DEPOT_IMAGES = 1000


def depot_images(cache: pathlib.Path | None) -> list[str]:
    """Every image published in the depot listing, as `pkg:version[--build]`.

    ⛔ PARSED FROM THE `href`, NOT FROM THE LINK TEXT. The listing is an nginx autoindex, and its
    link text is the DECODED name (`10x_bamtofastq:1.4.1`) while only the href carries `%3A`. A
    regex looking for `%3A` in the text therefore matched 0 of 124,750 entries when this was
    measured (2026-09-08), which does not fail: `resolve_container` simply finds nothing and every
    wrapper is refused with "no published biocontainer for <pkg>=<version>" -- a WRONG REFUSAL,
    the exact failure the header says made xml_to_udt.py untrustworthy. The href is also the more
    complete of the two: nginx truncates long link text with `..`, so the text yields 121,219
    names against the href's 124,750.
    """
    if cache and cache.exists():
        names = [n for n in cache.read_text().splitlines() if n.strip()]
        # ⛔ A TRUNCATED CACHE IS INDISTINGUISHABLE FROM AN EMPTY DEPOT, AND BOTH REFUSE. An
        # interrupted write leaves a 0-byte or half-written file that this used to trust
        # unconditionally, turning every later conversion into "no published biocontainer" with no
        # way to tell that wrong refusal from a real one. Cheap to check, so check it.
        if len(names) < MIN_DEPOT_IMAGES or not any(":" in n for n in names):
            raise Refusal(f"{cache} holds {len(names)} entries, which is not a depot listing "
                          f"(expected >{MIN_DEPOT_IMAGES:,} of the form pkg:version). It is "
                          f"truncated or was written by an interrupted run -- delete it and let "
                          f"this refetch, rather than reading refusals off a broken cache.")
        return names
    with urllib.request.urlopen(DEPOT, timeout=120) as fh:
        body = fh.read().decode("utf-8", "replace")
    names = [m.group(1).replace("%3A", ":")
             for m in re.finditer(r'href="([A-Za-z0-9_.\-]+%3A[^"]+)"', body)]
    if len(names) < MIN_DEPOT_IMAGES:
        raise Refusal(f"the depot listing parsed to {len(names)} images, which cannot be right. "
                      f"Its markup has changed and this regex no longer matches -- fix the parse "
                      f"rather than letting every conversion refuse for a fabricated reason.")
    if cache:
        cache.write_text("\n".join(names))
    return names


def build_number(name: str) -> int:
    """The `_N` build counter at the end of a biocontainer tag, or -1 if it has none."""
    m = re.search(r"_(\d+)$", name)
    return int(m.group(1)) if m else -1


def resolve_container(pkg: str, version: str, images: list[str]) -> str:
    """The biocontainer for one package, looked up in the depot listing rather than assembled.

    ⚠ The build suffix (`--h0b57e2e_0`) is not derivable from the requirement, which is exactly why
    this looks it up instead of formatting a string.

    ⛔ "NEWEST BUILD WINS" MEANS THE BUILD NUMBER, AND `sorted()[-1]` DOES NOT READ IT. The tag is
    `--<hash>_<build>`, so a lexical maximum sorts on the HASH first and on the build as a string
    second: `_2` beats `_10`, and a higher build under an alphabetically earlier hash loses
    outright. Measured over the whole depot listing (124,750 images, 2026-09-08), lexical selection
    picks a lower build number for 3,414 of 74,237 package:version groups -- 4.6%. A rebuild
    published precisely to fix a broken image is exactly the case that gets ignored, and the stale
    choice is then frozen into a provenance stamp as if it had been verified.

    ⚠ TIES ARE BROKEN LEXICALLY AND NOT REFUSED, DELIBERATELY. 4,910 of those groups (6.6%) publish
    more than one hash at the highest build number, so refusing the ambiguity would refuse one
    conversion in fifteen for something that is normal. The tie-break is the tag, which makes the
    choice arbitrary but reproducible.

    ⚠ AND WHAT IS VERIFIED IS PRESENCE IN THE DEPOT, WHICH IS THE SINGULARITY MIRROR, while the
    string emitted and stamped is a `quay.io/biocontainers` reference. They are the same
    BioContainers build under two distributions and agree in practice, but this is a mirror lookup
    and not a query against quay -- so a name present in one and absent from the other yields a
    stamped image that fails at run time with `manifest unknown`, which no linter can see.
    """
    # ⚠ TWO TAG SHAPES, AND MISSING THE SECOND IS A FALSE REFUSAL. Most biocontainers carry a build
    # suffix (`ucsc-chainstitchid:482--h0b57e2e_0`), but some are published under the bare
    # version -- `python:3.12` is, and it is the image four of this repo's UDTs already use. A
    # lookup that only matched `pkg:version--*` reported "no published biocontainer for
    # python=3.12" for tools whose container is sitting in the depot, which is the same class of
    # wrong answer that made xml_to_udt.py's refusals untrustworthy.
    suffixed = [n for n in images if n.startswith(f"{pkg}:{version}--")]
    if suffixed:
        return f"quay.io/biocontainers/{max(suffixed, key=lambda n: (build_number(n), n))}"
    if f"{pkg}:{version}" in images:
        return f"quay.io/biocontainers/{pkg}:{version}"
    raise Refusal(f"no published biocontainer for {pkg}={version} (looked for {pkg}:{version} and "
                  f"{pkg}:{version}--* among {len(images):,} depot images)")


def command_text(doc) -> str:
    node = doc.root.find("command")
    if node is None:
        raise Refusal("the wrapper has no <command>")
    return "".join(node.itertext())


def check_translatable(doc, cmd: str) -> list[str]:
    """Refuse every shape this converter cannot carry. Returns the data-reference names it may."""
    if re.search(r"^\s*#(if|else|elif|end|for|set|silent|import|def)\b", cmd, re.MULTILINE):
        d = re.search(r"^\s*#(\w+)", cmd, re.MULTILINE).group(1)
        raise Refusal(f"the command uses the Cheetah directive `#{d}`; in a shell_command a leading "
                      f"`#` is a COMMENT, so translating it away would silently drop the logic")
    names = []
    for ref in tool_cheetah_references(doc.root):
        raw = ref.name.strip("${}")
        if raw in PORTABLE_REFS:
            continue
        if raw.startswith("__"):
            raise Refusal(f"the command reads `{ref.name}`; a container has no directory beside the "
                          f"wrapper, so the helper it names would not exist at run time")
        if "element_identifier" in raw:
            raise Refusal(f"the command reads `{ref.name}`; a job receives PATHS and identifiers do "
                          f"not reach it by any route -- restructure the tool to run per element")
        names.append(raw.split(".")[0])
    return names


def param_name(p) -> str:
    """A param's Galaxy name, which is derived from `argument=` when `name=` is absent.

    ⛔ `p.get("name")` IS `None` FOR `<param argument="--x" type="data">`, AND THE CONSEQUENCE WAS
    A REFUSAL THAT BLAMED THE WRONG THING. Galaxy derives the name from the argument by stripping
    the dashes and mapping the rest to underscores, so `--out-idx` is `$out_idx`. With the name
    read as None, `translate()` matched nothing and the wrapper was refused with "these references
    survived translation and would reach the shell as literals: ['$idx', ...]" -- which names the
    symptom and points nowhere near `argument=`. tools/odgi/paths.xml is exactly this shape.
    """
    name = p.get("name")
    if name:
        return name
    arg = (p.get("argument") or "").lstrip("-").replace("-", "_")
    if not arg:
        raise Refusal("a <param> has neither `name` nor `argument`, so it has no Galaxy name to "
                      "translate; fix the wrapper")
    return arg


def params(doc) -> list[dict]:
    out = []
    section = doc.root.find("inputs")
    if section is None:
        raise Refusal("the wrapper has no <inputs>")
    # ⚠ SCOPED TO <inputs>, because <tests> carries <param> elements too -- and a test param has no
    # `type`, so an unscoped scan refuses every wrapper that HAS tests, for a reason that is false.
    for p in section.iter("param"):
        name = param_name(p)
        if p.get("type") != "data":
            raise Refusal(f"parameter `{name}` is type={p.get('type')!r}; this converter "
                          f"carries data inputs only, so port it by hand")
        out.append({"name": name, "format": p.get("format", "data"),
                    "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
    return out


def outputs(doc) -> list[dict]:
    out = []
    section = doc.root.find("outputs")
    if section is None:
        raise Refusal("the wrapper has no <outputs>")
    # ⛔ A <collection> OUTPUT WAS DROPPED WITHOUT A WORD, WHICH IS THE ONE THING THIS CONVERTER
    # PROMISES NOT TO DO. Scanning only for `data` meant a wrapper whose sole output is a
    # collection generated a UDT with an EMPTY `outputs:` list -- a tool that runs and produces
    # nothing. Eight wrappers here are that shape (masking_table, multiz_fold,
    # phase_e_graph_edges among them). The UDT equivalent is `discover_datasets`, which this does
    # not attempt; refusing by name is at least honest about that.
    for c in section.iter("collection"):
        raise Refusal(f"output `{c.get('name')}` is a <collection>, and this converter emits only "
                      f"`data` outputs -- a UDT needs a discovery pattern (`type: collection` plus "
                      f"`discover_datasets`) that has to be written by hand. Left implicit, this "
                      f"output would vanish and the tool would produce nothing.")
    for d in section.iter("data"):
        out.append({"name": d.get("name"), "format": d.get("format", "data"),
                    "label": d.get("label", "")})
    if not out:
        raise Refusal("the wrapper declares no <data> output, so the generated tool would claim "
                      "nothing it writes")
    return out


def translate(cmd: str, ins: list[dict], outs: list[dict]) -> tuple[str, dict[str, str]]:
    """Rewrite a `&&`-joined command into a shell_command, and name each output's work-dir file.

    ⚠ ONLY `&&`-JOINED LINES. Galaxy runs `<command>` through a shell that stops at the first
    failure because the wrapper chains with `&&`; a `shell_command` is an ordinary script, so the
    chain becomes `set -e` plus one statement per line. Any other joiner (`;`, `|` at a line end,
    a bare newline between statements) changes what a failure does, so it is refused.
    """
    body = cmd.strip()
    if re.search(r";\s*$", body, re.MULTILINE):
        raise Refusal("the command joins statements with `;`, which does not stop at a failure the "
                      "way the `&&` chain does; port it by hand")
    stmts = [s.strip() for s in re.split(r"&&\s*\n?", body) if s.strip()]
    workfiles = {o["name"]: f"{o['name']}.dat" for o in outs}
    # ⛔ LONGEST NAME FIRST, BECAUSE A PREFIX COLLISION CORRUPTS THE LONGER REFERENCE AND THE
    # LEFTOVER GUARD THEN CANNOT SEE IT. With outputs `output` and `output_gz`, substituting in
    # declaration order rewrote `$output_gz` as `output.dat_gz`: the `$` is consumed, so the guard
    # below finds NOTHING LEFT and the conversion succeeds. The generated tool writes
    # `output_gz.dat`, which its `from_work_dir` claims and nothing produces, while a stray
    # `output.dat_gz` is written and claimed by no output -- a green conversion, a green
    # registration, and one empty dataset at run time. Two wrappers here collide this way
    # (fasta_concat and pansn_rename: `output`/`output_gz`) and pggb has `output_lay`/
    # `output_layout_png`; all three are currently shielded only by unrelated refusals.
    subs = sorted(([(f"${i['name']}", f"$(inputs.{i['name']}.path)") for i in ins]
                   + [(f"${o['name']}", workfiles[o["name"]]) for o in outs]),
                  key=lambda kv: -len(kv[0]))
    lines = []
    for s in stmts:
        for ref, repl in subs:
            s = s.replace(f"'{ref}'", f"'{repl}'" if repl.startswith("$(") else repl)
            s = s.replace(ref, repl)
        lines.append(s)
    joined = "\n".join(lines)
    # ⚠ THE BRACED FORM IS MATCHED WHOLE, up to its closing brace, so a default expansion like
    # `${GALAXY_SLOTS:-1}` is reported as itself. Stopping at the first non-word character printed
    # `${GALAXY_SLOTS` in the refusal -- an unbalanced fragment that does not appear in the file the
    # reader is about to open.
    left = re.findall(r"(\\?)(\$\{[^}]*\}|\$[A-Za-z_][\w.]*)", joined)
    # ⚠ AN ESCAPED SHELL VARIABLE IS NOT AN UNTRANSLATED GALAXY REFERENCE, AND SAYING SO IS THE
    # WHOLE POINT. `\${GALAXY_SLOTS:-1}` is Cheetah-escaped -- it is a SHELL variable, idiomatic in
    # every threaded wrapper here (tools/odgi/paths.xml carries one) -- and reporting it as a
    # reference that "would reach the shell as literals" states the opposite of the truth. The
    # refusal is still right, for a different reason: `${NAME}` in a `shell_command` is fatal,
    # measured, while bare `$NAME` works and GALAXY_SLOTS is exported. So it is named separately,
    # with the fix attached.
    braced_shell = sorted({m[1] for m in left if m[0] == "\\" and m[1].startswith("${")})
    if braced_shell:
        raise Refusal(f"the command uses the braced shell form {braced_shell}, which is FATAL in a "
                      f"UDT shell_command (measured). Rewrite it unbraced -- `$GALAXY_SLOTS` works "
                      f"and is exported -- then convert.")
    galaxy_refs = sorted({m[1] for m in left
                          if m[0] != "\\" and m[1].strip("${}") not in PORTABLE_REFS})
    if galaxy_refs:
        raise Refusal(f"these references survived translation and would reach the shell as literals: "
                      f"{galaxy_refs}")
    # ⛔ THE CHEETAH ESCAPE MUST COME OFF, OR THE VARIABLE ARRIVES AS ITS OWN NAME. `\$GALAXY_SLOTS`
    # in a classic <command> is Cheetah being told to leave the dollar alone, so the shell receives
    # `$GALAXY_SLOTS` and expands it. A shell_command is not templated, so the backslash survives
    # into bash -- where `\$` is a LITERAL dollar (measured: `echo \$FOO` prints `$FOO`) -- and the
    # tool is handed the eight characters of the variable's name instead of the thread count. That
    # is not a crash; it is a tool silently running with a garbage argument.
    return joined.replace("\\$", "$"), workfiles


#: What is left of a Galaxy label idiom once its `${...}` parts are removed. `${tool.name} on
#: ${on_string}` -- the single most common output label in Galaxy -- leaves the bare word "on",
#: which is truthy, so the `or fallback` below could never fire for it and five wrappers here
#: (fasta_concat, pansn_rename and three under vg/) would have shipped an output labelled `on`.
LABEL_REMNANTS = {"on", "of", "and", "in", "for", "with", "from"}


def clean_label(label: str, fallback: str) -> str:
    """A UDT label is literal text, so the wrapper's Cheetah template must not survive into it.

    `${tool.name} on ${on_string}: stitched chains` is Galaxy templating that a UDT does not
    evaluate -- it would ship to the tool form verbatim -- and the unquoted `: ` would also make
    the YAML a mapping. Keep the human half after the last colon, and fall back when the template
    was the whole label.
    """
    text = re.sub(r"\$\{[^}]*\}", "", label)
    text = text.split(":")[-1].strip(" :")
    text = re.sub(r"\s+", " ", text).strip()
    if not text or text.lower() in LABEL_REMNANTS:
        return fallback
    return text


def repo_relative(xml: pathlib.Path) -> str:
    """The wrapper's path relative to the repository root.

    ⛔ SPLITTING ON A CHECKOUT DIRECTORY NAME IS NOT A PATH CALCULATION. This read
    `str(xml).split("brc-tools-fork/")[-1]`, so from a clone named anything else -- and with an
    absolute argument, which is the normal way to invoke a script -- the whole absolute path went
    into the stamp. check_udt_provenance.py then joins it onto ROOT, which for an absolute path
    yields that same absolute path: the check silently validates a file OUTSIDE the repository on
    the machine that generated it, and reports "source ... no longer exists" everywhere else.
    """
    root = pathlib.Path(__file__).resolve().parents[1]
    try:
        return str(xml.resolve().relative_to(root))
    except ValueError:
        raise Refusal(f"{xml} is outside this repository ({root}), so there is no path that a "
                      f"stamp could record for it") from None


def yaml_block(text: str, indent: str) -> str:
    return "\n".join(indent + ln if ln.strip() else "" for ln in text.splitlines())


def convert(xml: pathlib.Path, images: list[str]) -> tuple[str, str]:
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
        raise Refusal("the wrapper declares no conda requirement, so there is no image to resolve")
    pkg, version = reqs[0]
    container = resolve_container(pkg, version, images)

    cmd = command_text(doc)
    check_translatable(doc, cmd)
    ins, outs = params(doc), outputs(doc)
    shell, workfiles = translate(cmd, ins, outs)

    tool_id = doc.root.get("id")
    udt_id = "brc-" + re.sub(r"(?<!^)(?=[A-Z])", "-", tool_id).replace("_", "-").lower()
    desc = (doc.root.findtext("description") or "").strip()
    helpnode = doc.root.find("help")
    helptext = ("".join(helpnode.itertext()).strip() if helpnode is not None else "")

    macros = xml.parent / "macros.xml"
    prov = {
        "source": repo_relative(xml),
        "tool_id": tool_id,
        # ⛔ THE WHOLE WRAPPER IS HASHED, BECAUSE THE COMMAND IS NOT THE ONLY THING COPIED FROM IT.
        # This document also carries the wrapper's <description>, every input's format/label/help,
        # every output's FORMAT and label, and the entire <help> body -- none of which the two
        # narrow hashes below cover. Changing `<data name="output" format="chain">` to
        # `format="tabular"` left check_udt_provenance.py reporting "0 stale" while the generated
        # UDT still declared `format: chain`, demonstrated on this very wrapper. The narrow hashes
        # are kept because they LOCALISE a change once the file hash has detected it.
        "tool_sha256": _sha(xml.read_text()),
        "command_sha256": _sha(cmd),
        "requirements": f"{pkg}={version}",
        "macros_sha256": _sha(macros.read_text()) if macros.exists() else "(no macros.xml)",
        "container": container,
        "converter": f"scripts/udt_convert.py {CONVERTER_VERSION}",
    }
    head = ["# ⛔ GENERATED by scripts/udt_convert.py (fork-side) -- do not hand-edit.",
            "# Edit the wrapper under tools/ and regenerate. The stamp below is what",
            "# scripts/check_udt_provenance.py verifies, using the standard library only.",
            "#",
            "# provenance:"]
    head += [f"#   {k}: {v}" for k, v in prov.items()]

    doc_lines = [
        *head,
        "class: GalaxyUserTool",
        f"id: {udt_id}",
        'version: "0.1.0"',
        f"name: {tool_id} (BRC UDT)",
        f"description: {desc}" if desc else "",
        f"container: {container}",
        "shell_command: |",
        yaml_block("set -e\n" + shell, "  "),
        "inputs:",
    ]
    for i in ins:
        doc_lines += [f"  - name: {i['name']}", "    type: data", f"    format: {i['format']}"]
        if i["label"]:
            doc_lines.append(f"    label: {i['label']}")
        if i["help"]:
            doc_lines.append(f"    help: {i['help']}")
    doc_lines.append("outputs:")
    for o in outs:
        doc_lines += [f"  - name: {o['name']}", "    type: data", f"    format: {o['format']}",
                      f"    from_work_dir: {workfiles[o['name']]}"]
        if o["label"]:
            doc_lines.append(f"    label: {clean_label(o['label'], o['name'])}")
    if helptext:
        doc_lines += ["help:", "  format: markdown", "  content: |",
                      yaml_block(helptext, "    ")]
    return udt_id, "\n".join(ln for ln in doc_lines if ln is not None) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xml", nargs="+", type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, help="write <id>.gxtool.yml here (default: stdout)")
    ap.add_argument("--depot-cache", type=pathlib.Path,
                    help="file holding a cached depot listing (written if absent)")
    a = ap.parse_args()

    images = depot_images(a.depot_cache)
    rc = 0
    for xml in a.xml:
        try:
            udt_id, text = convert(xml, images)
        except Refusal as e:
            print(f"REFUSING {xml}: {e}", file=sys.stderr)
            rc = 1
            continue
        if a.out:
            # named for the UDT id, not the directory, so `udt/` reads as one set:
            # brc-chain-stitch-id -> chain_stitch_id.gxtool.yml
            path = a.out / (udt_id.removeprefix("brc-").replace("-", "_") + ".gxtool.yml")
            path.write_text(text)
            print(f"wrote {path} ({udt_id})")
        else:
            print(text)
    return rc


if __name__ == "__main__":
    sys.exit(main())
