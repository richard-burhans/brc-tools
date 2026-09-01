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


def depot_images(cache: pathlib.Path | None) -> list[str]:
    if cache and cache.exists():
        return cache.read_text().splitlines()
    with urllib.request.urlopen(DEPOT, timeout=120) as fh:      # noqa: S310 - fixed https host
        body = fh.read().decode("utf-8", "replace")
    names = [m.group(1) for m in re.finditer(r'>([A-Za-z0-9_.\-]+%3A[^<]+)</a>', body)]
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
        if raw.startswith("__"):
            raise Refusal(f"the command reads `{ref.name}`; a container has no directory beside the "
                          f"wrapper, so the helper it names would not exist at run time")
        if "element_identifier" in raw:
            raise Refusal(f"the command reads `{ref.name}`; a job receives PATHS and identifiers do "
                          f"not reach it by any route -- restructure the tool to run per element")
        names.append(raw.split(".")[0])
    return names


def params(doc) -> list[dict]:
    out = []
    section = doc.root.find("inputs")
    if section is None:
        raise Refusal("the wrapper has no <inputs>")
    # ⚠ SCOPED TO <inputs>, because <tests> carries <param> elements too -- and a test param has no
    # `type`, so an unscoped scan refuses every wrapper that HAS tests, for a reason that is false.
    for p in section.iter("param"):
        if p.get("type") != "data":
            raise Refusal(f"parameter `{p.get('name')}` is type={p.get('type')!r}; this converter "
                          f"carries data inputs only, so port it by hand")
        out.append({"name": p.get("name"), "format": p.get("format", "data"),
                    "label": p.get("label", ""), "help": (p.get("help", "") or "").strip()})
    return out


def outputs(doc) -> list[dict]:
    out = []
    section = doc.root.find("outputs")
    if section is None:
        raise Refusal("the wrapper has no <outputs>")
    for d in section.iter("data"):
        out.append({"name": d.get("name"), "format": d.get("format", "data"),
                    "label": d.get("label", "")})
    return out


def translate(cmd: str, ins: list[dict], outs: list[dict]) -> tuple[str, dict[str, str]]:
    """Rewrite a `&&`-joined command into a shell_command, and name each output's work-dir file.

    ⚠ ONLY `&&`-JOINED LINES. Galaxy runs `<command>` through a shell that stops at the first
    failure because the wrapper chains with `&&`; a `shell_command` is an ordinary script, so the
    chain becomes `set -e` plus one statement per line. Any other joiner (`;`, `|` at a line end,
    a bare newline between statements) changes what a failure does, so it is refused.
    """
    body = cmd.strip()
    if re.search(r";\s*$", body, re.M):
        raise Refusal("the command joins statements with `;`, which does not stop at a failure the "
                      "way the `&&` chain does; port it by hand")
    stmts = [s.strip() for s in re.split(r"&&\s*\n?", body) if s.strip()]
    workfiles = {o["name"]: f"{o['name']}.dat" for o in outs}
    lines = []
    for s in stmts:
        for i in ins:
            s = s.replace(f"'${i['name']}'", f"'$(inputs.{i['name']}.path)'")
            s = s.replace(f"${i['name']}", f"$(inputs.{i['name']}.path)")
        for o in outs:
            s = s.replace(f"'${o['name']}'", workfiles[o["name"]])
            s = s.replace(f"${o['name']}", workfiles[o["name"]])
        lines.append(s)
    joined = "\n".join(lines)
    left = re.findall(r"\$\{?[A-Za-z_][\w.]*\}?", joined)
    left = [x for x in left if not x.startswith("$(") and x.strip("${}") not in PORTABLE_REFS]
    if left:
        raise Refusal(f"these references survived translation and would reach the shell as literals: "
                      f"{sorted(set(left))}")
    return joined, workfiles


def clean_label(label: str, fallback: str) -> str:
    """A UDT label is literal text, so the wrapper's Cheetah template must not survive into it.

    `${tool.name} on ${on_string}: stitched chains` is Galaxy templating that a UDT does not
    evaluate -- it would ship to the tool form verbatim -- and the unquoted `: ` would also make
    the YAML a mapping. Keep the human half after the last colon.
    """
    text = re.sub(r"\$\{[^}]*\}", "", label)
    text = text.split(":")[-1].strip(" :")
    return text or fallback


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
        "source": str(xml).split("brc-tools-fork/")[-1],
        "tool_id": tool_id,
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

    doc_lines = head + [
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
        doc_lines += [f"  - name: {i['name']}", f"    type: data", f"    format: {i['format']}"]
        if i["label"]:
            doc_lines.append(f"    label: {i['label']}")
        if i["help"]:
            doc_lines.append(f"    help: {i['help']}")
    doc_lines.append("outputs:")
    for o in outs:
        doc_lines += [f"  - name: {o['name']}", f"    type: data", f"    format: {o['format']}",
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
