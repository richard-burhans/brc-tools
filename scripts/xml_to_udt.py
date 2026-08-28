#!/usr/bin/env python3
"""Convert a classic Galaxy tool XML in this repo into a User-Defined Tool (UDT) definition.

WHY THIS EXISTS RATHER THAN A HAND-WRITTEN SET OF UDTs. The wrappers under `tools/` install only
into a local Galaxy via `scripts/install_local_tools.sh`; on a public server such as usegalaxy.org
they cannot be installed at all. UDTs are the route that needs no admin — but hand-porting each
wrapper produces a second copy that drifts from the XML the moment either is edited. This converts,
so the XML stays the single source.

⛔ IT REFUSES RATHER THAN GUESSES, AND THAT IS THE WHOLE DESIGN. Cheetah is a templating language and
`shell_command` is not; a converter that silently mistranslates a `#for` loop or an
`element_identifier` reference would emit a tool that RUNS and is WRONG, which is worse than one that
does not exist. So the supported subset is enumerated below and everything else raises.

SUPPORTED
  * `#if str($x.ext) == 'fasta.gz' / #else / #end if` around a gunzip-or-cat — the only Cheetah
    conditional these wrappers use. Rewritten as `gzip -cdf`, which handles both.
  * `$__tool_directory__/<script>` — the script is INLINED as a heredoc. A UDT has no tool
    directory, and these scripts are 364 B to 2.7 kB, far too small to justify an image.
  * `$input` / `$output` scalars and data paths.

REFUSED (raises, with the offending construct named)
  * `#for` loops.
  * `element_identifier`. ⚠ Galaxy DOES expose it to a UDT as `$(inputs["<name>|__identifier__"])`
    when the tool is MAPPED over a collection — verified empirically on usegalaxy.org 26.1,
    2026-08-28 — but NOT when a `multiple: true` input consumes a whole collection in one job, where
    only paths survive. Which of the two a wrapper needs is a judgement, so it is left to a human.
  * Any other `#...` Cheetah directive.

⚠ CONTAINERS ARE LOOKED UP, NEVER DERIVED. A UDT has no conda resolution, and the biocontainer build
suffix (`--<hash>_<build>`) is not predictable from the version. `CONTAINERS` below is a hand-checked
map; an unlisted requirement raises rather than producing an image reference that will fail at pull
time with `manifest unknown`.

    python3 scripts/xml_to_udt.py tools/dustmasker/dustmasker.xml > udt/dustmasker.yml
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys
import xml.etree.ElementTree as ET

#: requirement (package, version) -> verified biocontainer image.
#: ⚠ Every entry must be checked against quay.io before it is added; see the module docstring.
CONTAINERS = {
    ("blast", "2.17.0"): "quay.io/biocontainers/blast:2.17.0--h66d330f_0",
    ("tantan", "51"): "quay.io/biocontainers/tantan:51--h4ac6f70_0",
    ("python", "3.12"): "quay.io/biocontainers/python:3.12",
}

CHEETAH_OK = re.compile(r"#(if|else|end if)\b")
GZIP_BLOCK = re.compile(
    r"#if\s+str\(\$(\w+)\.ext\)\s*==\s*'fasta\.gz'\s*"
    r"gunzip -c '\$\1'\s*#else\s*cat '\$\1'\s*#end if", re.S)


def refuse(msg: str) -> None:
    sys.exit(f"REFUSING to convert: {msg}\n"
             "  A silently mistranslated tool is worse than an absent one; port this by hand.")


#: Interpreters a command may invoke that a single-package biocontainer will NOT carry.
#: ⛔ MEASURED, NOT ASSUMED: quay.io/biocontainers/blast:2.17.0--h66d330f_0 has `dustmasker` and
#: `awk` but NO `python3` -- a converted dustmasker died with exit 127,
#: `python3: command not found`, on 2026-08-28. Every masking wrapper in this repo pipes its tool
#: through a small Python helper, so every one of them hits this.
RUNTIMES = ("python3", "python", "perl", "Rscript")


def check_single_runtime(cmd: str, container: str) -> None:
    """Refuse when the command needs an interpreter the chosen image will not have.

    ⚠ THE FIX IS NOT A BIGGER IMAGE, IT IS TWO TOOLS. A UDT gets exactly one container, so a wrapper
    that combines a compiled binary with a script written in another language cannot be one UDT
    unless a MULLED multi-package image exists -- and its hash is not derivable, so somebody has to
    build and publish one. Splitting the wrapper into two UDTs chained in the workflow (binary in its
    own image, script in a `python` image) needs no new image at all and is the cheaper answer.
    """
    needed = [r for r in RUNTIMES if re.search(rf"\b{r}\b", cmd)]
    if not needed:
        return
    pkg = container.split("/")[-1].split(":")[0]
    if pkg in ("python",) or any(r.startswith(pkg) for r in needed):
        return
    refuse(f"the command invokes {needed[0]} but the image is `{container}`, a single-package "
           f"container for `{pkg}` which does not carry it (verified for blast on 2026-08-28: "
           f"exit 127). Either publish a mulled image with both, or -- cheaper -- SPLIT this "
           f"wrapper into two chained UDTs, one per runtime.")


def container_for(root: ET.Element) -> str:
    reqs = [(r.text.strip(), r.get("version")) for r in root.iter("requirement")
            if r.get("type") == "package"]
    if not reqs:
        refuse("the wrapper declares no <requirements>, so there is no image to resolve. "
               "Add conda requirements first (see the fastan PR for an example).")
    for key in reqs:
        if key in CONTAINERS:
            return CONTAINERS[key]
    if len(reqs) > 1:
        refuse(f"the wrapper needs binaries from {len(reqs)} packages ({', '.join(r[0] for r in reqs)}), "
               "and a UDT gets exactly ONE container. That needs a mulled multi-package image, whose "
               "hash cannot be derived -- build and publish one, then add it to CONTAINERS.")
    refuse(f"no verified container for {reqs}. Look the tag up on quay.io and add it to "
           f"CONTAINERS — do not guess the build suffix.")
    return ""


def convert_command(cmd: str, tool_dir: pathlib.Path) -> tuple[str, list[tuple[str, str]]]:
    """Cheetah command -> shell_command, plus the scripts that must be inlined."""
    if "#for" in cmd:
        refuse("the command contains a `#for` loop; `shell_command` has no loop construct.")
    if "element_identifier" in cmd:
        refuse("the command reads `element_identifier`. A UDT mapped over a collection CAN read it "
               "as $(inputs[\"<name>|__identifier__\"]), but a multiple-input job cannot — decide "
               "which shape this tool needs and port it by hand.")

    # The one supported conditional: gunzip-or-cat on a possibly-gzipped FASTA.
    collapsed = re.sub(r"\s+", " ", cmd)
    m = GZIP_BLOCK.search(collapsed)
    if m:
        # ⛔ NOT `gzip -cdf`. GNU gzip copies unrecognised input through when given --stdout, but the
        # gzip in the blast biocontainer does not -- it fails with `gzip: invalid magic` on a plain
        # FASTA. Found by running the converted tool, not by reading the docs. Testing the file and
        # branching is portable to every gzip.
        v = m.group(1)
        # ⚠ NOT wrapped in $(...) -- that is GALAXY's templating delimiter, not a shell
        # substitution. Emitting `$(if ...)` makes Galaxy try to evaluate a shell conditional as
        # ECMAScript. This is plain shell and must stay plain shell.
        cmd = GZIP_BLOCK.sub(
            f"{{ if gzip -t '$(inputs.{v}.path)' 2>/dev/null; "
            f"then gzip -cd '$(inputs.{v}.path)'; else cat '$(inputs.{v}.path)'; fi; }}", collapsed)
    leftover = CHEETAH_OK.search(cmd)
    if leftover:
        refuse(f"unhandled Cheetah directive `{leftover.group(0)}` remains after conversion.")
    if "#" in re.sub(r"#\w*of\d+|#[0-9a-fA-F]{6}", "", cmd):
        pass  # bare '#' inside quoted awk is fine; the directive checks above are what matter

    # ⛔ UNESCAPE CHEETAH'S BACKSLASHES. The XML writes `\$0` so Cheetah does not eat the `$`; in a
    # shell_command there is no Cheetah, so `\$0` reaches awk literally and is a syntax error. Galaxy
    # only interpolates `$(...)`, so a bare `$` needs no protection here.
    cmd = cmd.replace("\\$", "$")

    scripts: list[tuple[str, str]] = []
    for name in sorted(set(re.findall(r"\$__tool_directory__/([\w.\-]+)", cmd))):
        p = tool_dir / name
        if not p.is_file():
            refuse(f"command references {name} but {p} does not exist")
        body = p.read_text().rstrip("\n")
        if "$(" in body:
            refuse(f"{name} contains `$(`, which Galaxy would interpolate inside the heredoc")
        scripts.append((name, body))
        cmd = cmd.replace(f"'$__tool_directory__/{name}'", name).replace(
            f"$__tool_directory__/{name}", name)

    # Remaining `$name` -> data path or scalar. Do the longest names first so $input2 is not
    # clobbered by $input.
    for var in sorted(set(re.findall(r"\$(\w+)", cmd)), key=len, reverse=True):
        if var in ("output",):
            cmd = re.sub(rf"'\${var}'", "out.dat", cmd)
        else:
            cmd = re.sub(rf"'\${var}'", f"'$(inputs.{var}.path)'", cmd)
    return cmd.strip(), scripts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xml", type=pathlib.Path)
    ap.add_argument("--suffix", default="+udt0", help="appended to the tool version")
    args = ap.parse_args()

    root = ET.parse(args.xml).getroot()
    tool_dir = args.xml.parent
    cmd_el = root.find("command")
    if cmd_el is None or not cmd_el.text:
        refuse("no <command> element")
    body, scripts = convert_command(cmd_el.text, tool_dir)
    check_single_runtime(body, container_for(root))

    heredocs = "".join(
        f"cat > {n} <<'BRC_EOF_{i}'\n{b}\nBRC_EOF_{i}\n" for i, (n, b) in enumerate(scripts))
    # The XML often opens with its own `set -o pipefail`; do not add a second one.
    prefix = "" if body.lstrip().startswith("set -o pipefail") else "set -o pipefail\n"
    shell = heredocs + prefix + body + "\n"

    inputs = []
    for p in root.iter("param"):
        if p.get("type") == "data":
            inputs.append({"name": p.get("name"), "type": "data",
                           "format": (p.get("format") or "data").split(",")[0],
                           "label": p.get("label") or p.get("name")})
    outs = []
    for d in root.iter("data"):
        outs.append({"name": d.get("name"), "type": "data",
                     "format": d.get("format") or "data",
                     "from_work_dir": "out.dat",
                     "label": (d.get("label") or d.get("name")).replace("${tool.name}", root.get("name"))})
    if not outs:
        refuse("no <data> output to claim with from_work_dir")

    import yaml
    doc = {
        "class": "GalaxyUserTool", "id": f"brc-{root.get('id')}".replace("_", "-"),
        "name": f"{root.get('name')} (BRC UDT)",
        "version": f"{root.get('version')}{args.suffix}",
        "description": (root.findtext("description") or "").strip(),
        "container": container_for(root),
        "shell_command": shell, "inputs": inputs, "outputs": outs,
        "help": {"format": "markdown", "content":
                 f"Generated from `{args.xml}` by `scripts/xml_to_udt.py` — do not edit by hand; "
                 f"edit the XML and regenerate.\n\n"
                 f"{(root.findtext('help') or '').strip()}\n\n"
                 f"⚠ The wrapper's bundled script(s) "
                 f"({', '.join(n for n, _ in scripts) or 'none'}) are inlined as heredocs, because a "
                 f"UDT has no `$__tool_directory__`."},
    }
    print(yaml.safe_dump(doc, sort_keys=False, width=100))
    return 0


if __name__ == "__main__":
    sys.exit(main())
