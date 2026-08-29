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
does not exist.

⛔ THE REFUSAL IS ENFORCED BY A COMPLETENESS ASSERTION, NOT BY A LIST OF KNOWN-BAD CONSTRUCTS.
`assert_fully_translated()` rejects the result if ANY `$name` / `${name}` / `#directive` survives
conversion. This is the load-bearing check and it is deliberately the last thing that runs: an
enumeration of what to refuse can only ever cover the constructs somebody thought of, and the three
that were missed all produced tools that ran and were wrong ---

  * `#set $in_ext = ...` was not in the refusal list, so it passed through into `shell_command`,
    where a leading `#` is a SHELL COMMENT. `${in_ext}` was then undefined and
    `ln -sf … 'input.${in_ext}'` created a file called `input.` (tools/longdust, tools/sdust).
  * Only SINGLE-QUOTED `'$var'` was rewritten, so `longdust -k$k -w$w -t$t` reached the container
    verbatim and every flag ran with an empty argument.
  * Only an output literally named `output` was recognised, so tools/sdust's `out_bed` became
    `$(inputs.out_bed.path)` --- the tool's OUTPUT addressed as an INPUT path that does not exist.

None of the three raised anything. An assertion that nothing untranslated survives catches all of
them, and the next one too.

SUPPORTED
  * `#if str($x.ext) == 'fasta.gz' / #else / #end if` around a gunzip-or-cat — the only Cheetah
    conditional these wrappers use. Rewritten as a portable test-and-branch.
  * `$__tool_directory__/<script>` — the script is INLINED as a heredoc. A UDT has no tool
    directory, and these scripts are 364 B to 2.7 kB, far too small to justify an image.
  * `<param type="data">` -> `$(inputs.<name>.path)`.
  * `<param>` scalars (integer/float/text/select) -> `$(inputs.<name>)`.
  * `<param type="boolean" truevalue= falsevalue=>` -> `$(inputs.<name> ? '<t>' : '<f>')`, which is
    the UDT idiom; `truevalue`/`falsevalue` have no meaning in a UDT.
  * `<data>` outputs -> a work-dir filename claimed by `from_work_dir`.

REFUSED (raises, naming the offending construct)
  * `#for` loops, `#set`, and every other Cheetah directive.
  * `element_identifier`. ⚠ Galaxy DOES expose it to a UDT as `$(inputs["<name>|__identifier__"])`
    when the tool is MAPPED over a collection — verified empirically on usegalaxy.org 26.1,
    2026-08-28 — but NOT when a `multiple: true` input consumes a whole collection in one job, where
    only paths survive. Which of the two a wrapper needs is a judgement, so it is left to a human.
  * Anything at all that the translation did not consume (see the assertion above).

⛔ A UDT ON usegalaxy.org GETS ONE CORE, AND THERE IS NO KNOB. Measured 2026-08-28 across five
controlled runs of one byte-identical probe: with no `requirements` block the job runs and reports
`GALAXY_SLOTS=1`, `nproc` 1, `sched_getaffinity` 1; adding a `resource` requirement makes it
UNRUNNABLE -- "No destinations are available to fulfill request: user_defined-.*" -- and that holds
for `cores_min: 4`, `cores_min: 2`, `ram_min: 8192` and even `cores_min: 1`. It is not a ceiling
being clipped, it is the requirement being refused. `export GALAXY_SLOTS=8` inside the command does
change the number the tool reads, and changes nothing else: `nproc` stays 1, so the tool would
oversubscribe a single core. Threading flags in a converted wrapper are therefore cosmetic on the
public server, and a wrapper's own `${VAR:-default}` is escaped through verbatim so that
the SAME yaml stays correct on a Galaxy that does allocate cores (this repo's own
`galaxy_config_job_conf.xml` hands every tool 16 slots and 60 GB).

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

#: Cheetah directive keywords. Only the `#if`/`#else`/`#end if` gzip idiom is handled; every other
#: one is a refusal. ⛔ `#set` is in this list because its ABSENCE was the first silent
#: mistranslation (see the module docstring) --- but the completeness assertion, not this list, is
#: what makes the refusal trustworthy.
CHEETAH_DIRECTIVES = (
    "if", "else", "elif", "end", "for", "while", "repeat", "unless", "set", "echo", "silent",
    "import", "from", "def", "block", "include", "raw", "try", "except", "finally", "pass",
    "break", "continue", "attr", "compiler", "call", "filter", "assert", "del", "return",
)
CHEETAH_RE = re.compile(rf"#\s*({'|'.join(CHEETAH_DIRECTIVES)})\b")

GZIP_BLOCK = re.compile(
    r"#if\s+str\(\$(\w+)\.ext\)\s*==\s*'fasta\.gz'\s*"
    r"gunzip -c '\$\1'\s*#else\s*cat '\$\1'\s*#end if", re.DOTALL)

#: OUR OWN emissions only -- `$(inputs...)`. ⚠ Deliberately NOT a general `\$\([^)]*\)`: blanking
#: every `$( )` let a wrapper's genuine shell command substitution through the completeness check,
#: which contradicts this module's own rule for inlined scripts (`convert_command` refuses `$(` in a
#: heredoc because Galaxy interpolates it). Galaxy claims `$( )` for templating; a shell
#: substitution in a `shell_command` is not shell, it is a Galaxy expression that will not evaluate.
GALAXY_EXPR = re.compile(r"\$\(inputs[^)]*\)")
#: A surviving template reference: `$name` or `${name}`, starting with a letter or underscore.
#: ⚠ `$0`/`$1` are NOT matched --- those are awk positionals, legitimate after the un-escaping
#: below, and flagging them would refuse every wrapper that pipes through awk.
SURVIVING_VAR = re.compile(r"\$\{?([A-Za-z_]\w*)")

#: Environment variables the UDT job container exports, which pass through a `shell_command`
#: unchanged and are NOT untranslated templates.
#:
#: ⛔ MEASURED ON usegalaxy.org, NOT ASSUMED. `udt/env_probe.gxtool.yml` dumps `printenv` from
#: inside a real UDT job; run it again if this list is ever in doubt. As of 2026-08-28 (job
#: 79420986, an Apptainer container on the main cluster) the exported set was exactly these plus
#: Apptainer/Singularity bookkeeping, and `GALAXY_SLOTS` was 1, matching `nproc`.
#:
#: ⚠ `USER` and `SHELL` are DELIBERATELY ABSENT, for two DIFFERENT reasons, and the distinction is
#: the point. Neither appears in `printenv`. `$USER` therefore expands to the empty string. `$SHELL`
#: does NOT -- it read `/bin/bash` -- but that value comes from the invoking bash, not from Galaxy,
#: so it is not a platform guarantee and does not belong on a list of things Galaxy provides. Both
#: were on this list until the probe was run, purely because they look like variables that ought to
#: exist.
SHELL_ENV_OK = frozenset({
    "TMPDIR", "TMP", "TEMP", "HOME", "PWD", "PATH", "LANG", "LC_ALL",
    "GALAXY_SLOTS", "GALAXY_MEMORY_MB", "GALAXY_MEMORY_MB_PER_SLOT",
    "_GALAXY_JOB_TMP_DIR", "_GALAXY_JOB_HOME_DIR",
})

#: `${...}` in a `shell_command`.
#:
#: ⛔ AN EARLIER VERSION OF THIS FILE SAID THE BRACE FORM WAS FATAL. IT IS NOT, AND THE MISTAKE WAS
#: EXPENSIVE -- it drove a refusal rule and a whole rewrite pass that were not needed. `${...}` is
#: Galaxy's own ECMAScript FUNCTION-BODY expression form, the sibling of `$(...)`:
#: `${ return 'HELLO'; }` evaluates to HELLO. `${HOME}` fails only because `HOME` is not a valid JS
#: function body -- a parse error, surfacing as "Error occurred while building command line for
#: tool" with both job streams empty, which is what made it look like a ban.
#:
#: ⚠ A SHELL BRACE EXPANSION NEEDS ESCAPING, NOT REFUSING. Measured on usegalaxy.org:
#:     \${GALAXY_SLOTS:-4}  ->  1          \${HOME}  ->  /home/g2main
#: so the wrapper's own `:-default` survives verbatim and needs no reconstruction.
BRACE_EXPANSION = re.compile(r"(?<!\\)\$\{")

#: A `${ ... }` whose body reads as JavaScript rather than a shell variable. Escaping one of these
#: would break a legitimate Galaxy expression, so the converter refuses rather than guessing.
JS_FUNCTION_BODY = re.compile(r"\$\{[^}]*\breturn\b[^}]*\}")


def escape_shell_braces(cmd: str) -> str:
    """Escape every unescaped `${` so Galaxy hands it to the shell instead of parsing it as JS."""
    out, i = [], 0
    while i < len(cmd):
        if cmd.startswith("${", i) and (i == 0 or cmd[i - 1] != "\\"):
            out.append("\\${")
            i += 2
        else:
            out.append(cmd[i])
            i += 1
    return "".join(out)



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
    if pkg == "python" or any(r.startswith(pkg) for r in needed):
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


def param_name(el: ET.Element) -> str | None:
    """The name Galaxy will address a param by.

    ⚠ `name` IS OPTIONAL. tools/longdust declares every tunable as `<param argument="-k" …>` with no
    `name`, and Galaxy derives `k` from the argument. A converter that reads only `name` sees nine
    unnamed params, translates none of them, and emits `longdust -k$k -w$w …` verbatim.
    """
    if el.get("name"):
        return el.get("name")
    arg = el.get("argument")
    if arg:
        return arg.lstrip("-").replace("-", "_")
    return None


def collect_io(root: ET.Element) -> tuple[dict[str, ET.Element], dict[str, ET.Element]]:
    """Params and outputs, scoped to `<inputs>`/`<outputs>`.

    ⚠ SCOPED DELIBERATELY. `root.iter("param")` also walks `<test>` blocks, where params carry a
    `name` and a `value` but no `type` -- so a test fixture can shadow a real param, and the count
    of "inputs" silently includes rows that are not inputs at all.
    """
    params: dict[str, ET.Element] = {}
    inputs_el = root.find("inputs")
    if inputs_el is not None:
        for el in inputs_el.iter("param"):
            name = param_name(el)
            if name:
                params[name] = el
    outputs: dict[str, ET.Element] = {}
    outputs_el = root.find("outputs")
    if outputs_el is not None:
        for el in outputs_el.iter("data"):
            if el.get("name"):
                outputs[el.get("name")] = el
    return params, outputs


def clean_label(raw: str, tool_name: str) -> str:
    """Render an output label, refusing any `${...}` the UDT cannot resolve.

    ⛔ THE COMMAND IS NOT THE ONLY PLACE `${...}` IS FATAL. Every masking wrapper here writes
    `label="${tool.name} on ${on_string}: ..."`, and only `${tool.name}` was ever rewritten -- so
    `${on_string}` reached the emitted YAML verbatim. `assert_fully_translated()` never saw it
    because it inspects the command alone. Galaxy claims the brace form for its own templating (see
    BRACE_EXPANSION), so at best the label reads as literal garbage and at worst the tool will not
    register, with both job streams empty and nothing naming the cause.

    `${on_string}` is dropped rather than translated: it names the inputs a classic job ran on, and
    a UDT has no equivalent. Anything else raises rather than being guessed at.
    """
    label = raw.replace("${tool.name}", tool_name)
    # "X on ${on_string}: Y" -> "X: Y"; a bare "${on_string}" just goes.
    label = re.sub(r"\s+on\s+\$\{on_string\}", "", label)
    label = label.replace("${on_string}", "").strip()
    leftover = re.search(r"\$\{[^}]*\}", label)
    if leftover:
        refuse(f"the output label {raw!r} contains `{leftover.group(0)}`, which a UDT cannot "
               f"resolve -- the brace form is Galaxy's own templating and is fatal in a UDT. "
               f"Simplify the label in the XML, or port this tool by hand.")
    return label


def workdir_file(name: str, el: ET.Element) -> str:
    """The work-dir filename a UDT output is claimed from.

    Honours the XML's own `from_work_dir` when it declares one -- tools/build_trackdb writes
    `selection_strict.bed` and never redirects to `$output`, so hardcoding one filename for every
    output both loses that name and gives every output of a multi-output wrapper the same source.
    """
    return el.get("from_work_dir") or f"{name}.dat"


def substitute(cmd: str, params: dict[str, ET.Element], outputs: dict[str, ET.Element]) -> str:
    """Replace every `$name`/`${name}` the XML declares. Longest first, so `$input2` survives `$input`."""
    for name in sorted(set(params) | set(outputs), key=len, reverse=True):
        if name in outputs:
            repl = workdir_file(name, outputs[name])
            quoted = f"'{repl}'"
        else:
            el = params[name]
            ptype = el.get("type")
            if ptype == "data":
                quoted = repl = f"'$(inputs.{name}.path)'"
            elif ptype == "boolean":
                tv, fv = el.get("truevalue", ""), el.get("falsevalue", "")
                if "'" in tv or "'" in fv:
                    refuse(f"boolean param `{name}` has a quote in its truevalue/falsevalue, which "
                           f"cannot be embedded in the UDT ternary; port this one by hand.")
                quoted = repl = f"$(inputs.{name} ? '{tv}' : '{fv}')"
            else:
                # ⚠ THE QUOTED FORM KEEPS ITS QUOTES. Only the boolean ternary above may drop them,
                # because an empty `falsevalue` has to vanish from the command line entirely. For a
                # text/select/number param the quotes in the XML are load-bearing: a value of
                # `-r 0.01` would word-split into two argv entries without them, and one containing
                # `;` or a backtick would be arbitrary shell in the job container.
                quoted = f"'$(inputs.{name})'"
                repl = f"$(inputs.{name})"
        # Quoted forms first, so the quotes are consumed rather than left wrapping the expression.
        cmd = cmd.replace(f"'${name}'", quoted).replace(f"'${{{name}}}'", quoted)
        cmd = re.sub(rf"\$\{{{re.escape(name)}\}}", repl, cmd)
        cmd = re.sub(rf"\${re.escape(name)}\b", repl, cmd)
    return cmd


def assert_fully_translated(cmd: str, extra_ok: set[str] | None = None) -> None:
    """⛔ THE LOAD-BEARING CHECK. Refuse if anything untranslated survived.

    Every silent mistranslation this converter has produced was a construct that no refusal rule
    named. This inverts that: instead of enumerating what is forbidden, require that the output
    contain nothing that still looks like a template. See the module docstring for the three.
    """
    directive = CHEETAH_RE.search(cmd)
    if directive:
        line = next((ln.strip() for ln in cmd.splitlines() if directive.group(0) in ln), "")
        refuse(f"the Cheetah directive `{directive.group(0)}` survived conversion — in a "
               f"`shell_command` a leading `#` is a COMMENT, so this would be silently dropped and "
               f"any variable it defines would be empty.\n"
               f"  offending line: {line}\n"
               f"  `shell_command` has no templating; restructure the XML or port this by hand.")
    # Blank out our own `$(inputs...)` emissions before looking for leftovers.
    residue = GALAXY_EXPR.sub("", cmd)
    stray = re.search(r"\$\([^)]*\)", residue)
    if stray:
        refuse(f"`{stray.group(0)}` is a shell command substitution, and Galaxy claims `$( )` for "
               f"its own templating -- it will be interpreted as an expression, not run by the "
               f"shell. This module already refuses `$(` inside an inlined script for the same "
               f"reason. Rewrite it with backticks, which pass through untouched, or port by hand.")
    brace = BRACE_EXPANSION.search(residue)
    if brace:  # escape_shell_braces() ran first, so an unescaped ${ here is this converter's bug
        inner = brace.group(0)[2:-1]
        bare = inner.split(":")[0].split("-")[0] or "VAR"
        refuse(f"`{brace.group(0)}` uses shell brace expansion, which Galaxy claims for its own "
               f"templating. The job then fails with \"Error occurred while building command line\" "
               f"and BOTH streams empty — measured on usegalaxy.org 2026-08-28.\n"
               f"  Write it bare instead: `${bare}`. The bare form reaches the shell untouched and "
               f"the container does export it (GALAXY_SLOTS read 1).\n"
               f"  ⚠ A `:-default` cannot survive this rewrite; drop it, or port by hand.")
    for match in SURVIVING_VAR.finditer(residue):
        var = match.group(1)
        if var in SHELL_ENV_OK or var in (extra_ok or ()):
            continue  # exported by the job container, or defined by our own preamble
        if var.startswith("GALAXY_"):
            refuse(f"`${var}` survived conversion. `udt/env_probe.gxtool.yml` measured which "
                   f"GALAXY_* variables a UDT container exports and this is not one of them, so it "
                   f"would expand to the empty string. Re-run the probe if you think that has "
                   f"changed, or substitute a literal in the XML.")
        refuse(f"`${var}` survived conversion — it matches no `<param>` or `<data>` in "
               f"this XML, so it would reach the container as an undefined shell variable and expand "
               f"to the empty string. Declare it, or port this tool by hand.")


def convert_command(cmd: str, tool_dir: pathlib.Path, params: dict[str, ET.Element],
                    outputs: dict[str, ET.Element]) -> tuple[str, list[tuple[str, str]]]:
    """Cheetah command -> shell_command, plus the scripts that must be inlined."""
    if "element_identifier" in cmd:
        refuse("the command reads `element_identifier`. A UDT mapped over a collection CAN read it "
               "as $(inputs[\"<name>|__identifier__\"]), but a multiple-input job cannot — decide "
               "which shape this tool needs and port it by hand.")

    # The one supported conditional: gunzip-or-cat on a possibly-gzipped FASTA.
    #
    # ⛔ DO NOT COLLAPSE WHITESPACE TO MAKE THIS MATCH. An earlier version ran
    # `re.sub(r"\s+", " ", cmd)` over the WHOLE command and then substituted into the collapsed
    # copy, which silently joined every newline-separated statement into one line: a wrapper whose
    # command continued `> in.fa` / `tantan in.fa > '$output'` on later lines emitted
    # `{ ... fi; } > in.fa tantan in.fa > 'output.dat'`, one command in which `tantan` and its
    # arguments are extra operands of the brace group. It exits 0 and writes an empty output, and
    # `assert_fully_translated()` CANNOT SEE IT because nothing untranslated survives -- exactly the
    # class of failure this module exists to prevent.
    #
    # The collapse was never necessary: GZIP_BLOCK already spans newlines via `\s*`/`\s+`, so it
    # matches the original text directly. Verified against tools/dustmasker/dustmasker.xml.
    m = GZIP_BLOCK.search(cmd)
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
            f"then gzip -cd '$(inputs.{v}.path)'; else cat '$(inputs.{v}.path)'; fi; }}", cmd)

    # ⛔ UNESCAPE CHEETAH'S BACKSLASHES. The XML writes `\$0` so Cheetah does not eat the `$`; in a
    # shell_command there is no Cheetah, so `\$0` reaches awk literally and is a syntax error. Galaxy
    # only interpolates `$(...)`, so a bare `$` needs no protection here.
    cmd = cmd.replace("\\$", "$")

    # ⚠ ESCAPE shell brace expansions; do not refuse them and do not rebuild their defaults. A
    # wrapper's `${GALAXY_SLOTS:-4}` is shell, and `\${...}` reaches the shell intact, default and
    # all. Refuse only what looks like a real Galaxy expression, since escaping that would break it.
    js = JS_FUNCTION_BODY.search(cmd)
    if js:
        refuse(f"the command contains what reads as a Galaxy expression, {js.group(0)[:60]!r}. "
               f"This converter escapes `${{...}}` as shell and will not guess which of the two was "
               f"meant. Port this tool by hand.")
    cmd = escape_shell_braces(cmd)

    scripts: list[tuple[str, str]] = []
    for name in sorted(set(re.findall(r"\$__tool_directory__/([\w.\-]+)", cmd))):
        p = tool_dir / name
        if not p.is_file():
            refuse(f"command references {name} but {p} does not exist")
        body = p.read_text(encoding="utf-8").rstrip("\n")
        if "$(" in body:
            refuse(f"{name} contains `$(`, which Galaxy would interpolate inside the heredoc")
        scripts.append((name, body))
        cmd = cmd.replace(f"'$__tool_directory__/{name}'", name).replace(
            f"$__tool_directory__/{name}", name)

    cmd = substitute(cmd, params, outputs)
    assert_fully_translated(cmd)
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
    params, out_els = collect_io(root)
    if not out_els:
        refuse("no <data> output to claim with from_work_dir")
    body, scripts = convert_command(cmd_el.text, tool_dir, params, out_els)
    check_single_runtime(body, container_for(root))

    heredocs = "".join(
        f"cat > {n} <<'BRC_EOF_{i}'\n{b}\nBRC_EOF_{i}\n" for i, (n, b) in enumerate(scripts))
    # The XML often opens with its own `set -o pipefail`; do not add a second one.
    prefix = "" if body.lstrip().startswith("set -o pipefail") else "set -o pipefail\n"
    shell = heredocs + prefix + body + "\n"

    inputs = []
    for name, el in params.items():
        ptype = el.get("type")
        if ptype == "data":
            # ⚠ The WHOLE format list, not the first entry. `format="fasta,fasta.gz"` split to
            # `fasta` makes a gzipped genome unselectable in the form -- and every masking wrapper
            # here accepts both.
            inputs.append({"name": name, "type": "data",
                           "format": el.get("format") or "data",
                           "label": el.get("label") or name})
            continue
        spec = {"name": name, "type": ptype or "text", "label": el.get("label") or name}
        if ptype == "boolean":
            # `checked` is the XML spelling; a UDT boolean carries a plain `value`.
            spec["value"] = el.get("checked", "false") == "true"
        else:
            # Emit numbers as numbers. A quoted `value: '7'` on an integer param is a string in the
            # generated YAML, and the form then rejects its own default.
            cast = {"integer": int, "float": float}.get(ptype)
            for key in ("value", "min", "max"):
                raw = el.get(key)
                if raw is None:
                    continue
                spec[key] = cast(raw) if cast else raw
        opts = list(el.iter("option"))
        if opts:
            spec["options"] = [{"label": (o.text or o.get("value")).strip(),
                                "value": o.get("value"),
                                "selected": o.get("selected") == "true"} for o in opts]
        inputs.append(spec)

    outs = [{"name": name, "type": "data",
             "format": el.get("format") or "data",
             "from_work_dir": workdir_file(name, el),
             "label": clean_label(el.get("label") or name, root.get("name") or "")}
            for name, el in out_els.items()]

    # ⚠ DEFERRED DELIBERATELY, not an oversight. PyYAML is needed only to EMIT, and every refusal
    # path above exits before reaching here. The repo installs it explicitly in CI
    # (.github/workflows/pages.yml) but the ambient python3 on a dev box generally does not have it,
    # and for the wrappers in this repo REFUSING is the common outcome -- all five masking tools
    # refuse today. Importing at module scope turns every one of those clean refusals into a
    # ModuleNotFoundError traceback, which is a worse tool for the case that actually happens.
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
