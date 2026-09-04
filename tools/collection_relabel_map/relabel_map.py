#!/usr/bin/env python3
"""Emit `{a}_{b}\t{a}.{b}` for every ordered pair of collection element identifiers.
The cross-product cells are named `A_B` (underscore join); downstream Phase E expects
`A.B`. This 2-col TSV drives WF-C's __RELABEL_FROM_FILE__ step."""
import argparse, itertools
ap = argparse.ArgumentParser()
# ⛔ `--ids-file` EXISTS BECAUSE A UDT CANNOT BUILD `--ids`. The UDT template used a shell
# command substitution wrapped around a Galaxy expression to turn the identifier FILE into a
# space-separated string. Galaxy owns that syntax in a shell_command and interpolates it ANYWHERE,
# heredocs included, so the outer substitution is consumed as a Galaxy expression, the command line
# cannot be built at all, and the job fails with exit=None and an empty command before anything
# runs. Measured on BOTH usegalaxy.org (26.1) and a local 25.0: "Error occurred while building
# command line". Those two tools had therefore never run on any Galaxy.
#
# ⚠ scripts/build_softmask_udts.py already refuses that construct -- its FORBIDDEN constant names
# the two-character opening sequence -- and build_inventory_udts.py did not. Passing the PATH and
# reading the file here needs no substitution at all, so nothing can be interpolated.
#
# ⚠ AND THIS COMMENT MUST NOT CONTAIN THE SEQUENCE IT DESCRIBES. This script is INLINED into the
# UDT's heredoc, so a literal example here would be interpolated exactly like real code -- and the
# generator's own guard rejected an earlier draft of this comment for precisely that reason.
ap.add_argument("--ids", help="space-separated element identifiers")
ap.add_argument("--ids-file", help="file of element identifiers, one per line")
ap.add_argument("--out", required=True)
a = ap.parse_args()
# ⚠ `is None`, NOT falsiness. The classic XML wrappers always pass --ids, and an EMPTY collection
# renders `--ids ''` -- which is a valid request for zero identifiers, and used to produce a
# zero-row file. Testing truthiness turned that into "one of --ids or --ids-file is required" and
# exit 2, a regression against the merged tool for the one input most likely to be automated.
if a.ids is None and a.ids_file is None:
    ap.error("one of --ids or --ids-file is required")
# ⛔ BOTH IS AN ERROR, NOT A PRECEDENCE. Silently preferring one meant `--ids "" --ids-file real.txt`
# produced a ZERO-ROW file from a fully populated list -- and a zero-row self-pair list removes
# nothing while a zero-row relabel map renames nothing, both without a word. Refuse the ambiguity.
if a.ids is not None and a.ids_file is not None:
    ap.error("give --ids or --ids-file, not both -- which one wins is not something to guess at")
# ⚠ `is not None` HERE TOO. `--ids ''` is falsy but PRESENT, and testing truthiness sent it to
# the file branch with ids_file=None -- a TypeError instead of the zero-row file the empty-string
# case is supposed to produce. The guard above and this selection must agree on what "given" means.
# ⚠ utf-8-SIG, NOT utf-8. A BOM on a hand-made identifier file survives into the first name --
# `\ufeffcs10` -- which then matches no collection element, silently no-opping exactly one row.
ids = ([x for x in a.ids.split() if x] if a.ids is not None
       else [x.strip() for x in open(a.ids_file, encoding="utf-8-sig") if x.strip()])
# ⛔ A TAB IN AN IDENTIFIER BREAKS THE COLUMN CONTRACT DOWNSTREAM. relabel_map emits
# `{a}_{b}<TAB>{a}.{b}`, so a tab inside a name yields a THREE-column row and
# `__RELABEL_FROM_FILE__` reads the wrong field. --ids could never carry one (it splits on
# whitespace); --ids-file can, so the check belongs here.
_bad = [i for i in ids if "\t" in i]
if _bad:
    raise SystemExit(f"identifier(s) contain a tab, which breaks the output's column contract: "
                     f"{_bad[:3]}")
with open(a.out, "w") as f:
    for x, y in itertools.product(ids, ids):
        f.write(f"{x}_{y}\t{x}.{y}\n")
