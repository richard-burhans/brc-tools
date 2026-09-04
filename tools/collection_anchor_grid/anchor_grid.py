#!/usr/bin/env python3
"""Select and rename the WF-C chains that WF-C2's TOGA2 pass needs.

WF-C emits one cleaned chain per ordered strain pair, keyed `{target}.{query}`
(56 cells for 8 strains). WF-C2's projection grid is anchors x strains with the
anchor self-cells dropped, keyed `{anchor}_{query}` (21 cells for 3 anchors x 8
strains). TOGA2's `--chain_file` is the reference-to-query chain and its
reference is the anchor, so the chain we want for cell `{anchor}_{query}` is the
one whose TARGET is the anchor: `{anchor}.{query}`.

Bridging the two therefore needs a selection list and a rename map:

  keep.txt     `{anchor}.{query}`                       -> __FILTER_FROM_FILE__
  relabel.tsv  `{anchor}.{query}<TAB>{anchor}_{query}`  -> __RELABEL_FROM_FILE__
  order.txt    `{anchor}_{query}`                       -> __SORTLIST__ (sort_type: file)

The third one is not optional. Galaxy pairs collections in a map-over by
POSITION, not by element identifier. The other grids come out of
__CROSS_PRODUCT_FLAT__ and are therefore in anchor-collection order, while a
chain collection filtered out of WF-C's 56 keeps WF-C's own (alphabetical)
order. Without re-sorting, every cell is handed another cell's chain, TOGA2
reports "Processed 0 chains" and exits -- and because the identifiers are
correct it all looks right. order.txt is the cross-product order, so sorting
the chains by it puts the two sides back in step.

All three are pure functions of the two collections' element identifiers, so
nothing has to be hand-authored per panel.
"""
import argparse


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    # THE `-file` FORMS EXIST BECAUSE A UDT CANNOT BUILD THE OTHER ONES. The XML tool turns a
    # collection into a space-separated string with Cheetah, which a User-Defined Tool has no
    # equivalent of: its templating owns the two-character sequence that opens a shell command
    # substitution and interpolates it ANYWHERE in a shell_command, heredocs included, so wrapping
    # one around a Galaxy expression to read the file inline consumes the outer form as an
    # expression and the job dies at command-build with an empty command line and exit=None.
    # Measured on usegalaxy.org (26.1) and a local 25.0; two sibling tools shipped that way and had
    # never run on any Galaxy. Passing the PATH and reading it here needs no substitution at all.
    #
    # THIS COMMENT MUST NOT CONTAIN THE SEQUENCE IT DESCRIBES: the file is inlined into the UDT's
    # heredoc verbatim, so a literal example here would be interpolated exactly like real code, and
    # the generator's guard rejects the file outright when one appears.
    p.add_argument("--anchors", help="Space-separated anchor element identifiers.")
    p.add_argument("--anchors-file", help="File of anchor element identifiers, one per line.")
    p.add_argument("--strains", help="Space-separated strain (query) element identifiers.")
    p.add_argument("--strains-file", help="File of strain element identifiers, one per line.")
    p.add_argument("--keep", required=True, help="Output: ids to keep, one per line.")
    p.add_argument("--relabel", required=True, help="Output: 2-column rename map.")
    p.add_argument("--order", required=True,
                   help="Output: cross-product order, one {anchor}_{query} per line.")
    a = p.parse_args(argv)
    for flag in ("anchors", "strains"):
        inline, from_file = getattr(a, flag), getattr(a, flag + "_file")
        if inline is None and from_file is None:
            p.error(f"one of --{flag} or --{flag}-file is required")
        # BOTH IS AN ERROR, NOT A PRECEDENCE. Silently preferring one meant an empty string beat a
        # populated file and produced a zero-row grid -- which selects no chains and renames
        # nothing, without a word. `is None`, not falsiness: an EMPTY collection legitimately
        # renders as an empty string, and that is a request for zero identifiers, not an omission.
        if inline is not None and from_file is not None:
            p.error(f"give --{flag} or --{flag}-file, not both -- which one wins is not something "
                    f"to guess at")
    return a


def identifiers(inline, path, what):
    """Element identifiers from whichever form was given, validated.

    utf-8-SIG, NOT utf-8: a BOM on a hand-made file survives into the first name and then matches
    no collection element, silently dropping exactly one row. And a TAB in an identifier breaks the
    column contract of the relabel map below -- `{a}.{b}<TAB>{a}_{b}` becomes a three-column row
    and __RELABEL_FROM_FILE__ reads the wrong field.
    """
    ids = ([x for x in inline.split() if x] if inline is not None
           else [x.strip() for x in open(path, encoding="utf-8-sig") if x.strip()])
    bad = [i for i in ids if "\t" in i]
    if bad:
        raise SystemExit(f"{what} identifier(s) contain a tab, which breaks the relabel map's "
                         f"column contract: {bad[:3]}")
    if len(set(ids)) != len(ids):
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        raise SystemExit(f"duplicate {what} identifier(s) {dupes[:3]}: the grid is a cross product, "
                         f"so a repeat produces duplicate cell names and the collection operations "
                         f"downstream silently keep only one of them.")
    return ids


def grid(anchors, strains):
    """Ordered (anchor, query) pairs, anchor self-cells dropped."""
    return [(a, q) for a in anchors for q in strains if q != a]


def main(argv=None):
    args = parse_args(argv)
    anchors = identifiers(args.anchors, args.anchors_file, "anchor")
    strains = identifiers(args.strains, args.strains_file, "strain")
    if not anchors:
        raise SystemExit("no anchor element identifiers")
    if not strains:
        raise SystemExit("no strain element identifiers")
    missing = [a for a in anchors if a not in strains]
    if missing:
        # Anchors are a subset of the panel; if one is absent the chain it needs
        # was never produced by WF-C, and TOGA2 would fail on a missing element.
        raise SystemExit(f"anchors absent from the strain collection: {missing}")

    pairs = grid(anchors, strains)
    with open(args.keep, "w") as kf, open(args.relabel, "w") as rf, open(args.order, "w") as of:
        for a, q in pairs:
            kf.write(f"{a}.{q}\n")
            rf.write(f"{a}.{q}\t{a}_{q}\n")
            of.write(f"{a}_{q}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
