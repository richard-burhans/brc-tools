#!/usr/bin/env python3
"""PyEuk distance + clustering driver for the Galaxy tool ``cyclospora_pyeuk``.

Reads a CDC-style haplotype data sheet (HDS: rows = specimens, columns =
haplotype names, ``X`` = present, empty = absent, an entirely empty locus block
= NOT CALLED), computes the binary KING-wIBS distance matrix with
``PyEukDistanceEngine.compute_revised_wibs_matrix`` and partitions the
specimens with ``CyclosporaClusterFinder.find_clusters``.

Two calls, no reimplementation: everything numeric is PyEuk's, pinned to
commit d8e45a81 (package ``__version__`` 2.1.0, binary KING-wIBS engine). This
file only marshals files in and out and prints diagnostics.

Label-free is the default: ``find_clusters(matrix, None)`` selects k from the
dendrogram merge-height gap knee and never sees an epidemiological label. Pass
``--gold`` only when a genuine gold-standard file exists; it switches PyEuk to
the supervised threshold calibration path.
"""

import argparse
import os
import sys

import pandas as pd

from cyclospora_pyeuk.clustering import CyclosporaClusterFinder
from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine

# The engine's own default. It parameterises the Bayesian distance only; the
# revised KING-wIBS matrix this tool computes does not read it. Kept at the
# reference value rather than exposed as a knob that would do nothing.
EPSILON = 0.3072


def die(msg):
    sys.stderr.write("cyclospora_pyeuk: %s\n" % msg)
    sys.exit(1)


def read_sheet(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:  # noqa: BLE001 - surface the parse failure verbatim
        die("could not parse the haplotype sheet as TSV: %s" % exc)

    df.columns = [str(c).strip() for c in df.columns]

    if df.shape[1] < 2:
        die(
            "the haplotype sheet has %d column(s); expected 'Seq_ID' plus at "
            "least one haplotype column" % df.shape[1]
        )
    if "Seq_ID" not in df.columns:
        die(
            "the haplotype sheet has no 'Seq_ID' column (first column is %r). "
            "PyEuk keys every specimen on 'Seq_ID'; build the sheet with "
            "cyclospora_hds_sheet or rename the column." % df.columns[0]
        )

    df["Seq_ID"] = df["Seq_ID"].astype(str).str.strip()
    df = df[(df["Seq_ID"] != "") & (df["Seq_ID"].str.lower() != "nan")].copy()
    if df.empty:
        die("the haplotype sheet has no specimen rows")

    dup = df["Seq_ID"][df["Seq_ID"].duplicated()].unique().tolist()
    if dup:
        die("duplicate Seq_ID values in the haplotype sheet: %s" % ", ".join(dup[:10]))

    marker_cols = [c for c in df.columns if c != "Seq_ID"]
    values = set()
    for c in marker_cols:
        values.update(v for v in df[c].dropna().astype(str).unique())
    odd = sorted(v for v in values if v.strip() not in ("X", ""))
    if odd:
        sys.stderr.write(
            "cyclospora_pyeuk: warning: %d marker value(s) other than 'X' or "
            "empty present and treated as ABSENT: %s\n"
            % (len(odd), ", ".join(repr(v) for v in odd[:10]))
        )

    return df, marker_cols


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sheet", required=True, help="CDC HDS wide sheet (TSV)")
    ap.add_argument("--matrix-out", required=True, help="output wIBS distance matrix (TSV)")
    ap.add_argument("--clusters-out", required=True, help="output cluster assignments (TSV)")
    ap.add_argument("--gold", default=None, help="optional gold-standard labels; omit for label-free")
    ap.add_argument("--min-completeness", type=float, default=0.10)
    ap.add_argument("--stringency", type=float, default=95.0)
    ap.add_argument("--robust", choices=["true", "false"], default="true")
    ap.add_argument("--default-threshold", type=float, default=0.05)
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=50)
    ap.add_argument("--report-excluded", choices=["true", "false"], default="true")
    ap.add_argument("--float-format", default="%.10g")
    args = ap.parse_args()

    if args.gold is not None and not os.path.exists(args.gold):
        die("gold-standard file not found: %s" % args.gold)
    if args.k_min < 2:
        die("--k-min must be >= 2; PyEuk guards against a k=1 single-cluster collapse")
    if args.k_max < args.k_min:
        die("--k-max (%d) is below --k-min (%d)" % (args.k_max, args.k_min))

    df, marker_cols = read_sheet(args.sheet)
    all_ids = df["Seq_ID"].tolist()
    mode = "supervised" if args.gold else "label-free"

    print("[cyclospora_pyeuk] mode              : %s" % mode)
    print("[cyclospora_pyeuk] specimens in sheet: %d" % len(all_ids))
    print("[cyclospora_pyeuk] haplotype columns : %d" % len(marker_cols))
    sys.stdout.flush()

    engine = PyEukDistanceEngine(epsilon=EPSILON, min_completeness=args.min_completeness)
    matrix = engine.compute_revised_wibs_matrix(df)
    sys.stdout.flush()

    if matrix.shape[0] == 0:
        die(
            "no specimen passed the completeness filter (--min-completeness %g). "
            "Every specimen has a smaller fraction of called haplotype columns "
            "than that; lower the threshold or check the sheet." % args.min_completeness
        )

    matrix.to_csv(args.matrix_out, sep="\t", index_label="Seq_ID",
                  float_format=args.float_format)

    finder = CyclosporaClusterFinder(
        stringency=args.stringency,
        robust=(args.robust == "true"),
        default_threshold=args.default_threshold,
    )
    report_excluded = args.report_excluded == "true"
    clusters, k, threshold = finder.find_clusters(
        matrix,
        gold_file_path=args.gold,
        k_min=args.k_min,
        k_max=args.k_max,
        output_dir=os.path.join(os.getcwd(), "pyeuk_clusters"),
        all_input_ids=all_ids if report_excluded else None,
    )
    sys.stdout.flush()

    # find_clusters returns early for a single specimen, before it can append
    # the excluded ones. Keep the promise of the option in that corner too.
    if report_excluded:
        missing = [s for s in all_ids if s not in set(clusters["Seq_ID"])]
        if missing:
            clusters = pd.concat(
                [clusters, pd.DataFrame({"Seq_ID": missing, "Assigned_cluster": -1})],
                ignore_index=True,
            )

    clusters.to_csv(args.clusters_out, sep="\t", index=False)

    assigned = clusters[clusters["Assigned_cluster"] != -1]
    sizes = assigned["Assigned_cluster"].value_counts().sort_index().to_dict()
    excluded = int((clusters["Assigned_cluster"] == -1).sum())
    print("[cyclospora_pyeuk] specimens clustered : %d" % len(assigned))
    print("[cyclospora_pyeuk] specimens excluded  : %d%s"
          % (excluded, " (reported as cluster -1)" if report_excluded else " (dropped)"))
    print("[cyclospora_pyeuk] k                   : %d" % k)
    print("[cyclospora_pyeuk] threshold           : %.6g" % threshold)
    print("[cyclospora_pyeuk] cluster sizes       : %s"
          % ", ".join("%s=%d" % (c, n) for c, n in sizes.items()))


if __name__ == "__main__":
    main()
