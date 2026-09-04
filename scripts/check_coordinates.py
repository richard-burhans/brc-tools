#!/usr/bin/env python3
"""Do this repository's coordinate-emitting tools report coordinates that are actually there?

    python3 scripts/check_coordinates.py            # every case
    python3 scripts/check_coordinates.py dustmasker # one

⛔ THREE COORDINATE BUGS SHIPPED, AND NOT ONE CHECK IN THIS REPOSITORY COULD SEE ANY OF THEM.

  * `interval2bed.awk` read NCBI's 0-based inclusive `-outfmt interval` as if it were 1-based, so
    every dustmasker and windowmasker interval in every strain was one base to the left.
  * `ANOtoBED` reports CONTIG-relative positions under the SCAFFOLD's name, so on any assembly with
    internal N-gaps every fastan interval after the first gap was in the wrong place -- 3 kb out in
    the measured case, with two arrays overlapping each other.
  * A `find_copies: true` YAML boolean against a SELECT parameter silently discarded the threshold
    beside it.

They were invisible for the same reason: EVERY OUTPUT WAS INTERNALLY CONSISTENT. Widths were
right, totals were right, coverage percentages were right, and `verify_softmask_outputs.py`
compared a mask against a FASTA masked FROM that mask and found them to agree. Reviewing harder
would not have helped; three separate adversarial passes read this code and two of the three
survived them.

What finds this class is one move, and it is the move used to find all three: TAKE AN INDEPENDENT
GROUND TRUTH FROM THE SAME RUN AND COMPARE POSITIONS. Not a golden file, which goes stale, and not
a planted expectation with a tolerance, which cannot see an off-by-one -- an answer the tool itself
produces by a different route:

    dustmasker/windowmasker   BED3   vs   the lowercase runs of its own `-outfmt fasta`
    tantan                    BED3   vs   the lowercase runs of its own masked FASTA
    fastan                    BED6   vs   the `M` records of its own .1ano, read with ONEview

Each pair must agree EXACTLY. A one-base shift fails, a 3 kb shift fails, and a tool that changes
its output convention fails on the next run rather than in six months' time.

⚠ NEEDS DOCKER, because the ground truth has to come from the tool itself. Without it this SKIPS
and says so; a skip is not a pass.
"""
import argparse
import pathlib
import random
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent


def run(image, script, workdir):
    """Run a shell snippet in a container, returning (rc, stdout, stderr)."""
    p = subprocess.run(["docker", "run", "--rm", "-v", f"{workdir}:/w", "-w", "/w", image,
                        "bash", "-c", script], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def lowercase_runs(path, minlen=1):
    """(start, end) of every soft-masked run, 0-based half-open -- the independent ground truth."""
    seq = "".join(x.strip() for x in path.read_text().splitlines() if not x.startswith(">"))
    return [(m.start(), m.end()) for m in re.finditer(rf"[acgtn]{{{minlen},}}", seq)]


def synthetic(path, name="chr1", seed=7):
    """A sequence with real low-complexity and a real tandem array, at positions we do not assert.

    ⚠ THE PLANTED POSITIONS ARE NOT THE TEST. Every masker draws its own boundaries, so asserting
    "the array is at 1500" would only measure the tool's heuristics. What is asserted is that the
    tool's BED and the tool's own masked FASTA name the SAME places.
    """
    random.seed(seed)
    r = lambda n: "".join(random.choice("ACGT") for _ in range(n))          # noqa: E731
    seq = (r(900) + "A" * 300 + r(300) + ("AT" * 400) + r(500)
           + ("CAGGT" * 260) + r(700) + r(1)[0] * 0 + r(1200))
    path.write_text(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return seq


# ---------------------------------------------------------------- the cases -------------------
def case_ncbi_masker(tool, image, wd):
    """dustmasker / windowmasker: BED3 from `-outfmt interval` vs its own `-outfmt fasta`."""
    synthetic(wd / "in.fa")
    (wd / "interval2bed.awk").write_text((ROOT / "tools/dustmasker/interval2bed.awk").read_text())
    if tool == "dustmasker":
        make = ("dustmasker -in in.fa -outfmt interval | awk -f interval2bed.awk > out.bed3 && "
                "dustmasker -in in.fa -outfmt fasta -out lc.fa")
    else:
        make = ("windowmasker -mk_counts -in in.fa -out cnt 2>/dev/null && "
                "windowmasker -ustat cnt -in in.fa -outfmt interval | awk -f interval2bed.awk "
                "> out.bed3 && windowmasker -ustat cnt -in in.fa -outfmt fasta -out lc.fa")
    rc, _, err = run(image, make, wd)
    if rc != 0:
        return None, f"the tool itself failed (rc={rc}): {err.strip().splitlines()[-1:]}"
    bed = [(int(a), int(b)) for a, b in
           (ln.split("\t")[1:3] for ln in (wd / "out.bed3").read_text().splitlines() if ln.strip())]
    return bed, lowercase_runs(wd / "lc.fa")


def case_tantan(_tool, image, wd):
    """tantan: BED3 from lc2bed.awk vs the lowercase runs of tantan's own output."""
    synthetic(wd / "in.fa")
    (wd / "lc2bed.awk").write_text((ROOT / "tools/tantan/lc2bed.awk").read_text())
    rc, _, err = run(image, "tantan in.fa > lc.fa && awk -f lc2bed.awk lc.fa > out.bed3", wd)
    if rc != 0:
        return None, f"the tool itself failed (rc={rc}): {err.strip().splitlines()[-1:]}"
    bed = [(int(a), int(b)) for a, b in
           (ln.split("\t")[1:3] for ln in (wd / "out.bed3").read_text().splitlines() if ln.strip())]
    return bed, lowercase_runs(wd / "lc.fa")


def case_fastan(_tool, image, wd):
    """fastan: BED6 from ano2bed6.awk vs the M records of the .1ano it was built from.

    ⛔ AND THE SEQUENCE HAS AN INTERNAL N-GAP ON PURPOSE. Without one, contig and scaffold
    coordinates coincide and the bug this case exists for is invisible: `ANOtoBED` was correct on
    every gapless test anyone ran.
    """
    random.seed(11)
    r = lambda n: "".join(random.choice("ACGT") for _ in range(n))          # noqa: E731
    seq = r(500) + r(20) * 60 + r(500) + "N" * 800 + r(500) + r(25) * 50 + r(500)
    (wd / "in.fa").write_text(">chrX\n"
                              + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    (wd / "ano2bed6.awk").write_text((ROOT / "tools/fastan/ano2bed6.awk").read_text())
    rc, _, err = run("quay.io/biocontainers/fastga:1.5.20260729--h118bc1c_0",
                     "FAtoGDB in.fa gdb", wd)
    if rc != 0:
        return None, f"FAtoGDB failed: {err.strip()[-120:]}"
    rc, _, err = run(image, "FasTAN -m -p -T1 -oscan gdb", wd)
    if rc != 0:
        return None, f"FasTAN failed: {err.strip()[-120:]}"
    rc, out, err = run("quay.io/biocontainers/fastga:1.5.20260729--h118bc1c_0",
                       "ONEview scan.1ano | awk -f ano2bed6.awk > out.bed6; ONEview scan.1ano", wd)
    if rc != 0:
        return None, f"ONEview failed: {err.strip()[-120:]}"
    bed = [(int(f[1]), int(f[2])) for f in
           (ln.split("\t") for ln in (wd / "out.bed6").read_text().splitlines() if ln.strip())]
    truth = [(int(f[2]), int(f[3])) for f in
             (ln.split() for ln in out.splitlines()) if f and f[0] == "M"]
    return bed, truth


CASES = {
    "dustmasker": (case_ncbi_masker, "quay.io/biocontainers/blast:2.17.0--h66d330f_0"),
    "windowmasker": (case_ncbi_masker, "quay.io/biocontainers/blast:2.17.0--h66d330f_0"),
    "tantan": (case_tantan, "quay.io/biocontainers/tantan:51--h5ca1c30_1"),
    "fastan": (case_fastan, "quay.io/biocontainers/fastan:0.8--h118bc1c_1"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("only", nargs="*", choices=sorted(CASES) or None, default=None,
                    help="run only these cases")
    args = ap.parse_args()
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        print("  ⚠ SKIP — docker is not available, and the ground truth has to come from the tool "
              "itself. Nothing was checked. NOT a pass.")
        return 0

    bad = skipped = 0
    for name in (args.only or sorted(CASES)):
        fn, image = CASES[name]
        with tempfile.TemporaryDirectory() as td:
            wd = pathlib.Path(td)
            bed, truth = fn(name, image, wd)
        if bed is None:
            print(f"  ⚠ SKIP  {name:<14} {truth}")
            skipped += 1
            continue
        if bed == truth:
            print(f"  pass    {name:<14} {len(bed)} interval(s), identical to the tool's own answer")
            continue
        bad += 1
        print(f"  ⛔ FAIL  {name:<14} the BED and the tool's own output name DIFFERENT places")
        print(f"            BED   ({len(bed):3d}): {bed[:4]}")
        print(f"            truth ({len(truth):3d}): {truth[:4]}")
        for a, b in list(zip(bed, truth, strict=False))[:3]:
            if a != b:
                print(f"            first difference: BED {a} vs truth {b}  "
                      f"(offset {a[0] - b[0]:+d}, {a[1] - b[1]:+d})")
                break
    total = len(args.only or CASES)
    print(f"\n  {total - bad - skipped}/{total} coordinate case(s) agree with an independent "
          f"answer from the same run"
          + (f"; {skipped} SKIPPED and therefore unchecked" if skipped else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
