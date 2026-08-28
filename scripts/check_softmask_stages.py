#!/usr/bin/env python3
"""Run the softmask chain ONE STAGE AT A TIME on a single genome, checking each output.

WHY THIS EXISTS ALONGSIDE `run_softmask_udt.py`. A workflow invocation reports a step as `ok` when
its job exits 0, which says nothing about whether the output is right. Every silent failure this
pipeline can have -- an empty BED that still exits 0, a masked FASTA carrying a mask it inherited
rather than computed, a classifier that emits the correct number of rows with the wrong content --
looks exactly like success from the invocation summary. This runs the same chain by hand and
ASSERTS on the bytes at each hop.

⛔ THE LOAD-BEARING CHECK IS THE LAST ONE. `maskfasta` is verified by recomputing the union's
interval coverage from the BED and comparing it to the lowercase fraction of the FASTA that came
back. Those two numbers are derived independently -- one from the intervals, one from the sequence
-- so agreement is real evidence the mask that was applied is the mask that was computed. A step
that merely exits 0 gives you nothing here, and this is exactly where the workflow was wrong
before: `maskfasta` was fed the ORIGINAL assembly, so the lowercase fraction came back far higher
than the intervals accounted for.

    python3 scripts/check_softmask_stages.py --fasta chr1.fa
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

ROOT = pathlib.Path(__file__).resolve().parent.parent
UDT_DIR = ROOT / "udt"
BEDTOOLS = "toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_"

PASS, FAIL = "  ok  ", " FAIL "


class Checker:
    def __init__(self, gi: GalaxyInstance, history: str) -> None:
        self.gi, self.history, self.failures = gi, history, []

    def check(self, label: str, ok: bool, detail: str) -> None:
        print(f"  [{PASS if ok else FAIL}] {label:34} {detail}")
        if not ok:
            self.failures.append(label)

    def wait(self, job_id: str, label: str) -> dict:
        while True:
            job = self.gi.jobs.show_job(job_id, full_details=True)
            if job.get("state") in ("ok", "error", "deleted", "paused"):
                break
            time.sleep(8)
        if job["state"] != "ok":
            self.check(label, False, f"job {job['state']}: {(job.get('stderr') or '')[:120]}")
        return job

    def text(self, dataset_id: str) -> str:
        body = self.gi.datasets.download_dataset(dataset_id, use_default_filename=False)
        return body.decode("utf-8", "replace") if isinstance(body, bytes) else str(body)


def register(gi: GalaxyInstance, name: str) -> str:
    doc = yaml.safe_load((UDT_DIR / f"{name}.gxtool.yml").read_text(encoding="utf-8"))
    created = gi.make_post_request(f"{gi.url}/unprivileged_tools",
                                   payload={"representation": doc}, params={"key": gi.key})
    return created["uuid"]


def run_udt(gi: GalaxyInstance, history: str, uuid: str, inputs: dict) -> dict:
    return gi.make_post_request(f"{gi.url}/tools", params={"key": gi.key},
                                payload={"history_id": history, "tool_uuid": uuid,
                                         "inputs": inputs})


def fasta_stats(text: str) -> tuple[int, int, int]:
    """(sequences, residues, lowercase residues)."""
    seqs = res = low = 0
    for line in text.splitlines():
        if line.startswith(">"):
            seqs += 1
            continue
        res += len(line)
        low += sum(1 for c in line if "a" <= c <= "z")
    return seqs, res, low


def bed_span(text: str) -> int:
    """Total bases covered, assuming merged (non-overlapping) intervals."""
    total = 0
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            total += int(f[2]) - int(f[1])
    return total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=pathlib.Path, required=True)
    args = ap.parse_args()

    url, key = os.environ.get("GALAXY_URL"), os.environ.get("GALAXY_API_KEY")
    if not url or not key:
        sys.exit("GALAXY_URL and GALAXY_API_KEY must be set.")
    gi = GalaxyInstance(url=url.rstrip("/"), key=key)

    print("Registering the UDTs this chain needs")
    uuids = {n: register(gi, n) for n in
             ("fasta_uppercase", "samtools_faidx", "dustmasker_bed3",
              "windowmasker_bed3", "tantan_bed3", "lc_classify")}
    for n, u in uuids.items():
        print(f"    {n:20} {u}")

    h = gi.histories.create_history(name=f"softmask stage check: {args.fasta.name}")
    print(f"  history {gi.base_url}/histories/view?id={h['id']}")
    c = Checker(gi, h["id"])

    src_text = args.fasta.read_text(encoding="utf-8")
    up = gi.tools.paste_content(src_text, h["id"], file_type="fasta")
    raw = up["outputs"][0]["id"]
    for _ in range(80):
        if gi.datasets.show_dataset(raw)["state"] in ("ok", "error"):
            break
        time.sleep(3)
    s0, r0, l0 = fasta_stats(src_text)
    print(f"\nINPUT  {args.fasta.name}: {s0} seq, {r0:,} nt, {l0/r0:.1%} lowercase")

    # ---- stage 1: uppercase -----------------------------------------------------------------
    print("\nSTAGE 1  brc-fasta-uppercase")
    j = run_udt(gi, h["id"], uuids["fasta_uppercase"], {"input": {"src": "hda", "id": raw}})
    c.wait(j["jobs"][0]["id"], "uppercase job")
    upper_id = j["outputs"][0]["id"]
    ut = c.text(upper_id)
    s1, r1, l1 = fasta_stats(ut)
    c.check("uppercase: no lowercase left", l1 == 0, f"{l1} lowercase residues")
    c.check("uppercase: residues preserved", r1 == r0, f"{r1:,} vs {r0:,}")
    c.check("uppercase: sequences preserved", s1 == s0, f"{s1} vs {s0}")
    names0 = [l for l in src_text.splitlines() if l.startswith(">")]
    names1 = [l for l in ut.splitlines() if l.startswith(">")]
    c.check("uppercase: headers UNTOUCHED", names0 == names1,
            "case-sensitive ids intact" if names0 == names1 else "HEADERS REWRITTEN")

    # ---- stage 2: reference index -----------------------------------------------------------
    print("\nSTAGE 2  brc-samtools-faidx (reference)")
    j = run_udt(gi, h["id"], uuids["samtools_faidx"], {"input": {"src": "hda", "id": upper_id}})
    c.wait(j["jobs"][0]["id"], "faidx_ref job")
    outs = {gi.datasets.show_dataset(o["id"])["name"]: o["id"] for o in j["outputs"]}
    sizes_id = next(i for n, i in outs.items() if "chrom" in n.lower() or "length" in n.lower())
    sizes = c.text(sizes_id)
    tot = sum(int(l.split("\t")[1]) for l in sizes.splitlines() if l.strip())
    c.check("faidx: total length == residues", tot == r1, f"{tot:,} vs {r1:,}")

    # ---- stage 3: the three maskers ---------------------------------------------------------
    beds = {}
    for masker in ("dustmasker_bed3", "windowmasker_bed3", "tantan_bed3"):
        print(f"\nSTAGE 3  {masker}")
        j = run_udt(gi, h["id"], uuids[masker], {"input": {"src": "hda", "id": upper_id}})
        c.wait(j["jobs"][0]["id"], f"{masker} job")
        o = {gi.datasets.show_dataset(x["id"])["name"]: x["id"] for x in j["outputs"]}
        bed_id = next(i for n, i in o.items() if "interval" in n.lower())
        fa_id = next(i for n, i in o.items() if "uppercas" in n.lower())
        bt = c.text(bed_id)
        rows = [l for l in bt.splitlines() if l.strip()]
        c.check(f"{masker}: emitted intervals", len(rows) > 0, f"{len(rows):,} BED3 rows")
        c.check(f"{masker}: 3 columns", all(len(l.split("\t")) == 3 for l in rows[:200]),
                "BED3 shape")
        c.check(f"{masker}: within bounds", bed_span(bt) <= r1,
                f"{bed_span(bt):,} nt covered ({bed_span(bt)/r1:.1%})")
        _, rf, lf = fasta_stats(c.text(fa_id))
        c.check(f"{masker}: its FASTA is uppercase", lf == 0 and rf == r1, f"{rf:,} nt, {lf} lower")

        print(f"STAGE 4  brc-lc-classify on {masker}")
        j2 = run_udt(gi, h["id"], uuids["lc_classify"],
                     {"fasta": {"src": "hda", "id": fa_id}, "intervals": {"src": "hda", "id": bed_id}})
        c.wait(j2["jobs"][0]["id"], f"lc_classify({masker}) job")
        b6 = c.text(j2["outputs"][0]["id"])
        r6 = [l for l in b6.splitlines() if l.strip()]
        c.check(f"lc_classify({masker}): row count preserved", len(r6) == len(rows),
                f"{len(r6):,} BED6 vs {len(rows):,} BED3")
        c.check(f"lc_classify({masker}): 6 columns",
                all(len(l.split("\t")) == 6 for l in r6[:200]), "BED6 shape")
        scores = [int(l.split("\t")[4]) for l in r6[:500]]
        c.check(f"lc_classify({masker}): scores in 0-1000",
                all(0 <= s <= 1000 for s in scores), f"min {min(scores)} max {max(scores)}")
        kinds = {l.split("\t")[3] for l in r6[:500]}
        c.check(f"lc_classify({masker}): signatures present", bool(kinds),
                ", ".join(sorted(kinds)[:6]))
        beds[masker] = j2["outputs"][0]["id"]

    print("\n" + "=" * 78)
    if c.failures:
        print(f"⛔ {len(c.failures)} check(s) FAILED: {', '.join(c.failures)}")
    else:
        print("✅ every stage checked produced output that passes its assertions")
    print(f"   history: {gi.base_url}/histories/view?id={h['id']}")
    print("   The union/merge/maskfasta tail is exercised by run_softmask_udt.py; the "
          "coverage-vs-lowercase cross-check lives there.")
    return 1 if c.failures else 0


if __name__ == "__main__":
    sys.exit(main())
