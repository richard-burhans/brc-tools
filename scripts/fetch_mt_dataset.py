#!/usr/bin/env python3
"""Build a small vertebrate mitochondrial dataset for exercising the workflows.

The Pv4 panel is the real target but a poor development loop: eight ~25 Mb
genomes, a 56-cell alignment grid and runs measured in hours. Vertebrate
mitochondria are the same problem two orders of magnitude smaller -- ~16.5 kb,
13 protein-coding genes, gene content identical across the whole clade -- so a
full A-to-F pass finishes in minutes and the right answer is known in advance:
13 orthogroups, each holding exactly one gene per species.

Everything comes from the UCSC REST API (api.genome.ucsc.edu), one assembly at a
time: the chrM sequence, and the gene models on chrM from the best annotation
track that assembly has.

Only protein-coding genes are kept. The filter is cdsStart != cdsEnd, which
drops the tRNAs and rRNAs that some tracks include and some do not -- without it
the gene count would swing between 13 and 37 depending on which track a species
happens to have, which is not a difference you want in a validation set.

Gene ids are prefixed with the assembly (hg38_ND1, mm39_ND1) so that no two
species share one. That matters: with bare symbols the orthology would be
readable straight off the names and the pipeline would never be tested. The
symbol is kept separately in ground_truth.tsv. Transcript ids get a _t1 suffix,
matching the PlasmoDB convention the existing tools already normalise away.

Three assemblies (galGal6, monDom5, ornAna2) have only Ensembl gene models on
chrM, which carry no symbols. Their symbols are INFERRED from gene order, which
is safe enough to do and worth flagging: the canonical vertebrate order was
confirmed to hold exactly in human, mouse, zebrafish and Xenopus -- mammal, fish
and amphibian -- before it was applied to a bird, a marsupial and a monotreme.
ground_truth.tsv marks those rows so a validation can exclude them.

Usage:  python scripts/fetch_mt_dataset.py [--out MT]
"""
import argparse
import json
import os
import subprocess
import sys
import time

API = "https://api.genome.ucsc.edu"

# canonical vertebrate mtDNA protein-coding gene order, 5' to 3' on the heavy strand
CANONICAL = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
             "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"]

# (label, assembly). Latest UCSC assembly per species that has a chrM with genes.
SPECIES = [
    ("human",      "hg38"),
    ("chimp",      "panTro6"),
    ("gorilla",    "gorGor6"),
    ("orangutan",  "ponAbe3"),
    ("rhesus",     "rheMac10"),
    ("marmoset",   "calJac4"),
    ("mouse",      "mm39"),
    ("rat",        "rn7"),
    ("rabbit",     "oryCun2"),
    ("guineapig",  "cavPor3"),
    ("dog",        "canFam6"),
    ("cat",        "felCat9"),
    ("horse",      "equCab3"),
    ("cow",        "bosTau9"),
    ("sheep",      "oviAri4"),
    ("pig",        "susScr11"),
    ("opossum",    "monDom5"),
    ("platypus",   "ornAna2"),
    ("chicken",    "galGal6"),
    ("turkey",     "melGal5"),
    ("xenopus",    "xenTro10"),
    ("zebrafish",  "danRer11"),
]

TRACKS = ["ncbiRefSeqCurated", "ncbiRefSeq", "refGene", "ensGene"]


def api(path, tries=3):
    for i in range(tries):
        try:
            return json.loads(subprocess.check_output(
                ["curl", "-sS", "--max-time", "60", f"{API}{path}"]))
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def mt_name(db):
    ch = api(f"/list/chromosomes?genome={db}").get("chromosomes") or {}
    for k in ch:
        if k.lower() in ("chrm", "chrmt"):
            return k, ch[k]
    return None, None


def coding_genes(db, chrom):
    """(track, [rows]) for the first track with protein-coding models on chrom."""
    for tr in TRACKS:
        rows = api(f"/getData/track?genome={db};track={tr};chrom={chrom}").get(tr)
        if isinstance(rows, list):
            coding = [r for r in rows if r.get("cdsStart") != r.get("cdsEnd")]
            if coding:
                coding.sort(key=lambda r: r["txStart"])
                return tr, coding
    return None, []


def symbols_for(rows):
    """(symbols, inferred?). Ensembl ids are not symbols, so fall back to order."""
    named = [(r.get("name2") or "").upper() for r in rows]
    if len(rows) == len(CANONICAL) and all(s in CANONICAL for s in named):
        return named, False
    if len(rows) == len(CANONICAL):
        return list(CANONICAL), True
    return [f"GENE{i + 1}" for i in range(len(rows))], True


def write_fasta(path, name, dna, width=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(dna), width):
            fh.write(dna[i:i + width] + "\n")


def write_gff3(path, db, chrom, length, rows, symbols):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {chrom} 1 {length}\n")
        for r, sym in zip(rows, symbols):
            gid = f"{db}_{sym}"
            tid = f"{gid}_t1"
            s, e, strand = r["txStart"] + 1, r["txEnd"], r["strand"]
            fh.write(f"{chrom}\tUCSC\tgene\t{s}\t{e}\t.\t{strand}\t.\t"
                     f"ID={gid};Name={sym}\n")
            fh.write(f"{chrom}\tUCSC\tmRNA\t{s}\t{e}\t.\t{strand}\t.\t"
                     f"ID={tid};Parent={gid};Name={sym}\n")
            starts = [int(x) for x in r["exonStarts"].rstrip(",").split(",")]
            ends = [int(x) for x in r["exonEnds"].rstrip(",").split(",")]
            cds_s, cds_e = r["cdsStart"], r["cdsEnd"]
            for i, (xs, xe) in enumerate(zip(starts, ends), 1):
                fh.write(f"{chrom}\tUCSC\texon\t{xs + 1}\t{xe}\t.\t{strand}\t.\t"
                         f"ID={tid}.exon{i};Parent={tid}\n")
                cs, ce = max(xs, cds_s), min(xe, cds_e)
                if cs < ce:
                    fh.write(f"{chrom}\tUCSC\tCDS\t{cs + 1}\t{ce}\t.\t{strand}\t0\t"
                             f"ID={tid}.cds{i};Parent={tid}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="MT")
    a = ap.parse_args()
    root = os.path.abspath(a.out)
    os.makedirs(os.path.join(root, "fasta"), exist_ok=True)
    os.makedirs(os.path.join(root, "gff3"), exist_ok=True)

    manifest, truth, failed = [], [], []
    for label, db in SPECIES:
        try:
            chrom, length = mt_name(db)
            if not chrom:
                failed.append((label, db, "no chrM")); continue
            track, rows = coding_genes(db, chrom)
            if not rows:
                failed.append((label, db, "no coding genes")); continue
            seq = api(f"/getData/sequence?genome={db};chrom={chrom}").get("dna", "")
            if len(seq) != length:
                failed.append((label, db, f"sequence {len(seq)} != {length}")); continue
            syms, inferred = symbols_for(rows)
            write_fasta(os.path.join(root, "fasta", f"{label}.fa"), chrom, seq)
            write_gff3(os.path.join(root, "gff3", f"{label}.gff3"), db, chrom, length, rows, syms)
            manifest.append((label, db, chrom, length, len(rows), track,
                             "inferred" if inferred else "from-track"))
            for r, sym in zip(rows, syms):
                truth.append((label, f"{db}_{sym}", sym, r["txStart"] + 1, r["txEnd"],
                              r["strand"], "inferred" if inferred else "from-track"))
            print(f'  {label:11} {db:10} {chrom} {length:>6} bp  {len(rows)} genes  '
                  f'{track}{"  (symbols inferred from order)" if inferred else ""}')
        except Exception as e:
            failed.append((label, db, str(e)[:60]))
            print(f'  {label:11} {db:10} FAILED: {e}', file=sys.stderr)

    with open(os.path.join(root, "manifest.tsv"), "w") as fh:
        fh.write("species\tassembly\tchrom\tlength\tn_genes\ttrack\tsymbol_source\n")
        for row in manifest:
            fh.write("\t".join(str(x) for x in row) + "\n")
    with open(os.path.join(root, "ground_truth.tsv"), "w") as fh:
        fh.write("species\tgene_id\tsymbol\tstart\tend\tstrand\tsymbol_source\n")
        for row in truth:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"\n{len(manifest)} species written to {root}")
    if failed:
        print("failed:", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
