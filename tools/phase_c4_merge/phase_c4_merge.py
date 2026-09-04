#!/usr/bin/env python3
"""
Phase C.4 — merge Liftoff-clean + TOGA2/CESAR2 outputs into unified annotation.

Rules:
  - Genes in liftoff_clean.gff3      → source=liftoff, intactness=I
  - Genes in TOGA2 query_annotation  → source=cesar2, intactness from loss_summary.tsv
  - Reference genes absent from both  → source=none, intactness=M (missing)

Outputs per query:
  {Q}.annotation.gff3        — unified GFF (Liftoff entries + CESAR2 BED-converted GFF)
  {Q}.classification.tsv     — flat per-gene table:
        reference_gene_id, query_gene_id, source, intactness,
        query_chrom, query_start, query_end, query_strand, orthology_class

Ported from /media/anton/data/sandbox/Pv4/v3/scripts/phase_c4_merge.py

TOGA2 v2.0.x output schema (verified against the toga2:local image source,
/opt/TOGA2/src/python/modules/constants.py + filter_loss_file.py +
finalise_orthology_files.py). All four files land directly in the run's
--output dir (toga_out/). This script consumes:

  loss_summary.tsv             header `level<TAB>entry<TAB>status` (3 cols).
                               level in {PROJECTION, TRANSCRIPT, GENE};
                               PROJECTION `entry` = `<transcript>#<chain>`
                               (may carry `#paralog`/`#retro` suffix);
                               status in the loss alphabet
                               {FI,I,PI,UL,M,L,PG,PP,N}.
  orthology_classification.tsv header
                               `t_gene<TAB>t_transcript<TAB>q_gene<TAB>q_transcript<TAB>orthology_class`;
                               orthology_class in
                               {one2one,one2many,many2one,many2many,one2zero}.
                               q_transcript = projection id `<transcript>#<chain>`.
  query_annotation.bed         BED12; col4 = projection id `<transcript>#<chain>`.
  query_genes.bed              BED (>=4 col); col4 = final query gene id (== q_gene).

The TOGA1 ports of these loaders already match the TOGA2 column layout
(level/entry/status; t_gene/.../orthology_class; projection `#` separator), so
no field renames were needed on the TOGA1->TOGA2 modernization. Gene-level
query coords come from query_genes.bed (keyed on q_gene); intactness is matched
from loss_summary via the `<t_transcript>#` projection-id prefix.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def parse_gff_attributes(s):
    d = {}
    for kv in s.strip().rstrip(';').split(';'):
        kv = kv.strip()
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d


def load_liftoff_clean(gff_path):
    """Return dict reference_gene_id -> list of GFF lines."""
    genes: dict = {}
    current_ref_id = None
    GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}
    if not Path(gff_path).exists():
        return genes
    with open(gff_path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip():
                continue
            fields = ln.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype in GENE_TYPES:
                attrs = parse_gff_attributes(fields[8])
                gid = attrs.get('ID', '')
                # ⛔ ASK LIFTOFF, DO NOT GUESS FROM THE NAME -- see the long note on
                # normalize_gene_id in phase_c2_triage.py, which this is the second copy of.
                # Stripping any `_<1-2 digits>` suffix on sight cannot tell an extra copy from a
                # reference gene named `SICAvar_1`, and it made this function merge two real genes
                # under one invented key: three genes in, four classification rows out, two of them
                # claiming genes were never projected that the GFF3 beside them reports as intact.
                # `extra_copy_number` is written by Liftoff on every gene it projects and survives
                # into liftoff_clean.gff3, because triage streams the original lines through.
                ref_id = gid
                try:
                    copies = int(str(attrs.get('extra_copy_number', '0')).strip() or '0')
                except ValueError:
                    copies = 0
                if copies > 0 and '_' in gid:
                    parts = gid.rsplit('_', 1)
                    if len(parts[1]) <= 2 and parts[1].isdigit() and not parts[0].endswith('_'):
                        ref_id = parts[0]
                current_ref_id = ref_id
                genes.setdefault(ref_id, []).append(ln)
            elif current_ref_id is not None:
                genes[current_ref_id].append(ln)
    return genes


def load_toga2_loss_summary(path):
    if not Path(path).exists():
        return {}
    status = {}
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip('\n')
            if not ln or ln.startswith('level'):
                continue
            fields = ln.split('\t')
            if len(fields) >= 3 and fields[0] == 'PROJECTION':
                status[fields[1]] = fields[2]
    return status


def load_toga2_orthology(path):
    if not Path(path).exists():
        return {}
    out: dict = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            t_gene = row.get('t_gene')
            if not t_gene:
                continue
            out.setdefault(t_gene, []).append({
                'q_gene': row.get('q_gene'),
                'q_tx': row.get('q_transcript'),
                'class': row.get('orthology_class'),
                't_tx': row.get('t_transcript'),
            })
    return out


def load_toga2_query_bed(path):
    if not Path(path).exists():
        return {}
    out = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith('#'):
                continue
            fields = ln.rstrip('\n').split('\t')
            if len(fields) >= 4:
                out[fields[3]] = ln
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Phase C.4: merge Liftoff-clean + TOGA2 outputs into unified annotation")
    ap.add_argument('--query', required=True, help='Query strain name (e.g. Pk_ANKA)')
    ap.add_argument('--triage-dir', required=True, help='work/02b_triage/{anchor}-as-ref/{Q}/')
    ap.add_argument('--toga-dir', required=True, help='work/02c_toga/{anchor}-as-ref/{Q}/')
    ap.add_argument('--out-dir', required=True, help='work/02d_merged/{anchor}-as-ref/')
    ap.add_argument('--ref-bed', required=True, help='inputs/annotations/{anchor}.bed12')
    args = ap.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Load Liftoff-clean GFF
    liftoff_clean = load_liftoff_clean(f"{args.triage_dir}/liftoff_clean.gff3")

    # 2) Load TOGA2 outputs
    loss = load_toga2_loss_summary(f"{args.toga_dir}/loss_summary.tsv")
    ortho = load_toga2_orthology(f"{args.toga_dir}/orthology_classification.tsv")
    query_bed = load_toga2_query_bed(f"{args.toga_dir}/query_annotation.bed")
    query_genes_bed = load_toga2_query_bed(f"{args.toga_dir}/query_genes.bed")

    # 3) Reference gene set from BED
    ref_genes: set = set()
    with open(args.ref_bed) as f:
        for ln in f:
            fields = ln.rstrip('\n').split('\t')
            if len(fields) >= 4:
                ref_genes.add(fields[3])

    # 4) Build classification rows
    rows = []
    seen_ref: set = set()

    # 4a) Liftoff clean
    for ref_id, lines in liftoff_clean.items():
        seen_ref.add(ref_id)
        gene_line = lines[0]
        fields = gene_line.rstrip('\n').split('\t')
        if len(fields) < 9:
            continue
        chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
        attrs = parse_gff_attributes(fields[8])
        rows.append({
            'reference_gene_id': ref_id,
            'query_gene_id': attrs.get('ID', ''),
            'source': 'liftoff',
            'intactness': 'I',
            'query_chrom': chrom,
            'query_start': start,
            'query_end': end,
            'query_strand': strand,
            'orthology_class': 'liftoff_clean',
        })

    # 4b) TOGA2 / CESAR2
    for ref_gene_id, projections in ortho.items():
        for p in projections:
            q_gene = p['q_gene']
            if q_gene in (None, '', 'None'):
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': '',
                    'source': 'cesar2',
                    'intactness': 'L',
                    'query_chrom': '',
                    'query_start': '',
                    'query_end': '',
                    'query_strand': '',
                    'orthology_class': p.get('class', 'one2zero'),
                })
                seen_ref.add(ref_gene_id)
                continue
            tx = p.get('t_tx', '')
            status_keys = [k for k in loss if k.startswith(tx + '#')]
            status = loss[status_keys[0]] if status_keys else '?'
            bed_line = query_bed.get(q_gene, '') or query_genes_bed.get(q_gene, '')
            if bed_line:
                f2 = bed_line.rstrip('\n').split('\t')
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': q_gene,
                    'source': 'cesar2',
                    'intactness': status,
                    'query_chrom': f2[0] if len(f2) > 0 else '',
                    'query_start': f2[1] if len(f2) > 1 else '',
                    'query_end': f2[2] if len(f2) > 2 else '',
                    'query_strand': f2[5] if len(f2) > 5 else '',
                    'orthology_class': p.get('class', ''),
                })
            else:
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': q_gene,
                    'source': 'cesar2',
                    'intactness': status,
                    'query_chrom': '',
                    'query_start': '',
                    'query_end': '',
                    'query_strand': '',
                    'orthology_class': p.get('class', ''),
                })
            seen_ref.add(ref_gene_id)

    # 4c) Missing (neither Liftoff nor TOGA2)
    for ref_id in ref_genes - seen_ref:
        rows.append({
            'reference_gene_id': ref_id,
            'query_gene_id': '',
            'source': 'none',
            'intactness': 'M',
            'query_chrom': '',
            'query_start': '',
            'query_end': '',
            'query_strand': '',
            'orthology_class': 'unprojected',
        })

    # 5) Write classification.tsv
    out_cls = outdir / f"{args.query}.classification.tsv"
    fieldnames = ['reference_gene_id', 'query_gene_id', 'source', 'intactness',
                  'query_chrom', 'query_start', 'query_end', 'query_strand', 'orthology_class']
    with open(out_cls, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # 6) Write merged GFF
    out_gff = outdir / f"{args.query}.annotation.gff3"
    n_lo = n_cs = 0

    # Build lookup dicts for TOGA2 gene decoration
    q_to_ref: dict = {}
    for rg, projs in ortho.items():
        for p in projs:
            qg = p.get('q_gene')
            if qg and qg != 'None':
                q_to_ref[qg] = (rg, p.get('class', ''), p.get('t_tx', ''))

    tx_to_status: dict = {}
    for k, v in loss.items():
        if '#' in k:
            tx_prefix = k.rsplit('#', 1)[0]
            if tx_prefix not in tx_to_status:
                tx_to_status[tx_prefix] = v

    GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}

    with open(out_gff, 'w') as f:
        f.write('##gff-version 3\n')
        f.write(f'# Phase C.4 merged annotation for {args.query}\n')
        f.write('# source=liftoff: triage-clean Liftoff projection (intactness=I)\n')
        f.write('# source=cesar2: TOGA2/CESAR2 fallback; intactness from TOGA2 loss_summary\n')

        for ref_id, lines in liftoff_clean.items():
            for ln in lines:
                fields = ln.rstrip('\n').split('\t')
                if len(fields) >= 9:
                    fields[8] = fields[8].rstrip(';') + ';source=liftoff;intactness=I'
                    f.write('\t'.join(fields) + '\n')
                    if fields[2] in GENE_TYPES:
                        n_lo += 1

        for q_gene, bed_ln in query_genes_bed.items():
            fields = bed_ln.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            chrom, start, end, name, score, strand = fields[:6]
            ref_gene_id, orth_class, t_tx = q_to_ref.get(q_gene, ('', '', ''))
            intactness = tx_to_status.get(t_tx, '?')
            attrs = (f"ID={name};reference_gene_id={ref_gene_id}"
                     f";source=cesar2;intactness={intactness}"
                     f";orthology_class={orth_class}")
            f.write(f"{chrom}\tTOGA2\tprotein_coding_gene\t{int(start)+1}\t{end}"
                    f"\t.\t{strand}\t.\t{attrs}\n")
            n_cs += 1

    src_counts = Counter(r['source'] for r in rows)
    intact_counts = Counter(r['intactness'] for r in rows)
    print(f"[{args.query}] {len(rows)} rows; sources={dict(src_counts)} "
          f"intactness={dict(intact_counts)}")
    print(f"  GFF: {n_lo} liftoff + {n_cs} CESAR2 genes → {out_gff}")
    print(f"  Classification: {out_cls}")


if __name__ == '__main__':
    main()
