# MT — a vertebrate mitochondrial dataset for exercising the workflows

A development-scale stand-in for the Pv4 panel. Eight ~25 Mb *P. vivax* genomes with
a 56-cell alignment grid take hours per pass, which is a poor loop for polishing a
workflow. Vertebrate mitochondria are the same problem two orders of magnitude
smaller — **22 species, ~16.5 kb each, 13 protein-coding genes** — so a full pass
finishes in minutes.

The point is that **the right answer is known in advance**: mitochondrial gene
content is identical across vertebrates, so a correct run should produce exactly
**13 orthogroups, each holding one gene per species, all labelled CORE-1:1**.
Anything else is a bug in the pipeline rather than a fact about the data.

## What is here

```
MT/
├── fasta/<species>.fa      one record, the chrM sequence
├── gff3/<species>.gff3     gene / mRNA / exon / CDS for the 13 protein-coding genes
├── manifest.tsv            species, assembly, chrom, length, gene count, source track
└── ground_truth.tsv        every gene with its symbol — the expected orthogroups
```

Total: 560 kB. Everything came from the UCSC REST API; `scripts/fetch_mt_dataset.py`
regenerates it.

22 species spanning ~430 My of divergence: six primates, mouse, rat, rabbit, guinea
pig, dog, cat, horse, cow, sheep, pig, opossum, platypus, chicken, turkey, Xenopus,
zebrafish.

## Suggested use

Human and mouse as anchors, the rest as queries — mirroring the Pv4 layout of three
anchors against a panel. That gives a 2 × 21 projection grid and a 22-strain
all-against-all alignment, both small enough to iterate on.

`ground_truth.tsv` maps every gene id to its symbol, so validating a run is a join:
group by symbol and check each group has 22 members and one orthogroup id.

## Three things to know about the data

**Gene ids are prefixed with the assembly** (`hg38_ND1`, `mm39_ND1`) so no two
species share one. Without that the orthology would be readable straight off the
names and the pipeline would never actually be tested. Transcript ids carry a `_t1`
suffix, matching the PlasmoDB convention the existing tools already normalise away.

**Only protein-coding genes are included.** The filter is `cdsStart != cdsEnd`,
which drops tRNAs and rRNAs. It is applied because some UCSC tracks include them and
some do not — without it the gene count would swing between 13 and 37 depending on
which track a species happens to have.

**Many CDS lengths are not divisible by three, and that is correct.** Vertebrate
mitochondrial genes frequently end in an incomplete stop codon (`T` or `TA`)
completed by polyadenylation. Across the 286 genes here, 64% are divisible by three,
31% are one base short and 5% are two short — human `ND2`, `COX3`, `ND3`, `ND4` and
`CYTB` among them. Anything translating these needs to truncate to a codon boundary;
`group_cds_by_og` already does.

## What this dataset cannot test

**Mitochondrial genes have no introns.** 285 of the 286 genes here are single-exon.
The one exception is turkey `ND3`, whose two "exons" are separated by a 1 bp gap —
the avian ND3 +1 frameshift, a single inserted nucleotide that is translationally
bypassed, not an intron.

Several things follow, and they bound what a green run here proves:

- **SpliceAI is irrelevant.** It predicts splice sites; there are none. Running the
  projection with or without it should make no difference, so this dataset cannot
  serve as a SpliceAI evaluation.
- **Splice-site and exon-boundary logic is never exercised.** Multi-exon CESAR
  placement, intron-aware alignment, and the intactness grades that depend on exon
  structure all go untested.
- **It cannot reproduce the `PVPAM_000040100_t1` class of bug.** That failure is a
  terminal 18 bp exon in a two-exon gene never receiving a CESAR placement. Nothing
  here has a terminal exon in that sense.

So this validates the *mechanics* — collection shapes, the alignment grid, chain
building, orthogroup assembly, the wiring between workflows — which is what most of
the recent breakage has been. It does not validate annotation transfer as biology.
Pv4 remains the only test for that, and *Plasmodium* genes do have introns, which is
also why Pv4 rather than this is the right place to evaluate SpliceAI.

## One caveat on the ground truth

Three assemblies — `galGal6` (chicken), `monDom5` (opossum), `ornAna2` (platypus) —
have only Ensembl gene models on chrM, which carry no gene symbols. Their symbols are
**inferred from gene order** and are marked `inferred` in both `manifest.tsv` and
`ground_truth.tsv`.

That inference is sound but worth stating: the canonical vertebrate order
(ND1, ND2, COX1, COX2, ATP8, ATP6, COX3, ND3, ND4L, ND4, ND5, ND6, CYTB) was checked
against human, mouse, zebrafish and Xenopus — mammal, fish and amphibian — and held
exactly in all four before being applied to a bird, a marsupial and a monotreme.
A strict validation can exclude those three species and still have 19.
