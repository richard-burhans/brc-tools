<!--
Prose for workflows/pipeline_io_map.html. Edit this file, then regenerate:

    python workflows/gen_pipeline_io_map.py

One `## <id>` section per workflow, where <id> matches an entry in VETTED in the
generator. Each needs a `### summary` (one line, shown in the Port detail
section) and a `### description` (the plain-language paragraph shown above the
diagram). Blank lines start a new paragraph. Inline **bold**, *italic*, `code`
and [links](https://example.org) work.
-->

# Workflow descriptions


## A

<!-- WF-A inventory -->

### summary

Panel inventory: sourmash sketch/compare for pairwise similarity, BUSCO for per-strain
completeness, and the derived panel bookkeeping (sizes, self-pairs, relabel map) that
WF-C needs — generated in-workflow rather than hand-authored.


### description

Before aligning anything, you want to know what is actually in the panel: how similar
the genomes are to each other, and whether any of them is missing a chunk of its gene
content. This workflow answers both. It takes the genome sequences and their protein
sets. It compares every genome against every other with sourmash, which comes to a
similarity score from shared k-mers without doing an alignment, so it is fast even on
whole genomes. Separately it runs BUSCO, which asks how many of the several hundred
genes expected to be present in exactly one copy in this clade can be found in each
strain. What comes out is an all-against-all similarity matrix with a clustered heatmap
and a tree, a completeness score per strain, and one QC report holding both. It also
writes the small bookkeeping files the alignment step needs later (chromosome lengths,
the list of genome pairs), so nobody has to make those by hand. The point is to catch a
bad assembly or an unexpected outlier here, before spending days of compute on it. On
the Pv4 panel it did exactly that: PvSY56 sits at about 0.24 similarity to everything
else while the rest sit at 0.63 to 0.70, and MHC087 came back only 87% complete.

One more piece of bookkeeping happens here: for each anchor -- the few strains whose gene
annotation you trust well enough to project onto the others -- it reads the curated GFF3
and writes out the two files the annotation-transfer workflow will need, a BED12 of the
protein-coding transcripts and a table saying which gene each transcript belongs to. These
were previously built by hand outside Galaxy, which meant the panel's provenance had a gap
in it right where the reference annotation entered the pipeline.



### step:assemblies

The genomes themselves, one per strain. Everything downstream is keyed by these element
identifiers, so the names chosen here become collection ids, pair ids (`A_B`) and eventually
track names in the genome browser.

### step:proteomes

One protein FASTA per strain, in the same order as the assemblies. BUSCO is run in protein
mode rather than genome mode, which is faster and avoids BUSCO doing its own gene prediction.
Normally these come from `gffread` on each strain's own annotation.

### step:sourmash_sketch

Reduces each genome to a **signature**: a small, fixed-size sample of the k-mers it contains,
chosen by hashing. Two genomes that share a lot of sequence share a lot of hashes, so
signatures can be compared instead of the genomes themselves. A 29 Mb genome collapses to a
few hundred kB, which is what makes the all-against-all comparison cheap.

### step:sourmash_compare

Compares every signature against every other and writes the similarity matrix, plus a
clustered heatmap and a dendrogram. Values are similarity: 1.0 is identical.
This is where an outlier becomes obvious — on the Pv4 panel PvSY56 sits near 0.24 against
everything while the rest sit at 0.63 to 0.70.

### step:busco

Asks how many of the single-copy genes expected across this clade can actually be found in
each strain's proteome. It is a completeness check on the **annotation**: a low score means genes are missing or
unannotated, which matters because a strain that is short on genes will also be short on
orthologues later.

### step:panel_faidx

Indexes each assembly with `samtools faidx`. The index is wanted for its second column,
the length of every sequence — the actual sequence retrieval it normally enables is not
used here.

### step:panel_sizes

Cuts the name and length columns out of each `.fai` to produce a plain
`chromosome<TAB>length` table. The UCSC chain and net tools in Phase C require this for
both genomes of every pair, and they will not run without it. Lengths are unaffected by
soft-masking, so this can be built from the raw assemblies.

### step:panel_self_pairs

Writes one `X_X` line per strain, e.g. `PvW1_PvW1`.

Phase C aligns every genome against every other by taking the cross-product of the genome
collection with itself. That grid includes the diagonal — each genome paired with itself —
which is meaningless to align. This file is the list of cells to drop, fed to Galaxy's
`__FILTER_FROM_FILE__`. For 8 strains it removes 8 of the 64 cells, leaving 56.

It is a file because Galaxy's collection operations filter by
matching element identifiers against file contents, and the identifiers are only known once
the panel is chosen.





### step:anchor_prep

Reads each anchor's curated GFF3 and writes the two annotation files the projection
workflow (WF-C2) consumes. `gffread --bed` turns the mRNA/CDS structure into BED12 --
one line per transcript, with the exon block sizes and offsets packed into columns 11
and 12 -- and the step keeps only protein-coding transcripts, keyed by transcript id in
column 4. Alongside it, the isoforms table is built by walking every `mRNA` feature and
reading two attributes out of column 9:

    ID=PVW1_000005000_t1;Parent=PVW1_000005000;...

`Parent` is the gene, `ID` is the transcript, so that line becomes:

    PVW1_000005000<TAB>PVW1_000005000_t1

TOGA2 wants this because it classifies each transcript separately, then has to decide the
fate of the *gene* -- and the only way it knows which projections belong to the same gene
is this table. For PvW1 that is 6,075 transcripts over 6,075 genes; PAM 6,462; PvSY56
5,328 (this panel has one isoform per gene, so the counts match, but the format allows
many transcripts per gene).

The reason both files come from one step on one GFF3 is that their transcript ids have to
agree exactly. Build the BED12 from one annotation release and the isoforms table from
another and TOGA2 will project transcripts it cannot map back to any gene, and drop them.

### step:panel_relabel_map

Writes a two-column table mapping `A_B` to `A.B` for every ordered pair, 64 rows for 8
strains:

    PvW1_PvP01<TAB>PvW1.PvP01
    PvW1_PAM<TAB>PvW1.PAM

**Why it is needed.** When Galaxy builds the all-against-all grid, it names each cell by
joining the two element identifiers with an underscore: `PvW1_PvP01`. It cannot use a dot,
because the dot is reserved as a separator in collection identifiers. But Phase E — the
orthology step much further downstream — expects pair ids in the dotted form `PvW1.PvP01`,
which is the convention the reciprocal-best chain files use.

So this is a rename table, applied by `__RELABEL_FROM_FILE__` after the alignment grid is
built. Without it the chains come out named in a form Phase E does not recognise, and the
orthology step silently finds nothing to join on.

It is generated here, from the panel's own element identifiers, rather than hand-written,
so that a new panel needs no manual bookkeeping.

### step:multiqc_report

Folds the per-strain BUSCO summaries and the sourmash heatmap into a single HTML report, so
the panel can be reviewed in one place instead of opening a dozen datasets. This is the
artefact to look at before committing compute to Phase B and C.

Two outputs. The **report** is a self-contained interactive HTML page, a couple of megabytes,
so it is not reproduced here. The **stats table** carries the numbers behind it, and is worth
reading directly:

    Sample   busco-complete  ..._single_copy  ..._duplicated  fragmented  missing

The denominator is the size of the BUSCO lineage set — 446 genes for `apicomplexa_odb10` —
so `complete` divided by 446 is the completeness percentage. `single_copy` versus
`duplicated` splits the complete genes by copy number; a high duplicated count can mean a
genuine expansion or an assembly that has not collapsed its haplotypes. `fragmented` genes
were found only in part, and `missing` ones not at all.

On the Pv4 panel this is where two strains stand out: MHC087 at 388/446 (87.0 %) with 18
fragmented and 40 missing, and PvSY56 at 401/446 (89.9 %) with 17 fragmented and 28 missing.
Everything else sits at 443–446. Those two are the strains to be careful with downstream, and
PvSY56 is the same genome the sourmash matrix flags as a similarity outlier — two independent
measurements agreeing that it is the odd one out.

## B

<!-- WF-B softmask -->

### summary

Soft-masking: four maintained low-complexity maskers run in parallel, each emitting a
content-annotated BED6 hub track; their union soft-masks the assembly.


### description

Genomes are full of repetitive and low-complexity sequence: homopolymer runs, short
tandem repeats, satellite arrays. Align without handling it and those regions generate
enormous numbers of meaningless matches, because a stretch of AT repeats matches every
other stretch of AT repeats in the genome. This workflow finds that sequence and
lowercases it. Lowercasing rather than cutting is deliberate: the bases stay where they
are, so every coordinate downstream is unchanged and nothing is lost, while alignment
tools that honour soft-masking can avoid seeding matches inside them. It takes the raw
genomes and runs four independent repeat finders over each one, then merges everything
any of them flagged and lowercases the union. Four tools rather than one because they
disagree with each other, and taking the union is the conservative choice. What comes
out is the soft-masked genome set that every later step uses, one browser track per tool
showing what it flagged and what kind of repeat it is, and a table of how much of each
genome each tool masked.



### step:assemblies

The raw panel genomes, same collection WF-A takes. `uppercase` maps over it and everything
after that maps over `uppercase`'s output, so all eight strains are masked independently and
in parallel.

### step:uppercase

Rewrites every sequence line in upper case, so that what this workflow publishes carries this
workflow's mask and nothing else. Assemblies frequently arrive already soft-masked — PvP01 in
the captured run does — and while each masker below upper-cases its own working copy anyway
(the BED tracks were always de novo), `maskfasta -soft` lower-cases the union *on top of*
whatever was lower-case already. Without this step the published FASTA is the union of two
masks with no way to tell them apart, and the tracks beside it describe only one of them.
Header lines are left byte-for-byte alone, because scaffold names are case-sensitive and
renaming them would break every downstream join.

### step:dustmasker

NCBI's low-complexity finder, the one BLAST uses. It looks for stretches whose composition
is skewed enough to align to almost anything. Our wrapper adds the fourth and fifth BED
columns: instead of just "this region is low-complexity", each interval is labelled with the
repeat unit it is made of — `polyA`, `(AT)n`, `lc` for anything without a clean period — and
scored by how pure that repeat is, out of 1000. That is what makes these usable as a browser track.

### step:windowmasker

Counts how often each k-mer occurs across the whole genome and masks the ones that are
over-represented. Unlike dustmasker it needs no prior repeat library: it learns what is
repetitive in *this* genome. It is consistently the most aggressive of the four here,
masking around 29 % of a *P. vivax* genome against dustmasker's 13 %.

### step:tantan

Finds tandem repeats using a probabilistic model, so it catches decaying repeats that a
strict periodicity test misses. Sits between dustmasker and windowmasker in aggressiveness.

### step:fastan

Finds tandem arrays and reports the size of the repeating unit. The most conservative of the
four, typically masking under 5 %, and the only one that reports a strand.

### step:union_cat

Concatenates the four BED files into one. No merging happens yet, so the intervals overlap
heavily — for one genome this is roughly 520,000 rows built from four sets of about 110,000
each. The four callers disagree, and taking the union rather than an intersection is the
conservative choice: anything any tool considers repetitive gets masked.

### step:sort_bed

Sorts the concatenated intervals by chromosome and start. `bedtools merge` requires sorted
input and will produce silently wrong output otherwise, so this is a correctness step, not
tidiness.

### step:merge_bed

Collapses the overlapping intervals into disjoint ones. The 520,000 concatenated rows become
about 249,000 merged regions. Note the output is plain BED3 — the repeat-unit labels and
scores are dropped here, because once four callers' intervals are fused the label no longer
belongs to any single call. The labelled per-tool files are kept separately for the browser.

### step:maskfasta

The actual masking. `bedtools maskfasta -soft` lowercases every base inside a merged interval
and leaves everything else untouched. Nothing is deleted and no coordinate moves, which is
why every downstream phase can use these genomes interchangeably with the raw ones.

### step:faidx

Indexes the soft-masked genomes. Downstream tools that read a region out of a FASTA need
this, and Phase C needs the sequence lengths it records.

### step:sizes_cut

Cuts the name and length columns from each `.fai` into a `chromosome<TAB>length` table.
Masking does not change any length, so these match the values WF-A derived from the raw
assemblies.

### step:gcov_dust

`bedtools genomecov` over the dustmasker intervals, reporting what fraction of each sequence
is covered. The columns are chromosome, depth, bases at that depth, sequence length, and the
fraction — so the `depth 1` row is the masked fraction and `depth 0` is the untouched
remainder. There is one of these per masker plus one for the union, and they are what the
summary table is built from.

### step:gcov_window

Coverage for windowmasker, same shape as the others.

### step:gcov_tantan

Coverage for tantan.

### step:gcov_fastan

Coverage for fastan.

### step:gcov_union

Coverage for the merged union — the number that actually matters, since it is the fraction of
the genome that ends up lowercased.

### step:masking_table

Reduces all of that to one row per strain and one column per masker, as a percentage:

    Sample  dustmasker  windowmasker  tantan  fastan  union
    PvP01        13.63         28.63   15.64    4.18   39.01

Read across a row to compare the callers, down a column to compare strains. The union is
always the largest and is not the sum, because the callers overlap heavily. This is the
table to check when a genome behaves oddly downstream: an unusually high union means an
assembly full of repeats, and an unusually low one can mean masking silently failed.

### step:masking_multiqc

Wraps the table into a MultiQC report so the masking can be reviewed alongside the other QC
without reading a TSV.

## C

<!-- WF-C align_chain -->

### summary

Pairwise alignment and chaining across the full ordered strain grid: KegAlign/LASTZ →
axtChain → chainNet → netChainSubset → reciprocal-best, one-click via native collection
operations.


### description

To compare two genomes you first have to establish which piece of one corresponds to
which piece of the other. This workflow aligns every genome against every other, all 56
ordered pairs for 8 genomes, and then assembles the thousands of short local alignments
into chains and nets. A chain is a run of alignments in consistent order and
orientation; the netting step picks the best chain for each region and discards the
rest. That turns a haystack of local hits into a statement about large-scale
correspondence, so a real inversion or translocation shows up as structure rather than
disappearing into noise. It takes the soft-masked genomes and their chromosome lengths.
What comes out is, for each ordered pair, the cleaned chains, the raw pairwise
alignments, and a reciprocal-best set, which keeps only the cases where two regions are
each other's best match. That reciprocal condition matters: it is the difference between
two regions being genuinely the same locus in two strains, and one of them being a
paralog elsewhere in the genome. This is by far the most expensive step, and the
alignment runs on a GPU.



### step:masked_fastas

The soft-masked genomes from WF-B. Masking matters here: the aligner is meant to skip the
lowercased repeat regions when seeding, which is what stops a stretch of AT repeats matching
every other stretch of AT repeats in the genome.

### step:sizes

Per-strain `chromosome<TAB>length` tables from WF-A. The UCSC chain and net tools need the
lengths of both genomes in a pair and will not run without them.

### step:xprod_fa

Builds the alignment grid. Step by step:

**What goes in.** One collection, `masked_fastas`, with 8 elements — `PvP01`, `PvW1`, `PAM`,
`PvSY56`, `Sal-I`, `PvT01`, `PvC01`, `MHC087`.

**What comes out.** **Two collections of 64 elements each**, `output_a` and `output_b`.

**They share their identifiers.** Both are keyed `PvP01_PvP01`, `PvP01_PvW1`, `PvP01_PAM`, …
The identifier names the *cell of the grid*; which genome sits inside it depends on which of
the two collections you look in.

**Each element is a whole genome.** Every one of the 64 elements on each side is a complete
soft-masked FASTA of about 29 MB, so each side of the grid is ~1.84 GB.

**The contents differ, and the FASTA headers show it.** The accession prefixes identify the
genome unambiguously: PvP01 contigs start `LT6356…`, PvW1 `CAJZCX…`, PAM `CASCJQ…`.

    cell PvP01_PvW1
       output_a -> tgt_fa   >LT635626.1  >LT635627.1  >LT635612.2         PvP01
       output_b -> qry_fa   >CAJZCX010000001.1  >CAJZCX010000002.1  …     PvW1

    cell PvP01_PAM
       output_a -> tgt_fa   >LT635626.1  >LT635627.1  >LT635612.2         PvP01 again
       output_b -> qry_fa   >CASCJQ010000001.1  >CASCJQ010000002.1  …     PAM

    cell PvW1_PvP01
       output_a -> tgt_fa   >CAJZCX010000001.1  …                         now PvW1 is target
       output_b -> qry_fa   >LT635626.1  …                                and PvP01 is query

Two things are visible there. `output_a` is identical across the first two cells, because
PvP01 is the target in all 8 of its cells — that is the "each genome repeated" pattern in the
actual data. And `PvW1_PvP01` is `PvP01_PvW1` with the two sides exchanged, which is why both
directions are kept as separate cells: chaining is not symmetric, since the net is built with
respect to whichever genome is the target.

{{figure:cross_product}}

**Then the diagonal goes.** `tgt_fa` filters `output_a` and `qry_fa` filters `output_b`,
both against the same `self_pairs` list, so both drop the same 8 cells:
`PvP01_PvP01`, `PvW1_PvW1`, `PAM_PAM`, `PvSY56_PvSY56`, `Sal-I_Sal-I`, `PvT01_PvT01`,
`PvC01_PvC01`, `MHC087_MHC087`. **64 in, 56 out, 8 discarded** — and you can see both
collections in the discarded output, identical.

**How it is consumed.** `kegalign` maps over `tgt_fa` and `qry_fa` together. For element *i*
it takes the target from one and the query from the other. Cell `PvP01_PvW1` therefore
aligns PvP01 against PvW1, and the name of the cell tells you which way round it is.

**Why this is worth understanding.** The pairing is carried entirely by *position*. Nothing
in the data says these two collections belong together — Galaxy simply hands a tool the *i*-th
element of each. If some later step introduces a collection whose order differs, every cell is
silently handed another cell's data, with correct-looking identifiers throughout. That is not
hypothetical: it happened in WF-C2, where chains filtered out of this workflow arrived in
alphabetical order while the grid was in panel order, and every projection got the wrong
chain.

### step:xprod_sz

The same cross-product over the `sizes` collection, so every cell also has the two size files
its pair needs, in the same order.

### step:tgt_fa

Drops the diagonal. The cross-product includes each genome paired with itself, which is
meaningless to align, so `__FILTER_FROM_FILE__` removes the 8 self-cells using the
`self_pairs` list from WF-A: **64 cells in, 56 out, 8 discarded**. This is the target side.

### step:qry_fa

The query side of the same filtered grid, kept aligned with the target side.

### step:tgt_sz

Target-side chromosome sizes, filtered identically so the sizes stay in step with the FASTAs.

### step:qry_sz

Query-side chromosome sizes.

### step:kegalign

The GPU aligner. It finds the seed matches between the two genomes of a pair — the expensive
part of the whole phase, and the reason Phase C is GPU-gated.

### step:batched_lastz

Turns the seeds into actual gapped local alignments and writes them as **axt**, one file per
pair, tens to hundreds of megabytes each. These are raw local alignments: many of them,
unordered, and overlapping. Everything that follows is about imposing structure on them.

### step:axtchain

Assembles the local alignments into **chains**: runs of alignments in consistent order and
orientation, allowing for gaps. A chain is the claim that this stretch of the target and that
stretch of the query are the same locus, read the same way. The header preserves the scoring
matrix and gap penalties used, so a chain file is self-documenting.

### step:chainsort_clean

Sorts chains by score, highest first, which is what the netting step expects.

### step:chainprenet

Removes chains that are already fully covered by something better, so the netting step is not
wasting work on redundant candidates.

### step:chainnet

Builds the **net**: for every region of the target, picks the single best chain, then fills the
gaps inside it with the next best, recursively. This is what turns a pile of overlapping
chains into one coherent statement about correspondence, where an inversion or a translocation
shows up as structure rather than disappearing into noise. It emits a net for each direction.

### step:netchainsubset

Pulls back out the subset of chains that the net actually kept. The output is a chain file
again, but now non-redundant.

### step:chainstitchid

Joins chain fragments that belong together but were broken apart, giving each surviving chain
a stable id.

### step:chainstitchid_clean

Joins chain fragments that belong together but were broken apart by the netting, and gives
each surviving chain a stable id. This is the cleaned chain set for the pair.

### step:relabel_cleaned

Renames the cells from `A_B` to `A.B` using the map from WF-A, because Phase E expects the
dotted form. Compare the element identifier here (`PvP01.PvW1`) with the one two steps
earlier (`PvP01_PvW1`) — the data is unchanged, only the label.

### step:rb_swap1

Starts the **reciprocal-best** branch. It transposes the chain so the roles of target and
query are exchanged: look at the first line of the output against the same line in
`relabel_cleaned` and the two coordinate blocks have swapped places. Running the whole netting
procedure again from the query's point of view is what makes the result reciprocal.

### step:rb_sort1

Sorts the swapped chains for the second netting pass.

### step:rb_net

Nets again, now from the query side.

### step:rb_subset

Takes the chains that survived the second net.

### step:rb_stitch

Stitches the surviving fragments, as on the forward pass.

### step:rb_swap2

Swaps target and query back, so the result is expressed in the original orientation again.

### step:rb_sort2

Final sort of the reciprocal-best chains.

### step:relabel_rbest

Renames to the dotted form, as for the cleaned chains.

**What reciprocal-best buys.** A plain net says "this is the best match for this target
region". Reciprocal-best says the two regions are each other's best match, in both directions.
That is the difference between two regions genuinely being the same locus in two strains, and
one of them being a paralogue somewhere else in the genome — which is why Phase E builds its
orthology graph from these rather than from the cleaned chains.

## C2

<!-- WF-C2 project_annotations -->

### summary

Annotation projection: each anchor's curated genes are lifted onto every other strain
with Liftoff, then triaged and merged into a per-pair classification.


### description

Most genomes in a panel have no gene annotation worth trusting. This workflow borrows
one. Each anchor genome comes with a curated set of gene models, and Liftoff carries
those genes onto every other genome: it takes the sequence of each annotated gene out of
the anchor, aligns that sequence against the target genome, and places the gene wherever
it finds the best match, keeping the exon structure intact.

Every projected gene is then examined, because a gene landing somewhere is not the same
as a gene surviving: the check asks whether it still has an intact reading frame, or
whether it is interrupted by a premature stop, shifted out of frame, or missing exons.

Genes that fail that check are not discarded. They go to a second pass, **TOGA2**, which
re-projects them with CESAR2 using the whole-genome chain alignments from WF-C. This is
the only reason WF-C feeds WF-C2 at all: Liftoff does its own gene-by-gene alignment
internally and needs nothing from WF-C, so pass 1 alone would leave the two halves of
Phase C independent. TOGA2 is worth the cost because it grades the outcome:
intact, partially intact, lost, or lost in a way the alignment cannot resolve. On the first verified pair that turned 4,707 usable gene calls
into roughly 8,900, and separated 1,876 genes that are genuinely absent from the far
larger set that had simply never been looked at properly.

It takes the anchor genomes with their gene models, the genomes being annotated, the
soft-masked genomes for both sides, a gene-to-transcript table per anchor, and WF-C's
cleaned chains. The unmasked sequence is what Liftoff aligns against; the soft-masked
sequence is what the intactness check and CESAR2 use.

What comes out, for each anchor and query combination, is a GFF3 of the projected genes
and a table classifying each one, tagged with which of the two passes produced it. That
classification is the raw evidence the orthology step downstream uses to decide which
genes are shared across the whole panel, and it weights a CESAR2 call slightly above a
Liftoff one because it is better evidenced.

One limitation to be clear about: this is projection. A gene present in a query genome but
absent from every anchor cannot be discovered by either pass.


### step:anchor_assemblies

The three genomes with curated annotation — PvW1, PAM, PvSY56 — unmasked. Liftoff uses these
as its reference coordinate system.

### step:anchor_gene_gff3s

The curated gene models for each anchor. Two hard requirements: the seqids must match the
anchor FASTA headers, and gene-level types must be renamed to `gene`, because Liftoff's
default mode finds zero `gene` features in a native `protein_coding_gene` / `ncRNA_gene` /
`pseudogene` GFF3 and silently projects nothing.

### step:anchor_bed12s

BED12 transcript models per anchor, used as the reference bed by triage, merge and TOGA2.
Column 4 holds **transcript** ids (`PVW1_000005000_t1`), matching `anchor_isoforms`.

### step:assemblies

The whole panel, unmasked. These are Liftoff's targets — the genomes receiving genes.

### step:query_masked

The soft-masked panel from WF-B. The intactness check and CESAR2 read these; Liftoff aligns
against the unmasked copies.

### step:anchor_masked

Soft-masked versions of the three anchors. TOGA2 wants soft-masked sequence on both sides,
so the anchor side is needed here even though Liftoff uses the unmasked copy.

### step:anchor_isoforms

A gene-to-transcript table per anchor, grouping the BED12 transcripts into genes. TOGA2 needs
it to tie projections back to genes; the transcript ids must match column 4 of the BED12
exactly.

### step:cleaned_chains

WF-C's cleaned chains for all 56 ordered strain pairs, keyed `target.query`. Only the TOGA2
pass reads these — Liftoff does its own gene-level alignment — so this input exists solely
for pass 2.

### step:gen_anchor_self_pairs

Emits `A_A` for each anchor. The projection grid is anchors × strains, so it contains three
cells where an anchor would project onto itself. This list is what removes them, and it is
derived from the anchor element identifiers rather than supplied as a file.

### step:p_grid_fa

Crosses the three anchors with all eight strains: **24 cells**, keyed `anchor_query`. Same
mechanics as WF-C's grid, except the two inputs are different collections, so the diagonal
appears only where an anchor is also a panel member — which is why the self-pair list is
generated from the anchors rather than assumed.

{{figure:cross_product}}

### step:p_grid_gff

The same 24-cell grid over the anchor GFF3s, so each cell carries the annotation matching its
anchor.

### step:p_grid_bed

The same grid over the anchor BED12s.

### step:p_grid_ancmask

The same grid over the soft-masked anchors, for TOGA2's reference side.

### step:p_grid_iso

The same grid over the isoforms tables.

### step:p_anc_fa

Drops the three anchor self-cells: **24 in, 21 out, 3 discarded**. Every parallel grid is
filtered against the same list, so all of them stay aligned cell for cell.

### step:p_qry_fa

Query-side genomes for the 21 cells.

### step:p_anc_gff

Anchor annotations for the 21 cells.

### step:p_anc_bed

Anchor BED12s for the 21 cells.

### step:p_qry_masked

Soft-masked query genomes for the 21 cells.

### step:p_anc_masked

Soft-masked anchor genomes for the 21 cells.

### step:p_anc_iso

Anchor isoforms tables for the 21 cells.

### step:gen_chain_grid

Bridges WF-C's chain collection onto this grid. WF-C emits one chain per ordered strain pair
keyed `{target}.{query}`; this grid is keyed `{anchor}_{query}`. TOGA2's `--chain_file` is
reference-to-query with the anchor as reference, so cell `{anchor}_{query}` needs chain
`{anchor}.{query}`.

It emits three small files, all derived from the two collections' element identifiers:

- **keep** — the 21 `{anchor}.{query}` ids to select
- **relabel** — `{anchor}.{query}` → `{anchor}_{query}`
- **order** — the same 21 ids in cross-product order

### step:p_chain_keep

Selects the 21 chains this grid needs out of WF-C's 56: **21 kept, 35 discarded**. The
discarded ones are the pairs whose target is not an anchor.

### step:p_chain

Renames the selected chains from the dotted form to the underscore form, so their identifiers
match the rest of the grid.

### step:p_chain_sort

Re-sorts the chains into cross-product order, and it is not optional.

Galaxy pairs collections in a map-over by **position**. Every other grid here comes from
`__CROSS_PRODUCT_FLAT__` and is in anchor-collection order; a collection filtered out of
WF-C's 56 keeps WF-C's own alphabetical order. Before this step existed the two disagreed,
and every cell was handed another cell's chain:

    ref_bed12   position 0   PvW1_PvP01
    chain       position 0   PAM_MHC087     <- wrong cell

TOGA2 then reported `Processed 0 chains` and exited in a minute or two instead of eighty.
Nothing looked wrong: the identifiers were all correct and the interface showed the right
names throughout. Only the positions were mismatched. A single-cell test cannot catch this,
because with one element position 0 is trivially correct.

### step:p_liftoff

Pass 1. Pulls each annotated gene's sequence out of the anchor, aligns it against the target
genome, and places the gene at the best match, preserving exon structure. It emits the lifted
annotation and a list of features it could not map.

### step:p_triage

Examines every projected gene, because a gene landing somewhere is not the same as a gene
surviving. Eight rules check the reading frame, sequence identity, reference coverage, copy
number, partial mappings, splice sites, subtelomeric position, and known variant-antigen
families. Genes that pass go forward as clean; the rest are flagged.

`needs_cesar2.bed` is **empty in every cell**, and that is a known defect rather than a
result: triage looks the flagged genes up in the anchor BED12 by gene id, but that file is
transcript-keyed as TOGA2 requires, so nothing ever matches. It is harmless here, because the
C.4 redesign has TOGA2 classify the full anchor annotation instead of a triage subset, so
nothing reads the file.

### step:p_toga2

Pass 2. Re-projects with CESAR2 over WF-C's chains, and grades every gene rather than passing
or failing it: `FI` fully intact, `I` intact, `PI` partially intact, `UL` uncertain loss,
`L` lost. This is where genes Liftoff could not place cleanly are recovered.

Expensive — 47 to 178 minutes per cell — and 6 of the 21 cells failed here on an upstream
TOGA2 v2.0.8 defect (hillerlab/TOGA2#41) triggered by a single anchor transcript.

### step:p_merge

Folds both passes into the final annotation, tagging each call `source=liftoff` or
`source=cesar2` with its intactness class. The GFF3 is the deliverable; the classification
table is the evidence Phase E consumes.

## I

<!-- WF-I multiz -->

### summary

Builds one multi-way alignment per hinge strain: convert each pairwise alignment to MAF,
then fold every query onto the hinge's coordinates closest-first, so each hinge ends up
with a single file holding all eight strains.

### description

Everything up to this point aligns genomes two at a time. That is enough to ask whether a
region is shared between a particular pair, but not to ask what a region looks like across
the whole panel at once — for that you need every strain laid out against a common
coordinate system, one column per position. That is what this workflow produces.

Pick one strain as the reference — the hinge. Take its alignment against each of the other
seven, convert each to MAF, and then merge them one at a time so that every new strain is
threaded through what is already aligned. The result is a single file for that hinge in
which each block is a stack of rows, one per strain that could be aligned there. Repeat for
every strain in turn, because a region missing from one strain's view may be present in
another's, and there is no single correct reference.

The order of the merge matters more than it looks. Each new strain has to thread through
the alignment built so far, so folding a distant strain in early costs you blocks that
never recover. The queries are therefore folded closest-first, using the sourmash
similarities from WF-A — which is why the inventory step's similarity matrix is an input
here and not just a QC figure.

On the Pv4 panel the eight hinges come to 15.1 GB of alignment. Reading the PvW1 file: 94,948
blocks and 195,329 rows, with all eight strains present. PvSY56 contributes the fewest rows
of any strain, 6,186, which is what you would expect from a genome sitting at 0.24 similarity
to the rest of the panel while they sit at 0.63 to 0.70 among themselves.

### step:axt_to_maf

Converts each pairwise alignment from AXT to MAF, which is the format multiz reads. This
runs once per (hinge, query) pair — 56 jobs for eight strains — and the two `.sizes` inputs
are what let it write the chromosome lengths into each block. The hinge is always the target
(top) row, which is the property the fold depends on.

The prefix parameters are left empty here, so the rows come out carrying bare contig
accessions such as `QMFC01000014.1`. `multiz_fold` fixes the names on the way in; see below
for why that matters.


### step:multiz_fold

Takes one hinge's seven pairwise MAFs and merges them into a single multi-way MAF, folding
the closest query first. The fold is iterative — `multiz prev.maf next.maf 1` — starting
from the closest query's alignment and threading each next-closest onto the running result.
The header of the output records the order it used.

Two things this step has to get right, both of which were broken until 2026-08-05:

**It has to know which strain is the hinge**, because that is how it looks up the row of
similarities in `compare.csv` that determines the fold order. The name arrives as the
contents of a one-line dataset rather than as a text box: Galaxy cannot route a map-over
element identifier into a text parameter, and connecting a collection to one silently
substitutes the internal representation of the collection element instead of its name.

**It has to give every row a strain-qualified name.** multiz decides which species a row
belongs to by reading the text before the first `.` in the row's name, so rows named
`QMFC01000014.1` read as one species per contig. This step therefore rewrites the names to
`strain.contig` as it stages each pairwise MAF — the target row takes the hinge's name, the
query row takes the collection element's. This failure was silent: multiz exits successfully
and writes a structurally normal file whose rows simply cannot be told apart.

## E

<!-- WF-E consensus orthology -->

### summary

Groups genes into orthologous families across the panel: reciprocal-best alignments say
which loci correspond, the projections say which genes sit there, and connected components
of that evidence become the orthogroups.

### description

Up to this point every comparison has been between two genomes at a time. This workflow
turns those pairwise statements into panel-wide ones. If PvW1's gene X corresponds to
PAM's gene Y, and PAM's Y corresponds to PvT01's Z, then X, Y and Z are the same gene seen
in three strains — and the set of all such linked genes is an orthogroup.

Two kinds of evidence go in. The reciprocal-best chains from the alignment step say that
two stretches of two genomes are each other's best match, which is the classic way to
call orthology from whole-genome alignment; intersected with each strain's gene
annotation, that yields a gene-to-gene link. The projections add a second layer: where an
annotation was transferred onto a strain, the transferred gene tells you a gene is present
even when the alignment alone was ambiguous, and TOGA2 additionally says whether that gene
is intact or degraded.

The output is one row per orthogroup with a column per strain, labelled by how the group is
distributed: present once in every strain, present everywhere but duplicated somewhere,
a larger family, or restricted to a few strains. On the Pv4 panel 4,286 of 5,516 groups are
single-copy in all seven usable strains, which is what you would expect from strains of the
same species.

Earlier versions of this workflow also took a pangenome graph as a third evidence source.
It was removed: on the run that was checked, the graph contributed zero edges, and running
without it produces the same groups gene for gene.

### step:cls_ok

Drops projection cells whose job failed. The projection grid is not always complete — 15 of
21 cells on the run captured here — and a single failed element makes every downstream
collection operation refuse to run.

Note the distinction that matters in practice: this removes cells that *failed*. A cell
still **paused**, waiting on a failure further upstream, leaves the whole collection
unresolved and nothing after it can start at all. Those have to be cancelled, or the failing
cells fixed and re-run, before this workflow can proceed.

### step:cls_nested

Regroups the projection tables. They arrive as one flat list with names like
`PvW1_PvP01`, and the consensus step reads its tables from a directory tree that encodes the
reference and the target strain as two separate levels. This step splits the name on the
first underscore and uses the halves as those two levels, which is what lets the projection
workflow's output be connected directly rather than rebuilt by hand.

### step:gene_bed

Reduces each strain's annotation to a plain list of gene positions — chromosome, start, end,
gene id. Everything downstream works with these rather than the full annotation, and the
same collection is used three times: to place the alignment evidence, and to tell a
transferred gene apart from the native gene already annotated at that position.

### step:rbest_overlap

Intersects the reciprocal-best alignments with the gene positions. For each pair of strains
it walks the alignment blocks, projects each gene of the first strain through them, and asks
which gene of the second strain the projection lands on. A link is recorded when the two
genes cover each other by at least the required fraction, so partial or spurious overlaps do
not create one.

The result is close to a one-to-one matching — on the run captured here every gene had
exactly one partner per strain bar a couple of dozen exceptions — which is what makes the
grouping step defensible.

### columns:consensus.ortholog_table

orthogroup_id -- Identifier for this group, assigned in this run. It is not stable across runs: regroup the panel and OG000123 will be a different set of genes.
label -- How the group is distributed across the panel. CORE-1:1 = present exactly once in every strain. CORE-VAR = present in every strain, more than once somewhere. FAMILY = three or more copies in some strain. LINEAGE-SPECIFIC = confined to two strains or fewer. PARTIAL = anything else, i.e. present in some strains and absent from others.
n_strains -- How many strains have at least one gene in this group.
max_copies -- The largest number of genes any single strain contributes. This is what the label is computed from, so a bug that double-counts genes shows up here first.
clique -- How much of the mutual agreement this group could have does it actually have, from 0 to 1. In a real orthogroup every member should have found every other member independently. Scored against what the evidence can provide rather than against every possible pair, because the alignment evidence is one-to-one and each gene can link to at most one gene per other strain. A low value means the group is held together by a few links rather than by consensus -- sort ascending to find groups that were chained together and may be two families spliced into one.
density -- The raw fraction of member pairs that are linked, from 0 to 1. Unlike clique this is not comparable between groups of different copy number, which is why it is not the number to threshold on, but it is what exposes a group that has swollen: OG000001 in this run holds 822 copies in one strain and scores 0.629 clique against 0.0007 density. Sort ascending to find swollen groups.
{strain} -- One column per strain, holding that strain's gene ids in this group, comma-separated, or `-` when the strain has no gene here. Where several genes at the same position collapsed into one entry they are joined with `|`.

### columns:rbest_overlap.rbest_edges

strain_a -- The strain the gene was projected FROM.
gene_a -- Its gene id in that strain.
strain_b -- The strain the gene was projected ONTO.
gene_b -- The gene found there.
overlap_a -- The fraction of gene_a covered by the alignment, 0 to 1.
overlap_b -- The fraction of gene_b covered, 0 to 1. Both must clear the threshold for the row to exist, which is what makes these links reciprocal rather than one-sided.

### columns:cls_ok, cls_nested, c4_classifications

reference_gene_id -- The gene in the reference (anchor) strain that was projected.
query_gene_id -- What it was called in the target strain. For a transferred annotation this is the REFERENCE gene's name, not the target's own -- which is why the consensus step has to match it back onto the native gene by position.
source -- Which tool produced the call. `liftoff` = a whole-gene lift. `cesar2` = TOGA2's exon-aware placement, used to rescue genes Liftoff could not place. `none` = the reference gene was not found in this target at all.
intactness -- The state of the projected gene. `I` = intact. For TOGA2 calls the finer grades apply, best to worst: FI (fully intact), I, PI (partially intact), UL (uncertain loss), L (lost), M (missing), PM, PG, PP. `M` accompanies source=none.
query_chrom -- Contig in the target genome where the gene was placed.
query_start -- Start coordinate there.
query_end -- End coordinate there.
query_strand -- Strand, + or -.
orthology_class -- TOGA2's own orthology call for the projection: one2one, one2many, many2one or many2many. Liftoff rows carry `liftoff_clean`, and unplaced genes `unprojected`.

### step:consensus

Merges all the evidence into orthogroups and labels them.

The subtlety this step exists to handle is that a transferred gene keeps the *source*
strain's name. Project PAM's gene onto PvC01 and the result is still called `PVPAM_…`,
so it never matches the `PVC01_…` gene already annotated at that spot. Left alone, every
gene enters the graph twice, every group looks like it has two copies per strain, and the
labels — which are computed from copy number — become meaningless. Each transferred gene is
therefore matched back onto the native gene at the same position; one that lands on nothing
is a gene the annotation missed, and is dropped rather than counted as an extra copy.

Two numbers come out alongside each group as a warning. `clique` asks how much of the
mutual agreement the group could have does it actually have: in a real orthogroup every
member should have found every other member independently, and a low value means the group
is held together by a few links rather than by consensus. `density` is the raw fraction of
member pairs that are linked, which is what exposes a group that has swollen to thousands of
genes. Both are reported rather than acted on, because groups are only ever merged and never
split, so a wrong link cannot be undone after the fact — sort by `density` to find swollen
groups and by `clique` to find thin ones.

## COMMON

<!-- column meanings shared by more than one workflow -->

### columns:classifications, p_merge.classification_tsv, c4_classifications, cls_ok, cls_nested

reference_gene_id -- The gene in the anchor strain that was projected.
query_gene_id -- What it is called in the query strain. For a transferred annotation this is the ANCHOR gene's name, since the projection carries its source's identifier.
source -- `liftoff` for a whole-gene lift, `cesar2` for TOGA2's exon-aware placement, `none` if the gene was not found in this query at all.
intactness -- The state of the projected gene: `I` intact, or for TOGA2 calls the finer grades FI, I, PI, UL, L, M, PM, PG, PP from best to worst. `M` accompanies source=none.
query_chrom -- Contig in the query genome where the gene was placed.
query_start -- Start coordinate.
query_end -- End coordinate.
query_strand -- Strand, + or -.
orthology_class -- TOGA2's own orthology call: one2one, one2many, many2one, many2many. Liftoff rows carry `liftoff_clean` and unplaced genes `unprojected`.

### columns:panel_sizes, sizes, target_sizes, query_sizes

chrom -- Contig or chromosome name, exactly as it appears in that strain's assembly FASTA. If an annotation for the same strain uses different names, nothing will intersect -- that is why Sal-I contributes no orthology evidence.
length -- Its length in bases. Alignment tools need this to write coordinates that other tools can interpret.

### columns:panel_relabel_map, relabel_map

pair_underscore -- An ordered strain pair written `A_B`, which is how alignment outputs are keyed.
pair_dot -- The same pair written `A.B`, which is what the orthology step expects. The file exists only to rename one convention into the other.

### columns:panel_self_pairs, self_pairs

self_pair -- One `X_X` per strain. These are the diagonal of the all-against-all grid -- a genome aligned to itself -- and the list exists so those cells can be dropped before any alignment runs.

### columns:anchor_prep.isoforms, anchor_isoforms

gene -- Gene identifier.
transcript -- One transcript of that gene. A gene with several transcripts gets several rows. TOGA2 classifies transcripts individually and needs this table to decide the fate of the gene they belong to.

