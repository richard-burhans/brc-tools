# Phase C — `align_chain_project/`

Phase C used to be one combined workflow (`align_chain_project.gxwf.yml`) that
did **two independent things** with no data dependency between them. It is now
**split into two workflows** in this directory:

| Workflow | File | Does |
|---|---|---|
| **WF-C** `align_chain` | `align_chain.gxwf.yml` | pairwise alignment → UCSC chain pipeline (cleaned + reciprocal-best chains) |
| **WF-C2** `project_annotations` | `project_annotations.gxwf.yml` | anchor gene-annotation projection (Liftoff → triage → merge) |

They share no steps; run them independently (re-project without realigning, or
vice-versa). Both **run one-click** — a single invocation with collection
inputs, no per-pair driving and no custom enumeration helpers — using Galaxy's
native collection-operation tools. Both were validated **bit-identical** to the
earlier per-pair driving on the Pv4 test panel (cleaned `PvP01.PvW1`: 11726
lines, chain score 293706235; merged `PvP01→PvW1`: 856 mRNA), within KegAlign GPU
run-to-run variance of the `run_all.sh` ground truth.

## How the map-over works (native built-ins, no helpers)

The historical gxformat2 blocker was that a `list:paired` slot could not be
addressed into the two separate `data` inputs of axtChain / chainNet, and the
anchor×query grid could not be co-fanned. Both are solved natively:

- **`__CROSS_PRODUCT_FLAT__`** takes two `list`s and emits **two aligned flat
  lists** (`output_a`, `output_b`). Feed `output_a`→target and `output_b`→query:
  Galaxy's element-wise (dot-product) matching then maps the two-input tool over
  every pair. This *is* the paired-slot→two-input binding, expressed natively.
- **`__FILTER_FROM_FILE__`** (`remove_if_present`) drops the diagonal/self pairs
  (chains 25→20; projection 15→12).
- **`__RELABEL_FROM_FILE__`** rewrites the `A_B` pair ids to `A.B` (the
  `join_identifier` select offers only `_ : -`, and Phase E parses `{a}.{b}` on
  dots).

## WF-C — `align_chain`

`cross_product_flat(masked_fastas, masked_fastas)` → 25 ordered pairs → filter
self → 20. KegAlign(target,query) → batched_lastz → axtChain → chainSort →
chainPreNet → chainNet → netChainSubset → chainStitchId = **20 directed cleaned
chains** (both directions come straight from the cross product). The
reciprocal-best branch (swap → sort → chainNet with swapped sizes → subset →
stitch → swap → sort) yields the rbest chains. `relabel_from_file` rewrites both
to `A.B`.

| Input | Type | Source |
|---|---|---|
| `masked_fastas` | list | WF-B softmasked FASTAs (id=strain) |
| `sizes` | list | WF-A `.sizes` (id=strain) |
| `self_pairs` | txt | strain self-pair ids `X_X` to exclude (WF-A) |
| `relabel_map` | tabular | `A_B<TAB>A.B` id map for Phase E (WF-A) |

**Outputs:** `cleaned_chains` (20, id `A.B`), `rbest_chains` (id `A.B`),
`pairwise_axt`.

> **Note — `rbest_chains` emits 20 directed, not 10 unordered.** The cross
> product yields both `A.B` and `B.A`, so a reciprocal-best chain is computed for
> each direction. Harmless downstream — Phase E `rbest_overlap` keys on the
> `{a}.{b}` stem and the union-find dedupes the redundant direction — but for
> exact GT parity, add a `__FILTER_FROM_FILE__` reducing the rbest branch to an
> `A<B` id set.

## WF-C2 — `project_annotations` (Liftoff + TOGA2 rescue)

`cross_product_flat(anchor_assemblies, assemblies)` builds the anchor×query grid
(3×5 = 15 cells; plus parallel grids for the anchor GFF/BED12 and query
softmasked FASTA, all keyed `anchor_query` and mutually aligned) → filter 3 self
→ 12. The anchor self-cell id list (`A_A`) is generated **internally** by
`collection_self_pairs` on the anchor element identifiers (no input file).
Liftoff → phase_c2_triage → phase_c4_merge map over the grid. Genes that fail
triage are **not** dropped: `toga2` re-projects them with CESAR2 over WF-C's
cleaned chains, and the merge folds both passes (`use_toga: yes`), tagging each
call `source=liftoff` or `source=cesar2` with a TOGA2 intactness class.

`collection_anchor_grid` bridges the two workflows. WF-C emits one cleaned chain
per ordered strain pair keyed `{target}.{query}`; this grid is keyed
`{anchor}_{query}`. TOGA2's `--chain_file` is reference-to-query with the anchor
as reference, so cell `{anchor}_{query}` needs chain `{anchor}.{query}`. The tool
emits the selection list and rename map from the two collections' element
identifiers — no hand-authored config per panel.

**Pass 1 alone needs nothing from WF-C**: Liftoff aligns each gene itself. The
`cleaned_chains` input exists solely for the rescue pass.

`query_name` is a constant label: per-element scalars derived from a map-over
identifier are not expressible in gxformat2, but the per-cell query *data* is
correct via id-aligned collections, and Phase E keys on the collection element
id, not the internal label.

| Input | Type | Source |
|---|---|---|
| `anchor_assemblies` | list | anchor unmasked FASTAs (Liftoff reference) |
| `anchor_gene_gff3s` | list | anchor `gene.gff3` (gene-level types renamed to `gene`; Liftoff default mode finds zero `gene` features in the native `protein_coding_gene`/`ncRNA_gene`/`pseudogene` GFF3) |
| `anchor_bed12s` | list | anchor BED12 (triage/merge ref-bed) |
| `assemblies` | list | query unmasked FASTAs (Liftoff target = query) |
| `query_masked` | list | query softmasked FASTAs (triage query-fasta, TOGA2 `--query_2bit`; WF-B) |
| `anchor_masked` | list | anchor softmasked FASTAs (TOGA2 `--ref_2bit`; WF-B) |
| `anchor_isoforms` | list | anchor gene→transcript tables (TOGA2 `--isoform_file`; `anchor_prep`) |
| `cleaned_chains` | list | cleaned chains for every ordered pair, id `target.query` (WF-C) |

**Outputs:** `merged_annotations` (12, id `anchor_query`), `classifications` (12).

> **TOGA2 note.** The rescue pass is wired in and on. It is container-only —
> upstream publishes no bioconda package and no biocontainer — and now runs
> from a pullable image, `quay.io/galaxy/toga2` (`tools/toga2/container/`),
> rather than a SIF each deployment stages itself. The pin moved from v2.0.8 to
> **v2.0.9**, which is a bug fix, not housekeeping: through v2.0.8 TOGA2's own
> `execute_joblist.nf` did not parse under Nextflow 25.10+. The wrapper's CLI
> mapping and output schema were verified against v2.0.8, so re-run the
> comparison before treating v2.0.9 output as interchangeable. It is GPU-independent but CPU/IO-heavy:
> **82 minutes for a single cell** on 16 cores, driving its own Nextflow
> pipeline with the local executor inside one allocation.
>
> Verified on PvW1→PvP01: 16,874 classification rows against 10,782 for Liftoff
> alone — 4,707 `liftoff`/`I` plus 6,092 `cesar2` (FI 3,906, L 1,876, UL 170,
> I 135, PI 5). The full 21-cell grid has not yet been run with it enabled.

## History

The earlier port used three bespoke helper tools (`__pair_strains__`,
`__pair_sizes__`, `__cross_product__`) plus per-pair API staging for the
slot-binding and anchor×query gaps. Those were replaced by the native
collection-operation built-ins. The combined `align_chain_project.gxwf.yml` was
then split into the two independent workflows above (the chaining and projection
halves never shared data).
