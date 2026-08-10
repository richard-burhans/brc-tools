# Comparative & Population Genomics Program — Coordination Plan

Working doc capturing the coordination story across the group's parallel efforts.
Goal: one coherent program, no redundant wrappers or near-identical workflows, clean IWC submissions.

Status: DRAFT, 2026-07-30

Hosting context: ties into **BRC Analytics** — workflows target IWC and will be hosted/runnable on the BRC site, with the group's own runs and outputs viewable by others. Provenance via IWC workflow invocations. The BRC site is planned to gain a **pangenome view with orthogroup awareness** — the orthology ledger (§4.1) is the data contract that view consumes.

---

## 1. From assets you know to suites we're proposing

You know one or more of these assets. Nobody knows the suite names yet — we just coined them here. This table is the translation:

| Asset you know | → Suite | What it is | What changes |
|---|---|---|---|
| **Crypto protocol** (v2.1c docx) — Module 1 read-level popgen, Module 2 MK divergence | → **PLAIG** | Read-level population genomics suite for haploid pathogens | New: generalize from *Cryptosporidium* study design into reusable Galaxy workflows |
| **Crypto protocol** — Module 3 orthology/synteny/ancestry/gene-family layer (§8.1–8.4, 8.6–8.8) | → **ANCOR** | Cross-species comparative suite: OrthoFinder → synteny → ancestral reconstruction → gene-family evolution | Existing designs; add ledger output + assembly QC gate |
| **Crypto protocol** — Module 3 HyPhy layer (§8.5) | → **CAPHEINE** (extend) | Shared selection-analysis engine (HyPhy suite) | Already in IWC but viral-tuned; extend for eukaryotes + full method coverage |
| **Pangenome workflows** (`brc-tools`, `pangenome-helpers`) — PGGB, Liftoff/TOGA2, consensus orthogroups, per-OG BUSTED, UCSC hubs | → **PANTEON** | Within-species pangenome pipeline | Has a working prototype (unlike PLAIG/ANCOR which are plan-only); add ledger sidecars, reroute selection through CAPHEINE |
| **MalariaGEN Plasmodium proxies** ([github.com/asgiraldoc/MalariaGEN-Pf8-Data-Retrieval-Proxy](https://github.com/asgiraldoc/MalariaGEN-Pf8-Data-Retrieval-Proxy)) — FastAPI over MalariaGEN Zarr releases: gene/sample-filtered FASTA/VCF extraction, 10 interactive pop-gen workflows (PCoA, Fws, NJ tree, PLINK2 PCA, Fst, diversity, selection scans, haplotype network, allele-frequency views), IGV.js genome browser. Deployed as **one instance per species per port** (Pf8:9000, Pv5:9001, Pk1:9050, Pmb1:9051, Poc1:9052, Pow1:9053) — likely includes lab-generated data alongside public MalariaGEN releases | → **PLAIG** (analysis workflows) + **BRC site** (interactive views) + **UCSC** (browser) | De-facto requirements document for malaria pop-gen across the genus; the per-port-per-species deployment is exactly the sprawl BRC + PLAIG consolidates (one site, all organisms; one workflow family, parameterized per species) | Workflows translate to Galaxy (`popgen-analysis` toggles); interactive views become BRC site features; browser becomes UCSC track hubs. See §12 and [`plaig-design.md`](plaig-design.md) gap analysis |

**Key insight:** ANCOR ≈ crypto Module 3's orthology/synteny/ancestry/gene-family layer, generalized. Module 3's HyPhy layer is CAPHEINE's job, not ANCOR's. PANTEON is the within-species layer the crypto protocol lacks. CAPHEINE sits underneath both PANTEON and ANCOR as the shared selection engine.

> Suite names are working titles — tentative, changeable, but better than generic letters.

Companion docs (details live there, decisions live here): **[`plaig-design.md`](plaig-design.md)** (PLAIG design), **[`ancor-reference.md`](ancor-reference.md)** (ANCOR reference), **[`panteon-reference.md`](panteon-reference.md)** (PANTEON reference), **[`capheine-updates.md`](capheine-updates.md)** (CAPHEINE change spec).

## 2. The three-timescale story (program narrative)

One umbrella program, three suites at three strata of the same question:

- **PLAIG — Pathogen Lineage Analysis & Inference of Genotypes** (new, from crypto Modules 1–2) — read-level population genomics: reads → haploid VCF → polymorphism, structure, recombination, sweeps, demography, MK divergence
- **PANTEON — PANgenome Trees & Evolutionary Orthology Networks** (existing pangenome) — within-species assemblies: graph variation, PAV/CNV, consensus orthogroups, per-OG selection
- **ANCOR — Analysis of Networked Core Orthologs for Reconstruction** (existing designs) — cross-species assemblies: synteny, ancestral reconstruction, gene-family evolution

CAPHEINE sits underneath PANTEON and ANCOR as the shared HyPhy engine.

Each stratum is the others' control: read-level allele frequencies are immune to assembly-collapse artifacts; assembly PAV is immune to capture/WGA coverage bias; cross-species divergence contextualizes both.

### 2.1 What needs to change in PANTEON

PANTEON is the only timescale suite with a working prototype (PLAIG and ANCOR are design-only). Changes are integration, not rebuild: add GARD + SRV to the selection step, emit orthology-ledger-format sidecars from the consensus step, add an assembly QC gate, and parameterize hardcoded *P. vivax* assumptions. Full change list with workflow-level detail: **[`panteon-reference.md`](panteon-reference.md)**.

> Specific workflow names and step descriptions throughout this doc and its companions are tentative — actual workflow structure may change as planning continues.

## 3. CAPHEINE — extend, don't fork

Decision: **all HyPhy work routes through CAPHEINE.** No new HyPhy workflow family. CAPHEINE is currently viral-tuned; extending it for eukaryotic/pathogen reuse requires: mandatory GARD pre-check (crypto §8.5), full method coverage (add aBSREL/FUBAR/SLAC to existing MEME/PRIME/BUSTED/FEL/CFEL/RELAX), eukaryotic input handling (internal stops, multi-exon CDS), pluggable aligner stage (cawlign/MACSE/MAFFT behind a uniform output contract), and a ledger-based input adapter so per-orthogroup inputs from PANTEON or ANCOR feed in directly.

*Full change spec, pluggable preprocessing design, tool availability, and build sequencing: **[`capheine-updates.md`](capheine-updates.md)**.*

## 4. Shared foundation (the redundancy sinks)

### 4.1 Orthology ledger — common format for BRC integration

Both PANTEON (within-species) and ANCOR (cross-species) produce orthogroup assignments. They won't be compared directly — different taxon scopes — but both should emit the **same TSV format** so downstream consumers (CAPHEINE, the BRC site's planned orthogroup-aware view) can ingest either without custom adapters.

**Common emission format (the contract):**

```
orthogroup_id  label (CORE-1:1 | CORE-VAR | FAMILY)  <per-taxon gene id columns...>
```

- PANTEON `consensus` already emits ~95% of this (within-species, graph-based).
- ANCOR `orthofinder-core` emits it with one formatting step (cross-species).
- Label rules (presence fraction, copy-number cutoffs, tandem-collapse handling) defined once, versioned — "FAMILY" means the same thing regardless of engine.
- Consumed by: CAFE5 (count matrix), synteny/AGORA (CORE-1:1 filter), CAPHEINE (per-OG sequence extraction), crypto Module 2 (MK ortholog set).
- Provenance is tracked by Galaxy workflow invocations — no separate bookkeeping needed.

- Satisfies crypto protocol §4.3 "Unified Orthology Ledger".

### 4.2 Shared infrastructure (cross-cutting, not suite-specific)

These are preparatory and connective components — QC gates, data-format contracts, and coordinate-transfer machinery that enable the suites but aren't themselves the science:

- **`assembly-qc-gate`** — BUSCO ≥85% + contig N50 check. Intake step 0 for PANTEON and ANCOR. New, small, no existing coverage.
- **`read-qc-decontam`** — reuse existing IWC read-QC workflows + crypto's competitive-mapping pass 2 as an optional extension.
- **`haploid-read-to-vcf`** — PLAIG's all-sites haploid variant calling. Justified vs existing IWC workflows (those are variants-only, per-sample, no joint genotyping). Details: **[`plaig-design.md`](plaig-design.md)**.
- **Annotation harmonization** — PANTEON's Liftoff/TOGA2/triage layer; ANCOR gains an optional entry point calling it.
- **Liftover/chains** — produced once by PANTEON (`align_chain_project`, `vcf_projection`), consumed by crypto §4.3.

## 5. PLAIG — new suite from crypto Modules 1–2

PLAIG generalizes the crypto protocol's read-level popgen (Module 1) and MK divergence (Module 2) into a reusable Galaxy suite for any haploid-dominant eukaryotic pathogen. Four generic subworkflows (read QC, variant calling, population analysis, MK divergence) + a thin organism-specific parent. The population analysis subworkflow bundles structure, diversity, recombination, selection scans, and demography as toggleable steps sharing one VCF input. Largest user base, worst current Galaxy tooling story (pixy, selscan, pyrho, dadi-cli all need wrappers).

*Full design (scope, workflow decomposition, config bundle, caller policy, open questions): **[`plaig-design.md`](plaig-design.md)**.*

## 6. Wrapper ownership — no overlap

| Owner | Wrappers |
|---|---|
| ANCOR | MCScanX, CAFE5, AGORA, jcvi |
| PLAIG | pixy, selscan, pyrho, MACSE or PRANK, dadi-cli |
| Shared | `orthofinder-core` (design stage; core dependency OrthoFinder is IUC-wrapped), CAPHEINE extensions, `assembly-qc-gate` |

## 7. Generalization across organisms

Principle: the architecture is organism-agnostic; the real boundaries are **ploidy and reproductive mode, not taxonomy**. Config-level porting covers any haploid-dominant eukaryotic pathogen.

### 7.1 Already general (no changes needed)

CAPHEINE (coding sequences in, selection out), ANCOR (any clade with decent assemblies), the orthology ledger format, `assembly-qc-gate`, read QC/trimming, and UCSC hubs/chains/PGGB machinery are all organism-neutral by design.

### 7.2 Porting to a new organism = workflow parameters, not code changes

Each new organism should require only workflow parameters (reference assembly, gene-family lists, outgroups, masks, etc.) — no code changes. The principle: **only parameterize what can't be ascertained from the data itself.** Things like ploidy, reference strain, and decontamination panel must be supplied; things like gene counts and assembly stats can be computed. Per-suite parameter details: **[`plaig-design.md`](plaig-design.md)**, **[`panteon-reference.md`](panteon-reference.md)**.

### 7.3 Suite-specific generality notes

- **PLAIG**: constraint is a **haploid-dominant life stage**. Clean fit: Plasmodium, Toxoplasma, Babesia, Theileria, Trypanosoma, haploid fungi (Cryptococcus, Aspergillus). Edge case: **Leishmania** — haploid calling is standard but aneuploidy is first-class biology; promote the depth/CNV module from provisional to primary there. Diploid/polyploid organisms are out of scope for v1 — see [`plaig-design.md`](plaig-design.md) for what a diploid extension would require.
- **PANTEON**: class boundary is paralogy-rich/subtelomeric biology, not apicomplexa — T. brucei (VSG), Leishmania (amastins), fungal adhesin families all fit the triage design.
- **Out of scope by design**: bacteria (mature separate ecosystem: snippy/roary etc.), diploid pathogens and helminths (need a genuinely different QC/popgen layer), Giardia (ploidy weirdness). **State this scope explicitly in each suite README** — preempts "does this work for Candida?" IWC review questions.

### 7.4 Adjacent efforts (non-competing, credited borrowing)

Two nf-core surveillance pipelines tile adjacent space — no functional overlap with this program, and both are Nextflow-only (no Galaxy-side conflict):

- **`nf-core/pathogensurveillance`** — broad-shallow triage: manifest-driven, kmer species ID → best-reference mapping → SNP phylogeny in context of NCBI-fetched references → poppr-style assignment → HTML report. Answers "what is this sample and where does it sit?" Overlaps only PLAIG's front half (assignment/mapping); none of the inference layer (diversity, LD, sweeps, MK, tiers, mixed infection). Natural flow: pathogensurveillance triage → PLAIG deep inference. **Source of two borrowed patterns** (§7.2): reference usage flags and NCBI-query reference fetching.
- **`PlasmoGenEpi/plasmodiumdrugres`** (nf-core-bound) — deep-narrow surveillance reporting: PMO/allele-table input → prevalence of known resistance alleles and multilocus haplotypes per user-defined population. Downstream of amplicon calling, locus-targeted. Their "population assignment" = our manifest covariates. Sibling repos of note: `recombuddy` (simulates P. falciparum polyclonality/relatedness — candidate validation harness for PLAIG's mixed-infection handling) and the **PMO format** (candidate PLAIG export, §5).

Positioning: pathogensurveillance = broad-shallow identification triage; plasmodiumdrugres = deep-narrow surveillance reporting; this program = deep-broad evolutionary inference on the Galaxy/IWC side. They tile the space.

### 7.5 Community demand (euk-pathogen priority, for effort sequencing)

1. **PLAIG** — largest user base: drug-resistance surveillance, transmission/outbreak tracking, vaccine escape (malaria-scale communities); also the worst current Galaxy tooling story (ToolShed gaps, §8).
2. **CAPHEINE (extended)** — "which genes are under selection" is the universal second question; cheap extension, already branded in IWC.
3. **PANTEON** — PAV/antigenic-family demand rising; per-species power-user analysis.
4. **ANCOR** — smallest audience (methods-paper territory); CAFE5 component has broadest standalone demand.

## 8. ToolShed audit summary (2026-07-30)

Available (IUC unless noted): fastp, MultiQC, Kraken2 (+DMs), bwa_mem2, minimap2, samtools/bcftools, GATK4, PLINK (1.9-era), OrthoFinder, IQ-TREE, ASTRAL, **full HyPhy suite** (gard, slac, fel, meme, fubar, busted, absrel, relax, cfel), DIAMOND, BUSCO, CrossMap, mosdepth, quast, sra_tools, seqkit, sourmash; ADMIXTURE (dereeper, vet before use).

Missing (need wrappers): pixy, selscan, pyrho, LDhat, SweeD/OmegaPlus, dadi-cli, stairwayplot2, MACSE, PRANK, fineSTRUCTURE, MCScanX, SynNet, AGORA, CAFE5, EDirect.

Best-practice swaps adopted: BWA-MEM → BWA-MEM2; LDhat → pyrho-primary.

## 9. Documentation policy — "different tools for similar jobs"

Rule: where two suites use different tools for a similar job, the choice must be documented and justified in one place (this doc + workflow READMEs), so IWC reviewers see intentionality, not sprawl.

Current justified divergences:

| Job | PANTEON choice | ANCOR choice | Justification |
|---|---|---|---|
| Orthology | Graph consensus (PGGB + rbest edges) | OrthoFinder | Within-species paralogy resolution vs cross-species orthogrouping; ledger unifies outputs |
| Synteny | chain/liftover projections | MCScanX/jcvi | Within-species liftover vs cross-species block detection |
| Codon alignment | (via CAPHEINE) | (via CAPHEINE) | single engine |
| Selection | CAPHEINE BUSTED profile | CAPHEINE full suite | same workflow, different profiles — the model case |

## 10. Open questions

- [ ] Umbrella program name (suite brands set: ANCOR, CAPHEINE, PANTEON, PLAIG)
- [ ] CAFE5 ultrametric tree: require user-provided, or add TreePL/r8s dating step? (crypto stance: relative branch lengths unless calibrations defensible)
- [ ] Species tree for ANCOR/CAFE5/AGORA: OrthoFinder STAG tree vs IQ-TREE+ASTRAL (crypto prefers latter, gCF/sCF)
- [ ] fineSTRUCTURE worth wrapping? (painful deps; ADMIXTURE+PCA+LD may suffice)
- [x] ~~Does PANTEON route its existing BUSTED step through CAPHEINE?~~ Resolved §2.1: yes, reroute WF-H after CAPHEINE extension lands.
- [ ] hyphy-cln vs MACSE frameshift/stop retention — decide whether hyphy-cln needs a codon-tolerant mode or whether frameshift OGs bypass it with documented justification
- [ ] UCSC assembly hub output: PANTEON produces hubs (WF-K). Should PLAIG and ANCOR also emit hubs for visualization on the BRC site?
- [ ] Popgen results contract (§11.2): define now alongside the orthology ledger, or defer until PLAIG v1 workflows exist?
- [ ] BRC gene pages (§11.3): which assemblies get gene pages (reference only? all?), and what is the workflow launch list for v1?

## 11. BRC Analytics integration — division of labor across UCSC, Galaxy, and BRC

The Pf8 proxy review (2026-08-05) clarified how features from an interactive analysis site decompose across the three hosting planes BRC Analytics already partners with. The principle: **Galaxy computes, UCSC visualizes, BRC navigates and renders precomputed results.** No plane should duplicate another's role.

### 11.1 What goes where

| Capability | Plane | Rationale |
|---|---|---|
| **Variant visualization** (track hubs, locus browsing) | **UCSC** | BRC already deep-links to UCSC (`ucscBrowserUrl`, `UcscTrack` types). Precompute variant/diversity tracks as bgzipped VCF/bigBed/bigWig once; UCSC serves byte-ranges. Replaces the proxy's IGV.js + PLINK2-track pipeline. |
| **Population-genetics compute** (PCA, Fst, Fws, selection scans, NJ trees, haplotype networks) | **Galaxy** | PLAIG `popgen-analysis` toggles. Launched from BRC via existing workflow-landing machinery (`galaxy-api.ts`, `WorkflowInputsView`). |
| **FASTA/VCF extraction** (genotype-aware, per-sample) | **Galaxy** (or a slimmed proxy API) | Can't be static — injecting 33k samples' genotypes is real compute. Either a Galaxy workflow or the proxy's extraction API linked from gene pages. |
| **Gene search & resolution** | **BRC site** | Static search index over annotation (the proxy's `GENE_NAME_MAP` is just a GFF lookup). Natural fit for BRC gene pages. |
| **Metadata dashboard / sample filtering** | **BRC site** | Faceted catalog tables are BRC's core competency (findable-ui). Pf8 samples could become a catalog entity type. |
| **Precomputed result views** (PCoA with filtering, Fws-by-population, allele-frequency heatmaps) | **BRC site** | Galaxy workflows emit versioned result artifacts; BRC renders them interactively — same pattern as the orthology ledger (§4.1). |

### 11.2 The popgen results contract

The orthology ledger (§4.1) already established the pattern: workflows emit a versioned TSV/JSON contract; the BRC site consumes it for interactive views. A parallel **popgen results contract** would let PLAIG workflows produce:

- PCoA coordinates + variance proportions (NPZ/TSV)
- Fws-by-population tables (TSV)
- Allele-frequency-by-population/time tables (TSV)
- Pairwise Fst matrices (TSV)

…once, and the BRC site renders them with filtering/exploration — replacing the proxy's entire reason for existing as a live API. This is the within-species analog of the orthology ledger: a data contract between Galaxy outputs and BRC views.

### 11.3 Gene pages as the connective tissue

Genes/loci are the natural join key across all three planes. BRC's entity model currently stops at organisms → assemblies. **Gene pages** — even minimal ones — give every feature a stable anchor. A v1 gene page need only contain:

- Gene ID + coordinates (chrom, start, end, strand)
- Organism and assembly name
- Functional annotation (where available)
- Sequence file downloads (CDS, protein, genomic)
- A list of workflows runnable with this gene as input (e.g. kmindex query, lexicap, AlphaFold structure prediction)
- **UCSC deep link** to the gene's locus (derived from coordinates + the assembly's UCSC hub URL), giving immediate biological context — surrounding genes, existing tracks, etc.

No precomputed per-gene summaries (variant counts, diversity stats) are assumed for v1 — those can come later if the popgen results contract (§11.2) makes them cheap to generate. The immediate value is the **launch point**: each workflow button deep-links into Galaxy with the gene/locus prefilled, and the UCSC link provides the visual context. Without gene pages, BRC can only link at assembly granularity and the user re-does gene-hunting inside each partner tool.

A note on what BRC renders vs. what Galaxy computes: the precomputed result views in §11.1 (PCoA, Fws, freq tables) are not produced by a separate BRC compute path — they are the **outputs of Galaxy workflow runs**, surfaced on the BRC site via linked Galaxy user accounts. BRC is the navigation and rendering layer; Galaxy remains the sole compute plane. This means the popgen results contract (§11.2) is also a provenance contract: every result view traces back to a specific Galaxy invocation.

#### 11.3.1 Future: variant data on gene pages (proxy parity)

To fully replicate the Pf8 proxy's gene-level experience, a gene page would eventually need to present a **slice of the variant data** for that locus — e.g. a per-sample variant table, SNP-injected reference sequence, or a mini-VCF download filtered to the gene's coordinates and the user's sample filters. This is genotype-aware compute (not static annotation), so it requires either:

- A Galaxy workflow (VCF slice + filter → download), launched from the gene page with the locus prefilled, or
- A thin backend API (the proxy's extraction endpoint, or a successor) that the gene page calls directly.

This is explicitly **out of scope for v1 gene pages**. It depends on the popgen results contract (§11.2) or a persistent extraction service existing first. Listed here so the v1 design doesn't accidentally preclude it — e.g. the gene page's coordinate data should be structured in a way that a future "Extract variants" button can consume.

### 11.4 What this means for the suites

- **PLAIG** should emit UCSC track hubs (open question §10) and popgen results contract artifacts alongside its analysis outputs.
- **PANTEON** already produces UCSC hubs (WF-K); its orthology ledger sidecars are the model for PLAIG's results contract.
- The proxy's `analysis/` and `browser/` packages become redundant once PLAIG workflows + UCSC hubs + BRC result views exist. The proxy's extraction API may persist as a thin backend or be replaced by a Galaxy workflow.

## 12. Next steps

1. [ ] Group sign-off on this doc + orthology ledger spec (one page)
2. [ ] PLAIG: start `haploid-read-to-vcf` (all existing tools, highest reuse)
3. [ ] CAPHEINE: GARD + eukaryotic-scope PR (itemized in §3)
4. [ ] ANCOR: proceed per existing design doc (Option B modular), add `assembly-qc-gate` + ledger output format
5. [ ] PANTEON change list per §2.1 (GARD+SRV in WF-H first; ledger sidecars in WF-E)
6. [ ] Wrapper queue per §6 ownership
7. [ ] Popgen results contract spec (§11.2) — define output shapes for PCoA/Fws/freq/Fst alongside the orthology ledger
8. [ ] BRC gene pages: scope entity model, search index, and v1 workflow launch list (§11.3)

---

### Appendix: source documents

Companion reference docs (linked throughout this doc):
- [`plaig-design.md`](plaig-design.md) — PLAIG (crypto Modules 1–2) Galaxy suite design
- [`ancor-reference.md`](ancor-reference.md) — ANCOR consolidated reference
- [`panteon-reference.md`](panteon-reference.md) — PANTEON pangenome reference
- [`capheine-updates.md`](capheine-updates.md) — CAPHEINE required changes

Primary sources:
- Crypto protocol: `Cryptosporidium_Protocol_v2_1c_Evolutionary_Genetics.docx`
- ANCOR designs: `ancor-workflow-designs.md`, `ANCOR_workflows.md`, `ANCOR_workflows_detailed.md`, `ancor-core.ga`
- Pangenome: `brc-tools/libs/pangenome-helpers/`
- CAPHEINE: `galaxyproject/iwc` → `workflows/comparative_genomics/hyphy/`; Nextflow sibling `veg/CAPHEINE`
- MalariaGEN Plasmodium proxies: [github.com/asgiraldoc/MalariaGEN-Pf8-Data-Retrieval-Proxy](https://github.com/asgiraldoc/MalariaGEN-Pf8-Data-Retrieval-Proxy) — reviewed 2026-08-05; deployed as 6 per-species instances (Pf8, Pv5, Pk1, Pmb1, Poc1, Pow1); gap analysis in [`plaig-design.md`](plaig-design.md), BRC integration notes in §11
