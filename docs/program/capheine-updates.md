# CAPHEINE Updates — Reference Doc

*Companion to `comparative-genomics-program.md` (§3, §3.1). What CAPHEINE is today, and the specific changes required for it to serve as the program's single selection-analysis engine.*

Status: in IWC (maintained by this group); updates at design stage, 2026-07-30

## What CAPHEINE is today

**CAPHEINE — Comprehensive Automated Pipeline using HyPhy for Evolutionary Inference with NExtflow** — Galaxy implementation of the HyPhy comparative-genomics pipeline, in IWC at `workflows/comparative_genomics/hyphy/`. Nextflow sibling: `veg/CAPHEINE`. Currently **viral-tuned** (authored for FASTQ viral genomes, HIV-1 test data).

| Workflow | Contents |
|---|---|
| `hyphy-preprocessing` | reference-anchored codon alignment: gffread → remove terminal stops → **cawlign** → cleanup → **hyphy-cln** → IQ-TREE |
| `hyphy-core` | **MEME, PRIME, BUSTED, FEL** |
| `hyphy-compare` | **CFEL, RELAX** |
| `capheine-core-and-compare` | orchestrator over the three subworkflows |

Inputs today: reference GTF + FASTA + per-sample unaligned FASTAs. HyPhy pinned 2.5.96.

## Why it must change (program drivers)

1. **Crypto protocol Module 3 (§8.5)** mandates the full HyPhy suite with mandatory GARD — CAPHEINE covers roughly half the method list and has no GARD step.
2. **Pangenome WF-H reroute** — PANTEON's per-OG BUSTED should become a CAPHEINE invocation (program doc §2.1 item 1), but per-OG pangenome inputs have no single reference, so cawlign-only preprocessing can't serve them.
3. **Eukaryotic scope** — internal stop codons, multi-exon CDS, frameshifts: current preprocessing assumes clean viral CDS.

## Required changes

1. **GARD pre-check** — GARD on every alignment, partition-specific topologies for recombining alignments (crypto §8.5 + Appendix B; CAPHEINE's own README already warns recombination misleads estimates). Add to `hyphy-preprocessing` or as an optional profile. `hyphy_gard` already IUC-wrapped. **Must be post-alignment and common to all aligner paths** — not bypassable by aligner choice.
2. **Method coverage** — crypto §8.5.1 requires: SLAC, BUSTED (SRV enabled), aBSREL, FUBAR (+ FEL follow-up), MEME, RELAX, Contrast-FEL. Have: MEME/PRIME/BUSTED/FEL + CFEL/RELAX. **Add: aBSREL, FUBAR, SLAC** (all IUC-wrapped). Also confirm the BUSTED wrapper exposes `--srv` (program doc §2.1: pinned 2.5.96+galaxy0 does not — bump or patch).
3. **Eukaryotic input handling** — strict alignment QC gates for internal stops / multi-exon CDS; fail loudly, never silently estimate.
4. **Pluggable preprocessing (aligner choice)** — see design below.
5. **Ledger-based input adapter** — thin front-end subworkflow: orthology ledger (per-OG sequence sets from PANTEON / crypto Module 3) → per-OG sequence collections → CAPHEINE core inputs, via the normalization layer.
6. **Branch-set parameterization** — `hyphy-compare` foreground labeling (regex/list) already exists; publish crypto §8.5.2 contrasts (human-infecting vs restricted, anthroponotic vs zoonotic, intestinal vs gastric, ruminant-associated) as **named preset lists** — documentation, no new workflow.

## Pluggable preprocessing design

Uniform contract: every path takes normalized input, emits codon-alignment collection (FASTA) + gene-tree collection (nwk) + QC report; `hyphy-core`/`hyphy-compare`/GARD ingest unchanged. Cleanup (hyphy-cln), IQ-TREE, and GARD stay post-alignment and shared.

| Aligner | Input regime | Caveat |
|---|---|---|
| cawlign | reference-anchored (many samples vs one reference gene set; current viral case) | requires a reference |
| MACSE | de novo per-OG codon MSA; tolerates internal stops/frameshifts — **eukaryotic/pangenome default** | slower; **not yet wrapped** (needs IUC wrapper) |
| PRANK | crypto protocol's named alternative to MACSE for the same de novo regime (§7.2.1: "MACSE or PRANK in codon mode") — a second tool under the MACSE branch, not a fourth path | not yet wrapped |
| MAFFT | clean conserved in-frame CDS; fastest | not codon-aware: aa-align + back-translate only (pal2nal), or gate to verified in-frame inputs; frame-breaks silently. Note: PANTEON WF-F already uses exactly this sanctioned path (MAFFT L-INS-I + pal2nal + trimAL) |

Normalization front-end maps input classes to the contract: (a) GTF+genome+per-sample FASTAs (current gffread path); (b) per-OG CDS collections from the ledger; (c) protein-only (back-translate); (d) pre-aligned (skip alignment — crypto Module 2 case). QC severity is aligner-dependent: internal stops tolerated on MACSE/cawlign, hard-filtered on MAFFT.

### Open items before building

- **VERIFY**: does `hyphy-cln`'s filtering hold for MACSE output containing documented frameshifts/internal stops? Most likely silent-sequence-drop seam. Test with a frameshift-containing OG, confirm retention; decide whether hyphy-cln needs a codon-tolerant mode or frameshift OGs bypass it (with justification). (Program doc §10 open question.)
- **Galaxy mechanics**: Option A — single workflow, `aligner` parameter + `when`-conditional branches + expression-tool output coalescing (less trodden in IWC review); Option B — three aligner-specific preprocessing subworkflows + thin dispatch parent (matches existing modular style). **Default: Option B** unless conditional execution is confirmed acceptable to IWC reviewers.

## Tool availability

- Already wrapped (IUC): full HyPhy suite — gard, slac, fel, meme, fubar, busted, absrel, relax, cfel (+ cawlign, hyphy_cln, gffread, IQ-TREE as preprocessing deps).
- **Needs new wrappers**: MACSE and/or PRANK (de novo codon alignment). Wrapper ownership: PLAIG (program doc §6).

## Sequencing

1. GARD + method-coverage + euk-QC PR to CAPHEINE (independent of aligner work).
2. MACSE/PRANK wrapper → pluggable preprocessing (Option B) → VERIFY hyphy-cln seam.
3. Ledger input adapter (after PANTEON WF-E ledger sidecars exist).
4. **Only then** reroute PANTEON WF-H as a thin CAPHEINE parent — do not churn WF-H twice (program doc §2.1).

Result: pangenome per-OG BUSTED and crypto Module 3's full suite both execute as CAPHEINE invocations with different inputs/profiles — one workflow family to maintain, one to defend in IWC review.

## Sources

- IWC: `galaxyproject/iwc` → `workflows/comparative_genomics/hyphy/` (README, `hyphy-preprocessing.ga` step list)
- Crypto protocol v2.1c §7.2.1 (alignment), §8.5/8.5.1/8.5.2 (HyPhy suite + contrasts), Appendix B (GARD failure mode)
- Program doc §3/§3.1 (decisions), §2.1 item 1 (WF-H reroute), §10 (open questions)
