# BRC Pathogen Comparative Genomics Program

One umbrella program, four suites at three evolutionary timescales — from read-level population genetics to cross-species comparative genomics — coordinated as a single effort for Galaxy/IWC delivery on BRC Analytics.

## The suites

Three suites at three evolutionary timescales, plus a shared engine underneath:

| Suite | Full name | Timescale | Status |
|---|---|---|---|
| **PLAIG** | Pathogen Lineage Analysis & Inference of Genotypes | recent (read-level) | design — from crypto protocol Modules 1–2 |
| **PANTEON** | PANgenome Trees & Evolutionary Orthology Networks | within-species | working prototype (the only timescale suite with one) |
| **ANCOR** | Analysis of Networked Core Orthologs for Reconstruction | cross-species | design stage |
| **CAPHEINE** | Comprehensive Automated Pipeline using HyPhy for Evolutionary Inference with NExtflow | any | in IWC, maintained by this group |

CAPHEINE is the outlier — it's not a timescale suite but the **shared selection-analysis engine** that sits underneath both PANTEON and ANCOR. PANTEON's per-orthogroup BUSTED and ANCOR's cross-species HyPhy work both route through it.

> Suite names are working titles — tentative, changeable, but better than generic letters.

## Start here

The **[Coordination Plan](program/comparative-genomics-program.md)** translates assets the group already knows (the crypto protocol, the pangenome workflows, the HyPhy suite) into the proposed suite structure. It's the entry point — companion docs linked from it contain the detailed designs.

## Documents

- [Coordination Plan](program/comparative-genomics-program.md) — the main summary and decision record
- [PLAIG Design](program/plaig-design.md) — read-level popgen suite (crypto Modules 1–2)
- [PANTEON Reference](program/panteon-reference.md) — within-species pangenome pipeline
- [ANCOR Reference](program/ancor-reference.md) — cross-species orthology/synteny/gene-family suite
- [CAPHEINE Updates](program/capheine-updates.md) — HyPhy suite extension spec
