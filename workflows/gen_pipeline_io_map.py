#!/usr/bin/env python3
"""Generate workflows/pipeline_io_map.html — a step-level flow diagram of the
workflows vetted in the Pv4 clean re-run (currently A, B, C, C2), plus the port
detail behind it.

Each workflow is drawn as its own DAG flowing top to bottom: workflow inputs on
the first row, tool steps ranked by dependency depth, workflow outputs on the
last. Clusters stack in pipeline order and inter-workflow edges run down a
right-hand lane from one workflow's output port into the next one's input port,
so the whole vetted stretch reads as a single graph and grows downward -- not
sideways -- as phases are added.

Where gen_pipeline_dataflow.py is the A->K overview (workflow-level boxes only),
this goes one level down: the actual steps, and what each port IS.

The graph, ports and docs are parsed live from the workflow files; the wiring is
imported from gen_pipeline_dataflow.EDGES, so the two documents cannot disagree.
Only what can't be derived from a workflow file is curated here: external input
provenance, observed run evidence, and per-port shapes.

Add a workflow to VETTED as it passes the clean re-run.

Run:  python workflows/gen_pipeline_io_map.py   (needs PyYAML)
"""
import html
import json
import os
import re
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_pipeline_dataflow import EDGES, WF_META, VETTED  # single source of truth

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Prose lives here, not in this file -- see load_descriptions().
DESCRIPTIONS = "workflows/workflow_descriptions.md"
# Per-step example data, frozen from real invocations by workflows/capture_examples.py.
# The Pages build has no Galaxy access, so these must be committed, not fetched.
EXAMPLES_DIR = "workflows/examples"

# ---- palette -------------------------------------------------------------
# Taken from BRC Analytics: site-config/brc-analytics/local/theme/constants.ts
# (primary) over the findable-ui base palette it builds on
# (DataBiosphere/findable-ui, src/theme/common/palette.ts). Light theme.
BRC = {
    "primary": "#28285B", "primaryDark": "#1F1F47",
    "info": "#00729C", "infoLight": "#97D6EA", "infoLightest": "#F2FAFC",
    "success": "#287555", "successLight": "#AEE9D1", "successLightest": "#F1F8F5",
    "warning": "#B54708", "warningLight": "#FFD79D", "warningLightest": "#FFFAEB",
    "alert": "#B42318", "alertLightest": "#FFF4F4",
    "caution": "#956F00", "cautionLight": "#FFEB78",
    "ink": "#212B36", "inkLight": "#637381",
    "smokeDark": "#C4CDD5", "smoke": "#E1E3E5",
    "smokeLight": "#F6F6F7", "smokeLightest": "#FAFBFB",
    "bg": "#F6F6F7", "white": "#FFFFFF",
}

# Per-workflow accent, drawn from the same palette rather than the brighter set
# gen_pipeline_dataflow.py uses for the A->K overview. success/alert are kept
# out of it: they mean "verified" and "gap" everywhere else in this page.
WFCOL = {"A": BRC["primary"], "B": BRC["info"], "C": BRC["warning"],
         "C2": BRC["caution"], "E": BRC["ink"], "I": BRC["success"]}

# Workflows that have passed the clean re-run, in pipeline order. They are
# stacked top-to-bottom in this order and each one's own graph also flows
# downward, so the canvas grows vertically as phases are added rather than
# sideways — the page scrolls the way a long pipeline actually reads.
# Producers must precede consumers so every cross-workflow edge points down.



EXTERNAL = {
    ("A", "assemblies"): "Staged panel genomes — $PV4_SSD/pv4_full/inputs/assemblies/{strain}.fa",
    ("A", "proteomes"): "gffread-derived protein FASTAs; PvP01 + Sal-I use PlasmoDB curated sets",
    ("A", "busco_lineage"): "Literal string — apicomplexa_odb10 for Pv4",
    ("B", "assemblies"): "Same staged panel genomes as WF-A",
    ("C2", "anchor_assemblies"): "Staged genomes, restricted to the 3 curated anchors",
    ("C2", "anchor_gene_gff3s"): "PlasmoDB-68 gene GFF3 for the 3 anchors",
    ("C2", "anchor_bed12s"): "PlasmoDB-68 BED12 for the 3 anchors (triage/merge reference bed)",
    ("C2", "assemblies"): "Staged panel genomes, unmasked (Liftoff target)",
    ("C2", "anchor_isoforms"): "anchor_prep output — gene<TAB>transcript per anchor (TOGA2 --isoform_file)",
}

SHAPE = {
    ("A", "in", "assemblies"): "list[8] · id=strain",
    ("A", "in", "proteomes"): "list[8] · id=strain",
    ("A", "in", "busco_lineage"): "string",
    ("A", "out", "similarity_matrix"): "1 CSV · 8x8",
    ("A", "out", "signatures"): "list[8]",
    ("A", "out", "busco_summaries"): "list[8]",
    ("A", "out", "sourmash_heatmap"): "1 PNG",
    ("A", "out", "sourmash_dendrogram"): "1 PNG",
    ("A", "out", "qc_report"): "1 HTML",
    ("A", "out", "sizes"): "list[8] · id=strain",
    ("A", "out", "self_pairs"): "1 txt · 8 rows",
    ("A", "out", "relabel_map"): "1 tabular · 64 rows",
    ("B", "in", "assemblies"): "list[8] · id=strain",
    ("B", "out", "softmasked_fasta"): "list[8] · id=strain",
    ("B", "out", "fasta_index"): "list[8]",
    ("B", "out", "dustmasker_bed"): "list[8] · BED6",
    ("B", "out", "windowmasker_bed"): "list[8] · BED6",
    ("B", "out", "tantan_bed"): "list[8] · BED6",
    ("B", "out", "fastan_bed"): "list[8] · BED6",
    ("B", "out", "masking_table"): "1 tabular · 8x4",
    ("B", "out", "masking_report"): "1 HTML",
    ("C", "in", "masked_fastas"): "list[8] · id=strain",
    ("C", "in", "sizes"): "list[8] · id=strain",
    ("C", "in", "self_pairs"): "1 txt · 8 rows",
    ("C", "in", "relabel_map"): "1 tabular · 64 rows",
    ("C", "out", "cleaned_chains"): "list[56] · id=A.B",
    ("C", "out", "rbest_chains"): "list[56] · id=A.B",
    ("C", "out", "pairwise_axt"): "list[56] · id=A_B",
    ("C2", "in", "anchor_assemblies"): "list[3] · id=anchor",
    ("C2", "in", "anchor_gene_gff3s"): "list[3] · id=anchor",
    ("C2", "in", "anchor_bed12s"): "list[3] · id=anchor",
    ("C2", "in", "assemblies"): "list[8] · id=strain",
    ("C2", "in", "query_masked"): "list[8] · id=strain",
    ("C2", "in", "anchor_masked"): "list[3] · id=anchor",
    ("C2", "in", "anchor_isoforms"): "list[3] · id=anchor",
    ("C2", "in", "cleaned_chains"): "list[56] · id=target.query",
    ("C2", "out", "merged_annotations"): "list[21] · id=anchor_query",
    ("C2", "out", "classifications"): "list[21] · id=anchor_query",
}

FALLBACK_DOC = {
    ("C", "out", "cleaned_chains"): "Per-pair chainNet/netChainSubset chains, relabelled to Phase E ids (A.B).",
    ("C", "out", "rbest_chains"): "Reciprocal-best chains per ordered pair, relabelled to A.B — the Phase E orthology signal.",
    ("C", "out", "pairwise_axt"): "Raw pairwise alignments (axt) off KegAlign/batched LASTZ; feeds the WF-I multiz fold.",
    ("C2", "out", "merged_annotations"): "Per-pair projected annotation GFF3 (Liftoff, merged against the anchor BED12 reference).",
    ("C2", "out", "classifications"): "Per-pair classification TSV: reference_gene_id, query_gene_id, source, intactness, query coords, orthology_class.",
}

STATUS = {
    "A": dict(jobs="39 / 39 ok", when="2026-08-05", history="9a7c6ce629ee1eba",
              invocation="2f5c1a0fa13e5b9d",
              note="Workflow version 6 in the editor (the API reports version 5, "
                   "zero-indexed), 13 steps. This run is the first to include anchor_prep, "
                   "which derives the anchor BED12 and isoforms table in-workflow; its three "
                   "isoforms tables are identical to the files previously staged by hand."),
    "B": dict(jobs="122 / 122 ok", when="2026-06-17", history="f4af09719299d4f6",
              invocation="75ed0631e40c6c7b",
              note="Workflow version 11 in the editor (the API reports version 10, "
                   "zero-indexed), 18 steps. 8 / 8 genomes soft-masked after the rewrite to "
                   "four maintained maskers; sdust and longdust were dropped. This is the only "
                   "invocation whose steps match the committed workflow -- an earlier 72-job "
                   "run in history feb4114d4f866b50 was a 10-step version with neither the "
                   "coverage nor the masking-table steps."),
    "C": dict(jobs="848 / 848 ok", when="2026-06-18", history="11862b69df84527c",
              invocation="b38588d5a72103c1",
              note="Workflow version 7 in the editor (API version 6), 27 steps, matching the "
                   "committed file exactly. Full ordered grid: 56 pairs from 8 strains, "
                   "KegAlign on GPU plus batched LASTZ. Previously recorded as 500 jobs -- "
                   "that was /api/jobs returning its default page of 500."),
    "E": dict(jobs="12 / 12 ok", when="2026-08-05", history="6bba8e615983c09d",
              invocation="feb4114d4f866b50",
              note="The first WF-E run on TOGA2-era projections, and the first with the "
                   "consensus tool resolving projected genes onto native ones. 4,286 of 5,516 "
                   "orthogroups are CORE-1:1 (77.7%); the previous run, on Liftoff-only "
                   "evidence and without that resolution, reported 21 of 5,817 (0.4%). This "
                   "version no longer takes a pangenome graph."),
    "I": dict(jobs="64 / 64 ok", when="2026-08-05", history="cb183a1d519b63a2",
              invocation="e37f25bfb58457fb",
              note="The first WF-I run to produce output at all. The two earlier "
                   "invocations (2026-06-13) report state 'completed', which is misleading: "
                   "every one of their 8 MAFs was in error. Two bugs, both fixed here -- the "
                   "hinge name reached multiz_fold as the repr of a collection element rather "
                   "than a string, and the MAF src names carried bare contig accessions, which "
                   "multiz reads as one species per contig."),
    "C2": dict(jobs="15 / 21 cells", when="2026-08-05", history="5fadef28a25ac561",
               invocation="301fbf6783f66c6e",
               note="Workflow version 2 in the editor (API version 1), 29 steps. Full grid with "
                    "the TOGA2 rescue pass on: 15 of 21 cells complete -- PvW1 7/7, PvSY56 7/7, "
                    "PAM 1/7. The six PAM failures are an upstream TOGA2 v2.0.8 defect "
                    "(hillerlab/TOGA2#41), triggered by one transcript and reproducing on main. "
                    "Two complete anchors is what Phase E needs. Across the completed cells, "
                    "41,919 liftoff plus 86,400 cesar2 classifications; usable Phase E evidence "
                    "goes from 41,919 to 101,435 edges. The grid ran on the Aug 4 version, whose "
                    "steps are identical to version 2."),
}

CAVEATS = [
    ("WF-C and WF-C2 outputs carry no <code>doc:</code> strings",
     "All five output ports across the two Phase C workflows are undocumented in the "
     "<code>.gxwf.yml</code>. The descriptions here are reconstructed from the verified run and live "
     "in <code>FALLBACK_DOC</code> in this generator — they should move into the workflow files "
     "before IWC submission."),
    ("The TOGA2 container is not reproducible",
     "Upstream's Apptainer def installs Nextflow from <code>get.nextflow.io</code>, i.e. whatever is "
     "newest at <b>build</b> time. A build made today ships a Nextflow whose strict parser rejects "
     "TOGA2 2.0.8's own <code>execute_joblist.nf</code>, so rebuilding the image on a different day "
     "gives a different result. The wrapper works around it with "
     "<code>NXF_SYNTAX_PARSER=v1</code>, but that parser is deprecated upstream; the durable fix is "
     "to pin <code>NXF_VER</code> in the def and rebuild."),
    ("The TOGA2 pass has only been proven on one cell",
     "Pass 2 is verified end to end on PvW1&rarr;PvP01 (82 minutes). The full 21-cell grid has not "
     "been run with TOGA2 enabled. Expect it to be dominated by the PvSY56 anchor, whose genes fall "
     "through to CESAR2 at 75&ndash;83&nbsp;% against roughly 21&nbsp;% for PvW1 and PAM."),
    ("Phase C.2 / C.4 TSVs are written with CRLF line endings",
     "<code>phase_c2_triage.py</code> and <code>phase_c4_merge.py</code> pass <code>newline=''</code> "
     "but leave <code>csv</code>'s default <code>lineterminator='\\r\\n'</code> in place. Harmless for "
     "<code>phase_e_consensus</code> (universal newlines) but wrong for a published artifact: the last "
     "column, <code>orthology_class</code>, carries a trailing CR for anything splitting on <code>\\n</code>."),
    ("<code>phase_c2_triage</code> writes an empty <code>needs_cesar2.bed</code>",
     "It looks flagged genes up in the anchor BED12 by gene id, but that file is transcript-keyed "
     "(<code>PVW1_000005000_t1</code>) as TOGA2 requires — zero overlap, so the BED is always empty. "
     "Harmless today: the C.4 redesign has TOGA2 classify the full anchor annotation rather than a "
     "triage subset, so nothing consumes it. Worth removing or fixing when that tool is next touched."),
    ("WF-B accepts any list collection",
     "<code>assemblies</code> is typed as a bare <code>collection_type: list</code>, so the run form "
     "accepts a non-FASTA collection and fails at the first tool rather than at selection time."),
]

# ---- the "choosing your inputs" guide ------------------------------------
# The Pv4 panel as staged, measured off inputs/assemblies/*.fa. P. vivax has 14
# chromosomes, so sequence COUNT is what separates a chromosome-level assembly
# from a scaffold-level one here -- N50 does not (every P. vivax assembly picks
# up near-chromosome scaffolds, so they all land at 1.5-2.1 Mb).
PANEL = [
    ("PvSY56", 14, 23.8, 2.09, "anchor"),
    ("PvW1", 19, 29.0, 2.12, "anchor"),
    ("PAM", 28, 29.4, 2.00, "anchor"),
    ("MHC087", 126, 29.2, 2.11, "query"),
    ("PvP01", 242, 29.0, 1.76, "query"),
    ("PvT01", 374, 29.0, 1.56, "query"),
    ("PvC01", 425, 30.2, 1.59, "query"),
    ("Sal-I", 2747, 27.0, 1.68, "query"),
]

ANCHOR_RULES = [
    ("Curated annotation that matches the assembly you staged",
     "The hard gate. The GFF3 seqids must match the FASTA headers exactly — an annotation built "
     "against a different assembly build of the same strain is unusable without chromosome "
     "reconciliation."),
    ("Contiguity",
     "Liftoff projects genes across a chain, so a fragmented donor breaks gene models at scaffold "
     "boundaries. Prefer chromosome-level assemblies."),
    ("Cost is linear in anchors",
     "The projection grid is anchors × strains, so each extra anchor adds N−1 Liftoff → triage → "
     "merge chains."),
    ("Two or three anchors buy independent evidence",
     "Phase E builds its consensus orthology graph from every anchor's classifications, so a gene "
     "supported by two anchors is stronger than one. One anchor works; three cross-validate."),
]

# How to choose / produce each externally staged input. Keyed like EXTERNAL.
INPUT_GUIDE = {
    ("A", "assemblies"): (
        "Every genome in the panel, unmasked. This defines the panel: everything downstream is "
        "keyed by these element identifiers, so pick the strain names carefully — they become "
        "collection ids, pair ids (<code>A_B</code>), and hub track names.",
        "8 P. vivax genomes fetched from NCBI by the accessions in <code>species.conf</code>. "
        "Element identifier = strain name (<code>PvP01</code>, <code>Sal-I</code>, …)."),
    ("A", "proteomes"): (
        "One protein FASTA per strain, parallel to <code>assemblies</code>. BUSCO runs in "
        "<code>-m prot</code> mode. Normally <code>gffread -y</code> from each strain's own "
        "annotation.",
        "Generated with gffread from the native GFFs — except PvP01 and Sal-I, which were swapped "
        "to the PlasmoDB curated sets because their NCBI assembly and annotation coordinates "
        "disagree."),
    ("A", "busco_lineage"): (
        "The smallest BUSCO clade that still contains your species. Too broad and completeness is "
        "meaningless; too narrow and the dataset may not exist. Must be provisioned by the data "
        "manager first — it does not populate the dropdown by clade alone.",
        "<code>apicomplexa_odb10</code> (n=446). Observed: PvP01/PvT01/Sal-I 100 %, PvC01/PvW1 "
        "99.6 %, PAM 99.3 %, PvSY56 89.9 %, MHC087 87.0 %."),
    ("B", "assemblies"): (
        "The same panel genomes as WF-A — unmasked. WF-B soft-masks them; every later phase "
        "consumes the masked output rather than these.",
        "Identical collection to WF-A's."),
    ("C2", "anchor_assemblies"): (
        "The subset of <code>assemblies</code> you have curated annotation for — the gene donors. "
        "See the anchor rules above.",
        "PvW1, PAM, PvSY56 — the only chromosome-level assemblies in the panel."),
    ("C2", "anchor_gene_gff3s"): (
        "The curated gene GFF3 for each anchor. Two hard requirements: seqids must match the "
        "anchor FASTA headers, and gene-level feature types must be renamed to <code>gene</code> — "
        "Liftoff's default mode finds zero <code>gene</code> features in a native "
        "<code>protein_coding_gene</code> / <code>ncRNA_gene</code> / <code>pseudogene</code> GFF3 "
        "and silently projects nothing.",
        "PlasmoDB-68 <code>*.gene.gff3</code> for the three anchors."),
    ("C2", "anchor_bed12s"): (
        "BED12 transcript models per anchor, used by triage and merge as the reference bed. "
        "Derive it from the same GFF3 you use above so the two can never drift apart.",
        "<code>gff3ToGenePred anchor.gff3 stdout | genePredToBed stdin anchor.bed12</code> "
        "(UCSC kent), or PlasmoDB's own BED."),
    ("C2", "anchor_masked"): (
        "The soft-masked version of each anchor genome. TOGA2 wants soft-masked sequence on "
        "<b>both</b> sides, so the anchor side is needed here even though Liftoff aligns against "
        "the unmasked copy. Same WF-B output as <code>query_masked</code>, restricted to anchors.",
        "PvW1, PAM, PvSY56 from WF-B's <code>softmasked_fasta</code>."),
    ("C2", "anchor_isoforms"): (
        "A two-column gene-to-transcript table per anchor, grouping the BED12 transcripts into "
        "genes. The transcript ids must match column 4 of <code>anchor_bed12s</code> exactly, or "
        "TOGA2 cannot tie projections back to genes.",
        "<code>anchor_prep</code> output; 6,075 transcripts for PvW1, matching its BED12 one for one."),
    ("C2", "cleaned_chains"): (
        "WF-C's cleaned chains for every ordered strain pair. Only the TOGA2 pass needs these; "
        "Liftoff does its own alignment. The workflow selects the anchor-to-query subset and "
        "renames it to match the grid, so pass the whole collection.",
        "56 chains for 8 strains, filtered in-workflow to the 21 anchor cells."),
    ("C2", "assemblies"): (
        "The query axis of the projection grid: the genomes receiving genes. Pass the whole panel — "
        "the anchor self-cells are dropped internally, so anchors receive each other's genes "
        "without projecting onto themselves.",
        "All 8 strains → 3 × 8 = 24 grid cells, filtered to 21."),
}

DERIVED_NOTE = (
    "Three inputs that used to be hand-authored are now produced in-workflow and need no staging: "
    "<code>sizes</code>, <code>self_pairs</code> and <code>relabel_map</code> come out of WF-A and "
    "wire straight into WF-C, and WF-C2 generates its own anchor self-pair list from the anchor "
    "element identifiers. The old <code>execution/cluster/gen_wfc_config.py</code> is only needed "
    "if you drive WF-C standalone without running WF-A first.")

# ---- layout constants ----------------------------------------------------
# Flow is top-to-bottom: a dependency rank is a ROW, and the nodes sharing a rank
# spread across it. GY is therefore the rank-to-rank (vertical) gap.
NW, NH = 152, 36           # node box
GX, GY = 20, 40            # within-rank (x) and rank-to-rank (y) gaps
PADX, PADY = 20, 42        # cluster padding (top pad leaves room for the header)
CLUST_GAP_Y = 74           # vertical space between workflow clusters
MAX_ROW = 5                # widest a rank may be drawn; longer ranks wrap onto
                           # extra rows so the canvas never needs sideways scrolling
GUTTER = 46                # right-hand lane where cross-workflow edges run
LANE = 15                  # per-edge offset inside the gutter


def esc(s):
    return html.escape(str(s or ""))


def md_inline(s):
    """Escape, then apply the inline markdown the descriptions file allows."""
    s = html.escape(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<i>\1</i>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    return s


def md_block(text, tag="p"):
    """Blank-line-separated paragraphs -> <p>...</p>.

    A paragraph whose lines are all indented four spaces becomes a <pre> instead,
    with its whitespace preserved -- these carry column-aligned samples (FASTA
    headers, TSV rows) where the alignment is the point.
    """
    out = []
    for para in re.split(r"\n\s*\n", text.strip()):
        if not para.strip():
            continue
        lines = para.split("\n")
        if all(l.startswith("    ") or not l.strip() for l in lines):
            body = "\n".join(l[4:] for l in lines).rstrip()
            out.append("<pre>" + html.escape(body) + "</pre>")
        else:
            out.append(f"<{tag}>{md_inline(' '.join(para.split()))}</{tag}>")
    return "".join(out)


def load_examples():
    """{wf_id: captured-run dict} for every workflow with a committed example file."""
    out = {}
    d = os.path.join(ROOT, EXAMPLES_DIR)
    if not os.path.isdir(d):
        return out
    for fn in sorted(os.listdir(d)):
        if not fn.startswith("wf_") or not fn.endswith("_examples.json"):
            continue
        with open(os.path.join(d, fn)) as fh:
            data = json.load(fh)
        out[data["workflow"]] = data
    return out


def parse_columns(block):
    """`name -- meaning` lines -> [(name, meaning)]; a tabular sample is a wall of
    text without them, and a header row alone does not say what a column means."""
    cols = []
    for line in (block or "").split("\n"):
        line = line.strip().lstrip("-").strip()
        if not line:
            continue
        m = re.match(r"^`?([A-Za-z0-9_.<>{} -]+?)`?\s+(?:--|\u2014)\s+(.+)$", line)
        if m:
            cols.append((m.group(1).strip(), m.group(2).strip()))
        elif cols:
            cols[-1] = (cols[-1][0], cols[-1][1] + " " + line)     # continuation
    return cols


def load_descriptions(path=DESCRIPTIONS):
    """Parse workflow_descriptions.md -> {wf_id: {"summary": .., "description": ..}}.

    The prose lives in markdown rather than in this file so it can be edited
    without touching Python. Missing or empty sections are a hard error: a
    silently blank description would publish to the live site unnoticed.
    """
    raw = open(os.path.join(ROOT, path)).read()
    raw = re.sub(r"<!--.*?-->", "", raw, flags=re.S)      # strip comments
    out, wid, field = {}, None, None
    for line in raw.split("\n"):
        m2, m3 = re.match(r"##\s+(\S+)\s*$", line), re.match(r"###\s+(.+?)\s*$", line)
        if m2 and not m3:
            wid, field = m2.group(1), None
            out.setdefault(wid, {})
        elif m3 and wid:
            # one heading may name several nodes: "### columns:cls_ok, cls_nested"
            head = m3.group(1).lower()
            if ":" in head:
                kind, names = head.split(":", 1)
                field = [f"{kind}:{n.strip()}" for n in names.split(",") if n.strip()]
            else:
                field = [head]                   # "summary", "description"
            for f in field:
                out[wid][f] = ""
        elif wid and field is not None:
            for f in field:
                out[wid][f] += line + "\n"

    missing = []
    for w in VETTED:
        for f in ("summary", "description"):
            if not (out.get(w, {}).get(f) or "").strip():
                missing.append(f"{w}.{f}")
    if missing:
        raise SystemExit(f"{path}: missing or empty section(s): {', '.join(missing)}\n"
                         f"every workflow in VETTED needs '### summary' and '### description'.")
    return out


def short_tool(tid):
    if not tid:
        return ""
    if tid.startswith("__") and tid.endswith("__"):
        return tid.strip("_").lower().replace("_", " ")
    if "/repos/" in tid:
        tid = tid.split("/")[-2]
    return tid


def parse_wf(wid, path):
    """Parse a .gxwf.yml into nodes + intra-workflow edges."""
    d = yaml.safe_load(open(os.path.join(ROOT, path)))
    inputs = d.get("inputs") or {}
    steps = d.get("steps") or {}
    outputs = d.get("outputs") or {}

    nodes, edges = {}, []

    def nid(kind, name):
        return f"{wid}|{kind}|{name}"

    for name, v in inputs.items():
        v = v or {}
        kind = f"collection[{v['collection_type']}]" if v.get("collection_type") else v.get("type", "data")
        nodes[nid("in", name)] = dict(wf=wid, kind="in", name=name, sub=kind,
                                      doc=(v.get("doc") or "").strip())
    for name, s in steps.items():
        s = s or {}
        nodes[nid("st", name)] = dict(wf=wid, kind="st", name=name,
                                      sub=short_tool(s.get("tool_id")), doc="")
    for name, o in outputs.items():
        o = o or {}
        nodes[nid("out", name)] = dict(wf=wid, kind="out", name=name,
                                       sub=SHAPE.get((wid, "out", name), ""),
                                       doc=(o.get("doc") or "").strip())

    def resolve(src):
        """'step/out' or 'input_label' -> node id (or None)."""
        if isinstance(src, dict):
            src = src.get("source", "")
        if not isinstance(src, str) or not src:
            return None
        head = src.split("/", 1)[0]
        if head in steps:
            return nid("st", head)
        if head in inputs:
            return nid("in", head)
        if src in inputs:
            return nid("in", src)
        return None

    for name, s in steps.items():
        for _port, src in ((s or {}).get("in") or {}).items():
            a = resolve(src)
            if a:
                edges.append((a, nid("st", name)))
    for name, o in outputs.items():
        a = resolve((o or {}).get("outputSource"))
        if a:
            edges.append((a, nid("out", name)))

    return nodes, edges


def layout(nodes, edges):
    """Longest-path ranking + barycenter ordering. Mutates nodes with x/y/rank."""
    preds, succs = {n: [] for n in nodes}, {n: [] for n in nodes}
    for a, b in edges:
        if a in nodes and b in nodes:
            succs[a].append(b)
            preds[b].append(a)

    rank, seen = {}, set()

    def rk(n):
        if n in rank:
            return rank[n]
        if n in seen:              # defensive: a cycle would otherwise recurse forever
            return 0
        seen.add(n)
        r = 0 if not preds[n] else 1 + max(rk(p) for p in preds[n])
        rank[n] = r
        return r

    for n in nodes:
        rk(n)
    # workflow inputs always sit at rank 0; outputs always at the far right
    for n, d in nodes.items():
        if d["kind"] == "in":
            rank[n] = 0
    maxr = max(rank.values()) if rank else 0
    for n, d in nodes.items():
        if d["kind"] == "out":
            rank[n] = maxr

    order = {}
    for n in nodes:
        order.setdefault(rank[n], []).append(n)

    # barycenter sweeps to cut edge crossings
    for it in range(6):
        keys = sorted(order) if it % 2 == 0 else sorted(order, reverse=True)
        for r in keys:
            ref = order.get(r - 1 if it % 2 == 0 else r + 1, [])
            pos = {n: i for i, n in enumerate(ref)}
            rel = preds if it % 2 == 0 else succs
            cur = {n: i for i, n in enumerate(order[r])}

            def bary(n):
                vals = [pos[m] for m in rel[n] if m in pos]
                return sum(vals) / len(vals) if vals else cur[n]

            order[r] = sorted(order[r], key=lambda n: (bary(n), nodes[n]["name"]))

    # A rank wider than MAX_ROW is drawn as several stacked rows. Ranks stay
    # logical (edges are unaffected); only the drawing wraps, which keeps the
    # canvas narrow enough to read without scrolling sideways.
    rows = []                       # [(rank, [nodes]), ...] in draw order
    for r in sorted(order):
        row = order[r]
        for i in range(0, len(row), MAX_ROW):
            rows.append((r, row[i:i + MAX_ROW]))

    wide = max((len(chunk) for _, chunk in rows), default=1)
    for ri, (r, chunk) in enumerate(rows):
        off = (wide - len(chunk)) * (NW + GX) / 2
        for i, n in enumerate(chunk):
            nodes[n]["rank"] = r
            nodes[n]["x"] = PADX + off + i * (NW + GX)
            nodes[n]["y"] = PADY + ri * (NH + GY)

    w = PADX * 2 + wide * NW + (wide - 1) * GX
    h = PADY + PADY // 2 + len(rows) * NH + (len(rows) - 1) * GY
    return w, h




def figure_cross_product():
    """Explain __CROSS_PRODUCT_FLAT__: it emits two aligned collections.

    Drawn with a 3-genome panel so the whole 3x3 grid fits, and the diagonal that
    the self-pairs filter removes is marked.
    """
    G = ["PvW1", "PAM", "Sal-I"]
    cw, ch, gap = 74, 26, 6
    x0, y0 = 14, 40
    ROW = 52                      # rank pitch: box height plus room for the caption
    S = []
    W = x0 + 9 * (cw + gap) + x0
    H = 372
    S.append(f'<svg class="fig" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
             f'xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMin meet">')

    def box(x, y, w, label, cls=""):
        S.append(f'<g class="fb {cls}"><rect x="{x}" y="{y}" width="{w}" height="{ch}" rx="5"/>'
                 f'<text x="{x + w/2:.0f}" y="{y + 17}" text-anchor="middle">{esc(label)}</text></g>')

    def cap(x, y, t, anchor="start"):   # y is the baseline
        S.append(f'<text class="fcap" x="{x}" y="{y}" text-anchor="{anchor}">{esc(t)}</text>')

    # the two input collections
    cap(x0, y0 - 12, "one input collection, shown here with 3 genomes instead of 8")
    for i, g in enumerate(G):
        box(x0 + i * (cw + gap), y0, cw, g, "in")
    cap(x0 + 4 * (cw + gap), y0 - 12, "\u2026 crossed with itself")
    for i, g in enumerate(G):
        box(x0 + (4 + i) * (cw + gap), y0, cw, g, "in")

    # the flattened outputs: every ordered combination, as two aligned lists
    ya, yb, yid = 140, 140 + ROW, 140 + 2 * ROW
    cap(x0, ya - 12, "output_a  \u2192  tgt_fa      each genome repeated, 9 elements")
    cap(x0, yb - 12, "output_b  \u2192  qry_fa      all genomes cycled, 9 elements")
    cap(x0, yid - 12, "identifier \u2014 the SAME on both collections; it names the cell, "
                       "not the genome in it")
    k = 0
    for a in G:
        for b in G:
            x = x0 + k * (cw + gap)
            diag = "diag" if a == b else ""
            box(x, ya, cw, a, "out " + diag)
            box(x, yb, cw, b, "out " + diag)
            box(x, yid, cw, f"{a}_{b}", "id " + diag)
            k += 1

    # arrows from inputs down into the flattened rows
    for cx, ylow in ((x0 + 1.5*(cw+gap), ya), (x0 + 5.5*(cw+gap), ya)):
        S.append(f'<path class="farrow" d="M{cx:.0f},{y0+ch+4} L{cx:.0f},{ylow-20}"/>')

    ytext = yid + ch + 26
    cap(x0, ytext, "Read down a column: the third column is output_a=PvW1 against "
                   "output_b=Sal-I, identifier PvW1_Sal-I.")
    cap(x0, ytext + 19, "Shaded = a genome paired with itself. __FILTER_FROM_FILE__ drops those "
                        "from both collections, leaving 6 of 9 (56 of 64 for the real panel).")
    cap(x0, ytext + 38, "The pairing is carried by POSITION alone. Reorder either collection and "
                        "every cell silently gets the wrong partner.")
    S.append("</svg>")
    return "".join(S)


FIGURES = {"cross_product": figure_cross_product}


def expand_figures(html_text):
    """Replace {{figure:name}} tokens left by md_block with generated SVG."""
    def sub(m):
        fn = FIGURES.get(m.group(1))
        return fn() if fn else m.group(0)
    return re.sub(r"\{\{figure:([a-z_]+)\}\}", sub, html_text)


EX_KEYS = {}   # {wf_id: set(node names with an example)}; filled in main()


def solo_svg(wid, g, allnodes):
    """Render one workflow's DAG on its own, for its detail tab."""
    W = g["w"] + 2 * PADX
    H = g["h"] + 2 * PADY
    c = WFCOL[wid]
    S = [f'<svg class="solo" id="solo-{wid}" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'preserveAspectRatio="xMidYMin meet" '
         f'xmlns="http://www.w3.org/2000/svg">',
         f'<defs><marker id="sar{wid}" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" '
         f'markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="{c}"/></marker></defs>']

    def xy(n):
        d = g["nodes"][n]
        return d["x"] + PADX, d["y"] + PADY

    for a, b in g["edges"]:
        if a not in g["nodes"] or b not in g["nodes"]:
            continue
        x1, y1 = xy(a); x2, y2 = xy(b)
        x1 += NW / 2; y1 += NH
        x2 += NW / 2
        dy = max(16, (y2 - y1) * 0.5)
        S.append(f'<path class="e" data-a="{esc(a)}" data-b="{esc(b)}" stroke="{c}" '
                 f'marker-end="url(#sar{wid})" fill="none" '
                 f'd="M{x1:.0f},{y1:.0f} C{x1:.0f},{y1+dy:.0f} {x2:.0f},{y2-dy:.0f} {x2:.0f},{y2:.0f}"/>')

    for n, d in g["nodes"].items():
        x, y = xy(n)
        nm = d["name"] if len(d["name"]) <= 21 else d["name"][:20] + "\u2026"
        sub = d["sub"] if len(d["sub"]) <= 24 else d["sub"][:23] + "\u2026"
        has = ' has-ex' if d["name"] in EX_KEYS.get(wid, ()) else ''
        S.append(
            f'<g class="n {d["kind"]}{has}" data-wf="{wid}" data-node="{esc(d["name"])}" '
            f'data-id="{esc(n)}" data-kind="{d["kind"]}" style="--c:{c}" '
            f'transform="translate({x:.0f},{y:.0f})">'
            f'<title>{esc(d["name"])}</title>'
            f'<rect width="{NW}" height="{NH}" rx="7"/>'
            f'<text class="nn" x="9" y="14">{esc(nm)}</text>'
            f'<text class="ns" x="9" y="26">{esc(sub)}</text></g>')
    S.append("</svg>")
    return "".join(S)


def main():
    desc = load_descriptions()
    graphs = {}
    for wid in VETTED:
        path = WF_META[wid][0]
        nodes, edges = parse_wf(wid, path)
        w, h = layout(nodes, edges)
        graphs[wid] = dict(path=path, title=WF_META[wid][3], ph=WF_META[wid][1],
                           nodes=nodes, edges=edges, w=w, h=h)

    # stack clusters top-to-bottom in pipeline order, centred on a common axis
    body_w = max(g["w"] for g in graphs.values())
    y = 0
    for wid in VETTED:
        g = graphs[wid]
        g["ox"] = (body_w - g["w"]) / 2
        g["oy"] = y
        y += g["h"] + CLUST_GAP_Y

    CH = y - CLUST_GAP_Y

    # absolute node coords
    allnodes = {}
    for wid, g in graphs.items():
        for n, d in g["nodes"].items():
            d = dict(d)
            d["ax"], d["ay"] = d["x"] + g["ox"], d["y"] + g["oy"]
            allnodes[n] = d

    # cross-workflow edges, mapped onto real port nodes
    cross = []
    for f, fo, t, ti in EDGES:
        if f in VETTED and t in VETTED:
            a, b = f"{f}|out|{fo}", f"{t}|in|{ti}"
            if a in allnodes and b in allnodes:
                cross.append((a, b))

    # widen the canvas for the gutter the cross-workflow edges run down
    CW = body_w + GUTTER + LANE * max(1, len(cross))

    def path(a, b):
        """Intra-workflow edge: straight down, bottom of source to top of target."""
        x1, y1 = a["ax"] + NW / 2, a["ay"] + NH
        x2, y2 = b["ax"] + NW / 2, b["ay"]
        dy = max(16, (y2 - y1) * 0.5)
        return f"M{x1:.0f},{y1:.0f} C{x1:.0f},{y1 + dy:.0f} {x2:.0f},{y2 - dy:.0f} {x2:.0f},{y2:.0f}"

    def cross_path(a, b, lane):
        """Cross-workflow edge: out to a right-hand lane, down past the intervening
        clusters, then back in. Keeps these edges off the workflow boxes entirely."""
        x1, y1 = a["ax"] + NW / 2, a["ay"] + NH
        x2, y2 = b["ax"] + NW / 2, b["ay"]
        gx = body_w + GUTTER + lane * LANE
        return (f"M{x1:.0f},{y1:.0f} "
                f"C{x1:.0f},{y1 + 26:.0f} {gx:.0f},{y1 + 6:.0f} {gx:.0f},{y1 + 40:.0f} "
                f"L{gx:.0f},{y2 - 40:.0f} "
                f"C{gx:.0f},{y2 - 6:.0f} {x2:.0f},{y2 - 26:.0f} {x2:.0f},{y2:.0f}")

    S = []
    S.append(f'<svg id="g" viewBox="0 0 {CW} {CH}" width="{CW}" height="{CH}" '
             f'xmlns="http://www.w3.org/2000/svg">')
    S.append('<defs>')
    for wid in VETTED:
        S.append(f'<marker id="ar{wid}" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" '
                 f'markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="{WFCOL[wid]}"/></marker>')
    S.append('<marker id="argy" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" '
             f'orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="{BRC["inkLight"]}"/></marker></defs>')

    # cluster frames
    for wid, g in graphs.items():
        c = WFCOL[wid]
        S.append(f'<g class="clu" data-wf="{wid}">'
                 f'<rect x="{g["ox"]}" y="{g["oy"]}" width="{g["w"]}" height="{g["h"]}" rx="12" '
                 f'class="cbox" style="--c:{c}"/>'
                 f'<rect x="{g["ox"]}" y="{g["oy"]}" width="{g["w"]}" height="26" rx="12" fill="{c}"/>'
                 f'<rect x="{g["ox"]}" y="{g["oy"]+14}" width="{g["w"]}" height="12" fill="{c}"/>'
                 f'<text x="{g["ox"]+12}" y="{g["oy"]+18}" class="ctitle">{esc(g["title"])}</text>'
                 f'<text x="{g["ox"]+g["w"]-12}" y="{g["oy"]+18}" class="cmeta" text-anchor="end">'
                 f'{esc(STATUS[wid]["jobs"])}</text></g>')

    # intra edges, then cross edges on top
    for wid, g in graphs.items():
        for a, b in g["edges"]:
            if a in allnodes and b in allnodes:
                S.append(f'<path class="e" data-a="{a}" data-b="{b}" stroke="{WFCOL[wid]}" '
                         f'marker-end="url(#ar{wid})" d="{path(allnodes[a], allnodes[b])}"/>')
    for i, (a, b) in enumerate(cross):
        S.append(f'<path class="e x" data-a="{a}" data-b="{b}" stroke="{BRC["inkLight"]}" '
                 f'marker-end="url(#argy)" d="{cross_path(allnodes[a], allnodes[b], i)}"/>')

    # nodes
    for n, d in allnodes.items():
        c = WFCOL[d["wf"]]
        cls = f'n {d["kind"]}'
        nm = d["name"] if len(d["name"]) <= 21 else d["name"][:20] + "…"
        sub = d["sub"] if len(d["sub"]) <= 24 else d["sub"][:23] + "…"
        tip = d["doc"] or d["sub"] or d["name"]
        S.append(
            f'<g class="{cls}" id="{esc(n)}" data-wf="{d["wf"]}" style="--c:{c}" '
            f'transform="translate({d["ax"]:.0f},{d["ay"]:.0f})">'
            f'<title>{esc(d["name"])} — {esc(tip)}</title>'
            f'<rect width="{NW}" height="{NH}" rx="7"/>'
            f'<text class="nn" x="9" y="14">{esc(nm)}</text>'
            f'<text class="ns" x="9" y="26">{esc(sub)}</text></g>')
    S.append("</svg>")

    # ---- examples ---------------------------------------------------------
    examples = load_examples()
    # resolve which node each captured example belongs to, and pre-render the
    # payload the detail panel shows when that node is clicked
    ex_payload = {}
    stale = {}    # node id -> the invocation it is missing from (step added later)
    for wid, cap in examples.items():
        if wid not in graphs:
            continue
        EX_KEYS[wid] = set(cap["steps"])
        outsrc = {}   # workflow output name -> (step, output name)
        d = yaml.safe_load(open(os.path.join(ROOT, WF_META[wid][0])))
        # parameter inputs never carry a dataset, so their absence from a capture
        # says nothing; data/collection inputs are a different matter
        datain = {n for n, i in (d.get("inputs") or {}).items()
                  if str((i or {}).get("type", "data")) in ("data", "collection")}
        for oname, o in (d.get("outputs") or {}).items():
            src = str((o or {}).get("outputSource", ""))
            if "/" in src:
                st, on = src.split("/", 1)
                outsrc[oname] = (st, on)
                if st in cap["steps"]:
                    EX_KEYS[wid].add(oname)
        for node, dd in graphs[wid]["nodes"].items():
            name, kind = dd["name"], dd["kind"]
            step, want = (name, None)
            if kind == "out" and name in outsrc:
                step, want = outsrc[name]
            entry = cap["steps"].get(step)
            if not entry:
                if kind in ("st", "out") or (kind == "in" and name in datain):
                    # in the .gxwf.yml but absent from the capture: either the step
                    # postdates that run, or it produced nothing worth sampling.
                    # The panel states the fact rather than guessing which.
                    stale[node] = f'{cap["invocation"]} ({cap["when"]})'
                continue
            outs = [o for o in entry["outputs"] if want is None or o["name"] == want]
            if outs:
                # prose for this step, if workflow_descriptions.md carries one
                prose = (desc.get(wid, {}).get(f"step:{step}") or "").strip()
                dw = desc.get(wid, {})
                for o in outs:
                    common = desc.get("COMMON", {})
                    blk = (dw.get(f"columns:{step}.{o['name']}")
                           or dw.get(f"columns:{o['name']}")
                           or dw.get(f"columns:{step}")
                           or common.get(f"columns:{step}.{o['name']}")
                           or common.get(f"columns:{o['name']}")
                           or common.get(f"columns:{step}"))
                    cols = parse_columns(blk)
                    if cols:
                        o["columns"] = cols
                ex_payload[node] = {"step": step, "outputs": outs,
                                    "why": expand_figures(md_block(prose)) if prose else ""}

    # ---- HTML body ------------------------------------------------------
    P = []
    A = P.append
    A('<header><div><h1>Pv4 pangenome pipeline</h1>'
      '<div class="sub">How to choose the inputs, and what every step of the vetted workflows '
      'does with them. Worked throughout on the 8-strain <i>Plasmodium vivax</i> panel.</div></div>'
      f'<div class="score"><b>{len(VETTED)}</b> of 12 verified</div>'
      '<a class="backlink" href="pipeline_dataflow.html">&larr; pipeline overview</a>'
      '</header>')

    tabbed = [w for w in VETTED if w in examples]
    A('<nav><button class="tab active" data-tab="overview">Overview</button>'
      + "".join(f'<button class="tab" data-tab="wf-{w}" style="--c:{WFCOL[w]}">'
                f'WF-{w}</button>' for w in tabbed)
      + '</nav>')
    A('<div class="tabpane active" id="tab-overview">')

    # ---- getting started -------------------------------------------------
    A('<div class="wrap"><section id="start"><h2>Choosing your inputs</h2>'
      '<p class="note">Everything the pipeline needs is staged up front as a handful of Galaxy '
      'collections. Most are obvious — the genomes you want in the panel. The one decision that '
      'actually needs thought is which genomes serve as <b>anchors</b>.</p>')

    A('<blockquote><b>anchor vs query.</b> An <i>anchor</i> is a genome that already has a trusted, '
      'curated gene annotation; it is the <b>donor</b> of genes. A <i>query</i> is any genome you '
      'want those genes projected onto. WF-C aligns genomes to each anchor and WF-C2 lifts the '
      'anchor\'s genes across that alignment onto the queries. So: <b>anchors = the genomes you '
      'have annotation for; queries = everything you want annotated.</b></blockquote>'
      '<p class="note">Anchors are queries too. The grid is anchors &times; <i>all</i> strains with '
      'the self-cells removed, so anchors receive each other\'s genes.</p>')

    A('<h3>Picking anchors</h3><ol class="rules">')
    for t, b in ANCHOR_RULES:
        A(f'<li><b>{t}.</b> {b}</li>')
    A('</ol>')

    A('<h3>Worked example — the <i>P. vivax</i> panel</h3>'
      '<p class="note"><i>P. vivax</i> has 14 chromosomes, so sequence <b>count</b> is what '
      'separates a chromosome-level assembly from a scaffold-level one. N50 does not discriminate '
      'here — every assembly in the panel picks up near-chromosome scaffolds and lands at '
      '1.5&ndash;2.1 Mb.</p>'
      '<table class="panel"><thead><tr><th>strain</th><th>sequences</th><th>total</th><th>N50</th>'
      '<th>role</th></tr></thead><tbody>')
    for s, n, tot, n50, role in PANEL:
        cls = ' class="anc"' if role == "anchor" else ""
        tag = (f'<span class="pill" style="--c:{WFCOL["C2"]}">anchor</span>' if role == "anchor"
               else '<span class="mut">query</span>')
        A(f'<tr{cls}><td><code>{esc(s)}</code></td><td><b>{n:,}</b></td><td class="mut">{tot} Mb</td>'
          f'<td class="mut">{n50} Mb</td><td>{tag}</td></tr>')
    A('</tbody></table>'
      '<p class="note">The three anchors are the only chromosome-level assemblies in the panel — '
      'everything else is 126 to 2,747 scaffolds — and the only three with curated PlasmoDB-68 '
      'GFF3 + BED12 staged.</p>')

    A('<div class="warn"><b>The trap:</b> the two genomes a <i>P. vivax</i> person reaches for '
      'first are the worst anchors in this panel. <b>Sal-I</b> is the classic Sal-1 reference, but '
      'here it is 2,747 scaffolds <i>and</i> its native GFF3 has <b>zero seqid overlap</b> with the '
      'staged assembly — of all 8 strains it is the only one that fails outright, and it would need '
      'chromosome reconciliation before it could donate anything. <b>PvP01</b> is the PlasmoDB '
      'reference but sits at 242 scaffolds here, and its assembly/annotation coordinates already '
      'disagree. Neither is used as an anchor.</div>')

    A('<h3>Every staged input, and how to produce it</h3>'
      '<table><thead><tr><th>input</th><th>workflow</th><th>how to choose it</th>'
      '<th><i>P. vivax</i> example</th></tr></thead><tbody>')
    seen = set()
    for wid in VETTED:
        for name, _kind, _doc in [(d["name"], d["sub"], d["doc"])
                                  for d in graphs[wid]["nodes"].values() if d["kind"] == "in"]:
            if (wid, name) not in INPUT_GUIDE or (wid, name) in seen:
                continue
            seen.add((wid, name))
            how, pv = INPUT_GUIDE[(wid, name)]
            A(f'<tr><td><code class="port">{esc(name)}</code></td>'
              f'<td><span class="pill" style="--c:{WFCOL[wid]}">{wid}</span></td>'
              f'<td>{how}</td><td class="mut">{pv}</td></tr>')
    A('</tbody></table>')
    A(f'<div class="tip">{DERIVED_NOTE}</div>')
    A('</section></div>')

    # ---- graph -----------------------------------------------------------
    A('<div class="wrap"><section id="graph"><h2>What each workflow does</h2>'
      '<p class="note">Four workflows have been run and verified end to end. Here is what each one '
      'is for, in plain terms, before the step-by-step diagram below.</p><div class="why-grid">')
    for wid in VETTED:
        g = graphs[wid]
        A(f'<div class="why" style="--c:{WFCOL[wid]}">'
          f'<div class="whyhdr"><span class="badge">{esc(g["ph"])}</span>'
          f'<b>{esc(g["title"])}</b>'
          f'<span class="ok">&#10003; {esc(STATUS[wid]["jobs"])}</span></div>'
          f'{md_block(desc[wid]["description"])}</div>')
    A('</div>'
      '<h3>The diagram</h3>'
      '<p class="note">Every step of those four workflows, wired input&rarr;output straight from '
      'the workflow files. Flow runs top to bottom; dashed grey edges in the right-hand lane '
      'cross workflow boundaries.</p></section></div>')

    A('<div class="legend"><span><i class="sw in"></i>workflow input</span>'
      '<span><i class="sw st"></i>tool step</span>'
      '<span><i class="sw out"></i>workflow output</span>'
      '<span><i class="ln"></i>internal edge</span>'
      '<span><i class="ln x"></i>cross-workflow edge</span>'
      '<span class="tools"><span class="tl">focus</span>'
      + "".join(f'<button class="fc" data-wf="{w}" style="--c:{WFCOL[w]}">{w}</button>' for w in VETTED)
      + '<button id="zf">all</button><button id="zo">&minus;</button><button id="zi">+</button></span>'
      '<span class="hint">click a node to trace up- and downstream &middot; drag to pan</span></div>')

    A('<div id="canvas">' + "".join(S) + '</div>')

    A('<div class="wrap">')

    # per-workflow port detail
    def edges_from(wid, port):
        return [(t, ti) for f, fo, t, ti in EDGES if f == wid and fo == port]

    def edges_into(wid, port):
        return [(f, fo) for f, fo, t, ti in EDGES if t == wid and ti == port]

    A('<section id="ports"><h2>Port detail</h2><p class="note">The contract behind the input row at the top and '
      'the output row at the bottom of every cluster above. Ports and docs are parsed from the workflow files; '
      'shapes are what the verified run actually produced.</p></section>')

    for wid in VETTED:
        g = graphs[wid]
        st = STATUS[wid]
        c = WFCOL[wid]
        ins = [(d["name"], d["sub"], d["doc"]) for d in g["nodes"].values() if d["kind"] == "in"]
        outs = [(d["name"], d["doc"]) for d in g["nodes"].values() if d["kind"] == "out"]

        A(f'<section class="wf"><div class="wfhdr" style="--c:{c}">'
          f'<span class="badge">{esc(g["ph"])}</span><h2>{esc(g["title"])}</h2>'
          f'<span class="ok">&#10003; {esc(st["jobs"])}</span>'
          f'<code class="path">{esc(g["path"])}</code></div>')
        A(f'<div class="what">{md_block(desc[wid]["summary"])}</div>')
        ev = (f'History <code>{esc(st["history"])}</code>'
              + (f' &middot; invocation <code>{esc(st["invocation"])}</code>' if st.get("invocation") else "")
              + f' &middot; {esc(st["when"])}')
        A(f'<div class="evidence"><b>Verified:</b> {esc(st["jobs"])} &middot; {ev}<br>'
          f'<span class="mut">{esc(st["note"])}</span></div>')

        A('<h3>Inputs</h3><table><thead><tr><th>port</th><th>kind</th><th>shape</th>'
          '<th>filled by</th><th>description</th></tr></thead><tbody>')
        for name, kind, doc in ins:
            src = edges_into(wid, name)
            origin = (" ".join(f'<span class="pill" style="--c:{WFCOL.get(f, BRC['inkLight'])}">{f}</span> <code>{esc(fo)}</code>'
                               for f, fo in src) if src else
                      f'<span class="pill ext">staged</span> <span class="mut">{esc(EXTERNAL.get((wid, name), ""))}</span>')
            A(f'<tr><td><code class="port">{esc(name)}</code></td><td class="mut">{esc(kind)}</td>'
              f'<td class="mut">{esc(SHAPE.get((wid, "in", name), ""))}</td><td>{origin}</td>'
              f'<td>{esc(doc)}</td></tr>')
        A('</tbody></table>')

        A('<h3>Outputs</h3><table><thead><tr><th>port</th><th>shape</th><th>consumed by</th>'
          '<th>description</th></tr></thead><tbody>')
        for name, doc in outs:
            dests = edges_from(wid, name)
            cons = (" ".join(f'<span class="pill{"" if t in VETTED else " todo"}" style="--c:{WFCOL.get(t, BRC['inkLight'])}">'
                             f'{t}</span> <code>{esc(ti)}</code>' for t, ti in dests) if dests else
                    '<span class="mut">terminal &mdash; deliverable</span>')
            body = doc or FALLBACK_DOC.get((wid, "out", name), "")
            gap = "" if doc else ' <span class="gap">no doc</span>'
            A(f'<tr><td><code class="port">{esc(name)}</code></td>'
              f'<td class="mut">{esc(SHAPE.get((wid, "out", name), ""))}</td>'
              f'<td>{cons}</td><td>{esc(body)}{gap}</td></tr>')
        A('</tbody></table></section>')

    # outbound contracts
    outbound = [(f, fo, t, ti) for f, fo, t, ti in EDGES if f in VETTED and t not in VETTED]
    A('<section id="downstream"><h2>Downstream contracts</h2><p class="note">Ports the eight remaining workflows '
      'depend on. Changing one of these invalidates work already banked.</p>'
      '<table><thead><tr><th>from</th><th>output port</th><th>consumed by</th><th>shape</th></tr>'
      '</thead><tbody>')
    for f, fo, t, ti in outbound:
        A(f'<tr><td><span class="pill" style="--c:{WFCOL.get(f, BRC['inkLight'])}">{f}</span></td>'
          f'<td><code>{esc(fo)}</code></td>'
          f'<td><span class="pill todo">{t}</span> <code>{esc(ti)}</code></td>'
          f'<td class="mut">{esc(SHAPE.get((f, "out", fo), ""))}</td></tr>')
    A('</tbody></table></section>')

    A('<section id="gaps"><h2>Known gaps in these four</h2><div class="caveats">')
    for title, body in CAVEATS:
        A(f'<div class="cav"><b>{title}</b><div class="mut">{body}</div></div>')
    A('</div></section></div>')

    A('</div>')   # end overview pane

    for wid in tabbed:
        g = graphs[wid]; cap = examples[wid]; c = WFCOL[wid]
        st = STATUS[wid]
        A(f'<div class="tabpane" id="tab-wf-{wid}"><div class="wrap">')
        A(f'<section><div class="wfhdr" style="--c:{c}">'
          f'<span class="badge">{esc(g["ph"])}</span><h2>{esc(g["title"])}</h2>'
          f'<span class="ok">&#10003; {esc(st["jobs"])}</span></div>')
        A(md_block(desc[wid]["description"]))
        A('</section>')

        A(f'<section><h3>Steps</h3><p class="note">Click any step to see what it actually '
          f'produced. The samples come from invocation <code>{esc(cap["invocation"])}</code> '
          f'({esc(cap["when"])}, {esc(cap["state"])}) &mdash; a real run on the 8-strain panel, '
          f'not synthetic data. Nodes with a sample are outlined; the rest carried no output '
          f'of their own.</p>')
        A(f'<div class="solowrap">{solo_svg(wid, g, allnodes)}</div>')
        A(f'<div class="expanel" id="ex-{wid}"><div class="exempty">Select a step above.</div></div>')
        A('</section></div></div>')

    A('<footer>Step graphs, ports and docs parsed from <code>workflows/*/*.gxwf.yml</code>; '
      'cross-workflow wiring imported from <code>gen_pipeline_dataflow.EDGES</code> so this and the '
      '<a href="pipeline_dataflow.html">pipeline overview</a> cannot disagree. That overview is the '
      'front page and shows all twelve workflows; this page carries the six that have been '
      'reviewed, one tab each. Run evidence read off '
      'the live Galaxy. Generated by <code>workflows/gen_pipeline_io_map.py</code> &mdash; do not edit '
      'by hand.</footer>')

    clusters = {wid: dict(x=g["ox"], y=g["oy"], w=g["w"], h=g["h"], title=g["title"])
                for wid, g in graphs.items()}

    adj = {}
    for wid, g in graphs.items():
        for a, b in g["edges"]:
            adj.setdefault(a, []).append(b)
    for a, b in cross:
        adj.setdefault(a, []).append(b)
    rev = {}
    for a, bs in adj.items():
        for b in bs:
            rev.setdefault(b, []).append(a)

    out = os.path.join(ROOT, "workflows/pipeline_io_map.html")
    with open(out, "w") as fh:
        fh.write(TEMPLATE
                 .replace("__BODY__", "".join(P))
                 .replace("__ADJ__", json.dumps(adj))
                 .replace("__REV__", json.dumps(rev))
                 .replace("__CLUSTERS__", json.dumps(clusters))
                 .replace("__EXAMPLES__", json.dumps(ex_payload))
                 .replace("__STALE__", json.dumps(stale)))
    n = sum(len(g["nodes"]) for g in graphs.values())
    e = sum(len(g["edges"]) for g in graphs.values())
    print(f"wrote workflows/pipeline_io_map.html — {len(VETTED)} workflows, {n} nodes, "
          f"{e} internal edges, {len(cross)} cross-workflow edges, canvas {CW}x{CH}")


TEMPLATE = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Pv4 pangenome pipeline — inputs, workflows and data flow</title>
<style>
/* Palette: BRC Analytics (primary #28285B) over the findable-ui base it builds on. */
:root{
  --bg:#F6F6F7; --surface:#FFFFFF; --ink:#212B36; --mut:#637381;
  --line:#E1E3E5; --line2:#C4CDD5; --smoke:#FAFBFB;
  --primary:#28285B; --primary-dark:#1F1F47;
  --info:#00729C; --info-l:#97D6EA; --info-lt:#F2FAFC;
  --ok:#287555; --ok-l:#AEE9D1; --ok-lt:#F1F8F5;
  --warn:#B54708; --warn-l:#FFD79D; --warn-lt:#FFFAEB;
  --alert:#B42318; --alert-lt:#FFF4F4;
  --caution:#956F00; --caution-l:#FFEB78;
}
*{box-sizing:border-box}
body{margin:0;font:14px/1.6 -apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;background:var(--bg);color:var(--ink)}
header{padding:22px 32px;background:var(--primary);color:#fff;display:flex;gap:22px;align-items:center;flex-wrap:wrap}
header h1{font-size:22px;margin:0 0 4px;font-weight:600;letter-spacing:-.01em}
header .sub{color:rgba(255,255,255,.82);font-size:13.5px;max-width:80ch}
.score{margin-left:auto;font-size:12.5px;border:1px solid rgba(255,255,255,.3);border-radius:8px;padding:8px 14px;white-space:nowrap;color:rgba(255,255,255,.9)}
.score b{color:#fff;font-size:17px}
nav{position:sticky;top:0;z-index:30;background:var(--surface);border-bottom:1px solid var(--line);padding:0 32px;display:flex;gap:4px;flex-wrap:wrap;box-shadow:0 1px 2px rgba(33,43,54,.05)}
nav button.tab{background:none;border:0;border-bottom:2px solid transparent;color:var(--mut);font:inherit;font-size:13px;padding:13px 16px;cursor:pointer}
nav button.tab:hover{color:var(--primary);border-bottom-color:var(--line2)}
nav button.tab.active{color:var(--c,var(--primary));border-bottom-color:var(--c,var(--primary));font-weight:600}
.tabpane{display:none}
.tabpane.active{display:block}
.solowrap{overflow-x:hidden;overflow-y:auto;background:var(--surface);border:1px solid var(--line);border-radius:8px;
  background-image:radial-gradient(circle at 1px 1px,var(--line) 1px,transparent 0);background-size:22px 22px;
  max-height:74vh;padding:4px}
/* scale to the pane rather than scrolling sideways; the viewBox keeps it sharp */
svg.solo{display:block;margin:8px auto;max-width:100%;height:auto}
svg.solo .n rect{fill:var(--surface);stroke:var(--line2);stroke-width:1.2}
svg.solo .n.in rect{fill:var(--info-lt);stroke:var(--info)}
svg.solo .n.out rect{fill:var(--ok-lt);stroke:var(--ok)}
svg.solo .n.has-ex{cursor:pointer}
svg.solo .n.has-ex rect{stroke-width:2.2;stroke-dasharray:none}
svg.solo .n:not(.has-ex){opacity:.55}
svg.solo .n.has-ex:hover rect{stroke:var(--primary);stroke-width:2.8}
svg.solo .n.sel rect{stroke:var(--primary);stroke-width:3}
svg.solo .n.hot rect{stroke:var(--primary);stroke-width:2.2}
svg.solo .n.dim{opacity:.18}
svg.solo .e.hot{opacity:1;stroke-width:2.6}
svg.solo .e.dim{opacity:.08}
svg.solo .nn{fill:var(--ink);font:600 11.5px ui-monospace,SFMono-Regular,Menlo,monospace}
svg.solo .ns{fill:var(--mut);font:9.5px -apple-system,Segoe UI,Roboto,sans-serif}
svg.solo .e{stroke-width:1.5;opacity:.5}
.expanel{margin-top:14px;border:1px solid var(--line);border-radius:8px;background:var(--surface);padding:14px 16px;min-height:90px}
.exempty{color:var(--mut);font-size:13px;font-style:italic}
.backlink{margin-left:16px;font-size:12.5px;color:var(--mut);text-decoration:none;
  border:1px solid var(--line);border-radius:12px;padding:3px 10px;white-space:nowrap}
.backlink:hover{color:var(--ink);border-color:var(--line2);background:var(--surface)}
.excols{border-collapse:collapse;margin:10px 0 6px;font-size:12.5px;width:100%}
.excols th{text-align:left;font-weight:600;color:var(--mut);border-bottom:1px solid var(--line);
  padding:4px 8px 4px 0;font-size:11px;text-transform:uppercase;letter-spacing:.04em}
.excols td{padding:4px 8px 4px 0;vertical-align:top;border-bottom:1px solid var(--line)}
.excols td:first-child{color:var(--mut);width:1.6em;text-align:right}
.excols code{font-size:12px;white-space:nowrap}
.exwhy{font-size:13.5px;line-height:1.65;max-width:88ch;border-left:3px solid var(--info);
  background:var(--info-lt);border-radius:0 6px 6px 0;padding:11px 15px;margin-bottom:14px}
.exwhy p{margin:0 0 8px}
.exwhy p:last-child{margin-bottom:0}
svg.fig{display:block;max-width:100%;height:auto;margin:10px 0 4px;background:var(--surface);
  border:1px solid var(--line);border-radius:6px}
svg.fig .fb rect{fill:var(--surface);stroke:var(--line2);stroke-width:1}
svg.fig .fb text{font:10.5px ui-monospace,Menlo,monospace;fill:var(--ink)}
svg.fig .fb.in rect{fill:var(--info-lt);stroke:var(--info)}
svg.fig .fb.out rect{fill:var(--smoke)}
svg.fig .fb.id rect{fill:var(--surface);stroke-dasharray:3 2}
svg.fig .fb.diag rect{fill:var(--warn-lt);stroke:var(--warn)}
svg.fig .fcap{font:10.5px -apple-system,Segoe UI,Roboto,sans-serif;fill:var(--mut)}
svg.fig .farrow{stroke:var(--line2);stroke-width:1.4;fill:none;marker-end:none}
.exwhy pre{background:var(--surface);border:1px solid var(--line);border-radius:5px;
  padding:8px 10px;font:11.5px/1.5 ui-monospace,Menlo,monospace;overflow-x:auto;white-space:pre}
.exhead{display:flex;gap:10px;align-items:baseline;flex-wrap:wrap;margin-bottom:4px}
.exhead b{font-size:14px}
.exmeta{color:var(--mut);font-size:11.5px}
.exout{border-top:1px solid var(--line);padding-top:10px;margin-top:10px}
.exout:first-of-type{border-top:0;padding-top:0;margin-top:0}
.expre{background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:9px 11px;
  font:11.5px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;overflow-x:auto;white-space:pre;margin:6px 0 0}
.eximg{max-width:100%;border:1px solid var(--line);border-radius:6px;margin-top:6px;background:#fff}
.exids{color:var(--mut);font-size:11.5px;margin-top:5px}
.wrap{padding:26px 32px;max-width:1180px}
section{margin:0 0 30px}
h2{font-size:19px;margin:0 0 10px;font-weight:600;letter-spacing:-.01em}
h3{font-size:15px;margin:26px 0 8px;font-weight:600}
.note{color:var(--mut);font-size:13.5px;margin:0 0 12px;max-width:90ch}
blockquote{margin:14px 0;padding:14px 18px;background:var(--info-lt);border-left:3px solid var(--info);border-radius:0 8px 8px 0;font-size:13.5px;max-width:92ch}
ol.rules{margin:10px 0;padding-left:22px;max-width:92ch}
ol.rules li{margin-bottom:9px;color:var(--mut);font-size:13.5px}
ol.rules b{color:var(--ink)}
table{width:100%;border-collapse:collapse;font-size:13px;background:var(--surface);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:6px 0}
thead th{text-align:left;font-weight:600;color:var(--mut);font-size:11px;text-transform:uppercase;letter-spacing:.05em;padding:9px 13px;border-bottom:1px solid var(--line);background:var(--smoke)}
td{padding:9px 13px;border-bottom:1px solid var(--line);vertical-align:top}
tbody tr:last-child td{border-bottom:0}
table.panel tr.anc{background:var(--info-lt)}
table.panel td b{font-variant-numeric:tabular-nums}
code{font:12px/1.45 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--bg);border:1px solid var(--line);border-radius:4px;padding:1px 5px}
code.port{background:var(--info-lt);border-color:var(--info-l);color:var(--info)}
code.path{background:none;border:0;color:var(--mut);margin-left:auto;font-size:11.5px}
.mut{color:var(--mut)}
.pill{display:inline-block;background:var(--c,var(--mut));color:#fff;border-radius:5px;padding:0 8px;font-size:11px;font-weight:700;line-height:19px}
.pill.ext{background:var(--mut)}
.pill.todo{background:none;border:1px dashed var(--line2);color:var(--mut);font-weight:600}
.warn{background:var(--warn-lt);border:1px solid var(--warn-l);border-left:3px solid var(--warn);border-radius:0 8px 8px 0;padding:13px 17px;font-size:13.5px;margin:14px 0;max-width:92ch}
.tip{background:var(--ok-lt);border:1px solid var(--ok-l);border-left:3px solid var(--ok);border-radius:0 8px 8px 0;padding:13px 17px;font-size:13.5px;margin:14px 0;max-width:92ch}
.legend{display:flex;gap:18px;align-items:center;padding:10px 32px;color:var(--mut);font-size:12px;flex-wrap:wrap;background:var(--surface);border-top:1px solid var(--line);border-bottom:1px solid var(--line)}
.legend>span{display:flex;gap:7px;align-items:center}
.legend .sw{width:14px;height:12px;border-radius:3px;display:inline-block;border:1.6px solid var(--line2);background:var(--surface)}
.legend .sw.in{background:var(--info-lt);border-color:var(--info)}
.legend .sw.out{background:var(--ok-lt);border-color:var(--ok)}
.legend .ln{width:22px;border-top:2px solid var(--primary);display:inline-block}
.legend .ln.x{border-top:2px dashed var(--mut)}
.legend .hint{margin-left:auto;font-style:italic}
.tools{display:flex;gap:6px;align-items:center}
.tools .tl{font-size:11px;text-transform:uppercase;letter-spacing:.06em}
.tools button{background:var(--surface);color:var(--ink);border:1px solid var(--line2);border-radius:6px;padding:4px 10px;cursor:pointer;font-size:12px}
.tools button:hover{border-color:var(--mut)}
.tools .fc{border-color:var(--c);background:var(--c);color:#fff;font-weight:700;min-width:32px}
.tools .fc:hover{filter:brightness(1.12)}
#canvas{overflow:auto;background:var(--surface);
  background-image:radial-gradient(circle at 1px 1px,var(--line) 1px,transparent 0);
  background-size:22px 22px;cursor:grab;border-bottom:1px solid var(--line);max-height:80vh}
#canvas.drag{cursor:grabbing}
#canvas svg{display:block;margin:20px}
.cbox{fill:var(--smoke);stroke:var(--c);stroke-width:1.2;opacity:.95}
.ctitle{fill:#fff;font:600 12.5px -apple-system,Segoe UI,Roboto,sans-serif}
.cmeta{fill:rgba(255,255,255,.9);font:600 10.5px ui-monospace,Menlo,monospace}
.n rect{fill:var(--surface);stroke:var(--line2);stroke-width:1.2}
.n.in rect{fill:var(--info-lt);stroke:var(--info)}
.n.out rect{fill:var(--ok-lt);stroke:var(--ok)}
.n .nn{fill:var(--ink);font:600 11.5px ui-monospace,SFMono-Regular,Menlo,monospace}
.n .ns{fill:var(--mut);font:9.5px -apple-system,Segoe UI,Roboto,sans-serif}
.n{cursor:pointer}
.n:hover rect{stroke:var(--primary);stroke-width:2}
.e{fill:none;stroke-width:1.5;opacity:.55}
.e.x{stroke-dasharray:6 4;stroke-width:2;opacity:.85}
.dim{opacity:.1}
.n.dim rect,.n.dim text{opacity:.3}
.hot rect{stroke:var(--primary);stroke-width:2.4}
.e.hot{opacity:1;stroke-width:2.6}
.why-grid{display:grid;gap:14px;margin:14px 0 4px}
.why{background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--c);border-radius:0 8px 8px 0;padding:14px 18px;box-shadow:0 1px 2px rgba(33,43,54,.04)}
.whyhdr{display:flex;align-items:center;gap:10px;margin-bottom:7px;flex-wrap:wrap}
.whyhdr b{font-size:15px}
.why p{margin:0;font-size:13.5px;color:var(--ink);max-width:95ch;line-height:1.65}
.wf{border:1px solid var(--line);border-radius:10px;background:var(--surface);padding:0 18px 18px;box-shadow:0 1px 2px rgba(33,43,54,.04)}
.wfhdr{display:flex;align-items:center;gap:12px;padding:14px 0;border-bottom:1px solid var(--line);margin-bottom:12px;flex-wrap:wrap}
.wfhdr h2{margin:0;font-size:17px}
.badge{background:var(--c);color:#fff;border-radius:6px;padding:2px 10px;font-size:12px;font-weight:700}
.ok{color:var(--ok);font-size:11.5px;font-weight:600;background:var(--ok-lt);border:1px solid var(--ok-l);border-radius:6px;padding:2px 9px}
.what{color:var(--mut);font-size:13.5px;margin:0 0 12px;max-width:94ch}
.what p{margin:0 0 6px}
.what p:last-child{margin-bottom:0}
.evidence{background:var(--ok-lt);border:1px solid var(--ok-l);border-left:3px solid var(--ok);border-radius:0 6px 6px 0;padding:10px 14px;font-size:12.5px}
.gap{color:var(--alert);background:var(--alert-lt);border:1px solid var(--alert);border-radius:4px;padding:0 5px;font-size:10px;white-space:nowrap}
.caveats{display:grid;gap:11px}
.cav{background:var(--warn-lt);border:1px solid var(--warn-l);border-left:3px solid var(--warn);border-radius:0 8px 8px 0;padding:12px 16px;font-size:13.5px}
.cav b{display:block;margin-bottom:4px}
footer{color:var(--mut);font-size:12px;padding:18px 32px 34px;border-top:1px solid var(--line);background:var(--surface);max-width:none}
footer a{color:var(--info)}
@media (max-width:820px){table{display:block;overflow-x:auto}.wrap{padding:18px}header,nav,.legend,footer{padding-left:18px;padding-right:18px}}
</style></head><body>__BODY__
<script>
const ADJ=__ADJ__, REV=__REV__, CLUSTERS=__CLUSTERS__;
const svg=document.getElementById('g'), cv=document.getElementById('canvas');
const W=+svg.getAttribute('width'), H=+svg.getAttribute('height'), M=20;
let z=1;
function setz(v){z=Math.max(.2,Math.min(2.2,v));svg.setAttribute('width',W*z);svg.setAttribute('height',H*z);}
document.getElementById('zi').onclick=()=>setz(z*1.2);
document.getElementById('zo').onclick=()=>setz(z/1.2);
document.getElementById('zf').onclick=()=>{
  setz(Math.max(.34,Math.min(1,(cv.clientWidth-2*M-8)/W)));cv.scrollTo({left:0,top:0,behavior:'smooth'});
};
document.querySelectorAll('.fc').forEach(b=>b.onclick=()=>{
  const c=CLUSTERS[b.dataset.wf]; if(!c) return;
  setz(Math.max(.34,Math.min(1.4,(cv.clientWidth-2*M-30)/c.w)));
  cv.scrollTo({left:(c.x+M)*z-14,top:(c.y+M)*z-14,behavior:'smooth'});
});
function walk(start,map){const seen=new Set(),st=[start];while(st.length){const n=st.pop();for(const m of (map[n]||[])){if(!seen.has(m)){seen.add(m);st.push(m);}}}return seen;}
function clear(){document.querySelectorAll('.n,.e').forEach(e=>e.classList.remove('dim','hot'));}
function trace(id){
  clear();
  const keep=new Set([id,...walk(id,ADJ),...walk(id,REV)]);
  document.querySelectorAll('.n').forEach(n=>{n.classList.toggle('dim',!keep.has(n.id));if(n.id===id)n.classList.add('hot');});
  document.querySelectorAll('.e').forEach(e=>{
    const on=keep.has(e.dataset.a)&&keep.has(e.dataset.b);
    e.classList.toggle('dim',!on); e.classList.toggle('hot',on);
  });
}
document.querySelectorAll('.n').forEach(n=>n.addEventListener('click',ev=>{ev.stopPropagation();trace(n.id);}));
cv.addEventListener('click',clear);
let down=false,sx,sy,sl,stp;
cv.addEventListener('mousedown',e=>{if(e.target.closest('.n'))return;down=true;sx=e.clientX;sy=e.clientY;sl=cv.scrollLeft;stp=cv.scrollTop;cv.classList.add('drag');});
window.addEventListener('mouseup',()=>{down=false;cv.classList.remove('drag');});
window.addEventListener('mousemove',e=>{if(!down)return;cv.scrollLeft=sl-(e.clientX-sx);cv.scrollTop=stp-(e.clientY-sy);});
window.addEventListener('load',()=>document.getElementById('zf').click());

// ---- tabs ----
function showTab(name, scroll){
  const pane = document.getElementById('tab-'+name);
  const btn  = document.querySelector('nav .tab[data-tab="'+name+'"]');
  if(!pane || !btn) return false;
  document.querySelectorAll('nav .tab').forEach(x=>x.classList.remove('active'));
  document.querySelectorAll('.tabpane').forEach(x=>x.classList.remove('active'));
  btn.classList.add('active');
  pane.classList.add('active');
  if(scroll !== false) window.scrollTo({top:0});
  return true;
}
document.querySelectorAll('nav .tab').forEach(b=>b.onclick=()=>{
  showTab(b.dataset.tab, true);
  history.replaceState(null, '', b.dataset.tab === 'overview' ? location.pathname : '#tab-'+b.dataset.tab);
});
// A pane is display:none until its tab is active, so a deep link from the
// pipeline overview (#tab-wf-C2) would otherwise land on a hidden element.
function tabFromHash(){
  const m = /^#tab-(.+)$/.exec(location.hash || '');
  if(m) showTab(m[1], false);
}
tabFromHash();
window.addEventListener('hashchange', tabFromHash);

// ---- per-step examples ----
const EXAMPLES=__EXAMPLES__;
const STALE=__STALE__;
function esc(s){return String(s).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));}
function human(n){return n>=1e6?(n/1e6).toFixed(1)+' MB':n>=1e3?(n/1e3).toFixed(1)+' kB':n+' B';}
function renderExample(wf,node){
  const panel=document.getElementById('ex-'+wf); if(!panel) return;
  const ex=EXAMPLES[wf+'|st|'+node]||EXAMPLES[wf+'|in|'+node]||EXAMPLES[wf+'|out|'+node];
  if(!ex){panel.innerHTML='<div class="exempty">No output was captured for this step.</div>';return;}
  let h='';
  if(ex.why) h+='<div class="exwhy">'+ex.why+'</div>';
  ex.outputs.forEach(o=>{
    h+='<div class="exout"><div class="exhead"><b>'+esc(o.name)+'</b>';
    h+='<span class="exmeta">'+esc(o.ext||'?')+' · '+human(o.bytes||0);
    if(o.kind==='collection') h+=' · collection of '+o.elements+', showing <code>'+esc(o.shown)+'</code>';
    if(o.total_lines!==undefined) h+=' · '+o.total_lines.toLocaleString()+' lines';
    h+='</span></div>';
    if(o.columns){
      h+='<table class="excols"><thead><tr><th>#</th><th>column</th><th>meaning</th></tr></thead><tbody>';
      o.columns.forEach((c,i)=>{h+='<tr><td>'+(i+1)+'</td><td><code>'+esc(c[0])+'</code></td><td>'+esc(c[1])+'</td></tr>';});
      h+='</tbody></table>';
    }
    if(o.image) h+='<img class="eximg" src="'+o.image+'" alt="'+esc(o.name)+'">';
    else if(o.text) h+='<pre class="expre">'+esc(o.text)+(o.truncated?'\n…':'')+'</pre>';
    else if(o.note) h+='<div class="exmeta">'+esc(o.note)+'</div>';
    if(o.element_ids&&o.element_ids.length>1)
      h+='<div class="exids">elements: '+o.element_ids.map(esc).join(', ')+'</div>';
    h+='</div>';
  });
  panel.innerHTML=h;
}
// clicking any node in a per-workflow diagram traces the path through it and,
// when a sample was captured for that step, shows the sample too
document.querySelectorAll('svg.solo .n').forEach(n=>n.addEventListener('click',()=>{
  const wf=n.dataset.wf, svg=document.getElementById('solo-'+wf);
  const id=n.dataset.id;
  const keep=new Set([id,...walk(id,ADJ),...walk(id,REV)]);
  svg.querySelectorAll('.n').forEach(x=>{
    x.classList.remove('sel','dim','hot');
    if(x.dataset.id===id) x.classList.add('sel');
    else if(!keep.has(x.dataset.id)) x.classList.add('dim');
    else x.classList.add('hot');
  });
  svg.querySelectorAll('.e').forEach(e=>{
    const on=keep.has(e.dataset.a)&&keep.has(e.dataset.b);
    e.classList.toggle('dim',!on); e.classList.toggle('hot',on);
  });
  if(n.classList.contains('has-ex')) renderExample(wf,n.dataset.node);
  else document.getElementById('ex-'+wf).innerHTML='<div class="exempty">'+(STALE[id]
    ? 'No sample for this node in invocation '+esc(STALE[id])+'. Re-run workflows/capture_examples.py against a newer invocation to add one.'
    : 'Path highlighted. This step produced no output of its own.')+'</div>';
}));
// click the background of a diagram to clear the trace
document.querySelectorAll('svg.solo').forEach(svg=>svg.addEventListener('click',ev=>{
  if(ev.target.closest('.n')) return;
  svg.querySelectorAll('.n,.e').forEach(x=>x.classList.remove('sel','dim','hot'));
}));
</script>
</body></html>'''


if __name__ == "__main__":
    main()
