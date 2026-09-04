#!/usr/bin/env python3
"""Stage the real WF-A input set on usegalaxy.org: 6 assemblies, 6 proteomes, 2 anchors.

⚠ SIX, NOT SIX, AND THE DIFFERENCE WAS MY OWN OVER-RESTRICTION. The first version took only the
PANEL-eligible assemblies (chromosome-level or better), which drops the two contig-level Jamaican
Lion genomes. That filter exists for MASKING -- WindowMasker estimates its frequency model from the
assembly, so fragmentation degrades the mask. WF-A is sourmash, BUSCO and bookkeeping; none of it
cares about contiguity. All six NCBI-annotated assemblies have a proteome AND an anchor-usable GFF3,
so all six belong here.

⚠ TWO OF THE SIX ARE NOT IN THE STAGED PANEL COLLECTION, for that same reason, so their FASTAs are
fetched from NCBI server-side rather than copied.

⛔ THE PROTEOMES ARE THE `_primary` COPIES, NOT WHAT NCBI SHIPS. RefSeq publishes every isoform, the
GenBank submitter annotations here publish one protein per gene, and BUSCO's duplication figure is
read directly off that: feeding cs10's 33,674 isoform records beside the T2T pair's 31,109 one-per-
gene records would report cs10 as massively duplicated for a reason that is not biology.
"""
from __future__ import annotations

import json
import os
import pathlib
import sys
import time
import urllib.request

#: Where the `_primary` proteome and anchor GFF3 files live. ⚠ NOT IN THIS REPOSITORY -- they are
#: several GB of NCBI downloads, so this points outside it and $WFA_PROTEOMES overrides. The script
#: refuses rather than staging a partial panel if a file is missing, because a collection that is
#: silently short is far worse than one that was never built.
PROTEOMES = pathlib.Path(os.environ.get("WFA_PROTEOMES", "proteomes")).expanduser()

#: Where the resulting collection ids are written, for the workflow driver to read.
OUT_DIR = pathlib.Path(os.environ.get("WFA_OUT_DIR", ".")).expanduser()
RAW_HDCA = "cab4808ec6fe5c51"          # the staged 19-genome raw panel

#: element identifier -> assembly_name its proteome/GFF3 files are keyed by
SIX = {
    "cs10_NCBI_RefSeq_softmasked": "cs10",
    "ASM2916894v1": "ASM2916894v1",
    "T2T_GCA_054642775.1_IMPD": "ASM5464277v1",
    "T2T_GCA_054642815.1_IMPD": "ASM5464281v1",
    "JL_Father": "JL_Father",
    "JL_Mother": "JL_Mother",
}
#: the two the panel collection does not hold; fetched by URL instead of copied
BY_URL = {
    "cs10_NCBI_RefSeq_softmasked": ("GCF_900626175.2", "cs10"),
    "ASM2916894v1": ("GCF_029168945.1", "ASM2916894v1"),
    "T2T_GCA_054642775.1_IMPD": ("GCA_054642775.1", "ASM5464277v1"),
    "T2T_GCA_054642815.1_IMPD": ("GCA_054642815.1", "ASM5464281v1"),
    "JL_Father": ("GCA_013030025.1", "JL_Father"),
    "JL_Mother": ("GCA_012923435.1", "JL_Mother"),
}

#: ⚠ Copying from the pre-staged panel is an optimisation available on ONE server. Everything has a
#: URL, so a server that holds nothing pre-staged fetches all six and gets identical bytes.
COPY_IF_PRESENT = os.environ.get("WFA_SERVER", "") == ""
#: ⚠ ANCHORS ARE A SUBSET, AND THE WORKFLOW SAYS SO ("shorter than assemblies"). The two RefSeq
#: annotations are the curated ones; the T2T pair are submitter annotations, fine as targets but a
#: different provenance to mix into a reference set.
ANCHORS = {"cs10_NCBI_RefSeq_softmasked": "cs10", "ASM2916894v1": "ASM2916894v1"}


def creds() -> tuple[str, str]:
    """Server chosen by $WFA_SERVER: "" for usegalaxy.org, "_3" for laila.

    ⚠ ONE SCRIPT, TWO SERVERS, because the two staging runs must produce the SAME inputs. Forking it
    per server is how the collections quietly drift apart and a difference in the RESULT gets
    blamed on the workflow.
    """
    sfx = os.environ.get("WFA_SERVER", "")
    u = os.environ.get(f"GALAXY_URL{sfx}", "").rstrip("/")
    k = os.environ.get(f"GALAXY_API_KEY{sfx}", "")
    if not (u and k):
        sys.exit(f"set GALAXY_URL{sfx} and GALAXY_API_KEY{sfx}")
    print(f"  server {u}")
    return u, k


URL, KEY = creds()


def api(path, payload=None):
    r = urllib.request.Request(URL + path, method="POST" if payload is not None else "GET",
                               headers={"x-api-key": KEY, "content-type": "application/json"},
                               data=json.dumps(payload).encode() if payload is not None else None)
    with urllib.request.urlopen(r, timeout=900) as f:
        return json.load(f)


def upload(history: str, path: pathlib.Path, name: str, ext: str) -> str:
    """Upload a LOCAL file through /api/tools/fetch as multipart.

    ⛔ NOT bioblend's upload_file, AND NOT `upload1`. Both fail here, for different reasons worth
    recording:

      * bioblend uses the TUS protocol on any Galaxy >= 22.01 with no way to opt out, and the TUS
        endpoint Galaxy hands back is built from ITS OWN configured base URL -- `localhost:8080` on
        laila. Behind an SSH tunnel that address is the sandbox's own localhost, not Galaxy's, so
        every chunk goes nowhere. `TusUploadFailed ... Max retries exceeded`.
      * `upload1` is the classic upload tool and simply is not in laila's panel, which carries 63
        tools. `Tool not found`, 400014.

    /api/tools/fetch takes the file as a multipart part and the plan as a JSON string, and is
    present on both servers.

    ⚠ `trust_env=False` for the tunnel. host.docker.internal must NOT go through the sandbox's
    outbound proxy; usegalaxy.org must.
    """
    import requests
    local = "host.docker.internal" in URL or "localhost" in URL
    s = requests.Session()
    s.headers.update({"x-api-key": KEY})
    s.trust_env = not local
    targets = [{"destination": {"type": "hdas"},
                "elements": [{"src": "files", "name": name, "ext": ext,
                              "to_posix_lines": False, "space_to_tab": False}]}]
    last = None
    for attempt in range(5):
        try:
            with path.open("rb") as fh:
                r = s.post(f"{URL}/api/tools/fetch",
                           data={"history_id": history, "targets": json.dumps(targets)},
                           files={"files_0|file_data": (name, fh)}, timeout=3600)
            r.raise_for_status()
            outs = r.json().get("outputs") or []
            if not outs:
                raise RuntimeError(f"fetch returned no output: {r.text[:200]}")
            return outs[0]["id"]
        except Exception as exc:                                    # noqa: BLE001
            last = exc
            print(f"      upload {name} attempt {attempt+1} failed "
                  f"({type(exc).__name__}: {str(exc)[:90]}); retrying", flush=True)
            time.sleep(10 * (attempt + 1))
    raise SystemExit(f"upload {name}: gave up after 5 attempts -- {type(last).__name__}: {last}")


def ncbi_fasta(acc: str, name: str) -> str:
    prefix, num = acc.split("_")
    num = num.split(".")[0]
    safe = name.replace(" ", "_")
    return (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{num[0:3]}/{num[3:6]}/{num[6:9]}/"
            f"{acc}_{safe}/{acc}_{safe}_genomic.fna.gz")


def wait(ids, label):
    while True:
        st = {i: api(f"/api/datasets/{i}")["state"] for i in ids}
        bad = [i for i, s in st.items() if s == "error"]
        if bad:
            sys.exit(f"{label}: error on {bad}")
        left = [i for i, s in st.items() if s not in ("ok", "empty")]
        if not left:
            return
        print(f"    {label}: {len(ids)-len(left)}/{len(ids)} ready", flush=True)
        time.sleep(20)


def collection(history, name, pairs):
    return api(f"/api/histories/{history}/contents", {
        "type": "dataset_collection", "collection_type": "list", "name": name,
        "element_identifiers": [{"src": "hda", "id": i, "name": n} for n, i in pairs]})


def main() -> int:
    # ⚠ NAMED FROM len(SIX), NOT A LITERAL. This said "(4 genomes)" for as long as the panel had
    # six in it -- the count was left behind when the two Jamaican Lion genomes were added, and a
    # history whose label disagrees with its contents is read as the contents being wrong.
    hist = api("/api/histories",
               {"name": f"WF-A UDT edition — cannabis test ({len(SIX)} genomes)"})["id"]
    print(f"  history {hist}")

    have = {}
    if COPY_IF_PRESENT:
        raw = api(f"/api/dataset_collections/{RAW_HDCA}?instance_type=history")
        have = {e["element_identifier"]: e["object"]["id"] for e in raw["elements"]}
    missing = [k for k in SIX if k not in have and k not in BY_URL]
    if missing:
        sys.exit(f"raw panel is missing {missing} and no URL is declared for them")

    asm = []
    for ident in SIX:
        if ident in have:
            d = api(f"/api/histories/{hist}/contents",
                    {"source": "hda", "content": have[ident], "type": "dataset"})
            asm.append((ident, d["id"]))
    print(f"  copied {len(asm)} assemblies from the staged panel")

    # ⚠ SERVER-SIDE FETCH, NOT AN UPLOAD. These two are not in the panel collection and their FASTAs
    # are not on this machine; Galaxy pulls them from NCBI directly, which costs no tunnel traffic.
    fetch = [{"src": "url", "url": ncbi_fasta(acc, nm), "name": ident, "ext": "fasta.gz",
              "to_posix_lines": False, "space_to_tab": False}
             for ident, (acc, nm) in BY_URL.items() if ident not in have]
    if fetch:
        r = api("/api/tools/fetch", {"history_id": hist,
                                     "targets": [{"destination": {"type": "hdas"},
                                                  "elements": fetch}]})
        for o in r.get("outputs", []):
            asm.append((o["name"], o["id"]))
        print(f"  fetching {len(fetch)} assembly FASTA(s) from NCBI")

    prot, anch = [], []
    for ident, key in SIX.items():
        f = next(PROTEOMES.glob(f"*_{key}_protein_primary.faa.gz"))
        prot.append((ident, upload(hist, f, f"{ident}.faa.gz", "fasta.gz")))
        print(f"    uploaded proteome {ident}")
    for ident, key in ANCHORS.items():
        f = next(PROTEOMES.glob(f"*_{key}_genomic.gff.gz"))
        anch.append((ident, upload(hist, f, f"{ident}.gff3", "gff3")))
        print(f"    uploaded anchor   {ident}")

    wait([i for _, i in asm + prot + anch], "staging")

    c_asm = collection(hist, "assemblies", asm)
    c_prot = collection(hist, "proteomes", prot)
    c_anch = collection(hist, "anchor_gene_gff3s", anch)
    out = {"server": URL, "history": hist, "assemblies": c_asm["id"], "proteomes": c_prot["id"],
           "anchor_gene_gff3s": c_anch["id"]}
    print(f"\n  assemblies        {c_asm['id']}  ({len(asm)})")
    print(f"  proteomes         {c_prot['id']}  ({len(prot)})")
    print(f"  anchor_gene_gff3s {c_anch['id']}  ({len(anch)})")
    dest = OUT_DIR / f"wfa_inputs{os.environ.get('WFA_SERVER', '')}.json"
    dest.write_text(json.dumps(out, indent=1))
    print(f"  wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
