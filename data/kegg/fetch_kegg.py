#!/usr/bin/env python3
"""Snapshot of the KEGG data used by PACFM 2 (KEGG REST API, academic use; https://www.kegg.jp/kegg/rest/).

Writes into data/kegg/<date>/:
  ko_pathway.tsv            KO -> pathway (ko numbers)              link/pathway/ko
  pathway_hierarchy.tsv     pathway -> level A / level B / name     br08901
  ko_module.tsv             KO -> module                            link/module/ko
  module_definition.tsv     module -> name, class, DEFINITION       get/<modules>
  module_pathway.tsv        module -> pathway                       link/pathway/module
  organisms.tsv             KEGG organism list (prokaryotes kept)   list/organism
  genome_ko/<org>.tsv       gene -> KO for the selected genomes     link/ko/<org>
  genomes_selected.tsv      one complete genome per genus (bacteria + archaea), deterministic choice
Run: python data/kegg/fetch_kegg.py [n_genomes]
"""
import re, sys, time, urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

BASE = "https://rest.kegg.jp"
OUT = Path(__file__).resolve().parent / str(date.today()); OUT.mkdir(parents=True, exist_ok=True)


def get(path, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(f"{BASE}/{path}", timeout=120) as r:
                return r.read().decode()
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def save(name, text):
    (OUT / name).write_text(text); print(name, len(text.splitlines()), "lines", flush=True)


def main(n_genomes=300):
    done = lambda f: (OUT / f).exists() and (OUT / f).stat().st_size > 0
    if done("module_pathway.tsv"):
        print("reference tables present, going to genomes", flush=True)
        return genomes(n_genomes)
    t = get("link/pathway/ko")
    rows = [l.split("\t") for l in t.splitlines() if "\tpath:ko" in l]
    save("ko_pathway.tsv", "ko\tpathway\n" + "\n".join(f"{a[3:]}\t{b[7:]}" for a, b in rows))
    br = get("get/br:br08901"); A = B = None; h = []
    for line in br.splitlines():
        if line.startswith("A"): A = re.sub(r"<.*?>", "", line[1:].strip())
        elif line.startswith("B"): B = line[1:].strip()
        elif line.startswith("C"):
            m = re.match(r"C\s+(\d{5})\s+(.*)", line)
            if m: h.append((m.group(1), A, B, m.group(2)))
    save("pathway_hierarchy.tsv", "pathway\tlevel_a\tlevel_b\tname\n" + "\n".join("\t".join(x) for x in h))
    t = get("link/module/ko")
    save("ko_module.tsv", "ko\tmodule\n" + "\n".join(f"{a[3:]}\t{b[3:]}" for a, b in (l.split("\t") for l in t.splitlines())))
    mods = [l.split("\t")[0].replace("md:", "") for l in get("list/module").splitlines()]
    defs = []
    for i in range(0, len(mods), 10):
        txt = get("get/" + "+".join(mods[i:i + 10]))
        for entry in txt.split("///"):
            m = re.search(r"^ENTRY\s+(M\d+)", entry, re.M)
            if not m:
                continue
            name = re.search(r"^NAME\s+(.*)$", entry, re.M)
            cls = re.search(r"^CLASS\s+(.*)$", entry, re.M)
            d = re.search(r"^DEFINITION\s+(.*?)(?=^\S)", entry, re.M | re.S)
            defs.append((m.group(1), name.group(1).strip() if name else "", cls.group(1).strip() if cls else "",
                         " ".join(d.group(1).split()) if d else ""))
        time.sleep(0.4)
    save("module_definition.tsv", "module\tname\tclass\tdefinition\n" + "\n".join("\t".join(x) for x in defs))
    t = get("link/pathway/module")
    save("module_pathway.tsv", "module\tpathway\n" + "\n".join(f"{a[3:]}\t{b[8:]}" for a, b in (l.split("\t") for l in t.splitlines()) if b.startswith("path:map")))
    return genomes(n_genomes)


def genomes(n_genomes):
    # list/organism no longer answers (HTTP 400, 2026-10): organisms and lineages from the BRITE hierarchy br08601
    rows, path = [], {}
    for line in get("get/br:br08601").splitlines():
        if not line or line[0] not in "ABCDEFGH":
            continue
        level = "ABCDEFGH".index(line[0]); text = re.sub(r"<.*?>", "", line[1:].strip())
        m = re.match(r"^([a-z]{3,5})\s{2,}(.+)$", text)
        if m and level >= 3:
            rows.append((m.group(1), m.group(2), ";".join(path[i] for i in sorted(path) if i < level)))
        else:
            path = {i: v for i, v in path.items() if i < level}; path[level] = re.sub(r"\s*\(\d+\)$", "", text)
    o = pd.DataFrame(rows, columns=["org", "name", "lineage"])
    o = o[o.lineage.str.startswith("Prokaryotes")]
    o.to_csv(OUT / "organisms.tsv", sep="\t", index=False)
    print(f"prokaryotic organisms: {len(o):,}", flush=True)
    o["genus"] = o.name.str.split().str[0]
    sel = o.sort_values(["genus", "org"]).drop_duplicates("genus")
    sel = sel.sample(n=min(n_genomes, len(sel)), random_state=1).sort_values("org")
    sel.to_csv(OUT / "genomes_selected.tsv", sep="\t", index=False)
    gd = OUT / "genome_ko"; gd.mkdir(exist_ok=True)
    for k, org in enumerate(sel.org, 1):
        f = gd / f"{org}.tsv"
        if not f.exists():
            t = get(f"link/ko/{org}")
            f.write_text("\n".join(f"{a}\t{b[3:]}" for a, b in (l.split("\t") for l in t.splitlines())))
            time.sleep(0.4)
        if k % 50 == 0:
            print(f"genomes {k}/{len(sel)}", flush=True)
    print("done:", OUT)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
