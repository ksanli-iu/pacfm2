#!/usr/bin/env python3
"""Validate the KEGG module parser of pacfm2.core against KEGG's own calls: each module entry lists the organisms in
which KEGG considers the module complete (COMPLETE field). For a random set of modules and the 300 benchmark genomes,
compare parser completeness (= 1) with membership in KEGG's COMPLETE list.
Run: python analyses/v0_module_parser_vs_kegg.py [n_modules]
"""
import re, sys, time, urllib.request
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from pacfm2.core import Kegg, module_completeness  # noqa: E402
KDIR = sorted((ROOT / "data" / "kegg").glob("20*"))[-1]


def complete_orgs(mod):
    t = urllib.request.urlopen(f"https://rest.kegg.jp/get/{mod}", timeout=120).read().decode()
    m = re.search(r"^COMPLETE\s+(.*?)(?=^\S)", t, re.M | re.S)
    return {l.split()[0] for l in m.group(1).splitlines() if l.strip()} if m else set()


def main(n=60):
    kegg = Kegg.load(KDIR)
    genomes = {f.stem: set(pd.read_csv(f, sep="\t", header=None)[1]) for f in sorted((KDIR / "genome_ko").glob("*.tsv"))}
    comp = {g: module_completeness(kegg, kos) for g, kos in genomes.items()}
    mods = kegg.module_def.sample(n=n, random_state=2).module
    rows = []
    for m in mods:
        ko = complete_orgs(m); time.sleep(0.3)
        for g in genomes:
            rows.append(dict(module=m, genome=g, parser=bool(comp[g].get(m, 0) == 1.0), kegg=g in ko))
    d = pd.DataFrame(rows); d.to_csv(ROOT / "analyses" / "out" / "v0_parser_vs_kegg.tsv", sep="\t", index=False)
    t = pd.crosstab(d.parser, d.kegg); print(t)
    agree = (d.parser == d.kegg).mean(); print(f"agreement {100 * agree:.2f} % over {len(d):,} module x genome pairs")
    bad = d[d.parser != d.kegg].groupby("module").size().sort_values(ascending=False)
    print("modules with disagreements:"); print(bad.head(10).to_string())


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
