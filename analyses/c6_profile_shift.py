#!/usr/bin/env python3
"""Case studies, profiling side: per sample, the share of the pathway signal assigned to Human Diseases and Organismal
Systems pathways (absent from bacterial communities, except bacterial infectious disease and antimicrobial
resistance) under each counting scheme, and the pathways whose share changes most between naive summation and
fraccov. Output: analyses/out/c6_profile_shift.tsv, c6_top_changes_<case>.tsv
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pacfm2.core import Kegg, pathway_abundance  # noqa: E402

kegg = Kegg.load(sorted((ROOT / "data" / "kegg").glob("20*"))[-1])
hz = kegg.hierarchy.set_index("pathway")
absent = set(hz.index[hz.level_a.isin(["Human Diseases", "Organismal Systems"])
                      & ~hz.level_b.str.contains("Infectious disease: bacterial|Drug resistance: antimicrobial", na=False)])
rows = []
for case in ("case_piglets", "case_crc"):
    X = pd.read_csv(ROOT / "data" / case / "ko_table.tsv", sep="\t", index_col=0)
    prof = {}
    for name, scheme, kw in [("naive", "naive", {}), ("cutoff2", "cutoff", {"cutoff": 2}), ("fractional", "fractional", {}),
                             ("em", "em", {}), ("fraccov", "fraccov", {}), ("emcov", "emcov", {})]:
        A = pathway_abundance(X, kegg, scheme, **kw)
        A = A.div(A.sum(axis=0), axis=1)
        prof[name] = A
        sh = A[A.index.isin(absent)].sum()
        rows.append(dict(case=case, scheme=name, median_absent_share=sh.median(), min=sh.min(), max=sh.max()))
    n, f = prof["naive"].mean(axis=1), prof["fraccov"].reindex(prof["naive"].index).fillna(0).mean(axis=1)
    top = pd.DataFrame({"name": hz.name.reindex(n.index), "level_a": hz.level_a.reindex(n.index), "naive": n, "fraccov": f})
    top["ratio"] = top.fraccov / top.naive
    top = top[top.naive > 0.002].sort_values("ratio")
    top.to_csv(ROOT / "analyses" / "out" / f"c6_top_changes_{case}.tsv", sep="\t")
    print(case, "most reduced:\n", top.head(8).round(4).to_string(), "\nmost increased:\n", top.tail(6).round(4).to_string())
t = pd.DataFrame(rows)
t.to_csv(ROOT / "analyses" / "out" / "c6_profile_shift.tsv", sep="\t", index=False)
print(t.round(4).to_string())
