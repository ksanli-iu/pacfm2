#!/usr/bin/env python3
"""Case studies: compare two groups of metagenomes under every counting scheme, with and without the evidence gate.

For each scheme, pathway abundances are converted to relative abundances per sample and compared with the Mann-Whitney
U test (two-sided) and Benjamini-Hochberg adjustment. The evidence gate (as in the web app and the benchmark): a
pathway-specific KO detected in at least half of the samples, or median KO coverage >= 0.5.
Outputs (analyses/out/<case>_*):
  _pathways.tsv   pathway x scheme: q-values, relative abundances, evidence and gate
  _summary.tsv    per scheme: significant pathways, of which in Human Diseases / Organismal Systems, naive-only
Run: python analyses/c4_case_study.py <case> <group_a> <group_b>
     e.g. c4_case_study.py case_piglets RWA conventional ; c4_case_study.py case_crc CRC control
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pacfm2.core import Kegg, pathway_abundance, pathway_evidence  # noqa: E402

SCHEMES = [("naive", "naive", {}), ("cutoff1", "cutoff", {"cutoff": 1}), ("cutoff2", "cutoff", {"cutoff": 2}),
           ("fractional", "fractional", {}), ("em", "em", {}), ("fraccov", "fraccov", {}), ("emcov", "emcov", {})]
EUK = {"Human Diseases", "Organismal Systems"}


def main(case, ga, gb):
    d = ROOT / "data" / case
    X = pd.read_csv(d / "ko_table.tsv", sep="\t", index_col=0)
    g = pd.read_csv(d / "groups.tsv", sep="\t").set_index(pd.read_csv(d / "groups.tsv", sep="\t").columns[0]).iloc[:, 0]
    g = g[g.isin([ga, gb])]
    X = X[g.index.intersection(X.columns)]
    g = g.reindex(X.columns)
    kegg = Kegg.load(sorted((ROOT / "data" / "kegg").glob("20*"))[-1])
    hz = kegg.hierarchy.set_index("pathway")
    a, b = (g == ga).values, (g == gb).values
    print(f"{case}: {X.shape[0]:,} KOs, {a.sum()} {ga} vs {b.sum()} {gb}; "
          f"{X.index.isin(kegg.n_assoc.index).mean():.1%} of KOs in KEGG pathways")
    ev = pathway_evidence(X, kegg)
    gate = ((ev["unique_coverage"].fillna(0) > 0).mean(axis=1) >= 0.5) | (ev["coverage"].median(axis=1) >= 0.5)
    res = {}
    for name, scheme, kw in SCHEMES:
        A = pathway_abundance(X, kegg, scheme, **kw)
        A = A.div(A.sum(axis=0), axis=1)
        p = A.apply(lambda v: mannwhitneyu(v[a], v[b]).pvalue if v.std() > 0 else 1.0, axis=1)
        q = pd.Series(multipletests(p, method="fdr_bh")[1], index=A.index)
        res[name] = pd.DataFrame({"q": q, f"mean_{ga}": A.loc[:, a].mean(axis=1), f"mean_{gb}": A.loc[:, b].mean(axis=1)})
    out = pd.concat(res, axis=1)
    out.columns = [f"{s}_{c}" for s, c in out.columns]
    out = out.join(hz[["level_a", "level_b", "name"]])
    out["median_coverage"] = ev["coverage"].median(axis=1)
    out["specific_ko_detected_share"] = (ev["unique_coverage"].fillna(0) > 0).mean(axis=1)
    out["gate"] = gate.reindex(out.index).fillna(False)
    out["fraccov_gated_sig"] = (out.fraccov_q < 0.05) & out.gate
    out["naive_only"] = (out.naive_q < 0.05) & ~out.fraccov_gated_sig
    out.to_csv(ROOT / "analyses" / "out" / f"{case}_pathways.tsv", sep="\t")
    rows = []
    for name, *_ in SCHEMES + [("fraccov+gate",)]:
        sig = out.fraccov_gated_sig if name == "fraccov+gate" else out[f"{name}_q"] < 0.05
        rows.append(dict(scheme=name, significant=int(sig.sum()), absent_category=int((sig & out.level_a.isin(EUK)).sum()),
                         metabolism=int((sig & (out.level_a == "Metabolism")).sum())))
    s = pd.DataFrame(rows).set_index("scheme")
    s.to_csv(ROOT / "analyses" / "out" / f"{case}_summary.tsv", sep="\t")
    print(s.to_string())
    show = out[out.naive_q < 0.05].sort_values("naive_q")
    print("\nsignificant with naive summation (q < 0.05), first 25:")
    print(show[["name", "level_a", "naive_q", "fraccov_q", "gate", "naive_only"]].head(25).to_string())
    print("\nfraccov + gate, first 25:")
    print(out[out.fraccov_gated_sig].sort_values("fraccov_q")[["name", "level_a", "naive_q", "fraccov_q", "fraccov_mean_" + ga, "fraccov_mean_" + gb]].head(25).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
