#!/usr/bin/env python3
"""Demo data for the web app: 2 x 8 simulated communities from the KEGG benchmark genomes (synthetic, shareable).
Group B is enriched (x8) for genomes in which a methanogenesis module (M00567, CO2 => methane) is complete.
Output: webapp/pacfm-web/public/demo_ko_table.tsv (KO x sample counts), demo_groups.tsv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from pacfm2.core import Kegg, module_completeness  # noqa: E402
KDIR = sorted((ROOT / "data" / "kegg").glob("20*"))[-1]
kegg = Kegg.load(KDIR)
C = pd.DataFrame({f.stem: pd.read_csv(f, sep="\t", header=None)[1].value_counts() for f in sorted((KDIR / "genome_ko").glob("*.tsv"))}).fillna(0)
meth = [g for g in C.columns if module_completeness(kegg, set(C.index[C[g] > 0])).get("M00567", 0) == 1.0]
print("genomes with complete M00567:", len(meth))
rng = np.random.default_rng(7); cols, grp = {}, []
for i in range(16):
    gs = list(rng.choice(C.columns, 25, replace=False))
    if i >= 8: gs = list(dict.fromkeys(gs[:22] + list(rng.choice(meth, 3, replace=False))))
    w = pd.Series(rng.lognormal(0, 1.5, len(gs)), index=gs)
    if i >= 8: w[w.index.isin(meth)] *= 8
    p = (C[gs] * (w / w.sum())).sum(axis=1); p = p[p > 0] / p.sum()
    name = f"{'A' if i < 8 else 'B'}{i % 8 + 1}"
    cols[name] = pd.Series(rng.multinomial(1_000_000, p.values), index=p.index); grp.append((name, "A" if i < 8 else "B"))
X = pd.DataFrame(cols).fillna(0).astype(int); X = X[X.sum(axis=1) > 0]
pub = ROOT / "webapp" / "pacfm-web" / "public"
X.to_csv(pub / "demo_ko_table.tsv", sep="\t", index_label="KO")
pd.DataFrame(grp, columns=["sample", "group"]).to_csv(pub / "demo_groups.tsv", sep="\t", index=False)
print(X.shape, "->", pub / "demo_ko_table.tsv")
