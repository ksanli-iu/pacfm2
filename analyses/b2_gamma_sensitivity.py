#!/usr/bin/env python3
"""Sensitivity of the coverage-weighted schemes to the exponent gamma (tests 1 and 2 of b1_insilico_benchmark.py:
spurious share and presence AUC, same communities and seed). Output: analyses/out/b2_gamma.tsv
Run: python analyses/b2_gamma_sensitivity.py [replicates]
"""
import sys

import numpy as np
import pandas as pd

import b1_insilico_benchmark as b1


def main(reps=20):
    kegg = b1.Kegg.load(b1.KDIR)
    C, T = b1.genome_tables(kegg)
    size = kegg.ko_pathway.groupby("pathway").ko.nunique()
    hz = kegg.hierarchy
    euk = set(hz.pathway[hz.level_a.isin(b1.EUK) & ~hz.level_b.str.contains("Infectious disease: bacterial|Drug resistance: antimicrobial", na=False)])
    rng = np.random.default_rng(1)
    rows = []
    for r in range(reps):
        gs = list(rng.choice(C.columns, b1.N_GENOMES, replace=False))
        X, _ = b1.community(C, rng, gs)
        X = X.to_frame("c")
        present = T[gs].any(axis=1)
        for scheme in ("fraccov", "emcov"):
            for gamma in (0.5, 1, 2, 4, 8):
                A = b1.pathway_abundance(X, kegg, scheme, gamma=gamma).c
                lab = present.reindex(A.index).dropna()
                score = (A / size.reindex(A.index)).reindex(lab.index).fillna(0)
                rows.append(dict(rep=r, scheme=scheme, gamma=gamma, spurious_share=A[A.index.isin(euk)].sum() / A.sum(),
                                 auc=b1.roc_auc_score(lab.astype(int), score)))
    t = pd.DataFrame(rows)
    t.to_csv(b1.OUT / "b2_gamma.tsv", sep="\t", index=False)
    print(t.groupby(["scheme", "gamma"])[["spurious_share", "auc"]].median().round(4).to_string())


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
