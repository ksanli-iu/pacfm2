#!/usr/bin/env python3
"""In silico benchmark of pathway counting schemes (PACFM 2).

Communities are assembled from KEGG prokaryotic reference genomes (data/kegg/<date>/genome_ko), so the true metabolic
content of every community is known.
  KO counts      a_k = sum_g w_g * copies_gk, sampled as a multinomial of DEPTH gene hits (sequencing noise)
  truth          a pathway is present in a genome when at least one of its KEGG modules is complete (all steps); for
                 pathways without modules, truth is undefined except for the Human Diseases / Organismal Systems
                 pathways, which are absent from prokaryotic communities by definition
Schemes: naive, cutoff c in {1, 2, 3, 10} (PACFM 2016), fractional, em, fraccov, emcov.
Evidence-gated testing (test 3 only): abundances from fraccov, and a pathway is called differential only when it also
passes an evidence gate computed from the KO table: gate_unique = at least one pathway-unique KO detected in >= half of
the samples; gate_cov = median KO coverage >= 0.5; gate_any = either.

Tests
  1 spurious signal   share of total pathway abundance assigned to Human Diseases + Organismal Systems pathways
  2 presence          per community, AUC of size-normalised pathway abundance for present vs absent pathways
                      (module-linked pathways only)
  3 differential      2 x 10 communities; group B enriches genomes carrying target modules; true pathway abundance =
                      sum_g w_g * [pathway present in g]; Mann-Whitney per pathway, BH q < 0.05; precision, recall,
                      false positives among Human Diseases / Organismal Systems pathways
Outputs: analyses/out/b1_*.tsv
Run: python analyses/b1_insilico_benchmark.py [replicates]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def roc_auc_score(y, score):
    """AUC = Mann-Whitney U / (n_pos * n_neg)"""
    y = np.asarray(y).astype(bool); score = np.asarray(score, float)
    return mannwhitneyu(score[y], score[~y], alternative="two-sided").statistic / (y.sum() * (~y).sum())

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pacfm2.core import pathway_evidence, Kegg, pathway_abundance, module_completeness  # noqa: E402

KDIR = sorted((ROOT / "data" / "kegg").glob("20*"))[-1]
OUT = ROOT / "analyses" / "out"
SCHEMES = [("naive", {}), ("cutoff1", {"cutoff": 1}), ("cutoff2", {"cutoff": 2}), ("cutoff3", {"cutoff": 3}),
           ("cutoff10", {"cutoff": 10}), ("fractional", {}), ("em", {}), ("fraccov", {}), ("emcov", {})]
EUK = {"Human Diseases", "Organismal Systems"}
DEPTH, N_GENOMES, N_PER_GROUP = 2_000_000, 25, 10


def scheme_abundance(X, kegg, name, kw):
    s = "cutoff" if name.startswith("cutoff") else name
    return pathway_abundance(X, kegg, s, **kw)


def genome_tables(kegg):
    copies = {}
    for f in sorted((KDIR / "genome_ko").glob("*.tsv")):
        g = pd.read_csv(f, sep="\t", header=None, names=["gene", "ko"])
        copies[f.stem] = g.ko.value_counts()
    C = pd.DataFrame(copies).fillna(0)                       # KO x genome
    modpath = kegg.module_pathway.groupby("pathway").module.apply(set)
    truth = {}
    for g in C.columns:
        comp = module_completeness(kegg, set(C.index[C[g] > 0]))
        complete = set(comp.index[comp == 1.0])
        truth[g] = {p: bool(ms & complete) for p, ms in modpath.items()}
    T = pd.DataFrame(truth).fillna(False)                    # pathway (module-linked) x genome
    return C, T


def community(C, rng, genomes, boost=None, boost_factor=8.0):
    w = pd.Series(rng.lognormal(0, 1.5, len(genomes)), index=genomes)
    if boost is not None:
        w[w.index.isin(boost)] *= boost_factor
    w /= w.sum()
    p = (C[genomes] * w).sum(axis=1)
    p = p[p > 0] / p.sum()
    counts = rng.multinomial(DEPTH, p.values)
    return pd.Series(counts, index=p.index), w


def main(reps=20):
    kegg = Kegg.load(KDIR)
    C, T = genome_tables(kegg)
    T.to_csv(OUT / "b1_genome_pathway_truth.tsv", sep="\t")
    size = kegg.ko_pathway.groupby("pathway").ko.nunique()
    # truly absent from prokaryotic communities: Organismal Systems and Human Diseases, except the categories that hold
    # bacterial genes (bacterial infectious disease: pathogenicity factors; antimicrobial drug resistance)
    hz = kegg.hierarchy
    euk = set(hz.pathway[hz.level_a.isin(EUK) & ~hz.level_b.str.contains("Infectious disease: bacterial|Drug resistance: antimicrobial", na=False)])
    rng = np.random.default_rng(1)
    rows1, rows2, rows3 = [], [], []
    for r in range(reps):
        # ---- tests 1 and 2: single communities
        gs = list(rng.choice(C.columns, N_GENOMES, replace=False))
        X, w = community(C, rng, gs)
        X = X.to_frame("c")
        present = T[gs].any(axis=1)
        for name, kw in SCHEMES:
            A = scheme_abundance(X, kegg, name, kw).c
            rows1.append(dict(rep=r, scheme=name, spurious_share=A[A.index.isin(euk)].sum() / A.sum()))
            lab = present.reindex(A.index).dropna()
            if lab.nunique() == 2:
                score = (A / size.reindex(A.index)).reindex(lab.index).fillna(0)
                rows2.append(dict(rep=r, scheme=name, auc=roc_auc_score(lab.astype(int), score)))
        # ---- test 3: differential abundance. Cohort-like design: both groups share one core community (30 genomes,
        # fixed base abundances); each sample varies its abundances (log-normal, sd 0.5); in group B the genomes that
        # carry any of 3 target modules are 4x more abundant.
        target = rng.choice(T.index[(T.mean(axis=1) > 0.1) & (T.mean(axis=1) < 0.4)], 3, replace=False)
        core = list(rng.choice(C.columns, 30, replace=False))
        carriers = [g for g in core if T.loc[target, g].any()] if all(g in T.columns for g in core) else []
        base = pd.Series(rng.lognormal(0, 1.5, len(core)), index=core)
        Xs, Ws, grp = [], [], []
        for i in range(2 * N_PER_GROUP):
            w = base * rng.lognormal(0, 0.5, len(core))
            if i >= N_PER_GROUP:
                w[w.index.isin(carriers)] *= 4
            w = w / w.sum()
            pk = (C[core] * w).sum(axis=1); pk = pk[pk > 0] / pk.sum()
            Xs.append(pd.Series(rng.multinomial(DEPTH, pk.values), index=pk.index).rename(i)); Ws.append(w.rename(i))
            grp.append(int(i >= N_PER_GROUP))
        X = pd.concat(Xs, axis=1).fillna(0)
        W = pd.concat(Ws, axis=1).fillna(0)
        truth_ab = T.reindex(columns=W.index, fill_value=False).astype(float) @ W     # pathway x community
        truth_ab = truth_ab.div(truth_ab.sum(axis=0), axis=1)
        grp = np.array(grp)

        def de(M):
            p = M.apply(lambda v: mannwhitneyu(v[grp == 0], v[grp == 1]).pvalue if v.std() > 0 else 1.0, axis=1)
            return pd.Series(multipletests(p, method="fdr_bh")[1], index=M.index) < 0.05
        true_de = de(truth_ab)
        ev = pathway_evidence(X, kegg)
        g_unique = (ev["unique_coverage"].fillna(0) > 0).mean(axis=1) >= 0.5
        g_cov = ev["coverage"].median(axis=1) >= 0.5
        gates = {"gate_unique": g_unique, "gate_cov": g_cov, "gate_any": g_unique | g_cov}
        calls = {}
        for name, kw in SCHEMES:
            A = scheme_abundance(X, kegg, name, kw)
            calls[name] = de(A.div(A.sum(axis=0), axis=1))
        for gname, gate in gates.items():
            calls[f"fraccov+{gname}"] = calls["fraccov"] & gate.reindex(calls["fraccov"].index, fill_value=False)
        for name, called in calls.items():
            mod = called.index.intersection(true_de.index)
            tp = int((called[mod] & true_de[mod]).sum()); fp = int((called[mod] & ~true_de[mod]).sum())
            fn = int((~called[mod] & true_de[mod]).sum())
            rows3.append(dict(rep=r, scheme=name, tp=tp, fp=fp, fn=fn,
                              precision=tp / (tp + fp) if tp + fp else np.nan, recall=tp / (tp + fn) if tp + fn else np.nan,
                              euk_false_positives=int(called[called.index.isin(euk)].sum())))
        print(f"rep {r + 1}/{reps}", flush=True)
    t1, t2, t3 = pd.DataFrame(rows1), pd.DataFrame(rows2), pd.DataFrame(rows3)
    t1.to_csv(OUT / "b1_spurious.tsv", sep="\t", index=False); t2.to_csv(OUT / "b1_presence_auc.tsv", sep="\t", index=False)
    t3.to_csv(OUT / "b1_differential.tsv", sep="\t", index=False)
    s = pd.concat([t1.groupby("scheme").spurious_share.median(), t2.groupby("scheme").auc.median(),
                   t3.groupby("scheme")[["precision", "recall", "euk_false_positives"]].median()], axis=1)
    s.to_csv(ROOT / "results" / "b1_benchmark_summary.tsv", sep="\t")
    print(s.round(3).to_string())
    tot = t3.groupby("scheme")[["tp", "fp", "fn", "euk_false_positives"]].sum()
    tot["precision"] = tot.tp / (tot.tp + tot.fp); tot["recall"] = tot.tp / (tot.tp + tot.fn)
    tot.to_csv(ROOT / "results" / "b1_differential_totals.tsv", sep="\t")
    print(tot.round(3).to_string())


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
