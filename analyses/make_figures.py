#!/usr/bin/env python3
"""Figures 1 and 2 and Table 2 of the PACFM 2 manuscript.
Run: python analyses/make_figures.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pacfm2.core import Kegg  # noqa: E402

FIG, TAB, OUT = ROOT / "figures", ROOT / "results", ROOT / "analyses" / "out"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
INK, MUTED, BLUE, ORANGE, GREEN = "#0b0b0b", "#898781", "#2a78d6", "#eb6834", "#1baf7a"

ORDER_A = ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes",
           "Organismal Systems", "Human Diseases"]
SHORT_A = ["Metabolism", "Genetic\ninform.", "Environ.\ninform.", "Cellular\nprocesses",
           "Organismal\nsystems", "Human\ndiseases"]
SCHEMES = ["naive", "cutoff1", "cutoff2", "cutoff3", "cutoff10", "fractional", "em", "fraccov", "emcov"]
LABEL = {"naive": "naive", "cutoff1": "cutoff 1", "cutoff2": "cutoff 2", "cutoff3": "cutoff 3", "cutoff10": "cutoff 10",
         "fractional": "fractional", "em": "EM", "fraccov": "fraccov", "emcov": "emcov",
         "fraccov+gate_any": "fraccov + gate", "fraccov+gate_unique": "fraccov + gate (specific KO)",
         "fraccov+gate_cov": "fraccov + gate (coverage)"}


def colour(s):
    return MUTED if s == "naive" else ORANGE if s.startswith("cutoff") else GREEN if "gate" in s else BLUE


def fig1():
    k = Kegg.load(sorted((ROOT / "data" / "kegg").glob("20*"))[-1])
    s = k.specificity()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1, 1.5]})
    n = k.n_assoc
    cnt = n.clip(upper=21).value_counts().sort_index()
    a.bar(cnt.index, cnt.values, color=BLUE, width=0.8)
    a.set_yscale("log"); a.set_xticks([1, 5, 10, 15, 21]); a.set_xticklabels(["1", "5", "10", "15", ">20"])
    a.set_xlabel("Pathways per KO"); a.set_ylabel("KOs")
    a.set_title(f"A  {100 * (n > 1).mean():.1f} % of {len(n):,} KOs are shared", loc="left", fontsize=8)
    rng = np.random.default_rng(1)
    for i, lvl in enumerate(ORDER_A):
        v = s.loc[s.level_a == lvl, "pct_unique"].values
        b.scatter(i + rng.uniform(-0.28, 0.28, len(v)), v, s=8, color=BLUE, alpha=0.55, linewidths=0)
        b.plot([i - 0.32, i + 0.32], [np.median(v)] * 2, color=INK, lw=1.5)
    b.set_xticks(range(len(ORDER_A))); b.set_xticklabels(SHORT_A, fontsize=7)
    b.set_ylabel("Pathway-specific KOs (%)")
    b.set_title("B  Specific KOs per pathway (line: median)", loc="left", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "fig1_specificity.png", dpi=300); fig.savefig(FIG / "fig1_specificity.pdf")


def fig2():
    sp = pd.read_csv(OUT / "b1_spurious.tsv", sep="\t"); auc = pd.read_csv(OUT / "b1_presence_auc.tsv", sep="\t")
    tot = pd.read_csv(TAB / "b1_differential_totals.tsv", sep="\t", index_col=0)
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.2, 2.9))
    for ax, t, col, ttl in [(a, sp, "spurious_share", "A  Signal in pathways absent\n    from prokaryotes (%)"),
                            (b, auc, "auc", "B  Present vs absent\n    pathways (AUC)")]:
        data = [t.loc[t.scheme == s, col].values * (100 if col == "spurious_share" else 1) for s in SCHEMES]
        bp = ax.boxplot(data, orientation="horizontal", widths=0.6, patch_artist=True, medianprops={"color": INK}, flierprops={"markersize": 2})
        for p, s in zip(bp["boxes"], SCHEMES):
            p.set_facecolor(colour(s)); p.set_alpha(0.6); p.set_edgecolor("none")
        ax.set_yticks(range(1, len(SCHEMES) + 1)); ax.set_yticklabels([LABEL[s] for s in SCHEMES]); ax.invert_yaxis()
        ax.set_title(ttl, loc="left", fontsize=8)
    b.set_yticklabels([])
    show = SCHEMES + ["fraccov+gate_any"]
    for s in show:
        r = tot.loc[s]
        c.scatter(r.euk_false_positives, r.recall, s=14 + 300 * (r.precision - 0.2), color=colour(s), alpha=0.85, linewidths=0)
        off = {"fractional": (-4, 5), "fraccov": (-4, -9), "naive": (-4, 5), "cutoff10": (-4, -10)}.get(s, (5, 2))
        c.annotate(LABEL[s], (r.euk_false_positives, r.recall), xytext=off, textcoords="offset points", fontsize=6.5,
                   ha="right" if s in ("fractional", "fraccov", "naive", "cutoff10") else "left")
    c.set_xlabel("Significant pathways absent\nfrom prokaryotes (20 replicates)"); c.set_ylabel("Recall")
    c.set_title("C  Differential test\n    (point size: precision)", loc="left", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "fig2_benchmark.png", dpi=300); fig.savefig(FIG / "fig2_benchmark.pdf")


def table2():
    med = pd.read_csv(TAB / "b1_benchmark_summary.tsv", sep="\t", index_col=0)
    tot = pd.read_csv(TAB / "b1_differential_totals.tsv", sep="\t", index_col=0)
    rows = []
    for s in SCHEMES + ["fraccov+gate_any", "fraccov+gate_unique", "fraccov+gate_cov"]:
        m, t = med.loc[s], tot.loc[s]
        sp = "" if pd.isna(m.spurious_share) else f"{100 * m.spurious_share:.1f}"
        au = "" if pd.isna(m.auc) else f"{m.auc:.3f}"
        rows.append(f"| {LABEL[s]} | {sp} | {au} | {t.precision:.3f} | {t.recall:.3f} | {int(t.euk_false_positives)} |")
    head = ("| Scheme | Spurious share (%) | Presence AUC | Precision | Recall | Significant absent pathways |\n"
            "|---|---:|---:|---:|---:|---:|")
    (TAB / "t2_benchmark.md").write_text(head + "\n" + "\n".join(rows) + "\n")
    print((TAB / "t2_benchmark.md").read_text())


if __name__ == "__main__":
    fig1(); fig2(); table2()
