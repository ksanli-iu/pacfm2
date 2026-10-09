#!/usr/bin/env python3
"""Case study 1, step 3: KO x sample abundance table for the 20 piglet metagenomes (PRJNA737271).

The KO of each MMseqs2 cluster representative (KofamScan, c2_piglets_kofam.sh) is passed on to all members of its
cluster; the abundance of a gene is the k-mer depth of its contig (MEGAHIT 'multi='), and the abundance of a KO in a
sample is the sum over its genes. Genes with more than one KO assignment count once for each.
Output: data/case_piglets/ko_table.tsv (KO x run), data/case_piglets/groups.tsv (run, group)
"""
import gzip
import re
from pathlib import Path

import pandas as pd

import os
W = Path(os.environ.get("PACFM_WORK", "work")) / "piglets"
A = Path(os.environ["AMRSURV_PROJECT"])   # AMRsurv run: metadata/expansion_samples.tsv, results_expansion/
OUT = Path(__file__).resolve().parents[1] / "data" / "case_piglets"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(A / "metadata" / "expansion_samples.tsv", sep="\t")
meta = meta[meta.bioproject == "PRJNA737271"]
meta[["run", "group"]].to_csv(OUT / "groups.tsv", sep="\t", index=False)

rep_ko = pd.read_csv(W / "rep_ko.tsv", sep="\t", header=None, names=["rep", "ko"])
clu = pd.read_csv(W / "clu_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
genes = clu.merge(rep_ko, on="rep")[["member", "ko"]]
print(f"cluster members with a KO: {genes.member.nunique():,} of {clu.member.nunique():,}")
genes["run"] = genes.member.str.split("|").str[0]
genes["contig"] = genes.member.str.split("|").str[1].str.replace(r"_\d+$", "", regex=True)

depth = []
pat = re.compile(r"^>(\S+).*?multi=([\d.]+)")
for run in meta.run:
    with gzip.open(A / "results_expansion" / "Assembly" / "MEGAHIT" / f"MEGAHIT-{run}.contigs.fa.gz", "rt") as fh:
        for line in fh:
            if line[0] == ">":
                m = pat.match(line)
                depth.append((run, m.group(1), float(m.group(2))))
depth = pd.DataFrame(depth, columns=["run", "contig", "depth"])
genes = genes.merge(depth, on=["run", "contig"], how="left")
print(f"genes without contig depth: {genes.depth.isna().sum():,}")
table = genes.pivot_table(index="ko", columns="run", values="depth", aggfunc="sum", fill_value=0)
table.index.name = "ko"
table.to_csv(OUT / "ko_table.tsv", sep="\t")
print(f"KO table: {table.shape[0]:,} KOs x {table.shape[1]} samples -> {OUT / 'ko_table.tsv'}")
