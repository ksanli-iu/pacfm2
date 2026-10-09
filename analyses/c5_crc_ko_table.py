#!/usr/bin/env python3
"""Case study 2, step 1: KO x sample table for the colorectal cancer cohort of Feng et al. 2015 (FengQ_2015).

Input: HUMAnN 3 unstratified UniRef90 gene family abundances (relative abundance) from curatedMetagenomicData
(ExperimentHub EH5525, exported to $PACFM_WORK/feng/gene_families_unstratified.tsv.gz) and the curated
sample metadata (sampleMetadata of curatedMetagenomicData). UniRef90 families are mapped to KOs with the HUMAnN utility
mapping map_ko_uniref90 (v201901b); a family mapped to several KOs counts for each.
Output: data/case_crc/ko_table.tsv (KO x sample), data/case_crc/groups.tsv (sample, group: CRC / control; adenomas
are left out of the comparison but kept in the table)
"""
import gzip
from pathlib import Path

import pandas as pd

import os
W = Path(os.environ.get("PACFM_WORK", "work")) / "feng"
MAP = Path(os.environ.get("HUMANN_KO_MAP", "map_ko_uniref90.txt.gz"))   # HUMAnN utility mapping v201901b
OUT = Path(__file__).resolve().parents[1] / "data" / "case_crc"
OUT.mkdir(parents=True, exist_ok=True)

pairs = []
with gzip.open(MAP, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        pairs += [(u, f[0]) for u in f[1:]]
m = pd.DataFrame(pairs, columns=["family", "ko"])
gf = pd.read_csv(W / "gene_families_unstratified.tsv.gz", sep="\t")
print(f"gene family records: {len(gf):,}; families mapped to a KO: {gf.family.isin(m.family).mean():.1%}")
x = gf.merge(m, on="family")
table = x.pivot_table(index="ko", columns="sample", values="value", aggfunc="sum", fill_value=0)
table.to_csv(OUT / "ko_table.tsv", sep="\t")
meta = pd.read_csv(W / "FengQ_2015_metadata.tsv", sep="\t")
meta = meta[meta.sample_id.isin(table.columns)]
meta.rename(columns={"sample_id": "sample", "study_condition": "group"})[["sample", "group"]].to_csv(OUT / "groups.tsv", sep="\t", index=False)
print(f"KO table: {table.shape[0]:,} KOs x {table.shape[1]} samples; groups: {meta.study_condition.value_counts().to_dict()}")
print("KO share of total gene family abundance per sample (median): %.1f %%" % (100 * table.sum().median()))
