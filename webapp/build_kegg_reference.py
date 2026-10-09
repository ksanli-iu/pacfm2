#!/usr/bin/env python3
"""Build the KEGG reference file for the PACFM 2 web app from a KEGG snapshot (data/kegg/<date>/).
Output: webapp/pacfm-web/public/kegg_reference.json  {"version", "kp": [[ko, pathway], ...], "h": [[pathway, A, B, name], ...]}
KEGG data are for academic use; a public deployment should ship this file only where KEGG's terms allow, otherwise users
build it themselves with this script (it uses only the KEGG REST API).
Run: python webapp/build_kegg_reference.py [snapshot_dir]
"""
import json, re, sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
snap = Path(sys.argv[1]) if len(sys.argv) > 1 else sorted((ROOT / "data" / "kegg").glob("20*"))[-1]
kp = pd.read_csv(snap / "ko_pathway.tsv", sep="\t", dtype=str)
h = pd.read_csv(snap / "pathway_hierarchy.tsv", sep="\t", dtype=str)
kp = kp[~kp.pathway.str.match(r"^01[12]\d\d$") & kp.pathway.isin(h.pathway)].drop_duplicates()
out = ROOT / "webapp" / "pacfm-web" / "public" / "kegg_reference.json"
out.write_text(json.dumps({"version": f"KEGG REST snapshot {snap.name}", "kp": kp.values.tolist(),
                           "h": h[h.pathway.isin(kp.pathway)].fillna("").values.tolist()}, separators=(",", ":")))
print(out, f"{out.stat().st_size / 1e6:.2f} MB", kp.pathway.nunique(), "pathways", kp.ko.nunique(), "KOs")
