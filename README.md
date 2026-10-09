# PACFM 2 — shared-enzyme-aware pathway analysis for functional metagenomics

Most KEGG orthologues (KOs) take part in several pathways. Summing the KOs of each pathway ("naive" profiling) credits
widely shared enzymes to every pathway that contains them, including pathways that cannot operate in a microbial
community (e.g. human disease pathways in a gut or soil metagenome). PACFM 2 quantifies and corrects this bias:

- **six counting schemes**: naive, association cutoff (PACFM 2016), fractional, evidence-proportional EM, and
  coverage-weighted fractional / EM (recommended: `fraccov`);
- **evidence per pathway**: KO coverage, coverage of pathway-specific KOs, share of abundance from specific KOs;
- **evidence-gated group tests**: Mann–Whitney U + Benjamini–Hochberg, reporting a pathway only when it is supported by
  its own KOs (or ≥ 50 % KO coverage), and flagging results that are significant under naive summation only;
- **web application**: radial view of the KEGG hierarchy with pathway bars, coloured sheets joining pathways that share
  KOs, in-figure legend, full-screen and full-width views, SVG/TSV export.

Paper: Şanlı K. *Shared enzymes inflate pathway profiles of metagenomes: specificity-aware counting, evidence gating and
interactive visualization with PACFM 2* (submitted). Archived release: Zenodo DOI to be added.

## Python package

```bash
pip install .            # numpy, pandas
```

```python
import pandas as pd
from pacfm2.core import Kegg, pathway_abundance, pathway_evidence

kegg = Kegg.load("data/kegg/2026-10-07")              # built with data/kegg/fetch_kegg.py
ko = pd.read_csv("ko_table.tsv", sep="\t", index_col=0)   # KO x sample
profile = pathway_abundance(ko, kegg, scheme="fraccov")   # pathway x sample
evidence = pathway_evidence(ko, kegg)                     # coverage, unique_coverage, unique_share
```

## KEGG data

KEGG data are subject to the [KEGG licence](https://www.kegg.jp/kegg/legal.html) and are **not distributed** here.
Build your own snapshot through the KEGG REST API:

```bash
python data/kegg/fetch_kegg.py                         # -> data/kegg/<date>/ (KO-pathway links, hierarchy, modules, 300 genomes)
python webapp/build_kegg_reference.py                  # -> webapp/pacfm-web/public/kegg_reference.json (for the web app)
```

## Web application

```bash
cd webapp/pacfm-web
npm install
npm run dev              # local development
npm run build            # static site in dist/ (copy kegg_reference.json next to index.html, or load it in the app)
```

The app reads KO tables in plain tabular form, HUMAnN gene family / KO tables, eggNOG-mapper annotations and KofamScan
detail output; in the standalone build, all data stay in the browser. `?demo=1` loads the demonstration data.
`webapp/test/` checks that the TypeScript and Python implementations give identical results.

## Reproducing the paper

| step | script | output |
|---|---|---|
| module parser vs KEGG | `analyses/v0_module_parser_vs_kegg.py` | `results/v0_parser_vs_kegg.tsv` |
| in silico benchmark | `analyses/b1_insilico_benchmark.py` | `results/b1_*` |
| sensitivity to γ | `analyses/b2_gamma_sensitivity.py` | `results/b2_gamma.tsv` |
| piglet case study | `analyses/c1_piglets_cluster.sh`, `c2_piglets_kofam.sh`, `c3_piglets_ko_table.py` | KO table |
| colorectal cancer case study | `analyses/c5_crc_ko_table.py` | KO table |
| group comparisons | `analyses/c4_case_study.py <case> <group_a> <group_b>` | `results/case_*` |
| profile shifts | `analyses/c6_profile_shift.py` | `results/c6_*` |
| figures | `analyses/make_figures.py` | `figures/` |

The case-study scripts read inputs from environment variables: `AMRSURV_PROJECT` (an AMRsurv run with MEGAHIT contigs
and Prodigal proteins), `KOFAM_DB` (KOfam profiles and ko_list), `HUMANN_KO_MAP` (HUMAnN utility mapping
`map_ko_uniref90`), `PACFM_WORK` (scratch folder) and `TOOLS_BIN` (MMseqs2, KofamScan).

## Licence

GPL-3.0 (see `LICENSE`).
