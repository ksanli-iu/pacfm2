#!/usr/bin/env bash
# Case study 1: pool the Prodigal proteins of the 20 piglet metagenomes (PRJNA737271; antibiotic-free vs conventional,
# assembled with AMRsurv in PROJ0022), cluster them with MMseqs2 (90 % identity, 80 % coverage) so that KofamScan only
# annotates one representative per cluster; KOs are passed on to all cluster members afterwards.
set -euo pipefail
P=$(cd "$(dirname "$0")/.." && pwd); A=${AMRSURV_PROJECT:?set AMRSURV_PROJECT to the AMRsurv run folder}
W=${PACFM_WORK:-$P/work}/piglets; E=${TOOLS_BIN:-$(dirname "$(command -v mmseqs)")}
mkdir -p "$W"
awk -F'\t' '$2=="PRJNA737271"{print $1}' $A/metadata/expansion_samples.tsv > "$W/runs.txt"
: > "$W/all.faa"
while read -r r; do
  zcat "$A/results_expansion/prodigal/MEGAHIT-${r}_prodigal.faa.gz" | sed "s/^>\(\S*\).*/>${r}|\1/; /^[^>]/s/\*$//" >> "$W/all.faa"
done < "$W/runs.txt"
grep -c '>' "$W/all.faa"
"$E/mmseqs" easy-linclust "$W/all.faa" "$W/clu" "$W/tmp" --min-seq-id 0.9 -c 0.8 --cov-mode 1 --threads 48 > "$W/linclust.log" 2>&1
grep -c '>' "$W/clu_rep_seq.fasta"
