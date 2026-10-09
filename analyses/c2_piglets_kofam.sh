#!/usr/bin/env bash
# Case study 1, step 2: KofamScan on the MMseqs2 cluster representatives (c1_piglets_cluster.sh), in chunks of 100,000
# proteins (adaptive thresholds, mapper output = only assignments above threshold), then KO x sample table by
# c3_piglets_ko_table.py. Timing test: 20,000 proteins in 51 s on 48 cores (~2.3 h for 3.16 M).
set -euo pipefail
P=$(cd "$(dirname "$0")/.." && pwd)
W=${PACFM_WORK:-$P/work}/piglets; K=${KOFAM_DB:?set KOFAM_DB to the KOfam folder (profiles/, ko_list)}; E=${TOOLS_BIN:-$(dirname "$(command -v exec_annotation)")}
export PATH=$E:$PATH
mkdir -p "$W/kofam_chunks"
cd "$W/kofam_chunks"
[ -s chunk_000.faa ] || awk 'BEGIN{n=-1} /^>/{n++; if (n % 100000 == 0) {if (f) close(f); f=sprintf("chunk_%03d.faa", n/100000)}} {print > f}' "$W/clu_rep_seq.fasta"
for f in chunk_*.faa; do
  o=${f%.faa}.ko.tsv
  [ -s "$o.done" ] && continue
  rm -rf tmp_$f
  exec_annotation -p $K/profiles/prokaryote.hal -k $K/ko_list --cpu 48 --tmp-dir tmp_$f -f mapper -o "$o" "$f" > /dev/null
  rm -rf tmp_$f; echo ok > "$o.done"; echo "$(date +%T) $f" >&2
done
cat chunk_*.ko.tsv | awk -F'\t' 'NF >= 2 && $2 != ""' > "$W/rep_ko.tsv"
echo "representatives with KO: $(wc -l < "$W/rep_ko.tsv")" >&2
python "$P/analyses/c3_piglets_ko_table.py"
echo C2_DONE >&2
