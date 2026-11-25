#!/bin/bash
# Plumbing + sanity test of the IBD pipeline on simulated data (needs python3 with numpy, and R with data.table + igraph):
#   bash tests/run_ibd_test.sh [workdir]
# Truth: Mod_Zoroastrian = 100% Steppe_BA ; Mod_Central_Iranian = 70% Armenia_BA + 30% Anatolia_N ; Mod_Turk = 50% Anatolia_N + 50% Han_like.
set -eo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
W=${1:-$(mktemp -d)}; mkdir -p "$W"
python3 "$HERE/tests/simulate_ibd.py" "$W/sim"
# cluster names depend on the clustering result, so run the clustering first and read the names back
CONDA_ENV="" REPO_DIR=$HERE META=$W/sim/sample_meta.tsv IBD_FILE=$W/sim/fake.final.with_cM.tsv.gz OUTDIR=$W/out CHUNK_LINES=30000 \
  STEPS="samples pairs cluster" bash "$HERE/run_ibd_ancestry.sh"
pool=$(cut -f1 "$W/out/ibd_cluster_summary.tsv" | tail -n +2)
echo "$pool" > "$W/pool.txt"
nm() { echo "$pool" | grep "^$1"; }
printf 'model_id\ttarget\tsources\nz\tMod_Zoroastrian\t%s\ncentral\tMod_Central_Iranian\t%s,%s\nturk\tMod_Turk\t%s,%s\n' \
  "$(nm Steppe)" "$(nm Armenia)" "$(nm Anatolia)" "$(nm Anatolia)" "$(nm Han)" > "$W/models.tsv"
CONDA_ENV="" REPO_DIR=$HERE OUTDIR=$W/out MODELS=$W/models.tsv SOURCE_POOL=$W/pool.txt STEPS="paint nnls" bash "$HERE/run_ibd_ancestry.sh"
echo; echo "=== explicit models, PROFILE_NORM=none (unbiased) and sum (default, biased towards sources with more total IBD) ==="
for pn in none sum; do echo "-- $pn"; awk -F'\t' '$NF=="explicit" || 1' "$W/out/nnls_$pn/ibd_nnls_models.tsv" | cut -f1,2,5,6,7 | grep -E '^(z|central|turk|model_id)'; done
