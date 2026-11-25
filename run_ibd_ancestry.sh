#!/bin/bash
#SBATCH --job-name=ibd_ancestry
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=ibd_ancestry_%j.out
# =============================================================================
# IBD-based clustering + supervised ancestry modelling (IBDseq output -> Leiden clusters -> painting profiles -> NNLS)
#
#   cd .../ancient_modern_all/qpWave-Adm && sbatch run_ibd_ancestry.sh
#
# Starts from the files you already made:  *.paper_filtered.final.with_cM.tsv.gz  (LOD >= 3, >= 2 cM, excess-long-IBD regions removed)
#
# Steps (STEPS="samples pairs cluster paint nnls", default = all):
#   samples  fam_meta/sample_meta.tsv -> sample table (ancient / modern, dates, coordinates, modern target groups)
#   pairs    stream the 3.8 GB IBD table once -> total cM per (pair, chromosome)                [scripts/ibd_01_pair_totals.R]
#   cluster  hierarchical Leiden clustering of the ANCIENT individuals, spatiotemporal collapse    [scripts/ibd_02_cluster.R]
#   paint    painting tables: IBD of every individual with every donor group, per chromosome        [scripts/ibd_03_painting.R]
#   nnls     NNLS ancestry proportions + leave-one-chromosome-out jackknife SEs                     [scripts/ibd_04_nnls.R]
#
# Recommended: first  STEPS="samples pairs cluster" sbatch run_ibd_ancestry.sh , look at $OUTDIR/ibd_cluster_summary.tsv,
# put source clusters in config/ibd_source_pool.tsv / config/ibd_models.tsv, then  STEPS="paint nnls" sbatch run_ibd_ancestry.sh
# Every setting can be overridden from the environment (sbatch --export=ALL ...); see README_ibd.md.
# =============================================================================
ARGS=("$@")
set --
set -eo pipefail

DATA_DIR=${DATA_DIR:-/scratch/nme3923/IranGen/modern_anc_unfilt_PLUS54/merged_datasets/ancient_modern_all}
IBD_DIR=${IBD_DIR:-$DATA_DIR}
IBD_FILE=${IBD_FILE:-$IBD_DIR/ancient.modern.iran.filtered.LDpruned_Allentoft.bgz.paper_filtered.final.with_cM.tsv.gz}
CONDA_ENV=${CONDA_ENV-/scratch/nme3923/conda-envs/glimpse2}      # set CONDA_ENV="" to skip conda activation
CONDA_ACTIVATE=${CONDA_ACTIVATE:-/share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate}
REPO_DIR=${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}
[[ -f $REPO_DIR/scripts/ibd_01_pair_totals.R ]] || REPO_DIR=$DATA_DIR/qpWave-Adm
OUTDIR=${OUTDIR:-$DATA_DIR/ibd_ancestry_out}
STEPS=${STEPS:-"samples pairs cluster paint nnls"}
META=${META:-$REPO_DIR/fam_meta/sample_meta.tsv}

# ---- step 1 (pairs) ----
export MIN_CM=${MIN_CM:-0}                   # extra minimum segment length on top of the file's own 2 cM filter
export ID_SUB=${ID_SUB:-}                    # regex stripped from IDs in the IBD file (e.g. '_[12]$' if IBDseq IDs carry a haplotype suffix)
export CHUNK_LINES=${CHUNK_LINES:-4000000}
# COL_ID1 COL_ID2 COL_CHR COL_LEN_CM COL_START_CM COL_END_CM: only needed if the header is not auto-detected
# ---- step 2 (cluster) ----
export R_START=${R_START:-0.5} R_STEP=${R_STEP:-0.5} R_MAX=${R_MAX:-3}
export MIN_COMMUNITY=${MIN_COMMUNITY:-1}     # smallest community a Leiden split may create
export MAX_DATE_SPAN=${MAX_DATE_SPAN:-1000}  # years: a node is collapsed into one cluster if its dated members span <= this ...
export MAX_KM=${MAX_KM:-500}                 # ... and lie within this distance of each other
export OVERRIDES=${OVERRIDES:-$REPO_DIR/config/ibd_cluster_overrides.tsv}
export SEED=${SEED:-1}
# ---- step 3 (paint) ----
export DONOR_MIN_N=${DONOR_MIN_N:-2}         # groups with fewer members are not used as donors
# ---- step 4 (nnls) ----
export MODELS=${MODELS:-$REPO_DIR/config/ibd_models.tsv}
export SOURCE_POOL=${SOURCE_POOL:-$REPO_DIR/config/ibd_source_pool.tsv}
export MAX_SOURCES=${MAX_SOURCES:-3} MIN_WEIGHT=${MIN_WEIGHT:-0.02} TOP_N=${TOP_N:-5}
export SIZE_NORM=${SIZE_NORM:-1} EXCLUDE_SOURCES=${EXCLUDE_SOURCES:-1} METHOD=${METHOD:-lsei}
PROFILE_NORMS=${PROFILE_NORMS:-"sum none"}   # run the models with both profile normalisations (see README_ibd.md, "Normalisation")
export NCORES=${NCORES:-${SLURM_CPUS_PER_TASK:-4}}

# ------------------------------- environment --------------------------------
if [[ -n $CONDA_ENV ]]; then
  set +eu
  source ~/.bashrc
  source "$CONDA_ACTIVATE"
  conda activate "$CONDA_ENV" || { echo "ERROR: cannot activate conda env $CONDA_ENV"; exit 1; }
  set -e
fi
set -u
command -v Rscript >/dev/null || { echo "ERROR: Rscript not found -> bash setup_ibd_packages.sh"; exit 1; }
Rscript -e 'for (p in c("data.table","igraph")) if (!requireNamespace(p, quietly=TRUE)) quit(status=1)' \
  || { echo "ERROR: R packages data.table / igraph missing -> run on a login node: bash setup_ibd_packages.sh"; exit 1; }
mkdir -p "$OUTDIR"
export OUTDIR SAMPLES=$OUTDIR/ibd_samples.tsv PAIR_RDS=$OUTDIR/ibd_pair_chr.rds GROUPS_TSV=$OUTDIR/ibd_groups.tsv PAINT_RDS=$OUTDIR/painting.rds IBD_FILE
echo "[$(date +%T)] steps: $STEPS   out: $OUTDIR"
{ echo "date: $(date -Is)"; echo "git_commit: $(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || echo NA)"; echo "ibd_file: $IBD_FILE"
  echo "MIN_CM=$MIN_CM R=$R_START/$R_STEP/$R_MAX MIN_COMMUNITY=$MIN_COMMUNITY MAX_DATE_SPAN=$MAX_DATE_SPAN MAX_KM=$MAX_KM DONOR_MIN_N=$DONOR_MIN_N"
  echo "MAX_SOURCES=$MAX_SOURCES MIN_WEIGHT=$MIN_WEIGHT SIZE_NORM=$SIZE_NORM EXCLUDE_SOURCES=$EXCLUDE_SOURCES METHOD=$METHOD PROFILE_NORMS=$PROFILE_NORMS"
} >> "$OUTDIR/run_info.txt"
has() { [[ " $STEPS " == *" $1 "* ]]; }
SC=$REPO_DIR/scripts

if has samples; then
  echo "[$(date +%T)] samples"
  python3 -I "$SC/ibd_make_samples.py" --meta "$META" --out "$SAMPLES"
fi
if has pairs; then
  [[ -s $IBD_FILE ]] || { echo "ERROR: IBD file not found: $IBD_FILE (set IBD_FILE or IBD_DIR)"; exit 1; }
  echo "[$(date +%T)] pairs  (streams $IBD_FILE)"
  Rscript "$SC/ibd_01_pair_totals.R"
fi
if has cluster; then
  echo "[$(date +%T)] cluster"
  Rscript "$SC/ibd_02_cluster.R"
fi
if has paint; then
  echo "[$(date +%T)] paint"
  Rscript "$SC/ibd_03_painting.R"
fi
if has nnls; then
  nmod=$( { grep -v '^\s*#' "$MODELS" 2>/dev/null || true; } | grep -c . || true )
  npool=$( { grep -v '^\s*#' "$SOURCE_POOL" 2>/dev/null || true; } | grep -c . || true )
  if (( nmod > 1 )); then :; else export MODELS=""; fi
  export ENUMERATE=${ENUMERATE:-$(( npool > 0 ? 1 : 0 ))}
  if [[ -z $MODELS && $ENUMERATE != 1 ]]; then
    echo "NOTE: no models in $MODELS and no source pool in $SOURCE_POOL -> nothing to fit. Look at $OUTDIR/ibd_cluster_summary.tsv, fill in config/ibd_models.tsv or config/ibd_source_pool.tsv, then STEPS=\"nnls\"."
  else
    for pn in $PROFILE_NORMS; do
      echo "[$(date +%T)] nnls (PROFILE_NORM=$pn)"
      PROFILE_NORM=$pn OUTDIR="$OUTDIR/nnls_$pn" Rscript "$SC/ibd_04_nnls.R"
    done
  fi
fi
echo "[$(date +%T)] done. Results: $OUTDIR"
ls -1 "$OUTDIR"/*.tsv "$OUTDIR"/nnls_*/*.tsv 2>/dev/null || true
