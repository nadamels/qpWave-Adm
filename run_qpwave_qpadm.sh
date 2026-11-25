#!/bin/bash
#SBATCH --job-name=qpwave_qpadm
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=qpwave_qpadm_%j.out
# =============================================================================
# qpWave -> source-set screening -> qpAdm -> plots, for modern Iranian sub-populations
#
#   cd .../ancient_modern_all/qpWave-Adm && sbatch run_qpwave_qpadm.sh
#
# Every setting below can be overridden from the environment, e.g.
#   MAX_SOURCES=4 PTHRESH=0.01 sbatch --export=ALL run_qpwave_qpadm.sh
# Edit config/pop_pools.tsv to change which populations are sources / fixed outgroups.
# One-time install of ADMIXTOOLS 2 (needs internet -> login node):  bash setup_admixtools.sh
# =============================================================================
# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --
set -eo pipefail

# ------------------------------- settings -----------------------------------
DATA_DIR=${DATA_DIR:-/scratch/nme3923/IranGen/modern_anc_unfilt_PLUS54/merged_datasets/ancient_modern_all}
BED_PREFIX=${BED_PREFIX:-all_anc_mod.run3}   # all SNPs (preferred for f2 blocks). LD-pruned alternative: all_anc_mod.run3.LDpruned
CONDA_ENV=${CONDA_ENV-/scratch/nme3923/conda-envs/glimpse2}   # set CONDA_ENV="" to skip conda activation
CONDA_ACTIVATE=${CONDA_ACTIVATE:-/share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate}
REPO_DIR=${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}
[[ -f $REPO_DIR/scripts/qpwave_qpadm.R ]] || REPO_DIR=$DATA_DIR/qpWave-Adm
# DESIGN=default      : rotating qpAdm, p >= 0.05 (config/pop_pools.tsv)
# DESIGN=distal       : ancient pre-Copper-Age components (Anatolia_N, Iran_N, WHG, EEHG, WSHG, Onge, Han, PPN, Natufian) as sources, rotated between left and
#                       right (Narasimhan 2019 SM S4.2.4); p >= 0.01, every weight +-2 SE inside [0,1], up to 4 sources (config/pop_pools_distal.tsv)
# DESIGN=narasimhan   : paper-style (Narasimhan et al. 2019, SM S4.3): small FIXED right set of pre-Copper-Age populations, sources never rotated
#                       onto the right, p >= 0.01, every weight +-2 SE inside [0,1], <= 3 sources, model competition (config/pop_pools_narasimhan.tsv)
DESIGN=${DESIGN:-default}
case $DESIGN in
  default)    POOLS_FILE=${POOLS_FILE:-$REPO_DIR/config/pop_pools.tsv}
              : "${ROTATE:=1}" "${CI_Z:=0}" "${COMPETE:=0}"
              OUTDIR=${OUTDIR:-$DATA_DIR/qpwave_qpadm_out} ;;
  narasimhan) POOLS_FILE=${POOLS_FILE:-$REPO_DIR/config/pop_pools_narasimhan.tsv}
              : "${ROTATE:=0}" "${CI_Z:=2}" "${COMPETE:=1}" "${PTHRESH:=0.01}" "${MAX_SOURCES:=3}" "${MIN_WEIGHT:=0}"
              OUTDIR=${OUTDIR:-$DATA_DIR/qpwave_qpadm_out_narasimhan} ;;
  distal)     POOLS_FILE=${POOLS_FILE:-$REPO_DIR/config/pop_pools_distal.tsv}
              : "${ROTATE:=1}" "${CI_Z:=2}" "${COMPETE:=0}" "${PTHRESH:=0.01}" "${MAX_SOURCES:=4}" "${MIN_WEIGHT:=0}" "${RUN_HOMOGENEITY:=0}"
              OUTDIR=${OUTDIR:-$DATA_DIR/qpwave_qpadm_out_distal} ;;
  *) echo "ERROR: DESIGN must be 'default', 'narasimhan' or 'distal'"; exit 1 ;;
esac
export ROTATE CI_Z COMPETE
# Step 0 inputs (ADMIXTURE K=8 + PCA of the LD-pruned merged data; row order = the fam given here)
RUN_EXPLORE=${RUN_EXPLORE:-1}
EXPLORE_FAM=${EXPLORE_FAM:-$REPO_DIR/fam_meta/all_anc_mod.run3.LDpruned.pca.fam}
EXPLORE_Q=${EXPLORE_Q:-$REPO_DIR/fam_meta/all_anc_mod.run3.LDpruned.8.Q}
EXPLORE_EIGENVEC=${EXPLORE_EIGENVEC:-$REPO_DIR/fam_meta/all_anc_mod.run3.LDpruned.pca.eigenvec}

export PTHRESH=${PTHRESH:-0.05}              # qpWave / qpAdm significance level
export MAX_SOURCES=${MAX_SOURCES:-3}         # largest source set tried
export MIN_SOURCES=${MIN_SOURCES:-2}         # smallest model that may be selected (1-source fits are listed separately)
export MIN_WEIGHT=${MIN_WEIGHT:-0.05}        # each source must contribute >= this fraction, else it is not a real n-way model
export TOP_N=${TOP_N:-3}                     # competing models kept per target (SEs + plots)
export WEIGHT_TOL=${WEIGHT_TOL:-0}           # allow weights in [-tol, 1+tol]
export FORCE_F2=${FORCE_F2:-0}               # 1 = recompute f2 blocks
export NCORES=${NCORES:-${SLURM_CPUS_PER_TASK:-8}}
export MAXMEM_MB=${MAXMEM_MB:-$(( ${SLURM_MEM_PER_NODE:-64000} * 6 / 10 ))}
export PLINK_PREFIX=$OUTDIR/plink/data OUTDIR GROUP_SUMMARY=$OUTDIR/group_summary.tsv
export TMPDIR=${TMPDIR:-$OUTDIR/tmp}; mkdir -p "$TMPDIR"   # temp files on scratch, not in home
export OMP_NUM_THREADS=1                     # parallelism is across models, not inside one fit

# ------------------------------- environment --------------------------------
if [[ -n $CONDA_ENV ]]; then
  set +eu                                      # .bashrc / conda activate are not strict-mode safe
  source ~/.bashrc
  source "$CONDA_ACTIVATE"
  conda activate "$CONDA_ENV" || { echo "ERROR: cannot activate conda env $CONDA_ENV"; exit 1; }
  set -e
fi
set -u
command -v Rscript >/dev/null  || { echo "ERROR: Rscript not found in env '${CONDA_ENV:-<current>}' -> run: bash setup_admixtools.sh"; exit 1; }
command -v python3 >/dev/null  || { echo "ERROR: python3 not found"; exit 1; }
Rscript -e 'suppressMessages(library(admixtools))' 2>/dev/null \
  || { echo "ERROR: R package 'admixtools' (ADMIXTOOLS 2) is not installed -> run on a login node: bash setup_admixtools.sh"; exit 1; }

for f in "$DATA_DIR/$BED_PREFIX".{bed,bim,fam} "$REPO_DIR/fam_meta/sample_meta.tsv" "$POOLS_FILE"; do
  [[ -s $f ]] || { echo "ERROR: missing input $f"; exit 1; }
done
mkdir -p "$OUTDIR/plink"
echo "[$(date +%T)] data: $DATA_DIR/$BED_PREFIX   out: $OUTDIR   cores: $NCORES   maxmem: ${MAXMEM_MB}MB"

# ---------------- provenance: the candidate pool is frozen BEFORE qpAdm is looked at ----------
cp "$POOLS_FILE" "$OUTDIR/pop_pools.used.tsv"
{
  echo "date: $(date -Is)"; echo "git_commit: $(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || echo NA)"
  echo "git_dirty: $(git -C "$REPO_DIR" status --porcelain config scripts 2>/dev/null | wc -l) modified file(s) in config/ scripts/"
  echo "pop_pools_md5: $(md5sum "$OUTDIR/pop_pools.used.tsv" | cut -d' ' -f1)"
  echo "data: $DATA_DIR/$BED_PREFIX"
  echo "DESIGN=$DESIGN POOLS_FILE=$POOLS_FILE ROTATE=$ROTATE CI_Z=$CI_Z COMPETE=$COMPETE"
  echo "PTHRESH=$PTHRESH MIN_SOURCES=$MIN_SOURCES MAX_SOURCES=$MAX_SOURCES MIN_WEIGHT=$MIN_WEIGHT WEIGHT_TOL=$WEIGHT_TOL TOP_N=$TOP_N"
} > "$OUTDIR/run_info.txt"

# ---------------- 0. ADMIXTURE / PCA exploration (hypothesis generation only; never an input to qpAdm) ----
if [[ $RUN_EXPLORE == 1 && -s $EXPLORE_Q && -s $EXPLORE_EIGENVEC && -s $EXPLORE_FAM ]]; then
  echo "[$(date +%T)] step 0/3  ADMIXTURE + PCA exploration (candidate coverage check)"
  python3 -I "$REPO_DIR/scripts/explore_admixture_pca.py" --meta "$REPO_DIR/fam_meta/sample_meta.tsv" --fam "$EXPLORE_FAM" \
      --pools "$POOLS_FILE" --q "$EXPLORE_Q" --eigenvec "$EXPLORE_EIGENVEC" --outdir "$OUTDIR/explore" \
    && OUTDIR="$OUTDIR/explore" Rscript "$REPO_DIR/scripts/plot_admixture_pca.R" \
    || echo "WARNING: exploration step failed; continuing with qpWave/qpAdm (it does not depend on it)"
else
  echo "[$(date +%T)] step 0/3  skipped (RUN_EXPLORE=$RUN_EXPLORE or inputs missing)"
fi

# ---------------- 1. population labels: FID = target group / source group -----
echo "[$(date +%T)] step 1/3  building population-labelled .fam"
python3 -I "$REPO_DIR/scripts/make_pops.py" \
  --meta  "$REPO_DIR/fam_meta/sample_meta.tsv" \
  --fam   "$DATA_DIR/$BED_PREFIX.fam" \
  --pools "$POOLS_FILE" \
  --sublabels \
  --outdir "$OUTDIR"
ln -sf "$DATA_DIR/$BED_PREFIX.bed" "$OUTDIR/plink/data.bed"     # original genotypes are never modified
ln -sf "$DATA_DIR/$BED_PREFIX.bim" "$OUTDIR/plink/data.bim"
ln -sf "$DATA_DIR/$BED_PREFIX.bed" "$OUTDIR/plink_sub/data.bed"       # sub-label labelling for the homogeneity QC
ln -sf "$DATA_DIR/$BED_PREFIX.bim" "$OUTDIR/plink_sub/data.bim"

# ---------------- 2. f2 blocks -> qpWave -> source screen -> qpAdm -----------
echo "[$(date +%T)] step 2/3  qpWave + qpAdm"
Rscript "$REPO_DIR/scripts/qpwave_qpadm.R"

# ---------------- 2b. homogeneity QC of the pooled groups (sub-label pairs, qpWave rank 0) ----------------
if [[ ${RUN_HOMOGENEITY:-1} == 1 ]]; then
  echo "[$(date +%T)] step 2b/3 homogeneity QC of pooled groups"
  PLINK_PREFIX_SUB="$OUTDIR/plink_sub/data" Rscript "$REPO_DIR/scripts/homogeneity_test.R" \
    || echo "WARNING: homogeneity QC failed; the qpAdm results above do not depend on it"
fi

# ---------------- 3. plots ---------------------------------------------------
echo "[$(date +%T)] step 3/3  plots"
Rscript "$REPO_DIR/scripts/plot_qpadm.R"

echo "[$(date +%T)] done. Results: $OUTDIR"
ls -1 "$OUTDIR"/*.tsv "$OUTDIR"/*.png 2>/dev/null
