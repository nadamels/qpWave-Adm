#!/bin/bash
# One-time install of ADMIXTOOLS 2 (R package) into the conda env. Run on a LOGIN node (needs internet).
#   bash setup_admixtools.sh            # installs into /scratch/nme3923/conda-envs/glimpse2; override with CONDA_ENV=...
#
# The source is fetched with plain `git clone` (github.com), NOT remotes::install_github, which needs api.github.com
# (often blocked on HPC). Fallbacks: codeload tarball, or a copy you made yourself:
#   ADMIXTOOLS_SRC=/path/to/admixtools-master.tar.gz   (or .zip, or an unpacked directory)  bash setup_admixtools.sh
ARGS=("$@"); set --                      # conda's activate rejects positional arguments
# Keep EVERYTHING off the (small) home directory: conda package cache, temp files and R packages go to scratch.
SCRATCH_BASE=${SCRATCH_BASE:-/scratch/nme3923}
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-$SCRATCH_BASE/conda-pkgs}
export TMPDIR=${TMPDIR:-$SCRATCH_BASE/tmp}
mkdir -p "$CONDA_PKGS_DIRS" "$TMPDIR"
export R_LIBS_USER=/nonexistent          # never create/use ~/R; packages are installed into the conda env (on scratch)
CONDA_ENV=${CONDA_ENV-/scratch/nme3923/conda-envs/glimpse2}
CONDA_ACTIVATE=${CONDA_ACTIVATE:-/share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate}
if [[ -n $CONDA_ENV ]]; then
  source ~/.bashrc
  source "$CONDA_ACTIVATE"
  conda activate "$CONDA_ENV"
fi
set -eo pipefail
INSTALLER=$(command -v mamba || command -v conda)

# ---- 1. R + C++17 toolchain + binary builds of the dependencies (no source compilation, no CRAN needed for these)
$INSTALLER install -y -c conda-forge r-base r-remotes r-rcpp r-rcpparmadillo r-igraph r-dplyr r-tidyr r-readr \
  r-purrr r-ggplot2 r-tibble r-scales r-stringr r-stringi gxx_linux-64 make
EXTRA=(r-abind r-cli r-crayon r-cubelyr r-data.tree r-digest r-doparallel r-foreach r-furrr r-future r-jsonlite
       r-lobstr r-magrittr r-parallelly r-pgenlibr r-quadprog r-rlang)
$INSTALLER install -y -c conda-forge "${EXTRA[@]}" || {
  echo "NOTE: not all extra packages are on conda-forge; trying them one by one (the rest will come from CRAN)"
  for p in "${EXTRA[@]}"; do $INSTALLER install -y -c conda-forge "$p" || echo "  - $p not available from conda-forge"; done
}

# ---- 2. get the ADMIXTOOLS 2 source (no api.github.com)
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
SRC_DIR=
if [[ -n ${ADMIXTOOLS_SRC:-} ]]; then
  case $ADMIXTOOLS_SRC in
    *.tar.gz|*.tgz) tar -xzf "$ADMIXTOOLS_SRC" -C "$WORK" ;;
    *.zip)          unzip -q "$ADMIXTOOLS_SRC" -d "$WORK" ;;
    *)              cp -r "$ADMIXTOOLS_SRC" "$WORK/admixtools-src" ;;
  esac
elif git clone --depth 1 https://github.com/uqrmaie1/admixtools.git "$WORK/admixtools-src"; then
  echo "cloned with git"
elif curl -fL --retry 2 -o "$WORK/a.tgz" https://codeload.github.com/uqrmaie1/admixtools/tar.gz/refs/heads/master \
     || wget -q -O "$WORK/a.tgz" https://codeload.github.com/uqrmaie1/admixtools/tar.gz/refs/heads/master; then
  tar -xzf "$WORK/a.tgz" -C "$WORK"
else
  cat <<MSG
ERROR: this node cannot reach github.com. On a machine that can, download
  https://github.com/uqrmaie1/admixtools/archive/refs/heads/master.zip
copy it to the cluster (scp) and run:
  ADMIXTOOLS_SRC=/path/to/master.zip bash setup_admixtools.sh
MSG
  exit 1
fi
SRC_DIR=$(find "$WORK" -maxdepth 2 -name DESCRIPTION -printf '%h\n' | head -1)
[[ -n $SRC_DIR ]] || { echo "ERROR: no R package (DESCRIPTION) found in the downloaded source"; exit 1; }
echo "ADMIXTOOLS 2 source: $SRC_DIR (commit: $(git -C "$SRC_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown))"

# ---- 3. install; anything still missing is fetched from CRAN (needs cran access), otherwise it is reported
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
src <- commandArgs(TRUE)[1]
imp <- trimws(sub("\\(.*", "", strsplit(read.dcf(file.path(src, "DESCRIPTION"), "Imports")[1], ",")[[1]]))
imp <- setdiff(imp[nzchar(imp)], c("stats", "utils", "graphics", "grDevices", "methods"))
miss <- imp[!vapply(imp, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  cat("Installing missing dependencies from CRAN:", paste(miss, collapse = ", "), "\n")
  tryCatch(install.packages(miss, lib = .Library), error = function(e) NULL)
  miss <- imp[!vapply(imp, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Could not install: ", paste(miss, collapse = ", "),
                         "  -> install them with conda (conda install -c conda-forge r-<name>) and rerun")
}
install.packages(src, repos = NULL, type = "source", lib = .Library)   # .Library = R library inside the conda env
library(admixtools); cat("ADMIXTOOLS 2 installed, version", as.character(packageVersion("admixtools")), "\n")' "$SRC_DIR"
