#!/bin/bash
# One-time install of the R packages used by run_ibd_ancestry.sh into the conda env. Run on a LOGIN node (needs internet).
#   bash setup_ibd_packages.sh            # CONDA_ENV=... to override
ARGS=("$@"); set --
CONDA_ENV=${CONDA_ENV-/scratch/nme3923/conda-envs/glimpse2}
CONDA_ACTIVATE=${CONDA_ACTIVATE:-/share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate}
if [[ -n $CONDA_ENV ]]; then
  source ~/.bashrc
  source "$CONDA_ACTIVATE"
  conda activate "$CONDA_ENV"
fi
set -eo pipefail
INSTALLER=$(command -v mamba || command -v conda)
$INSTALLER install -y -c conda-forge r-base r-data.table r-igraph r-limsolve r-quadprog r-ggplot2 r-remotes r-rcpp gxx_linux-64 make
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
if (!requireNamespace("leidenAlg", quietly = TRUE)) {
  ok <- tryCatch({ install.packages("leidenAlg"); requireNamespace("leidenAlg", quietly = TRUE) }, error = function(e) FALSE)
  if (!ok) remotes::install_github("kharchenkolab/leidenAlg", upgrade = "never")
}
for (p in c("data.table", "igraph", "limSolve", "leidenAlg", "ggplot2"))
  cat(sprintf("%-12s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING"))'
