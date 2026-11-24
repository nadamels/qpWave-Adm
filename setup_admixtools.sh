#!/bin/bash
# One-time install of ADMIXTOOLS 2 (R package) into the conda env. Run on a LOGIN node (needs internet).
#   bash setup_admixtools.sh            # installs into /scratch/nme3923/conda-envs/glimpse2; override with CONDA_ENV=...
ARGS=("$@"); set --                      # conda's activate rejects positional arguments
CONDA_ENV=${CONDA_ENV-/scratch/nme3923/conda-envs/glimpse2}
CONDA_ACTIVATE=${CONDA_ACTIVATE:-/share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate}
if [[ -n $CONDA_ENV ]]; then
  source ~/.bashrc
  source "$CONDA_ACTIVATE"
  conda activate "$CONDA_ENV"
fi
set -eo pipefail
INSTALLER=$(command -v mamba || command -v conda)

# R + a C++17 toolchain + binary builds of the heavy dependencies (everything else comes from CRAN below)
$INSTALLER install -y -c conda-forge r-base r-remotes r-rcpp r-rcpparmadillo r-igraph r-dplyr r-tidyr r-readr \
  r-purrr r-ggplot2 r-tibble r-scales r-stringr r-stringi gxx_linux-64 make

Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
remotes::install_github("uqrmaie1/admixtools", dependencies = TRUE, upgrade = "never")
library(admixtools); cat("ADMIXTOOLS 2 installed, version", as.character(packageVersion("admixtools")), "\n")'
