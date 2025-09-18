#!/bin/bash
set -euo pipefail
# This script is used to generate or update requirements/lock files for R and Python packages
# Requires that the `pip-tools` python package is installed.
# Before running, make sure that the renv.lock file and installed libaries are
# consistent with renv::snapshot() or renv::restore()

# To upgrade packages, set the `UPGRADE_PY` environment variable as follows:
# UPGRADE_PY=1 bash make-requirements.sh

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

UPGRADE_PY=${UPGRADE_PY:-0}

if [ "$UPGRADE_PY" -eq 0 ]; then
  UPGRADE_PY_FLAG=""
else
  UPGRADE_PY_FLAG="--upgrade"
fi

# Python package lists
pip-compile --no-annotate --strip-extras $UPGRADE_PY_FLAG --output-file=requirements.txt requirements.in
pip-compile --no-annotate --strip-extras $UPGRADE_PY_FLAG --output-file=requirements_anndata.txt requirements_anndata.in
pip-compile --no-annotate --strip-extras $UPGRADE_PY_FLAG --output-file=requirements_scvi.txt requirements_scvi.in

# slim lockfile
Rscript scripts/make-lockfile.R -f renv_slim.lock -p batchelor,rtracklayer,scDblFinder

# reports lockfile
Rscript scripts/make-lockfile.R -f renv_reports.lock -p ComplexHeatmap,ggforce,kableExtra,rmarkdown,ggmap,DT,patchwork,ggridges

# zellkonverter lockfile
Rscript scripts/make-lockfile.R -f renv_zellkonverter.lock -p zellkonverter

# Seurat lockfile
Rscript scripts/make-lockfile.R -f renv_seurat.lock -p Seurat

# infercnv lockfile
Rscript scripts/make-lockfile.R -f renv_infercnv.lock -p infercnv
