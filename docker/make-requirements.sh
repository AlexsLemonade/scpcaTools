#!/bin/bash
set -euo pipefail
# This script is used to generate requirements files for R and Python packages
# Requires that the `pip-tools` python package is installed.
# Before running, make sure that the renv.lock file and installed libaries are
# consistent with renv::snapshot() or renv::restore()

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Python package lists
pip-compile --no-annotate --strip-extras --output-file=requirements.txt requirements.in
pip-compile --no-annotate --strip-extras --output-file=requirements_anndata.txt requirements_anndata.in
pip-compile --no-annotate --strip-extras --output-file=requirements_scvi.txt requirements_scvi.in

# slim lockfile
Rscript scripts/make-lockfile.R -f renv_slim.lock

# reports lockfile
Rscript scripts/make-lockfile.R -f renv_reports.lock -p ComplexHeatmap,ggforce,kableExtra,rmarkdown

# zellkonverter lockfile
Rscript scripts/make-lockfile.R -f renv_zellkonverter.lock -p zellkonverter

# Seurat lockfile
Rscript scripts/make-lockfile.R -f renv_seurat.lock -p Seurat
