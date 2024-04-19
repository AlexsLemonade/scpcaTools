#!/bin/bash
set -euo pipefail
# This script is used to generate requirements files for R and Python packages
# Requires that the `pip-tools` python package is installed.

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Python package lists
pip-compile --no-annotate --strip-extras --output-file=requirements.txt requirements.in
pip-compile --no-annotate --strip-extras --output-file=requirements_anndata.txt requirements_anndata.in

# slim lockfile
Rscript scripts/make-lockfile.R -f renv_slim.lock

# zellkonverter lockfile
Rscript scripts/make-lockfile.R -f renv_zellkonverter.lock -p zellkonverter

# Seurat lockfile
Rscript scripts/make-lockfile.R -f renv_seurat.lock -p Seurat
