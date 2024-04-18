#!/bin/bash
set -euo pipefail
# This script is used to generate requirements files for R and python packages
# Requires that the `pip-tools` python package is installed.

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Python package lists
pip-compile --no-annotate --output-file=requirements.txt requirements.in
pip-compile --no-annotate --output-file=requirements_slim.txt requirements_slim.in

# move up a directory for R scripts to capture the package files
cd ..
Rscript - <<EOF
lockfile <- "docker/renv_slim.lock"
# get dependencies of scpcaTools
renv::snapshot(lockfile = lockfile, type = "explicit")

# additional dependencies for added packages
added_packages <- c(
  "tidyverse"
)
tools::package_dependencies(added_packages) |>
  unlist() |>
  unique() |>
  renv::record(lockfile = lockfile)
EOF
