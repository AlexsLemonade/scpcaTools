#!/usr/bin/env Rscript

# This script is used to create or update a lockfile based on the current renv environment,
# including scpcaTools dependencies and any additional packages specified by the user.
# This means that it does _not_ look at what files are used by scripts in the project, except
# for the package DESCRIPTION file.
# The main `renv.lock` file is not updated by this script, but must be up to date with the
# current environment before running this script.

library(optparse)

option_list <- list(
  make_option(
    c("-f", "--lockfile"),
    type = "character",
    default = NULL,
    help = "The lockfile to create."
  ),
  make_option(
    c("-p", "--packages"),
    type = "character",
    default = "",
    help = paste(
      "Additional packages to include in the lockfile beyond scpcaTools requirements.",
      "Comma, semicolon, or space separated."
    )
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# change to project root to get renv (renv::load does not set bioc mirrors)
current_dir <- getwd()
setwd(here::here())
source(".Rprofile")

# Check that lockfile is up to date with currently installed packages
if (!renv::status()$synchronized) {
  stop(paste(
    "renv.lock file does not match installed packages.",
    "Please run `renv::restore()` or `renv::snapshot()` from within R before running this script."
  ))
}

# get dependencies of scpcaTools
minimal_packages <- renv::dependencies("DESCRIPTION")$Package

# additional dependencies for added packages
added_packages <- opts$packages |>
  stringr::str_split_1("[,;\\s]+") |>
  stringr::str_subset(".+") |> # remove empty strings
  c("optparse", "scater", "scran", "tidyverse") # always include these

all_packages <- c(minimal_packages, added_packages)

# go back to working directory and snapshot
setwd(current_dir)
renv::snapshot(lockfile = opts$lockfile, packages = all_packages)
