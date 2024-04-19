#!/usr/bin/env Rscript

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
setwd(current_dir)

# get dependencies of scpcaTools
renv::snapshot(lockfile = opts$lockfile, type = "explicit")

# additional dependencies for added packages
added_packages <- opts$packages |>
  stringr::str_split_1("[,;\\s]+") |>
  c("tidyverse") |> # always include tidyverse
  tools::package_dependencies() |>
  unlist() |>
  unique()

renv::record(added_packages, lockfile = opts$lockfile)
