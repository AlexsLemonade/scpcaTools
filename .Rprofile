# Use PPM for linux (not Mac or Windows) (CRAN by default)
if (! Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) {
  options(renv.config.repos.override = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
}

# activate renv
source("renv/activate.R")
