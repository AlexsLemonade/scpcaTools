source("renv/activate.R")

# Use CRAN only for Mac(Darwin) and Windows for binary installs
if (Sys.info()[["sysname"]] == "Darwin") {
  options(renv.config.repos.override = c(CRAN = "https://cloud.r-project.org"))
} else if (Sys.info()[["sysname"]] == "Windows") {
  options(renv.config.repos.override = c(CRAN = "https://cloud.r-project.org"))
}
