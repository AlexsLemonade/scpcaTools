
# scpcaTools: Useful tools for analysis of single-cell RNA seq counts data

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/AlexsLemonade/scpcaTools/workflows/R-CMD-check/badge.svg)](https://github.com/AlexsLemonade/scpcaTools/actions)
  <!-- badges: end --> 

The `scpcaTools` package contains a set of tools for working with single-cell and single-nuclei RNA-seq counts data.
Mainly, this package can work with data that has been produced using Alevin, Alevin-Fry, Cellranger, or Kallisto. 
Counts matrices from either of these pre-processing tools can be imported into R using the `import_quant_data` function to return a `SingleCellExperiment` object to be used for downstream analysis. 

Currently, the `import_quant_data` function can support single-nuclei RNA seq counts data that has been aligned to a pre-mRNA index by using the `intron_mode=TRUE` and `which_counts=intron` options. There is also support to import data processed using Alevin-fry's USA mode by using the `usa_mode=TRUE` option. 

## Installation 

The package can be installed from github with:

```r
remotes::install_github("AlexsLemonade/scpcaTools")
```

## Testing

To test the import functions in this package, you must also have installed the [`scpcaData` package](https://github.com/AlexsLemonade/scpcaData).
This can be done by running `remotes::install_github("AlexsLemonade/scpcaData")`

## Docker image and `renv`

This repository also includes a Dockerfile and associated scripts for building a Docker image that includes the `scpcaTools` package and some additional packages used in other ScPCA repositories, notably [scpca-nf](https://github.com/AlexsLemonade/scpca-nf).
To support this, we use [`renv`](https://rstudio.github.io/renv/index.html) to track the versioned set of packages that will be installed in that Docker image.

Packages that are required for the main `scpcaTools` package should be included in the `renv.lock` file as they are installed and used.
To keep this up to date, `renv::snapshot()` should be run periodically during development(before submitting PRs), which should add any packages that are used in `scpcaTools` scripts and notebooks to `renv.lock`.  
When checking out a branch, `renv::restore()` can be used to keep local package installations in sync.

Packages that are not required for the `scpcaTools` scripts directly, but that should be installed in the Docker image can be added to `docker/dependencies.R` to be sure they are part of the image.

Conversely, if there are packages that should _not_ be added to the docker image (such as `scpcaData`), they should be set as ignored in `renv/settings.dcf`. 

Note that `renv` is not used when testing the package in CI, so package versions may be different from those specified in `renv.lock`.



