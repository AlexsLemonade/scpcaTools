
# scpcaTools: Useful tools for analysis of single-cell RNA seq counts data

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/AlexsLemonade/scpcaTools/workflows/R-CMD-check/badge.svg)](https://github.com/AlexsLemonade/scpcaTools/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

The `scpcaTools` package contains a set of tools for working with single-cell and single-nuclei RNA-seq counts data.
Mainly, this package can work with data that has been produced using Alevin, Alevin-Fry, Cell Ranger, or Kallisto.
Counts matrices from either of these pre-processing tools can be imported into R using the `import_quant_data` function to return a `SingleCellExperiment` object to be used for downstream analysis.

Currently, by default the `import_quant_data` function will return a `SingleCellExperiment` object with two assays, the `counts` assay will contain the total counts from both unspliced and spliced reads and the `spliced` assay will contain only reads from the spliced reads.
Input data must be aligned to a reference containing both spliced and unspliced transcripts, otherwise the `include_unspliced=FALSE` option can be used to return just the spliced reads only as the main `counts` assay.
Note that this option is only for reading in data from Alevin, Alevin-fry, and Kallisto and is not used for Cell Ranger.

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
To keep this up to date, `renv::snapshot()` should be run periodically during development (before submitting PRs), which should add any packages that are used in `scpcaTools` scripts and notebooks to `renv.lock`.
When checking out a branch, `renv::restore()` can be used to keep local package installations in sync.

Packages that are not required for the `scpcaTools` scripts directly, but that should be installed in the Docker image can be added to `docker/dependencies.R` to be sure they are part of the image.

Conversely, if there are packages that should _not_ be added to the docker image (such as `scpcaData`), they should be set as ignored in `renv/settings.dcf`.

Note that `renv` is not used when testing the package in CI, so package versions may be different from those specified in `renv.lock`.

## Python package installation and tracking with `renv`

Python packages are also tracked with `renv`.
These are present mostly for use in the Docker image, as the `scpcaTools` package itself does not require any Python packages beyond those installed via [`basilisk`-based packages](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html).
It is important to note that some of the Python functionality may not be exposed in the RStudio GUI, so you will need to use functions in the console and/or terminal to manage these packages.

### Setup

To get started, you will want to open an R session in the project directory and run `renv::restore()`.
You will prompted first to update any R packages (if required), then for the Python enviroment updates.
This will create a `.venv` directory and populate it with the required Python packages, as specified in `requirements.txt`.

To use the python environment from the terminal, you may need to activate it.
This activation is done automatically by RStudio, but for other terminals you will need to run `source .venv/bin/activate`.
Once the environment is activated, you can run `python` to start a Python session with the required packages available.

### Adding Python packages

To add new Python packages, you can use `pip install <package>` from the terminal.
This will _not_ add the packages to `requirements.txt`; to record packages to `requirements.txt`, you should again run `renv::snapshot()` from the R console, which will update both `renv.lock` and `requirements.txt` as required.

Note that `renv::status()` may not report that the `requirements.txt` file is out of date, even if you have installed new packages.

### Docker installation

Any packages in `requirements.txt` will automatically be added to the Docker image when it is built, so you should not need to do anything special beyond the above steps to add them to the docker image.





