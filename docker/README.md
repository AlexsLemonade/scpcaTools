This folder contains a Dockerfile to build an image with R, scpcaTools and all of its dependencies.

The image is built from a versioned Rocker image, with additional R packages and other utilities added as needed.

It can be built with the following command run in this directory:

```
docker buildx build . -t scpcatools --platform linux/amd64
```

Note that this image does not include RStudio, in attempt to make the image smaller.

This image is built automatically on every update to `scpcaTools` (`main`) and pushed to `ghcr.io/alexslemonade/scpca-tools:edge`.

Release versions of `scpcaTools` will be tagged as `ghcr.io/alexslemonade/scpcatools:latest` and tagged by version number.

The main image includes all of the recommended packages for `scpcaTools` (aside from those only required for development), as well as the Python `scanpy` and `scvi-tools` packages and dependencies.

There are also smaller images that can be built, which include only the R packages required for `scpcaTools`, with particular additions for specific purposes. These include:

- `scpcatools_slim.Dockerfile`: A smaller image that includes only the required R packages for `scpcaTools`, with no additional packages.
- `scpcatools_scvi.Dockerfile`: The slim package with the addition of the `zellkonverter` package for reading and writing `AnnData` objects, and the Python `anndata` and `scvi` packages.
- `scpcatools_reports.Dockerfile`: The slim package with the addition `rmarkdown` and `ComplexHeatmap` and associated packages for generating reports.
- `scpcatools_seurat.Dockerfile`: The slim package with the addition `Seurat`.

These images are built automatically on every update to `scpcaTools` (`main`) and pushed to `ghcr.io/alexslemonade/scpcatool-slims:edge` and `ghcr.io/alexslemonade/scpcatools-scvi:edge`, etc.

## Generating R and Python dependency files

Dependency lock files, with the exception of `renv.lock`, are built using the `make-requirements.sh` script.
This script depends on `pip-tools` and should be run using Python 3.10 to match the default version of Python in the Docker images.
For convenience, a conda `environment.yml` file is included in this directory that can be used to create an `scpcatools-dev` environment with the necessary Python version and packages.

This script will generate the `requirements_*.txt` and `renv_*.lock` files for each set of Python packages to be installed within the various Docker images.
- The requirements files are based on `requirements_*.in` files that specify the high-level package requirements for each image.
- The `renv_*.lock` files are strict subsets of the main `renv.lock` file that is use for development and for the full image.
  - The main `renv.lock` file should be updated manually using `renv::snapshot()` from within the `scpcaTools` project.
  - After updating `renv.lock`, the `make-requirements.sh` script should be run to update the `renv_*.lock` files.
