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
- `scpcatools_anndata.Dockerfile`: The slim package with the addition of the `zellkonverter` package for reading and writing `AnnData` objects, and the Python `anndata` package.

These images are built automatically on every update to `scpcaTools` (`main`) and pushed to `ghcr.io/alexslemonade/scpcatools:edge-slim` and `ghcr.io/alexslemonade/scpcatools:edge-anndata`, respectively.

## Setting R and Python dependencies

Dependency lock files are built using the `make-requirements.sh` script.
This script depends on `pip-tools` and should be run using Python 3.10 to match the default version of Python in the Docker images.
For convenience, a conda `environment.yml` file is included in this directory that can be used to create an `scpcatools-dev` environment with the necessary Python version and packages.

This script will generate the `requirements.txt` and `renv.lock` files for each set of Python packages to be installed within the various Docker images.
These are based on `requirements.in` files that specify the high-level package requirements for each image.
