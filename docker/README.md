This folder contains a Dockerfile to build an image with R, scpcaTools and all of its dependencies.

The image is built from a versioned Rocker image, with additional R packages and other utilities added as needed.

It can be built with the following command run in this directory:

```
docker buildx build . -t scpca-tools --platform linux/amd64
```

Note that this image does not include RStudio, in attempt to make the image smaller. 

This image is built automatically on every update to `scpcaTools` (`main`) and pushed to `ghcr.io/alexslemonade/scpca-tools:edge`.

Release versions of `scpcaTools` will be tagged as `ghcr.io/alexslemonade/scpca-tools:latest` and tagged by version number.
