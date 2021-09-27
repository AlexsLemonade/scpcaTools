#!/bin/bash
set -e

NCPUS=${NCPUS:-1}

apt-get update && apt-get install -y --no-install-recommends \
  apt-utils \
  dialog

# Add curl, bzip2 and dev libs
apt-get -y --no-install-recommends install \
    awscli \
    bzip2 \
    curl \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libgdal-dev \
    libgit2-dev \
    libglpk-dev \
    liblzma-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libudunits2-dev \
    libxml2-dev \
    libxtst6 \
    pandoc \
  && rm -rf /var/lib/apt/lists/*

#### R packages
###############

# Set CRAN mirror to allow access to more recent matrixStats
CRAN="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-08-27"
echo "options(repos = c(CRAN = '${CRAN}'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site
# update all to that repo
Rscript -e "update.packages(ask = FALSE)"

# this comes first since bioconductor packages are
# dependent on updated matrixStats
install2.r --error --skipinstalled -n $NCPUS \
    BiocManager \
    DBI \
    devtools \
    here \
    matrixStats \
    optparse \
    rmarkdown \
    rprojroot \
    RSQLite \
    tidyverse


##########################
# Install bioconductor packages
Rscript -e "withCallingHandlers(
  BiocManager::install(c( \
    'AnnotationHub', \
    'Biostrings', \
    'bluster', \
    'BSgenome', \
    'DropletUtils', \
    'eisaR', \
    'ensembldb', \
    'fishpond', \
    'GenomicFeatures', \
    'miQC', \
    'scater', \
    'scran', \
    'scuttle', \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'tximport'), \
    update = FALSE), \
  warning = function(w) stop(w))"

rm -rf /tmp/downloaded_packages
rm -rf /tmp/Rtmp*
