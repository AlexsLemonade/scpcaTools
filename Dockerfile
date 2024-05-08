##########################
# Build the slim image target first --------------------------------------------
FROM bioconductor/r-ver:3.19 AS slim
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.title "scpcatools-slim"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

#### R packages
# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('renv'))"

COPY docker/renv_slim.lock renv.lock
# restore renv and remove cache files
RUN Rscript -e 'renv::restore()'\
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*


# Install scpcaTools package (& test loading)
COPY . scpcaTools
RUN R CMD INSTALL scpcaTools/
RUN Rscript -e "require(scpcaTools)"

# restore slim renv
COPY docker/renv_slim.lock renv.lock

##########################
# Add Seurat support target ----------------------------------------------------
FROM slim AS seurat
LABEL org.opencontainers.image.title "scpcatools-seurat"

COPY docker/renv_seurat.lock renv.lock
RUN Rscript -e 'renv::restore()'\
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*


##########################
# Add Reports support target ---------------------------------------------------
FROM slim AS reports
LABEL org.opencontainers.image.title "scpcatools-reports"

RUN apt-get -y update &&  \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --no-install-recommends -y \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

COPY docker/renv_reports.lock renv.lock
RUN Rscript -e 'renv::restore()'\
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*


##########################
# Add zellkonverter/anndata support target -------------------------------------
FROM slim AS anndata
LABEL org.opencontainers.image.title "scpcatools-anndata"

COPY docker/renv_zellkonverter.lock renv.lock
RUN Rscript -e 'renv::restore()'\
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc); \
  basilisk.utils::cleanConda()"

#### Python packages
COPY docker/requirements_anndata.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


##########################
# Add scvi support target ------------------------------------------------------
FROM anndata AS scvi
LABEL org.opencontainers.image.title "scpcatools-scvi"

COPY docker/requirements_scvi.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


##########################
# Full image with all dependencies ---------------------------------------------
FROM scvi AS full
LABEL org.opencontainers.image.title "scpcatools"

RUN apt-get -y update &&  \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --no-install-recommends -y \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

COPY renv.lock renv.lock
RUN Rscript -e 'renv::restore()'\
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

COPY docker/requirements.txt requirements.txt

RUN pip install --no-cache-dir -r requirements.txt
