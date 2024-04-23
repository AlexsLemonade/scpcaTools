##########################
# Build the slim image first --------------------------------------------------
FROM bioconductor/r-ver:3.18 AS slim
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.title "scpcatools-slim"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

#### R packages
# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('remotes', 'renv'))"

WORKDIR /usr/local/renv
COPY renv_slim.lock renv.lock
# restore renv and remove cache files
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# bust cache if needed
ADD "https://api.github.com/repos/AlexsLemonade/scpcaTools/commits?per_page=1" latest_commit
# Install scpcaTools package (& test loading)
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools', upgrade = 'never'); \
  require(scpcaTools)"

WORKDIR /home

##########################
# Add Seurat support
FROM slim AS seurat
LABEL org.opencontainers.image.title "scpcatools-anndata"

WORKDIR /usr/local/renv
COPY renv_seurat.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

WORKDIR /home

##########################
# Add Reports support
FROM slim AS reports
LABEL org.opencontainers.image.title "scpcatools-anndata"

WORKDIR /usr/local/renv
COPY renv_reports.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

WORKDIR /home


##########################
# Add zellkonverter/anndata support
FROM slim AS anndata
LABEL org.opencontainers.image.title "scpcatools-anndata"

WORKDIR /usr/local/renv
COPY renv_zellkonverter.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc)"

#### Python packages
COPY requirements_anndata.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /home

##########################
# Add scvi support
FROM anndata AS scvi
LABEL org.opencontainers.image.title "scpcatools-scvi"

WORKDIR /usr/local/renv
COPY requirements_scvi.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /home

##########################
# Full image with all dependencies
FROM scvi AS full
LABEL org.opencontainers.image.title "scpcatools"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

COPY requirements.txt requirements.txt

RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /home
