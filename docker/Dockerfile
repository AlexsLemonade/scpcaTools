##########################
# Build the slim image first --------------------------------------------------
FROM bioconductor/r-ver:3.18 AS slim
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.title "scpcatools-slim"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

#### R packages
# Use renv for R packages
RUN Rscript -e "install.packages(c('remotes', 'renv')); renv::init(bare=TRUE)"


COPY renv_slim.lock renv.lock
# restore renv and remove cache files
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv/source && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# bust cache if needed
ADD "https://api.github.com/repos/AlexsLemonade/scpcaTools/commits?per_page=1" latest_commit
# Install scpcaTools package (& test loading)
RUN Rscript -e "renv::install('AlexsLemonade/scpcaTools'); \
  require(scpcaTools)"


##########################
# Add Seurat support
FROM slim AS seurat
LABEL org.opencontainers.image.title "scpcatools-anndata"

COPY renv_seurat.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv/source && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*



##########################
# Add reports support
FROM slim AS reports
LABEL org.opencontainers.image.title "scpcatools-anndata"

COPY renv_reports.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv/source && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*



##########################
# Add zellkonverter/anndata support
FROM slim AS anndata
LABEL org.opencontainers.image.title "scpcatools-anndata"

COPY renv_zellkonverter.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv/source && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc)"

#### Python packages
COPY requirements_anndata.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


##########################
# Add scvi support
FROM anndata AS scvi
LABEL org.opencontainers.image.title "scpcatools-scvi"

COPY requirements_scvi.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


##########################
# Full image with all dependencies
FROM scvi AS full
LABEL org.opencontainers.image.title "scpcatools"

COPY renv.lock renv.lock
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv/source && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

COPY requirements.txt requirements.txt

RUN pip install --no-cache-dir -r requirements.txt
