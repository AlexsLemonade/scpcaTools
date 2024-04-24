FROM bioconductor/r-ver:3.18
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools
LABEL org.opencontainers.image.title "scpcatools-anndata"

#### R packages
# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('remotes', 'renv'))"

WORKDIR /usr/local/renv
COPY renv_zellkonverter.lock renv.lock
# restore renv and remove cache files
RUN Rscript -e "renv::restore()" && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
  basilisk::basiliskStop(proc); \
  basilisk.utils::cleanConda()"

#### Python packages
COPY requirements_scvi.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

##########################
# bust cache if needed
ADD "https://api.github.com/repos/AlexsLemonade/scpcaTools/commits?per_page=1" latest_commit
# Install scpcaTools package (& test loading)
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools', upgrade = 'never'); \
  require(scpcaTools)"

# set final workdir for commands
WORKDIR /home
