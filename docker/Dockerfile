FROM rocker/r-ver:4.3.2
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

COPY scripts/install_scpca_deps.sh .

RUN bash ./install_scpca_deps.sh

#### R packages
# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('remotes', 'renv'))"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
# restore with PPM repo set for binary installs
RUN Rscript -e "renv::consent(provided = TRUE); \
      renv::restore(repos = 'https://packagemanager.posit.co/cran/latest')" && \
      rm -rf ~/.local/share/renv && \
      rm -rf /tmp/downloaded_packages && \
      rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
ENV BASILISK_EXTERNAL_DIR /usr/local/renv/basilisk
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
                basilisk::basiliskStop(proc)"

#### Python packages
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt && pip cache purge


##########################
# bust cache if needed
ADD "https://api.github.com/repos/AlexsLemonade/scpcaTools/commits?per_page=1" latest_commit
# Install scpcaTools package (& test loading)
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools', upgrade = 'never'); \
  require(scpcaTools)"

# set final workdir for commands
WORKDIR /home
