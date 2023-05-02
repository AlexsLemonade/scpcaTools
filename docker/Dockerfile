FROM rocker/r-ver:4.2.3
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

COPY scripts/install_scpca_deps.sh .

RUN bash ./install_scpca_deps.sh

# Use renv for R packages
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages(c('remotes', 'renv'))"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN Rscript -e "renv::consent(provided = TRUE); \
      renv::restore()" && \
      rm -rf ~/.local/share/renv && \
      rm -rf /tmp/downloaded_packages && \
      rm -rf /tmp/Rtmp*

# Complete installation of zellkonverter conda env
RUN Rscript -e "proc <- basilisk::basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata'); \
                basilisk::basiliskStop(proc)"

##########################
# bust cache if needed
ADD "https://api.github.com/repos/AlexsLemonade/scpcaTools/commits?per_page=1" latest_commit
# Install scpcaTools package (& test loading)
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools'); \
  require(scpcaTools)"

# set final workdir for commands
WORKDIR /home
