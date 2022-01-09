FROM rocker/r-ver:4.1.1
LABEL maintainer="ccdl@alexslemonade.org"
LABEL org.opencontainers.image.source https://github.com/AlexsLemonade/scpcaTools

COPY scripts/install_scpca_deps.sh .

RUN bash ./install_scpca_deps.sh

# Use renv for R packages
ENV RENV_VERSION 0.15.0
ENV RENV_CONFIG_CACHE_ENABLED FALSE
RUN Rscript -e "install.packages('remotes'); \
      remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /usr/local/renv
COPY renv.lock renv.lock
RUN Rscript -e "renv::consent(provided = TRUE); \
      renv::restore()" && \
      rm -rf ~/.local/share/renv && \
      rm -rf /tmp/downloaded_packages && \
      rm -rf /tmp/Rtmp*

##########################
# Install scpcaTools package (& test loading)
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools'); \
  require(scpcaTools)"

# set final workdir for commands
WORKDIR /home
