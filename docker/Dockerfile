FROM rocker/r-ver:4.0.5
LABEL maintainer="ccdl@alexslemonade.org"

COPY scripts/install_scpca_deps.sh .

RUN bash ./install_scpca_deps.sh

##########################
# Install scpcaTools package
RUN Rscript -e "remotes::install_github('AlexsLemonade/scpcaTools')"

# set final workdir for commands
WORKDIR /home
