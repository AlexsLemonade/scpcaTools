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
    libfribidi-dev \
    libgdal-dev \
    libgit2-dev \
    libglpk-dev \
    libharfbuzz-dev \
    liblzma-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libudunits2-dev \
    libxml2-dev \
    libxtst6 \
    libxt-dev \
    pandoc \
    python3-pip

rm -rf /var/lib/apt/lists/*
