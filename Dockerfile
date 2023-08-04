# Base: https://hub.docker.com/r/rocker/tidyverse/tags (Ubuntu 22.04 with R packages)
# Purpose: Image for CI Pipeline to test this repo
# Build: docker build -t gitlab.spang-lab.de:4687/grw28475/metabodecon:1.1.0 .
# Run: docker run -it --rm -v "$(pwd):/root/metabodecon" -w /root/metabodecon gitlab.spang-lab.de:4687/grw28475/metabodecon:1.1.0
# Push:
#   docker login gitlab.spang-lab.de:4687
#   docker push gitlab.spang-lab.de:4687/grw28475/metabodecon:1.1.0
FROM rocker/tidyverse:4.3

# APT Dependencies
# libgsl27: dependency of R package speaq
# qpdf: required by R CMD check
RUN apt-get update && \
    apt-get install -y --no-install-recommends libgsl27 qpdf &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# R Dependencies
COPY . /root/metabodecon
WORKDIR /root/metabodecon
RUN R -e "devtools::install_deps()"
WORKDIR /root
RUN rm -rf /root/metabodecon

# Default start command
CMD /bin/bash
