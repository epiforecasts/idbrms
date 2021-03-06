FROM  rocker/geospatial:latest


RUN apt-get update -y && \
    apt-get install -y \
    libudunits2-dev \
    libgdal-dev \
    libqpdf-dev \
    libmagick++-dev \
    xdg-utils \
    && apt-get clean
    
## Copy files to working directory of server
ADD . /home/rstudio/idbrms

## Set working directory to be this folder
WORKDIR /home/rstudio/idbrms

## Install missing packages
RUN Rscript -e "devtools::install_dev_deps()"

## Install the local version of idbrms
Run R CMD INSTALL --no-multiarch --with-keep.source .
