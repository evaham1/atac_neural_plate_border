FROM  evahamrud/sc_multi_omic-tidyverse_4.1.2:latest

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with hichipper"

# Ubuntu dependencies
RUN   apt-get install bedtools libssl-dev libcurl4-openssl-dev libxml2-dev

# R packages
RUN   R -e "install.packages('data.table')" && \
      R -e "install.packages('devtools')" && \
      R -e "install.packages('foreach')" && \
      R -e "install.packages('ggplot2')" && \
      R -e "install.packages('knitr')" && \
      R -e "install.packages('networkD3')" && \
      R -e "install.packages('readr')" && \
      R -e "install.packages('reshape2')" && \
      R -e "install.packages('BiocManager')" && \
      R -e "BiocManager::install('diffloop')"

# Install hichipper
RUN   conda install -c anaconda pip \
      pip install hichipper

