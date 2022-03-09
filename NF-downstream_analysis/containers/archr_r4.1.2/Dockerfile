FROM  rocker/tidyverse:4.1.2

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR"

## Non-R stuff, lifted from Alex's seurat Docker image
# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      libhdf5-dev \
      libboost-filesystem1.67-dev \
      libxtst6 \
      libxt-dev \
      libfontconfig1-dev \
      libfribidi-dev \
      libgsl0-dev \
      libbz2-dev \
      libclang-dev \
      bzip2-devel.x86_64

## For ArchR:
# Install ArchR
RUN   R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"

# when ran in tidyverse docker image got error: ERROR: dependencies ‘motifmatchr’, ‘chromVAR’, ‘Rsamtools’ are not available for package ‘ArchR’
RUN   R -e "BiocManager::install('motifmatchr')"


# Install ArchR dependencies
RUN   R -e "ArchR::installExtraPackages()"
# ArchR bug means need to install hexbin to make feature plots
RUN   R -e "install.packages('hexbin')"

## For Annotation: 
# Install GenomicFeatures
RUN   R -e "BiocManager::install('GenomicFeatures')"
# Install OrgDb 
RUN   R -e "if (!requireNamespace('org.Gg.eg.db', quietly = TRUE)){BiocManager::install('org.Gg.eg.db')}"
# Install UCSC genome (will have to change to ENSEMBL later)
RUN   R -e "if (!requireNamespace('BSgenome.Ggallus.UCSC.galGal6', quietly = TRUE)){BiocManager::install('BSgenome.Ggallus.UCSC.galGal6')}"
