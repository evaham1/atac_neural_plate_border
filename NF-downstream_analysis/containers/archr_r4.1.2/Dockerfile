FROM  rocker/tidyverse:4.1.2

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2 and ArchR"

## For ArchR:
# Install ArchR
RUN   R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
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
