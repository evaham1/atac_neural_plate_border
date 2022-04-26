FROM  evahamrud/sc_multi_omic-tidyverse_4.1.2:latest

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image with R4.1.2, MACS2 and ArchR"

# Installing MACS2 peak calling:
RUN   conda install -c anaconda pip && \
      pip install MACS2

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

## Other packages:
# Install other packages that are needed
RUN   R -e "install.packages('gridExtra')" && \
      R -e "install.packages('getopt')" && \
      R -e "install.packages('optparse')" && \
      R -e "install.packages('pheatmap')" && \
      R -e "install.packages('patchwork')" && \
      R -e "install.packages('future')" && \
      R -e "install.packages('future.apply')" && \
      R -e "install.packages('clustree')"