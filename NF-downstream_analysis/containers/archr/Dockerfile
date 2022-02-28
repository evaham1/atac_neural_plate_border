FROM  evahamrud/sc_multi_omic-schelper:latest

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Signac and scHelper and ArchR"

# Install ArchR
RUN   R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
# Install GenomicFeatures - required for annotation
RUN   R -e "BiocManager::install('GenomicFeatures')"
# Install ArchR dependencies
RUN   R -e "ArchR::installExtraPackages()"
# Install OrgDb
RUN   R -e "if (!requireNamespace('org.Gg.eg.db', quietly = TRUE)){BiocManager::install('org.Gg.eg.db')}"
# Install UCSC genome (will have to change to ENSEMBL later)
RUN   R -e "if (!requireNamespace('BSgenome.Ggallus.UCSC.galGal6', quietly = TRUE)){BiocManager::install('BSgenome.Ggallus.UCSC.galGal6')}"
