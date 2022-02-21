FROM  evahamrud/sc_multi_omic-schelper:latest

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Signac and scHelper and ArchR"

# Install ArchR
RUN   R -e "devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())"
# Install GenomicFeatures - required for annotation
RUN   R -e "BiocManager::install("GenomicFeatures")"