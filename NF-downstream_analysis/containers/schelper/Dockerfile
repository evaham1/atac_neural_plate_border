FROM  evahamrud/sc_multi_omic-r_signac:v2

LABEL authors="eva.hamrud@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Signac and scHelper"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE)"