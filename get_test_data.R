#!/usr/bin/env Rscript 

### Extract example data to run the tests
library(scPred)
library(Seurat)
library(tidyr)

# extract reference and query datasets
reference = scPred::pbmc_1
query = scPred::pbmc_2

# add processing steps done by Seurat 
reference <- reference %>%
  NormalizeData(normalization.method="RC", scale.factor = 1e6) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# write data 
saveRDS(reference, file = "post_install_tests/reference_pbmc.rds")
saveRDS(query, file = "post_install_tests/query_pbmc.rds")

