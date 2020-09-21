#!/usr/bin/env Rscript 

### Extract example data to run the tests
library(scPred)

# extract reference and query datasets
reference = scPred::pbmc_1
query = scPred::pbmc_2

# write data 
saveRDS(reference, file = "post_install_tests/reference_pbmc.rds")
saveRDS(query, file = "post_install_tests/query_pbmc.rds")

