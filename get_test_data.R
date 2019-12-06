#!/usr/bin/env Rscript 

### Extract example data to run the tests

download.file("https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/pollen.rds", 
               destfile = "post_install_tests/pollen.rds")

if(!file.exists("post_install_tests/pollen.rds")) stop("Test input file does not exist.")
suppressPackageStartupMessages(require("SingleCellExperiment"))
pollen = readRDS("post_install_tests/pollen.rds")
saveRDS(pollen, file = "post_install_tests/pollen_cpm.rds")
