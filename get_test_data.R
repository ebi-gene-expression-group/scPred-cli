#!/usr/bin/env Rscript 
### Extract example data for tesing purpose

download.file("https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/pollen.rds", 
               destfile = "post_install_tests/pollen.rds")

if(!file.exists("post_install_tests/pollen.rds")) stop("Test input file does not exist.")
suppressPackageStartupMessages(require("SingleCellExperiment"))
pollen = readRDS("post_install_tests/pollen.rds")
#pollen_counts = normcounts(pollen)
#normcounts(pollen) = apply(normcounts(pollen), 2, function(x) (x/sum(x))*1000000)
#pollen_metadata = as.data.frame(colData(pollen))

saveRDS(pollen, file = "post_install_tests/pollen_cpm.rds")
#write.csv(pollen_metadata, file = "post_install_tests/pollen_metadata.txt", sep = "\t")
