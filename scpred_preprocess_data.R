#!/usr/bin/env Rscript 


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(caret))
suppressPackageStartupMessages(require(SingleCellExperiment))

### Extract matrix data and/or training labels from SCE object
### assume that the counts have been normalised previously, i.e. the 'normcounts' slot is present     

# argument parsing 
option_list = list(
    make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object in .rds format'
  ),
    make_option(
    c("-t", "--normalised-counts-slot"),
    action = "store",
    default = "normcounts",
    type = 'character',
    help = 'Name of the slot with normalised counts matrix in SCE object. Default: normcounts'
  ),
    make_option(
    c("-m", "--output-matrix-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output matrix object in .rds format'
  ),
    make_option(
    c("-l", "--output-labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the metadata file with cell type labels in text format'
  )
)

opt = wsc_parse_args(option_list, mandatory=c("input_sce_object", "output_matrix_object"))
sce = readRDS(opt$input_sce_object)

# extract matrix and labels
if(opt$normalised_counts_slot %in% names(assays(sce))){
    matrix = as.matrix(assays(sce)[[opt$normalised_counts_slot]])
} else{
    stop("Specified counts slot not found in SCE object")
}

matrix = apply(matrix, 2, function(x) (x/sum(x))*1000000)
saveRDS(matrix, file=opt$output_matrix_object)
if(!is.na(opt$output_labels)){
    labels = as.data.frame(colData(sce))
    write.csv(labels, file=opt$output_labels)
}
