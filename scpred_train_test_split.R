#!/usr/bin/env Rscript


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(caret))
suppressPackageStartupMessages(require(SingleCellExperiment))

### Split the dataset into train/test
### training/test labels need to be supplied in ColData slot of SCE object   

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
    c("-m", "--training-matrix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output path for training matrix in .rds format'
  ),

    make_option(
    c("-n", "--test-matrix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output path for test matrix in .rds format'
  ),

    make_option(
    c("-c", "--cell-types-column"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Column name for assigned cell types'
  ),

    make_option(
    c("-o", "--training-labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output path for training labels in text format'
  ),
    make_option(
    c("-p", "--test-labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output path for test labels in text format'
  ),
    make_option(
    c("-s", "--random-seed"),
    action = "store",
    default = 123,
    type = 'numeric',
    help = 'Seed for random number generation'
  ), 
    make_option(
    c("-r", "--training-ratio"),
    action = "store",
    default = 0.70,
    type = 'numeric',
    help = 'Proportion of training/testing dataset split'
  ) 
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object",
                                                "training_matrix",
                                                "test_matrix",
                                                "training_labels",
                                                "cell_types_column",
                                                "test_labels"))

set.seed(opt$random_seed)
print(opt$input_sce_object)
# preprocess data from SCE object 
sce_object = readRDS(opt$input_sce_object)
sce_counts = as.matrix(assays(sce_object)[[opt$normalised_counts_slot]])

sce_counts_cpm = apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
sce_metadata = as.data.frame(colData(sce_object))

# run data partitioning 
i = createDataPartition(sce_metadata[, opt$cell_types_column], p = opt$training_ratio, list = FALSE)
train_data =  sce_counts_cpm[, i]
train_labels = sce_metadata[i, , drop = FALSE]

test_data = sce_counts_cpm[, -i]
test_labels = sce_metadata[-i, , drop = FALSE]

saveRDS(train_data, file = opt$training_matrix)
saveRDS(test_data, file = opt$test_matrix)
write.csv(train_labels, file = opt$training_labels)
write.csv(test_labels, file = opt$test_labels)
