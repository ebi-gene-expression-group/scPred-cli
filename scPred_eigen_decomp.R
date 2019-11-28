#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))

# Calculate n first pricipal components and apply log-transform to the matrix if specified; initialise scPred object 

option_list = list(
    make_option(
        c("-i", "--training-matrix"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the training matrix in .rds format'
  ),
    make_option(
        c("-l", "--log-transform"), 
        action = "store",
        default = TRUE,
        type = 'logical',
        help = 'Should log-transform be performed on the matrix? Default: TRUE'
  ),

    make_option(
        c("-t", "--training-labels"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the training labels (metadata) in text format'
  ),
    make_option(
        c("-n", "--principal-comps"), 
        action = "store",
        default = 10,
        type = 'numeric',
        help = 'Number of pricipal components for eigenvalue decomposition'
  ),
    make_option(
        c("-s", "--random-seed"), 
        action = "store",
        default = 123,
        type = 'numeric',
        help = 'Seed for random number generator'
  ),
    make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for the scPred object'
   )
)

opt = wsc_parse_args(option_list, mandatory = c("training_matrix", 
                                                "training_labels",
                                                "output_path"))

set.seed(opt$random_seed)
labels = read.csv(opt$training_labels, header=TRUE)
row.names(labels) = labels[,1]
labels = labels[,-1]
mat = readRDS(opt$training_matrix)
# create scPred object 
scp = eigenDecompose(mat, 
                     n = opt$principal_comps, 
                     seed = opt$random_seed, 
                     pseudo = opt$log_transform)

scPred::metadata(scp) = data.frame(labels)
saveRDS(scp, opt$output_path)
