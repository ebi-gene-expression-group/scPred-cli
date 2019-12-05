#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))

# Make cell type predicitons using trained model
# this script can be used both for evaluation of model performance on test data and obtaining predictions 
# on new data 

option_list = list(
    make_option(
        c("-i", "--input-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input object of scPred or seurat class in .rds format'
  ),
    make_option(
        c("-p", "--pred-data"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input prediction matrix in .rds format'
  ),
    make_option(
        c("-s", "--test-labels"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the test labels for evalutation of model performance'
  ),
    make_option(
        c("-r", "--cell-types-column"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Column name of true labels in provided metadata file'
  ),

    make_option(
        c("-l", "--threshold-level"), 
        action = "store",
        default = 0.9,
        type = 'numeric',
        help = 'Classification threshold value'
  ),
     make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for values predicted by the model'
  ),
     make_option(
        c("-a", "--plot-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for prediction probabilities histograms'
  ), 
     make_option(
        c("-b", "--confusion-table"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for confusion table in text format'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_object",
                                                "pred_data", 
                                                "output_path"))

scp = readRDS(opt$input_object)
pred_data = readRDS(opt$pred_data)

# get predictions 
scp = scPredict(scp, newData = as.matrix(pred_data), threshold = opt$threshold_level)
predictions = getPredictions(scp)
write.csv(predictions, file=opt$output_path, row.names = TRUE, col.names = TRUE)

# if test labels supplied, run model performance evaluation block 
if(!is.na(opt$test_labels)){
    test_labels = read.csv(opt$test_labels, header = TRUE)
    row.names(test_labels) = test_labels[,1]
    test_labels = test_labels[,-1]
    scp@predMeta = test_labels
    conf_table = crossTab(scp, true = opt$cell_types_column)
    write.table(conf_table, file = opt$confusion_table, sep="\t")
    png(opt$plot_path)
    print(plotPredProbs(scp))
    dev.off()
}

