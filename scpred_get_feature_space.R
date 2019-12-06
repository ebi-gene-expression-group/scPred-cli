#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))

# Select principal components that will be used as features for training the model
option_list = list(
    make_option(
        c("-i", "--input-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input object of scPred or seurat class in .rds format'
  ),
    make_option(
        c("-p", "--prediction-column"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Name of the metadata column that contains training labels'
  ),
     make_option(
        c("-v", "--explained-var-limit"), 
        action = "store",
        default = 0.01,
        type = 'numeric',
        help = 'Threshold to filter principal components based on variance explained'
  ),
     make_option(
        c("-c", "--correction-method"), 
        action = "store",
        default = 'fdr',
        type = 'character',
        help = 'Multiple testing correction method. Default: fdr'
  ),
     make_option(
        c("-s", "--significance-threshold"), 
        action = "store",
        default = 0.05,
        type = 'numeric',
        help = 'Significance threshold for principal components explaining class identity'
  ),
      make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for the modified scPred object in .rds format'
  ), 
      make_option(
        c("-e", "--eigenvalue-plot-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for eigenvalue plot for principal components in .png format'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_object",
                                                "prediction_column",
                                                "output_path"))
scp = readRDS(opt$input_object)
scp = getFeatureSpace(scp, 
                      pVar = opt$prediction_column, 
                      varLim = opt$explained_var_limit, 
                      correction = opt$correction_method, 
                      sig = opt$significance_threshold)

saveRDS(scp, opt$output_path)
# if plot path supplied, produce eigenvalue plots 
if(!is.na(opt$eigenvalue_plot_path)){
    png(opt$eigenvalue_plot_path)
    print(plotEigen(scp, group = opt$prediction_column))
    dev.off() 
}
