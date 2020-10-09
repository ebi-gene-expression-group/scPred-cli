#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))
suppressPackageStartupMessages(require(Seurat))

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
        c("-n", "--normalise-data"),
        action = "store",
        default = FALSE,
        type = 'logical',
        help = 'Should the predicted expression data be normalised? Default: False'
  ),
    make_option(
        c("-m", "--normalisation-method"),
        action = "store",
        default = "RC",
        type = 'character',
        help = 'If --normalise-data specified, what normalisation method to use? Default: LogNormalize
                NB: normalisation method must be identical to that used for reference data'
  ),
    make_option(
        c("-s", "--scale-factor"),
        action = "store",
        default = 1e6,
        type = 'numeric',
        help = 'If --normalise-data specified, what scale factor should be applied? 
                Note: for CPM normalisation, select 1e6'
  ),
    make_option(
        c("-l", "--threshold-level"), 
        action = "store",
        default = 0.8,
        type = 'numeric',
        help = 'Classification threshold value'
  ),
    make_option(
        c("-x", "--max-iter-harmony"), 
        action = "store",
        default = 20,
        type = 'numeric',
        help = 'Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step'
  ),
    make_option(
        c("-r", "--recompute-alignment"),
        action = "store",
        default = TRUE,
        type = 'logical',
        help = 'Recompute alignment? Useful if scPredict() has already been run. Default: TRUE'
  ),
    make_option(
        c("-k", "--reference-scaling"),
        action = "store",
        default = TRUE,
        type = 'logical',
        help = 'Scale new dataset based on means and stdevs from reference dataset before alignment. Default: TRUE'
  ),
    make_option(
        c("-e", "--random-seed"), 
        action = "store",
        default = 66,
        type = 'numeric',
        help = 'Random number generator seed'
  ),
     make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for Seurat object holding predicted values'
  ),
     make_option(
        c("-a", "--plot-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for prediction probabilities histograms in .png format'
  ) 
)

opt = wsc_parse_args(option_list, mandatory = c("input_object",
                                                "pred_data", 
                                                "output_path"))

ref_data = readRDS(opt$input_object)
pred_data = readRDS(opt$pred_data)


# normalise query data, if specified
if(opt$normalise_data){
    pred_data = NormalizeData(object = pred_data,
                              normalization.method = opt$normalisation_method,
                              scale.factor = opt$scale_factor)
}

# get predictions 
pred_data = scPredict(reference = ref_data,
                      new = pred_data,
                      threshold = opt$threshold_level,
                      max.iter.harmony = opt$max_iter_harmony,
                      recompute_alignment = opt$recompute_alignment,
                      reference_scaling = opt$reference_scaling,
                      seed = opt$random_seed)
saveRDS(pred_data, opt$output_path)

# generate prediction plot 
if(!is.na(opt$plot_path)){
    png(opt$plot_path)
    print(DimPlot(pred_data, group.by = "scpred_prediction", reduction = "scpred"))
    dev.off()
}
