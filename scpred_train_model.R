#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(doParallel))

# Use principal component-projected data and selected features to train a specified classification model

option_list = list(
    make_option(
        c("-i", "--input-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input object of Seurat class in .rds format'
  ),
    make_option(
        c("-f", "--train-id"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'ID of the training dataset (optional)'
  ),
    make_option(
        c("-m", "--model"), 
        action = "store",
        default = 'svmRadial',
        type = 'character',
        help = 'Model type used for training. Must be one of the models supported by Caret package. 
                Default: svmRadial'
  ),
    make_option(
        c("-r", "--resample-method"), 
        action = "store",
        default = 'cv',
        type = 'character',
        help = 'Resampling method used for model fit evaluation'
  ),
    make_option(
        c("-n", "--iter-num"), 
        action = "store",
        default = 5,
        type = 'numeric',
        help = 'Number of resampling iterations. Default: 5'
  ),
    make_option(
        c("-s", "--random-seed"), 
        action = "store",
        default = 123,
        type = 'numeric',
        help = 'Random seed'
  ),
    make_option(
        c("-p", "--allow-parallel"), 
        action = "store",
        default = TRUE,
        type = 'logical',
        help = 'Should parallel processing be allowed? Default: TRUE'
  ),
    make_option(
        c("-c", "--num-cores"), 
        action = "store",
        default = 1,
        type = 'numeric',
        help = 'For parallel processing, how many cores should be used?'
  ),
    make_option(
        c("-t", "--tune-length"), 
        action = "store",
        default = 3,
        type = 'numeric',
        help = 'An integer denoting the amount of granularity in the tuning parameter grid'
  ),
    make_option(
        c("-a", "--metric"), 
        action = "store",
        default = "ROC,PR,Accuracy,Kappa",
        type = 'character',
        help = "Performance metric to be used to select best model"
  ), 
    make_option(
        c("-e", "--preprocess"), 
        action = "store",
        default = "center,scale",
        type = 'character',
        help = "A string vector that defines a pre-processing of the predictor data. Enter values as comma-separated string. Current possibilities are 
                'BoxCox', 'YeoJohnson', 'expoTrans', 'center', 'scale', 'range', 'knnImpute', 'bagImpute', 'medianImpute' 
                'pca', 'ica' and 'spatialSign'. The default is 'center' and 'scale'."
  ), 
    make_option(
        c("-u", "--return-data"), 
        action = "store",
        default = FALSE,
        type = 'logical',
        help = 'If TRUE, training data is returned within scPred object. Default: FALSE'
  ),
     make_option(
        c("-v", "--save-predictions"), 
        action = "store",
        default = "final",
        type = 'character',
        help = "Specifies the set of hold-out predictions for each resample that should be
                returned. Values can be either 'all', 'final' or 'none'."
  ), 
     make_option(
        c("-y", "--reclassify"), 
        action = "store",
        default = NULL,
        type = 'character',
        help = "Cell types to reclassify using a different model"
  ), 
    make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for the output scPred object in .rds format'
  ), 
    make_option(
        c("-g", "--get-scpred"), 
        action = "store",
        default = FALSE,
        type = 'logical',
        help = 'Should scpred object be extracted from Seurat object after model training? Default: FALSE'
  ),
    make_option(
        c("-d", "--train-probs-plot"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for training probabilities plot in .png format'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_object", "output_path"))
preprocess = wsc_split_string(opt$preprocess, ",")
metric = wsc_split_string(opt$metric, ",")

if(!is.null(opt$reclassify)){
    cells_to_reclassify = wsc_split_string(opt$preprocess, ",")
} else {
    cells_to_reclassify = NULL
}


data_seurat = readRDS(opt$input_object)


# model training step 
clust = makePSOCKcluster(opt$num_cores)
registerDoParallel(clust)
classifier = trainModel(data_seurat, 
                 seed = opt$random_seed, 
                 model = opt$model,
                 resampleMethod = opt$resample_method, 
                 number = opt$iter_num, 
                 allowParallel = opt$allow_parallel,
                 preProcess = preprocess,
                 tuneLength = opt$tune_length,
                 metric = metric,
                 returnData = opt$return_data,
                 savePredictions = opt$save_predictions,
                 reclassify = NULL
                 )
stopCluster(clust)

if(opt$get_scpred){
    classifier = get_scpred(classifier)
}

# plot class probs
if(!is.na(opt$train_probs_plot)){
    png(opt$train_probs_plot)
    print(plot_probabilities(classifier))
    dev.off()
}

# add dataset field to the object 
if(!is.na(opt$train_id)){
    attributes(classifier)$dataset = opt$train_id
    } else{
        attributes(classifier)$dataset = NA
    }

saveRDS(classifier, opt$output_path)
