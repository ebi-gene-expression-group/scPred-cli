#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))

# Use principal component-projected data and selected features to train a specified classification model

option_list = list(
    make_option(
        c("-i", "--input-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input object of scPred or seurat class in .rds format'
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
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for the output scPred object in .rds format'
  ), 
    make_option(
        c("-t", "--training-results"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path for training step results object in .rds format'
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
scp = readRDS(opt$input_object)
# model training step 
scp = trainModel(scp, 
                 seed = opt$random_seed, 
                 model = opt$model,
                 resampleMethod = opt$resample_method, 
                 number = opt$iter_num, 
                 allowParallel = opt$allow_parallel)
# obtain training results 
if(!is.na(opt$training_results)){
    res = getTrainResults(scp)
    saveRDS(res, opt$training_results)
}

# plot class probs
if(!is.na(opt$train_probs_plot)){
    png(opt$train_probs_plot)
    print(plotTrainProbs(scp))
    dev.off()
}

# add dataset field to the object 
if(!is.na(opt$train_id)){
    attributes(scp)$dataset = opt$train_id
    } else{
        attributes(scp)$dataset = NA
    }

saveRDS(scp, opt$output_path)
