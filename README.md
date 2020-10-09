# scPred-scripts
This is a collection of R scripts to allow the execution of the [scPred](https://github.com/powellgenomicslab/scPred) cell-type classification method.

## Installation
The simplest way to install package is to use conda. It is recommended to install into clean environment to avoid dependency conflicts.
```
conda create -n scpred-cli 
conda activate scpred-cli 
conda install scpred-cli 
```

### Testing
To run post-install tests, you will need to install [bats](https://github.com/sstephenson/bats) testing system into the environment you created:
```
conda install bats
```
Then, run `scpred_post_install_tests.sh`

## Commands

Currently wrapped scPred functions are described below. Each script has usage instructions available via --help, consult function documentation in scPred for further details.

### Extract test data
Input to the workflow will be a serialized single-cell experiment object. A test dataset can be downloaded by running the following script.

```
get_test_data.R
```

### Get feature space
Select principal components that will be used as features for training the model.

```
scpred_get_feature_space.R\
    --input-object < Path to the input object of scPred or seurat class in .rds format >\
    --prediction-column < Name of the metadata column that contains training labels >\
    --correction method < Multiple testing correction method. Default: fdr >\
    --significance-threshold < Significance threshold for principal components explaining class identity >\ 
    --reduction-parameter < Name of reduction in Seurat objet to be used to determine the feature space. Default: "pca" >\
    --output-path < Path for the modified scPred object in .rds format >
```

### Train classification model
Use principal component-projected data and selected features to train a specified classification model.

```
scpred_train_model.R\
    --input-object < Path to the input object of scPred or seurat class in .rds format >\ 
    --train-idf <Path to the training data IDF file (optional)>\
    --model < Model type used for training. Must be one of the models supported by Caret package >\
    --resample-method < Resampling method used for model fit evaluation >\
    --iter-num < umber of resampling iterations. Default: 5 >\
    --allow-parallel < Should parallel processing be allowed? Default: TRUE >\
    --num-cores < For parallel processing, how many cores should be used? >\
    --tune-length < An integer denoting the amount of granularity in the tuning parameter grid >\
    --metric < Performance metric to be used to select best model >\
    --preprocess < A string vector that defines a pre-processing of the predictor data >\
    --reclassify < Cell types to reclassify using a different model >\
    --output-path < Path for the output scPred object in .rds format >\
    --get-scpred < Should scpred object be extracted from Seurat object after model training? Default: FALSE >\
    --train-probs-plot < Path for training probabilities plot in .png format >
```
### Obtain predictions using trained model
Make cell type predictions using trained model this script can be used both for evaluation of model performance on test data and obtaining predictions on new data.

```
scpred_predict.R\ 
    --input-object < Path to the input object of scPred or seurat class in .rds format >\ 
    --pred-data < Path to the input prediction object in .rds format >\
    --normalise-data < Should the predicted expression data be normalised? Default: False >\
    --normalisation-method < If --normalise-data specified, what normalisation method to use? Default: LogNormalize
                             NB: normalisation method must be identical to that used for reference data >\
    --scale-factor < If --normalise-data specified, what scale factor should be applied? >\
    --threshold-level < Classification threshold value >\
    --output-path < Output path for values predicted by the model in text format >\ 
    --plot-path < Output path for prediction probabilities histograms in .png format >\ 

```
### Get standardised output table 
This script allows to get predicted labels in a standardised format, simplifying downstream analyses. 
```
scpred_get_std_output.R\
    --predictions-object <Path to the predictions object in .rds format>\
    --get-scores <Boolean: should the prediction scores be included? default: FALSE>\
    --classifier <Path to the classifier object in .rds format>\
    --output-table <Path to the final output file in text format>\
```
