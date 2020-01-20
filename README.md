# scPred-scripts
This is a collection of R scripts to allow the execution of the [scPred](https://github.com/powellgenomicslab/scPred) cell-type classification method.

## Commands

Currently wrapped scPred functions are described below. Each script has usage instructions available via --help, consult function documentation in scPred for further details.

### Extract test data
Input to the workflow will be a serialized single-cell experiment object. A test dataset can be downloaded from [here](https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/pollen.rds) or equivalently running the following script.

```
get_test_data.R
```

### Split dataset into train and test 
Input dataset is pre-processed (CPM normalization) and partitioned with a size of the training subset specified in **training-ratio** (70% by default).

```
scpred_train_test_split.R 
--input-sce-object < Path to the input SCE object in .rds format > \
--normalised-counts-slot < Name of the slot with normalised counts matrix in SCE object. Default: normcounts > \ 
--training-matrix < Output path for training matrix in .rds format > \ 
--test-matrix < Output path for test matrix in .rds format > \ 
--cell-types-column < Column name for assigned cell types > \
--training-labels < Output path for training labels in text format > \  
--test-labels < Output path for test labels in text format > \
--training-ratio < Proportion of training/testing dataset split > \ 
```

### Eigenvalue-decomposition of training matrix
Calculate n first principal components and apply log-transformation to the matrix if specified; initialize object of scPred class. 

```
scpred_eigen_decomp.R 
--training-matrix < Path to the training matrix in .rds format > \
--log-transform < Should log-transform be performed on the matrix? Default: TRUE > \ 
--training-labels < Path to the training labels (metadata) in text format > \ 
--principal-comps < Number of pricipal components for eigenvalue decomposition > \ 
--output-path < Output path for the scPred object in .rds format > \
```

### Get feature space
Select principal components that will be used as features for training the model.

```
scpred_get_feature_space.R 
--input-object < Path to the input object of scPred or seurat class in .rds format > \
--prediction-column < Name of the metadata column that contains training labels> \
--explained-var-limit < Threshold to filter principal components based on variance explained > \
--output-path < Path for the modified scPred object in .rds format > \
--eigenvalue-plot-path < Path for eigenvalue plot for principal components in .png format >
```

### Train classification model
Use principal component-projected data and selected features to train a specified classification model.

```
scpred_train_model.R 
--input-object < Path to the input object of scPred or seurat class in .rds format > \ 
--output-path < Path for the output scPred object in .rds format > \
--train-probs-plot < Path for training probabilities plot in .png format > \
```
### Obtain predictions using trained model
Make cell type predictions using trained model this script can be used both for evaluation of model performance on test data and obtaining predictions on new data.

```
scpred_predict.R 
--input-object < Path to the input object of scPred or seurat class in .rds format > \ 
--pred-data < Path to the input prediction matrix in .rds format >\
--test-labels < Path to the test labels file for evalutation of model performance in text format > \ 
--cell-types-column < Column name of true labels in provided metadata file > \ 
--output-path < Output path for values predicted by the model in text format > \ 
--plot-path < Output path for prediction probabilities histograms in .png format > \ 
--confusion-table < Output path for confusion table in text format > \
```
