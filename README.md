# scPred-cli
Command-line interface for the [scPred](https://github.com/IMB-Computational-Genomics-Lab/scPred) cell type classifaction package.

### Installation 


To run post-install tests, enter the following command: `scpred_post_install_tests.sh 'test'` 


### Commands 
A full set of parameters and their description can be obtained by running `<command>.R --help`. 

#### Train/test split 
Split the expression matrix into training and testing datasets.
```
scpred_train_test_split.R\
          --input-sce-object <path to the input object of SingleCellExperiment class in .rds format>\
          --training-matrix <output path for training matrix>\
          --test-matrix <output path for testing matrix>\
          --training_labels <output path for training labels>\
          --test-labels <output path for test labels>\
          --cell-types-column <column name for assigned cell types in the metadata slot of SCE object>
```

#### Eigenvalue decomposition of training matrix
Calculate n first principal components from the data matrix and build an object of scPred class.
```
scpred_eigen_decomp.R\
          --training-matrix <path to the training matrix in .rds format>\
          --log-transform <should log-transform be performed on the matrix data? Default: TRUE>\
          --training-labels <path to the training labels (metadata) in text format>\
          --output-path <output path for the scPred object in .rds format>
```

#### Feature selection 
Get representative features from the set of principal components.
```
scpred_get_feature_space.R\
          --input-object <path to the input object of scPred or seurat class in .rds format>\
          --prediction-column <name of the metadata column that contains training labels>\
          --output-path <path for the modified scPred object in .rds format>\
          --eigenvalue-plot-path <path for eigenvalue plot for principal components in .png format>
```

#### Train the model 
Use principal component-projected data to train a specified classification model.
```
scpred_train_model.R\
          --input-object <path to the input object of scPred of seurat class in .rds format>\
          --model <model type used for training. Must be one of the models supported by Caret package. Default: svmRadial>\
          --output-path <path for the output scPred object in .rds format>\
          --train-probs-plot <path for training probabilities plot in .png format>
```

#### Obtain predictions 
Use the trained model either to predict cell types for unseen data or evaluate performance on test data.
```
scpred_predict.R\
          --input-object <path to the input object of scPred or seurat class in .rds format>\
          --pred-data <path to the input object for prediction matrix in .rds format>\
          --test-labels <optional: path to the test labels for evalutation of model performance>\
          --cell-types-column <column name of true labels in provided labels (metadata) file>\
          --output-path <output path for predicted values in .txt format>\
          --plot-path <output path for prediction probabilities histograms in .png format>\
          --confusion-table <output path for confusion table in text format>

```





