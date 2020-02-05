#!/usr/bin/env bats 

# download test sce object from the link provided in package docs
@test "extract test data" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$test_sce" ]; then
        skip "$test_sce exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $test_sce && get_test_data.R
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$test_sce" ]
  
}

@test "train/test split" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$train_matrix" ]; then
        skip "$train_matrix exists and use_existing_outputs is set to 'true'"
    fi
    run rm -f $train_matrix $test_matrix $train_metadata $test_metadata\
    && scpred_train_test_split.R\
                        --input-sce-object $test_sce\
                        --normalised-counts-slot $normalised_counts_slot\
                        --training-matrix $train_matrix\
                        --test-matrix $test_matrix\
                        --cell-types-column $cell_types_column\
                        --training-labels $train_metadata\
                        --test-labels $test_metadata\
                        --training-ratio $training_ratio

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$train_matrix" ]
}

@test "eigenvalue-decomposition of training matrix" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scPred_object" ]; then
        skip "$scPred_object exists and use_existing_outputs is set to 'true'"
    fi
    
    run rm -f $scPred_object && scpred_eigen_decomp.R\
                                    --training-matrix $train_matrix\
                                    --log-transform $log_trainsform\
                                    --training-labels $train_metadata\
                                    --principal-comps $n_pr_components\
                                    --output-path $scPred_object

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scPred_object" ]
}

@test "get feature space" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scPred_feat_space" ]; then
        skip "$scPred_feat_space exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scPred_feat_space $eigenvalue_plot &&\
                                 scpred_get_feature_space.R\
                                        --input-object $scPred_object\
                                        --prediction-column $cell_types_column\
                                        --explained-var-limit $explained_var_limit\
                                        --output-path $scPred_feat_space\
                                        --eigenvalue-plot-path $eigenvalue_plot

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scPred_feat_space" ]
}

@test "train classification model" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scPred_trained" ]; then
        skip "$scPred_trained exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scPred_trained $training_results $train_probs_plot &&\
                               scpred_train_model.R\
                                    --input-object $scPred_feat_space\
                                    --output-path $scPred_trained\
                                    --train-probs-plot $train_probs_plot

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scPred_trained" ]
    
}

@test "obtain predictions using trained model" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$predictions_output" ]; then
        skip "$predictions_output exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $predictions_output $predict_probs_plot $confusion_table &&\
                            scpred_predict.R\
                                    --input-object $scPred_trained\
                                    --pred-data $test_matrix\
                                    --test-labels $test_metadata\
                                    --cell-types-column $cell_types_column\
                                    --output-path $predictions_output\
                                    --plot-path $predict_probs_plot\
                                    --confusion-table $confusion_table

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scPred_trained" ]
}


@test "Get standard output table" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scpred_output_tbl" ]; then
        skip "$scpred_output_tbl exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scpred_output_tbl && scpred_get_std_output.R\
                                        --predictions-file $predictions_output\
                                        --get-scores\
                                        --output-table $scpred_output_tbl

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scpred_output_tbl" ]
}
