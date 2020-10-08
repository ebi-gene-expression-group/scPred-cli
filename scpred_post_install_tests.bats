#!/usr/bin/env bats 

# download test sce object from the link provided in package docs
@test "extract test data" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$ref_seurat" ]; then
        skip "$test_sce exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $test_sce && get_test_data.R
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$ref_seurat" ]
    [ -f "$query_seurat" ]
  
}

@test "get feature space" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scPred_feat_space" ]; then
        skip "$scPred_feat_space exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scPred_feat_space  &&\
                                 scpred_get_feature_space.R\
                                        --input-object $ref_seurat\
                                        --output-path $scPred_feat_space\

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
                                    --train-id $train_id\
				    --num-cores $num_cores\
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
                                    --pred-data $query_seurat\
                                    --normalise-data $norm_query_data\
                                    --output-path $predictions_output\
                                    --plot-path $predict_probs_plot\

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$predictions_output" ]
}

@test "Get standard output table" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scpred_output_tbl" ]; then
        skip "$scpred_output_tbl exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scpred_output_tbl && scpred_get_std_output.R\
                                        --predictions-object $predictions_output\
                                        --get-scores $get_scores\
                                        --classifier $scPred_trained\
                                        --output-table $scpred_output_tbl

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scpred_output_tbl" ]
}
