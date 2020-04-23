#!/usr/bin/env bash 
script_name=$0

# This is a test script designed to test that everything works in the various
# accessory scripts in this package. Parameters used have absolutely NO
# relation to best practice and this should not be taken as a sensible
# parameterisation for a workflow.

function usage {
    echo "usage: scPred_post_install_tests.sh [action] [use_existing_outputs]"
    echo "  - action: what action to take, 'test' or 'clean'"
    echo "  - use_existing_outputs, 'true' or 'false'"
    exit 1
}

action=${1:-'test'}
use_existing_outputs=${2:-'false'}

if [ "$action" != 'test' ] && [ "$action" != 'clean' ]; then
    echo "Invalid action"
    usage
fi

if [ "$use_existing_outputs" != 'true' ] && [ "$use_existing_outputs" != 'false' ]; then
    echo "Invalid value ($use_existing_outputs) for 'use_existing_outputs'"
    usage
fi

test_working_dir=`pwd`/'post_install_tests'
output_dir=$test_working_dir/outputs

# Clean up if specified
if [ "$action" = 'clean' ]; then
    echo "Cleaning up $test_working_dir ..."
    rm -rf $test_working_dir
    exit 0
elif [ "$action" != 'test' ]; then
    echo "Invalid action '$action' supplied"
    exit 1
fi 

# Initialise directories
mkdir -p $test_working_dir
mkdir -p $output_dir

################################################################################
# List tool outputs/inputs & parameters 
################################################################################
export test_sce=$test_working_dir/'pollen_cpm.rds'
export train_id='E-ENAD-16'
export train_matrix=$output_dir/'train_matrix.mtx'
export test_matrix=$output_dir/'test_matrix.mtx'
export train_metadata=$output_dir/'train_metadata.txt'
export test_metadata=$output_dir/'test_metadata.txt'
export scPred_object=$output_dir/'scPred_object.rds'
export scPred_feat_space=$output_dir/'scPred_feat_space.rds'
export eigenvalue_plot=$output_dir/'eigenvalue_plot.png'
export scPred_trained=$output_dir/'scPred_trained.rds'
export training_results=$output_dir/'training_results.rds'
export train_probs_plot=$output_dir/'train_probs_plot.png'
export predict_probs_plot=$output_dir/'prediction_probs_plot.png'
export predictions_output=$output_dir/'predictions_table.txt'
export confusion_table=$output_dir/'confusion_table.txt'
export scpred_output_tbl=$output_dir/'scpred_output_tbl.txt'

### Workflow parameters
export cell_types_column='cell_type2'
export normalised_counts_slot='normcounts'
export training_ratio=0.7
export log_trainsform='TRUE'
export n_pr_components=10
export explained_var_limit=0.01
export correction_method='fdr'
export model='svmRadial'
export resample_method='cv'
export iter_num=5
export classification_threshold=0.9

################################################################################
# Test individual scripts
################################################################################

# Make the script options available to the tests so we can skip tests e.g.
# where one of a chain has completed successfullly.

export use_existing_outputs

# Derive the tests file name from the script name
tests_file="${script_name%.*}".bats

# Execute the bats tests
$tests_file
