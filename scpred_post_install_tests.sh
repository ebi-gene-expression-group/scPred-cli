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
export ref_seurat=$test_working_dir/'reference_pbmc.rds'
export query_seurat=$test_working_dir/'query_pbmc.rds'
export scPred_feat_space=$output_dir/'scPred_feat_space.rds'
export scPred_trained=$output_dir/'scPred_trained.rds'

export train_probs_plot=$output_dir/'train_probs_plot.png'
export predict_probs_plot=$output_dir/'prediction_probs_plot.png'
export predictions_output=$output_dir/'predicted_data.rds'
export scpred_output_tbl=$output_dir/'scpred_output_tbl.txt'

### Workflow parameters
export train_id='E-ENAD-16'
export norm_query_data="TRUE"
export get_scores="TRUE"
export num_cores=2

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
