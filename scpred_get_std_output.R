#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

### create final workflow output in standard format
option_list = list(
    make_option(
        c("-i", "--predictions-file"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the predictions file in text format'
    ),
    make_option(
        c("-b", "--cell-id-col"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Cell id column name. If not provided, it is assumed cell ids are represented by index'
    ),
    make_option(
        c("-s", "--get-scores"),
        action = 'store_true',
        default = FALSE,
        type = 'logical',
        help = 'Should score column be added? Default: FALSE'
    ),
    make_option(
        c("-o", "--output-table"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the final output file in text format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("predictions_file", "output_table"))
data = read.csv(opt$predictions_file, row.names = 1)
if(!is.na(opt$cell_id_col)){
    barcodes = as.character(data[, opt$cell_id_col])
    } else{
        barcodes = row.names(data)
    }

output = cbind(cell_id=barcodes, pred_label=as.character(data[, "predClass"]))
if(opt$get_scores){
    scores = data[,-ncol(data)]
    .get_max = function(vec){
        vec = as.numeric(vec)
        return(max(vec))
    }
    scores = apply(scores, 1, .get_max)
    output = cbind(output, scores)
    col_names = c('cell_id', "predicted_label", "score")
} else{
    col_names = c('cell_id', "predicted_label")
}

colnames(output) = col_names
# get rid of the goddamn dots
output[, "predicted_label"] = sapply(output[, "predicted_label"], function(x) x = gsub(pattern = '.', replacement = ' ', x, fixed = TRUE))
write.table(output, file = opt$output_table, sep="\t", row.names=FALSE)
