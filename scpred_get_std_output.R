#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

### create final workflow output in standard format
option_list = list(
    make_option(
        c("-i", "--predictions-object"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the Seurat predictions object in .rds format'
    ),
    make_option(
        c("-s", "--get-scores"),
        action = 'store_true',
        default = TRUE,
        type = 'logical',
        help = 'Should score column be added? Default: TRUE'
    ),
    make_option(
        c("-k", "--classifier"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the classifier object in .rds format (Optional; required to add dataset of origin to output table)'
    ),
    make_option(
        c("-o", "--output-table"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the final output file in text format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("predictions_object", "output_table"))
data = readRDS(predictions_object)
pred_labels = data$scpred_prediction
barcodes = names(pred_labels)

output = cbind(barcodes, unname(pred_labels))
if(opt$get_scores){
    scores = unname(data$scpred_max)
    output = cbind(output, scores)
    col_names = c('cell_id', "predicted_label", "score")
} else{
    col_names = c('cell_id', "predicted_label")
}
colnames(output) = col_names

# get rid of the goddamn dots
output[, "predicted_label"] = sapply(output[, "predicted_label"], function(x) x = gsub(pattern = '.', replacement = ' ', x, fixed = TRUE))
# add metadata if classifier is specified 
append = FALSE
if(!is.na(opt$classifier)){
    append = TRUE
    cl = readRDS(opt$classifier)
    dataset = attributes(cl)$dataset
    system(paste("echo '# tool scpred' >", opt$output_table))
    system(paste("echo '# dataset'", dataset, ">>", opt$output_table))
}
write.table(output, file = opt$output_table, sep="\t", row.names=FALSE, append=append)
