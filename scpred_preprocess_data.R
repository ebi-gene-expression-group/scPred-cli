#!/usr/bin/env Rscript 


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scPred))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(dplyr))

### Pre-process the initial Seurat object
### The following steps are performed:
###     - Normalise data (if specified)
###     - Find relevant feature genes
###     - Scale data (subtract mean expression and divide by the standard deviation)
###     - Run PCA
###     - Generate UMAP plot (if specified)


# argument parsing 
option_list = list(
    make_option(
    c("-i", "--raw-seurat-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input Seurat object in .rds format'
  ),
    make_option(
    c("-n", "--normalise-data"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Should the expression data be normalised? Default: False'
  ),
    make_option(
    c("-m", "--normalisation-method"),
    action = "store",
    default = "LogNormalize",
    type = 'character',
    help = 'If --normalise-data specified, what normalisation method to use? Default: LogNormalize'
  ),
    make_option(
    c("-s", "--scale-factor"),
    action = "store",
    default = 10000,
    type = 'numeric',
    help = 'If --normalise-data specified, what scale factor should be applied? 
            Note: for PCM normalisation, select 1e6'
  ),
    make_option(
    c("-t", "--selection-method"),
    action = "store",
    default = "vst",
    type = 'character',
    help = 'Selection method to find feature genes'
  ),
    make_option(
    c("-m", "--processed-seurat-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the processed Seurat object in .rds format'
  ),
    make_option(
    c("-d", "--umap-plot-dim"),
    action = "store",
    default = 30,
    type = 'numeric',
    help = 'Number of dims for UMAP plot. Default: 30'
  ),
    make_option(
    c("-c", "--cell-type-col"),
    action = "store",
    default = "cell_type",
    type = 'character',
    help = 'Name of cell type column in object metadata. Default: "cell_type"'
  ),
    make_option(
    c("-p", "--umap-plot-output-path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the UMAP plot file in .png format'
  )
)

opt = wsc_parse_args(option_list, mandatory=c("raw_seurat_object", "processed_seurat_object"))
data_seurat = readRDS(opt$raw_seurat_object)

# normalise data, if specified
if(opt$normalise_data){
    data_seurat = NormalizeData(object = data_seurat,
                                normalization.method = opt$normalisation_method,
                                scale.factor = opt$scale_factor)
}


# find feature genes; scale; run PCA
data_seurat = data_seurat %>%
          FindVariableFeatures(selection.method = opt$selection_method) %>%
          ScaleData() %>%
          RunPCA()

# generate UMAP plot, if specified
if(!is.na(opt$umap_plot_output_path)){
    data_seurat = RunUMAP(object = data_seurat, dims = 1:opt$umap_plot_dim)
    png(opt$umap_plot_output_path)
    print(DimPlot(data_seurat, group.by = opt$cell_type_col, label = TRUE, repel = TRUE))
    dev.off()
}




