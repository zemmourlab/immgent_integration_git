# David Zemmour
# R
# usage: Rscript limma_sample_seuratobject.R [path_to_seurat_object] [output_dir] [so_file_name]

options(max.print=1000)

# Parse arguments
args = commandArgs(TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript limma_sample_seuratobject.R [path_to_seurat_object] [output_dir] [so_file_name]")
}
path_to_seurat_object = args[1] 
output_dir = args[2] 
so_file_name = args[3] 

# path_to_wd = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4"
# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"
# output_dir = "DGE_limma/20241214"
# so_file_name = "igt1_96_CD4_20241113_sampled.Rds"

# Validate inputs
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Folder created:", output_dir, "\n")
} else {
    cat("Folder already exists:", output_dir, "\n")
}

if (!file.exists(path_to_seurat_object)) {
    stop("The specified Seurat object file does not exist: ", path_to_seurat_object)
}

# setwd(output_dir)

message("loading R libraries and custom functions")
libs = c("Seurat","BPCells", "dplyr") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

message("Loading Seurat object")
so_orig = readRDS(file = sprintf("%s", path_to_seurat_object))

message("Sample cells so that no more than 20 cell per cluster per sample")
sampled_cells = so_orig@meta.data %>%
    group_by(IGTHT, annotation_level2) %>%
    group_modify(~ slice_sample(.x, n = min(20, nrow(.x)), replace = FALSE)) %>%
    ungroup() %>%
    pull(cellID)
if (length(sampled_cells) == 0) {
    stop("No cells were sampled. Check the input Seurat object and metadata.")
}
cat("Number of sampled cells:", length(sampled_cells), "\n")
write.table(sampled_cells, file = sprintf("%s/sampledcells_IGTHT_annotationlevel2.csv", "."), row.names = F, col.names = F, sep = ",", quote = F)
so = so_orig[,sampled_cells]
saveRDS(so, file = sprintf("%s/%s",output_dir,so_file_name))
message("Saved sampled Seurat object to:", sprintf("%s/%s",output_dir,so_file_name))
