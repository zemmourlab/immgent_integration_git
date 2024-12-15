# David Zemmour
# R
# usage: Rscript limma_sample_seuratobject.R [path_to_seurat_object] [output_dir] [so_file_name]

options(max.print=1000)

args = commandArgs(TRUE)
path_to_seurat_object = args[1] 
output_dir = args[2] 
so_file_name = args[3] 

# path_to_wd = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4"
# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"
# output_dir = "DGE_limma/20241214"
# so_file_name = "igt1_96_CD4_20241113_sampled.Rds"

folder_path = sprintf("%s/", output_dir)
if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
} else {
    cat("Folder already exists:", folder_path, "\n")
}

setwd(output_dir)

message("loading R libraries and custom functions")
libs = c("Seurat","BPCells", "dplyr") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

message("loading seurat object")
so_orig = readRDS(file = sprintf("%s", path_to_seurat_object))

message("Sample cells so that no more than 20 cell per cluster per sample")
sampled_cells = so_orig@meta.data %>%
    group_by(IGTHT, annotation_level2) %>%
    group_modify(~ slice_sample(.x, n = min(20, nrow(.x)), replace = FALSE)) %>%
    ungroup() %>%
    pull(colnames)
length(sampled_data)
write.table(sampled_cells, file = sprintf("%s/sampledcells_IGTHT_annotationlevel2.csv", "."), row.names = F, col.names = F, sep = ",", quote = F)
so = so_orig[,sampled_cells]
saveRDS(so, file = so_file_name)
