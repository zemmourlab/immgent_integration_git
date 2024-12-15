# David Zemmour
# R
# usage: Rscript limma_fit_level2.IGTHT.R [path_to_seurat_object] [path_to_tmm_object] [output_dir] [fit_file_name]

options(max.print=1000)

args = commandArgs(TRUE)
path_to_seurat_object = args[1]
path_to_tmm_object = args[2]
output_dir = args[3] 
tmm_file_name = args[4]
# prefix2 = args[4]

# path_to_wd = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4"
# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"
# prefix = gsub(pattern = ".Rds", replacement = "", x = path_to_seurat_object)

folder_path = sprintf("%s/", output_dir)
if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
} else {
    cat("Folder already exists:", folder_path, "\n")
}

setwd(output_dir)

message("loading R libraries and custom functions")
libs = c("limma", "edgeR", "Seurat","BPCells", "dplyr", "rlang","reshape2") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

message("loading seurat object")
so = readRDS(file = path_to_seurat_object)

message("loading tmm object")
tmm = readRDS(file = path_to_tmm_object)

message("Fitting")
design = data.frame(so@meta.data[,c("annotation_level2","IGTHT")]) 
design$annotation_level2.IGTHT = paste(design$annotation_level2, design$IGTHT, sep = ".")
all(rownames(design) == colnames(so))
design = model.matrix(~ 0 + annotation_level2.IGTHT, data=design) 
colnames(design) = gsub("annotation_level2.IGTHT", "", colnames(design))
#cat("ncol=",ncol(design),"rank=", qr(design)$rank,"\n")
fit = lmFit(tmm, design = design)
saveRDS(fit, file = fit_file_name)