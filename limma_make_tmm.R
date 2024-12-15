# David Zemmour
# R
# usage: Rscript limma_make_tmm.R [path_to_seurat_object] [output_dir] [tmm_file_name]

options(max.print=1000)

args = commandArgs(TRUE)
path_to_seurat_object = args[1] 
output_dir = args[2] 
tmm_file_name = args[3]
# prefix2 = args[4]

# path_to_wd = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4"
# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"

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
so_orig = readRDS(file = path_to_seurat_object)

message("Remove TCR, mt, ribo, Gm and Rik genes and genes that are not expressed")
tcr_genes = rownames(so)[grepl(x = rownames(so), pattern = "Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|Trgc|Trdv|Trdj|Trdc")]
gm_rik_genes = rownames(so)[grepl(x = rownames(so), pattern = "Gm|Rik$")]
ribo_genes = rownames(so)[grepl(x = rownames(so), pattern = "Rpl|Rps|Mrpl|Mrps|Rsl")]
mt_genes = rownames(so)[grepl(x = rownames(so), pattern = "^mt-")]
genes_not_expressed = colSums(so[["RNA"]]$counts) == 0
genes_to_keep = !(rownames(so) %in% c(genes_not_expressed, tcr_genes, gm_rik_genes, ribo_genes, mt_genes))
length(which(genes_to_keep))

#Make tmm objet for limma
count = so[["RNA"]]$counts
dge = DGEList(count)
dge = calcNormFactors(dge)
tmm = new("EList")
message("TMM normalization")
tmm$E = edgeR::cpm(dge, log = TRUE, prior.count = 0.1)

saveRDS(tmm, file = tmm_file_name)