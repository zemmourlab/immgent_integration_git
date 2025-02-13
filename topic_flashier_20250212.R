# David Zemmour
# R
# Usage: Rscript topic_flashier_20250203.R [path_to_seurat_object] [output_dir] [backfit (True/False)]
# LOOK FOR EDIT to find where to edit the script! 

options(max.print=1000)
options(expressions = 50000)

# Parse arguments
args = commandArgs(TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript topic_flashier_20250203.R [path_to_seurat_object] [output_dir] [backfit True/False]")
}
path_to_seurat_object = args[1]
output_dir = args[2] 
backfit_option = tolower(args[3])

# Validate and convert backfit_option
if (backfit_option %in% c("true", "t")) {
    backfit_option=  TRUE
} else if (backfit_option %in% c("false", "f")) {
    backfit_option = FALSE
} else {
    stop("Error: backfit should be 'True' or 'False'.")
}

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

message("loading R libraries and custom functions")
libs = c("fastTopics", "flashier", "Matrix", "Seurat","BPCells") 

sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

message("loading seurat object")
so = readRDS(file = path_to_seurat_object)

message("log1p with scale factor as the average rowSums(counts)")
norm_fac = mean(so$nCount_RNA)
so = NormalizeData(so, assay = "RNA", normalization.method = "LogNormalize", scale.factor = norm_fac)
counts = so[['RNA']]$counts
shifted_log_counts = so[['RNA']]$data

#transpose!
counts = t(counts)
shifted_log_counts = t(shifted_log_counts)

message("Removing TCR, mt, ribo, Gm, and Rik genes and genes that are not expressed")
tcr_genes = grepl(x = rownames(so), pattern = "Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|Trgc|Trdv|Trdj|Trdc")
gm_rik_genes = grepl(x = rownames(so), pattern = "Gm|Rik$|\\-ps$")
ribo_genes = grepl(x = rownames(so), pattern = "Rpl|Rps|Mrpl|Mrps|Rsl")
mt_genes = grepl(x = rownames(so), pattern = "^mt-")
genes_not_expressed = rowSums(so[["RNA"]]$counts) == 0
genes_to_keep = !(tcr_genes | gm_rik_genes | ribo_genes | mt_genes | genes_not_expressed)
cat("Number of genes to keep:", sum(genes_to_keep), "\n")

shifted_log_counts = shifted_log_counts[,genes_to_keep]
counts = counts[,genes_to_keep]

message("Variance regularization")
n = nrow(counts)
x = rpois(1e7, 1/n)
s1 = sd(log(x + 1))

message("Fitting")
fit = flash(shifted_log_counts, 
            ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
            var_type = 2, 
            S = s1,
            backfit = backfit_option)

saveRDS(fit, file = sprintf("%s/fit.Rds", output_dir))




