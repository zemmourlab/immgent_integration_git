### Pre-process paper data ###

library(dplyr)
library(Seurat)
library(SeuratObject)
library(MASS)
library(readxl)
library(DropletUtils)

so_path = "N:/CBDM_Lab/Odhran/scVerse/ImmgenT_Workshop/unconventional/Tgd/"
so_name = "gd_holygrail_filtered_final.rds"
batche_of_interest = "batch1_oct"

# for seurat object
so <- readRDS(paste0(so_path, so_name))
View(so@meta.data)
so <- so[,so$batch ==batch_of_interest]
query_df = as.data.frame(so@assays$RNA@counts)
genes_query = rownames(so)
cells_query = colnames(so)

# make genes matrix
sc <- readRDS("N:/CBDM_Lab/scRNAseq/ImmgenT/IGT92/dataset_clean.Rds")
genes <- rownames(sc)
query_df = query_df[rownames(query_df) %in% genes,]
gene_diff = length(genes) - length(rownames(query_df))
missing_genes <- setdiff(genes, rownames(query_df))

zero_rows <- as.data.frame(matrix(0, nrow = gene_diff, ncol = length(cells_query)))
colnames(zero_rows) <- colnames(query_df)
rownames(zero_rows) <- missing_genes
query_df <- as.data.frame(query_df)
query_df_all = rbind(query_df, zero_rows)

query_df_reorder = query_df_all[genes,]

Query_matrix_finished = as(query_df_reorder, "sparseMatrix")

meta = so@meta.data
write.table(meta, paste0(so_path, "cell_metadata.csv"), quote = F, row.names = T, col.names = NA, sep = ",")

write10xCounts(
  paste0(so_path, "matrix"),
  Query_matrix_finished,
  barcodes = cells_query,
  gene.id = genes,
  gene.symbol = genes,
  overwrite = TRUE
)
