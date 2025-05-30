### Transfer labels - David's code



library(RcppAnnoy)
library(Seurat)
#install.packages("fastTopics")
library(ggplot2)

libs = c("fgsea","fastTopics", "flashier", "Matrix", "Seurat","BPCells", "ggplot2", "dplyr", "reshape2", "ggrastr", "RColorBrewer", "pals", "scales", "pheatmap") 
#libs = c("matrixStats", "gplots","ggplot2", "reshape2", "scales", "gridExtra", "dplyr", "RColorBrewer", "grid", "Rtsne", "limma",   "RColorBrewer", "pheatmap", "Seurat",  "ggrastr", "ggbeeswarm") #pwr, preprocessCore, "readxl", vegan", "genefilter" "scde" "cellrangerRkit", "Rsamtools", "GenomicRanges", "GenomicAlignments", "VGAM", "WGCNA",rafalib, ‘UsingR’,MAST "ggrastr",
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))


#Gradients
library(RColorBrewer)
ColorRamp = rev(colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100))
ColorRamp = rev(colorRampPalette(c('red','white','blue'))(20))
ColorRamp = rev(rainbow(10, end = 4/6))
library(viridis)
ColorRamp = rev(viridis(100))
ColorRamp = rev(cividis(100))
#image(1:100, 1, as.matrix(1:100), col = ColorRamp, main = "Color Palette", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")

#mycol = c("black","blue","red","darkgreen", "orange","purple", "cyan", "greenyellow",    "salmon", "magenta","pink", "tan", "brown") # c("magenta", "red",  "darkgreen", "cyan", "blue", "blue", "blue", "black", "blue", "blue", "blue", "blue", "orange")

#Large Color palette
library("pals")
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal1 = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal1 = mypal1[-4]
mypal = c(glasbey(), polychrome(), mypal1)
names(mypal) = NULL
#pal.bands(mypal, main="Colormap suggestions")


# Read in so with reference and query together
# sKIP IF CREATED ALREADY
so_integrated <- readRDS("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/igt1_96_withtotalvi20250109.Rds")
so_integrated@meta.data <- so_integrated@meta.data[, c("IGT", "IGTHT")]
IGT_97_104 <- readRDS("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/IGT97_104/IGT97_104_seurat_object.Rds")
IGT_97_104@meta.data <- IGT_97_104@meta.data[, c("IGT", "IGTHT")]
RNA_so <- merge(so_integrated, IGT_97_104)

# use previous RNA_so object with transfer labels
RNA_so <- readRDS("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/IGT97_104/IGT1_104_seurat_object_PostQC_20250505.Rds")
annotations_20250226 <- read.csv("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/annotation_table_20250226.csv", row.names = 1)
RNA_so$annotation_level1_20250226 <- "NA"
RNA_so$annotation_level1_20250226[rownames(annotations_20250226)] <- annotations_20250226$level1
RNA_so$annotation_level2_group_20250226 <- "NA"
RNA_so$annotation_level2_group_20250226[rownames(annotations_20250226)] <- annotations_20250226$level2.group
RNA_so$annotation_level2_20250226 <- "NA"
RNA_so$annotation_level2_20250226[rownames(annotations_20250226)] <- annotations_20250226$level2

# subset RNA_so and work with CD4 only
RNA_so <- readRDS("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/IGT97_104/Seurat_PostQC/seuratobject_merged_RNA_withmetadata_withreductions_transferlabels.Rds")
CD4 <- RNA_so[, RNA_so$level1 %in% "CD4"]

# Create dimreduc for latent space and use this for the transfer labels 
# Already in seurat object from previous integration
latent <- read.csv("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/IGT97_104/Seurat_PostQC_AllLatent/latent.csv", row.names = 1)
#rownames(latent) <- latent$X
#latent$X <- NULL
# columns must be labelled with key_ as a prefix
latent_matrix <- as.matrix(latent)
rownames(latent_matrix) <- colnames(CD4)
colnames(latent_matrix) <- c("latentSCVI_1", "latentSCVI_2", "latentSCVI_3", "latentSCVI_4", "latentSCVI_5", "latentSCVI_6", "latentSCVI_7", 
                                    "latentSCVI_8", "latentSCVI_9", "latentSCVI_10", "latentSCVI_11", "latentSCVI_12", "latentSCVI_13", "latentSCVI_14", 
                                    "latentSCVI_15", "latentSCVI_16", "latentSCVI_17", "latentSCVI_18", "latentSCVI_19", "latentSCVI_20", "latentSCVI_21", 
                                    "latentSCVI_22", "latentSCVI_23", "latentSCVI_24", "latentSCVI_25", "latentSCVI_26", "latentSCVI_27", "latentSCVI_28", 
                                    "latentSCVI_29", "latentSCVI_30")
latent_dimreduc <- CreateDimReducObject(latent_matrix, key = "latentSCVI_", assay = "RNA")
CD4@reductions$latent_TOTALVI <- latent_dimreduc

# Complete for level 2, then use the cluster memberships to transfer the group information
CD4$level2_new_latent_all_20250226 = as.character(CD4$annotation_level2_20250226)
# CD4$level2_new[CD4$annotation_level2_parent1 == "proliferating"] = "CD8AB_prolif"
#CD4$level2_new[CD4$annotation_level2_parent1 == "miniverse"] = "CD8AB_miniverse"
# CD4$level2_new[CD4$annotation_level2_parent1 == "preT"] = "CD8AB_preT"
#table(CD4$annotation_level2_20250226, CD4$ref_query)
#table(CD4$level2_new, CD4$annotation_level2_parent1)

CD4@meta.data$is_ref <- "other"
CD4$is_ref[!(colnames(CD4) %in% rownames(annotations_20250226))] = "query"
CD4$is_ref[rownames(annotations_20250226)] = "ref"
#CD4$is_ref[grepl("miniverse", CD4$cluster_query)] = "other" #"proliferating", "preT"
#CD4$is_ref[grepl("prolif", CD4$cluster_query)] = "proliferating" #"proliferating"
cell_embeddings_ref_latent = CD4[["latent_TOTALVI"]]@cell.embeddings[CD4$is_ref == "ref",]
cell_embeddings_query_latent = CD4[["latent_TOTALVI"]]@cell.embeddings[CD4$is_ref == "query",]

# Use the Annoy algorithm to find nearest neighbors between query (3k) and ref (10K) datasets
n_neighbors = 30  # Number of nearest neighbors to find k =30

# Create Annoy index for ref dataset
annoy_index = new(AnnoyAngular, ncol(cell_embeddings_ref_latent)) ##use cosine distance Angular
for (i in 1:nrow(cell_embeddings_ref_latent)) {
  annoy_index$addItem(i - 1, cell_embeddings_ref_latent[i, ])
}
annoy_index$build(10)  # Build the index with 10 trees

# Find nearest ref neighbors for each cell in query dataset
nn_indices <- t(sapply(1:nrow(cell_embeddings_query_latent), function(i) {
  annoy_index$getNNsByVector(cell_embeddings_query_latent[i, ], n_neighbors)
}))

# nn_indices gives you the indices of nearest neighbors in the ref dataset
# the rows are cells from query dataset, columns are the 30 nearest cells in the ref dataset
dim(nn_indices)
head(nn_indices)

labels_ref = as.character(CD4$annotation_level2_20250226[CD4$is_ref == "ref"])
# Transfer labels based on majority vote from nearest neighbors
transfer_labels <- apply(nn_indices, 1, function(neighbors) {
  # Get labels for the nearest neighbors
  neighbor_labels <- labels_ref[neighbors + 1]  # Add 1 for R's 1-based index
  
  # Return the most common label (majority vote)
  most_common_label <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  return(most_common_label)
})

# Now, transfer_labels contains the predicted labels for the 3k PBMC dataset
CD4$level2_new_latent_all_20250226[CD4$is_ref == "query"] = transfer_labels
head(transfer_labels)
table(transfer_labels)
#table(CD4$level2_predictions[CD4$is_ref == "query"])

#table(CD4$level2_predictions, CD4$level2_new_latent_all_20250226)
table(CD4$is_ref, CD4$level2_new_latent_all_20250226)
#table(CD4$annotation_level2_20241216, CD4$level2_new)

Cluster_annotation_table <- table(CD4$level2_new_latent_all_20250226[CD4$is_ref == "query"], CD4$IGTHT[CD4$is_ref == "query"])
write.table(Cluster_annotation_table, "Seurat_PostQC_AllLatent/Transfer_label_annotation_level2_20250226_table.txt", sep = ",", quote = F, col.names = T, row.names = T)
write.csv(Cluster_annotation_table, "Seurat_PostQC_AllLatent/Transfer_label_annotation_level2_20250226_table.csv", quote = F, col.names = T)

# Read in created so
CD4 <- readRDS("N:/CBDM_Lab/Odhran/scVerse/ImmgenT_freeze_20250109/IGT97_104/Seurat_PostQC/seuratobject_merged_RNA_withmetadata_withreductions_transferlabels.Rds")

## Add scvi mde embedding 
# Read in co-ords
mde_CD4 <- read.csv("Trial_TOTALVI_PostQC/mde_incremental.csv", row.names = 1)
mde_CD4 <- mde_CD4[colnames(CD4) ,]
mde_CD4_AllLatent <- read.csv("Trial_TOTALVI_PostQC_AllLatent/mde_incremental.csv", row.names = 1)
mde_CD4_AllLatent <- mde_CD4_AllLatent[colnames(CD4) ,]
View(mde_CD4)
umap_python <- read.csv("Trial_TOTALVI_PostQC/umap_python.csv", row.names = 1)
View(umap_python)
# Take only co-ordinates for cell embeddings
mde_coords <- mde_CD4[,1:2]
mde_coords_AllLatent <- mde_CD4_AllLatent[,1:2]
python_coords <- umap_python[,1:2]

# columns must be labelled with key_ as a prefix
# pyMDE
mde_coords_matrix <- as.matrix(mde_coords)
rownames(mde_coords_matrix) <- colnames(CD4)
colnames(mde_coords_matrix) <- c("mdeTOTALVI_1", "mdeTOTALVI_2")
mde_dimreduc <- CreateDimReducObject(mde_coords_matrix, key = "mdeTOTALVI_", assay = "RNA")
CD4@reductions$mde_incremental_TOTALVI <- mde_dimreduc

# pyMDEAllLatent
mde_coords_AllLatent_matrix <- as.matrix(mde_coords_AllLatent)
rownames(mde_coords_AllLatent_matrix) <- colnames(CD4)
colnames(mde_coords_AllLatent_matrix) <- c("mdeTOTALVIAllLatent_1", "mdeTOTALVIAllLatent_2")
mde_dimreduc_AllLatent <- CreateDimReducObject(mde_coords_AllLatent_matrix, key = "mdeTOTALVIAllLatent_", assay = "RNA")
CD4@reductions$mde_incremental_TOTALVI_AllLatent <- mde_dimreduc_AllLatent

# python - not necessary for all data so
python_coords_matrix <- as.matrix(python_coords)
rownames(python_coords_matrix) <- colnames(CD4)
colnames(python_coords_matrix) <- c("pythonTOTALVICD8_1", "pythonTOTALVICD8_2")
python_dimreduc <- CreateDimReducObject(python_coords_matrix, key = "pythonTOTALVICD8_", assay = "RNA")
CD4@reductions$python_TOTALVI_CD8 <- python_dimreduc



# Check clusters and re-annotate if necessary
# FeaturePlot(CD4, reduction = "mde_incremental",slot = "count", features = "Stmn1", order = T, raster = T, raster.dpi = c(1024,1024))

prefix_latent = "latent_TOTALVI"
n_latent = ncol(CD4[[prefix_latent]]@cell.embeddings)
CD4 = FindNeighbors(CD4, dims = 1:n_latent, reduction = prefix_latent)
CD4 = FindClusters(CD4, reCD4lution = 1, cluster.name = sprintf("Cluster_%s_Res1",prefix_latent))

# Colour palette
mypal_level2 <- setNames(mypal, sort(unique(CD4$level2)))
names(mypal_level2)[50] <- "query"
mypal_level2["query"] <- "gray"


# pyMDE
label_transfer <- DimPlot(CD4, reduction = "mde_incremental_TOTALVI", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("transfer labels from annotation_level2_20250226 using transfer RcppAnnoy (SCVI latent space, all)") +
  xlim(-2.5, 6) + ylim(-2.5,6)
#IGT97_104_cluster_query <- DimPlot(CD4, reduction = "mde_incremental_TOTALVI", group.by = "cluster_query", raster = F) + scale_color_manual(values = mypal_level2)+
# xlim(-2.5, 6) + ylim(-2.5,6)
IGT97_104_only_query <- DimPlot(CD4[, CD4$is_ref %in% "query"], reduction = "mde_incremental_TOTALVI", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("Transfer label result - only query cells") +
  xlim(-2.5, 6) + ylim(-2.5,6)
IGT1_96_only_query <- DimPlot(CD4[, CD4$is_ref %in% "ref"], reduction = "mde_incremental_TOTALVI", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("Original level 2 labels - IGT1_96") +
  xlim(-2.5, 6) + ylim(-2.5,6)

png("Seurat_PostQC/Transfer_label_latent_space_TOTALVI_20250226.png", width = 900, height = 900)
label_transfer
dev.off()

png("Seurat_PostQC/IGT97_104_only_cluster_mde_20250226.png", width = 1000, height = 900)
IGT97_104_only_query
dev.off()

#png("Seurat_PostQC/IGT97_104_cluster_query.png", width = 1000, height = 900)
#IGT97_104_cluster_query
#dev.off()

png("Seurat_PostQC/IGT1_96_cluster_query_20250226.png", width = 1000, height = 900)
IGT1_96_only_query
dev.off()

pdf("Seurat_PostQC/Transfer_label_David's_code_query_TOTALVI_20250226.pdf", 20, 20, useDingbats = F)
DimPlot(CD4, reduction = "mde_incremental_TOTALVI", group.by = "is_ref", raster = F, raster.dpi = c(512,512), label = T) +
  xlim(-2.5, 6) + ylim(-2.5,6)
#DimPlot(CD4, reduction = "mde_incremental_TOTALVI", group.by = sprintf("Cluster_%s_Res1",prefix_latent), raster = T, raster.dpi = c(512,512), label = T) + scale_color_manual(values = mypal)
IGT1_96_only_query
IGT97_104_only_query
label_transfer
dev.off()

# pyMDE All Latent
label_transfer_AllLatent <- DimPlot(CD4, reduction = "mde_incremental_TOTALVI_AllLatent", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("transfer labels from annotation_level2_20250226 using transfer RcppAnnoy (SCVI latent space, all)") +
  xlim(-2.5, 6) + ylim(-2.5,6)
#IGT97_104_cluster_query_AllLatent <- DimPlot(CD4, reduction = "mde_incremental_TOTALVI_AllLatent", group.by = "cluster_query", raster = F) + scale_color_manual(values = mypal_level2) +
# xlim(-2.5, 6) + ylim(-2.5,6)
IGT97_104_only_query_AllLatent <- DimPlot(CD4[, CD4$is_ref %in% "query"], reduction = "mde_incremental_TOTALVI_AllLatent", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("Transfer label result - only query cells") +
  xlim(-2.5, 6) + ylim(-2.5,6)
IGT1_96_only_query_AllLatent <- DimPlot(CD4[, CD4$is_ref %in% "ref"], reduction = "mde_incremental_TOTALVI_AllLatent", group.by = "level2", raster = T, raster.dpi = c(512,512)) + scale_color_manual(values = mypal_level2) + ggtitle("Original level 2 labels - IGT1_96") +
  xlim(-2.5, 6) + ylim(-2.5,6)

png("Seurat_PostQC_AllLatent/Transfer_label_latent_space_TOTALVI_20250226.png", width = 900, height = 900)
label_transfer_AllLatent
dev.off()

png("Seurat_PostQC_AllLatent/IGT97_104_only_cluster_mde_20250226.png", width = 1000, height = 900)
IGT97_104_only_query_AllLatent
dev.off()

#png("Seurat_PostQC_AllLatent/IGT97_104_cluster_query.png", width = 1000, height = 900)
#IGT97_104_cluster_query_AllLatent
#dev.off()

png("Seurat_PostQC_AllLatent/IGT1_96_cluster_query_20250226.png", width = 1000, height = 900)
IGT1_96_only_query_AllLatent
dev.off()

pdf("Seurat_PostQC_AllLatent/Transfer_label_David's_code_query_TOTALVI_20250226.pdf", 20, 20, useDingbats = F)
DimPlot(CD4, reduction = "mde_incremental_TOTALVI_AllLatent", group.by = "is_ref", raster = F, raster.dpi = c(512,512), label = T) +
  xlim(-2.5, 6) + ylim(-2.5,6)
#DimPlot(CD4, reduction = "mde_incremental_TOTALVI_AllLatent", group.by = sprintf("Cluster_%s_Res1",prefix_latent), raster = T, raster.dpi = c(512,512), label = T) + scale_color_manual(values = mypal)
IGT1_96_only_query_AllLatent
IGT97_104_only_query_AllLatent
label_transfer_AllLatent
dev.off()

saveRDS(CD4, "Seurat_PostQC_AllLatent/seuratobject_merged_RNA_withmetadata_withreductions_transferlabels.Rds")


## Annotate level2.group by cluster memberships
table(CD4$annotation_level2, CD4$annotation_level2_group_20250226)
resting_clusters <- unique(CD4$annotation_level2_20250226[CD4$annotation_level2_group_20250226 %in% "resting"])
activated_clusters <- unique(CD4$annotation_level2_20250226[CD4$annotation_level2_group_20250226 %in% "activated"])
CD4$level2_group_transfer_20250226 <- CD4$annotation_level2_group_20250226
CD4$level2_group_transfer_20250226[CD4$is_ref %in% "query" & CD4$level2_new_latent_all_20250226 %in% resting_clusters] <- "resting"
CD4$level2_group_transfer_20250226[CD4$is_ref %in% "query" & CD4$level2_new_latent_all_20250226 %in% activated_clusters] <- "activated"
CD4$level2_group_transfer_20250226[CD4$is_ref %in% "query" & CD4$level2_new_latent_all_20250226 %in% "CD4_miniverse"] <- "miniverse"
CD4$level2_group_transfer_20250226[CD4$is_ref %in% "query" & CD4$level2_new_latent_all_20250226 %in% "CD4_prolif"] <- "proliferating"
table(CD4$level2_new_latent_all_20250226, CD4$level2_group_transfer_20250226)

## Save transfer label annotation
annotation_level2_with_IGT97_104 <- as.data.frame(CD4$level2_new_latent_all_20250226)
colnames(annotation_level2_with_IGT97_104) <- "annotation_level2_20250226"
annotation_level2_with_IGT97_104$annotation_level2_group_20250226 <- CD4$level2_group_transfer_20250226
write.csv(annotation_level2_with_IGT97_104, "annotation_level2_20250226_IGT1_104_CD4.csv", quote = F, col.names = T)

# Save mde co-ordinates
mde_incremental_TOTALVI_coords <- as.data.frame(CD4@reductions$mde_incremental_TOTALVI_AllLatent@cell.embeddings)
colnames(mde_incremental_TOTALVI_coords) <- paste0("level2_", colnames(mde_incremental_TOTALVI_coords))
write.csv(mde_incremental_TOTALVI_coords, "mde_incremental_TOTALVI_coords_IGT1_104_CD4.csv", quote = F, col.names = T)


#CD4$level2_new[CD4$Cluster_totalvi_20241008_rmIGTsample_Res1 %in% "25" & CD4[["mde_incremental"]]@cell.embeddings[,1] < 1 ] = "CD8AB_cl25"
#CD4$level2_new[CD4$Cluster_totalvi_20241008_rmIGTsample_Res1 %in% "15" & CD4[["mde_incremental"]]@cell.embeddings[,1] > 0 ] = "CD8AB_cl26"
#CD4$level2_new[CD4$Cluster_totalvi_20241008_rmIGTsample_Res1 %in% "13" ] = "CD8AB_cl27"
#CD4$level2_new[CD4$Cluster_totalvi_20241008_rmIGTsample_Res1 %in% "20" ] = "CD8AB_cl28"
#CD4$level2_new[CD4$level2_new  %in% c("CD8AB_cl19")] = "CD8AB_prolif"

