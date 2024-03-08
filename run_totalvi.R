# David Zemmour
# R
# usage: Rscript run_totalvi.R /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration [seurat Rds] [prefix for output]

options(max.print=1000)
options(Seurat.object.assay.version = 'v5')

args = commandArgs(TRUE)
path_to_wd = args[1] 
seurat_object = args[2] 
prefix = args[3] 

#path_to_wd = "/scratch/midway3/zemmour/integration"
#seurat_object = "immgent_20231030.Rds" 
#prefix = "totalvi_igt1_56_20231030_allgenes"

setwd(path_to_wd)

message("loading R libraries and custom functions")
libs = c("Seurat", "ggplot2", "dplyr", "ggrastr", "scales", "reticulate", "sceasy", "SingleCellExperiment", "anndata") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

#My functions
ConvertS5toS3 = function(so, assay1 = "RNA", assay2 = "ADT") {
	    so.v3 = CreateAssayObject(counts = so[[assay1]]$counts)
    so.v3 = CreateSeuratObject(so.v3)
        so.v3[[assay2]] = CreateAssayObject(counts = so[[assay2]]$counts) #samples@assays$ADT$counts
        message("same colnames between 2 assays")
	    print(all(colnames(so.v3[[assay1]]$counts) == colnames(so.v3[[assay2]]$counts)))
	    message("same colnames between S3 and S5")
	        print(all(colnames(so.v3) == colnames(so)))
	        so.v3@meta.data = so@meta.data
		    return(so.v3)
}

message("Loading python library")
sc = tryCatch(import("scanpy", convert = FALSE), error = function(e) {message("scanpy error usually not a problem undefined symbol: cblas_cdotc_sub")})
scvi = import("scvi", convert = FALSE)
sys = import ("sys", convert = FALSE)
np = import("numpy", convert=FALSE)

message("Loading seurat object")
so = readRDS(seurat_object)
#prefix = "totalvi_igt1_48_20230706_allgenes"

so$colnames = colnames(so)

message("Check that ADT counts are numerical. 0 is cells when no ADT data, no NA. Range:")
print(range(so[["ADT"]]$counts))

message("Remove IGT with only 1 cells")
rm_igt = so@meta.data %>% group_by(IGT) %>% summarize(counts = n()) %>% filter(counts == 1) %>% pull(IGT)
print(rm_igt)
so = so[,!so$IGT %in% rm_igt]

message("Finding variable features across datasets (not necessary when integrating with all genes)")
#Find Variable features: not useful if running totalvi on all genes
#scvi-tools models require the raw counts; normalization just here becayse often used for other purposes
#Another important thing to keep in mind is highly-variable gene selection. While scVI and scANVI both accomodate for large gene sets in terms of runtime, we usually recommend filtering genes for best performance when the dataset has few number of cells. As a rule of thumb, performance starts to decrease when number of cells and number of genes are comparable. This point is emphasized in this comparative analysis of data integration algorithms for scRNA-seq data.
DefaultAssay(so) = "RNA"
so[["RNA"]] = split(so[["RNA"]], f = so$IGT)
so = NormalizeData(so)
so = FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
plot1 = VariableFeaturePlot(so)
plot2 = LabelPoints(plot = plot1, points = VariableFeatures(so), repel = TRUE)

pdf(file = sprintf("%s_HVG.pdf", prefix), 20, 20, useDingbats = F)
plot2
dev.off()

so = JoinLayers(so)

message("Converting Seurat V5 to V3")
so.v3 = ConvertS5toS3(so) #so[VariableFeatures(so)]is using totalvi with variable genes only
sce = Seurat::as.SingleCellExperiment(so.v3)

message("Converting to AnnData")
adata = convertFormat(sce, from = "sce", to= "anndata", main_layer = "counts", drop_single_values=FALSE)
adata_protein = convertFormat(altExp(sce), from="sce", to="anndata", main_layer="counts", drop_single_values=FALSE)
adata$obsm[["protein"]] = adata_protein$to_df()
print(adata) # Note generally in Python, dataset conventions are obs x var

# run setup_anndata: https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.setup_anndata

message("Setting up AnnData")
scvi$model$TOTALVI$setup_anndata(adata, protein_expression_obsm_key="protein", categorical_covariate_keys = list("IGT"))
message("SAYS NO GPU EVEN IF THERE IS!!")

message("Setting up the model")
n_latent = as.integer(30)
#model = scvi$model$TOTALVI(adata, n_latent = n_latent)
model = scvi$model$TOTALVI(adata, n_latent = n_latent, gene_likelihood = "nb")

print(model)

message("Training the model")
#max_epochs = min(round((20000 / as.numeric(as.character(adata$n_obs))) * 400), 400)
model$train() #use_gpu = T
# to specify the number of epochs when training:
# model$train(max_epochs = as.integer(400))

message("Saving the model")
model$save(sprintf("%s/", prefix), save_anndata = T) 
#model = scvi$model$TOTALVI$load(sprintf("%s/", prefix))

message("Getting latent representation")
latent = model$get_latent_representation()
latent = as.matrix(latent)
rownames(latent) = colnames(so)
so[["totalvi"]] = CreateDimReducObject(embeddings = latent, key = "totalvi_", assay = DefaultAssay(so))

message("UMAP, clustering")
#so = RunUMAP(so, dims = 1:10, reduction = "scvi", n.components = 2)
so = RunUMAP(so, dims = 1:n_latent, reduction = "totalvi", n.components = 2, reduction.name = "umap_totalvi")

so = FindNeighbors(so, dims = 1:n_latent, reduction = "totalvi")
so = FindClusters(so, resolution = 0.25, cluster.name = "ClusterTOTALVI_Res0.25")
so = FindClusters(so, resolution = 0.5, cluster.name = "ClusterTOTALVI_Res0.5")
so = FindClusters(so, resolution = 1, cluster.name = "ClusterTOTALVI_Res1")
so = FindClusters(so, resolution = 2, cluster.name = "ClusterTOTALVI_Res2")
so = FindClusters(so, resolution = 3, cluster.name = "ClusterTOTALVI_Res3")
so = FindClusters(so, resolution = 4, cluster.name = "ClusterTOTALVI_Res4")

#so = JoinLayers(so)
message("Saving regular seurat object V5")
saveRDS(so, file = sprintf("%s_so.Rds", prefix))

message("Making BPCells object")
library(BPCells)
write_matrix_dir(mat = so[["RNA"]]$counts, dir = sprintf("bpcells/%s_counts", prefix))
counts.mat = open_matrix_dir(dir = sprintf("bpcells/%s_counts", prefix))

so_bp = so
so_bp[["RNA"]]$counts = counts.mat

message("Saving Seurat BPCells object")
saveRDS(
    object = so_bp,
    file = sprintf("%s_BPCells_so.Rds", prefix),
    destdir = sprintf("bpcells/%s_counts", prefix)
)


