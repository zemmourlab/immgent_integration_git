# David Zemmour
# R
# usage: Rscript run_totalvi_plot.R /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration [prefix for output]

options(max.print=1000)
options(Seurat.object.assay.version = 'v5')

args = commandArgs(TRUE)
path_to_wd = args[1]
prefix = args[2] 

#path_to_wd = "/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration"
#seurat_object = "spleen_standards.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"
#prefix = "test1"

setwd(path_to_wd)

message("loading R libraries and custom functions")
libs = c("Seurat", "ggplot2", "dplyr", "ggrastr", "RColorBrewer", "pals", "scales") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal = mypal[-4]
parade = function(n) { return(Seurat::DiscretePalette(n, palette = "parade", shuffle = F)) }
length(glasbey())
length(polychrome())
mypal2 = c(glasbey(), polychrome(), mypal)
names(mypal2) = NULL

MyPlots = function (seurat_object = so, dim1 = so[["umap_unintegrated"]]@cell.embeddings[,1], dim2 = so[["umap_unintegrated"]]@cell.embeddings[,2], color_by = "spleen_standard", split_by1 = "IGT", split_by2 =  NULL, genes = c("Foxp3", "Il2ra"), cluster_key = "ClusterSCVI_Res", mypal = glasbey()) {
    
    so = seurat_object
    #so@meta.data[,"split_by"] = so@meta.data[,split_by]
    so@meta.data[,"color_by"] = factor(so@meta.data[,color_by])
    so@meta.data[,"split_by1"] = so@meta.data[,split_by1]
    so@meta.data[,"split_by2"] = so@meta.data[,split_by2]
    
    alpha = 0.5
    #sample_name_colors = color_palette[1:length(unique(so@meta.data[,color_by]))]
    #names(sample_name_colors) = levels(so$sample_name)
    #sample_name_colors2 = sample_name_colors
    #sample_name_colors2[grepl("WT", names(sample_name_colors))] = "grey"
    
    message("Plot 1: UMAP")
    plot1 = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) + 
        geom_point_rast(aes(dim1, dim2, color = color_by), alpha = I(alpha), raster.dpi = 100) +
        theme_bw() + scale_color_manual(values = mypal)
    print(plot1)
    
    message("Plot 2: UMAP split")
    tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)
    bkgrd = data.frame(dim1 = dim1, dim2 = dim2)
    
    p = ggplot(bkgrd) + geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, alpha = 0.2, raster.dpi = 50)
    p2 = geom_point_rast(data = tmp, aes(dim1, dim2, color = color_by), size = 1,  alpha = alpha) 
    
    if (is.null(split_by2)) {
        plot2 = p + p2 + scale_color_manual(values = mypal) + theme_bw() + facet_wrap(facets = vars(so@meta.data[,"split_by1"])) + ggtitle(label = sprintf("color: %s, grid: %s", color_by, split_by1))
    } else {
        plot2 = p + p2 + scale_color_manual(values = mypal)  + theme_bw() + facet_grid(rows = vars(so@meta.data[,"split_by1"]), cols = vars(so@meta.data[,"split_by2"])) + ggtitle(label = sprintf("color: %s, grid: %s x %s", color_by, split_by1, split_by2))
    }
    print(plot2)
    
    message("Plot 3: UMAP genes")
    
    for (g in genes) {
        print(g)
        tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2, size = so@assays$RNA$counts[rownames(so@assays$RNA$counts) %in% g,])
        bkgrd = data.frame(dim1 = dim1, dim2 = dim2)
        
        p = ggplot(bkgrd) + geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, raster.dpi = 50)
        p2 = geom_point(data = tmp, aes(dim1, dim2, color = size > 0, size = size,  alpha = size > 0)) 
        #p2 =  geom_point(data = tmp2, aes(dim1, dim2, color = size > 0, size = size), alpha = I(alpha))  + scale_color_manual
        
        
        if (is.null(split_by2)) {
            plot3 = p + p2  + 
                scale_color_manual(values = c("black", "red")) + 
                scale_alpha_manual(values = c(0.5,1))  + 
                theme_bw()  + 
                facet_wrap(facets = vars(so@meta.data[,"split_by1"])) + 
                ggtitle(label = sprintf("gene: %s, color: %s, grid: %s", g, color_by, split_by1))
        } else {
            plot3 = p + p2  + 
                scale_color_manual(values = c("black", "red")) + 
                scale_alpha_manual(values = c(0.5,1))  + 
                theme_bw() + 
                facet_grid(rows = vars(so@meta.data[,"split_by1"]), cols = vars(so@meta.data[,"split_by2"])) +
                ggtitle(label = sprintf("gene: %s, color: %s, grid: %s", g, color_by, split_by1))
        }
        print(plot3)
        
    }
    
}

message("loading seurat object")
so = readRDS(file = sprintf("%s_so.Rds", prefix))

message("Plotting")

pdf(sprintf("%s_UMAP.pdf", prefix), 20, 20, useDingbats = F)
MyPlots(seurat_object = so, 
        dim1 = so[["umap_totalvi"]]@cell.embeddings[,1], 
        dim2 = so[["umap_totalvi"]]@cell.embeddings[,2], 
        color_by = "IGT", 
        split_by1 = "is_spleen_standard",
        split_by2 =  NULL, 
        genes = NULL, 
        cluster_key = NULL, mypal = mypal2) #glasbey()

MyPlots(seurat_object = so, 
        dim1 = so[["umap_totalvi"]]@cell.embeddings[,1], 
        dim2 = so[["umap_totalvi"]]@cell.embeddings[,2], 
        color_by = "organ", 
        split_by1 = "is_spleen_standard",
        split_by2 =  NULL, 
        genes = NULL, 
        cluster_key = NULL, mypal = mypal2) #glasbey()
dev.off()

igt = as.character(unique(so$IGT))
x = seq(1, length(unique(so$IGT)), by = 5)
for (n in 1:length(x)) {
    i = x[n]
    j = ifelse(is.na(x[n+1]-1), yes = length(igt), no = x[n+1]-1)
    igt_to_plot = igt[i:j]
    print(sprintf("Plotting %s_Plots_%s_%s.pdf", prefix, igt[i],igt[j]))
    pdf(sprintf("%s_Plots_%s_%s.pdf", prefix, igt[i],igt[j]), 20, 30, useDingbats = F)
    MyPlots(seurat_object = so[, so$IGT %in% igt_to_plot], 
            dim1 = so[, so$IGT %in% igt_to_plot][["umap_totalvi"]]@cell.embeddings[,1], 
            dim2 = so[, so$IGT %in% igt_to_plot][["umap_totalvi"]]@cell.embeddings[,2], 
            color_by = "IGT", 
            split_by1 = "IGT",
            split_by2 =  "is_spleen_standard", 
            genes = c("Ms4a1", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Trgc1", "Itgax", "Itgam"), 
            cluster_key = NULL, mypal = mypal) #glasbey()
    dev.off()
}


#x = seq(1, length(unique(so$IGT)), by = 5)
#for (n in 1:length(x)) {
#    i = x[n]
#    j = x[n+1]-1
#    print(sprintf("Plotting %s_Plots_IGT%s_%s.pdf", prefix, i,j))
#    pdf(sprintf("%s_Plots_IGT%s_%s.pdf", prefix, i,j), 20, 30, useDingbats = F)
    
    #pdf(sprintf("%s_Plots.pdf", prefix), 20, 100, useDingbats = F)
#    MyPlots(seurat_object = so[, so$IGT %in% sprintf("IGT%s", i:j)], 
#            dim1 = so[, so$IGT %in% sprintf("IGT%s", i:j)][["umap_totalvi"]]@cell.embeddings[,1], 
#            dim2 = so[, so$IGT %in% sprintf("IGT%s", i:j)][["umap_totalvi"]]@cell.embeddings[,2], 
#            color_by = "IGT", 
#            split_by1 = "IGT",
#            split_by2 =  "is_spleen_standard", 
#            genes = c("Ms4a1", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Trgc1", "Itgax", "Itgam"), 
#            cluster_key = NULL, mypal = mypal) #glasbey()
#    dev.off()
#}




