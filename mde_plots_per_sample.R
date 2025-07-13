# David Zemmour
# R
# usage: Rscript mde_plots_per_sample.R [path_to_seurat_object] [output_dir]

# Parse arguments
args = commandArgs(TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript mde_plots_per_sample.R [path_to_seurat_object] [output_dir]")
}
path_to_seurat_object = args[1]
output_dir = args[2]

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
libs = c("gridExtra","Seurat", "ggplot2", "gridExtra", "scattermore", "dplyr", "reshape2", "RColorBrewer", "pals", "scales") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))
source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")


message("loading seurat object")
so_orig = readRDS(file = path_to_seurat_object)

message("Set color palette")
library("pals")
# n = 70
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# mypal1 = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
# mypal1 = mypal1[-4]
# mypal = c(glasbey(), polychrome(), mypal1)
# names(mypal) = NULL
# mypal_organ = setNames(mypal, unique(so_orig$organ_simplified))
# mypal_organ = mypal_organ[!is.na(names(mypal_organ))]
# mypal_level1 = setNames(mypal, unique(so_orig$annotation_level1))
# mypal_level1 = mypal_level1[!is.na(names(mypal_level1))]
# mypal_level2 = setNames(mypal, unique(so_orig$annotation_level2))
# mypal_level2 = mypal_level2[!is.na(names(mypal_level2))]
# mypal_level2["CD4_cl18"] = "blue"
# mypal_igt = setNames(mypal, unique(so_orig$IGT))
# mypal_igt = mypal_igt[!is.na(names(mypal_igt))]

mypal_list = readRDS("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/immgent_color_palettes/20250709/immgent_color_palettes_list_20250709.Rds")
#mypal_list = readRDS("/Users/david/Library/CloudStorage/GoogleDrive-david.puyraimond@gmail.com/My Drive/ImmgenT/analysis/immgent_integration_git/immgent_color_palettes/20250709/immgent_color_palettes_list_20250709.Rds")
mypal_level1 = mypal_list[["level1"]]
mypal_level2_all = mypal_list[["level2_all"]]
mypal_level2 = c(mypal_list[["CD4"]], mypal_list[["CD8"]], mypal_list[["Treg"]], mypal_list[["gdT"]], mypal_list[["nonconv"]], mypal_list[["CD8aa"]], mypal_list[["DN"]], mypal_list[["DP"]],mypal_list[["thymocyte"]])
mypal_organ = mypal_list[["organ_simplified"]]


message("Plotting")

so_orig$is_agspe = !is.na(so_orig$Ag_spe)

so = list()
so[["CD4"]] = so_orig[,so_orig$annotation_level1 == "CD4"]
so[["CD8"]] =  so_orig[,so_orig$annotation_level1 == "CD8"]
so[["Treg"]] =  so_orig[,so_orig$annotation_level1 == "Treg"]
so[["gdT"]] = so_orig[,so_orig$annotation_level1 == "gdT"]
so[["CD8aa"]] = so_orig[,so_orig$annotation_level1 == "CD8aa"]
so[["nonconv"]] = so_orig[,so_orig$annotation_level1 == "nonconv"]
so[["DN"]] = so_orig[,so_orig$annotation_level1 == "DN"]
so[["DP"]] = so_orig[,so_orig$annotation_level1 == "DP"]

# i="IGT37"
level1 = c("CD4", "CD8", "Treg", "gdT", "CD8aa", "nonconv", "DN", "DP")

message("Plots for Ag specific cells only")
igts = so_orig@meta.data %>% filter(is_agspe == T) %>% pull(IGT) %>% unique()

pdf(sprintf("%s/IGT1-96_MDEs_AgSpe.pdf", output_dir), width = 20, height = 20, useDingbats = F)
for (i in igts) {
    samples = so_orig@meta.data %>% filter(IGT == i) %>% pull(sample_code) %>% unique()
    for (s in samples) {
        print(s)
        ps2 = list()
        ps3 = list()
        for (l in level1) {
            # print(l)
            p = MyDimPlotHighlight(seurat_object = so[[l]], umap_to_plot = "mde_incremental", cells_to_highlight = colnames(so[[l]])[which(so[[l]]$sample_code == s)], highlight_column_name = "is_agspe", pixels = c(512, 512), mycols = c("black", "red"), title = sprintf("%s in %s", l, s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
            ps2[[l]] = p$plot2
            ps3[[l]] = p$plot3 #to get plot with labels
        }
        grid.arrange(grobs = ps3, ncol = 3, nrow = 3)
        grid.arrange(grobs = ps2, ncol = 3, nrow = 3)
        # dev.off()
    }
}
dev.off()

pdf(sprintf("%s/IGT1-96_MDEs_AgSpe2.pdf", output_dir), width = 20, height = 20, useDingbats = F)
ps2 = list()
ps3 = list()
for (l in level1) {
    print(l)
    p = MyDimPlotHighlight(seurat_object = so[[l]], umap_to_plot = "mde_incremental", cells_to_highlight = colnames(so[[l]])[which(so[[l]]$is_agspe == T)], highlight_column_name = "Ag_spe", pixels = c(512, 512), mycols = mypal, title = sprintf("Ag spe in %s", l), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
    ps2[[l]] = p$plot2
    ps3[[l]] = p$plot3 #to get plot with labels
}
grid.arrange(grobs = ps3, ncol = 3, nrow = 3)
grid.arrange(grobs = ps2, ncol = 3, nrow = 3)
dev.off()


message("One file per sample, with labels only")
for (i in unique(so_orig$IGT)) {
    samples = so_orig@meta.data %>% filter(IGT == i) %>% pull(sample_code) %>% unique()
    for (s in samples) {
        print(s)
        ensure_directory(sprintf("%s/tmp/", output_dir))
        pdf(sprintf("%s/tmp/%s_MDEsWithLabels.pdf", output_dir, s), width = 20, height = 20, useDingbats = F)
        message(sprintf("PDF IN %s/tmp/%s_MDEsWithLabels.pdf", output_dir, s))
        ps2 = list()
        ps3 = list()
        p = MyDimPlotHighlight(seurat_object = so_orig, umap_to_plot = "mde2_totalvi_20241006", cells_to_highlight = colnames(so_orig)[which(so_orig$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s",s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
        # ps2[["all"]] =  p$plot2
        ps3[["all"]] =  p$plot3
        for (l in level1) {
            # print(l)
            p = MyDimPlotHighlight(seurat_object = so[[l]], umap_to_plot = "mde_incremental", cells_to_highlight = colnames(so[[l]])[which(so[[l]]$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s in %s", l, s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
            # ps2[[l]] = p$plot2
            ps3[[l]] = p$plot3 #to get plot with labels
        }
        grid.arrange(grobs = ps3, ncol = 3, nrow = 3)
        # grid.arrange(grobs = ps2, ncol = 3, nrow = 3)
        dev.off()
    }
}

message("One file per sample in separate IGT folders, with and without labels")
for (i in unique(so_orig$IGT)) {
    samples = so_orig@meta.data %>% filter(IGT == i) %>% pull(sample_code) %>% unique()
    for (s in samples) {
        print(s)
        ensure_directory(sprintf("%s/%s/",output_dir, i))
        pdf(sprintf("%s/%s/%s_MDEs.pdf", output_dir, i, s), width = 20, height = 20, useDingbats = F)
        message(sprintf("PDF IN %s/%s/%s_MDEs.pdf", output_dir, i, s))
        ps2 = list()
        ps3 = list()
        p = MyDimPlotHighlight(seurat_object = so_orig, umap_to_plot = "mde2_totalvi_20241006", cells_to_highlight = colnames(so_orig)[which(so_orig$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s",s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
        ps2[["all"]] =  p$plot2
        ps3[["all"]] =  p$plot3
        for (l in level1) {
            # print(l)
            p = MyDimPlotHighlight(seurat_object = so[[l]], umap_to_plot = "mde_incremental", cells_to_highlight = colnames(so[[l]])[which(so[[l]]$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s in %s", l, s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
            ps2[[l]] = p$plot2
            ps3[[l]] = p$plot3 #to get plot with labels
        }
        grid.arrange(grobs = ps3, ncol = 3, nrow = 3)
        grid.arrange(grobs = ps2, ncol = 3, nrow = 3)
        dev.off()
    }
}

# message("One file with all plots")
# pdf(sprintf("%s/IGT1-96_MDEs",output_dir), width = 20, height = 20, useDingbats = F)
# for (i in unique(so_orig$IGT)) {
#     samples = so_orig@meta.data %>% filter(IGT == i) %>% pull(sample_code) %>% unique()
#     for (s in samples) {
#         print(s)
#         ensure_directory(sprintf("MDE_InEachSample/%s/", i))
#         # pdf(sprintf("MDE_InEachSample/%s/%s_MDEs.pdf", i, s), width = 20, height = 20, useDingbats = F)
#         ps2 = list()
#         ps3 = list()
#         p = MyDimPlotHighlight(seurat_object = so_orig, umap_to_plot = "mde2_totalvi_20241006", cells_to_highlight = colnames(so_orig)[which(so_orig$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s",s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
#         ps2[["all"]] =  p$plot2
#         ps3[["all"]] =  p$plot3
#         for (l in level1) {
#             # print(l)
#             p = MyDimPlotHighlight(seurat_object = so[[l]], umap_to_plot = "mde_incremental", cells_to_highlight = colnames(so[[l]])[which(so[[l]]$sample_code == s)], highlight_column_name = "annotation_level2", pixels = c(512, 512), mycols = mypal_level2, title = sprintf("%s in %s", l, s), highlight_size = 1, highlight_alpha = 1, print_plot1 = T, print_plot2 = T, labelclusters = T)
#             ps2[[l]] = p$plot2
#             ps3[[l]] = p$plot3 #to get plot with labels
#         }
#         grid.arrange(grobs = ps3, ncol = 3, nrow = 3)
#         # grid.arrange(grobs = ps2, ncol = 3, nrow = 3)
#         # dev.off()
#     }
# }
# dev.off()





