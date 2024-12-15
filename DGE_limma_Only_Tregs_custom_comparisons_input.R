## o2 Setup for DGE Limma Trend Signature calculation ##

# Rscript DGE_limma_Only_Tregs_custom_comparisons_input.R [cwd] [path to R script with functions (David)] [path to so] [subgroup] 

options(max.print=1000)
args = commandArgs(TRUE)
#setwd(args[1])
source(args[2])

libs = c("limma", "edgeR", "Seurat", "ggplot2", "dplyr", "rlang","reshape2", "ggrastr", "RColorBrewer", "pals", "scales", "pheatmap") # BPCells
#libs = c("matrixStats", "gplots","ggplot2", "reshape2", "scales", "gridExtra", "dplyr", "RColorBrewer", "grid", "Rtsne", "limma",   "RColorBrewer", "pheatmap", "Seurat",  "ggrastr", "ggbeeswarm") #pwr, preprocessCore, "readxl", vegan", "genefilter" "scde" "cellrangerRkit", "Rsamtools", "GenomicRanges", "GenomicAlignments", "VGAM", "WGCNA",rafalib, ‘UsingR’,MAST "ggrastr",
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

#options(Seurat.object.assay.version = 'v5')

#Gradients
library(RColorBrewer)
ColorRamp = rev(colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100))
ColorRamp = rev(colorRampPalette(c('red','white','blue'))(20))
ColorRamp = rev(rainbow(10, end = 4/6))
library(viridis)
#ColorRamp = rev(viridis(100))
#ColorRamp = rev(cividis(100))
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

#parade = function(n) { return(Seurat::DiscretePalette(n, palette = "parade", shuffle = F)) }


#pdf("mypal.pdf", 10, 10)
#pie(rep(1,n), col=unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))) )
#pie(rep(1,n), col=unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))), labels =  unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
#dev.off()

PlotsAfterIntegration = function (seurat_object = so, dim1 = so@reductions$umap@cell.embeddings[,1], dim2 = so@reductions$umap@cell.embeddings[,2], split_by = "IGT", genes = c("Foxp3", "Il2ra"), cluster_key = "ClusterSCVI_Res") {
  
  so = seurat_object
  so@meta.data[,"split_by"] = so@meta.data[,split_by]
  
  alpha = 0.5
  sample_name_colors = mypal[1:length(unique(so@meta.data[,split_by]))]
  names(sample_name_colors) = levels(so$sample_name)
  sample_name_colors2 = sample_name_colors
  sample_name_colors2[grepl("WT", names(sample_name_colors))] = "grey"
  
  message("Plot 1: UMAP")
  p1 = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) + geom_point_rast(aes(dim1, dim2, color = split_by), alpha = I(alpha), raster.dpi = 100) + theme_bw() + scale_color_manual(values = mypal)
  print(p1)
  
  message("Plot 2: UMAP split")
  tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)
  bkgrd = data.frame(dim1 = dim1, dim2 = dim2)
  
  p = ggplot(bkgrd) + geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, alpha = 0.2, raster.dpi = 50)
  p2 = geom_point(data = tmp, aes(dim1, dim2, color = split_by), size = 1,  alpha = alpha) 
  
  q = p + p2 + scale_color_manual(values = sample_name_colors)  + theme_bw() + facet_wrap(~ split_by)
  print(q)
  
  message("Plot 3: UMAP genes")
  
  for (g in genes) {
    print(g)
    tmp = data.frame(so@meta.data, dim1 = dim1, dim2 = dim2, size = so@assays$RNA$counts[rownames(so@assays$RNA$counts) %in% g,])
    bkgrd = data.frame(dim1 = dim1, dim2 = dim2)
    
    p = ggplot(bkgrd) + geom_point_rast(aes(dim1, dim2), color = "grey", size = 0.1, raster.dpi = 50)
    p2 = geom_point(data = tmp, aes(dim1, dim2, color = size > 0, size = size,  alpha = size > 0)) 
    #p2 =  geom_point(data = tmp2, aes(dim1, dim2, color = size > 0, size = size), alpha = I(alpha))  + scale_color_manual
    
    p3 = p + p2  + scale_color_manual(values = c("black", "red")) + scale_alpha_manual(values = c(0.5,1))  + theme_bw() + facet_wrap(~split_by) + ggtitle(g)
    print(p3)
  }
  
  
  message("Plot 4: UMAP Clusters") 
  for (i in colnames(so@meta.data)[grepl("cluster_key", colnames(so@meta.data))]) {
    print(i)
    p = ggplot(data.frame(so@meta.data, dim1 = dim1, dim2 = dim2)) + 
      geom_point_rast(aes(dim1, dim2), color = so@meta.data[,i], alpha = I(alpha), raster.dpi = 100) + scale_color_manual(values = mypal) + ggtitle(i) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p + facet_wrap(~split_by))
  }
  
}

GetGroups <- function(metadata, group1_filter, group2_filter, id_column) {
  require(dplyr)
  require(rlang)
  # Evaluate the filter expressions within the context of metadata
  group1 <- metadata %>% filter(!! group1_filter)
  group2 <- metadata %>% filter(!! group2_filter)
  
  # Pull the specified identifier column
  group1_ids <- group1 %>% pull(!! sym(id_column))
  group2_ids <- group2 %>% pull(!! sym(id_column))
  
  # Combine groups and add group identifiers
  group1$group <- "group1"
  group2$group <- "group2"
  combined_groups <- rbind(group1, group2)
  
  # Function to perform Chi-squared test
  test_balanced <- function(var) {
    obs <- table(combined_groups$group, combined_groups[[var]])
    test <- tryCatch({
      chisq.test(obs)
    }, error = function(e) NA)  # Return NA on error
    if (!is.na(test$p.value)) {
      return(test$p.value)
    } else {
      return(1)  # Return non-significant p-value if test could not be performed
    }
  }
  
  # Apply test to each relevant column and suppress warnings
  results <- suppressWarnings(
    lapply(colnames(metadata)[!colnames(metadata) %in% c("sample_id", id_column)], test_balanced)
  )
  names(results) <- colnames(metadata)[!colnames(metadata) %in% c("sample_id", id_column)]
  
  # Check for significant results
  if (any(na.omit(unlist(results) < 0.05))) {
    message("WARNING: IMBALANCED GROUPS FOR:")
    print(names(results)[which(unlist(results) < 0.05)])
  }
  
  return(list(group1_ids = group1_ids, group2_ids = group2_ids, results = results))
}

CreateComparisonName <- function(group1_filter, group2_filter) {
  require(rlang)
  require(stringr)
  # Convert the expressions to strings
  group1_name <- expr_text(group1_filter)
  group2_name <- expr_text(group2_filter)
  
  # Function to extract quoted substrings
  extract_quoted_parts <- function(name) {
    str_match_all(name, '"([^"]+)"')[[1]][,2]
  }
  
  # Extract quoted parts for each group
  group1_parts <- extract_quoted_parts(group1_name)
  group2_parts <- extract_quoted_parts(group2_name)
  
  # Concatenate the extracted parts
  group1_clean <- paste0(group1_parts, collapse = "")
  group2_clean <- paste0(group2_parts, collapse = "")
  
  # Create the comparison name
  comparison_name <- paste0(group1_clean, "_vs_", group2_clean)
  
  return(comparison_name)
}

CollapseDiff_limmatrend = function(l) {
  #l = lapply(l, function(x) x = data.frame(x, rownames = rownames(x)))
  l = lapply(l, function(x) x = data.frame(x, rownames = l[[1]]$SYMBOL))

  fc = data.frame(l[[1]][,c("rownames", "logFC")])
  for (i in 2:length(l)) { fc = merge(x = fc, y = l[[i]][,c("rownames", "logFC")], by.x = "rownames", by.y = "rownames", all = T) }
  
  mean = data.frame(l[[1]][,c("rownames", "AveExpr")])
  for (i in 2:length(l)) { mean = merge(x = mean, y = l[[i]][,c("rownames", "AveExpr")], by.x = "rownames", by.y = "rownames", all = T) }
  
  pv = data.frame(l[[1]][,c("rownames", "P.Value")])
  for (i in 2:length(l)) { pv = merge(x = pv, y = l[[i]][,c("rownames", "P.Value")], by.x = "rownames", by.y = "rownames", all = T) }
  
  qcl = data.frame(l[[1]][,c("rownames", "adj.P.Val")])
  for (i in 2:length(l)) { qcl = merge(x = qcl, y = l[[i]][,c("rownames", "adj.P.Val")], by.x = "rownames", by.y = "rownames", all = T) }
  
  # OC Edit - 24-12-02
  pct.1 = data.frame(l[[1]][,c("rownames", "pct.1")])
  for (i in 2:length(l)) { pct.1 = merge(x = pct.1, y = l[[i]][,c("rownames", "pct.1")], by.x = "rownames", by.y = "rownames", all = T) }
  
  pct.2 = data.frame(l[[1]][,c("rownames", "pct.2")])
  for (i in 2:length(l)) { pct.2 = merge(x = pct.2, y = l[[i]][,c("rownames", "pct.2")], by.x = "rownames", by.y = "rownames", all = T) }
  
  rownames(fc) = fc[,1]
  rownames(mean) = fc[,1]
  rownames(pv) = pv[,1]
  rownames(qcl) = qcl[,1]
  rownames(pct.1) = pct.1[,1]
  rownames(pct.2) = pct.2[,1]
  fc = fc[,-1]
  mean = mean[,-1]
  pv = pv[,-1]
  qcl = qcl[,-1]
  pct.1 = pct.1[,-1]
  pct.2 = pct.2[,-1]
  colnames(fc) = paste("logFC", names(l), sep = "_")
  colnames(mean) = paste("AveExpr", names(l), sep = "_")
  colnames(pv) = paste("P.Value", names(l), sep = "_")
  colnames(qcl) = paste("adj.P.Val", names(l), sep = "_")
  colnames(pct.1) = paste("pct.1", names(l), sep = "_")
  colnames(pct.2) = paste("pct.2", names(l), sep = "_")
  all(rownames(fc) == rownames(pv))
  all(rownames(fc) == rownames(mean))
  all(rownames(fc) == rownames(qcl))
  all(rownames(fc) == rownames(pct.1))
  all(rownames(fc) == rownames(pct.2))
  fc[is.na(fc)] = 0
  pv[is.na(pv)] = 1
  qcl[is.na(qcl)] = 1
  pct.1[is.na(pct.1)] = 0
  pct.2[is.na(pct.2)] = 0
  diff = data.frame(fc, mean, pv, qcl, pct.1, pct.2)
  which(is.na(diff))
  
  return(diff)
  
}

#########
Vplot = function(vplot, xlab = "FC", xlimits = c(0.1, 10), ylimits = c(10^-300,1)) {
  p = ggplot(data = vplot) + geom_point(aes(x = fc, y = pval), colour = "black", alpha = I(1), size = I(0.5)) +
    scale_x_continuous(trans = log_trans(10), limits = xlimits) + #breaks = c(1/50,1/40, 1/30, 1/20, 1/10, 1/5,1/2,1,2,5,10, 20, 30, 40, 50),labels =c("1/50","1/40", "1/30", "1/20", "1/10", "1/5","1/2","1","2","5","10", "20", "30", "40", "50")
    scale_y_continuous(trans = log_trans(10), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = ylimits) + #limits = c(10^-5.5, 1)
    #geom_hline(aes(yintercept = 0.05), linetype="dashed", color = "brown") +
    annotation_logticks(sides = "b") +
    xlab(xlab) +
    theme_bw() +
    ylab("p value") +
    theme(axis.text.x  = element_text(size=20,angle = 0, hjust = 0.5), axis.text.y  = element_text(size=20), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20))
  return(p)
}

VplotAddSig = function(p, vplot, y_text = 10^-100) {
  a = length(which(vplot$fc < 1))
  b = length(which(vplot$fc > 1))
  c_up = length(which(vplot$fc[vplot$sig == "up"] < 1))
  d_up = length(which(vplot$fc[vplot$sig == "up"] > 1))
  cont_table = rbind(total = c(round(a/(a+b)*(c_up+d_up)), round(b/(a+b)*(c_up+d_up))), geneset = c(c_up, d_up))
  colnames(cont_table) = c("FC<1", "FC>1")
  c = chisq.test(cont_table)
  pvalsigup = c$p.value
  
  a = length(which(vplot$fc < 1))
  b = length(which(vplot$fc > 1))
  c_down = length(which(vplot$fc[vplot$sig == "down"] < 1))
  d_down = length(which(vplot$fc[vplot$sig == "down"] > 1))
  cont_table = rbind(total = c(round(a/(a+b)*(c_down+d_down)), round(b/(a+b)*(c_down+d_down))), geneset = c(c_down, d_down))
  colnames(cont_table) = c("FC<1", "FC>1")
  c = chisq.test(cont_table)
  pvalsigdown = c$p.value
  
  q = p + geom_point_rast(aes(fc,pval, color = sig), data = vplot[vplot$sig %in% c("up", "down"), ], raster.dpi = 100) + scale_color_manual(values = c("blue", "red")) + ggtitle(label = sigtoplot) + 
    #annotate("text", label = sprintf("Expected: %s vs %s \n geneset: %s vs %s \n p = %s ", cont_table[1,1], cont_table[1,2], cont_table[2,1], cont_table[2,2], round(c$p.value, 6)), x = 1, y = 10^-6) + 
    #annotate("text", label = sprintf("sigup: p = %s \n sigdown: p = %s", round(pvalsigup, 6), round(pvalsigdown, 6)), x = 1, y = 10^-6) +
    annotate("text", label = sprintf("total = %s - %s \n sigup: p = %s, observed = %s - %s \n sigdown: p = %s, observed = %s - %s",a,b,round(pvalsigup, 6), c_up, d_up, round(pvalsigdown, 6), c_down, d_down), x = 1, y = y_text)
  return(q)
  
}

FCFCplot = function(vplot, xlab = "fc1", ylab = "fc2", main = "", printgeomtext = T, xlimits = c(0.1, 10),  ylimits = c(0.1, 10)) {
  p = ggplot(data = vplot) + geom_point(aes(x = fc1, y = fc2), colour = "black", alpha = I(1), size = I(0.5)) +
    scale_x_continuous(trans = log_trans(10), breaks = c(0.125,0.5, 1, 2,5, 10),  labels = c(0.125,0.5, 1, 2,5, 10), limits = xlimits) +
    scale_y_continuous(trans = log_trans(10), breaks = c(0.125,0.5, 1, 2,5, 10),  labels = c(0.125,0.5, 1, 2,5, 10), limits = ylimits) +
    geom_hline(aes(yintercept = 1), linetype="dashed", color = "brown") +
    geom_vline(aes(xintercept = 1), linetype="dashed", color = "brown") +
    #geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "brown") +
    annotation_logticks(sides = "bl") +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main) +
    theme_bw() +
    theme(axis.text.x  = element_text(size=15,angle = 0, hjust = 1), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20))
  return(p)
}


# Edited - 20241108_OC
so = readRDS(args[3])
so
DefaultAssay(so) <- "RNA"
head(rownames(so), 50)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)


########
VplotAddSig = function(p, vplot, y_text = 10^-100) {
  if (any(vplot$sig == "up")) {
    a = length(which(vplot$fc < 1))
    b = length(which(vplot$fc > 1))
    c_up = length(which(vplot$fc[vplot$sig == "up"] < 1))
    d_up = length(which(vplot$fc[vplot$sig == "up"] > 1))
    cont_table = rbind(total = c(round(a/(a+b)*(c_up+d_up)), round(b/(a+b)*(c_up+d_up))), geneset = c(c_up, d_up))
    colnames(cont_table) = c("FC<1", "FC>1")
    c = chisq.test(cont_table)
    pvalsigup = c$p.value
  }
  
  if (any(vplot$sig == "down")) {
    a = length(which(vplot$fc < 1))
    b = length(which(vplot$fc > 1))
    c_down = length(which(vplot$fc[vplot$sig == "down"] < 1))
    d_down = length(which(vplot$fc[vplot$sig == "down"] > 1))
    cont_table = rbind(total = c(round(a/(a+b)*(c_down+d_down)), round(b/(a+b)*(c_down+d_down))), geneset = c(c_down, d_down))
    colnames(cont_table) = c("FC<1", "FC>1")
    c = chisq.test(cont_table)
    pvalsigdown = c$p.value
  }
  
  
  q = p + geom_point_rast(aes(fc,pval, color = sig), data = vplot[vplot$sig %in% c("up", "down"), ], raster.dpi = 100) + scale_color_manual(values = c("blue", "red")) + ggtitle(label = sigtoplot) + 
    #annotate("text", label = sprintf("Expected: %s vs %s \n geneset: %s vs %s \n p = %s ", cont_table[1,1], cont_table[1,2], cont_table[2,1], cont_table[2,2], round(c$p.value, 6)), x = 1, y = 10^-6) + 
    #annotate("text", label = sprintf("sigup: p = %s \n sigdown: p = %s", round(pvalsigup, 6), round(pvalsigdown, 6)), x = 1, y = 10^-6) +
    annotate("text", label = sprintf("total = %s - %s \n sigup: p = %s, observed = %s - %s \n sigdown: p = %s, observed = %s - %s",a,b,round(pvalsigup, 6), c_up, d_up, round(pvalsigdown, 6), c_down, d_down), x = 1, y = y_text)
  return(q)
  
}


######

Vplot = function(vplot, xlab = "FC", xlimits = c(0.1, 10), ylimits = c(10^-300,1)) {
  p = ggplot(data = vplot) + geom_point_rast(aes(x = fc, y = pval), colour = "grey", alpha = I(1), size = I(1), raster.dpi = 100) +
    scale_x_continuous(trans = log_trans(10), limits = xlimits) + #breaks = c(1/50,1/40, 1/30, 1/20, 1/10, 1/5,1/2,1,2,5,10, 20, 30, 40, 50),labels =c("1/50","1/40", "1/30", "1/20", "1/10", "1/5","1/2","1","2","5","10", "20", "30", "40", "50")
    scale_y_continuous(trans = log_trans(10), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = ylimits) + #limits = c(10^-5.5, 1)
    #geom_hline(aes(yintercept = 0.05), linetype="dashed", color = "brown") +
    annotation_logticks(sides = "b") +
    xlab(xlab) +
    theme_bw() +
    ylab("p value") +
    theme(axis.text.x  = element_text(size=20,angle = 0, hjust = 0.5), axis.text.y  = element_text(size=20), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20))
  return(p)
}

###############################

# set variables
subgroup = args[4]


#Make tmm objet for limma
# Edited - 20241108_OC
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

count = so@assays$RNA@counts #so[["RNA"]]$counts
dge = DGEList(count)
dge = calcNormFactors(dge)
tmm = new("EList")
message("TMM normalization")
tmm$E = edgeR::cpm(dge, log = TRUE, prior.count = 0.1)
saveRDS(tmm, file = sprintf("igt1_96_%s_20241113_so_tmm_o2.Rds", subgroup))


#######################
# Edited - 20241108_OC

design = data.frame(so@meta.data[,c("annotation_level2","IGTHT")]) 
design$annotation_level2.IGTHT = paste(design$annotation_level2, design$IGTHT, sep = ".")
all(rownames(design) == colnames(so))
design = model.matrix(~ 0 + annotation_level2.IGTHT, data=design) #totalvi_igt1_56_allgenes_Treg_20240529_fit_onCluster_totalvi20240525rmigtsample_Res0.5Sample_idConcat.Rds
colnames(design) = gsub("annotation_level2.IGTHT", "", colnames(design))
cat("ncol=",ncol(design),"rank=", qr(design)$rank,"\n")
fit = lmFit(tmm, design = design)
saveRDS(fit, file = sprintf("igt1_96_%s_20241113_so_fit_annotation_level2IGTHT.Rds", subgroup))




# Custom group comparisons example
library(readxl)
comparisons <- read_xlsx(sprintf("%s_Group_comparisons.xlsx", subgroup))
comparisons[,1]
comparisons[,2]
comparisons[,3]
group_index <- c(which(!is.na(comparisons[,1])), (length(rownames(comparisons)) + 1 ))
gp1_index <- c(group_index[1], group_index[2] -1)
comparisons_custom <- comparisons[gp1_index[1]:gp1_index[2],]
custom_clusters <- c()
for (i in gp1_index[1]:gp1_index[2]) {
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
    }
    
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
    }
    
    custom_clusters <- c(custom_clusters, cl_1, cl_2)
    custom_clusters <- unique(custom_clusters)
}



gp2_index <- c(group_index[2], group_index[3] -1)
comparisons_resting <- comparisons[gp2_index[1]:gp2_index[2],]
resting_clusters <- c()
for (i in gp2_index[1]:gp2_index[2]) {
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
    }
    
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
    }
    
    resting_clusters <- c(resting_clusters, cl_1, cl_2)
    resting_clusters <- unique(resting_clusters)
}

gp3_index <- c(group_index[3], group_index[4] -1)
comparisons_activated <- comparisons[gp3_index[1]:gp3_index[2],]
activated_clusters <- c()
for (i in gp3_index[1]:gp3_index[2]) {
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
    }
    
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
    }
    
    activated_clusters <- c(activated_clusters, cl_1, cl_2)
    activated_clusters <- unique(activated_clusters)
}

gp4_index <- c(group_index[4], group_index[5] - 1)
comparisons_resting_activated <- comparisons[gp4_index[1]:gp4_index[2],]
resting_activated_clusters <- c()
for (i in gp4_index[1]:gp4_index[2]) {
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
    }
    
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
    }
    
    resting_activated__clusters <- c(resting_activated_clusters, cl_1, cl_2)
    resting_activated_clusters <- unique(resting_activated_clusters)
}



gp5_index <- c(group_index[5], group_index[6] -1 )
comparisons_one_v_all <- comparisons[gp5_index[1]:gp5_index[2],]
one_v_all_clusters <- c()
for (i in gp5_index[1]:gp5_index[2]) {
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
    }
    
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
    }
    
    one_v_all_clusters <- c(one_v_all_clusters, cl_1, cl_2)
    one_v_all_clusters <- unique(one_v_all_clusters)
}

list_of_comparisons <- c()
list_of_comparisons <- list(comparisons_custom, comparisons_resting, comparisons_activated, comparisons_resting_activated, comparisons_one_v_all)
list_of_clusters <- list(custom_clusters, resting_clusters, activated_clusters, resting_activated_clusters, one_v_all_clusters)


## Read in Fit
fit = readRDS(sprintf("igt1_96_%s_20241113_so_fit_annotation_level2IGTHT.Rds", subgroup))
design = fit$design


### Cusotm Comparisons - in all (not just healthy)
comparisons[,2]
comparisons[,3]

metadata=so@meta.data
metadata = data.frame(so@meta.data[,c("annotation_level2","organ_simplified","condition_broad", "condition_detailed", "sex", "IGT", "IGTHT")])
metadata$annotation_level2.IGTHT = paste(metadata$annotation_level2, metadata$IGTHT, sep = ".")
metadata = metadata %>% unique()
dim(metadata)
rownames(metadata) = metadata$annotation_level2.IGTHT

contrasts = c()
namescontrasts = c()
for (i in 1:length(rownames(comparisons))) {
    #group1_filter <- expr(annotation_level2 == cl & condition_broad == "healthy")
    if (grepl(",", comparisons[i,2][[1]])) {
        cl_1 <- strsplit(as.character(comparisons[i,2]), ", ")[[1]]
        cl1 <- gsub(", ", comparisons[i,2], replacement = "")
        cl1 <- gsub(sprintf("%s_", subgroup), cl1, replacement = "")
    } 
    else {
        cl_1 <- comparisons[i,2][[1]]
        cl1 <- gsub(sprintf("%s_", subgroup), cl_1, replacement = "")
    }
    if (grepl(",", comparisons[i,3][[1]])) {
        cl_2 <- strsplit(as.character(comparisons[i,3]), ", ")[[1]]
        cl2 <- gsub(", ", comparisons[i,3], replacement = "")
        cl2 <- gsub(sprintf("%s_", subgroup), cl2, replacement = "")
    } 
    else {
        cl_2 <- comparisons[i,3][[1]]
        cl2 <- cl2 <- gsub(sprintf("%s_", subgroup), cl_2, replacement = "")
    }
    
    group1_filter <- expr((annotation_level2 %in% cl_1)) # & condition_broad == "healthy")
    
    #group2_filter <- expr(!(annotation_level2 %in% cl) & condition_broad == "healthy")
    group2_filter <- expr((annotation_level2 %in% cl_2)) # & condition_broad == "healthy")
    
    groups = GetGroups(metadata, group1_filter,group2_filter, "annotation_level2.IGTHT")
    
    nmecontrast = sprintf("%s_vs_%s", cl1, cl2)
    # nmecontrast = CreateComparisonName(group1_filter, group2_filter)
    
    # contrasts = c()
    # namescontrasts = c()
    # # i = 1
    group1 = groups[[1]]
    group2 = groups[[2]]
    contrasts[i] = sprintf("(%s)/%s-(%s)/%s", paste(group1, collapse = "+"), length(group1), paste(group2, collapse = "+"), length(group2))
    namescontrasts[i] = nmecontrast
    # names(contrasts) = namescontrasts
    gene_symbol = rownames(so)
}

cont.matrix = makeContrasts( contrasts = contrasts, levels = design)
colnames(cont.matrix) = namescontrasts
# colnames(cont.matrix) = names(contrasts)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, trend = TRUE, robust = TRUE)
tt_list_all = list()
j = 1
for (i in colnames(cont.matrix)) {
  print(i)
  tt_list_all[[i]] = topTable(fit2,coef = i ,n = Inf, adjust.method = "BH", sort.by = "none")

  if (grepl(",", comparisons[j,2][[1]])) {
  cl_1 <- strsplit(as.character(comparisons[j,2]), ", ")[[1]]
  } else {
  cl_1 <- comparisons[j,2][[1]]
  }
  #cells.1 <- colnames(so[,so@meta.data$annotation_level2 %in% cl_1])
  #cells.1 <- WhichCells(so, idents = cl_1)
  #pct.1 <- round(
  #    x = rowSums(x = so[cells.1, drop = FALSE]) /
  #        length(x = cells.1),
  #    digits = 3
  #)
  
  so_cl1 <- so[,so$annotation_level2 %in% cl_1]
  Idents(so_cl1) <- "annotation_level2"
  expression_info = FetchData(so_cl1, vars = rownames(so_cl1@assays$RNA), cells = colnames(so_cl1))
  pct.1 <- apply(X = expression_info, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  pct.1 <- round(pct.1, 3)

  if (grepl(",", comparisons[j,3][[1]])) {
  cl_2 <- strsplit(as.character(comparisons[j,3]), ", ")[[1]]
  } else {
  cl_2 <- comparisons[j,3][[1]]
  }
  #cells.2 <- colnames(so[,so@meta.data$annotation_level2 %in% cl_2])
  #cells.2 <- WhichCells(so, idents = cl_2)
  #pct.2 <- round(
  #    x = rowSums(x = so[cells.2, drop = FALSE]) /
  #        length(x = cells.2),
  #    digits = 3
  #)

  #so_cl2 <- so[,so$annotation_level2 %in% cl_2]
  #Idents(so_cl2) <- "annotation_level2"
  #pct.2 <- Percent_Expressing(so_cl2, features = rownames(so))
  #pct.2 <- pct.2[,1]

  so_cl2 <- so[,so$annotation_level2 %in% cl_2]
  Idents(so_cl2) <- "annotation_level2"
  expression_info = FetchData(so_cl2, vars = rownames(so_cl2@assays$RNA), cells = colnames(so_cl2))
  pct.2 <- apply(X = expression_info, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  pct.2 <- round(pct.2, 3)

  tt_list_all[[i]]$pct.1 = pct.1
  tt_list_all[[i]]$pct.2 = pct.2
  tt_list_all[[i]]$SYMBOL = gene_symbol
  tt_list_all[[i]]$comparison <- i
  j = j+1
}


tmp = tt_list_all
names(tmp)

pdf("DiffExpression_volcano_All_Comparisons__customcomparisons_annotation_level2.IGTHT.pdf", width = 10, height = 10)
for (i in 1:length(tmp)) {
    print(i)
    vplot = data.frame(SYMBOL = tmp[[i]]$SYMBOL, fc = 2^tmp[[i]]$logFC, pval = tmp[[i]]$P.Value+10^-300, AveExpr = 2^(tmp[[i]]$AveExpr))
    p = Vplot(vplot = vplot, xlab = names(tmp)[i], xlimits = c(2^min(sapply(tmp, function(x) min(x$logFC, na.rm=TRUE))), 2^max(sapply(tmp, function(x) max(x$logFC, na.rm=TRUE)))), ylimits = c(min(vplot$pval), 1)) #c(min(vplot$pval[vplot$pval != 0]),1) # xlimits = c(min(vplot$fc[vplot$fc != 0]), max(vplot$fc))
    sub = (log2(vplot$fc) > 0.5 | log2(vplot$fc) < -0.5)  & vplot$pval < 0.05
    sub = vplot[sub,] %>% arrange(desc(log2(fc)))
    sub1 = sub %>% pull(SYMBOL) %>% head(25) #vplot$SYMBOL %in% as.character(sub2$SYMBOL)[1:200]
    sub2 = sub %>% pull(SYMBOL) %>% tail(25)
    #sub2 = vplot[sub,] %>% arrange(desc(abs(log2(fc))))
    #sub = vplot$SYMBOL %in% as.character(sub2$SYMBOL)[1:50]
    sub = vplot$SYMBOL %in% c(sub1, sub2)#as.character(sub2$SYMBOL)[1:50]
    library(ggrepel)
    print(p + geom_text_repel(data = vplot[sub,], aes(x = fc, y = pval, label = SYMBOL)))
    #print(p + geom_text(data = vplot[sub,], aes(x = fc, y = pval, label = SYMBOL)))
}
dev.off()

#Use to create whole list - we want two separate lists at this stage (UP and DOWN)

#for (i in 1:length(tt_list_all)) {
#    data <- tt_list_all[[i]] %>% arrange(desc(logFC), adj.P.Val) # %>% slice_head(n = 20)
#    data <- data[abs(data$logFC) > 0.5 & data$adj.P.Val < 0.05,]
#    colnames(data) <- paste(names(tt_list_all)[i], colnames(data), sep = ".")
#    #data$comparison <- names(tt_list_all)[i]
#    data <- data[,-c(2:3, 6)]
#    write.csv(data, paste0("Signature_genes_list_",subgroup,"_", names(tt_list_all)[i], ".csv"), quote = F)
#}

# Create a heatmap for each group of comparisons, with only these clusters involved - top 10 up-regulated genes from each comparison



############### SKIP FOR GENERAL CALCULATIONS - COULD BE USED FOR EXPERT CLUSTER COMPARISONS PROPOSED BY MEMBERS OF WORKING GROUPS ##########



## Custom Fit- Virus v healthy in the same clusters
#fit = readRDS("igt1_96_Treg_20241029_so_fit_annotation_level2IGTHT.Rds")
#design = fit$design

#metadata = data.frame(so@meta.data[,c("annotation_level2","organ_simplified","condition_broad", "condition_detailed", "sex", "IGT", "IGTHT")])
#metadata$annotation_level2.IGTHT = paste(metadata$annotation_level2, metadata$IGTHT, sep = ".")
#metadata = metadata %>% unique()
#dim(metadata)
#rownames(metadata) = metadata$annotation_level2.IGTHT

#contrasts = c()
#namescontrasts = c()
#for (cl in unique(metadata$annotation_level2)) {
#    group1_filter <- expr(annotation_level2 %in% cl & condition_broad == "healthy")
#    group2_filter <- expr(annotation_level2 %in% cl & condition_broad == "virus")
#    groups = GetGroups(metadata, group1_filter,group2_filter, "annotation_level2.IGTHT")
#    nmecontrast = sprintf("Healthy_vs_Virus_In%s", cl)
#    # nmecontrast = CreateComparisonName(group1_filter, group2_filter)
#    
#    # contrasts = c()
#    # namescontrasts = c()
#    # # i = 1
#    group1 = groups[[1]]
#    group2 = groups[[2]]
#    contrasts[cl] = sprintf("(%s)/%s-(%s)/%s", paste(group1, collapse = "+"), length(group1), paste(group2, collapse = "+"), length(group2))
#    namescontrasts[cl] = nmecontrast
#    # names(contrasts) = namescontrasts
#    gene_symbol = rownames(so[["RNA"]]$counts)
#}

#cont.matrix = makeContrasts( contrasts = contrasts, levels = design)
#colnames(cont.matrix) = namescontrasts
## colnames(cont.matrix) = names(contrasts)
#fit2 = contrasts.fit(fit, cont.matrix)
#fit2 = eBayes(fit2, trend = TRUE, robust = TRUE)
#tt_list_2 = list()
#for (i in colnames(cont.matrix)) {
#    print(i)
#    tt_list_2[[i]] = topTable(fit2,coef = i ,n = Inf, adjust.method = "BH", sort.by = "none")
#    tt_list_2[[i]]$SYMBOL = gene_symbol
#}

#tmp = tt_list_2
#names(tmp)

#pdf("CustomComparisonsTrial/Healthy_vs_Virus/DiffExpression_volcano_Healthy_vs_Virus_annotation_level2.IGTHT.pdf", width = 10, height = 10)
#for (i in 1:length(tmp)) {
#    print(i)
#    vplot = data.frame(SYMBOL = tmp[[i]]$SYMBOL, fc = 2^tmp[[i]]$logFC, pval = tmp[[i]]$P.Value+10^-300, AveExpr = 2^(tmp[[i]]$AveExpr))
#    p = Vplot(vplot = vplot, xlab = names(tmp)[i], xlimits = c(min(vplot$fc[vplot$fc != 0]), max(vplot$fc)), ylimits = c(min(vplot$pval), 1)) #c(min(vplot$pval[vplot$pval != 0]),1)
#    sub = (log2(vplot$fc) > 0.5 | log2(vplot$fc) < -0.5)  & vplot$pval < 0.05
#    sub = vplot[sub,] %>% arrange(desc(log2(fc)))
#    sub1 = sub %>% pull(SYMBOL) %>% head(25) #vplot$SYMBOL %in% as.character(sub2$SYMBOL)[1:200]
#    sub2 = sub %>% pull(SYMBOL) %>% tail(25)
#    #sub2 = vplot[sub,] %>% arrange(desc(abs(log2(fc))))
#    #sub = vplot$SYMBOL %in% as.character(sub2$SYMBOL)[1:50]
#    sub = vplot$SYMBOL %in% c(sub1, sub2)#as.character(sub2$SYMBOL)[1:50]
#    library(ggrepel)
#    print(p + geom_text_repel(data = vplot[sub,], aes(x = fc, y = pval, label = SYMBOL)))
#    #print(p + geom_text(data = vplot[sub,], aes(x = fc, y = pval, label = SYMBOL)))
#}
#dev.off()

#tt_list_2[[i]] %>% arrange(desc(logFC), adj.P.Val) %>% slice_head(n = 20) 
#tt_list_2[[i]] %>% arrange(logFC, adj.P.Val) %>% slice_head(n = 20) 

####### Gene lists ##########
#change tt_list to suit
#for (i in 1:length(tt_list)) {
#    data <- tt_list_healthy[[i]] %>% arrange(desc(logFC), adj.P.Val) # %>% slice_head(n = 20)
#    data <- data[abs(data$logFC) > 0.5 & data$adj.P.Val < 0.05,]
#    colnames(data) <- paste(names(tt_list_healthy)[i], colnames(data), sep = ".")
#    write.csv(data, paste0("CustomComparisonsTrial/Healthy_vs_Virus/Signature_genes_list_", names(tt_list_healthy)[i], ".csv"), quote = F)
#}



### Exporting gene lists
# Change tt_list for each comparison - in our case one for custom comparisons in Healthy and one for custom comparisons in All

################## SKIP HEALTHY ######################
# Healthy
# save tt_list

#saveRDS(tt_list_healthy, file = "Healthy/ttlist_healthy_20241113.Rds")

#Collapse tt_list
# generic tt_list from here

#tt_list = readRDS("Healthy/ttlist_healthy_20241113.Rds")
#tt = CollapseDiff_limmatrend(tt_list)

#rownames(tt) = do.call(rbind,strsplit(rownames(tt), split = "_"))[,2]

#tt$SYMBOL = rownames(tt)#fData(cds2)$SYMBOL[match(rownames(tt), rownames(cds2))]
#colnames(tt) = gsub(pattern = "AveExpr", replacement = "log2AveExpr", colnames(tt))
#colnames(tt) = gsub(pattern = "logFC", replacement = "log2FC", colnames(tt))

#tt[tt$SYMBOL == "FOXP3",]

#write.table(x = tt, file = "Healthy/ttlist_20241121.txt", sep = "\t", quote = F)

#tt = read.table(file = sprintf("%s/%s/ttlist_20240615.txt", prefix, prefix2), sep = "\t")

# Create 4 separate tt list tables for each of the columns found in a seurat object
# one for FC, one for p value etc ... pct.1, pct.2.....
# First, add pct.1 and pct.2 to the whole tt_list


#Make and export Gene signatures
#tt2 = read.table(file = "Healthy/ttlist_20241121.txt", sep = "\t")
#ttl = list(log2FC = tt2[,grepl("log2FC", colnames(tt2))],
#           log2AveExpr = tt2[,grepl("log2AveExpr", colnames(tt2))],
#           P.Value = tt2[,grepl("P.Value", colnames(tt2))],
#           adj.P.Val = tt2[,grepl("adj.P.Val", colnames(tt2))],
#           pct.1 = tt2[,grepl("pct.1", colnames(tt2))],
#           pct.2 = tt2[,grepl("pct.2", colnames(tt2))])

#ttl = lapply(ttl, function(x) cbind(SYMBOL = tt2$SYMBOL, x))

#ttl$SYMBOL = tt2$SYMBOL
#colnames = gsub(pattern = "^AveExpr_", replacement = "", colnames(ttl$AveExpr))

#for (i in 1:2) { print(colnames(ttl[[i]])) }
#for (i in 1:(length(ttl)-1)) {
#    colnames(ttl[[i]]) = gsub(pattern = "^log2FC_|^log2AveExpr_|^P.Value_|^adj.P.Val_|^pct.1_|^pct.2_", replacement = "", colnames(ttl[[i]]))
#}


#ttl$log2FC["Foxp3",ttl$adj.P.Val["Foxp3",] < 0.05]
#ttl$log2FC["Rorc",ttl$adj.P.Val["Rorc",] < 0.05]
#ttl$adj.P.Val["Rorc",ttl$adj.P.Val["Rorc",] < 0.05]

#Column_tt_list <- c("log2FC", "log2AveExpr", "PValue", "adjPValue", "pct1", "pct2")
#for (i in 1:length(ttl)) {
#    write.table(x = ttl[[i]], file = sprintf("Healthy/ttlist_%s_20241121.csv", Column_tt_list[i]), sep = ",", quote = F, col.names = NA, row.names = T)
#}


#th = 0.5

# gene_list = list()
# for (j in colnames(ttl$log2FC)) {
#     gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[j]), decreasing=T)]
#     gene_list[[sprintf("%s_up", j)]] = gene_ord[gene_ord %in% as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] > th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
#     gene_list[[sprintf("%s_down", j)]] = gene_ord[gene_ord %in%as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] < -th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
# }
# which(unlist(lapply(gene_list, function(x) any(x %in% "Ccr7"))))
# lapply(gene_list, length)

#gene_list = list()
#for (j in colnames(ttl$log2FC)) {
#    #gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[j]), decreasing=T)]
#    gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[[j]]), decreasing=T)]
#    up = gene_ord[gene_ord %in% as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] > th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
#    down = gene_ord[gene_ord %in%as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] < -th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
#    gene_list[[sprintf("%s_up", j)]] = data.frame(SYMBOL = up, log2FC = ttl$log2FC[j][up,], P.Value = ttl$P.Value[j][up,], adj.P.Val = ttl$adj.P.Val[j][up,], pct.1 = ttl$pct.1[j][up,], pct.2 = ttl$pct.2[j][up,])
#    gene_list[[sprintf("%s_down", j)]] = data.frame(SYMBOL = down, log2FC = ttl$log2FC[j][down,], P.Value = ttl$P.Value[j][down,], adj.P.Val = ttl$adj.P.Val[j][down,], pct.1 = ttl$pct.1[j][down,], pct.2 = ttl$pct.2[j][down,])
#}

#which(unlist(lapply(gene_list, function(x) any(x$SYMBOL %in% "Ccr7"))))
#lapply(gene_list, dim)

#gene_list_original = gene_list

#dir.create("DiffExpression_signatures")
#for (i in 1:length(gene_list)) {
#    print(i)
#    write.table(gene_list[[i]], sprintf("Healthy/DavidLists/%s.csv", names(gene_list)[i]), sep = ",", quote = F, row.names = F, col.names = T)
    # write.table(gene_list[[i]], sprintf("./DiffExpression_signatures/%s.txt", names(gene_list)[i]), sep = "\t", quote = F, row.names = F, col.names = F)
#}

#unique(so$annotation_level2)
#ttl$log2FC["Tbx21",grepl("InHealthy$", colnames(ttl$log2FC))]
#ttl$log2FC["Rorc",grepl("InHealthy$", colnames(ttl$log2FC))]
#ttl$log2FC["Gata3",grepl("InHealthy$", colnames(ttl$log2FC))]

# Produce other plots - subgroup heatmap etc
# Heatmap for each group with top 10 genes in each comparison (Only Up)
# Create a heatmap for each group of comparisons, with only these clusters involved - top 10 up-regulated genes from each comparison

#heatmap_genes <- c()
#genes <- c()
#for (i in 1:(length(group_index) - 1)) {
#    lower = (group_index[i] * 2) - 1
#    upper = ((group_index[i+1] -1) * 2)
#    comparison_gene_list <- gene_list[lower:upper]
#    for (k in 1:length(comparison_gene_list)) {
#        if (k %% 2 == 1) {
#            genes <- append(genes, head(comparison_gene_list[[k]]$SYMBOL, 10))
#            print(k)
#        }
#    }
#    
#    heatmap_genes[[length(heatmap_genes) + 1]] <- unique(genes)
#    genes <- c()
#}


# subset seurat object to only include clusters involved in each comparison

#DefaultAssay(so) <- 'RNA'
#so = NormalizeData(so, normalization.method =  "LogNormalize", scale.factor = 10000)

#so = ScaleData(so, verbose = T)

#library(forcats)
#i = 0
#for(group in list_of_comparisons) {
#    i = i+1
#    pdf(sprintf("Healthy/Heatmap_%s_top10_up_genes.pdf", gsub(" ", "_", list_of_comparisons[[i]]$Comparisons[1]), height = 15, width = 8))
#    print(DoHeatmap(so[, so$annotation_level2 %in% list_of_clusters[[i]]], features = heatmap_genes[[i]], group.by = "annotation_level2")) #, raster = F)
#    dev.off()
#}


#pdf("Healthy/Heatmap_top10_up_down_genes.pdf", height = 12, width = 8)
#DoHeatmap(so, features = heatmap_genes, group.by = "annotation_level2")
#dev.off()

#pdf("CustomComparisonsTrial/Ditto_Heatmap_top10_up_down_genes.pdf", height = 12, width = 8)
#dittoHeatmap(so, genes = heatmap_genes, annot.by = "annotation_level2")
#dev.off()


########################################## All #################################################
# save tt_list
saveRDS(tt_list_all, file = "ttlist_20241210.Rds")

#Collapse tt_list
# generic tt_list from here
tt_list = readRDS("ttlist_20241210.Rds")
tt = CollapseDiff_limmatrend(tt_list)
#rownames(tt) = do.call(rbind,strsplit(rownames(tt), split = "_"))[,2]
tt$SYMBOL = rownames(tt)#fData(cds2)$SYMBOL[match(rownames(tt), rownames(cds2))]
colnames(tt) = gsub(pattern = "AveExpr", replacement = "log2AveExpr", colnames(tt))
colnames(tt) = gsub(pattern = "logFC", replacement = "log2FC", colnames(tt))
#tt[tt$SYMBOL == "FOXP3",]

write.table(x = tt, file = "ttlist_20241210.txt", sep = "\t", quote = F)

#tt = read.table(file = sprintf("%s/%s/ttlist_20240615.txt", prefix, prefix2), sep = "\t")

# Create 4 separate tt list tables for each of the columns found in a seurat object
# one for FC, one for p value etc ... pct.1, pct.2.....
# First, add pct.1 and pct.2 to the whole tt_list


#Makr and export Gene signatures
tt2 = read.table(file = "ttlist_20241210.txt", sep = "\t")
ttl = list(log2FC = tt2[,grepl("log2FC", colnames(tt2))],
           log2AveExpr = tt2[,grepl("log2AveExpr", colnames(tt2))],
           P.Value = tt2[,grepl("P.Value", colnames(tt2))],
           adj.P.Val = tt2[,grepl("adj.P.Val", colnames(tt2))],
           pct.1 = tt2[,grepl("pct.1", colnames(tt2))],
           pct.2 = tt2[,grepl("pct.2", colnames(tt2))])
#ttl = lapply(ttl, function(x) cbind(SYMBOL = tt2$SYMBOL, x))
ttl$SYMBOL = tt2$SYMBOL
#colnames = gsub(pattern = "^AveExpr_", replacement = "", colnames(ttl$AveExpr))
for (i in 1:2) { print(colnames(ttl[[i]])) }
for (i in 1:(length(ttl)-1)) {
    colnames(ttl[[i]]) = gsub(pattern = "^log2FC_|^log2AveExpr_|^P.Value_|^adj.P.Val_|^pct.1_|^pct.2_", replacement = "", colnames(ttl[[i]]))
}

#ttl$log2FC["Foxp3",ttl$adj.P.Val["Foxp3",] < 0.05]
#ttl$log2FC["Rorc",ttl$adj.P.Val["Rorc",] < 0.05]
#ttl$adj.P.Val["Rorc",ttl$adj.P.Val["Rorc",] < 0.05]

Column_tt_list <- c("log2FC", "log2AveExpr", "PValue", "adjPValue", "pct1", "pct2")
for (i in 1:length(ttl)) {
    write.table(x = ttl[[i]], file = sprintf("ttlist_%s_20241210.csv", Column_tt_list[i]), sep = ",", quote = F, col.names = NA, row.names = T)
}


th = 0.5

# gene_list = list()
# for (j in colnames(ttl$log2FC)) {
#     gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[j]), decreasing=T)]
#     gene_list[[sprintf("%s_up", j)]] = gene_ord[gene_ord %in% as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] > th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
#     gene_list[[sprintf("%s_down", j)]] = gene_ord[gene_ord %in%as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] < -th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
# }
# which(unlist(lapply(gene_list, function(x) any(x %in% "Ccr7"))))
# lapply(gene_list, length)

gene_list = list()
for (j in colnames(ttl$log2FC)) {
    #gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[j]), decreasing=T)]
    gene_ord = ttl$SYMBOL[order(abs(ttl$log2FC[[j]]), decreasing=T)]
    up = gene_ord[gene_ord %in% as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] > th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
    down = gene_ord[gene_ord %in%as.character(ttl$SYMBOL[ttl$adj.P.Val[j] < 0.05 & ttl$log2FC[j] < -th & !grepl("^Rpl|^Rps",ttl$SYMBOL)])]
    gene_list[[sprintf("%s_up", j)]] = data.frame(SYMBOL = up, log2FC = ttl$log2FC[j][up,], P.Value = ttl$P.Value[j][up,], adj.P.Val = ttl$adj.P.Val[j][up,], pct.1 = ttl$pct.1[j][up,], pct.2 = ttl$pct.2[j][up,])
    gene_list[[sprintf("%s_down", j)]] = data.frame(SYMBOL = down, log2FC = ttl$log2FC[j][down,], P.Value = ttl$P.Value[j][down,], adj.P.Val = ttl$adj.P.Val[j][down,], pct.1 = ttl$pct.1[j][down,], pct.2 = ttl$pct.2[j][down,])
}

which(unlist(lapply(gene_list, function(x) any(x$SYMBOL %in% "Ccr7"))))
lapply(gene_list, dim)

#gene_list_original = gene_list

#dir.create("DiffExpression_signatures")
for (i in 1:length(gene_list)) {
    print(i)
    write.table(gene_list[[i]], paste0(subgroup, "_",names(gene_list)[i] ,".csv"), sep = ",", quote = F, row.names = F, col.names = T)
    # write.table(gene_list[[i]], sprintf("./DiffExpression_signatures/%s.txt", names(gene_list)[i]), sep = "\t", quote = F, row.names = F, col.names = F)
}

unique(so$annotation_level2)
#ttl$log2FC["Tbx21",grepl("InAll$", colnames(ttl$log2FC))]
#ttl$log2FC["Rorc",grepl("InAll$", colnames(ttl$log2FC))]
#ttl$log2FC["Gata3",grepl("InAll$", colnames(ttl$log2FC))]

# Produce other plots - subgroup heatmap etc
# Heatmap for each group with top 10 genes in each comparison (Only Up)
# Create a heatmap for each group of comparisons, with only these clusters involved - top 10 up-regulated genes from each comparison

heatmap_genes <- c()
genes <- c()
for (i in 1:(length(group_index) - 1)) {
    lower = (group_index[i] * 2) - 1
    upper = ((group_index[i+1] -1) * 2)
    comparison_gene_list <- gene_list[lower:upper]
    for (k in 1:length(comparison_gene_list)) {
        if (k %% 2 == 1) {
            genes <- append(genes, head(comparison_gene_list[[k]]$SYMBOL, 10))
            print(k)
        }
    }
    
    heatmap_genes[[length(heatmap_genes) + 1]] <- unique(genes)
    genes <- c()
}


# subset seurat object to only include clusters involved in each comparison
DefaultAssay(so) <- 'RNA'
#so = NormalizeData(so, normalization.method =  "LogNormalize", scale.factor = 10000)
so = ScaleData(so, verbose = T)

#library(forcats)
i = 0
for(group in list_of_comparisons) {
    i = i+1
    pdf(sprintf("Heatmap_%s_top10_up_genes.pdf", as.charachter(substitute(group)), height = 15, width = 8))
    print(DoHeatmap(so[, so$annotation_level2 %in% list_of_clusters[[i]]], features = heatmap_genes[[i]], group.by = "annotation_level2")) #, raster = F)
    dev.off()
}

