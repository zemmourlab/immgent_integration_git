#custom functions
#David Zemmour
#source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

extract_numeric = function(x) {
    as.numeric(gsub("\\D", "", x))
}

#UMAP plots for seurat
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
    
    if (is.null(split_by1)) {
        plot2 = p + p2 + scale_color_manual(values = mypal) + theme_bw() + ggtitle(label = sprintf("color: %s", color_by))
    } else if (is.null(split_by2)) {
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
        
        
        if (is.null(split_by1)) {
            plot3 = p + p2  + 
                scale_color_manual(values = c("black", "red")) + 
                scale_alpha_manual(values = c(0.5,1))  + 
                theme_bw()  + 
                #facet_wrap(facets = vars(so@meta.data[,"split_by1"])) + 
                ggtitle(label = sprintf("gene: %s", g))   
        } else if (is.null(split_by2)) {
            plot3 = p + p2  + 
                scale_color_manual(values = c("black", "red")) + 
                scale_alpha_manual(values = c(0.5,1))  + 
                theme_bw()  + 
                facet_wrap(facets = vars(so@meta.data[,"split_by1"])) + 
                ggtitle(label = sprintf("gene: %s, grid: %s", g, split_by1))
        } else {
            plot3 = p + p2  + 
                scale_color_manual(values = c("black", "red")) + 
                scale_alpha_manual(values = c(0.5,1))  + 
                theme_bw() + 
                facet_grid(rows = vars(so@meta.data[,"split_by1"]), cols = vars(so@meta.data[,"split_by2"])) +
                ggtitle(label = sprintf("gene: %s, grid: %s by ", g, color_by, split_by1, split_by2))
        }
        print(plot3)
        
    }
    
}

#Functions for integration_midway3.Rmd

ConvertS5toS3 = function(so, assay1 = "RNA", assay2 = "ADT") {
    so.v3 = CreateAssayObject(counts = so[[assay1]]$counts)
    so.v3 = CreateSeuratObject(so.v3)
    so.v3[[assay2]] = CreateAssayObject(counts = so[[assay2]]$counts) #samples@assays$ADT$counts
    #print(all(colnames(so.v3@assays$RNA$counts) == colnames(so.v3@assays$ADT$counts)))
    #all(colnames(so.v3) == colnames(so))
    so.v3@meta.data = so@meta.data
    return(so.v3)
}

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


###For milo
MyplotDAbeeswarm = function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
    require(ggbeeswarm)
    if (!is.null(group.by)) {
        if (!group.by %in% colnames(da.res)) {
            stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
                 group.by, ")?")
        }
        if (is.numeric(da.res[, group.by])) {
        }
        da.res <- mutate(da.res, group_by = da.res[, group.by])
    }
    else {
        da.res <- mutate(da.res, group_by = "g1")
    }
    if (!is.factor(da.res[, "group_by"])) {
        message("Converting group_by to factor...")
        da.res <- mutate(da.res, group_by = factor(group_by, 
                                                   levels = unique(group_by)))
    }
    if (!is.null(subset.nhoods)) {
        da.res <- da.res[subset.nhoods, ]
    }
    beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(PValue < 
                                                                          alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                             logFC)) + geom_quasirandom())
    pos_x <- beeswarm_pos$data[[1]]$x
    pos_y <- beeswarm_pos$data[[1]]$y
    n_groups <- unique(da.res$group_by) %>% length()
    da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                         1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 1, 
                                                                                logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                            levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y) %>% 
        ggplot(aes(pos_x, pos_y, color = logFC, size = PValue < 
                       alpha))  + scale_color_gradient2() + xlab(group.by) + ylab("Log Fold Change") + #guides(color = "none")
        scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                        seq(1, n_groups))) + geom_point() + coord_flip() + 
        theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}

### To run limma-trend DGE, VPlot and FCFCplot

run_limmatrend_contrasts_counfoundings = function(tmm = tmm, group = group, confoundings = confoundings, formula.mod.matrix = formula.mod.matrix, contrasts =contrasts, gene_symbol= gene_symbol) { #count = count, dge = dge, 
    suppressPackageStartupMessages(library(limma))
    suppressPackageStartupMessages(library(edgeR))
    message("limmatrend")
    #dge = DGEList(count, group = group)
    #dge = calcNormFactors(dge)
    design = model.matrix(formula(formula.mod.matrix))
    #colnames(design) = gsub("confoundings\\[\\[[0-9]\\]\\]", "", colnames(design))
    colnames(design) = gsub(paste(c("group", sprintf("confoundings\\[\\[%s\\]\\]", 1:length(confoundings))), collapse = "|"), "", colnames(design))
    cont.matrix = makeContrasts( contrasts = contrasts, levels = design)
    colnames(cont.matrix) = names(contrasts)
    #tmm = new("EList")
    #message("TMM normalization")
    #tmm$E = edgeR::cpm(dge, log = TRUE, prior.count = 0.1)
    message("lmFit")
    fit = lmFit(tmm, design = design)
    fit2 = contrasts.fit(fit, cont.matrix)
    fit2 = eBayes(fit2, trend = TRUE, robust = TRUE)
    
    tt_list = list()
    for (i in colnames(cont.matrix)) {
        print(i)
        tt_list[[i]] = topTable(fit2,coef = i ,n = Inf, adjust.method = "BH", sort.by = "none")
        tt_list[[i]]$SYMBOL = gene_symbol
    }
    
    return(tt_list)
}


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

CollapseDiff_limmatrend = function(l) {
    l = lapply(l, function(x) x = data.frame(x, rownames = rownames(x)))
    
    fc = data.frame(l[[1]][,c("rownames", "logFC")])
    for (i in 2:length(l)) { fc = merge(x = fc, y = l[[i]][,c("rownames", "logFC")], by.x = "rownames", by.y = "rownames", all = T) }
    
    mean = data.frame(l[[1]][,c("rownames", "AveExpr")])
    for (i in 2:length(l)) { mean = merge(x = mean, y = l[[i]][,c("rownames", "AveExpr")], by.x = "rownames", by.y = "rownames", all = T) }
    
    pv = data.frame(l[[1]][,c("rownames", "P.Value")])
    for (i in 2:length(l)) { pv = merge(x = pv, y = l[[i]][,c("rownames", "P.Value")], by.x = "rownames", by.y = "rownames", all = T) }
    
    qcl = data.frame(l[[1]][,c("rownames", "adj.P.Val")])
    for (i in 2:length(l)) { qcl = merge(x = qcl, y = l[[i]][,c("rownames", "adj.P.Val")], by.x = "rownames", by.y = "rownames", all = T) }
    
    rownames(fc) = fc[,1]
    rownames(mean) = fc[,1]
    rownames(pv) = pv[,1]
    rownames(qcl) = qcl[,1]
    fc = fc[,-1]
    mean = mean[,-1]
    pv = pv[,-1]
    qcl = qcl[,-1]
    colnames(fc) = paste("logFC", names(l), sep = "_")
    colnames(mean) = paste("AveExpr", names(l), sep = "_")
    colnames(pv) = paste("P.Value", names(l), sep = "_")
    colnames(qcl) = paste("adj.P.Val", names(l), sep = "_")
    all(rownames(fc) == rownames(pv))
    all(rownames(fc) == rownames(mean))
    all(rownames(fc) == rownames(qcl))
    fc[is.na(fc)] = 0
    pv[is.na(pv)] = 1
    qcl[is.na(qcl)] = 1
    diff = data.frame(fc, mean, pv, qcl)
    which(is.na(diff))
    
    return(diff)
    
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

#splitColors
SplitColors = function(pal = mypal[1:length(levels(cds$clustersCCA9))], splitvector = sapply(levels(cds$clustersCCA9), function(x) { length(which(grepl(x, levels(cds$clustersCCA9Split))))}) ) {
    newpal = c()
    for (c in 1:length(splitvector)){
        #print(c)
        n = splitvector[c]
        newcol = matrix(rep(col2rgb(pal[c], alpha = T), n), ncol = n)
        newcol[4,] = seq(50, 255, length.out = n)
        newpal = c(newpal, rgb(red = newcol[1,]/255, blue = newcol[2,]/255, green = newcol[3,]/255, alpha = newcol[4,]/255))
    }
    return(newpal)
}
sp = SplitColors
#pie(1:length(levels(cds$clustersCCA9Split)), col = sp)


#removeDuplicateColumns

removeDuplicateColumns <- function(df) {
    # Determine unique columns by comparing all pairs
    uniqueCols <- !duplicated(as.list(df))
    # Subset dataframe to keep only unique columns
    dfUnique <- df[, uniqueCols]
    print(colnames(df)[!uniqueCols])
    return(dfUnique)
}




