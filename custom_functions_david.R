#custom functions




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

#VPlot and FCFCplot

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
