```{r}
pdf("AlluvialPlot_level2_tissue3.pdf", 25, 8, useDingbats = F) #no more than 1000 cells per organ before subsetting gdT
df = so_orig@meta.data %>% filter(condition_broad == "healthy" & grepl(pattern = "allT|CD45p", x = target_cells_simplified)) %>% group_by(organ_simplified) %>% slice_sample(n = 1000) %>% ungroup() %>% filter(annotation_level1 %in% "nonconv" & & !(annotation_level2 %in% c("nonconv_prolif", "nonconv_miniverse"))) %>% group_by(annotation_level2, organ_simplified) %>% summarize(count = n()) %>% as.data.frame()


df$annotation_level2 = factor(df$annotation_level2, levels = rev(c("gdT_cl2","gdT_cl11","gdT_cl28","gdT_cl35","gdT_cl4","gdT_cl12","gdT_cl42","gdT_cl23","gdT_cl52","gdT_cl54","gdT_cl33","gdT_cl26","gdT_cl10","gdT_cl9","gdT_cl46","gdT_cl5","gdT_cl43","gdT_cl50","gdT_cl7")))
df$organ_simplified = factor(df$organ_simplified, levels = rev(c("small intestine","colon","mammary gland","liver","peritoneal cavity","lung","blood","spleen","LN","bone marrow","thymus","kidney","submandibular gland","placenta","uterus","skin")))
p = ggplot(df,
           aes(y = count, axis1 = organ_simplified, axis2 = annotation_level2)) +
    coord_flip() +
    geom_alluvium(aes(fill = annotation_level2), width = 1/5) +
    geom_stratum(width = 1/5, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle = 90) +
    # geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)),
    #                 size = 5,
    #                 min.segment.length = 0, # Ensures a line is drawn even if the label is close
    #                 box.padding = 0.5,      # Adjusts space around text before line is drawn
    #                 point.padding = 0.5,    # Adjusts space around the point (stratum)
    #                 max.overlaps = Inf) +   # Allows all labels to be drawn, even if they overlap initially
    scale_x_discrete(limits = c("organ_simplified", "annotation_level2"), expand = c(.05, .05)) +
    scale_fill_manual(values = mypal_level2) + # Or scale_fill_brewer
    ggtitle("Level2 across Tissues at baseline") +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(), # Removes major grid lines
        panel.grid.minor = element_blank(), # Removes minor grid lines
        axis.title = element_blank(),       # Removes axis titles (e.g., "count")
        axis.text = element_blank(),        # Removes axis labels/text (e.g., numbers on the axis)
        axis.ticks = element_blank()        # Removes axis tick marks
    )
print(p)
p = ggplot(df,
           aes(y = count, axis1 = organ_simplified, axis2 = annotation_level2)) +
    coord_flip() +
    geom_alluvium(aes(fill = organ_simplified), width = 1/5) +
    geom_stratum(width = 1/5, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle = 90) +
    # geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)),
    #                 size = 5,
    #                 min.segment.length = 0, # Ensures a line is drawn even if the label is close
    #                 box.padding = 0.5,      # Adjusts space around text before line is drawn
    #                 point.padding = 0.5,    # Adjusts space around the point (stratum)
    #                 max.overlaps = Inf) +   # Allows all labels to be drawn, even if they overlap
    scale_x_discrete(limits = c("organ_simplified", "annotation_level2"), expand = c(.05, .05)) +
    scale_fill_manual(values = mypal_organ) + # Or scale_fill_brewer
    ggtitle("Level2 across Tissues at baseline") +
    theme_minimal() +
    theme(
        panel.grid.major = element_blank(), # Removes major grid lines
        panel.grid.minor = element_blank(), # Removes minor grid lines
        axis.title = element_blank(),       # Removes axis titles (e.g., "count")
        axis.text = element_blank(),        # Removes axis labels/text (e.g., numbers on the axis)
        axis.ticks = element_blank()        # Removes axis tick marks
    )
print(p)
dev.off()

```



**Topic DGE**
    
    ```{r}

TopicDGEPlots = function(F1, L1, group1, group2, title_plot = "") {
    
    require(ggplot2)
    require(ggrepel)
    require(dplyr)
    require(scattermore)
    
    # Calculate the loadings
    loadings_group1 = colMeans(L1[group1,])
    loadings_group2 = colMeans(L1[group2,])
    loadings_groups = colMeans(L1[c(group1, group2),])
    
    # Mean genes for each group
    mean_genes_group1 = F1 %*% loadings_group1
    mean_genes_group2 = F1 %*% loadings_group2
    mean_genes = F1 %*% colMeans(L1[c(group1, group2),])
    
    # Fold Change between the two groups
    fc_loadings = loadings_group1 - loadings_group2
    fc_genes = F1 %*% fc_loadings %>% as.data.frame()  # take exp to have the actual FC
    
    # Create vplot for fold changes (MA plot of factors)
    vplot = data.frame(SYMBOL = names(fc_loadings), log2FC = fc_loadings / log(2), AveExpr = loadings_groups)
    
    # Maximum FC for symmetrical x-axis range
    max_fc = ceiling(max(abs(vplot$log2FC)))  # Get the max of the absolute FC values and round up
    top_genes = vplot %>%
        arrange(desc(abs(log2FC))) %>%
        head(50)
    
    # MA Plot (Factors)
    p1 = ggplot(data = vplot) + 
        geom_point(aes(x = log2FC, y = AveExpr), colour = "black", alpha = I(1), size = I(1)) +
        xlim(-max_fc, max_fc) +  # Set symmetrical x-axis range
        geom_text_repel(data = top_genes, aes(x = log2FC, y = AveExpr, label = SYMBOL), 
                        size = 3, color = "red", box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps  = 20) +
        theme_minimal() +
        labs(
            x = "Fold Change (log2)",
            y = "Average Expression",
            title = title_plot
        )
    
    # Create vplot for gene expression (MA plot of genes)
    vplot_genes = data.frame(SYMBOL = rownames(fc_genes), log2FC = fc_genes[,1] / log(2), AveExpr = mean_genes[,1])
    
    # Maximum FC for symmetrical x-axis range (for genes)
    max_fc_genes = ceiling(max(abs(vplot_genes$log2FC)))  # Get the max of the absolute FC values and round up
    top_genes_genes <- vplot_genes %>%
        arrange(desc(abs(log2FC))) %>%
        head(50)
    
    # MA Plot (Genes)
    p2 <- ggplot(data = vplot_genes) + 
        geom_scattermore(aes(x = log2FC, y = AveExpr), colour = "black", alpha = I(1), size = I(1), pixels = c(512, 512)) +
        xlim(-max_fc_genes, max_fc_genes) +  # Set symmetrical x-axis range
        geom_text_repel(data = top_genes_genes, aes(x = log2FC, y = AveExpr, label = SYMBOL), 
                        size = 3, color = "red", box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps  = 20) +
        theme_minimal() +
        labs(
            x = "Fold Change (log2)",
            y = "Average Expression",
            title = title_plot
        )
    
    # Return the plots and vplot tables
    list(p1 = p1, p2 = p2, diff_factors = vplot, diff_genes = vplot_genes)
}

F1 = read.table(sprintf("%s/%s/gene_factor_matrix.txt", prefix, prefix2), header = T, sep = "\t") %>% as.matrix()
L1 = read.table(sprintf("%s/%s/cell_factor_matrix.txt", prefix, prefix2), header = T, sep = "\t") %>% as.matrix()
mdata = read.table(sprintf("%s/%s/metadata_20250513_mde.txt", prefix, prefix2), header = T, sep = "\t")

```

```{r}
# Example usage

group1 = mdata_orig %>% filter(annotation_level2 %in% "nonconv_cl8" & organ_simplified == "uterus" & condition_broad == "healthy") %>% pull(cellID)
group2 = mdata_orig %>% filter(annotation_level2 %in% "nonconv_cl8" & organ_simplified == "skin" & condition_broad == "healthy") %>% pull(cellID)
results = TopicDGEPlots(F1, L1, group1, group2, title_plot = "nonconv_cl8 uterus vs skin")

group1 = mdata_orig %>% filter(annotation_level2 %in% "nonconv_cl8" & iNKT & condition_broad == "healthy") %>% pull(cellID)
group2 = mdata_orig %>% filter(annotation_level2 %in% "nonconv_cl8" & MAIT & condition_broad == "healthy") %>% pull(cellID)
results = TopicDGEPlots(F1, L1, group1, group2, title_plot = "nonconv_cl8 iNKT vs MAIT")



results$p1  # First plot (MA plot of factors)
results$p2  # Second plot (MA plot of genes)
results$diff_factors  # vplot for fold changes
results$diff_genes %>% arrange(desc(log2FC)) %>% head(50) # vplot for gene expression
results$diff_genes %>% arrange(desc(log2FC)) %>% tail(50) 

# Call the function with F1, L1, and group data
results = TopicDGEPlots(F1, L1, group1, group2,, title_plot = "Treg_cl8 vs Treg_cl2 in colon healthy")

# Access the plots and data
results$p1  # First plot (MA plot of factors)
results$p2  # Second plot (MA plot of genes)
results$diff_factors  # vplot for fold changes
results$diff_genes  # vplot for gene expression


```

