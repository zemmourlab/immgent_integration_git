# David Zemmour
# R
# usage: Rscript limma_constrasts_custom.R [path_to_seurat_object] [path_to_tmm_object] [path_to_fit_object] [output_dir] [prefix_file_name]
# LOOK FOR EDIT to find where to edit the script! Essentially in loading metadata and in making contrasts

options(max.print=1000)
options(expressions = 50000)

# Parse arguments
args = commandArgs(TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript limma_constrasts_custom.R [path_to_seurat_object] [path_to_tmm_object] [path_to_fit_object] [output_dir] [prefix_file_name]")
}
path_to_seurat_object = args[1]
path_to_tmm_object = args[2]
path_to_fit_object = args[3]
output_dir = args[4] 
prefix_file_name = args[5]

# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/Treg/igt1_96_Treg_20241216.Rds"
# path_to_tmm_object = "DGE_limma/20241217/igt1_96_Treg_20241216_tmm.Rds"
# path_to_fit_object = "DGE_limma/20241217/igt1_96_Treg_20241216_fit.Rds"
# output_dir = "DGE_limma"
# prefix_file_name = "in_Activated"

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

if (!file.exists(path_to_tmm_object)) {
    stop("The specified TMM object file does not exist: ", path_to_tmm_object)
}

if (!file.exists(path_to_fit_object)) {
    stop("The specified limma fit object file does not exist: ", path_to_fit_object)
}

message("loading R libraries and custom functions")
libs = c("limma", "edgeR", "Seurat","BPCells", "ggplot2","ggrepel", "dplyr", "rlang","reshape2", "RColorBrewer", "pals", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

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

message("loading seurat object")
so = readRDS(file = path_to_seurat_object)

message("loading tmm object")
tmm = readRDS(file = path_to_tmm_object)

message("loading fit object")
fit = readRDS(path_to_fit_object)
design = fit$design

message("loading metadata") ##EDIT
metadata = data.frame(so@meta.data[,c("annotation_level2","annotation_level2_parent1", "annotation_level2.IGTHT", "IGTHT","organ_simplified","condition_broad", "condition_detailed", "sex")])
metadata = metadata %>% unique()
dim(metadata)
rownames(metadata) = metadata$annotation_level2.IGTHT

message("creating constrats") ##EDIT
metadata = metadata[metadata$annotation_level2_parent1 %in% c("activated", "resting"),]
metadata$annotation_level2_parent1 = factor(metadata$annotation_level2_parent1, levels = c("activated", "resting"))

contrasts = c()
namescontrasts = c()
for (cl in levels(metadata$annotation_level2_parent1)[1]) {
    print(cl)
    group1_filter <- expr(annotation_level2_parent1 == cl) ##EDIT
    group2_filter <- expr(!(annotation_level2_parent1 %in% cl) ) #& level2_parent1 == "resting" ##EDIT
    groups = GetGroups(metadata, group1_filter,group2_filter, "annotation_level2.IGTHT")
    nmecontrast = "Activated_vs_Resting" ##EDIT
    # nmecontrast = sprintf("%s_vs_AllResting", cl) ##EDIT
    # nmecontrast = CreateComparisonName(group1_filter, group2_filter)
    
    # contrasts = c()
    # namescontrasts = c()
    # # i = 1
    group1 = groups[[1]]
    group2 = groups[[2]]
    contrasts[cl] = sprintf("(%s)/%s-(%s)/%s", paste(group1, collapse = "+"), length(group1), paste(group2, collapse = "+"), length(group2))
    namescontrasts[cl] = nmecontrast
    # names(contrasts) = namescontrasts
    gene_symbol = rownames(tmm$E)#rownames(so[["RNA"]]$counts)
}

message("contrasts.fit...")
cont.matrix = makeContrasts( contrasts = contrasts, levels = design)
colnames(cont.matrix) = namescontrasts
# colnames(cont.matrix) = names(contrasts)
fit2 = contrasts.fit(fit, cont.matrix)
message("eBayes...")
fit2 = eBayes(fit2, trend = TRUE, robust = TRUE)

tt_list = list()
for (i in colnames(cont.matrix)) {
    print(i)
    tt_list[[i]] = topTable(fit2,coef = i ,n = Inf, adjust.method = "BH", sort.by = "none")
    tt_list[[i]]$SYMBOL = gene_symbol
}

tmp = tt_list
names(tmp)

saveRDS(tt_list, file = sprintf("%s/ttlist_%s.Rds", output_dir, prefix_file_name))
message("Results saved to: ", sprintf("%s/ttlist_%s.pdf", output_dir, prefix_file_name))

pdf(sprintf("%s/DiffExpression_volcano_%s.pdf", output_dir, prefix_file_name), width = 10, height = 10)
for (i in 1:length(tmp)) {
    print(i)
    vplot = na.omit(data.frame(SYMBOL = tmp[[i]]$SYMBOL, fc = 2^tmp[[i]]$logFC, pval = tmp[[i]]$P.Value + 10^-300, AveExpr = 2^(tmp[[i]]$AveExpr))) #dealing with NA
    # vplot = data.frame(SYMBOL = tmp[[i]]$SYMBOL, fc = 2^tmp[[i]]$logFC, pval = tmp[[i]]$P.Value+10^-300, AveExpr = 2^(tmp[[i]]$AveExpr))
    p = Vplot(vplot = vplot, xlab = names(tmp)[i], xlimits = c(min(vplot$fc[vplot$fc != 0]), max(vplot$fc)), ylimits = c(min(vplot$pval), 1)) #c(min(vplot$pval[vplot$pval != 0]),1)
    sub = vplot %>%
        filter(abs(log2(fc)) > 0.5 & pval < 0.05) %>%
        arrange(desc(log2(fc))) %>%
        pull(SYMBOL)
    highlight_genes = unique(c(head(sub, 25), tail(sub, 25)))
    # sub = (log2(vplot$fc) > 0.5 | log2(vplot$fc) < -0.5)  & vplot$pval < 0.05
    # sub = vplot[sub,] %>% arrange(desc(log2(fc)))
    # sub1 = sub %>% pull(SYMBOL) %>% head(25) #vplot$SYMBOL %in% as.character(sub2$SYMBOL)[1:200]
    # sub2 = sub %>% pull(SYMBOL) %>% tail(25)
    highlight_filter = vplot$SYMBOL %in% highlight_genes
    # sub = vplot$SYMBOL %in% c(sub1, sub2)#as.character(sub2$SYMBOL)[1:50]
    library(ggrepel)
    print(p + geom_text_repel(data = vplot[highlight_filter,], aes(x = fc, y = pval, label = SYMBOL)))
    # print(p + geom_text(data = vplot[highlight_filter,], aes(x = fc, y = pval, label = SYMBOL)))
}
dev.off()
message("Plots saved to: ", sprintf("%s/DiffExpression_volcano_%s.pdf", output_dir, prefix_file_name))


