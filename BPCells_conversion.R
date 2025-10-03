# convert_to_bpcells.R
suppressPackageStartupMessages({
    library(Seurat)
})

# You need BPCells installed for the conversion step:
# remotes::install_github("bnprks/BPCells")

convert_seurat_to_bpcells <- function(
        in_rds,
        out_rds,
        assays = NULL,             # default: all assays in object
        convert_data = FALSE,      # also back the normalized "data" matrix (optional)
        out_parent = NULL,         # where to place BPCells dirs; default sits next to out_rds
        overwrite = TRUE
) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
        stop("BPCells is not installed. Install with remotes::install_github('bnprks/BPCells').")
    }
    
    obj <- readRDS(in_rds)
    if (is.null(assays)) assays <- Assays(obj)
    
    if (is.null(out_parent)) {
        out_parent <- file.path(
            dirname(out_rds),
            paste0(tools::file_path_sans_ext(basename(out_rds)), "_bpc")
        )
    }
    dir.create(out_parent, showWarnings = FALSE, recursive = TRUE)
    
    for (assay in assays) {
        if (!assay %in% Assays(obj)) next
        message("Converting assay: ", assay)
        a <- obj[[assay]]
        
        # Counts -> BPCells
        counts_mat <- try(Seurat::GetAssayData(a, slot = "counts"), silent = TRUE)
        if (!inherits(counts_mat, "try-error") && !is.null(counts_mat)) {
            if (!any(grepl("BPCells", class(counts_mat)))) {
                dir_counts <- file.path(out_parent, paste0(assay, "_counts"))
                BPCells::write_matrix_dir(counts_mat, dir = dir_counts, overwrite = overwrite)
                a <- Seurat::SetAssayData(a, slot = "counts", new.data = BPCells::open_matrix_dir(dir_counts))
            } else {
                message("  counts already BPCells-backed; skipping")
            }
        }
        
        # Normalized "data" -> BPCells (optional)
        if (isTRUE(convert_data)) {
            data_mat <- try(Seurat::GetAssayData(a, slot = "data"), silent = TRUE)
            if (!inherits(data_mat, "try-error") && !is.null(data_mat)) {
                if (!any(grepl("BPCells", class(data_mat)))) {
                    dir_data <- file.path(out_parent, paste0(assay, "_data"))
                    BPCells::write_matrix_dir(data_mat, dir = dir_data, overwrite = overwrite)
                    a <- Seurat::SetAssayData(a, slot = "data", new.data = BPCells::open_matrix_dir(dir_data))
                } else {
                    message("  data already BPCells-backed; skipping")
                }
            }
        }
        
        obj[[assay]] <- a
    }
    
    saveRDS(obj, out_rds)
    invisible(list(rds = out_rds, bpc_dir = out_parent))
}

# ---- Example usage ----
# convert_seurat_to_bpcells(
#   in_rds       = "my_big_seurat.rds",
#   out_rds      = "my_big_seurat_bpc.rds",
#   assays       = c("RNA","SCT"),       # or leave NULL for all assays
#   convert_data = FALSE                  # set TRUE to also back the 'data' matrix
# )

convert_seurat_to_bpcells(
  in_rds       = "/Users/david/Desktop/Treg_App/igt1_96_withtotalvi20250513_clean_Treg.Rds",
  out_rds      = "/Users/david/Desktop/Treg_App/igt1_96_withtotalvi20250513_clean_Treg_bpc.Rds",
  assays       = c("RNA"),       # or leave NULL for all assays
  convert_data = FALSE                  # set TRUE to also back the 'data' matrix
)
