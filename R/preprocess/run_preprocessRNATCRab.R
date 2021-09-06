
message("Merge seurat-objects with TCRab-info...")

mergeTCRtoSeurat <- function(seurat_object, tcr_df){

    ## Merge with TCRab data with seurat-object metadata
    # seurat_object@meta.data <- seurat_object@meta.data[,-1]
    # seurat_object$barcode <- colnames(seurat_object)
    # tcr_df$barcode_uniq <- gsub("\\-", "_", tcr_df$barcode_uniq)
  
    metadata_temp <- merge(seurat_object@meta.data, dplyr::select(tcr_df, -barcode), all.x = T, by.x = "barcode", by.y = "barcode_uniq")
    metadata_temp <- metadata_temp[match(colnames(seurat_object), metadata_temp$barcode), ]
    dim(metadata_temp)
    
    ## Add some meta data;

    # ## 1) Major: over 10 cells in clonotype
    major_clonotypes               <- unique(subset(metadata_temp, metadata_temp$frequency > 10)$new_clonotypes_id)
    metadata_temp$major_clonotypes <- metadata_temp$new_clonotypes_id
    metadata_temp$major_clonotypes[!metadata_temp$new_clonotypes_id %in% major_clonotypes] <- "minor"
     
     
    ## 2) Expanded: over 2 cells in clonotype
    expanded_clonotypes <- unique(subset(metadata_temp, metadata_temp$frequency > 2)$new_clonotypes_id)
    metadata_temp$expanded_clonotypes <- metadata_temp$new_clonotypes_id
    metadata_temp$expanded_clonotypes[!metadata_temp$new_clonotypes_id %in% expanded_clonotypes] <- "unexpanded"

    ## Add metadata into Seurat object; make sure that the colnames match
    rownames(metadata_temp) <- metadata_temp$barcode
    colnames(seurat_object) == rownames(metadata_temp)
    seurat_object@meta.data <- metadata_temp
    return(seurat_object)

}

tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt") 
lag3_seurat <- removeTCRabData(lag3_seurat)
lag3_seurat <- mergeTCRtoSeurat(seurat_object = lag3_seurat, tcr_df = tot_barcode)
saveRDS(lag3_seurat, "results/lag3_seurat_latest.rds")

lag3_seurat@meta.data$tcrab <- ifelse(is.na(lag3_seurat$new_clonotypes_id), "no", "yes")
table(lag3_seurat$tcrab, lag3_seurat$orig.ident)



