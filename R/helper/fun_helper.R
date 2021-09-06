
putLatentsSeurat <- function(seurat_object, latent){
  
  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)
  
  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)
  
  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)
  
  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))
  
  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}



getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){
  
  hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
  blueprint <- SingleR::BlueprintEncodeData()
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not 
  
  if(!is.null(sample)){
    
    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])
    
  }
  
  sce       <- as.SingleCellExperiment(seurat_object)
  
  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)
    
    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.blu$first.labels
      return(seurat_object)  
    }
    
    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }
    
  }
  
  
  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}

getSCviCode <- function(seurat_object, folder_name){
  
  #### Make script for scVI
  names <- unique(seurat_object$orig.ident) %>% gsub(pattern = "\\.", replacement = "\\_")
  names <- names %>% gsub(pattern = "\\-", replacement = "\\_")
  names <- ifelse(grepl('[0-9]', substr(names, 1, 1)), paste0("t", names), names)
  
  df <- list()
  i <- 1 
  
  for(name in names){
    df[[i]] <- paste0(name,  " = CsvDataset(filename='", folder_name, name, ".csv', save_path='', sep=',', new_n_genes=False)")
    i <- i + 1
  }
  
  a <- do.call(rbind, df) 
  
  df <- list()
  i <- 1 
  
  for(name in names){
    df[[i]] <- paste0(name, ".X,")
    i <- i + 1
  }
  
  b <- do.call(rbind, df)
  
  rbind(a, b)
  
}


## For GSEA
forGSEA <- function(de_df){
  
  ## Create .rnk files for GSEA
  rnk  <- de_df %>% arrange(desc(avg_logFC)) %>% dplyr::select(gene, avg_logFC)
  colnames(rnk) <- c('Name', 'metric')
  return(rnk)
  
}



preprocessSeurat <- function(orig_object, cells.to.use){

  clonality_genes <- getClonalityGenes(orig_object)  
  unwanted_genes <- getUnwantedGenes(orig_object)  
  
  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)
  
  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta
  
  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)
  
  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]
  
  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))
  
  ## Scale data
  object <- ScaleData(object, features = hvg)
  
  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))
  
  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  
  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df
  
  return(object)
  
}





plotClustering <- function(seurat_object){
  
  res <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]
  
  q <- NULL; i <- 1
  
  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }
  
  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()
  
}


getLatentClustering <- function(seurat_object){
  
  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)
  
}

readGliphFile <- function(filename){
  
  gliph_df   <- fread(filename, fill = T, stringsAsFactors = F)
  
  gliph_df$vb_score <- as.numeric(gliph_df$vb_score)
  gliph_df$Fisher_score <- as.numeric(gliph_df$Fisher_score)
  gliph_df$final_score <- as.numeric(gliph_df$final_score)
  gliph_df$number_unique_cdr3 <- as.numeric(gliph_df$number_unique_cdr3)
  
  gliph_df$number_subject <- as.numeric(gliph_df$number_subject)
  gliph_df$length_score <- as.numeric(gliph_df$length_score)
  gliph_df$cluster_size_score <- as.numeric(gliph_df$cluster_size_score)
  
  ## Remove everything after last numeric index
  max_index  <- max(as.numeric(gliph_df$index), na.rm = T)
  max_indexs <- which(gliph_df$index == max_index)
  gliph_df   <- gliph_df[1:max_indexs[length(max_indexs)], ]
  return(gliph_df)
  
}

readFACS <- function(path){
  
  # path = "facs/anna/lag3_tube4_final.txt"
  tube_raw <- fread(path)
  tubename <- substr(path, 16, 20)
  colnames(tube_raw) <- gsub("\\+", "pos", colnames(tube_raw))
  colnames(tube_raw) <- gsub("\\-", "neg", colnames(tube_raw))
  colnames(tube_raw) <- make.names(colnames(tube_raw)) %>% make.unique() 
  
  if("FM.ID" %in% colnames(tube_raw)){
    tube_raw <- tube_raw %>% dplyr::rename("V1" = "FM.ID")
  }
  
  ## Remove all other than mel naive or not?
  # tube_raw <- tube_raw %>% filter(Cohort == 1)
  # tube_raw <- tube_raw[,-c(2:4)]
  
  table(is.na(tube_raw))
  m <- tube_raw
  
  m[rowSums(is.na(m)) != ncol(m), ]
  
  tube_raw <- tube_raw %>% 
    mutate(timepoint = as.factor(Sample)) %>% 
    mutate(timepoint = plyr::revalue(timepoint, c("1" = "0m", "2" = "1m", "3" = "3m", "4" = "6m"))) %>% 
#    dplyr::rename(name = V1) %>% 
    dplyr::select(-Sample) %>% mutate(key = paste(V1, timepoint)) # %>% select(-name,-timepoint)
  
}


require(clusterProfiler)
require(org.Hs.eg.db)

if(me == "hru"){
  
  hallmark   <- read.gmt("/Users/hru/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
  
}


if(me == "janihuuh"){
  
  hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")
  
}




## == No need to modify

plotHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  
  if(table(enrich@result$p.adjust < 0.05) %>% length() > 1){
    heatplot(enrich)
  }
  
  else(NULL)
  
}

getHypergeometric <- function(genes_df, universe_df, term_df){
  
  require(clusterProfiler)
  
  if(nrow(genes_df) == 0) return(NULL)
  
  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)
  
  # out: df with enrichment results
  
  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
  enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(enrich)
  
}



getDEGbyCluster <- function(seurat_object, cluster){
  
  message(paste0("===== ", cluster, " ====="))
  
  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= 50) return(NULL)
  
  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$timepoint
  
  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()
  
  n1 <- n_df %>% filter(Var1 == 0) %>% pull(Freq) >= 50
  n2 <- n_df %>% filter(Var1 == 1) %>% pull(Freq) >= 50
  n3 <- n_df %>% filter(Var1 == 3) %>% pull(Freq) >= 50
  
  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  
  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL
  
  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "0", ident.2 = "1", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "2v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "0", ident.2 = "3", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "1", ident.2 = "3", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  
  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2) 
  
  if(!is.null(df)) df <- df %>% 
    # filter(p_val_adj < 0.05) %>% 
    mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)
  
}


removeModelNames <- function(epitopes){
  
  
  epitopes <- gsub("\\_comb\\_tra", "", epitopes)
  epitopes <- gsub("\\_comb\\_trb", "", epitopes)
  epitopes <- gsub("\\_comb\\_trab", "", epitopes)
  
  epitopes <- gsub("\\_comb\\_cdr3a", "", epitopes)
  epitopes <- gsub("\\_comb\\_cdr3b", "", epitopes)
  epitopes <- gsub("\\_comb\\_cdr3ab", "", epitopes)
  
  epitopes <- gsub("\\_cdr3ab", "", epitopes)
  epitopes <- gsub("\\_cdr3a", "", epitopes)
  epitopes <- gsub("\\_cdr3b", "", epitopes)
  epitopes <- gsub("\\_cdr3v", "", epitopes)
  
  epitopes <- gsub("\\_trab", "", epitopes)
  epitopes <- gsub("\\_tra", "", epitopes)
  epitopes <- gsub("\\_trb", "", epitopes)
  
  return(epitopes)
  
}


## Rename thresholds
renameEpitopes <- function(thresholds){
  
  thresholds <- plyr::revalue(thresholds, c(PKYVKQNTLKLAT_cdr3b = "IAV_HA_PKY_cdr3b",
                                            GILGFVFTL_cdr3b = "IAV_M1_GIL_cdr3b",
                                            GILGFVFTL_cdr3ab = "IAV_M1_GIL_cdr3ab",
                                            
                                            GLCTLVAML_cdr3ab = "EBV_BMLF1_GLC_cdr3ab",
                                            GLCTLVAML_cdr3b = "EBV_BMLF1_GLC_cdr3b",
                                            
                                            NLVPMVATV_cdr3b = "CMV_p65_NLV_cdr3b",
                                            TPRVTGGGAM_cdr3b = "CMV_p65_TPR_cdr3b",
                                            RAKFKQLL_cdr3b = "EBV_BZLF1_RAF_cdr3b",
                                            IPSINVHHY_cdr3b = "CMV_p65_IPS_cdr3b",
                                            YVLDHLIVV_cdr3b = "EBV_BRLF1_YVL_cdr3b",
                                            
                                            meloe1_cdr3b  = "MAA_MELOE1_FKY_cdr3b",
                                            
                                            PKYVKQNTLKLAT_cdr3b = "IAV_HA1_PKY_cdr3b",
                                            RPRGEVRFL_cdr3b = "HSV2_B7_RPR_cdr3b",
                                            
                                            melana_cdr3b = "MAA_MART1_ELA_cdr3b",
                                            
                                            
                                            ELAGIGILTV_trab = "MAA_MART1_ELA_trab",
                                            ELAGIGILTV_cdr3ab = "MAA_MART1_ELA_cdr3ab",
                                            ELAGIGILTV_tra_comb = "MAA_MART1_ELA_comb_tra",
                                            ELAGIGILTV_cdr3a_comb = "MAA_MART1_ELA_comb_cdr3a",
                                            ELAGIGILTV_trb_comb = "MAA_MART1_ELA_comb_trb",
                                            ELAGIGILTV_cdr3b_comb = "MAA_MART1_ELA_comb_cdr3b",
                                            
                                            AMFWSVPTV_trb = "MAA_TKT_AMF_trb",
                                            
                                            FLYNLLTRV_trb = "MAA_SEC24_FLY_trb",
                                            FLYNLLTRV_cdr3b = "MAA_SEC24_FLY_cdr3b",
                                            
                                            mart1_trb = "MAA_MART1_AAG_trb",
                                            mart1_cdr3b = "MAA_MART1_AAG_cdr3b"
  ), warn_missing = F)
  
  return(thresholds)
  
}



extractFACSfeature <-  function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\..*", "", str1)
}


reorderClusters <- function(cluster_vec){
  
  ## Get clusters in order
  clusters <- cluster_vec %>% unique() 
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)
  
}

extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}



plotLolliplot <- function(viz_df, timepoint_temp){
  
  df1 <- viz_df %>% group_by(cluster, timepoint, overall) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% filter(timepoint == timepoint_temp & overall == "R")
  df2 <- viz_df %>% group_by(cluster, timepoint, overall) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>% filter(timepoint == timepoint_temp & overall == "N")
  
  df_tot <- left_join(df1, df2, by = "cluster") %>% mutate(log2fc = log2(freq.x/freq.y)) %>%
    mutate(dir = ifelse(log2fc > 1, "up", "unsigf")) %>%
    mutate(dir = ifelse(log2fc < -1, "down", dir))
  
  max_y <- abs(max(df_tot$log2fc))
  
  # ggplot(df_tot, aes(cluster,log2fc,fill=dir)) + geom_bar(stat = "identity") + coord_flip() + ylim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon"))
  
  ggplot(df_tot, aes(log2fc, reorder(cluster, log2fc), fill=dir, size = n.x)) +
    geom_segment(aes(x = 0, xend = log2fc, y=reorder(cluster, log2fc), yend = reorder(cluster, log2fc)), color = "lightgrey", size = 0.5) +
    geom_point(shape = 21) + geom_vline(xintercept = 0)  + geom_vline(xintercept = -1, linetype = "dotted") + geom_vline(xintercept = 1, linetype = "dotted") +
    xlim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon")) + labs(y = "", size = "nCells", fill = "") + add_guide
  
}



plotDynamicsLolliplot <- function(viz_df, overall_temp){
  
  df1 <- viz_df %>% filter(timepoint == 0 & overall == overall_temp) %>% group_by(cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
  df2 <- viz_df %>% filter(timepoint == 1 & overall == overall_temp) %>% group_by(cluster) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
  
  df_tot <- left_join(df1, df2, by = "cluster") %>% mutate(log2fc = log2(freq.x/freq.y)) %>%
    mutate(dir = ifelse(log2fc > 1, "post treatment", "unsigf")) %>%
    mutate(dir = ifelse(log2fc < -1, "pre treatment", dir))
  
  max_y <- abs(max(df_tot$log2fc))
  
  # ggplot(df_tot, aes(cluster,log2fc,fill=dir)) + geom_bar(stat = "identity") + coord_flip() + ylim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon"))
  
  ggplot(df_tot, aes(log2fc, cluster, fill=dir, size = n.x)) +
    geom_segment(aes(x = 0, xend = log2fc, y=cluster, yend = cluster), color = "lightgrey", size = 0.5) +
    geom_point(shape = 21) + geom_vline(xintercept = 0)  + geom_vline(xintercept = -1, linetype = "dotted") + geom_vline(xintercept = 1, linetype = "dotted") +
    xlim(values = c(-max_y,max_y)) + scale_fill_manual(values = c("salmon", "dodgerblue", "lightgrey")) + labs(y = "", size = "nCells", fill = "") + add_guide
  
}



readFACS <- function(path){
  
  tube_raw <- fread(path)
  tubename <- substr(path, 16, 20)
  colnames(tube_raw) <- gsub("\\+", "pos", colnames(tube_raw))
  colnames(tube_raw) <- gsub("\\-", "neg", colnames(tube_raw))
  colnames(tube_raw) <- make.names(colnames(tube_raw)) %>% make.unique() 
  
  if("FM.ID" %in% colnames(tube_raw)){
    tube_raw <- tube_raw %>% dplyr::rename("V1" = "FM.ID")
  }
  
  ## Remove all other than mel naive or not?
  tube_raw <- tube_raw %>% filter(Cohort == 1)
  
  tube_raw <- tube_raw[,-c(2:4)]
  tube_raw <- tube_raw %>% 
    mutate(timepoint = as.factor(Sample)) %>% 
    mutate(timepoint = plyr::revalue(timepoint, c("1" = "0m", "2" = "1m", "3" = "3m", "4" = "6m"))) %>% 
    dplyr::rename(name = V1) %>% 
    select(-Sample) %>% mutate(key = paste(name, timepoint)) %>% select(-name,-timepoint)
  
}



vdjToGliph <- function(vdj_df){
  
  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj 
  
  df <- vdj_df %>% 
    dplyr::select("cdr3aa", "v", "j", "name", "freq") %>% 
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = freq) 
  
  return(df)
  
}


facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))
add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set1"))


extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}



getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}


plotLatentUmap <- function(viz_df, cluster){

  viz_df_temp <- data.frame(viz_df, "seurat_cluster" = cluster)
  nClusters   <- unique(cluster) %>% length

  ## Visualise
  umap_mean <- data.frame(aggregate(umap_1 ~ seurat_cluster, viz_df_temp, median), umap_2 = aggregate(umap_2 ~ seurat_cluster, viz_df_temp, median)[,2])


  ## Plot UMAPs with TCRGP predictions highlighted
  ggplot() +
    geom_point(data = viz_df_temp, aes(x = umap_1, y = umap_2, color = seurat_cluster), size = 0.8) +

    # stat_ellipse(data = viz_df_temp, geom = "polygon", aes(x = umap_1, y = umap_2, color = seurat_cluster, fill = seurat_cluster), alpha = 0.1, lty = "dotted") +
    ggrepel::geom_label_repel(data = umap_mean, aes(x = umap_1, y = umap_2, color = seurat_cluster, label = seurat_cluster), size = 5, color = "black") +

    theme_void() + theme(legend.position = "none") +
    scale_color_manual(values = getPalette(nClusters)) +
    scale_fill_manual(values = getPalette(nClusters)) + labs()


}


extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){

  strsplit(str1, "[//]")[[1]][4]

}


plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){

  ## Plot univariate violin plots with filter thresholds

  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable

  viz_df_temp <- viz_df %>% select(var_to_plot)

  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table

  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +

    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +

    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +

    labs(x = "", title = var_to_plot) + theme(legend.position = "none")

}






# getClusterPhenotypes <- function(clusters){
#
#   clusters <- plyr::revalue(clusters, replace = c(
#
#     "0"  = "00_NK_NK",
#     "1"  = "01_CD4_naive",
#     "2"  = "02_CD4_EM/Th1-like",
#     "3"  = "03_CD8_EM",
#     "4"  = "04_CD8_TRM",
#     "5"  = "05_CD4CD8_naive",
#     "6"  = "06_CD4CD8_naive",
#     "7"  = "07_CD4_TRM/Th2",
#     "8"  = "08_other_MAIT",
#     "9"  = "09_B",
#     "10" = "10_CD8_effector/exhausted",
#     "11" = "11_CD4_Treg",
#     "12" = "12_B_B/T",
#     "13" = "13_B_?",
#     "14" = "14_NK_?",
#     "15" = "15_B_?",
#     "16" = "16_CD8_?",
#     "17" = "17_monocyte_promonocyte",
#     "18" = "18_monocyte_monocyte",
#     "19" = "19_NK_?",
#     "20" = "20_other_HSC-like?",
#     "21" = "21_other_pDC",
#     "22" = "22_other_junk",
#     "23" = "23_other_plasma"))
#
#
#   return(clusters)
#
#
# }



getClusterPhenotypes <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "0 NK CD56dim",
    "1"  = "1 CD4 naive/CM",
    "2"  = "2 CD4 EM/Th1-like",
    "3"  = "3 CD8 EM",
    "4"  = "4 CD8 effector",
    "5"  = "5 CD8 naive",
    "6"  = "6 T-cell activated",    ## CD8+, cytotoxic markers, also CD4+ and clusters nearer to CD4+ cells. Huge cluster heterogeneity, possible three phenotypes in one cluster. Also exhaustion markers
    "7"  = "7 CD8 CM", ## CD8+, but no cytotoxic markers, clusters nearer to CD4 cells
    "8"  = "8 T-cell MAIT",
    "9"  = "9 B naive/resting",
    "10" = "10 NK adaptive",
    "11" = "11 CD4 treg",
    "12" = "12 B B/T",
    "13" = "13 B class-switched memory",
    "14" = "14 NK CD56bright",
    "15" = "15 B naive/resting CMV-specific", # female-specific
    "16" = "16 T-cell exhausted cycling",
    "17" = "17 monocyte CD16-",
    "18" = "18 monocyte CD16+",
    "19" = "19 NK cycling",
    "20" = "20 low quality",
    "21" = "21 monocyte pDC",
    "22" = "22 megakaryocytes",
    "23" = "23 B plasma"))

  return(clusters)

}


getClusterPhenotypesImmune <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK",
                                                    
                                                    "1"  = "1",
                                                    "2"  = "2",
                                                    "3"  = "3 CD4+ Tem",
                                                    "4"  = "4 NK",
                                                    "5"  = "5 Tregs",
                                                    "6"  = "6 CD8+ Tcm",
                                                    "7"  = "7 Class-switched memory B-cells",
                                                    "8"  = "8 Monocytes",
                                                    "9"  = "9 CD8+ Tcm",
                                                    "10" = "10 CD4+ Tem",
                                                    
                                                    "11" = "11 CD4+ Tem ",
                                                    "12" = "12 CD8+ Tcm",
                                                    "13" = "13",
                                                    "14" = "14",
                                                    "15" = "15",
                                                    "16" = "16 CD8+ Tcm",
                                                    "17" = "17",
                                                    "18" = "18 Monocytes",
                                                    "19" = "19 B-cell",
                                                    "20" = "20 CD4+ Tem",

                                                    "21" = "21 CD8+ Tcm",
                                                    "22" = "22 CD8+ Tcm",
                                                    "23" = "23 Monocytes",
                                                    "24" = "24 Tregs",
                                                    "25" = "25",
                                                    "26" = "26 Plasma cells",
                                                    "27" = "27 NK",
                                                    "28" = "28 CD8+ Tcm",
                                                    "29" = "29 CD8+ Tcm",
                                                    "30" = "30 Macrophages",
                                                    
                                                    "31" = "31 B-cell", 
                                                    "32" = "32 Epithelial cells", 
                                                    "33" = "33 CD8+ Tcm", 
                                                    "34" = "34 Memory B-cells", 
                                                    "35" = "35", 
                                                    "36" = "36 Monocytes", 
                                                    "37" = "37 CD4+ Tem", 
                                                    "38" = "38 Monocytes", 
                                                    "39" = "39 Macrophages", 
                                                    "40" = "40", 
                                                    
                                                    "41" = "41 CD8+ Tcm", 
                                                    "42" = "42 CD4+ Tem", 
                                                    "43" = "43", 
                                                    "44" = "44 Fibroblasts", 
                                                    "45" = "45 Monocytes",     
                                                    "46" = "46 Melanocytes", 
                                                    "47" = "47 mv Endothelial cells", 
                                                    "48" = "48 MEP", 
                                                    "49" = "49", 
                                                    "50" = "50"
                                                    
  ))
  
  return(clusters)
  
}


getClusterPhenotypesImmuneRP <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 CD4+ Tem",
                                                    
                                                    "1"  = "1 NK cells",
                                                    "2"  = "2 NK cells",
                                                    "3"  = "3 CD4+ Tem",
                                                    "4"  = "4 CD4+ T-cells",
                                                    "5"  = "5 Class-switched memory B-cells",
                                                    "6"  = "6 CD8+ Tcm",
                                                    "7"  = "7 CD8+ Tcm",
                                                    "8"  = "8 CD8+ Tcm",
                                                    "9"  = "9 Tregs",
                                                    "10" = "10 Monocytes",
                                                    
                                                    "11" = "11 CD8+ Tcm ",
                                                    "12" = "12 Monocytes",
                                                    "13" = "13 CD4+ T-cells",
                                                    "14" = "14 CD8+ Tcm",
                                                    "15" = "15 CD4+ Tem",
                                                    "16" = "16 CD4+ Tem",
                                                    "17" = "17 CD4+ Tem",
                                                    "18" = "18 CD8+ Tcm",
                                                    "19" = "19 Tregs",
                                                    "20" = "20 CD8+ Tcm ",
                                                    
                                                    "21" = "21 Plasma cells",
                                                    "22" = "22 CD8+ Tcm",
                                                    "23" = "23 Monocytes",
                                                    "24" = "24 B cells",
                                                    "25" = "25 Epithelial cells",
                                                    "26" = "26 Macrophages",
                                                    "27" = "27 CD4+ Tem",
                                                    "28" = "28 Fibroblasts",
                                                    "29" = "29 DC",
                                                    "30" = "30 Melanocytes",
                                                    
                                                    "31" = "31 MEP", 
                                                    "32" = "32 Class-switched memory B-cells"))
                                                    
  
  return(clusters)
  
}




getClusterPhenotypesTNK <- function(clusters){
  
  # clusters <- substr(clusters, 1, 2) %>% as.numeric() %>% as.character()
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0" = "0 NK cells",
    "1" = "1 CD4 EM/Th1-like", ## EM: ; Th1: 
    "2" = "2 CD4 CM",
    "3" = "3 CD8 EM",
    "4" = "4 CD8 effector",
    "5" = "5 T naive", ## 
    "6" = "6 CD4 activated", ## CD69, NRA42
    "7" = "7 MAIT",
    "8" = "8 CD4 treg",
    "9" = "9 CD8 effector/exhausted",
    "10" = "10 CD8 NK-like", ## FCGR3A. XCL1, XCL2, TBX21, HIF1
    "11" = "11 CD8 CM",
    "12" = "12 CD8",
    "13" = "13 CD4 CM",
    "14" = "14 NK NKT", 
    "15" = "15 NK cycling")) ## MKI67
  
  return(clusters)
  
}




getClusterCoarsePhenotype <- function(clusters){

  # func_temp <- function(str1){strsplit(clusters, "[ ]")[[1]][2]}
  # vec <- lapply(clusters, func_temp)
  # vec <- do.call(what = "c", vec)
  # return(vec)

  clusters <- substr(clusters, 1, 2) %>% as.numeric() %>% as.character()

  clusters <- plyr::revalue(clusters, replace = c(

    "0"  = "NK",
    "1"  = "CD4",
    "2"  = "CD4",
    "3"  = "CD8",
    "4"  = "CD8",
    "5"  = "CD8",
    "6"  = "CDCD8",
    "7"  = "CD8",
    "8"  = "other",
    "9"  = "B",
    "10" = "NK",
    "11" = "CD4",
    "12" = "B",
    "13" = "B",
    "14" = "NK",
    "15" = "B",
    "16" = "CDCD8",
    "17" = "other",
    "18" = "other",
    "19" = "NK",
    "20" = "CD4CD8",
    "21" = "other",
    "22" = "other",
    "23" = "other"))

  return(clusters)

}





## Helper function to break TCRb names into clinical data
breakName <- function(name_temp){
  
  ## Breaks the filename into relevant information
  
  # "2121_PR_PB_LAG3_R_0m_MNC.txt"
  
  name      <- substr(name_temp, 1, 4)
  response  <- substr(name_temp, 6, 7)
  type      <- substr(name_temp, 9, 10)
  regimen   <- substr(name_temp, 12, 15)
  overall   <- substr(name_temp, 17, 17)
  timepoint <- substr(name_temp, 19, 20)
  celltype  <- substr(name_temp, 22, 24)
  io_stat   <- substr(name_temp, 26, 27)
  
  
  for(i in 1:length(type)){
    if(type[i] == "PB"){type[i] = "Blood"}
    if(type[i] == "TX"){type[i] = "Tumor"}
  }
  
  for(i in 1:length(regimen)){
    if(regimen[i] == "CTLA"){regimen[i] = "antiCTLA4"}
    if(regimen[i] == "iPD1"){regimen[i] = "antiPD1"}
    if(regimen[i] == "NiIp"){regimen[i] = "antiPD1+antiCTLA4"}
    if(regimen[i] == "IpNi"){regimen[i] = "antiCTLA4+antiPD1"}
    if(regimen[i] == "LAG3"){regimen[i] = "antiPD1+antiLAG3"}
    if(regimen[i] == "NANA"){regimen[i] = "TreatmentNA"}
  }
  
  for(i in 1:length(overall)){
    if(overall[i] == "U"){overall[i] = "NA"}
  }
  
  for(i in 1:length(io_stat)){
    if(io_stat[i] == "pr"){io_stat[i] = "Prior.IO"}
    if(io_stat[i] == "nv"){io_stat[i] = "IO.naive"}
    
  }
  
  return(data.frame(name, response, type, regimen, overall, timepoint, celltype, io_stat, stringsAsFactors = F))
  
}




fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}





extractCoarsePhenotype <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }
  
  return(p)
  
}




removeTCRabData <- function(seurat_object){
  
  seurat_object@meta.data$cdr3s_nt <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$cdr3s_nt               <- NULL
  seurat_object@meta.data$tra_cdr3s_nt <- NULL
  seurat_object@meta.data$trb_cdr3s_nt <- NULL
  seurat_object@meta.data$cdr3s_aa <- NULL
  seurat_object@meta.data$tra_cdr3s_aa <- NULL
  seurat_object@meta.data$trb_cdr3s_aa <- NULL
  seurat_object@meta.data$clonotype_id <- NULL
  seurat_object@meta.data$frequency <- NULL
  seurat_object@meta.data$proportion <- NULL
  seurat_object@meta.data$cdr3s_aa.1 <- NULL
  seurat_object@meta.data$cdr3s_nt.1 <- NULL
  seurat_object@meta.data$chain_tra <- NULL
  seurat_object@meta.data$v_tra <- NULL
  seurat_object@meta.data$d_tra <- NULL
  seurat_object@meta.data$j_tra <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_tra <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_tra <- NULL
  seurat_object@meta.data$v_freq_tra <- NULL
  seurat_object@meta.data$d_freq_tra <- NULL
  seurat_object@meta.data$j_freq_tra <- NULL
  seurat_object@meta.data$chain_trb             <- NULL
  seurat_object@meta.data$v_trb <- NULL
  seurat_object@meta.data$d_trb <- NULL
  seurat_object@meta.data$j_trb <- NULL
  seurat_object@meta.data$cdr3s_nt_freq_trb <- NULL
  seurat_object@meta.data$cdr3s_aa_freq_trb <- NULL
  seurat_object@meta.data$v_freq_trb <- NULL
  seurat_object@meta.data$d_freq_trb <- NULL
  seurat_object@meta.data$j_freq_trb <- NULL
  seurat_object@meta.data$tcr_type <- NULL
  seurat_object@meta.data$patient <- NULL
  seurat_object@meta.data$new_clonotypes_id  <- NULL
  
  return(seurat_object)

}









## Olink categories
apoptosis_cell_killing <- c("CASP8", "CD40L", "FasL", "Gal9", "GZMA", "GZMB", "GZMH", "MMP7",
                            "TRAIL", "TWEAK", "TNFRSF12A", "TNFRSF21")
chemotaxis <- c("MCP1", "CCL3" , "CCL4"  , "MCP3" , "MCP2" , "MCP4" , "CCL17" , "CCL19" , "CCL20" , "CCL23" , "CXCL1" , "CXCL5"  , "CXCL9"  , "CXCL10" , "CXCL11")

metabolism_autophagy <- c("ADA", "CAIX", "HO1")

promote_tumor_immunity <- c("CD27" ,"CD40" ,"CD40L","CD70" ,"CD83" ,"CXCL9"  ,"CXCL10" ,"CXCL11" ,"CXCL13" ,"CRTAM" ,"CX3CL1"  ,"ICOSLG" ,"IFNG" ,"IL1" ,"IL2" ,"IL6" ,"IL7" ,"IL12RB1" ,"IL18" ,"NCR1" ,"CD244" ,"KLRD1" ,"CD4" ,"CD5" ,"CD8A" ,"CD28" ,"TNF" ,"TNFSF14"  ,"TNFRSF4"  ,"TNFRSF9" ,"IL15")

suppress_tumor_immunity <- c("MCP4","CCL17","CCL19","CCL20","CXCL1","CXCL5","CXCL11","CXCL13","Gal1","Gal9","ARG1","IL4","IL5","IL6","IL8","IL10","IL13","IL18","IL33","LAP","LAMP3","CSF1","MMP12","MMP7","MICA/B","PDL1","PDL2","PDCD1","CXCL12","CD4","CD5","TNF","MUC16","KIR3DL1","IL15","LAG3","IL15","MUC16","IL4","IL5","IL6","IL8","IL10","IL13","IL18","IL33")

vascular_tissue_remodelling <- c("ADGRG1","ANGPT1","TIE2","ANGPT2","CAIX","MCP1","MCP4","CCL23","CXCL1","CXCL5","CXCL9","CXCL10","CXCL11","DCN","FGF2","Gal1","Gal9","HGF","IL1","IL8","MMP12","NOS3","PGF","PDGF","PTN","EGF","CXCL12","TNF","TWEAK","TNFRSF12A","VEGFA","VEGFR2")
