

message("Create seurat-object from individual seurat-objects ...")

## Select non-tumor samples
folders        <- list.dirs("data/scRNAseq/", recursive = T) %>% grep(pattern = "filtered_gene_bc_matrices|filtered_feature_bc_matrix", value = T) %>% grep(pattern = "tumor", invert = T, value = T)
scrnaseq_files <- lapply(folders, function(x){print(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 0, min.features = 0)})
lag3_seurat    <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = lapply(folders, extractSeuratName) %>% unlist, min.cells = 0, min.genes = 0, do.normalize = F, do.scale = F, do.center = F)

## Basic QC
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

lag3_seurat  <- PercentageFeatureSet(lag3_seurat, pattern = "^MT-", col.name = "percent.mt")
lag3_seurat  <- PercentageFeatureSet(lag3_seurat, pattern = "^RP", col.name = "percent.ribo")
lag3_seurat  <- PercentageFeatureSet(lag3_seurat, features = cycle.genes, col.name = "percent.cycle")
lag3_seurat@meta.data$barcode <- colnames(lag3_seurat)

saveRDS(lag3_seurat, "results/lag3_seurat.rds")
