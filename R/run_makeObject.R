
message("Create seurat-object from individual seurat-objects ...")

getSeuratName  <- function(x){substr(x, 16, nchar(x))}

# getSeuratName(folders)
folders        <- list.dirs("data/scRNAseq/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
# lag3_seurat     <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = getSeuratName(folders))
lag3_seurat     <- scrnaseq_files[[1]]


## Basic QC
dir.create("results/qc/", showWarnings = F)
dir.create("results/qc/before_1/", showWarnings = F)
dir.create("results/qc/after_1/", showWarnings = F)

lag3_seurat  <- PercentageFeatureSet(lag3_seurat, pattern = "^MT-", col.name = "percent.mt")
lag3_seurat  <- PercentageFeatureSet(lag3_seurat, pattern = "^RP", col.name = "percent.ribo")
lag3_seurat  <- PercentageFeatureSet(lag3_seurat, features = cycle.genes, col.name = "percent.cycle")
lag3_seurat@meta.data$barcode   <- colnames(lag3_seurat)

lag3_seurat %>% plotQC(folder = "results/qc/before_1/")
lag3_seurat <- lag3_seurat %>% getQC()


## Get SingleR predictions; omit predictions from cell types rare than 10 cells
lag3_seurat             <- lag3_seurat %>% getSingler()
relevant_hpca_clusters <- lag3_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- lag3_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_blueprint_pred)

lag3_seurat$singler_hpca_pred      <- ifelse(lag3_seurat$singler_hpca_pred %in% relevant_hpca_clusters, lag3_seurat$singler_hpca_pred, "rare")
lag3_seurat$singler_blueprint_pred <- ifelse(lag3_seurat$singler_blueprint_pred %in% relevant_blue_clusters, lag3_seurat$singler_blueprint_pred, "rare")

## Get doublets
lag3_seurat <- lag3_seurat %>% getDoublets()
lag3_seurat <- subset(lag3_seurat, hybrid_doublet_score < 1.8)

## Get Seurat
clonality_genes <- getClonalityGenes(lag3_seurat)
unwanted_genes  <- getUnwantedGenes(lag3_seurat)

lag3_seurat <- lag3_seurat %>% preprocessSeurat(cells.to.use = colnames(lag3_seurat))
lag3_seurat <- lag3_seurat %>% getClustering()

# ## Get scVI
dir.create("results/scvi/", showWarnings = F)
dir.create("results/scvi/input_files/", showWarnings = F)
lag3_seurat$orig.ident <- gsub("\\/", "\\_", lag3_seurat$orig.ident)
lag3_seurat %>% getScviInput(folder = "results/scvi/input_files/")
#
## Get scVI results
latents    <- fread("results/scvi/results/lag3_latent.csv")
lag3_seurat <- lag3_seurat %>% putLatentsSeurat(latent = latents)
lag3_seurat <- lag3_seurat %>% getLatentClustering() %>% fixSeurat()

## Decide on clustering
lag3_seurat %>% plotClustering()
ggsave("results/qc/after_1/scatter_clustering.png", width = 5, height = 4)
saveRDS(lag3_seurat, "results/lag3_seurat.rds")
