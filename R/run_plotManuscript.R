
dir.create("results/manuscript", showWarnings = F)
dir.create("results/manuscript/figure1/", showWarnings = F)

lag3_seurat <- readRDS("results/lag3_seurat_latest.rds")

## Figure 1
Idents(lag3_seurat)    <- Idents(lag3_seurat) %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()
lag3_seurat$cluster    <- Idents(lag3_seurat) 
lag3_seurat$patient    <- substr(colnames(lag3_seurat), 1, 4)
lag3_seurat$cmv_status <- ifelse(lag3_seurat$patient %in% c("2348", "2350"), "neg", "pos") 
lag3_seurat$hla_a2     <- ifelse(lag3_seurat$patient %in% c("2083", "2342"), "neg", "pos") 

all_markers <- FindAllMarkers(lag3_seurat, test = "t", max.cells.per.ident = 1e3)
all_markers <- all_markers %>% filter(avg_log2FC > 0)
fwrite(all_markers, "results/manuscript/all_markers.txt", sep = "\t", quote = F, row.names = F)


DimPlot(lag3_seurat, label = T, repel = T) + theme_bw(base_size = 17) + theme(legend.position = "none") + scale_color_manual(values = getPalette(24)) + labs(x = "UMAP from latents 1", y = "UMAP from latents 2")
ggsave("results/manuscript/figure1/latent_umap.png", width = 6, height = 5)

DimPlot(lag3_seurat, label = F, repel = T) + theme_void(base_size = 17) + theme(legend.position = "none") + scale_color_manual(values = getPalette3(24)) + labs(x = "UMAP from latents 1", y = "UMAP from latents 2")
ggsave("results/manuscript/figure1/latent_umap2.png", width = 6, height = 5)

lag3_seurat$cluster <- Idents(lag3_seurat)
patient_df <- lag3_seurat@meta.data %>% group_by(orig.ident, patient, timepoint) %>% summarise(n = n()) %>% dplyr::select(-n) %>% mutate(prepost = ifelse(timepoint == 0, "pre", "post"))

lag3_seurat@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df, by = "orig.ident") %>% 
  ggplot(aes(timepoint, prop)) + geom_boxplot(outlier.shape = NA) + facet_grid(~cluster) + ggpubr::stat_compare_means(label = "p", method = "t")

lag3_seurat@meta.data %>% 
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df, by = "orig.ident") %>% 
  ggplot(aes(prepost, prop)) + geom_boxplot(outlier.shape = NA) + facet_grid(~cluster) + ggpubr::stat_compare_means(label = "p") + geom_jitter(size=0.5)


lag3_seurat@meta.data %>% 
  group_by(cluster, patient) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% #left_join(patient_df, by = "orig.ident") %>% 
  ggplot(aes(cluster, prop, fill = patient)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(10)) + ggpubr::rotate_x_text(0) + labs(x = "") + coord_flip() + 
  theme_classic(base_size = 17) +
  theme(legend.position = "top")
ggsave("results/manuscript/figure1/bar_prop.pdf", width = 7, height = 7)


viz_df <- lag3_seurat@meta.data %>% mutate(cluster = Idents(lag3_seurat))

plotLolliplot(viz_df, timepoint_temp = "1") + theme_classic(base_size = 17) + xlim(values = c(-5,5)) + labs(x = "log2(fold change)") + scale_fill_manual(values = c("dodgerblue", "lightgrey", "salmon"), guide=FALSE)
ggsave("results/manuscript/figure1/lolliplot_baseline.png", width = 8, height = 5)

p <- DotPlot(lag3_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/manuscript/figure1/dotplot_big_markers.png", width = 14, height = 6)

p <- VlnPlot(lag3_seurat, features = c("LAG3", "PDCD1"), pt.size = 1e-1000, cols = getPalette(24)) + labs(x = "")
ggsave(plot = p, "results/manuscript/figure1/vln_lag3.png", width = 12, height = 4)

p <- FeaturePlot(lag3_seurat, features = c(big_markers, "PRF1", "LAG3", "PDCD1", "MKI67"), ncol = 6, order = T, cols = c("gray90", "red3"), combine = TRUE) & NoAxes() 
ggsave(plot = p, "results/manuscript/figure1/features_lag3.png", width = 16, height = 6)

p <- FeaturePlot(lag3_seurat, features = c("LAG3", "PDCD1", "CTLA4", "PRF1", "GZMB", "MKI67"), ncol = 2, order = F, cols = c("gray90", "red3"), combine = TRUE) & NoAxes() 
ggsave(plot = p, "results/manuscript/figure1/features_lag3_manu.png", width = 5, height = 6)

df <- lag3_seurat@meta.data %>% 
  bind_cols(data.frame(lag3 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "LAG3", ])) %>% 
  bind_cols(data.frame(pdcd1 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "PDCD1", ]))

lag3pos_df <- df  %>% filter(timepoint == "0") %>% group_by(cluster, pos = lag3 > 0) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(pos == T)
pd1pos_df  <- df  %>% filter(timepoint == "0") %>% group_by(cluster, pos = pdcd1 > 0) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(pos == T)

lag3pos_df %>% left_join(pd1pos_df, by = c("cluster")) %>% 
  ggplot(aes(prop.x,prop.y,label=cluster)) + geom_point() + ggrepel::geom_text_repel()

lag3pos_df %>% ggplot(aes(reorder(cluster, prop), prop)) + geom_bar(stat = "identity") + coord_flip()


df <- lag3_seurat@meta.data %>% bind_cols(data.frame(lag3 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "LAG3", ])) %>% filter(timepoint == 0)
median_df <- df %>% filter(lag3 > 0) %>% group_by(cluster) %>% summarise(median = median(lag3), sd = sd(lag3)) %>% filter(!cluster %in% c("12 B B/T"))
cl.to.include <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n > 2e2) %>% filter(!cluster %in% c("12 B B/T"))

df %>% left_join(median_df) %>% filter(cluster %in% cl.to.include$cluster) %>% 
  ggplot(aes(reorder(cluster, median), log2(lag3), fill = median)) + geom_violin(draw_quantiles = 0.5, adjust = 1) + theme_classic(base_size = 17) + 
  # scale_y_log10() + 
  # geom_jitter(size = 0.5, alpha = 0.5) +
  labs(x = "", y = expression(paste("scaled ", italic("LAG3"), " expression"))) + scale_fill_distiller(palette ="RdBu", direction = -1) + coord_flip() + theme(legend.position = "none")
ggsave("results/manuscript/figure1/vln_lag3.pdf", width = 7, height = 7)


df %>% left_join(median_df) %>% filter(cluster %in% cl.to.include$cluster) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(cluster = ifelse(cluster == "15 B naive/resting CMV-specific", "15 B resting CMV-specific", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "13 B class-switched memory", "13 B class-switched", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "16 T-cell exhausted cycling", "16 T-cell exhausted cycling", cluster)) %>% 
  
  ggplot(aes(y = reorder(cluster, median), x = log2(lag3), fill = median)) + 
  ggridges::geom_density_ridges(scale = 3, alpha = 0.8, quantile_lines = TRUE, quantiles = 2, jittered_points = F,
                                position = ggridges::position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.5) +
  theme_classic(base_size = 19) + 
  labs(y = "", x = expression(paste("log2 scaled ", italic("LAG3"), " expression"))) + scale_fill_distiller(palette ="RdBu", direction = -1) + theme(legend.position = "none")
ggsave("results/manuscript/figure1/ridge_lag32.pdf", width = 7, height = 6)


df <- lag3_seurat@meta.data %>% bind_cols(data.frame(lag3 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "PDCD1", ])) %>% filter(timepoint == 0)
median_df <- df %>% filter(lag3 > 0) %>% group_by(cluster) %>% summarise(median = median(lag3))
cl.to.include <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n > 2e2) %>% filter(!cluster %in% c("12 B B/T"))
df %>% left_join(median_df) %>% filter(cluster %in% cl.to.include$cluster) %>% 
  ggplot(aes(reorder(cluster, median), lag3, fill = median)) + geom_violin(draw_quantiles = 0.5, adjust = 1) + theme_classic(base_size = 17) + scale_y_log10() + labs(x = "", y = "scaled PDCD1 expression") + scale_fill_distiller(palette ="RdBu", direction = -1) + coord_flip() + theme(legend.position = "top")
# ggpubr::rotate_x_text(45)
ggsave("results/manuscript/figure1/vln_pdcd1.pdf", width = 7, height = 7)

df %>% left_join(median_df) %>% filter(cluster %in% cl.to.include$cluster) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(cluster = ifelse(cluster == "15 B naive/resting CMV-specific", "15 B resting CMV-specific", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "13 B class-switched memory", "13 B class-switched", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "16 T-cell exhausted cycling", "16 T-cell exhausted cycling", cluster)) %>% 
  
#  mutate(cluster = ifelse(cluster == "15 B naive/resting CMV-specific", "15 B-naive/resting\nCMV-specific", cluster)) %>% 
  ggplot(aes(y = reorder(cluster, median), x = log2(lag3), fill = median)) + 
  ggridges::geom_density_ridges(scale = 3, alpha = 0.8, quantile_lines = TRUE, quantiles = 2, jittered_points = F,
                                position = ggridges::position_points_jitter(width = 0.05, height = 0),
                                point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.5) +
  theme_classic(base_size = 17) + 
  labs(y = "", x = expression(paste("log2 scaled ", italic("PDCD1"), " expression"))) + scale_fill_distiller(palette ="RdBu", direction = -1) + theme(legend.position = "none")
ggsave("results/manuscript/figure1/ridge_pdcd12.pdf", width = 7, height = 6)




## scatter of lag3 and pdcd1
df <- lag3_seurat@meta.data %>% bind_cols(data.frame(lag3 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "LAG3", ]),
                                                    pdcd1 = lag3_seurat@assays$RNA@scale.data[rownames(lag3_seurat@assays$RNA@scale.data) == "PDCD1", ]) %>% filter(timepoint == 0)
median_lag3_df  <- df %>% 
  # filter(lag3 > 0) %>% 
  group_by(cluster) %>% summarise(median_lag3 = mean(lag3), sd_lag3 = sd(lag3))
median_pdcd1_df <- df %>% 
  # filter(pdcd1 > 0) %>% 
  group_by(cluster) %>% summarise(median_pdcd1 = mean(pdcd1), sd_pdcd1 = sd(pdcd1))

cl.to.include <- df %>% group_by(cluster) %>% filter(timepoint == "0") %>% summarise(n = n()) %>% filter(n > 2e2) %>% filter(!cluster %in% c("12 B B/T"))

median_lag3_df %>% left_join(median_pdcd1_df) %>% filter(cluster %in% cl.to.include$cluster) %>% 
  
  filter(median_lag3 > 0 | median_pdcd1 > 0) %>% 
  # filter(median_pdcd1 > 0) %>% 
  ggplot(aes(x = median_lag3, y = median_pdcd1, label = cluster,
             xmin = median_lag3 - sd_lag3, xmax = median_lag3 + sd_lag3,
             ymin = median_pdcd1 - sd_pdcd1, ymax = median_pdcd1 + sd_pdcd1,
             )) + geom_errorbarh(color="gray90") + geom_pointrange(color="gray90") +
  geom_point(color="black") + ggrepel::geom_text_repel() 

lag3_seurat$RNA_snn_res.0.6
lag3_seurat2 <- lag3_seurat
Idents(lag3_seurat2) <- lag3_seurat2$RNA_snn_res.0.6
p <- FeaturePlot(lag3_seurat2, features = c("LAG3", "PDCD1"), pt.size = 1e-4, order = T, label = T, repel = T, cols = c("gray90", "red3")) & theme_void(base_size = 17) & theme(legend.position = "top") & theme(plot.title = element_text(face = "bold.italic"))
ggsave(plot = p, "results/manuscript/feature_lag3pdcd1.png", width = 8, height = 5)

# scale_fill_gradient2(midpoint = 3.5, direction = -1)

lag3_seurat@meta.data %>% bind_cols(data.frame(lag3 = lag3_seurat@assays$RNA@data[rownames(lag3_seurat@assays$RNA@data) == "LAG3", ])) %>% 
  ggplot(aes(cluster, lag3)) + geom_violin(draw_quantiles = 0.5) + coord_flip()


DotPlot(diet_rp_seurat, features = rev(unique(c("CD3E", "LAG3", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH"))), cols = "RdYlBu") + 
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")

## Supplementary figure 1
DimPlot(lag3_seurat, label = F, repel = T, group.by = "orig.ident", split.by = "orig.ident", ncol = 3) + 
  theme_bw(base_size = 17) + theme(legend.position = "none") + scale_color_manual(values = getPalette(18)) + labs(x = "UMAP from latents 1", y = "UMAP from latents 2")
ggsave("results/manuscript/figure1/latent_umap_orig_ident.png", width = 8, height = 18)

p <- FeaturePlot(lag3_seurat, features = signature_features, cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0, ncol = 8, pt.size = 0.5) + theme_bw() + labs(x = "UMAP 1 from latents", y = "UMAP 2 from latents") #+ labs(title = paste0("Zheng", signature_feature))
ggsave(plot = p, "results/manuscript/figure1/umap_meta.png", width = 24, height = 16)

DotPlot(lag3_seurat, features = rev(unique(c("NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRB1", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH"))), cols = "RdYlBu") + 
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave("results/manuscript/figure1/dotplot_dufva.pdf", width = 8, height = 5)




DimPlot(diet_rp_seurat, label = T, repel = T) + 
  theme_bw(base_size = 17) + theme(legend.position = "none") + scale_color_manual(values = getPalette(33)) + labs(x = "UMAP from latents 1", y = "UMAP from latents 2")

DotPlot(diet_rp_seurat, features = rev(unique(c("CD3E", "LAG3", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH"))), cols = "RdYlBu") + 
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")


DotPlot(diet_rp_seurat, features = rev(unique(c("NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH"))), cols = "RdYlBu") + 
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave("results/manuscript/figure1/dotplot_dufva.pdf", width = 8, height = 5)

## CMV plot
set.seed(111)
min_cells = 40e3
cmv_df <- lag3_seurat@meta.data %>% mutate(cluster = Idents(lag3_seurat)) %>% 
  group_by(cmv_status) %>% sample_n(min_cells) %>% 
  group_by(cluster, cmv_status) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) 

cluster_order <- cmv_df %>% filter(cmv_status == "pos") %>% arrange(prop) %>% pull(cluster) %>% as.character()

cmv_df %>% ungroup() %>% 
  #mutate(cluster = factor(cluster, levels = c("23 B plasma", cluster_order))) %>% 
  mutate(cluster = factor(cluster, levels = c(cluster_order))) %>% 
  filter(cluster != "20 low quality") %>% 
  
  # mutate(cluster = reorder(cluster, cluster_order)) %>% 
  ggplot(aes(cluster, prop, fill = cmv_status)) + geom_bar(stat = "identity") + coord_flip() + geom_hline(yintercept = 0.5, linetype = "dotted") +
  # scale_fill_manual(values = getPalette(4)) + 
  scale_fill_manual(values = c("darkblue", "darkred")) + 
  labs(x = "", fill = "CMV") + theme_classic(base_size = 17) + theme(legend.position = "top")
ggsave("results/manuscript/figure1/bar_cmv.pdf", width = 6, height = 6)




  ## -------------------- Figure 2
dir.create("results/manuscript/figure2/", showWarnings = F)


lag3_seurat@meta.data %>% mutate(cluster = Idents(lag3_seurat)) %>% 
  group_by(timepoint,overall,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% 
  
  ggplot(aes(timepoint, prop, color = overall)) + geom_path(lwd = 1.5) + scale_color_manual(values = c("lightgrey", "salmon")) +
  facet_wrap(~cluster, scales = "free_y") + facets_nice
ggsave("results/manuscript/figure2/path_all.pdf", width = 12, height = 8)


#### Analyze
DEG_cluster_df   <- read.delim("results/effectOnClusters/deg_cluster.txt")   %>% mutate(celltype = getClusterCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))
DEG_cluster_r_df <- read.delim("results/effectOnClusters/deg_r_cluster.txt") %>% mutate(celltype = getClusterCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))
DEG_cluster_n_df <- read.delim("results/effectOnClusters/deg_n_cluster.txt") %>% mutate(celltype = getClusterCoarsePhenotype(as.character(cluster)), cluster = getClusterPhenotypes(as.character(cluster)))


## Volcano plot on genes

## 2v1 total
volcano_df <- DEG_cluster_df %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_2v1.png", width = 7, height = 6)

## 3v1 total
volcano_df <- DEG_cluster_df %>% filter(timepoint == "3v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_3v1.png", width = 7, height = 6)



## 2v1 r
volcano_df <- DEG_cluster_r_df %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_2v1_r.png", width = 7, height = 6)

## 3v1 r
volcano_df <- DEG_cluster_r_df %>% filter(timepoint == "3v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_3v1_r.png", width = 7, height = 6)



## 2v1 n
volcano_df <- DEG_cluster_n_df %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_2v1_n.png", width = 7, height = 6)

## 3v1 n
volcano_df <- DEG_cluster_n_df %>% filter(timepoint == "3v1") %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/figure2/volcano_3v1_n.png", width = 7, height = 6)



## 3v1 r, 10 nk adaptive
volcano_df <- DEG_cluster_r_df %>% filter(timepoint == "3v1" | timepoint == "3v2") %>% mutate(cluster = reorderClusters(cluster)) %>% filter(cluster == "10 NK adaptive")
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-1,1)) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice +
#  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long | -log10(p_val) > 30), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
ggsave("results/manuscript/figure2/volcano_3v1_nk.png", width = 5, height = 4)



clonality_genes <- getClonalityGenes(lag3_seurat)
unwanted_genes  <- getUnwantedGenes(lag3_seurat)

volcano_df <- DEG_cluster_df %>% filter(timepoint == "3v1") %>% mutate(cluster = reorderClusters(cluster)) %>% filter(cluster == "10 NK adaptive")
volcano_df$avg_logFC <- -volcano_df$avg_logFC

ggplot(volcano_df, aes(avg_logFC, -log10(p_val_adj), color = avg_logFC)) + geom_point(alpha=0.5) + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC > 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5) + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC < 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5) + 
  # ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(avg_logFC, -log10(p_val_adj), label = gene), color = "black", fontface = 3) +
  
  geom_hline(yintercept = -log2(0.10), linetype = "dotted") +
  theme_classic(base_size = 18) + 
  xlim(values = c(-0.5,0.5)) +
  
  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC", y = "-log10(Padj)")
ggsave("results/manuscript/figure2/volcano_nk_3v1.png", width = 6, height = 5)



ggplot(volcano_df, aes(avg_logFC, -log10(p_val_adj), color = avg_logFC)) + geom_point(alpha=0.5) + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC > 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5, color = "black") + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC < 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5, color = "black") + 
  # ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(avg_logFC, -log10(p_val_adj), label = gene), color = "black", fontface = 3) +
  
  geom_hline(yintercept = -log2(0.10), linetype = "dotted") +
  theme_classic(base_size = 18) + 
  xlim(values = c(-0.5,0.5)) +
  
  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC", y = "-log10(Padj)")
ggsave("results/manuscript/figure2/volcano_nk_3v1_2.png", width = 6, height = 5)



volcano_df <- DEG_cluster_df %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster)) %>% filter(cluster == "10 NK adaptive")
volcano_df$avg_logFC <- -volcano_df$avg_logFC

ggplot(volcano_df, aes(avg_logFC, -log10(p_val_adj), color = avg_logFC)) + geom_point(alpha=0.5) + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC > 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5, color = "black") + 
  ggrepel::geom_text_repel(data = volcano_df %>% filter(p_val_adj < 10^-20 & avg_logFC < 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log10(p_val_adj), label = gene), fontface = "italic", size = 3.5, color = "black") + 
  # ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(avg_logFC, -log10(p_val_adj), label = gene), color = "black", fontface = 3) +
  
  geom_hline(yintercept = -log2(0.10), linetype = "dotted") +
  theme_classic(base_size = 18) + 
  xlim(values = c(-0.5,0.5)) +
  
  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC", y = "-log10(Padj)")
ggsave("results/manuscript/figure2/volcano_nk_2v1_2.png", width = 6, height = 5)









## Meta picture
DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% 
  ggplot(aes(reorder(cluster, n), n, fill = cluster, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(22)) +
  labs(x = "", y = "number of DEG as a function of time") + theme_classic(base_size = 12) + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/manuscript/figure2/bar_2v1_nDEG.pdf", width = 6, height = 4)

DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v2")  %>% 
  ggplot(aes(reorder(cluster, n), n, fill = cluster, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(22)) +
  labs(x = "", y = "number of DEG as a function of time") + theme_classic(base_size = 12) + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/manuscript/figure2/bar_3v2_nDEG.pdf", width = 6, height = 4)

DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v1")  %>% 
  ggplot(aes(reorder(cluster, n), n, fill = cluster, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values =  getPalette(22)) +
  labs(x = "", y = "number of DEG as a function of time") + theme_classic(base_size = 12) + facets_nice + geom_text() + theme(legend.position = "none")
ggsave("results/manuscript/figure2//bar_3v1_nDEG.pdf", width = 6, height = 4)



DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v1")  %>% 
  ggplot(aes(reorder(cluster, n), n, label = n)) + 
  geom_point(size = 5, color = "darkred") + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n)) +
  # geom_text(size = 4, color = "gray80") +
  coord_flip() + scale_fill_manual(values =  getPalette(22)) +
  theme_classic(base_size = 17) +
  labs(x = "", y = "number of DEG \n3 mo vs dg") + facets_nice + theme(legend.position = "none") + ggpubr::rotate_x_text(90)
ggsave("results/manuscript/figure2//lolli_3v1_nDEG.pdf", width = 7, height = 5)

DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "2v1")  %>% 
  ggplot(aes(reorder(cluster, n), n, label = n)) + 
  geom_point(size = 5, color = "darkgreen") + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n)) +
  # geom_text(size = 4, color = "gray80") +
  coord_flip() + scale_fill_manual(values =  getPalette(22)) +
  theme_classic(base_size = 17) +
  labs(x = "", y = "number of DEG \n3 mo vs dg") + facets_nice + theme(legend.position = "none") + ggpubr::rotate_x_text(90)
ggsave("results/manuscript/figure2//lolli_2v1_nDEG.pdf", width = 7, height = 5)


df_3v1 <- DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v1") %>% arrange(desc(n))


df <- DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>%  mutate(cluster = factor(cluster, levels = df_3v1$cluster))
df <- df[order(df$cluster), ]

df %>% 
  filter(!is.na(cluster)) %>% 
  ggplot(aes(cluster, n, label = n, color = timepoint)) + geom_path(aes(group=timepoint)) +
  geom_point(size = 5) + 
  # geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n)) +
  # geom_text(size = 4, color = "gray80") +
  scale_fill_manual(values =  getPalette(22)) +
  theme_classic(base_size = 17) +
  labs(x = "", y = "number of DEG") + ggpubr::rotate_x_text(90) + scale_color_manual(values = c("darkgreen", "darkred", "darkblue")) + theme(legend.position = "top") + labs(color = "")
ggsave("results/manuscript/figure2//lolli_nDEG.pdf", width = 7, height = 6)






DEG_cluster_df %>% group_by(timepoint, cluster, celltype) %>% summarise(n = n()) %>% filter(timepoint == "3v1")  %>% 
  mutate(cluster = ifelse(cluster == "15 B naive/resting CMV-specific", "15 B naive/resting\nCMV-specific", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "13 B class-switched memory", "13 B class-swtiched\nmemory", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "16 T-cell exhausted cycling", "16 T-cell\nexhausted cycling", cluster)) %>% 
  
  ggplot(aes(reorder(cluster, n), n, label = n)) + 
  geom_point(size = 5, color = "darkred") + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n)) +
  # geom_text(size = 4, color = "gray80") +
  coord_flip() + 
  theme_classic(base_size = 17) +
  labs(x = "", y = "number of DEG \n3 mo vs dg") + facets_nice + theme(legend.position = "none") + ggpubr::rotate_x_text(45)
ggsave("results/manuscript/figure2//lolli_3v1_nDEG2.pdf", width = 5.5, height = 7)




## combine to see which time point has bigger effect, 1mo or 3mo
DEG_cluster_2v1_summary <- DEG_cluster_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "2v1")
DEG_cluster_3v1_summary <- DEG_cluster_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "3v1")

deg_cluster_cummary <- merge(DEG_cluster_2v1_summary[,-1], DEG_cluster_3v1_summary[,-1], by = "cluster", all = T) 
colnames(deg_cluster_cummary)[2:3] <- c("t2v1", "t3v1") 

deg_cluster_cummary %>% mutate(log2fc = log2(t3v1/t2v1)) %>% filter(cluster %in% Idents(lag3_seurat)) %>% 
  mutate(log2fc = ifelse(is.na(log2fc), 4, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) > 5, 5, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) < -5, -5, log2fc)) %>% 
  
  ggplot(aes(reorder(cluster,log2fc),log2fc,fill=reorder(cluster,log2fc))) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + scale_fill_manual(values = getPalette(22)) + 
  theme_classic(base_size = 12) + theme(legend.position = "none") + ylim(values = c(-5,5)) + geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted") + labs(x = "")
ggsave("results/manuscript/figure2//bar_2v1_3v2_nDEG_log2fc.pdf", width = 6, height = 4)


lag3_seurat <- AddModuleScore(lag3_seurat, features = list(c("LAG3")), name = "LAG3")
lag3_seurat <- AddModuleScore(lag3_seurat, features = list(c("PDCD1")), name = "PDCD1")




df <- lag3_seurat@meta.data %>% mutate(cluster = Idents(lag3_seurat)) %>% 
  group_by(cluster) %>% select(PDCD11,LAG31)
df <- rowMeans(df[,2:3])
  
expression_df <- lag3_seurat@meta.data %>% mutate(cluster = Idents(lag3_seurat), lag3pdcd1 = df) %>% 
  filter(timepoint == "1") %>% 
  group_by(cluster) %>% 
  summarise(lag3_median = mean(LAG31), pdcd1_median = mean(PDCD11), tot_median = mean(lag3pdcd1))

deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t2v1, lag3_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3 expression at baseline")
ggsave("results/manuscript/figure2//scatter_2v1_deg_lag3.pdf", width = 6, height = 4)

deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t2v1, pdcd1_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3 expression at baseline")
ggsave("results/manuscript/figure2//scatter_2v1_deg_pdcd1.pdf", width = 6, height = 4)

deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t2v1, tot_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3+PDCD1 expression at baseline")
ggsave("results/manuscript/figure2//scatter_2v1_deg_lag3pdcd1.pdf", width = 6, height = 4)


deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t3v1, lag3_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3 expression at baseline")
ggsave("results/manuscript/figure2//scatter_3v1_deg_lag3.pdf", width = 6, height = 4)

deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t3v1, pdcd1_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3 expression at baseline")
ggsave("results/manuscript/figure2//scatter_3v1_deg_pdcd1.pdf", width = 6, height = 4)

deg_cluster_cummary %>% left_join(expression_df) %>% 
  ggplot(aes(t3v1, tot_median, label = cluster)) + geom_point() + ggrepel::geom_text_repel() + ggpubr::stat_cor() + 
  labs(x = "number of DEG as a function of time", y = "mean LAG3+PDCD1 expression at baseline")
ggsave("results/manuscript/figure2//scatter_3v1_deg_lag3pdcd1.pdf", width = 6, height = 4)


## combine

## 2v1
DEG_cluster_r_summary <- DEG_cluster_r_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "2v1")
DEG_cluster_n_summary <- DEG_cluster_n_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "2v1")

deg_cluster_cummary <- merge(DEG_cluster_r_summary[,-1], DEG_cluster_n_summary[,-1], by = "cluster", all = T) 
colnames(deg_cluster_cummary)[2:3] <- c("R", "N") 

deg_cluster_cummary <- deg_cluster_cummary %>% 
  mutate(N = ifelse(is.na(N), 0, N)) %>% 
  mutate(R = ifelse(is.na(R), 0, R)) %>% 
  mutate(log2fc = log2(R/N)) %>% 
  filter(cluster %in% Idents(lag3_seurat)) %>% 
  mutate(fill   = ifelse(log2fc > 0, "up in responders", "up in non-responders")) %>% 
  mutate(fill   = ifelse(abs(log2fc) < 1, "unsigf", fill)) %>% 
  mutate(fill   = factor(as.character(fill), levels = c("up in responders", "unsigf", "up in non-responders"))) 

deg_cluster_cummary %>% arrange(desc(log2fc))

deg_cluster_cummary %>% 
  filter(!cluster %in% c("20 low quality", "23 B plasma")) %>% 
  mutate(log2fc = ifelse(is.na(log2fc), 5, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) > 5, 5, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) < -5, -5, log2fc)) %>% 
  
  ggplot(aes(reorder(cluster,log2fc),log2fc,fill=fill)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + 
  scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + 
  labs(fill = "") +
  theme_classic(base_size = 12) + theme(legend.position = "top") +
  ylim(values = c(-5,5)) + geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted") + labs(x = "")
ggsave("results/manuscript/figure2//bar_r_v_n_2v1_nDEG_log2fc.pdf", width = 6, height = 4)


## 3v1
DEG_cluster_r_summary <- DEG_cluster_r_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "3v1")
DEG_cluster_n_summary <- DEG_cluster_n_df %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% filter(timepoint == "3v1")

deg_cluster_cummary <- merge(DEG_cluster_r_summary[,-1], DEG_cluster_n_summary[,-1], by = "cluster", all = T) 
colnames(deg_cluster_cummary)[2:3] <- c("R", "N") 

deg_cluster_cummary <- deg_cluster_cummary %>% 
  mutate(N = ifelse(is.na(N), 0, N)) %>% 
  mutate(R = ifelse(is.na(R), 0, R)) %>% 
  mutate(log2fc = log2(R/N)) %>% 
  filter(cluster %in% Idents(lag3_seurat)) %>% 
  mutate(fill   = ifelse(log2fc > 0, "up in responders", "up in non-responders")) %>% 
  mutate(fill   = ifelse(abs(log2fc) < 1, "unsigf", fill)) %>% 
  mutate(fill   = factor(as.character(fill), levels = c("up in responders", "unsigf", "up in non-responders"))) 

deg_cluster_cummary %>% arrange(desc(log2fc))

deg_cluster_cummary %>% 
  filter(!cluster %in% c("20 low quality", "23 B plasma")) %>% 
  mutate(log2fc = ifelse(is.na(log2fc), 5, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) > 5, 5, log2fc)) %>% 
  mutate(log2fc = ifelse((log2fc) < -5, -5, log2fc)) %>% 
  
  ggplot(aes(reorder(cluster,log2fc),log2fc,fill=fill)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + 
  scale_fill_manual(values = c("salmon", "lightgrey", "dodgerblue")) + 
  labs(fill = "") +
  theme_classic(base_size = 12) + theme(legend.position = "top") +
  ylim(values = c(-5,5)) + geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted") + labs(x = "")
ggsave("results/manuscript/figure2//bar_r_v_n_3v1_nDEG_log2fc.pdf", width = 6, height = 4)




## Genes
inhibitory_genes_r <- DEG_cluster_r_df %>% filter(gene %in% inhibitory_long) %>% 
  group_by(cluster, timepoint) %>% summarise(n = n()) %>% mutate(key = paste(cluster,timepoint))

inhibitory_genes_n <- DEG_cluster_n_df %>% filter(gene %in% inhibitory_long) %>% 
  group_by(cluster, timepoint) %>% summarise(n = n()) %>% mutate(key = paste(cluster,timepoint)) %>% ungroup %>% select(-cluster,-timepoint)

inhibitory_tot <- merge(inhibitory_genes_r, inhibitory_genes_n, by = "key", all.x = T) %>% dplyr::rename(n.r = n.x, n.n = n.y) 

inhibitory_tot %>% melt(id = c("key", "cluster", "timepoint")) %>% 
  filter(timepoint != "3v2") %>% 
  ggplot(aes(cluster,value,fill=variable)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + facet_wrap(~timepoint)

DEG_cluster_r_df %>% filter(gene %in% c("LAG3", "PDCD1"))
DEG_cluster_n_df %>% filter(gene %in% c("LAG3", "PDCD1"))


## Enrichemnts
enrichment_n_2v1_up <- fread("results/effectOnClusters/enrichment/meta_nonresponder_up_2v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))
enrichment_n_2v1_dn <- fread("results/effectOnClusters/enrichment/meta_nonresponder_down_2v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))

enrichment_r_2v1_up <- fread("results/effectOnClusters/enrichment/meta_responder_up_2v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))
enrichment_r_2v1_dn <- fread("results/effectOnClusters/enrichment/meta_responder_down_2v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))

enrichment_n_3v1_up <- fread("results/effectOnClusters/enrichment/meta_nonresponder_up_3v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))
enrichment_n_3v1_dn <- fread("results/effectOnClusters/enrichment/meta_nonresponder_down_3v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))

enrichment_r_3v1_up <- fread("results/effectOnClusters/enrichment/meta_responder_up_3v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))
enrichment_r_3v1_dn <- fread("results/effectOnClusters/enrichment/meta_responder_down_3v1.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypes()) %>% filter(p.adjust < 0.05) %>% mutate(key = paste(ID, cluster))


enrichment_r_2v1_up[enrichment_r_3v1_up$cluster == "10 NK adaptive", ] %>% View
enrichment_r_3v1_up[enrichment_r_3v1_up$cluster == "10 NK adaptive", ] %>% View
enrichment_n_3v1_up[enrichment_n_3v1_up$cluster == "10 NK adaptive", ]









##### Figure 3

sigf_means_r_0m <- fread("results/cellphonedb/out/r_0m/significant_means.txt")
sigf_means_r_1m <- fread("results/cellphonedb/out/r_1m/significant_means.txt")
sigf_means_r_3m <- fread("results/cellphonedb/out/r_3m/significant_means.txt")

sigf_means_n_0m <- fread("results/cellphonedb/out/n_0m/significant_means.txt")
sigf_means_n_1m <- fread("results/cellphonedb/out/n_1m/significant_means.txt")
sigf_means_n_3m <- fread("results/cellphonedb/out/n_3m/significant_means.txt")


## How many interactions are there?
sigf_means_r_0m_pairs <- getPairAmounts(sigf_means_r_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))
sigf_means_r_1m_pairs <- getPairAmounts(sigf_means_r_1m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))
sigf_means_r_3m_pairs <- getPairAmounts(sigf_means_r_3m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))

sigf_means_n_0m_pairs <- getPairAmounts(sigf_means_n_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))
sigf_means_n_1m_pairs <- getPairAmounts(sigf_means_n_1m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))
sigf_means_n_3m_pairs <- getPairAmounts(sigf_means_n_3m) %>% filterRedundantPairs() %>% mutate(cluster1 = getClusterPhenotypes(extractClusterNumber(cluster1)), cluster2 = getClusterPhenotypes(extractClusterNumber(cluster2)))

table(sigf_means_r_3m_pairs$cluster1)
table(sigf_means_r_0m_pairs$cluster2)


p <- plotDiffHeatmap(df1 = sigf_means_n_0m_pairs, df2 = sigf_means_r_0m_pairs)
ggsave(plot = p, "results/manuscript/figure3/heatmap_sigf_rvn_0.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_n_1m_pairs, df2 = sigf_means_r_1m_pairs)
ggsave(plot = p, "results/manuscript/figure3/heatmap_sigf_rvn_1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_n_3m_pairs, df2 = sigf_means_r_3m_pairs)
ggsave(plot = p, "results/manuscript/figure3/heatmap_sigf_rvn_3.png", width = 8, height = 6)


p <- plotDiffHeatmap(df1 = sigf_means_r_0m_pairs, df2 = sigf_means_r_3m_pairs)
ggsave(plot = p, "results/manuscript/figure3/heatmap_sigf_r_3v0.png", width = 8, height = 6)

p <- plotDiffHeatmap(sigf_means_n_0m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/manuscript/figure3/heatmap_sigf_n_3v0.png", width = 8, height = 6)


ab_test <- ab
colnames(ab_test) <- colnames(ab_test) %>% make.names() %>% make.unique()
colnames(ab_test) <- gsub("CMV", "CMV.IgG", colnames(ab_test))

ab_test %>% #filter(pos == "pos") %>% 
  select(overall, CMV.IgG.D0, CMV.IgG.1mo, CMV.IgG.3mo) %>% melt(id = "overall") %>% 
  ggplot(aes(overall,value,fill=overall)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + geom_jitter(size = 0.5) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = c("lightgrey", "salmon")) + theme(legend.position = "none") + facets_nice
ggsave("results/manuscript/figure3/box_cmv_ab.pdf", width = 5, height = 4)

ab_test %>% filter(pos == "pos") %>% 
  select(overall, CMV.IgG.D0, CMV.IgG.1mo, CMV.IgG.3mo) %>% melt(id = "overall") %>% 
  ggplot(aes(overall,value,fill=overall)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable) + geom_jitter(size = 0.5) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = c("lightgrey", "salmon")) + theme(legend.position = "none") + facets_nice
ggsave("results/manuscript/figure3/box_cmv_ab_pos.pdf", width = 5, height = 4)




###### Figure 4
dir.create("results/manuscript/figure4/", showWarnings = F)

div20k <- fread("results/cmv/tcrb_info.txt")

## Lag3 cohort, cmv pos vs neg
div_20k %>%
  filter(cohort == "melanoma" & timepoint == "0m" & regimen == "antiPD1+antiLAG3") %>% 
  melt(id = meta_columns) %>% 
  filter(!is.na(pos)) %>% 
  mutate(pos = paste(overall, pos)) %>% 
  mutate(pos = factor(pos, levels = c("R pos", "N pos", "R neg", "N neg"))) %>% 
  # filter(io_stat == "IO.naive") %>% 
  
  ggplot(aes(pos, value, fill = pos)) + 
  theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + 
  facet_wrap(~variable, scales = "free", ncol = 5) + 
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(4)) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice #+ ggpubr::stat_compare_means(ref = "R pos", label = "p.signif")
ggsave("results/manuscript/figure4/box_lag3_cmv_response_0m.pdf", width = 12, height = 5)


div_20k$cohort

div_20k %>%
  # filter(cohort == "melanoma" & timepoint == "0m" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>% 
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>% 
  
  select(clonality, inverseSimpsonIndex_mean, meta_columns) %>% 
  melt(id = meta_columns) %>% 
  filter(!is.na(pos)) %>% 
  filter(pos != "Unknown") %>% 
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>% 
  mutate(pos = paste(overall, pos)) %>% 
  mutate(pos = factor(pos, levels = c("R pos", "N pos", "healthy pos", "R neg", "N neg", "healthy neg"))) %>% 
  # filter(io_stat == "IO.naive") %>% 
  
  ggplot(aes(pos, value, fill = pos)) + 
  theme_classic(base_size = 12) +
  # geom_violin() + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.3) + 
  
  facet_wrap(~variable, scales = "free", ncol = 5) + 
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(6)) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice #+ ggpubr::stat_compare_means(ref = "R pos", label = "p.signif")
ggsave("results/manuscript/figure4/box_lag3_cmv_response_all.pdf", width = 5, height = 3)


div_20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3") %>% 
  melt(id = meta_columns) %>% 
  filter(!is.na(pos)) %>% 
  
  ggplot(aes(timepoint, value, fill = pos)) + 
  theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + 
  facet_wrap(~variable, scales = "free", ncol = 5) + 
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = c("lightgrey", "salmon")) + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice
ggsave("results/cmv/box_lag3_cmv_timepoint.pdf", width = 12, height = 5)












lag3_seurat@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>% 
  filter(cluster %in% cl.to.include$cluster) %>% 
  mutate(cluster = as.character(cluster)) %>% 
  mutate(cluster = ifelse(cluster == "15 B naive/resting CMV-specific", "15 B resting CMV-specific", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "13 B class-switched memory", "13 B class-switched", cluster)) %>% 
  mutate(cluster = ifelse(cluster == "16 T-cell exhausted cycling", "16 T-cell exhausted cycling", cluster)) %>% 
  
  ggplot(aes(timepoint,prop, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, scales = "free_y", strip.position="right", ncol = 6) + ggpubr::stat_compare_means(label = "p", label.y.npc = 0.9) + scale_fill_manual(values = c("red3", "darkolivegreen4", "darkolivegreen4")) +
  theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + theme(strip.text.y = element_text(size = 10))
ggsave("results/manuscript/figure1/box_lag3_timepoint.pdf", width = 12, height = 6)



table(lag3_seurat$cluster)

nk.cells.to.keep <- lag3_seurat@meta.data %>% filter(grepl(pattern = "NK", cluster)) %>% pull(barcode)
nk_seurat <- subset(lag3_seurat, cells = nk.cells.to.keep)
nk_seurat <- nk_seurat %>% getLatentUMAP()

nk_adaptive <- subset(lag3_seurat, idents = "10 NK adaptive")
nk_adaptive <- nk_adaptive %>% getLatentUMAP()

DimPlot(nk_seurat, label = T, repel = T, cols = getPalette3(4), label.size = 7) + theme_void(base_size = 17) + theme(legend.position = "none")
ggsave("results/manuscript/figure1/umap_nk.png", width = 6, height = 4)


FeaturePlot(nk_seurat, features = c("KLRC2", "LAG3"), order = T, cols = c("gray90", "red3"), label = T, repel = T) & NoAxes()
ggsave("results/manuscript/figure1/umap_adaptive_nk.png", width = 7, height = 3)

p <- FeaturePlot(nk_seurat, features = c("CD3E", "KLRC2", "LAG3", "TIGIT", "ZEB2", "CD2"), order = F, cols = c("gray90", "red3"), label = T, repel = T, ncol = 3, label.size = 5, pt.size = 1e-15) & NoAxes() & theme(title = element_text(face = "bold.italic"))
ggsave(plot = p, "results/manuscript/figure1/umap_adaptive_nk2.png", width = 10, height = 5)

p <- FeaturePlot(nk_seurat, features = c("CD3E", "KLRC2", "LAG3", "TIGIT", "PDCD1", "CD2"), order = T, cols = c("gray90", "red3"), label = F, repel = T, ncol = 3, label.size = 10, pt.size = 1e-1) & theme_void(base_size = 30) & theme(title = element_text(face = "bold.italic")) #+ theme(plot.margin = margin(rep(1e-5,4), "cm"), plot.margin)
ggsave(plot = p, "results/manuscript/figure1/umap_adaptive_nk2.1.png", width = 10*2, height = 5*2)


p <- FeaturePlot(nk_seurat, features = dufva_genes, order = T, cols = c("gray90", "red3"), label = F) & NoAxes()

nk_markers <- FindAllMarkers(nk_seurat, test.use = "t")
fwrite(nk_markers, "results/manuscript/nk_markers.txt", sep = "\t", quote = F, row.names = F)

nk_markers$avg_logFC
FeaturePlot(nk_adaptive, features = c("LAG3"), order = T)


nk_markers %>% filter(gene %in% inhibitory_markers)

top10 <- nk_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

DoHeatmap(subset(nk_seurat, downsample = 100), features = top10$gene, group.colors = getPalette(7), angle = 270, hjust = 1, draw.lines = T, raster = T) + guides(color = FALSE) + 
  # scale_fill_gradientn(colors = c("darkblue", "white", "darkred")) +
  # scale_fill_gradientn(colors = c("blue", "white", "red")) +
  # viridis::scale_fill_viridis() +
  theme(legend.position = "bottom") + theme(legend.position = "none")
ggsave("results/manuscript/figure1/heatmap_de_nk.png", width = 4, height = 11)
ggsave("results/manuscript/figure1/heatmap_de_nk.pdf", width = 4, height = 11)

DoHeatmap(subset(nk_seurat, downsample = 100), features = top10$gene, group.colors = getPalette(7), angle = 270, hjust = 1, draw.lines = T, raster = T) + guides(color = FALSE) + 
theme(legend.position = "bottom") + theme(legend.position = "top")
ggsave("results/manuscript/figure1/heatmap_de_nk2.pdf", width = 4, height = 11)

