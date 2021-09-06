
## Run slingshot to get pseudotime
library(slingshot)
library(SingleCellExperiment)

nk_sce   <- as.SingleCellExperiment(nk_seurat)
nk_sling <- slingshot(data = nk_sce, clusterLabels = 'cluster', reducedDim = "LATENT_UMAP", start.clus = "14 NK CD56bright")
nk_curve <- SlingshotDataSet(nk_sling)

cololors <- nk_sling$ident
levels(cololors) <- getPalette3(4)
curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
curve_cols <- c("black")

pdf("results/nk_sling_pca.pdf", width = 5, height = 5)
plot(reducedDims(nk_sling)$LATENT_UMAP, col = as.character(cololors), pch = 16, asp = 1, xlim=c(-5,5), cex = 0.3)
for(i in 1){
  lines(nk_curve@curves[[i]], lwd = 3, col = curve_cols[i])
}
dev.off()



png("results/nk_sling_pca.png", width = 5, height = 5, res = 1086, units = "in")
plot(reducedDims(nk_sling)$LATENT_UMAP, col = as.character(cololors), pch = 16, asp = 1, xlim=c(-5,5), cex = 0.3)
for(i in 1){
  lines(nk_curve@curves[[i]], lwd = 3, col = curve_cols[i])
}
dev.off()

p <- FeaturePlot(nk_seurat, features = c("CD3E", "LAG3", "KLRC2", "TIGIT", "PDCD1", "CD2"), ncol = 3, cols = c("gray90", "red3"), pt.size = 0.01, order = T) &
  theme_void(base_size = 15) &
  theme(legend.position = "none") & theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) #& theme(plot.title = element_text(family = "bold.italic"))
ggsave(plot = p, "results/manuscript/feature_nk_adaptive.png", width = 8, height = 4)



p <- FeaturePlot(nk_seurat, features = c("CD3E", "LAG3", "KLRC2", "TIGIT", "PDCD1", "CD2"), ncol = 3, cols = c("gray90", "red3"), pt.size = 0.01, order = T) &
  theme_void(base_size = 15) & theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) #& theme(plot.title = element_text(family = "bold.italic"))
ggsave(plot = p, "results/manuscript/feature_nk_adaptive2.png", width = 8, height = 4)

p <- FeaturePlot(nk_seurat, features = c("NCAM1", "FCGR3A", "LAG3", "KLRC2", "HAVCR2", "PDCD1"), ncol = 3, cols = c("gray90", "red3"), pt.size = 0.01, order = T) &
  xlim(c(-5,5)) &
  theme_void(base_size = 15) &
  theme(legend.position = "none") & theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) & theme(plot.title = element_text(face = "bold.italic"))
ggsave(plot = p, "results/manuscript/feature_nk_adaptive.png", width = 6, height = 4)
