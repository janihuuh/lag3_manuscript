
dir.create("results/manuscript/timepoint", showWarnings = F)

## Across patients
lag3_seurat$cluster <- Idents(lag3_seurat)

lag3_seurat@meta.data %>% group_by(orig.ident, timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  # filter(timepoint %in% c(0,3)) %>%
  ggplot(aes(timepoint,prop, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~cluster, scales = "free_y") + ggpubr::stat_compare_means(label = "p.format") + geom_jitter(size=0.5) + facets_nice +
  scale_fill_manual(values = getPalette(4)) + theme(legend.position = "none") + labs(x = "months")
ggsave("results/manuscript/timepoint/box_abundance.pdf", width = 9.5, height = 8)



## Get DEGs
deg_timepoints <- lapply(unique(Idents(lag3_seurat)), FUN = getDEGbyCluster, seurat_object = lag3_seurat) %>% rbindlist()
fwrite(deg_timepoints, "results/manuscript/timepoint/deg_timepoint2.txt", sep = "\t", quote = F, row.names = F)
# deg_timepoints <- fread("results/manuscript/timepoint/deg_timepoint2.txt")

deg_timepoints$cluster %>% extractCoarsePhenotype %>% do.call(what = "c")
deg_timepoints$p_val_adj
deg_timepoints %>%
  filter(p_val_adj < 0.05) %>%
  group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(celltype = cluster %>% extractCoarsePhenotype %>% do.call(what = "c")) %>%
  filter(timepoint != "3v2") %>%
  ggplot(aes(timepoint,n,group=cluster,color=cluster)) + geom_path(color="black") + geom_point(size = 3) + facet_wrap(~celltype) + facets_nice + labs(x = "") + scale_color_manual(values = getPalette(21))
ggsave("results/manuscript/timepoint/path_deg.pdf", width = 8, height = 5)

deg_timepoints %>%   filter(p_val_adj < 0.05) %>%
  group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(celltype = cluster %>% extractCoarsePhenotype %>% do.call(what = "c")) %>%
  filter(timepoint != "3v1") %>%
  ggplot(aes(reorder(cluster, n),n,color=celltype,label=n)) + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n), color = "lightgrey") +
  geom_point(size = 10) + geom_text(color = "white") + facets_nice + labs(x = "") + scale_color_manual(values = getPalette(6)) +
  coord_flip() + theme(legend.position = "none")
ggsave("results/manuscript/timepoint/bar_3v1_deg.pdf", width = 5, height = 4)

deg_timepoints %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(celltype = cluster %>% extractCoarsePhenotype %>% do.call(what = "c")) %>%
  filter(timepoint == "3v2") %>%
  ggplot(aes(reorder(cluster, n),n,color=celltype,label=n)) + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n), color = "lightgrey") +
  geom_point(size = 10) + geom_text(color = "white") + facets_nice + labs(x = "") + scale_color_manual(values = getPalette(6)) +
  coord_flip() + theme(legend.position = "none")
ggsave("results/manuscript/timepoint/bar_3v2_deg.pdf", width = 5, height = 4)

deg_timepoints %>% group_by(timepoint, cluster) %>% summarise(n = n()) %>% mutate(celltype = extractCoarsePhenotype(cluster)) %>%
  filter(timepoint == "2v1") %>%
  ggplot(aes(reorder(cluster, n),n,color=celltype,label=n)) + geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n), color = "lightgrey") +
  geom_point(size = 10) + geom_text(color = "white") + facets_nice + labs(x = "") + scale_color_manual(values = getPalette(6)) +
  coord_flip() + theme(legend.position = "none")
ggsave("results/manuscript/timepoint/bar_2v1_deg.pdf", width = 5, height = 4)



## Get enrichment
universe_df = rownames(lag3_seurat)

hallmark_df <-
  lapply(unique(deg_timepoints$timepoint), FUN = function(z){
    lapply(unique(deg_timepoints$direction), FUN = function(y){
      lapply(unique(deg_timepoints$cluster), FUN = function(x){
        p <- deg_timepoints %>% filter(cluster == x & timepoint == z & direction == y) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
        if(!is.null(p)){p %>% mutate(cluster = x, timepoint = z, dir = y)}
    }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist() %>% filter(p.adjust < 0.05)

fwrite(hallmark_df, "results/manuscript/timepoint/hallmark_pathways.txt", sep = "\t", quote = F, row.names = F)

hallmark_df <-
  lapply(unique(deg_timepoints_edit$timepoint), FUN = function(z){
    lapply(unique(deg_timepoints_edit$direction), FUN = function(y){
      lapply(unique(deg_timepoints_edit$cluster), FUN = function(x){
        p <- deg_timepoints_edit %>% filter(cluster == x & timepoint == z & direction == y) %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
        if(!is.null(p)){p %>% mutate(cluster = x, timepoint = z, dir = y)}
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist() %>% filter(p.adjust < 0.05)
fwrite(hallmark_df, "results/manuscript/timepoint/hallmark_pathways_updated.txt", sep = "\t", quote = F, row.names = F)





go_df <-
  lapply(unique(deg_timepoints_edit$timepoint), FUN = function(z){
    lapply(unique(deg_timepoints_edit$direction), FUN = function(y){
      lapply(unique(deg_timepoints_edit$cluster), FUN = function(x){
        p <- deg_timepoints_edit %>% filter(cluster == x & timepoint == z & direction == y) %>% getHypergeometric(universe_df = universe_df, term_df = go)
        if(!is.null(p)){p %>% mutate(cluster = x, timepoint = z, dir = y)}
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist() %>% filter(p.adjust < 0.05)
fwrite(go_df, "results/manuscript/timepoint/go_pathways_updated.txt", sep = "\t", quote = F, row.names = F)

go_df <-
  lapply(unique(deg_timepoints_edit2$timepoint), FUN = function(z){
    lapply(unique(deg_timepoints_edit2$direction), FUN = function(y){
      lapply(unique(deg_timepoints_edit2$cluster), FUN = function(x){
        p <- deg_timepoints_edit2 %>% filter(cluster == x & timepoint == z & direction == y) %>% getHypergeometric(universe_df = universe_df, term_df = go)
        if(!is.null(p)){p %>% mutate(cluster = x, timepoint = z, dir = y)}
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist() %>% filter(p.adjust < 0.05)
fwrite(go_df, "results/manuscript/timepoint/go_pathways_updated2.txt", sep = "\t", quote = F, row.names = F)
