
dir.create("results/manuscript_new/tcr/")

lag3_cd4_tcrseq <- fread("results/tcrb/diversity/lag3_cd4.diversity.aa.exact.txt")
lag3_cd8_tcrseq <- fread("results/tcrb/diversity/lag3_cd8.diversity.aa.exact.txt")

lag3_cd4_tcrseq$timepoint <- as.factor(paste0(substr(lag3_cd4_tcrseq$sample_id, 6,6), "mo"))
lag3_cd8_tcrseq$timepoint <- as.factor(paste0(substr(lag3_cd8_tcrseq$sample_id, 6,6), "mo"))

lag3_cd4_tcrseq$type <- "CD4+ LAG3+\nsorted cells"
lag3_cd8_tcrseq$type <- "CD8+ LAG3+\nsorted cells"

lag3_tcrseq <- rbind(lag3_cd4_tcrseq, lag3_cd8_tcrseq)

lag3_tcrseq %>% ggplot(aes(timepoint, 1-normalizedShannonWienerIndex_mean, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
  # ggpubr::stat_compare_means(label="p") +
  facet_wrap(~type, scales = "free_y") + labs(x = "", y = "clonality") +
  scale_fill_manual(values = getPalette5(8)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") + facets_nice + theme(legend.position = "none") + scale_fill_manual(values = c("darkolivegreen4", "darkolivegreen4")) + ggpubr::rotate_x_text(angle = 45)
ggsave("results/manuscript_new/tcr/box_lag3_cells.pdf", width = 3.5, height = 3.5)

lag3_tcrseq %>% ggplot(aes(type, 1-normalizedShannonWienerIndex_mean, fill = type)) + geom_boxplot() + geom_jitter() + ggpubr::stat_compare_means(label="p") + labs(x = "", y = "clonality") +
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + facets_nice + theme(legend.position = "none")
ggsave("results/manuscript_new/tcr/box_lag3_cells_diff.pdf", width = 5, height = 3.5)




## To gliph
lag3_cd8_chains <- lapply(list.files("data/tcrseq/lag3_sorted/cd8/", full.names = T), function(x) fread(x) %>% mutate(name = extractFileName(x))) %>% rbindlist() %>% as.data.frame()
lag3_cd8_chains_gliph <- lag3_cd8_chains %>% vdjToGliph()
fwrite(lag3_cd8_chains_gliph, "results/manuscript_new/tcr/lag3_cd8_input_gliph.txt", sep = "\t", quote = F, row.names = F)

lag3_cd4_chains <- lapply(list.files("data/tcrseq/lag3_sorted/cd4/", full.names = T), function(x) fread(x) %>% mutate(name = extractFileName(x))) %>% rbindlist() %>% as.data.frame()
lag3_cd4_chains_gliph <- lag3_cd4_chains %>% vdjToGliph()
fwrite(lag3_cd4_chains_gliph, "results/manuscript_new/tcr/lag3_cd4_input_gliph.txt", sep = "\t", quote = F, row.names = F)

lag3_cd8 <- readGliphFile("results/gliph/lag3_cd8_gliph_online.csv") %>% filter(vb_score < 0.05 & number_unique_cdr3 > 3)
lag3_cd8_patterns <- unique(lag3_cd8$pattern)

lag3_cd4 <- readGliphFile("results/gliph/lag3_cd4_gliph_online.csv") %>% filter(vb_score < 0.05 & number_unique_cdr3 > 3)
lag3_cd4_patterns <- unique(lag3_cd4$pattern)


lag3_cd8 %>% group_by(pattern) %>% summarise(n = n()) %>% filter(n>4) %>%
  ggplot(aes(reorder(pattern,n),n)) + geom_bar(stat="identity",fill="lightgrey") + coord_flip() + labs(x = "", y = "TCRs with motif")
ggsave("results/manuscript_new/tcr/bar_cd8_gliph.pdf", width = 3, height = 5)

lag3_cd8 %>% group_by(pattern) %>% summarise(n = n()) %>% filter(n>4) %>%
  ggplot(aes(reorder(pattern,n),n)) + geom_bar(stat="identity",fill="lightgrey") + labs(x = "", y = "TCRs with motif") + ggpubr::rotate_x_text(90)
ggsave("results/manuscript_new/tcr/bar_cd8_gliph2.pdf", width = 5, height = 3)


lag3_cd4 %>% group_by(pattern) %>% summarise(n = n()) %>% filter(n>4) %>%
  ggplot(aes(reorder(pattern,n),n)) + geom_bar(stat="identity",fill="lightgrey") + coord_flip() + labs(x = "", y = "TCRs with motif")
ggsave("results/manuscript_new/tcr/bar_cd4_gliph.pdf", width = 3, height = 4)







## CMV to
div20k <- fread("results/cmv/tcrb_info.txt")
meta_columns <- c("name", "response", "type", "regimen", "overall", "timepoint", "celltype", "io_stat", "cohort", "pos")

table(div20k$cohort)


div20k %>% filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3") %>% pull(name) %>% unique() %>% length()

## Lag3 cohort, cmv pos vs neg
div20k %>%
  filter(cohort == "melanoma" & timepoint == "0m" & regimen == "antiPD1+antiLAG3") %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  mutate(pos = paste(pos)) %>%
  mutate(pos = factor(pos, levels = c("pos", "N pos", "R neg", "N neg"))) %>%
  # filter(io_stat == "IO.naive") %>%

  ggplot(aes(pos, value, fill = pos)) +
  theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(4)) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice #+ ggpubr::stat_compare_means(ref = "R pos", label = "p.signif")
ggsave("results/manuscript/figure4/box_lag3_cmv_response_0m.pdf", width = 12, height = 5)

table(div20k$regimen)
table(div20k$cohort)

div20k %>%
  filter(regimen == "antiPD1+antiLAG3") %>%
  filter(timepoint != "6m") %>%
  mutate(prepost = ifelse(timepoint == "0m", "pre", "post")) %>%

  ggplot(aes(timepoint, clonality)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5, alpha = 0.5) + facet_wrap(~io_stat) + ggpubr::stat_compare_means(label = "p")


table(div20k$type)

div20k %>%
  mutate(io_stat = ifelse(io_stat == "IO.naive", "IO naive", io_stat)) %>%
  mutate(io_stat = ifelse(io_stat == "Prior.IO", "IO refractory", io_stat)) %>%
  filter(regimen == "antiPD1+antiLAG3") %>%
  filter(timepoint != "6m") %>%
  mutate(prepost = ifelse(timepoint == "0m", "pre", "post")) %>%
  mutate(prepost = factor(as.character(prepost), levels = c("pre", "post"))) %>%
  ggplot(aes(prepost, clonality, fill = prepost)) + geom_boxplot(outlier.shape = NA) +
  # geom_path(aes(group = name)) +
  geom_jitter(size = 0.5, alpha = 0.5) + facet_wrap(~io_stat) +
  theme_classic(base_size = 15) +
  # ggpubr::stat_compare_means(label = "p", method = "t") +
  scale_fill_manual(values = c("red3", "darkolivegreen4")) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45)
ggsave("results/manuscript/box_clonality_io_stat.pdf", width = 3, height = 3)




div_20k$cohort

div20k %>%
  # filter(cohort == "melanoma" & timepoint == "0m" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%

  dplyr::select(clonality, inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  # mutate(pos = factor(pos, levels = c("melanoma pos", "healthy pos", "melanoma neg", "healthy neg"))) %>%
  # filter(io_stat == "IO.naive") %>%

  ggplot(aes(pos, value, fill = pos)) +
  theme_classic(base_size = 12) +
  # geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.3) +

  facet_wrap(~variable, scales = "free", ncol = 5) +
  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette(6)) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice #+ ggpubr::stat_compare_means(ref = "R pos", label = "p.signif")
ggsave("results/manuscript/timepoint/box_lag3_cmv_response_all.pdf", width = 5, height = 3)





div20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%

  dplyr::select(clonality, inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  filter(pos %in% c("healthy neg", "healthy pos", "melanoma neg", "melanoma pos")) %>%
  mutate(pos = plyr::revalue(pos, replace = c("healthy neg" = "healthy CMV-", "healthy pos" = "healthy CMV+", "melanoma neg" = "SKCM CMV-", "melanoma pos" = "SKCM CMV+"))) %>%
  mutate(type = ifelse(grepl("healthy", x = pos), "healthy", "SKCM")) %>%
  # mutate(pos = factor(pos, levels = c("melanoma pos", "healthy pos", "melanoma neg", "healthy neg"))) %>%
  # filter(io_stat == "IO.naive") %>%

  ggplot(aes(pos, value, fill = pos)) +
  theme_classic(base_size = 12) +
  # geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.3) +

  facet_wrap(~variable, scales = "free", ncol = 5) +
  # ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette3(6)) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("melanoma pos", "healthy pos"), c("melanoma neg", "healthy neg")), step_increase = 0.05)
ggsave("results/manuscript_new/tcr/box_lag3_cmv_response_all.pdf", width = 5, height = 5)





div20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%
  dplyr::select(clonality, Simpson = inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  filter(pos %in% c("healthy neg", "healthy pos", "melanoma neg", "melanoma pos")) %>%
  mutate(pos = plyr::revalue(pos, replace = c("healthy neg" = "healthy CMV-", "healthy pos" = "healthy CMV+", "melanoma neg" = "SKCM CMV-", "melanoma pos" = "SKCM CMV+"))) %>%
  mutate(type = ifelse(grepl("healthy", x = pos), "healthy", "SKCM")) %>%
  mutate(positivity = ifelse(grepl("\\+", x = pos), "CMV+", "CMV-")) %>%

  # mutate(pos = factor(pos, levels = c("melanoma pos", "healthy pos", "melanoma neg", "healthy neg"))) %>%
  # filter(io_stat == "IO.naive") %>%

  ggplot(aes(positivity, value, fill = type)) +
  theme_classic(base_size = 15) +
  # geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, alpha = 0.1) +

  facet_wrap(~variable, scales = "free_y") +
  # ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette3(6)) + theme(legend.position = "top") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) +
  ggsignif::geom_signif(comparisons = list(c("melanoma pos", "healthy pos"), c("melanoma neg", "healthy neg")), step_increase = 0.05) + scale_fill_manual(values = c("blue3", "red3"))
ggsave("results/manuscript_new/tcr/box_lag3_cmv_response_all2.pdf", width = 4, height = 5)


a <- div20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%
  dplyr::select(clonality, Simpson = inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  filter(pos %in% c("healthy neg", "healthy pos", "melanoma neg", "melanoma pos")) %>%
  mutate(pos = plyr::revalue(pos, replace = c("healthy neg" = "healthy CMV-", "healthy pos" = "healthy CMV+", "melanoma neg" = "SKCM CMV-", "melanoma pos" = "SKCM CMV+"))) %>%
  mutate(type = ifelse(grepl("healthy", x = pos), "healthy", "SKCM")) %>%
  mutate(positivity = ifelse(grepl("\\+", x = pos), "CMV+", "CMV-"))

a <- a %>% filter(timepoint == "0m" | timepoint == "")
table(a$type) / 2
a$timepoint

div20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%
  dplyr::select(clonality, Simpson = inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  filter(pos %in% c("healthy neg", "healthy pos", "melanoma neg", "melanoma pos")) %>%
  mutate(pos = plyr::revalue(pos, replace = c("healthy neg" = "healthy CMV-", "healthy pos" = "healthy CMV+", "melanoma neg" = "SKCM CMV-", "melanoma pos" = "SKCM CMV+"))) %>%
  mutate(type = ifelse(grepl("healthy", x = pos), "healthy", "SKCM")) %>%
  mutate(positivity = ifelse(grepl("\\+", x = pos), "CMV+", "CMV-")) %>%
  filter(variable == "clonality") %>%

  filter(timepoint == "0m" | timepoint == "") %>%
  mutate(type = ifelse(type == "SKCM", "Melanoma", "Healthy")) %>%

  ggplot(aes(positivity, value, fill = type)) +
  theme_classic(base_size = 15) +

  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, alpha = 0.1) +

  ggpubr::stat_compare_means(label = "p.format") +
  scale_fill_manual(values = getPalette3(6)) + theme(legend.position = "top") + labs(x = "") + ggpubr::rotate_x_text(angle = 45) +
  ggsignif::geom_signif(comparisons = list(c("melanoma pos", "healthy pos"), c("melanoma neg", "healthy neg")), step_increase = 0.05) + scale_fill_manual(values = c("darkgray", "brown")) + labs(y = "clonality", fill = "")
ggsave("results/manuscript_new/tcr/box_lag3_cmv_response_all3.pdf", width = 3, height = 5)





df <- div20k %>%
  filter(cohort == "melanoma" & regimen == "antiPD1+antiLAG3" | cohort == "healthy") %>%

  dplyr::select(clonality, inverseSimpsonIndex_mean, meta_columns) %>%
  melt(id = meta_columns) %>%
  filter(!is.na(pos)) %>%
  filter(pos != "Unknown") %>%
  mutate(overall = ifelse(overall == "", "healthy", overall)) %>%
  mutate(pos = paste(cohort, pos)) %>%
  filter(pos %in% c("healthy neg", "healthy pos", "melanoma neg", "melanoma pos"))

dim(df)
table(df$cohort)/2
