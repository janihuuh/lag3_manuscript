
facs_tot <- fread("results/manuscript_new/facs/facs_tot.txt")

## Test wlcx for io naive
df <- krskl.df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>%
  filter(Cohort == 1)
df <- df[!duplicated(df), ]

krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]
krskl.df$timepoint <- ifelse(krskl.df$timepoint == "0m", "prior", "post")
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  # message(x)
  y <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(timepoint) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()
  if(n.vars == 2){

    wilcox.test(value~timepoint, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, timepoint == "prior")$value), median.y = median(subset(y, timepoint != "prior")$value),
                                                                                                                mean.x   = mean(subset(y, timepoint == "prior")$value), mean.y   = mean(subset(y, timepoint != "prior")$value), log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))

    }

}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_wlcx_p_io_naive.txt", sep = "\t", quote = F, row.names = F)
p.df <- fread("results/manuscript_new/facs/facs_wlcx_p_io_naive.txt")

ggplot(p.df, aes(x = log2fc.mean, y = -log2(p.adj), color = log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log2(0.10), linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.10), force = 3) + #, direction = "y", force = 0.5, nudge_x           = 0.15) +

  theme_classic(base_size = 18) +
  xlim(values = c(-4,4)) +

  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC of means", y = "-log2(Padj)")
ggsave("results/manuscript_new/facs/volcano_wlcx_io_naive.pdf", width = 6, height = 5)


p.df %>% filter(p.value < 0.05)
p.df %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% filter(p.adj < 0.1)
p.df %>% filter(!grepl(pattern = "NKT", x = variable)) %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% filter(p.adj < 0.1)



p.df %>% filter(p.adj < 0.05)

p.df <- fread("results/manuscript_new/facs/facs_wlcx_p_io_naive.txt")

p.df <- p.df %>%
  mutate(variable = gsub(pattern = "\\.", replacement = "\\ ", variable)) %>%
  mutate(variable = gsub(pattern = "pos", replacement = "+", variable)) %>%
  mutate(variable = gsub(pattern = "neg", replacement = "-", variable))

p.df <- p.df %>%
  mutate(variable = paste0(variable, "+")) %>%
  mutate(variable = gsub(pattern = "\\++", replacement = "\\+", variable)) %>%
  mutate(variable = gsub(pattern = "\\-\\+", replacement = "\\-", variable))

ggplot(p.df, aes(x = log2fc.mean, y = -log10(p.value), color = log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.05), force = 3, color = "black", size = 6) + #, direction = "y", force = 0.5, nudge_x           = 0.15) +

  theme_classic(base_size = 18) +
  xlim(values = c(-4,4)) +

  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC of means", y = "-log10(P)")
ggsave("results/manuscript_new/facs/volcano_wlcx_io_naive2.pdf", width = 6, height = 5)
fwrite(p.df, "results/manuscript_new/facs/facs_wlcx_p_io_naive.txt", sep = "\t", quote = F, row.names = F)

krskl.df <- krskl.df %>%
  mutate(variable = gsub(pattern = "\\.", replacement = "\\ ", variable)) %>%
  mutate(variable = gsub(pattern = "pos", replacement = "+", variable)) %>%
  mutate(variable = gsub(pattern = "neg", replacement = "-", variable)) %>%
  mutate(variable = paste0(variable, "+")) %>%
  mutate(variable = gsub(pattern = "\\++", replacement = "\\+", variable)) %>%
  mutate(variable = gsub(pattern = "\\-\\+", replacement = "\\-", variable))

krskl.df %>% filter(variable %in% subset(p.df, p.adj < 0.1)$variable) %>%
  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable, scales = "free_y", strip.position="right", ncol = 7) +
  ggpubr::stat_compare_means(label = "p", label.y.npc = 0.9) + scale_fill_manual(values = c("red3", "darkolivegreen4", "darkolivegreen4")) +
  theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + theme(strip.text.y = element_text(size = 10)) + ggpubr::rotate_x_text(angle = 45)
ggsave("results/manuscript_new/facs/box_io_naive_sigf.pdf", width = 12, height = 4)

krskl.df %>% filter(variable == "CD4 LAG3+") %>%
  mutate(timepoint = gsub("m", " month", x = timepoint)) %>%
  ggplot(aes(timepoint, value/100, fill = timepoint)) + geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label = "p", label.y.npc = 0.9) + scale_fill_manual(values = c("red3", "darkgreen", "darkgreen")) +
  theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + theme(strip.text.y = element_text(size = 10)) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "prop of repertoire", title = "CD4+\nLAG3+")
ggsave("results/manuscript_new/treg/box_cd4_lag3.pdf", width = 2, height = 4)




## Test wlcx for io resistant
df <- krskl.df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>%
  filter(Cohort == 2)
df <- df[!duplicated(df), ]

krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]
krskl.df$timepoint <- ifelse(krskl.df$timepoint == "0m", "prior", "post")
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  # message(x)
  y <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(timepoint) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()
  if(n.vars == 2){

    # t.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
    wilcox.test(value~timepoint, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, timepoint == "prior")$value), median.y = median(subset(y, timepoint != "prior")$value),
                                                                                                                mean.x   = mean(subset(y, timepoint == "prior")$value), mean.y   = mean(subset(y, timepoint != "prior")$value), log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))
    # coin::wilcoxsign_test(value~timepoint, data = y, distribution="exact") %>% broom::tidy() %>% mutate(variable = x)

  }

}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_wlcx_p_prior_io.txt", sep = "\t", quote = F, row.names = F)


ggplot(p.df, aes(x = log2fc.mean, y = -log2(p.adj), color = log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.10), force = 3) + #, direction = "y", force = 0.5, nudge_x           = 0.15) +

  theme_classic(base_size = 18) +
  xlim(values = c(-4,4)) +

  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC of means", y = "-log2(Padj)")
ggsave("results/manuscript_new/facs/volcano_wlcx_prior_io.pdf", width = 6, height = 5)



p.df <- fread("results/manuscript_new/facs/facs_wlcx_p_prior_io.txt")

p.df <- p.df %>%
  mutate(variable = gsub(pattern = "\\.", replacement = "\\ ", variable)) %>%
  mutate(variable = gsub(pattern = "pos", replacement = "+", variable)) %>%
  mutate(variable = gsub(pattern = "neg", replacement = "-", variable))

p.df <- p.df %>%
  mutate(variable = paste0(variable, "+")) %>%
  mutate(variable = gsub(pattern = "\\++", replacement = "\\+", variable)) %>%
  mutate(variable = gsub(pattern = "\\-\\+", replacement = "\\-", variable))

ggplot(p.df, aes(x = log2fc.mean, y = -log10(p.value), color = log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.05), force = 7, color = "black", size = 6) + #, direction = "y", force = 0.5, nudge_x           = 0.15) +

  theme_classic(base_size = 18) +
  xlim(values = c(-4,4)) +

  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC of means", y = "-log10(P)")
ggsave("results/manuscript_new/facs/volcano_wlcx_prior_io2.pdf", width = 6, height = 5)
fwrite(p.df, "results/manuscript_new/facs/facs_wlcx_p_prior_io.txt", sep = "\t", quote = F, row.names = F)
p.df <- fread("results/manuscript_new/facs/facs_wlcx_p_prior_io.txt")

krskl.df <- krskl.df %>%
  mutate(variable = gsub(pattern = "\\.", replacement = "\\ ", variable)) %>%
  mutate(variable = gsub(pattern = "pos", replacement = "+", variable)) %>%
  mutate(variable = gsub(pattern = "neg", replacement = "-", variable)) %>%
  mutate(variable = paste0(variable, "+")) %>%
  mutate(variable = gsub(pattern = "\\++", replacement = "\\+", variable)) %>%
  mutate(variable = gsub(pattern = "\\-\\+", replacement = "\\-", variable))

krskl.df %>% filter(variable %in% subset(p.df, p.adj < 0.1)$variable) %>%
  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~variable, scales = "free_y", strip.position="right", ncol = 5) +
  ggpubr::stat_compare_means(label = "p", label.y.npc = 0.9) + scale_fill_manual(values = c("red3", "darkolivegreen4", "darkolivegreen4")) +
  theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + theme(strip.text.y = element_text(size = 10)) + ggpubr::rotate_x_text(angle = 45)
ggsave("results/manuscript_new/facs/box_prior_io_sigf.pdf", width = 12, height = 4)











sigf_variables_naive <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)
p.df %>% filter(p.adj < 0.05)
krskl.df %>% filter(variable %in% sigf_variables_naive) %>% left_join(p.df %>% filter(p.adj < 0.1)) %>%
  mutate(timepoint = factor(as.character(timepoint), levels = c("prior", "post"))) %>%


  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) +
  # facet_grid(~reorder(variable, p.adj), scales = "free_y") +
  facet_grid(~variable) +

  # geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)[c(3,2,1)]) + geom_jitter(size = 0.5) + facets_nice +
  ggpubr::stat_compare_means(label="p") + theme(legend.position = "none") #+ labs("")
ggsave("results/manuscript_new/facs/box_krskl_p.pdf", width = 6, height = 5)







## Test krskl for io naive
df <- krskl.df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>%
  filter(Cohort == 1)
df <- df[!duplicated(df), ]

krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_naive <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)
p.df %>% filter(p.adj < 0.05)
krskl.df %>% filter(variable %in% sigf_variables_naive) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>%


  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~reorder(variable, p.adj), scales = "free_y", ncol = 2) +
  # geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(4)[c(3,2,1)]) + geom_jitter(size = 0.5) + facets_nice +
  ggpubr::stat_compare_means(label="p") + theme(legend.position = "none") #+ labs("")
ggsave("results/manuscript_new/facs/box_krskl_p.pdf", width = 6, height = 5)





## Test krskl for io prior
df <- krskl.df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>%
  filter(Cohort == 2)
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_prior.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_prior <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)

p.df %>% filter(p.adj < 0.01)

krskl.df %>% filter(variable %in% sigf_variables_prior) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>%

  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~reorder(variable, p.adj), scales = "free_y", ncol = 3) +
  # geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette5(4)[c(3,2,1)]) + geom_jitter(size = 0.5) + facets_nice +
  ggpubr::stat_compare_means(label="p") + theme(legend.position = "none") #+ labs("")
ggsave("results/manuscript_new/facs/fbox_krskl_p_prior.pdf", width = 8, height = 5)

sigf_variables_prior %in% sigf_variables_naive








## Test krskl for cmv
lag3_clin <- fread("data/facs/lag3_tot.csv")
colnames(lag3_clin) <- make.names(colnames(lag3_clin))
lag3_cmv <- lag3_clin %>% dplyr::select(FM,CMV.seropositive) %>% dplyr::rename(name=FM) #%>% mutate(name=as.character(name))
facs_tot <- facs_tot %>% left_join(lag3_cmv)

df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(!is.na(CMV.seropositive)) %>%
  filter(timepoint != "6m") %>%
  filter(CMV.seropositive == "CMVpos") %>%
  filter(Cohort == 1)
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_io_naive_cmvpos.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_naive_cmv_pos <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)

## cmv neg
df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(!is.na(CMV.seropositive)) %>%
  filter(timepoint != "6m") %>%
  filter(CMV.seropositive == "CMVneg") %>%
  filter(Cohort == 1)
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_io_naive_cmvneg.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_naive_cmv_neg <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)



## prior io cmv pos
df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(!is.na(CMV.seropositive)) %>%
  filter(timepoint != "6m") %>%
  filter(CMV.seropositive == "CMVpos") %>%
  filter(Cohort == 2)
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_prior_io_cmvpos.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_prior_cmv_pos <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)


## cmv neg
df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(!is.na(CMV.seropositive)) %>%
  filter(timepoint != "6m") %>%
  filter(CMV.seropositive == "CMVneg") %>%
  filter(Cohort == 2)
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_prior_io_cmvneg.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_prior_cmv_neg <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)


krskl.df %>% filter(variable %in% sigf_variables_prior_cmv_pos) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>%

  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~reorder(variable, p.adj), scales = "free_y", ncol = 3) +
  # geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette5(4)[c(3,2,1)]) + geom_jitter(size = 0.5) + facets_nice +
  ggpubr::stat_compare_means(label="p") + theme(legend.position = "none") #+ labs("")
ggsave("results/manuscript_new/facs/fbox_krskl_p_prior.pdf", width = 8, height = 5)






## Test krskl for cmv pos
df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(!is.na(CMV.seropositive)) %>%
  filter(timepoint != "6m") %>%
  filter(CMV.seropositive == "CMVpos")
df <- df[!duplicated(df), ]
krskl.df <- df %>% dplyr::select(timepoint,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_krskl_p_io_naive_cmvpos.txt", sep = "\t", quote = F, row.names = F)





## Test mann-whitney for io naive cmv pos vs neg at base line
df <- facs_tot %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>%
  filter(timepoint == "0m") %>%
  filter(Cohort == 1)
df <- df[!duplicated(df), ]

krskl.df <- df %>% dplyr::select(CMV.seropositive,CD56bright:GD.CX3CR1, CD3:CD8.EMRA) %>% melt(id="CMV.seropositive")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$CMV.seropositive)) == 2){
    wilcox.test(value~CMV.seropositive, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_cmv_io_naive_baseline.txt", sep = "\t", quote = F, row.names = F)

sigf_variables_naive_cmv <- p.df %>% filter(p.adj < 0.05)
sigf_variables_naive_cmv <- p.df %>% filter(p.adj < 0.05) %>% pull(variable)

levels <- c("CD4negCD8pos.CX3CR1", "CD4posCD8neg.CX3CR1", "NKTdim.CX3CR1", "CD4.CD57",
  "CD8.naive", "CD8.CD27", "NKTdim.CD27", "NKTbright.CD27")

krskl.df %>% filter(variable %in% sigf_variables_naive_cmv) %>% left_join(p.df %>% filter(p.adj < 0.05)) %>%
  filter(variable %in% levels) %>%
  mutate(variable = factor(as.character(variable), levels = levels)) %>%
  mutate(variable = gsub("CD8neg", "", variable)) %>%
  mutate(variable = gsub("CD4neg", "", variable)) %>%

  ggplot(aes(CMV.seropositive, value, fill = CMV.seropositive)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~reorder(variable, p.adj), scales = "free_y", ncol = 4) +
  # geom_vline(xintercept = 1.5, linetype = "dotted") + geom_vline(xintercept = 2.5, linetype = "dotted") +
  ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = c("lightgrey", "antiquewhite")) + geom_jitter(size = 0.5) + facets_nice +
  ggpubr::stat_compare_means(label="p") + theme(legend.position = "none") + labs(x="")
ggsave("results/manuscript_new/facs/box_krskl_p_cmv_io_naive.pdf", width = 8, height = 5)



## the same but with volcano
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  # message(x)
  y <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(CMV.seropositive) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()
  if(n.vars == 2){
      wilcox.test(value~CMV.seropositive, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, CMV.seropositive == "CMVpos")$value), median.y = median(subset(y, CMV.seropositive != "CMVpos")$value),
                                                                                                                mean.x   = mean(subset(y, CMV.seropositive == "CMVpos")$value), mean.y   = mean(subset(y, CMV.seropositive != "CMVpos")$value), log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))
  }

}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/facs/facs_cmv_io_naive_baseline2.txt", sep = "\t", quote = F, row.names = F)

p.df <- fread("results/manuscript_new/facs/facs_cmv_io_naive_baseline2.txt")

p.df %>%
  filter(!grepl("NKT", variable)) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>% filter(p.adj < 0.05)


ggplot(p.df, aes(x = -log2fc.mean, y = -log2(p.adj), color = -log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +
  # ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.05), force = 5, size = 3) +
  ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.05),
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20, size = 3
  ) +

  theme_classic(base_size = 18) +
  xlim(values = c(-4,4)) + ylim(c(0,6)) +

  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average log2FC of means", y = "-log2(Padj)")
ggsave("results/manuscript_new/facs/volcano_wlcx_cmv_baseline_io_naive.pdf", width = 6, height = 5)
