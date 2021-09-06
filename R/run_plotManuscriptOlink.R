
dir.create("results/manuscript_new/olink/")

## Test krskl for io naive
df <- olink_pca %>% dplyr::rename(timepoint=Timepoint)
df <- df[!duplicated(df), ]

krskl.df <- df %>% dplyr::select(timepoint,IL8:CSF1) %>% melt(id="timepoint")
krskl.df <- krskl.df[!duplicated(krskl.df), ]

p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(variable == x)
  if(length(unique(y$timepoint)) == 3){
    kruskal.test(value~timepoint, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/olink/olink_krskl_p_naive.txt", sep = "\t", quote = F, row.names = F)


olink_pca %>%
  mutate(timepoint = factor(as.character(Timepoint), levels = c("D0", "1mo", "3mo"))) %>%
#   plyr::revalue(timepoint, replace = c("D0" = "0mo")) %>%
  ggplot(aes(timepoint, PDCD1, fill = timepoint)) + geom_boxplot() + scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + ggpubr::stat_compare_means(label="p")
ggsave("results/manuscript_new/olink/box_io_naive.pdf", width = 5, height = 4)






############################################################

## Test mann-whitney
olink    <- fread("data/olink/olink_IO_naive.txt")
krskl.df <- olink %>% reshape::melt(id = "timepoint")

## the same but with volcano
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){

  message(x)
  y      <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(timepoint) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()

  if(n.vars == 2){
    wilcox.test(value~timepoint, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, timepoint == "pre")$value), median.y = median(subset(y, timepoint != "pre")$value),
                                                                                                                mean.x   = mean(subset(y, timepoint == "pre")$value), mean.y   = mean(subset(y, timepoint != "pre")$value),
                                                                                                                log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))
  }

}) %>% rbindlist() %>% arrange(p.value) %>% mutate(p.adj = p.adjust(p.value, method = "BH"))
fwrite(p.df, "results/manuscript_new/olink/p_olink_io_naive.txt", sep = "\t", quote = F, row.names = F)
# p.df <- fread("results/manuscript_new/olink/p_olink_io_naive.txt")

ggplot(p.df, aes(x = log2fc.mean, y = -log2(p.value), color = log2fc.mean, label = variable)) + geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +

  ggrepel::geom_text_repel(data = subset(p.df, p.value < 0.05),
                           nudge_x = .15,
                           force = 5,
                           box.padding = 0.5,
                           nudge_y = 1,
                           segment.curvature = -0.1,
                           segment.ncp = 3,
                           segment.angle = 20, size = 3, color = "black"
  ) +

  theme_classic(base_size = 18) + xlim(values = c(-0.75,0.75)) +

  # scale_colour_gradient(low = "darkblue", high = "darkred", limits = c(-0.2,0.8)) +
  scale_colour_gradient(low = "dodgerblue", high = "darkred", limits = c(-0.2,0.8)) +
  theme(legend.position = "right") + labs(x = "average log2FC of means", y = "-log2(P)",)
ggsave("results/manuscript_new/olink/volcano_olink_io_naive.pdf", width = 6, height = 4)





## Test mann-whitney
olink <- fread("data/olink/olink_IO_naive.txt") %>% dplyr::select(-NOS3, -FGF2, -CD28)
length(unique(olink$FM))

krskl.df <- olink %>%
  filter(!Timepoint %in% c("", "healthy")) %>%
  mutate(timepoint = factor(as.character(Timepoint), levels = c("healthy", "D0", "1mo", "3mo"))) %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>% dplyr::select(-FM,-Timepoint) %>%
  mutate(timepoint = ifelse(timepoint == "D0", "pre", "post"))

krskl.df <- krskl.df[!duplicated(krskl.df), ]
krskl.df <- krskl.df %>% melt(id = "timepoint")

## the same but with volcano
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  # message(x)
  y <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(timepoint) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()
  if(n.vars == 2){
    wilcox.test(value~timepoint, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, timepoint == "pre")$value), median.y = median(subset(y, timepoint != "pre")$value),
                                                                                                                       mean.x   = mean(subset(y, timepoint == "pre")$value), mean.y   = mean(subset(y, timepoint != "pre")$value), log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))
  }

}) %>% rbindlist() %>% arrange(p.value) %>% mutate(p.adj = p.adjust(p.value, method = "BH"))# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/olink/p_olink_io_naive.txt", sep = "\t", quote = F, row.names = F)
p.df <- fread("results/manuscript_new/olink/p_olink_io_naive.txt")

ggplot(p.df, aes(x = log2fc.mean, y = -log2(p.value), color = -log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +
  # ggrepel::geom_text_repel(data = subset(p.df, p.adj < 0.05), force = 5, size = 3) +
  ggrepel::geom_text_repel(data = subset(p.df, p.value < 0.05),
                           nudge_x = .15,
                           box.padding = 0.5,
                           nudge_y = 1,
                           segment.curvature = -0.1,
                           segment.ncp = 3,
                           segment.angle = 20, size = 3
  ) +

  theme_classic(base_size = 18) + xlim(values = c(-0.75,0.75)) +

  scale_colour_gradient(low = "darkblue", high = "darkred", limits = c(-1,1)) + theme(legend.position = "none") + labs(x = "average log2FC of means", y = "-log2(P)")
ggsave("results/manuscript_new/olink/volcano_olink_io_naive.pdf", width = 5, height = 4)



olink %>%
  mutate(timepoint = factor(as.character(Timepoint), levels = c("healthy", "D0", "1mo", "3mo"))) %>%

  filter(!is.na(timepoint)) %>% melt(id = c("Timepoint", "timepoint", "FM")) %>%
  filter(variable %in% subset(p.df, p.value < 0.05)$variable) %>%


  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5, alpha = 0.5) +
  theme_classic(base_size = 17) +
  scale_fill_manual(values = c("gray50", "red3", "darkolivegreen4", "darkolivegreen4")) + theme(legend.position = "none") +
  labs(x = "") + ggpubr::rotate_x_text(45) + ggpubr::stat_compare_means(label = "p", label.y.npc = 0.85) + facet_wrap(~variable, scales = "free_y", ncol = 4)
ggsave("results/manuscript_new/olink/box_p05_naive.pdf", width = 7, height = 9)













## Test mann-whitney for prior
olink_prior <- fread("data/olink/olink_prior_io.txt") %>% dplyr::select(-NOS3, -FGF2, -CD28)
table(olink_prior$Timepoint)

dim(olink_prior)
length(unique(olink_prior$FM))


length(unique(subset(olink_prior, Timepoint != "no")$FM))

olink_prior_clin <- olink_prior %>% dplyr::select(FM:Response)
olink_prior <- olink_prior %>% dplyr::select(-c(FM:Response))
olink_prior <- apply(olink_prior, 2, FUN = function(x) gsub(x, pattern = "\\,", replacement = "\\.") %>% as.character %>% as.numeric)  %>% as.data.frame()
olink <- olink_prior %>% bind_cols(dplyr::select(olink_prior_clin, FM, Timepoint))

krskl.df <- olink %>%
  filter(!Timepoint %in% c("", "healthy")) %>%
  mutate(timepoint = factor(as.character(Timepoint), levels = c("healthy", "D0", "1mo", "3mo"))) %>%
  filter(!is.na(timepoint)) %>%
  filter(timepoint != "6m") %>% dplyr::select(-FM,-Timepoint) %>%
  mutate(timepoint = ifelse(timepoint == "D0", "pre", "post"))

krskl.df <- krskl.df[!duplicated(krskl.df), ]
krskl.df <- krskl.df %>% melt(id = "timepoint")

## the same but with volcano
p.df <- lapply(unique(krskl.df$variable), FUN = function(x){
  # message(x)
  y <- krskl.df %>% filter(variable == x) %>% filter(!is.na(value))
  n.vars <- y %>% group_by(timepoint) %>% summarise(n=n(), verbose = T) %>% filter(n>1) %>% nrow()
  if(n.vars == 2){
    wilcox.test(value~timepoint, data = y, conf.int=TRUE) %>% broom::tidy() %>% mutate(variable = x) %>% mutate(median.x = median(subset(y, timepoint == "pre")$value), median.y = median(subset(y, timepoint != "pre")$value),
                                                                                                                mean.x   = mean(subset(y, timepoint == "pre")$value), mean.y   = mean(subset(y, timepoint != "pre")$value), log2fc.median = log2(median.y/median.x), log2fc.mean = log2(mean.y/mean.x))
  }

}) %>% rbindlist() %>% arrange(p.value) %>% mutate(p.adj = p.adjust(p.value, method = "BH"))# %>% left_join(var_df)
fwrite(p.df, "results/manuscript_new/olink/p_olink_prior_io.txt", sep = "\t", quote = F, row.names = F)
p.df <- fread("results/manuscript_new/olink/p_olink_prior_io.txt")

ggplot(p.df, aes(x = log2fc.mean, y = -log2(p.value), color = -log2fc.mean, label = variable)) + geom_point(size = 3) +
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(p.df, p.value < 0.1), force = 5, size = 3) +
  # ggrepel::geom_text_repel(data = subset(p.df, p.value < 0.1),
  #                          nudge_x = .15,
  #                          box.padding = 0.5,
  #                          nudge_y = 1,
  #                          segment.curvature = -0.1,
  #                          segment.ncp = 3,
  #                          segment.angle = 20, size = 3) +

  theme_classic(base_size = 18) + xlim(values = c(-0.2,0.2)) +

  scale_colour_gradient(low = "darkblue", high = "darkred", limits = c(-1,1)) + theme(legend.position = "none") + labs(x = "average log2FC of means", y = "-log2(P)")
ggsave("results/manuscript_new/olink/volcano_olink_prior_io.pdf", width = 5, height = 4)



olink %>% mutate(Timepoint = ifelse(Timepoint == "no", "healthy", Timepoint)) %>%
  mutate(timepoint = factor(as.character(Timepoint), levels = c("healthy", "D0", "1mo", "3mo"))) %>%
  filter(!is.na(timepoint)) %>% melt(id = c("Timepoint", "timepoint", "FM")) %>%
  filter(variable %in% subset(p.df, p.value < 0.10)$variable) %>%


  ggplot(aes(timepoint, value, fill = timepoint)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5, alpha = 0.5) +
  theme_classic(base_size = 17) +
  scale_fill_manual(values = c("gray50", "red3", "darkolivegreen4", "darkolivegreen4")) + theme(legend.position = "none") +
  labs(x = "") + ggpubr::rotate_x_text(45) + ggpubr::stat_compare_means(label = "p", label.y.npc = 0.85) + facet_wrap(~variable, scales = "free_y")
ggsave("results/manuscript_new/olink/box_p010_prior_io.pdf", width = 5, height = 4)
