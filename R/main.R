
me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/lag3/"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(data.table)
library(gridExtra)
library(RColorBrewer)

theme_set(theme_classic(base_size = 17))

## Run all fun_* codes
for(code in list.files("src/R/", "fun", full.names = T, recursive = T)){

  message(code)
  source(code)

}

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
cytolytic <- c("GZMA", "GZMH", "GZMM", "PRF1", "GNLY")
facets_nice2 <- theme(strip.text.x = element_text(size=15, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))
