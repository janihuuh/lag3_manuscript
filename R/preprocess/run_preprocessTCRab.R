
message("Preprocess TCRab...")

## Read in the sequencing results; the most important files are all_contig and clonotypes -files
contig_files    <- list.files("data/scRNAseq+TCRseq/", recursive = T, pattern = "all_contig_annotations.csv", full.names = T)
clonotype_files <- list.files("data/scRNAseq+TCRseq/", recursive = T, pattern = "clonotypes.csv", full.names = T)

contigs    <- lapply(contig_files, read.delim, sep = ",")
clonotypes <- lapply(clonotype_files, read.delim, sep = ",")

## Pre-process the data. The function produces two files, barcoded and clonotyped
for(i in 1:length(clonotypes)){

  name <- substr(clonotype_files[i], 1, nchar(clonotype_files[i]) - 15) %>% extractFileName()

  if(!name %in% gsub("_barcoded.txt", "", extractFileName(barcoded_files))){
  # if(name %in% c("2351_0", "2351_1")){

    message(name)
    preprocess10Xtcr(contig_file = contigs[[i]], 
                       clonotype_file = clonotypes[[i]], 
                       prefix = paste0("data/scRNAseq+TCRseq/preprocessed/", name))

  }

}


## Read in the just processed files
clonotyped_files <- list.files("data/scRNAseq+TCRseq/preprocessed/", full.names = T, pattern = "clonotyped")
barcoded_files   <- list.files("data/scRNAseq+TCRseq/preprocessed/", full.names = T, pattern = "barcoded")

clonotyped <- lapply(clonotyped_files, fread)
barcoded   <- lapply(barcoded_files, function(x){fread(x) %>% mutate(barcode_uniq = paste0(gsub("_barcoded.txt", "", extractFileName(x)), "_", barcode))}) %>% rbindlist()

## Make new clonotype id:s for each of the patient profiled
barcoded$patient <- extractName(barcoded$barcode_uniq)
barcoded         <- barcoded %>% mutate(patient = substr(barcode_uniq, 1, 4))
barcoded_per_pt  <- split(barcoded, f = barcoded$patient)

tot_barcode <- lapply(barcoded_per_pt, function(x){message(unique(x$patient)); newClonotype_df(x)}) %>% rbindlist()
tot_barcode <- tot_barcode %>%
  mutate(new_clonotypes_id = paste0(tot_barcode$patient, "_", tot_barcode$new_clonotypes_id)) %>%
  mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
tot_barcode <- tot_barcode %>% mutate(barcode_uniq = gsub("\\-", "_", barcode_uniq))


## Write down results
write.table(tot_barcode, "data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt",  sep = "\t", row.names = F, quote = F)
tot_barcode[!duplicated(tot_barcode$new_clonotypes_id), ] %>%
  write.table("data/scRNAseq+TCRseq/preprocessed/tot_uniq_tcrab.txt", sep = "\t", row.names = F, quote = F)


tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt") 



