
## Run the preprocessing pipeline
setwd("/csc/mustjoki/lag3/") ## FIMM

source("src/R/main.R")
source("src/R/preprocess/run_preprocessRNA.R")
source("src/R/preprocess/run_preprocessTCRab.R")
source("src/R/preprocess/run_preprocessRNATCRab.R")

message("Fin")
