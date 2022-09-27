# Copyright 2022, Michael Schneider, All rights reserved.
## collect QDNAseq objects or data frames for single cell profiles into single object
# two examples usage cases
# Rscript merge.R target.rds input_folder file_in_folder
# Rscript merge.R target.rds input_1.rds input_2.rds input_3.rds

require(BiocGenerics, quietly=TRUE)
require(gtools, quietly=TRUE)
require(dplyr, quietly=TRUE)
#require(scAbsolute)
source("~/scAbsolute/R/core.R")
args = commandArgs(trailingOnly=TRUE)
indexSort = TRUE
options(future.globals.maxSize= 4096*1024^2)

print(args)
rdsFiles = tail(args, n=-1)
if (length(args) >= 3){
  d = args[2]
  f = args[3]  
}else{
  d = args[2]
  f = NULL
}

if (dir.exists(d) && (is.null(f) || !file.exists(f))){
  # merge based on filename ending (f) in folder (d) to file args[1]
  setwd(d)
  rdsFiles = list.files(pattern=paste0(f, ".rds$"), recursive=FALSE)
}else{
  # merge all files in args[2:end] to file args[1]
  stopifnot(all(endsWith(rdsFiles, ".rds")))
  stopifnot(all(file.exists(rdsFiles)))
  
  if(indexSort){
    rdsFiles = sort(rdsFiles, decreasing=FALSE)
  }
}

if (is.data.frame(readRDS(rdsFiles[[1]]))) {
  rdsData = base::lapply(rdsFiles, readRDS)#, USE.NAMES = FALSE)
  mergedRDS = do.call(rbind, rdsData)
} else {
  rdsData = base::sapply(rdsFiles, readRDS, USE.NAMES = FALSE)  
  # alternative command, but issues with NaN values in pData: 
  #mergedRDS = do.call(BiocGenerics::combine, rdsData)
  mergedRDS = combineQDNASets(unlist(rdsData))
}

saveRDS(mergedRDS, file = args[1])
