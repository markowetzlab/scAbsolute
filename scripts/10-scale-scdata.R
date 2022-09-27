## scAbsolute - parallel run script

## Variables ####
args = commandArgs(trailingOnly=TRUE)

if (interactive()){
  rm(list=ls()) 
  
  # bamPath = "~/Data/project1/20-align/"
  # bamFile = "UID-10X-Andor-2020-48959-MKN-45_SLX-00000_000365_CTAATGGGTTGATCGT-1.bam"
  #bamPath = "~/Data/project1/20-align/UID-BUG/"
  bamPath = "~/Data/Figures/scAbsolute/rawdata/"
  bamFile = c("UID-DLP-SA1044_SLX-A96139A_000007_R03-C17.bam",
              "UID-DLP-SA1044_SLX-A96139A_000057_R04-C47.bam")

  RESULTPATH = "~/ExampleOutput.rds"
  
  # options
  binSize = 500
  
  minPloidy = 1.1
  maxPloidy = 8.0

  reticulate::use_condaenv(condaenv = "rstudio-server", conda = "~/.anaconda3/bin/conda")
  BASEDIR="~/scAbsolute"
}else{
  bamPath = args[1]
  bamFile = args[2]
  
  RESULTPATH = args[3]

  binSize = args[4]
  
  if(length(args) == 6){
    minPloidy = as.numeric(args[5])
    maxPloidy = as.numeric(args[6])
  }else{
    minPloidy = NULL
    maxPloidy = NULL
  }
  
  BASEDIR="/opt/scAbsolute"
}

# PACKAGE DEPENDENCIES
library(reticulate, quietly=TRUE, warn.conflicts = FALSE)
library(future.apply, quietly=TRUE, warn.conflicts = FALSE)
suppressWarnings(plan(multiprocess))
library(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)

library(tidyverse, quietly=TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts = FALSE)

library(Biobase, quietly=TRUE, warn.conflicts = FALSE)
library(BiocGenerics, quietly=TRUE, warn.conflicts = FALSE)
library(devtools, quietly=TRUE, warn.conflicts = FALSE)
library(digest, quietly=TRUE, warn.conflicts = FALSE)
library(IRanges, quietly=TRUE, warn.conflicts = FALSE)
library(MASS, quietly=TRUE, warn.conflicts = FALSE)
library(robustbase, quietly=TRUE, warn.conflicts = FALSE)
library(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
library(matrixStats, quietly=TRUE, warn.conflicts = FALSE)
print(paste0("Base directory: ", BASEDIR))

source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
source(file.path(BASEDIR, "R/scSegment.R"))
source(file.path(BASEDIR, "R/scAbsolute.R"))
source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/mean-variance.R"))
source(file.path(BASEDIR, "R/visualization.R"))
source(file.path(BASEDIR, "R/cellcycle.R"))
# END DEPENDENCIES

## PARAMETERS ====
set.seed(2020)
randomSeed = 2020
bamFile = do.call("c", (base::strsplit(bamFile, split=",")))

# these variables are script arguments
binSize = as.numeric(binSize)
# bin sizes currently supported by QDNAseq framework
stopifnot(binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000))

species = "Human"
genome = "hg19"

trimLength = 1
minLength = 10
globalModel = NULL
# optimization of segmentation
testStatistic = "NegBinMM"
splitPerChromosome=TRUE
optimizeSegmentation=FALSE

max_iterations=101
change_prob=1e-3
max_states=9

## Save data for later processing
if(!dir.exists(dirname(RESULTPATH))) dir.create(dirname(RESULTPATH), recursive = TRUE)
hmm_path = base::dirname(RESULTPATH)

# ploidy description, example how to specify this on a per cell level
# cellnames = c("UID-10X-COLO829_SLX-00000_000326_ATGTGTGCAAAGTAAC-1", "UID-10X-Fibroblast-cell_SLX-00000_000277_CACACCTTCACACGTA-1")
# NOTE that order has to be correct, minPloidy entry 1 is cell one, entry 2 is cell 2, ...
# minPloidy = list(1.7, 5.1); names(minPloidy) = cellnames
# maxPloidy = list(3.4, 6.9); names(maxPloidy) = cellnames
# Default ranges without any prior information
if(is.null(minPloidy)) minPloidy = 1.1
if(is.null(maxPloidy)) maxPloidy = 8.0

ploidyWindow = 0.1
ploidyRegion=c(paste0("chr", as.character(seq(1,22))), "chrX")
selectRegion=c(paste0("chr", as.character(seq(1,22))), "chrX", "chrY")
limitPloidy=16

method = "error" # model or error
## END PARAMETERS

# load correct data
filePaths = gsub("//", "/", c(paste(bamPath, bamFile, sep="/")))
filePaths = gsub("//", "/", filePaths)
stopifnot(all(file.exists(filePaths)) || dir.exists(bamPath))

if(all(file.exists(filePaths)) && all(utils::file_test("-f", filePaths))){
} else {
  if(any(is.null(bamFile) || length(bamFile) == 0 || bamFile == "")){
    print(paste0("Collecting all bam files in folder", bamPath))
    filePaths = file.path(bamPath, list.files(path=bamPath, pattern = "\\.bam$"))
  }
}
stopifnot(all(file.exists(filePaths)))


## run the core algorithm - using scAbsolute wrapper function
scaledCN = scAbsolute(filePaths, method=method, globalModel=globalModel, binSize=binSize, genome=genome,
           minLength=minLength, batchSize=1, testStatistic=testStatistic, limitPloidy=limitPloidy,
           minPloidy=minPloidy, maxPloidy=maxPloidy, ploidyWindow=ploidyWindow, ploidyRegion=ploidyRegion, selectRegion=selectRegion,
           splitPerChromosome=splitPerChromosome, optimizeSegmentation=optimizeSegmentation,
           max_iterations=max_iterations, hmm_path=hmm_path, change_prob=change_prob, max_states=max_states,
           randomSeed=randomSeed, outputPath=RESULTPATH, outputSegmentation=FALSE, debug=TRUE)

## Debugging and examples
# 
# readCounts = readData(filePaths, binSize=binSize, extendedBlacklisting=TRUE,
#                       filterChromosomes = c("MT"), genome="hg19")
# 
# 
# plotCounts(readCounts[100000:102954,2])
# 
# ## check copynumber profile
# cell = 1
# plotCopynumber(scaledCN[, cell], ylim=c(0, 15))
# 
# ## check fit
# cell = 1
# p = plotFit(scaledCN[,cell], scale=TRUE)
#
# ## write output to file
# gr = getSegTable(scaledCN)
# writeSegTable(scaledCN, filename="~/test.bed", trimLength=trimLength, minLength=minLength)

## store metadata ====

# NOTE add whole genome values (in case of restricted ploidy calculation via ploidyRegion)
scaledCN = computeRPC(scaledCN, rpc_value="rpc.all", ploidy_value="ploidy.all")

# compute gini coefficient (quality control)
scaledCN = computeGini(scaledCN)

# compute MAPD coefficient (quality control)
scaledCN = computeMAPD(scaledCN)

# compute information theoretic measures (quality control)
scaledCN = computeInfotheo(scaledCN)

# compute l2 segmentation error
scaledCN = computeL2(scaledCN)


## Include metadata
sample_names = rownames(pData(scaledCN))
description <- pData(scaledCN)

description$minLength = minLength
description$method = method
description$globalModel = globalModel
description$genome = genome

unpack <- function(object, x){
  if(!is.null(names(x))){
    l = c()
    for(n in colnames(object)){
      l = c(x[[n]], l)
    }
    return(l)
  }else{
    return(x)
  }
}
description$rpc.post = apply(scaledCN@assayData$calls[binsToUseInternal(scaledCN),,drop=FALSE], 2, sum) / apply(scaledCN@assayData$copynumber, 2, mean, na.rm=TRUE) / sum(binsToUseInternal(scaledCN))
description$ploidy.post = apply(scaledCN@assayData$copynumber[binsToUseInternal(scaledCN),,drop=FALSE], 2, mean, na.rm=TRUE)
description$changeProb = change_prob
description$minPloidy = unpack(scaledCN, minPloidy)
description$maxPloidy = unpack(scaledCN, maxPloidy)
description$ploidyWindow = ploidyWindow
description$ploidyRegion =  paste(ploidyRegion, collapse="-")
description$selectRegion =  paste(selectRegion, collapse="-")
description$limitPloidy = limitPloidy

description$testStatistic = testStatistic
description$optimizeSegmentation = optimizeSegmentation
description$splitPerChromosome = splitPerChromosome

description$n_bins = dim(scaledCN)[[1]]
description$n_effective_bins = dim(scaledCN)[[1]] - sum(!binsToUseInternal(scaledCN))

if(!is.null(globalModel)){
  description$globalModel_md5sum <- digest::digest(globalModel, algo="md5", serialize=F)
}else{
  description$globalModel_md5sum = NA
}

## Include quality control data
flagstat_description = readFlagstat(filePaths)

## Join descriptive data
if(!is.null(flagstat_description)){
  flagstat_description$name = as.character(flagstat_description$name)
  description = dplyr::left_join(description, flagstat_description, by="name")
  rownames(description) = sample_names
}
pData(scaledCN) = description

## Save data for later processing
if(!dir.exists(dirname(RESULTPATH))) dir.create(dirname(RESULTPATH), recursive = TRUE)
saveRDS(scaledCN, file=RESULTPATH)

## Reproducibility info
devtools::session_info()
