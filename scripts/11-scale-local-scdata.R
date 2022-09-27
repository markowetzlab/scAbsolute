# Copyright 2022, Michael Schneider, All rights reserved.
## scAbsolute - parallel run script

## Variables ####
args = commandArgs(trailingOnly=TRUE)

if (interactive()){
  rm(list=ls()) 
  
  # bamPath = "~/Data/project1/20-align/UID-10X-COLO829/SLX-00000/"
  # bamFile = c("UID-10X-COLO829_SLX-00000_000001_AAACCTGAGGCTACGA-1.bam",
  #             "UID-10X-COLO829_SLX-00000_000006_AAACGGGCAAAGTAAC-1.bam")
  bamPath = "~/Data/project1/20-align/UID-SIM-ploidy-Bulk-300/SLX-00000/"
  bamFile = c("UID-SIM-ploidy-Bulk-300_SLX-00000_000019_NNNNN.bam",
              "UID-SIM-ploidy-Bulk-300_SLX-00000_000089_NNNNN.bam")
    # "UID-10X-Spikein-1pc_SLX-00000_000515_CTCGAAAAGTCGATAA-1.bam",
    #           "UID-10X-Fibroblast-cell_SLX-00000_000914_TCATTTGTCTGGTGGC-1.bam",
    #           "UID-10X-COLO829_SLX-00000_000326_ATGTGTGCAAAGTAAC-1.bam",
    #           "UID-10X-MKN45-Gastric_SLX-00000_000493_ACGCCGACAAACTCAC-1.bam")
  # bamFile = "UID-10X-Spikein-1pc_SLX-00000_000515_CTCGAAAAGTCGATAA-1.bam"
  RESULTPATH = "~/UncertainDecember.rds"
  
  # options
  binSize = 500
  minLength = NULL
  
  minPloidy = 2.5
  maxPloidy = 3.5

  # reticulate::use_condaenv(condaenv = "rstudio-server", conda = "/home/schnei01/.anaconda3/bin/conda")
}else{
  bamPath = args[1]
  bamFile = args[2]
  
  RESULTPATH = args[3]

  binSize = args[4]
  
  if(length(args) == 5){
    minLength = args[5]  
  }else{
    minLength = NULL
  }
  
  if(length(args) == 6){
    minPloidy = as.numeric(args[5])
    maxPloidy = as.numeric(args[6])
  }else{
    minPloidy = NULL
    maxPloidy = NULL
  }
  
}

# TODO BIOCONDUCTOR
library(reticulate, quietly=TRUE, warn.conflicts = FALSE)
use_python("/opt/conda/bin/conda/python3")
# use_miniconda(condaenv = "base")
# reticulate::use_condaenv(condaenv = "base", conda = "/opt/conda/bin/conda")
# .libPaths( c( .libPaths(), "/home/schnei01/.R", "/home/schnei01/.local/share/r-miniconda/envs/r-reticulate") )
# .libPaths()
py_discover_config()

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

source("~/scAbsolute/data/changepoint/wrap_PELT.R")
source("~/scAbsolute/R/scSegment.R")
source("~/scAbsolute/R/scAbsolute.R")
source("~/scAbsolute/R/core.R")
source("~/scAbsolute/R/mean-variance.R")
source("~/scAbsolute/R/visualization.R")
source("~/scAbsolute/R/cellcycle.R")
# library(scAbsolute)
# END TODO BIOCONDUCTOR

## PARAMETERS ====
set.seed(2020)
randomSeed = 2020
bamFile = do.call("c", (base::strsplit(bamFile, split=",")))

# these variables are script arguments
binSize = as.numeric(binSize)
# bin sizes currently supported by QDNAseq framework
stopifnot(binSize %in% c(1, 5, 10, 15, 30, 50, 100, 500, 1000))

if(is.null(minLength)){
  idx = which(binSize == c(10, 15, 30, 50, 100, 500, 1000))
  minLengths = c(50, 40, 30, 25, 20, 15, 12)
  minLength = minLengths[[idx]]
}else{
  minLength = as.numeric(minLength)
}

species = "Human"
genome = "hg19"
# if(binSize <= 100){
#   splitPerChromosome=TRUE
# }else{
#   splitPerChromosome=FALSE
# }

trimLength = 1
# optimization of segmentation
testStatistic = "NegBinMM"
penalty="AIC"
splitPerChromosome=FALSE
optimizeSegmentation=FALSE

# ploidy description, example how to specify this on a per cell level
# cellnames = c("UID-10X-COLO829_SLX-00000_000326_ATGTGTGCAAAGTAAC-1", "UID-10X-Fibroblast-cell_SLX-00000_000277_CACACCTTCACACGTA-1")
# NOTE that order has to be correct, minPloidy entry 1 is cell one, entry 2 is cell 2, ...
# minPloidy = list(1.7, 5.1); names(minPloidy) = cellnames
# maxPloidy = list(3.4, 6.9); names(maxPloidy) = cellnames
# Default ranges without any prior information
if(is.null(minPloidy)) minPloidy = 1.1
if(is.null(maxPloidy)) maxPloidy = 12.0

ploidyWindow = 0.1
ploidyRegion=c(paste0("chr", as.character(seq(1,22))), "chrX")
selectRegion=c(paste0("chr", as.character(seq(1,22))), "chrX", "chrY")
limitPloidy=16

method = "error" # model or error
# method = "model"
globalModel = NULL
# MAKE SURE TO CHOOSE THE APPROPRIATE BIN SIZE (see publication for details)
# globalModel=paste0("~/scAbsolute/data/models/model.10X.", binSize, "kb.rds")
stopifnot(is.null(globalModel) || 
    (file.exists(globalModel) || file.exists(paste0("~/scAbsolute/data/models/", globalModel)) || 
     file.exists(normalizePath(paste0(.libPaths(), "/scAbsolute/data/models/", globalModel))) 
    && endsWith(globalModel, paste0(binSize, "kb\\.rds"))))
## END PARAMETERS

# load correct data
filePaths = gsub("//", "/", c(paste(bamPath, bamFile, sep="/")))
stopifnot(all(file.exists(filePaths)) || dir.exists(bamPath))
if(all(file.exists(filePaths))){
} else {
  if(any(bamFile == "" | length(bamFile) == 0)){
    filePaths = file.path(bamPath, list.files(path=bamPath, pattern = "\\.bam$"))
  }
}
stopifnot(all(file.exists(filePaths)))

# readCounts = readData(filePaths, binSize=binSize,
#                       species=species, genome=genome,
#                       filterChromosomes=c("MT"))
# 
# countsObject = readCounts


## run the core algorithm - using scAbsolute wrapper function
scaledCN = scAbsolute(filePaths, method=method, globalModel=globalModel,binSize=binSize, genome=genome,
           minLength=minLength, batchSize=1, testStatistic=testStatistic, limitPloidy=limitPloidy,
           minPloidy=minPloidy, maxPloidy=maxPloidy, ploidyWindow=ploidyWindow, ploidyRegion=ploidyRegion, selectRegion=selectRegion,
           splitPerChromosome=splitPerChromosome, optimizeSegmentation=optimizeSegmentation,
           penalty=penalty,
           randomSeed=randomSeed, outputPath=RESULTPATH, outputSegmentation=TRUE, debug=TRUE)


# single_cell_segmentation <- function(scaledCN, psi, 
                                      # max_state=12, truncate_cn=30, 
                                      # gc_correction=TRUE, debug=FALSE,
# n_samples=10000, randomSeed=2021)

psi = 0.05
max_state=10
truncate_cn=30
debug=FALSE

n_samples = 10000
randomSeed = 2021
# TODO clean up dependencies
source("~/scUnique/R/gcCorrection.R")
gc_correction = estimate_gc_correction(scaledCN)
# TODO move single cell segmenation to scAbsolute
source("~/scUnique/R/segmentation.R")


n_cells = dim(scaledCN)[[2]]
# TODO CHECK that element are available and within reasonable range
rpc = scaledCN@phenoData@data$rpc
# rpc = object@phenoData@data$rpc
alpha = scaledCN@phenoData@data$alpha_nb
# scaledCN@phenoData@data$omega
omega = rep(1, n_cells)
scaledCN@phenoData@data$zero_offset = 0.5
zero_offset = rep(0.5, n_cells)


## true changepoints
groundtruth <- readr::read_tsv("~/UID-SIM-ploidy-Bulk-300_SLX-00000_000089.groundtruth.bed", 
                               col_names = FALSE, col_types=list(col_character(), col_double(), 
                                                                 col_double(), col_double()))
library(BSgenome.Hsapiens.UCSC.hg19)
source("~/Iapetus/validation/00_validation-library.R")
binned.genome = bin_genome(binSize * (1+0))
valid = binsToUseInternal(scaledCN)
binned.genome[!valid] = NULL
SV.bed.binary <- bed_to_binary(groundtruth, binned.genome,
                               as.character(seq(1, 22)), windowSize = 1 * binSize * 1000)
cpt_true = which(SV.bed.binary == 1)
cpt_true = cpt_true[c(-100, diff(cpt_true)) == 1]




# function(cellindex){
#   
# }
# TODO remove
  cellindex = 2
  object = scaledCN[, cellindex]
  cell_gc_correction = gc_correction[, cellindex,drop=FALSE]

  logprob = compute_logprob(object, responsibilities=NULL,
                            correction = cell_gc_correction, 
                            max_state=max_state, truncate_cn=truncate_cn)

  stopifnot(dim(logprob)[[3]] == 1)
  lp = logprob[,,1]

  source("~/scAbsolute/data/changepoint/wrap_PELT.R")
  # TODO constrain to AIC to MBIC range
  # if(debug){
    cpts = wrap_state_space_segmentation(lp, penalty = "MBIC", method="CROPS", debug=TRUE)
  # }else{
  #   { sink("/dev/null"); cpts = wrap_state_space_segmentation(lp, penalty = "MBIC", method="CROPS", debug=FALSE); sink(); }
  # }

  ## segmentation based on psi threshold on type I error
  n_cpts = length(cpts)

  h_list = list()
  time_list = list()
  h_list_hits = list()
  # for each pair of neighboring segmentations, compute pvalues of new segments
  for(cutoff in seq(1,n_cpts-1)){
    start_time <- Sys.time()


    if(cutoff %% 10 == 1){
      print(cutoff)
    }
    
    # TODO remove
    # cutoff = 1
    # lp, h0_cpt, h1_cpt, rpc, h0, h1, calls
    h1_cpt = cpts[[cutoff]]
    h0_cpt = cpts[[cutoff+1]]
    # delta_cpt = setdiff(h1_cpt, h0_cpt)
    delta_cpt = h1_cpt
    hits = 0
    for(i1 in delta_cpt){
      if(any(abs(i1 - cpt_true) < 3))
        hits = hits + 1
    }
    
    # hits = 
    # delta_cpt %-% cpt_true
  
    # NOTE, we use raw value here
    reads = Biobase::assayDataElement(object, "calls")
    valid = binsToUseInternal(object)
    good_reads = reads[valid,1]
  
    cpt_to_seg = function(ct, lp){
      if(ct[length(ct)] != dim(lp)[[1]]){
        ct = c(ct, dim(lp)[[1]])
      }
      segments = rep(c(1:length(ct)), times=c(ct[1], diff(ct)))
      reg = S4Vectors::Rle(segments)
      ir = IRanges::ranges(reg)
      cnvalues = smooth_operator(1:(dim(lp)[[1]]), ir, function(xrt){return(which.max(apply(lp[xrt, ,drop=FALSE], 2, sum)) - 1)}, minLength = 0, trimLength=0)
      cnseq = rep(cnvalues, times=reg@lengths)
      return(cnseq)
    }
  
    h1_seq = cpt_to_seg(h1_cpt, lp)
    h0_seq = cpt_to_seg(h0_cpt, lp)
    stopifnot(length(h1_seq) == length(h0_seq))
    stopifnot(length(h1_seq) == length(good_reads))
  
    candidates = Rle(paste0(h1_seq != h0_seq, "-", h0_seq, "-", h1_seq))
    counter = 0

    rpc_i = rpc[cellindex]
    alpha_i = alpha[cellindex]
    omega_i = omega[cellindex]
    
    plist = c()
    invalid = c()
    counter = 0
    for(i in 1:length(candidates@values)){
      start = counter + 1
      end = counter + candidates@lengths[i]
      
      # if(i %% 10 == 1){
      #   print(i)
      # }
      # print(paste0(start, " - ", end))
      # start = 7
      # end = 8
      # i = 3
      
      if(base::startsWith(candidates@values[i], "TRUE")){
        # get alpha for this test
        observed_reads = as.vector(good_reads[start:end])
        valid_reads = !is.na(observed_reads)
        h0 = h0_seq[start:end][valid_reads]
        # h0 = table(h0)
        temp <- table(as.vector(h0))
        h0 = as.integer(names(temp)[temp == max(temp)])
        h1 = h1_seq[start:end][valid_reads]
        # print("----")
        # print(i)
        # print(h0)
        stopifnot(length(unique(h1)) == 1)
        h1 = h1[1]
        
        # TODO replace with gc correction
        mu_offset = as.vector(as.matrix(rep(0.0, end-start-sum(!valid_reads)+1), nrow=1))
        stopifnot(length(mu_offset) == length(observed_reads[valid_reads]))
        
        julia <- julia_setup()
        julia_source("/home/schnei01/scUnique/R/likelihoodRatioTestNegBin.jl")
        pvalue = julia_do.call(
          "likelihood_ratio_test_interface",
          list(as.integer(observed_reads[valid_reads]), 
               as.integer(h0), as.integer(h1), rpc_i, alpha_i, omega_i, 
               mu_offset, as.integer(n_samples), as.integer(randomSeed)),
          need_return = "R",
          show_value = FALSE
        )
        
        plist = c(plist, pvalue)
        invalid = c(invalid, i)
      }
      
      counter = counter + candidates@lengths[i]
    }
    
    
    values = unlist(lapply(strsplit(candidates@values, "-"), function(x) return(as.logical(x[[1]]))))
    h_list[[cutoff]] = plist
    h_list_hits[[cutoff]] = hits
    end_time <- Sys.time()
    time_list[[cutoff]] = difftime(end_time,start_time,units="mins")
    # ids = base::startsWith(candidates@values, "TRUE")
    # consider adjusting plist for multiple testing
    
  }

require(purrr)
df_p <- map_df(h_list, ~as.data.frame(.x), .id="cpt")
df_p$index = 1:dim(df_p)[[1]]
df_p$pvalue = df_p$.x
df_p$cpt = as.numeric(df_p$cpt)
df_p$pvalue.adj = p.adjust(df_p$pvalue, method="BH") #+ rnorm(mean=0.0, sd=0.001, n=length(df_p$pvalue))
# hits = data.frame(index=1:(n_cpts-1), hits=unlist(lapply(h_list_hits, sum)))
# dplyr::left_join(df_p, hits, by=)

fit = stats::loess(pvalue.adj ~ as.numeric(cpt), data=df_p, span=0.3)
preds = predict(fit, newdata=data.frame(cpt=df_p$cpt))#seq(1, n_cpts)))
which(pred < 0.01)
# max(1, min())
df_p$preds = preds

library(cowplot)
# %>% dplyr::filter(pvalue.adj < 0.05)
ggplot(data=df_p) + geom_point(aes(x=cpt, y=pvalue.adj)) + 
  geom_point(aes(x=cpt, y=preds), color="cyan") + 
  stat_smooth(aes(x=cpt, y=pvalue.adj), formula=y ~ x, span=0.3, method="loess", color="orange") + 
  geom_hline(yintercept = 0.01, color="red") #+ #coord_cartesian(ylim=c(0, 0.1))
  # coord_trans(y='log2') + theme_cowplot()

  

# values[invalid[p.adjust(plist, method="BY") >= 0.05]] = FALSE
cnv_profile = rep(NA, sum(valid))
choice = rep(values, candidates@lengths)
cnv_profile[!choice] = h0_seq[!choice]
cnv_profile[choice] = h1_seq[choice]

cnv_profile
Rle(cnv_profile)

qqnorm(p.adjust(unlist(h_list), method="BY"))

# quantile()
plot(log(unlist(lapply(h_list, mean, na.rm=TRUE))))





# 
# lrtests = t(sapply(seq(1, n_segments), function(index){
#   observed_reads = as.vector(reads[start_indexes[index]:end_indexes[index]])
#   valid_reads = !is.na(observed_reads)
#   # print(observed_reads[valid_reads])
#   # print(profile_background[!valid_reads])
#   # print(profile_background[start_indexes[index]:end_indexes[index]][valid_reads])
#   stopifnot(all(is.na(profile_background[start_indexes[index]:end_indexes[index]][!valid_reads])))
#   h0 = profile_background[start_indexes[index]:end_indexes[index]][valid_reads]
#   stopifnot(length(unique(h0)) == 1)
#   h0 = h0[1]
#   h1 = as.vector(copynumber[index])
#   h1_alt = profile[start_indexes[index]:end_indexes[index]][valid_reads]
#   stopifnot(length(unique(h1_alt)) == 1)
#   stopifnot(all(h1_alt == h1))
#   
#   if(h0 == h1){
#     return(c(NA, h0, h1, sum(valid_reads), mean(observed_reads[valid_reads]), NA))
#   }
#   
#   mu_offset = as.vector(cell_gc_correction[start_indexes[index]:end_indexes[index]])
#   
#   # TODO update from package installation
#   # TODO only load once
#   julia <- julia_setup()
#   julia_source("/home/schnei01/scUnique/R/likelihoodRatioTestNegBin.jl")
#   pvalue = julia_do.call(
#     "likelihood_ratio_test_interface",
#     list(as.integer(observed_reads[valid_reads]), 
#          as.integer(h0), as.integer(h1), rpc, alpha, omega, 
#          mu_offset[valid_reads], as.integer(n_samples), as.integer(randomSeed)),
#     need_return = "R",
#     show_value = FALSE
#   )
#   
#   return(c(pvalue, 
#            h0, h1, sum(valid_reads), 
#            mean(observed_reads[valid_reads]), mean(mu_offset[valid_reads])))
# }))
# 
# # output 



cpt = cpts[[20]]
segments = rep(c(1:length(cpt)), times=c(cpt[1], diff(cpt)))
reg = S4Vectors::Rle(segments)
ir = IRanges::ranges(reg)
cnvalues = smooth_operator(1:(dim(lp)[[1]]), ir, function(xrt){return(which.max(apply(lp[xrt, ,drop=FALSE], 2, sum)) - 1)}, minLength = 0, trimLength=0)
cnseq = rep(cnvalues, times=reg@lengths)
seg_space = assayDataElement(object, "segmented")
valid = binsToUseInternal(object)
seg_space[valid] = cnseq

object = assayDataElementReplace(object, "segmented", matrix(seg_space, ncol=1))
plotCopynumber(selectChromosomes(object, include=c("5", "6")))
plotCopynumber(object)


output = statespace_segmentation(logprob[,,1], pen.value=0, penalty = "NONE", method="CROPS")

Rle(output)

output = statespace_segmentation(logprob[,,1], pen.value=0, penalty = "AIC")
Rle(output)


# plotCopynumber(scaledCN)
# plotCounts(readCounts[,1])
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
