## scAbsolute - parallel run script

## Variables ####
args = commandArgs(trailingOnly=TRUE)

if (interactive()){
  rm(list=ls()) 
  
  bamPath = "~/Data/project1/20-align/"
  # bamFile = c("UID-10X-COLO829_SLX-00000_000326_ATGTGTGCAAAGTAAC-1.bam", "UID-10X-Fibroblast-cell_SLX-00000_000277_CACACCTTCACACGTA-1.bam")
  # bamFile = "UID-10X-COLO829_SLX-00000_000326_ATGTGTGCAAAGTAAC-1.bam"
  bamFile = "UID-JBL-CIOV1_SLX-18430_000022_SINCEL-31-Plate-95-C2-CGTACTAG-CCTAGAGT.bam"
  # bamFile = "UID-10X-Fibroblast-cell_SLX-00000_000914_TCATTTGTCTGGTGGC-1.bam"
  # bamFile = c("UID-PEO23-2D_SLX-18072_000127_GGTATTGAGTCGCCGT-1.bam",
  #             "UID-PEO23-2D_SLX-18072_000085_CTCTAATCATCCTGGG-1.bam")
  # c("UID-PEO23-2D_SLX-18072_000117_GCTTCCATCATATCTC-1.bam",
  #   "UID-PEO23-2D_SLX-18072_000106_GAGCAGATCGTTACCC-1.bam",
  #   "UID-PEO14-2D_SLX-18073_000280_CTCGAGGGTTCCCTAC-1.bam")
  
  RESULTPATH = "~/UncertainDecember.rds"
  
  # options
  binSize = 30
  
  minPloidy = 1.5
  maxPloidy = 2.5
  
  reticulate::use_condaenv(condaenv = "rstudio-server", conda = "/home/schnei01/.anaconda3/bin/conda")
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
  
}

# TODO BIOCONDUCTOR
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

species = "Human"
genome = "hg19"

trimLength = 1
minLength = 10
# optimization of segmentation
testStatistic = "NegBinMM"
splitPerChromosome=FALSE
optimizeSegmentation=FALSE

max_iterations=201
change_prob=1e-6
max_states=10

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
if(is.null(maxPloidy)) maxPloidy = 12.0

ploidyWindow = 0.1
ploidyRegion=c(paste0("chr", as.character(seq(1,22))), "chrX")
selectRegion=c(paste0("chr", as.character(seq(1,22))), "chrX", "chrY")
limitPloidy=16

method = "error" # model or error
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

readCounts = readData(filePaths, binSize=binSize, extendedBlacklisting=TRUE,
                      filterChromosomes = c("MT"), genome="hg19")






## Simulate GC impact on reads
valid = binsToUseInternal(readCounts)
all(readCounts@assayData$calls[valid,1] == readCounts@assayData$copynumber[valid,1])
all(readCounts@assayData$calls[valid,1] == readCounts@assayData$probdloss[valid,1])

gc = readCounts@featureData@data$gc[valid]
map = readCounts@featureData@data$map[valid]

states = rep(sample(c(1:9, 0.001), size=2000, replace=TRUE), each=500)[1:sum(valid)]
rpc_true = 7
covariate_base = scale(gc^2 + 0.05 * (map^2), scale=TRUE)

# covariate = (rnorm(covariate_base, mean=covariate_base, sd=0.2) + 1) / 27 * 7 * 100
covariate = exp(covariate_base / max(covariate_base)) #(exp(covariate_base))
summary(covariate)

mu_vec = states * rpc_true * covariate
ratio = covariate / (states * rpc_true)
table(mu_vec <= 0)
mu_vec[mu_vec <= 0] = 1
counts = rnbinom(mu_vec, size=100, mu=mu_vec)
newCalls = Biobase::assayDataElement(readCounts, "calls")
newSegments = Biobase::assayDataElement(readCounts, "calls")
newCalls[valid] = counts
newSegments[valid] = states * rpc_true
newReadCounts = readCounts
newReadCounts = Biobase::assayDataElementReplace(newReadCounts, "calls", newCalls)
newReadCounts = Biobase::assayDataElementReplace(newReadCounts, "copynumber", newCalls)
newReadCounts = Biobase::assayDataElementReplace(newReadCounts, "probdloss", newCalls)

## run the core algorithm - using scAbsolute wrapper function
scaledCN = scAbsolute(newReadCounts, method=method, globalModel=globalModel,binSize=binSize, genome=genome,
                      minLength=minLength, batchSize=1, testStatistic=testStatistic, limitPloidy=limitPloidy,
                      minPloidy=minPloidy, maxPloidy=maxPloidy, ploidyWindow=ploidyWindow, ploidyRegion=ploidyRegion, selectRegion=selectRegion,
                      splitPerChromosome=splitPerChromosome, optimizeSegmentation=optimizeSegmentation,
                      max_iterations=max_iterations, hmm_path=hmm_path, change_prob=change_prob, max_states=max_states,
                      randomSeed=randomSeed, outputPath=RESULTPATH, outputSegmentation=TRUE)

p1 = plotCopynumber(scaledCN, copynumber = FALSE, plotSegmentation = TRUE, correction = TRUE)
p2 = plotCopynumber(scaledCN, copynumber = FALSE, plotSegmentation = TRUE, correction = FALSE)

cowplot::plot_grid(p1, p2, nrow=2)
plotCounts(scaledCN[, 1],copynumber=FALSE,ylim=c(0,100))

debug = data.frame(gc=gc, raw=as.numeric(scaledCN[valid,1]@assayData$calls), 
                   ccor=as.numeric(scaledCN[valid,1]@assayData$probdloss))

q1 = ggplot(debug, aes(x=gc, y=raw, v=ccor)) + geom_point() + 
  stat_smooth(span=0.1, color="red", formula=y ~ x) + 
  stat_smooth(span=0.1, color="blue", formula=v ~ x) +
  theme_cowplot() + 
  geom_abline(slope=0.2, intercept=8, color="purple") + 
  geom_hline(yintercept=17, color="cyan")
q1




initialSegmentedCounts = segment(newReadCounts, penalty = "MBIC", pen.value = NULL, splitPerChromosome=TRUE)
# plotCopynumber(initialSegmentedCounts, ylim=c(0, 100), copynumber=FALSE,correction = FALSE, plotSegmentation = TRUE)

initialSegmentedCountsCorrected = initial_gc_correction(initialSegmentedCounts)
segmentedCounts = segment(initialSegmentedCountsCorrected, penalty = "MBIC", pen.value = NULL, splitPerChromosome=FALSE)
# plotCopynumber(segmentedCounts, ylim=c(0, 100), copynumber=FALSE,correction = TRUE, plotSegmentation = TRUE)

scaledCN = scAbsolute(newReadCounts, method=method, globalModel=globalModel,binSize=binSize, genome=genome,
                      minLength=minLength, batchSize=1, testStatistic=testStatistic, limitPloidy=limitPloidy,
                      minPloidy=1.2, maxPloidy=10.0, ploidyWindow=ploidyWindow, ploidyRegion=ploidyRegion, selectRegion=selectRegion,
                      splitPerChromosome=splitPerChromosome, optimizeSegmentation=optimizeSegmentation,
                      max_iterations=max_iterations, hmm_path=hmm_path, change_prob=change_prob, max_states=max_states,
                      randomSeed=randomSeed, outputPath=RESULTPATH, outputSegmentation=TRUE)
plotCopynumber(scaledCN, ylim=c(0, 50))

correctSegmentedCounts = Biobase::assayDataElementReplace(initialSegmentedCounts, "segmented", newSegments)
# plotCopynumber(correctSegmentedCounts, ylim=c(0, 100), copynumber=FALSE,correction = TRUE, plotSegmentation = TRUE)

gc_cor1 = estimate_gc_correction(initialSegmentedCounts)
gc_cor2 = estimate_gc_correction(segmentedCounts)
gc_cor3 = estimate_gc_correction(correctSegmentedCounts)
newScaledCN = Biobase::assayDataElementReplace(scaledCN, "segmented", scaledCN@assayData$copynumber)
gc_cor4 = estimate_gc_correction(newScaledCN)

# plot(gc_cor[valid,1], covariate)
cor(gc_cor1[valid,1], covariate)
cor(gc_cor2[valid,1], covariate)
cor(gc_cor3[valid,1], covariate)
cor(gc_cor4[valid,1], covariate)

table(round(states, digits=0) == round(scaledCN[valid,1]@assayData$segmented, digits=0))
table(round(states, digits=0) == round(scaledCN[valid,1]@assayData$copynumber, digits=0))

# all(correctSegmentedCounts[valid,1]@assayData$calls == correctSegmentedCounts[valid,1]@assayData$probdloss)

# ccor = as.numeric(correctSegmentedCounts[valid,1]@assayData$calls) * gc_cor[valid]
debug = data.frame(gc=gc, raw=as.numeric(correctSegmentedCounts[valid,1]@assayData$calls), 
                   ccor1=as.numeric(correctSegmentedCounts[valid,1]@assayData$calls) * (1/gc_cor1[valid]),
                   ccor2=as.numeric(correctSegmentedCounts[valid,1]@assayData$calls) * (1/gc_cor2[valid]),
                   ccor3=as.numeric(correctSegmentedCounts[valid,1]@assayData$calls) * (1/gc_cor3[valid]),
                   ccor=as.numeric(correctSegmentedCounts[valid,1]@assayData$probdloss))

q1 = ggplot(debug, aes(x=gc, y=raw, z1=ccor1, z2=ccor2, z3=ccor3, v=ccor)) + geom_point() + 
  stat_smooth(span=0.1, color="red", formula=y ~ x) + 
  stat_smooth(span=0.1, color="orange", formula=z1 ~ x) +
  stat_smooth(span=0.1, color="blue", formula=z2 ~ x) +
  stat_smooth(span=0.1, color="green", formula=z3 ~ x) +
  theme_cowplot() + 
  geom_abline(slope=0.2, intercept=8, color="purple") + 
  geom_hline(yintercept=17, color="cyan")
q1

