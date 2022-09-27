# script can be run independently, create artefact at the end of script (make sure to uncomment saveRDS command)
library(QDNAseq)
library(GenomicRanges)
require(S4Vectors)
rm(list=ls())
set.seed(2020)
exampleBamFile = "~/Data/replication-timing/UID-10X-Fibroblast-cell_SLX-00000_000001_AAACCTGAGCAGCCTC-1.bam"
binSize = 10

options(future.globals.maxSize= 4096*1024^2)
source("~/scAbsolute/R/core.R")
test = readData(bamFile, binSize=binSize)
# GRanges from readCounts object
info = rownames(test)
seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
stopifnot(length(seqName) == length(chunkStart))

gr = GRanges(seqName, IRanges(chunkStart, chunkEnd), mcols=data.frame(gc_baseline=test@featureData@data$gc))

saveRDS(gr, file=paste0("~/scAbsolute/data/replicationTiming/GenomicRanges-GRCh37-", binSize, "kb.RDS"))
