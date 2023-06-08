# script can be run independently, create artefact at the end of script (make sure to uncomment saveRDS command)
library(QDNAseq)
library(GenomicRanges)
require(S4Vectors)
rm(list=ls())

# reptime_path <- "~/scAbsolute/data/replicationTiming"
reptime_path <- "/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/data/replicationTiming" 

#Human
set.seed(2020)
exampleBamFile = "~/Data/replication-timing/UID-10X-Fibroblast-cell_SLX-00000_000001_AAACCTGAGCAGCCTC-1.bam"

options(future.globals.maxSize= 4096*1024^2)
source("~/scAbsolute/R/core.R")

for (binSize in c(10, 15, 30, 50, 100, 500, 1000)) { 
  #bins smaller than 10kb not supported
  
  print(paste0('Bin size: ', binSize))
  
  test = readData(exampleBamFile, binSize=binSize)
  # GRanges from readCounts object
  info = rownames(test)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
  chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
  chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
  stopifnot(length(seqName) == length(chunkStart))
  
  gr = GRanges(seqName, IRanges(chunkStart, chunkEnd), mcols=data.frame(gc_baseline=test@featureData@data$gc))
  
  saveRDS(gr, file=paste0(reptime_path, "/gr_data/GenomicRanges-GRCh37-", binSize, "kb.RDS"))

}


#Mouse
set.seed(2020)
exampleBamFile = "/home/adr44/rds/hpc-work/scAbsolute/github/scDNAseq-workflow/mouse_compatibility/data/align/SH171012_I_028.bam"

options(future.globals.maxSize= 4096*1024^2)
source("~/scAbsolute/R/core.R")
# source("/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/R/core.R")


for (binSize in c(10, 15, 30, 50, 100, 500, 1000)) {
  #bins smaller than 10kb not supported
  
  print(paste0('Bin size: ', binSize))

  test = readData(exampleBamFile, binSize=binSize, species="Mouse", genome="GRCm38")
  # GRanges from readCounts object
  info = rownames(test)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
  chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
  chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
  stopifnot(length(seqName) == length(chunkStart))
  
  gr = GRanges(seqName, IRanges(chunkStart, chunkEnd), mcols=data.frame(gc_baseline=test@featureData@data$gc))
  
  saveRDS(gr, file=paste0(reptime_path, "/gr_data/GenomicRanges-GRCm38-", binSize, "kb.RDS"))
}



