#Run post replicationInfo.R and createGR.R

rm(list=ls())
set.seed(2020)
require(rtracklayer)
require(GenomicRanges)
require(S4Vectors)
require(QDNAseq)
source("~/scAbsolute/R/core.R")

# reptime_path <- "~/scAbsolute/data/replicationTiming"
reptime_path <- "/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/data/replicationTiming" 

##Human
rep_rle <- readRDS(paste0(reptime_path, "/replicationInfo.RDS"))
repTime_list <- list()
for (binSize in c(10, 15, 30, 50, 100, 500, 1000)) {
  gr <- readRDS(paste0(reptime_path, "/gr_data/GenomicRanges-GRCh37-", binSize, "kb.RDS"))
  
  replicationTiming = binnedAverage(gr, rep_rle, "replicationTime", na.rm=TRUE)
  
  repTime_list[[binSize]] <- replicationTiming
}

saveRDS(repTime_list, file=paste0(reptime_path, "/replicationTiming_per_binSize.RDS"))


##Mouse
rep_rle <- readRDS(paste0(reptime_path, "/replicationInfo_mouse.RDS"))
repTime_list <- list()
for (binSize in c(10, 15, 30, 50, 100, 500, 1000)) {
  gr <- readRDS(paste0(reptime_path, "/gr_data/GenomicRanges-GRCm38-", binSize, "kb.RDS"))
  
  replicationTiming = binnedAverage(gr, rep_rle, "replicationTime", na.rm=TRUE)
  
  repTime_list[[binSize]] <- replicationTiming
}

saveRDS(repTime_list, file=paste0(reptime_path, "/replicationTiming_per_binSize_mouse.RDS"))