## Replication timing -> cell cycle prediction
# script can be run independently, create artefact at the end of script (make sure to uncomment saveRDS command)
rm(list=ls())
set.seed(2020)
require(rtracklayer)
require(GenomicRanges)
require(S4Vectors)
require(QDNAseq)

# source("~/scAbsolute/R/core.R")
source("/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/R/core.R")

# reptime_path <- "~/scAbsolute/data/replicationTiming"
reptime_path <- "/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/data/replicationTiming" 

##### Human replication information

## load replication timing tracks ====
# files = sub(pattern = "(.*)\\..*$", replacement = "\\1", list.files(path, pattern = "\\.bigWig$"))
# file source: http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeFsuRepliChip  

path = "~/Data/replication-timing/"
files = c("wgEncodeFsuRepliChipBg02esWaveSignalRep1",
          "wgEncodeFsuRepliChipBg02esWaveSignalRep2",
          "wgEncodeFsuRepliChipGm06990WaveSignalRep1",
          "wgEncodeFsuRepliChipGm06990WaveSignalRep2",
          "wgEncodeFsuRepliChipH1hescWaveSignalRep1",
          "wgEncodeFsuRepliChipH1hescWaveSignalRep2",
          "wgEncodeFsuRepliChipH1hescWaveSignalRep3",
          "wgEncodeFsuRepliChipH7esWaveSignalRep1",
          "wgEncodeFsuRepliChipH7esWaveSignalRep2",
          "wgEncodeFsuRepliChipH9esWaveSignalRep1",
          "wgEncodeFsuRepliChipHelas3WaveSignalRep1",
          "wgEncodeFsuRepliChipImr90WaveSignalRep1",
          "wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep1",
          "wgEncodeFsuRepliChipIpshfib2ips4WaveSignalRep2",
          "wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep1",
          "wgEncodeFsuRepliChipIpshfib2ips5WaveSignalRep2")
ordering =  c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
              "chr22", "chrX", "chrY")

replicationInfo = list()
for(f in files){
  print(paste0("File ", f))

  # load signal track
  signal <- rtracklayer::import(paste0(path, f, ".bigWig"))

  ## create bin level replication timing information
  signal_rle = mcolAsRleList(signal, "score")

  signal_rle = signal_rle[order(match(names(signal_rle), ordering))]
  signal_rle = signal_rle[names(signal_rle) %in% ordering]

  replicationInfo[[f]]=signal_rle
}

## Average information across data sets ====
files = names(replicationInfo)
avgSignal = list()

for(chr in 1:24){
  print(paste0("Chromosome ", chr))
  chr_length = sum(as.numeric(replicationInfo[[1]][[chr]]@lengths))
  vector = rep(NA, chr_length)
  support = rep(0, chr_length)
  for(f in files){
    tmp = S4Vectors::decode(replicationInfo[[f]][[chr]])
    idx = is.na(vector) & !is.na(tmp)
    vector = ifelse(idx, 0, vector) + tmp
    support = support + ifelse(!is.na(tmp), 1.0, 0)
  }
  vector = vector / as.numeric(ifelse(support == 0, NA, support))
  s = S4Vectors::Rle(vector)
  avgSignal[[gsub("chr", "", ordering[chr])]] = s
}

# convert to RleList object
replicationInfo = avgSignal
rep_gr = do.call("c", lapply(names(replicationInfo), function(xy){
  xz = replicationInfo[[xy]]
  grt = GRanges(xy,IRanges(cumsum(c(0,runLength(xz)[-nrun(xz)])),
                           width=runLength(xz)),
                values = runValue(xz))
  seqlevels(grt) <- names(replicationInfo)
  return(grt)
}))
rep_rle = GenomicRanges::mcolAsRleList(rep_gr[!is.na(rep_gr$values)], "values")

## save RleList to file for import into main programs
saveRDS(rep_rle, file=paste0(reptime_path, "/replicationInfo.RDS"))

## to be used with
# replicationTiming = binnedAverage(gr, rep_rle, "replicationTime", na.rm=TRUE)
## see replicationInfo_per_bin.R


##### Mouse replication information

# file source: http://hgdownload.soe.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeFsuRepliChip/ or http://genome.ucsc.edu/cgi-bin/hgFileUi?db=mm9&g=wgEncodeFsuRepliChip   ###mouse


path = "/home/adr44/rds/hpc-work/scAbsolute/github/scAbsolute/data/replicationTiming/repliChip_data/" #"~/Data/replication-timing/"
files = c("wgEncodeFsuRepliChipCh12FWaveSignalRep1",
          "wgEncodeFsuRepliChipCh12FWaveSignalRep2",
          "wgEncodeFsuRepliChipEpisc5MWaveSignalRep1",
          "wgEncodeFsuRepliChipEpisc5MWaveSignalRep2",
          "wgEncodeFsuRepliChipEpisc7FWaveSignalRep1",
          "wgEncodeFsuRepliChipEpisc7FWaveSignalRep2",
          "wgEncodeFsuRepliChipEs46cMDifff6dWaveSignalRep1",
          "wgEncodeFsuRepliChipEs46cMWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MDiffe3dWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MDiffe3dWaveSignalRep2",
          "wgEncodeFsuRepliChipEsd3MDiffe6dWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MDiffe6dWaveSignalRep2",
          "wgEncodeFsuRepliChipEsd3MDiffe9dWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MDiffe9dWaveSignalRep2",
          "wgEncodeFsuRepliChipEsd3MDiffg3dWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MDiffg3dWaveSignalRep2",
          "wgEncodeFsuRepliChipEsd3MWaveSignalRep1",
          "wgEncodeFsuRepliChipEsd3MWaveSignalRep2",
          "wgEncodeFsuRepliChipEsem5sUDiffhsoxmWaveSignalRep1",
          "wgEncodeFsuRepliChipEsem5sUDiffhsoxmWaveSignalRep2",
          "wgEncodeFsuRepliChipEsem5sUDiffhsoxpWaveSignalRep1",
          "wgEncodeFsuRepliChipEsem5sUDiffhsoxpWaveSignalRep2",
          "wgEncodeFsuRepliChipEstt2MDifff9dWaveSignalRep1",
          "wgEncodeFsuRepliChipEstt2MDifff9dWaveSignalRep2",
          "wgEncodeFsuRepliChipEstt2MWaveSignalRep1",
          "wgEncodeFsuRepliChipEstt2MWaveSignalRep2",
          "wgEncodeFsuRepliChipJ185aUWaveSignalRep1",
          "wgEncodeFsuRepliChipJ185aUWaveSignalRep2",
          "wgEncodeFsuRepliChipL1210FWaveSignalRep1",
          "wgEncodeFsuRepliChipL1210FWaveSignalRep2",
          "wgEncodeFsuRepliChipMefMWaveSignalRep1",
          "wgEncodeFsuRepliChipMelMWaveSignalRep1",
          "wgEncodeFsuRepliChipMelMWaveSignalRep2")
ordering =  c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
              "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")

replicationInfo = list()
for(f in files){
  print(paste0("File ", f))
  
  # load signal track
  signal <- rtracklayer::import(paste0(path, f, ".bigWig"))
  
  ## create bin level replication timing information
  signal_rle = mcolAsRleList(signal, "score")
  
  signal_rle = signal_rle[order(match(names(signal_rle), ordering))]
  signal_rle = signal_rle[names(signal_rle) %in% ordering]
  
  replicationInfo[[f]]=signal_rle
}

## Average information across data sets ====
files = names(replicationInfo)
avgSignal = list()

for(chr in 1:21){
  print(paste0("Chromosome ", chr))
  chr_length = sum(as.numeric(replicationInfo[[1]][[chr]]@lengths))
  vector = rep(NA, chr_length)
  support = rep(0, chr_length)
  for(f in files){
    tmp = S4Vectors::decode(replicationInfo[[f]][[chr]])
    idx = is.na(vector) & !is.na(tmp)
    vector = ifelse(idx, 0, vector) + tmp
    support = support + ifelse(!is.na(tmp), 1.0, 0)
  }
  vector = vector / as.numeric(ifelse(support == 0, NA, support))
  s = S4Vectors::Rle(vector)
  avgSignal[[gsub("chr", "", ordering[chr])]] = s
}

# convert to RleList object
replicationInfo = avgSignal
rep_gr = do.call("c", lapply(names(replicationInfo), function(xy){
  xz = replicationInfo[[xy]]
  grt = GRanges(xy,IRanges(cumsum(c(0,runLength(xz)[-nrun(xz)])),
                           width=runLength(xz)),
                values = runValue(xz))
  seqlevels(grt) <- names(replicationInfo)
  return(grt)
}))
rep_rle = GenomicRanges::mcolAsRleList(rep_gr[!is.na(rep_gr$values)], "values")

## save RleList to file for import into main programs
saveRDS(rep_rle, file=paste0(reptime_path, "/replicationInfo_mouse.RDS")) 




