## Replication timing -> cell cycle prediction
# script can be run independently, create artefact at the end of script (make sure to uncomment saveRDS command)
rm(list=ls())
set.seed(2020)
require(rtracklayer)
require(GenomicRanges)
require(S4Vectors)
require(QDNAseq)
source("~/scAbsolute/R/core.R")

## load replication timing tracks ====
path = "~/Data/replication-timing/"
# files = sub(pattern = "(.*)\\..*$", replacement = "\\1", list.files(path, pattern = "\\.bigWig$"))
# file source: http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeFsuRepliChip
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
saveRDS(rep_rle, file="~/scAbsolute/data/replicationTiming/replicationInfo.RDS")

## to be used with
# replicationTiming = binnedAverage(gr, rep_rle, "replicationTime", na.rm=TRUE)
