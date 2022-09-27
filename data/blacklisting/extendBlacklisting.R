## Script to identify regions that map badly or are not covered by existing centromer maps in hg19
# NOTE uncomment end in order to write to file and update the blacklisted regions
# NOTE the current map has a resolution of 10kb, for lower resolution we would require deeper sequencing of normal cells

# test = readRDS("~/Data/project1/keep_for_publication/extendedBlacklisting/GRCh38/UID-10X-HG38-BJ-CELLS_10.rds")

rm(list=ls()) 

binsToUse <- function(object){
           if ("use" %in% colnames(fData(object))) {
              return(fData(object)$use)
            } else if ("filter" %in% colnames(fData(object))) {
              fData(object)$use <- fData(object)$filter
              fData(object)$filter <- NULL
              return(fData(object)$use)
            }
            rep(TRUE, times=nrow(object))
}

median_read_depth <- function(object){
  
  reads = Biobase::assayDataElement(object, "copynumber")
  print(dim(reads))
  
  # NOTE we add a pseudocount here
  reads = reads + 1
  
  locations = rownames(object)
  median_reads = apply(reads, 2, median, na.rm=TRUE)
  
  # proportion of reads observed vs median reads of sample (diploid sample -> median_reads corresponds to ok bins)
  preads = t(t(reads) / median_reads)
  rm(median_reads)
  rm(reads)
  
  # median over cells
  position_reads = apply(preads, 1, median, na.rm=FALSE)
  # summary(position_reads)
  return(position_reads)
}

# blacklisting is done on files with original blacklisting by QDNAseq (hg19 build)
set.seed(2020)
library(QDNAseq, quietly = TRUE)
library(Biobase, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
library(cowplot, quietly = TRUE)

# percentage of unmapped regions in a bin
cutoff_svregions = 0.25 #(200bp region maximum)
lower_bound = 0.95
upper_bound = 1.05
lower_bound_x = 0.9
upper_bound_x = 1.1
# bound for cross datasets (DLP, CCL) - all chromosomes, larger range
lower_bound_lq = 0.8
upper_bound_lq = 1.2
read_cutoff = 600000
read_cutoff_x = 800000
isGRCh37=TRUE
isCROSS=TRUE
stopifnot(!isGRCh37 | (isCROSS & isGRCh37))
if(isGRCh37){
  sampleFile1 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh37/UID-10X-HG19-BJ-CELLS_10.rds"
  sampleFile2 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh37/UID-10X-HG19-BJ-NUCLEI_10.rds"
  if(isCROSS){
    sampleFile3 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh37/UID-DLP-SA928_10.rds"
    CN3 = readRDS(sampleFile3)
    sampleFile4 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh37/UID-CCL-DIPLOID-SUBSET3_10.rds"
    CN4 = readRDS(sampleFile4)
  }
}else{
  sampleFile1 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh38/UID-10X-HG38-BJ-CELLS_10.rds"
  sampleFile2 = "~/Data/Studies/SOURCES/blacklistingRegions/GRCh38/UID-10X-HG38-BJ-NUCLEI_10.rds"
  # NOTE TO REMOVE THE EXTERNAL BLACKLISTING SOURCE WHEN USING GRCh38! 
}

CN1 = readRDS(sampleFile1)
CN2 = readRDS(sampleFile2)
object = CN1

# create gr object based on binsize of segmentedCounts
info = rownames(object)
seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
stopifnot(length(seqName) == length(chunkStart))
gr = GRanges(seqName, IRanges(chunkStart, chunkEnd))

# load blacklisted regions (GRCh37)
blacklistedRegions = "~/scAbsolute/data/blacklisting/ceph18.b37.lumpy.exclude.2014-01-15.bed"
excluded_gr = rtracklayer::import(blacklistedRegions)
chromosomes = c(as.character(seq(1:22)), "X", "Y")
excluded_gr = keepSeqlevels(excluded_gr[seqnames(excluded_gr) %in% chromosomes], chromosomes)
excluded_gr$width = excluded_gr@ranges@width

rle_list = mcolAsRleList(excluded_gr, "width")
for(l in 1:length(rle_list)){
  x = rle_list[[l]]
  x@values=ifelse(is.na(x@values), 0, 1)
  rle_list[[l]] = x
}

bavg = binnedAverage(gr, rle_list, "percentage_excluded")
is_blacklisted = bavg$percentage_excluded > cutoff_svregions
is_na = is.na(CN1@assayData$copynumber[,1])
is_notvalid = !binsToUse(CN1)
## Analysis needs to be done separately for chr1:22 and chrX (chrY we don't have additional information)

# Chromsomes 1:22 ----

## Initial blacklist based on SV information and external (QDNAseq, Speedseq) information ====
gr_chr = gr[seqnames(gr) %in% as.character(seq(1, 22))]
idx_chr = gr %over% gr_chr
stopifnot(length(idx_chr) == dim(object)[[1]])

table(is_na[idx_chr], is_blacklisted[idx_chr])
if(isGRCh37){
  blacklist_external = (is_na | is_blacklisted) & idx_chr
}else{
  # NOTE: in case of GRCh38
  blacklist_external = (is_na) & idx_chr
}

## Additional blacklist based on empirical coverage 
# Identify bins that have unusual relative proportion of reads mapped to them -> based on CN1 ####

position_reads = median_read_depth(CN1[idx_chr, CN1@phenoData@data$used.reads > read_cutoff])
# Determine cutoff for additional bins that will be blacklisted
blacklist2 = position_reads < lower_bound | position_reads > upper_bound
table(blacklist2)

# visualize distribution of median percentage 
p1 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
  geom_vline(xintercept=c(lower_bound, upper_bound), color="red") + theme_minimal()


## Check whether blacklist2 is robust with independent data set (CN2) ####
## Identify bins that have unusual relative proportion of reads mapped to them

position_reads = median_read_depth(CN2[idx_chr, CN2@phenoData@data$used.reads > read_cutoff])
# Determine cutoff for additional bins that will be blacklisted
blacklist3 = position_reads < lower_bound | position_reads > upper_bound
table(blacklist3)

# visualize distribution of median percentage 
p2 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
  geom_vline(xintercept=c(lower_bound, upper_bound), color="red") + theme_minimal()
table(blacklist2, blacklist3)

if(isCROSS){
  position_reads = median_read_depth(CN3[idx_chr, CN3@phenoData@data$used.reads > read_cutoff])
  # Determine cutoff for additional bins that will be blacklisted
  blacklist4 = position_reads < lower_bound_lq | position_reads > upper_bound_lq
  table(blacklist4)
  
  p3 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
    geom_vline(xintercept=c(lower_bound_lq, upper_bound_lq), color="red") + theme_minimal()
  
  
  position_reads = median_read_depth(CN4[idx_chr, CN4@phenoData@data$used.reads > read_cutoff])
  # Determine cutoff for additional bins that will be blacklisted
  blacklist5 = position_reads < lower_bound_lq | position_reads > upper_bound_lq
  table(blacklist5)
  p4 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
    geom_vline(xintercept=c(lower_bound_lq, upper_bound_lq), color="red") + theme_minimal()
  
}



# Chromsome X and Y ----
gr_chr = gr[seqnames(gr) %in% c("X", "Y")]
idx_chr = gr %over% gr_chr
stopifnot(length(idx_chr) == dim(object)[[1]])

table(is_blacklisted[idx_chr])
blacklist_external_X = is_blacklisted & idx_chr


## Additional blacklist based on empirical coverage
## Identify bins that have unusual relative proportion of reads mapped to them -> based on CN1 ####
# Only select cells with higher read depth to adjust for 50% less reads on X chromosome
position_reads = median_read_depth(CN1[idx_chr, CN1@phenoData@data$used.reads > read_cutoff_x])
# Determine cutoff for additional bins that will be blacklisted
blacklist2_X = position_reads < lower_bound_x | position_reads > upper_bound_x
table(blacklist2_X)

# visualize distribution of median percentage 
q1 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
  geom_vline(xintercept=c(lower_bound_x, upper_bound_x), color="red") + theme_minimal()

## Check whether blacklist2 is robust with independent data set (CN2) ####
## Identify bins that have unusual relative proportion of reads mapped to them

position_reads = median_read_depth(CN2[idx_chr, CN2@phenoData@data$used.reads > read_cutoff_x])
# Determine cutoff for additional bins that will be blacklisted
blacklist3_X = position_reads < lower_bound_x | position_reads > upper_bound_x
table(blacklist3_X)

# visualize distribution of median percentage 
q2 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
  geom_vline(xintercept=c(lower_bound_x, upper_bound_x), color="red") + theme_minimal()


if(isCROSS){
  position_reads = median_read_depth(CN3[idx_chr, CN3@phenoData@data$used.reads > read_cutoff_x])
  # Determine cutoff for additional bins that will be blacklisted
  blacklist4_X = position_reads < lower_bound_lq | position_reads > upper_bound_lq
  table(blacklist4_X)
  
  q3 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
    geom_vline(xintercept=c(lower_bound_lq, upper_bound_lq), color="red") + theme_minimal()
  
  
  position_reads = median_read_depth(CN4[idx_chr, CN4@phenoData@data$used.reads > read_cutoff_x])
  # Determine cutoff for additional bins that will be blacklisted
  blacklist5_X = position_reads < lower_bound_lq | position_reads > upper_bound_lq
  table(blacklist5_X)
  
  q4 = ggplot(data=data.frame(x=position_reads)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 2.3)) + 
    geom_vline(xintercept=c(lower_bound_lq, upper_bound_lq), color="red") + theme_minimal()
  
}



## Create joined blacklist


## Make sure gene regions are kept

# # TODO check whether genes fall in NA regions
# # TODO replicate for non
# # TODO finish this part
# cancer_genes = readr::read_tsv(file="~/scAbsolute/data/blacklisting/NCG6_cancergenes.tsv", col_names=TRUE)
# gene_annotation = read_tsv("~/scAbsolute/data/blacklisting/gencode.v19.tsv", col_names = c("chromosome", "annotation", "start", "end", "gene_id", "gene_name"))
# gene_annotation = gene_annotation[gene_annotation$gene_name %in% cancer_genes$symbol, ]
# 
# gene_gr = GRanges(seqnames=gsub("chr", "", gene_annotation$chromosome), ranges = IRanges(start=gene_annotation$start, end=gene_annotation$end), 
#                   gene_id=gene_annotation$gene_id, gene_name=gene_annotation$gene_name)
# 
# gr_genes = GenomicRanges::findOverlaps()
# # gr = getSegTable(scaledCN, minLength = minLength)[[1]]
# 
# # gr_nan = createGR(scaledCN, variable=is.na(scaledCN@assayData$calls))
# a = gr_nan$`annotation.UID-10X-Spikein-1pc_SLX-00000_000515_CTCGAAAAGTCGATAA-1`
# gr_subset_nan = gr_nan[a]
# hits = GenomicRanges::findOverlaps(gene_gr, gr_subset_nan)
# # unique(subjectHits(hits))
# 
# a = gr_nan$`annotation.UID-10X-Spikein-1pc_SLX-00000_000515_CTCGAAAAGTCGATAA-1`[subjectHits(hits)]
# table(a)
# 
# a = gr_nan$`annotation.UID-10X-Spikein-1pc_SLX-00000_000515_CTCGAAAAGTCGATAA-1`[hits]
# table(a)

# table(blacklist2, blacklist5)
# table(blacklist2 & blacklist3, blacklist5 & blacklist4)
# table((blacklist2 & blacklist3) | (blacklist5 & blacklist4))
# 
# table(blacklist2_X & blacklist3_X & blacklist4_X & blacklist5_X)


blacklist = blacklist_external | blacklist_external_X
blacklist_reads = c((blacklist2 & blacklist3 & blacklist4 & blacklist5), 
                    (blacklist2_X & blacklist3_X & blacklist4_X & blacklist5_X)) #, rep(TRUE, length(gr[seqnames(gr) == "Y"])))
blacklist_reads = is.na(blacklist_reads) | blacklist_reads

stopifnot(length(blacklist) == length(blacklist_reads))
blacklist_all = blacklist | blacklist_reads
table(blacklist_all)

# original number of QDNAseq blacklisted regions
table(is_notvalid)

# The 0.5 peak is X chromosome
# pos = median_read_depth(CN1[!blacklist_all,])
# ggplot(data=data.frame(x=pos)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 1.3)) + 
#   geom_vline(xintercept=c(lower_bound, upper_bound), color="red") + theme_minimal()
# pos = median_read_depth(CN2[!blacklist_all,])
# ggplot(data=data.frame(x=pos)) + geom_histogram(aes(x=x), bins = 100) + xlim(c(0, 1.3)) + 
#   geom_vline(xintercept=c(lower_bound, upper_bound), color="red") + theme_minimal()

# distributions
# ggplot(data=data.frame(x=as.vector(CN1[!idx_chr, ]@assayData$calls))) + geom_histogram(aes(x=x), binwidth = 1) + theme_cowplot() + coord_cartesian(xlim=c(0, 10))
ggplot(data=data.frame(x=as.vector(CN1[idx_chr, CN1@phenoData@data$used.reads > read_cutoff_x]@assayData$calls))) + geom_histogram(aes(x=x), binwidth = 1) + theme_cowplot() + coord_cartesian(xlim=c(0, 10))
ggplot(data=data.frame(x=as.vector(CN2[idx_chr, CN2@phenoData@data$used.reads > read_cutoff_x]@assayData$calls))) + geom_histogram(aes(x=x), binwidth = 1) + theme_cowplot() + coord_cartesian(xlim=c(0, 10))



# Visualize distribution across genome
copynumber = matrix(data=NA, nrow=dim(CN1)[[1]], ncol=1)
fdat = Biobase::featureData(CN1[,1])
fdat[["use"]] = rep(TRUE, length(fdat[["use"]]))
copynumber[!blacklist_all] = pos
obj = new("QDNAseqCopyNumbers",
          bins=fdat,
          # copynumber=copynumber,
          copynumber = matrix(data=ifelse(blacklist_all, 1.2, 0.2),nrow=dim(CN1)[[1]], ncol=1),
          phenodata=Biobase::pData(CN1[,1]))

Biobase::assayDataElement(obj, "segmented") = matrix(data=ifelse(blacklist_all, 0, 2), ncol=1)

QDNAseq::plot(obj[150000:160000,1], title="median read ratio across genome", ylab="read ratio", ylim=c(0.0,1.5), logTransform=FALSE)
QDNAseq::plot(obj[300000:309000,1], title="median read ratio across genome", ylab="read ratio", ylim=c(0.0,1.5), logTransform=FALSE)

abline(h=0.8, col="red")
abline(h=1.2, col="red")

# plotCounts(obj[52700:52900,1])

# GRanges from readCounts object
info = rownames(object)
seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
stopifnot(length(seqName) == length(chunkStart))
# remove chrX and chrY from blacklist
# start = min(which(startsWith(names(blacklist), "X:")))
# blacklist[start:length(blacklist)] = FALSE
gr = GRanges(seqName, IRanges(chunkStart, chunkEnd))
gr$blacklist = ifelse(blacklist_all, 0, 1)

rl <- mcolAsRleList(gr, "blacklist")

# NOTE uncomment in order to overwrite existing file
# saveRDS(rl, file="~/scAbsolute/data/blacklisting/extendedBlacklisting-GRCh37.RDS")
# and for GRCh38
# saveRDS(rl, file="~/scAbsolute/data/blacklisting/extendedBlacklisting-GRCh38.RDS")
saveRDS(rl, file="~/scAbsolute/data/blacklisting/extendedBlacklisting-crosstechnologies-GRCh37.RDS")

## example on how to use the gr object ####
# library(BSgenome.Hsapiens.UCSC.hg19)
# tiles <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg19), tilewidth=100001)
# gr2 = do.call("c", tiles)
# chromosomes = paste0("chr", c(as.character(seq(1, 22)), "X", "Y"))
# gr2 = keepSeqlevels(gr2, chromosomes, pruning.mode = "tidy")
# # gr3 = renameSeqlevels(gr2, c(as.character(seq(1, 22)), "X", "Y"))
# 
# ag <- binnedAverage(gr2, rl, "avg")
# ag

## Analyse existing blacklisting and example ####
# readCounts = readRDS("~/Data/blacklistingRegions/dummyReads.RDS")
# valid.old = !readRDS("~/Data/blacklistingRegions/extendedBlacklisting_old.RDS")$mcols.blacklist
# valid.new = fData(readCounts)[["use"]]
# 
# table(valid.old, valid.new)
# 
# # copynumber = matrix(data=NA, nrow=dim(readCounts)[[1]], ncol=1)
# fdat = Biobase::featureData(readCounts[,1])
# fdat[["use"]] = rep(TRUE, length(fdat[["use"]]))
# obj = new("QDNAseqCopyNumbers",
#           bins=fdat,
#           copynumber = matrix(data=ifelse(valid.new, 1.2, 0.2),nrow=dim(readCounts)[[1]], ncol=1),
#           phenodata=Biobase::pData(readCounts[,1]))
# 
# # Biobase::assayDataElement(obj, "segmented") = matrix(data=ifelse(blacklist_all, 0, 2), ncol=1)
# 
# QDNAseq::plot(obj[50000:55000,1], title="median read ratio across genome", ylab="read ratio", ylim=c(0.0,1.5), logTransform=FALSE)
# QDNAseq::plot(obj[,1], title="median read ratio across genome", ylab="read ratio", ylim=c(0.0,1.5), logTransform=FALSE)
