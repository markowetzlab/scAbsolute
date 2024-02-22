# Copyright 2022, Michael Schneider, All rights reserved.
rm(list=ls())
BASEDIR="~/scAbsolute/"
set.seed(2021)
library(tidyverse, quietly=TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts = FALSE)
library(gghalves, quietly=TRUE, warn.conflicts = FALSE)
library(cowplot, quietly=TRUE, warn.conflicts = FALSE)
library(ggpubr, quietly=TRUE, warn.conflicts = FALSE)
library(ggbeeswarm, quietly=TRUE, warn.conflicts = FALSE)
library(rtracklayer, quietly=TRUE, warn.conflicts = FALSE)
source("~/scUnique/R/visualize.R")
source("~/scAbsolute/R/core.R")
source("~/scAbsolute/R/visualization.R")

binSizeValue=100
files = c(
  paste0("~/Data/project1/30-scale/readCoverage/", 
         paste0(c(
           "UID-10X-Fibroblast-cell",
           "UID-10X-Fibroblast-nuclei",
           "UID-10X-Spikein-1pc",
           "UID-CCL-DIPLOID",
           "UID-DLP-SA928",
           "UID-JBL-NA12878"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files, 
function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

# load cellcycle info
df_cellcycle_DLP = readr::read_csv("~/mean-variance-model/data/cellcycle-info-DLP.csv", col_types = "ccccdd")
df_cellcycle_JBL = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv", col_types = "ccccdd")
df_meta = dplyr::bind_rows(df_cellcycle_DLP, df_cellcycle_JBL) %>%
  dplyr::mutate(name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

cycling = df_meta %>% dplyr::filter(cellcycle %in% c("G1", "D")) %>% dplyr::pull(name)
df1 = pd %>% dplyr::filter(name %in% cycling | SLX=="SLX-00000" | SLX %in% c("SLX-18431","SLX-18434", "SLX-20747", "SLX-20748") | 
                             SLX %in% c(
                                        # "SLX-12655", "SLX-12658", "SLX-12659", "SLX-13903", "SLX-13906", "SLX-13908",
                                        # "SLX-13909", "SLX-15390", "SLX-16839", "SLX-16971", "SLX-16976", "SLX-16979",
                                        # "SLX-17692", "SLX-17728", "SLX-17747", "SLX-17751", "SLX-18615", "SLX-18616",
                                        # "SLX-18816", "SLX-18817", "SLX-18818", "SLX-18819", "SLX-18820", "SLX-19311",
                                        "SLX-19312", "SLX-19314", "SLX-19315"))

# ## cellcycle analysis with MEDICC => on cluster
CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
object1 = CN[, colnames(CN) %in% df1$name]
valid_bins = binsToUseInternal(CN)
n_bins = dim(CN)[[1]]
rnames = rownames(CN)

# # valid_50 
# t50 = readRDS("~/Data/project1/30-scale/50/scaling/UID-CIOV1_50.rds")
# valid_50 = binsToUseInternal(t50)
# t100 = readRDS("~/Data/project1/30-scale/100/scaling/UID-CIOV1_100.rds")
# valid_100 = binsToUseInternal(t100)

## Split data up in X, Y and rest -> only create filters for rest
valid = binsToUseInternal(CN)
x_data = startsWith(rownames(object1), "X:") #& valid
y_data = startsWith(rownames(object1), "Y:") #& valid
data = !startsWith(rownames(object1), "Y:") & !startsWith(rownames(object1), "X:") #& valid_50
index = data

result = list()
result[[1]]=list();result[[2]]=list();result[[3]]=list();result[[4]]=list();
counter = 0
gc = object1@featureData@data$gc[index]
map = object1@featureData@data$map[index]
for(SLX_set in list( list("SLX-00000"), list("SLX-A73044A", "SLX-A90553C"),
                     list("SLX-18431","SLX-18434","SLX-20747","SLX-20748"),
                     list(
                          "SLX-12655", "SLX-12658", "SLX-12659", "SLX-13903", "SLX-13906", "SLX-13908",
                          "SLX-13909", "SLX-15390", "SLX-16839", "SLX-16971", "SLX-16976", "SLX-16979",
                          "SLX-17692", "SLX-17728", "SLX-17747", "SLX-17751", "SLX-18615", "SLX-18616",
                          "SLX-18816", "SLX-18817", "SLX-18818", "SLX-18819", "SLX-18820", "SLX-19311",
                          "SLX-19312", "SLX-19314", "SLX-19315"))){
  print(SLX_set)
  counter = counter + 1
  subset = object1@phenoData@data %>% tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
    dplyr::filter(SLX %in% SLX_set) %>% dplyr::pull(name)
  
  if(length(subset) == 0){
    next
  }
  object = object1[, subset]
  rown = rownames(object)[index]
  print(dim(object))
  rpc = apply(object@assayData$probdloss[index,], 2, mean, na.rm=TRUE)
  ratio = (object@assayData$probdloss[index,] / rpc)
  r = apply(ratio, 1, median, na.rm=TRUE)
  
  ratio_X = (object@assayData$probdloss[x_data,] / rpc)
  ratio_Y = object@assayData$probdloss[y_data,]
  r_X = apply(ratio_X, 1, mean, na.rm=TRUE)
  r_Y = apply(ratio_Y, 1, function(x) sum(x != 0, na.rm=TRUE))
  
  hist(r, breaks=100)
  print(summary(r))
  
  result[[counter]] = list("SLX_set"=SLX_set, "r"=r, "rX"=r_X, "rY"=r_Y)
}

# rm(CN); rm(object1)

# Copied from http://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

n_items = 4
udat_X = dplyr::tibble(
  position=rep(rownames(object)[x_data], times=n_items),
  index=rep(seq(1, length(index))[x_data], times=n_items),
  gc=rep(object1@featureData@data$gc[x_data], times=n_items),
  map=rep(object1@featureData@data$map[x_data], times=n_items),
  qdnaseq_valid = rep(valid_bins[x_data], times=n_items),
  ratio=c(result[[1]][["rX"]], result[[2]][["rX"]], result[[3]][["rX"]], result[[4]][["rX"]]))
udat_X$dataset = base::rep(c("10X", "DLP", "JBL", "CCL"), each=sum(x_data))
datX = udat_X %>% dplyr::filter(ratio > 0) %>% dplyr::group_by(dataset) %>%
  dplyr::group_modify(~ { .x %>% dplyr::mutate(density = get_density(ratio, gc))}) %>% dplyr::ungroup()

udat_Y = dplyr::tibble(
  position=rep(rownames(object)[y_data], times=n_items),
  index=rep(seq(1, length(index))[y_data], times=n_items),
  gc=rep(object1@featureData@data$gc[y_data], times=n_items),
  map=rep(object1@featureData@data$map[y_data], times=n_items),
  qdnaseq_valid = rep(valid_bins[y_data], times=n_items),
  ratio=c(result[[1]][["rY"]], result[[2]][["rY"]], result[[3]][["rY"]], result[[4]][["rY"]]))
udat_Y$dataset = base::rep(c("10X", "DLP", "JBL", "CCL"), each=sum(y_data))
datY = udat_Y %>% dplyr::filter(ratio > 0) %>% dplyr::group_by(dataset) %>%
  dplyr::group_modify(~ { .x %>% dplyr::mutate(density = get_density(ratio, gc), avg=median(ratio))}) %>% dplyr::ungroup()


## X chromosome analysis

gg <- ggplot(data = datX %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL")), aes(x = gc, y = ratio)) +
  geom_point(cex = 0.5) + 
  ylab("Median of normalized reads per bin") + 
  facet_wrap(~dataset) +
  geom_density_2d(size = 0.8, n = 100, color="red", bins=50) + theme_cowplot() +
  geom_smooth(aes(x=gc, y=ratio), color="blue", model="gam")
gg


gb <- ggplot_build(gg)
contour_levels <- unique(gb[["data"]][[2]][["level"]])

datX = datX %>% dplyr::mutate(outlier = density < contour_levels[1] | ratio < 0.52)
g1 = ggplot(data=datX %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL"))) + 
  facet_wrap(~dataset) + 
  geom_point(aes(x=gc, y=ratio, color=outlier)) + 
  scale_colour_manual(values = c("FALSE" = "#00BFC4","TRUE" = "#F8766D")) + 
  theme_pubclean()
g1
#ggpubr::ggexport(g1, filename="~/scAbsolute/figures/Supplement-QC-Filtering-QC_sex_X.pdf")

decisionX = datX %>% mutate(outlier_density = density < contour_levels[1], outlier_ratio = ratio < 0.52)
outlier_10X = decisionX %>% dplyr::filter(dataset == "10X" & (outlier_density | outlier_ratio)) %>% dplyr::pull(position)
outlier_DLP = decisionX %>% dplyr::filter(dataset == "DLP" & (outlier_density | outlier_ratio)) %>% dplyr::pull(position)
outlier_JBL = decisionX %>% dplyr::filter(dataset == "JBL" & (outlier_density | outlier_ratio)) %>% dplyr::pull(position)
length(intersect(outlier_10X, outlier_DLP))
length(outlier_10X)
length(outlier_DLP)
length(outlier_JBL)
decisionX$valid = !(decisionX$position %in% union(union(outlier_10X, outlier_DLP), outlier_JBL))


## Y chromosome analysis

g2 = ggplot(data=datY %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL"))) + 
  facet_wrap(~dataset) +
  geom_point(aes(x=gc, y=ratio)) + 
  geom_point(aes(x=gc, y=avg), color="red") + 
  theme_pubclean()

decisionY = datY %>% mutate(outlier_ratio = (ratio < (avg - 40)) | (ratio > avg + 40))

g2 = ggplot(data=decisionY %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL"))) +
  geom_point(aes(x=gc, y=ratio, color=outlier_ratio)) + 
  scale_colour_manual(values = c("FALSE" = "#00BFC4","TRUE" = "#F8766D")) + 
  facet_wrap(~dataset) + theme_pubclean()
g2
#ggpubr::ggexport(g2, filename="~/scAbsolute/figures/Supplement_QC-Filtering-QC_sex_Y.pdf")


outlier_10X = decisionY %>% dplyr::filter(dataset == "10X" & (outlier_ratio)) %>% dplyr::pull(position)
outlier_DLP = decisionY %>% dplyr::filter(dataset == "DLP" & (outlier_ratio)) %>% dplyr::pull(position)
outlier_JBL = decisionY %>% dplyr::filter(dataset == "JBL" & (outlier_ratio)) %>% dplyr::pull(position)
length(intersect(outlier_10X, outlier_DLP))
length(outlier_10X)
length(outlier_DLP)
length(outlier_JBL)
decisionY$valid = !(decisionY$position %in% union(union(outlier_10X, outlier_DLP), outlier_JBL))


## Export files

blacklist_regions = gtools::mixedsort(unique(dplyr::bind_rows(decisionX, decisionY) %>% dplyr::filter(!valid, dataset %in% c("10X", "DLP")) %>% dplyr::pull(position)))
use = rep(TRUE, n_bins)
names(use) = rnames 
use[blacklist_regions] = FALSE
# only extending to X and Y areas - we don't have blacklist information here
use[index] = TRUE

# create GRanges object for blacklisting
gr = createGR(object[,1])
exclude_regions = gr[!use]
exclude_regions = GenomicRanges::reduce(exclude_regions)

## add centromere and telomere regions
centromers = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/centromere.tsv"),
                             col_names = c("chromosome", "start", "end", "arm", "acen"), col_types = "ciicc")
# custom filter extending regions around chromosome ends
conflicting_regions = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/extended_telomere_filter.tsv"),
                                      col_names = c("chromosome", "start", "end"),
                                      col_types = "cii")
subject1 = GRanges(seqnames=gsub("chr", "", centromers$chromosome),
                   ranges = IRanges(start=centromers$start+1, end=centromers$end), seqinfo = seqinfo(gr))
subject2 = GRanges(seqnames=gsub("chr", "", conflicting_regions$chromosome),
                   ranges = IRanges(start=conflicting_regions$start+1, end=conflicting_regions$end), seqinfo=seqinfo(gr))
subject = sort(c(exclude_regions, subject1, subject2))

subject = GenomicRanges::reduce(subject[subject@seqnames %in% c("X", "Y")])
sum(subject@ranges@width) / sum(gr[seqnames(gr) %in% c("X", "Y")]@ranges@width)


fig_QC_Sex1 = ggpubr::ggarrange(
  ggpubr::ggarrange(ggparagraph(text=" ", face = "italic", size = 16, color = "white"), gg + rremove("xlab") + rremove("ylab"),
                    ggparagraph(text=" ", face = "italic", size = 16, color = "white"), g1 + theme_cowplot() + theme(strip.background = element_blank(),strip.text.x = element_blank()) + rremove("legend") + rremove("xlab") + ylab("Median of normalized reads per bin") + rremove("ylab"),
                    ncol=2, nrow=2, labels=c(" ", "a", " ", "b"), widths = c(0.1, 4), hjust=-1.0, vjust=1.5, label.y=c(1, 1, 1, 1.07)),
                    g2 + theme_cowplot() + theme(strip.background = element_blank(),strip.text.x = element_blank()) + rremove("xlab") +
                                 labs(color="Outlier bin ") + ylab("Total reads") + xlab("GC content") + theme(legend.background = element_rect(fill = "white"),
                                                               legend.key = element_rect(fill = "white", color = NA),
                                                               legend.key.width = unit(0.5,"cm"),
                                                               legend.position = c(0.85, 0.8)
                                                               ),
                   labels=c("", "c"), label.x=0.03, label.y=1.05, ncol=1, nrow=2, heights = c(2, 1))
fig_QC_Sex = ggpubr::annotate_figure(fig_QC_Sex1,
  left = text_grob("Median of normalized reads", rot = 90, size=14, vjust=2.7, hjust=0.1, color="black"),
  bottom = text_grob("GC content", rot = 0, size=14, vjust=0.0, hjust=0.0, color="black"))
fig_QC_Sex

# NOTE: file requires manual export for some reason to keep proper aspect ratio
# cowplot::ggsave2(filename = "~/scAbsolute/figures/Supplement-QC-Sex.pdf", fig_QC_Sex, width=4.8, height = 2.7)
ggpubr::ggexport(fig_QC_Sex, filename = "~/scAbsolute/figures/Sup_QC-Sex_fix.pdf", width=9, height = 6)

stop("Manually uncomment file export to overwrite existing file!")
# rtracklayer::export.bed(subject, "~/scAbsolute/data/blacklisting/final_exclude_regions_hg19_sex.bed")
