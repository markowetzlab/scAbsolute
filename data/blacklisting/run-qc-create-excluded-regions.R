# Copyright 2022, Michael Schneider, All rights reserved.
## Reduce overall noise level further by identifying outliers (on diploid data and three technologies)

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
                                        "SLX-17692", "SLX-17728", "SLX-17747", "SLX-17751", "SLX-18615", "SLX-18616",
                                        "SLX-18816", "SLX-18817", "SLX-18818", "SLX-18819", "SLX-18820", "SLX-19311",
                                        "SLX-19312", "SLX-19314", "SLX-19315") |
                             SLX %in% c("SLX-A90689B", "SLX-A90689C", "SLX-A90694B"))

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
x_data = startsWith(rownames(object1), "X:") & valid
y_data = startsWith(rownames(object1), "Y:") & valid
data = !startsWith(rownames(object1), "Y:") & !startsWith(rownames(object1), "X:") #& valid_50
index = data

result = list()
result[[1]]=list();result[[2]]=list();result[[3]]=list();result[[4]]=list();
counter = 0
gc = object1@featureData@data$gc[index]
map = object1@featureData@data$map[index]
for(SLX_set in list( list("SLX-00000"), list("SLX-A73044A", "SLX-A90553C", "SLX-A90689B", "SLX-A90689C", "SLX-A90694B"),
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

  hist(r, breaks=100)
  print(summary(r))
  
  result[[counter]] = list("SLX_set"=SLX_set, "r"=r)
}

# rm(CN); rm(object1)

n_items = 4
udat = dplyr::tibble(
  position=rep(rownames(object)[index], times=n_items),
  index=rep(seq(1, length(index))[index], times=n_items),
  gc=rep(gc, times=n_items), 
  map=rep(map, times=n_items),
  qdnaseq_valid = valid_bins[index],
  ratio=c(result[[1]][["r"]], result[[2]][["r"]], result[[3]][["r"]], result[[4]][["r"]]))
udat$dataset = base::rep(c("10X", "DLP", "JBL", "CCL"), each=length(gc))

dat = udat %>% dplyr::group_by(dataset) %>%
  dplyr::group_modify(~ { .x %>% dplyr::mutate(extrema = (ratio > 4.0 | ratio < 0.10 | map < 70)) %>% mutate(m=mean(ratio[!extrema]),
    m_l = mean(ratio[!extrema], na.rm=TRUE) - 1.96 * sd(ratio[!extrema], na.rm=TRUE),
    m_u = mean(ratio[!extrema], na.rm=TRUE) + 1.96 * sd(ratio[!extrema], na.rm=TRUE))}) %>% dplyr::ungroup()
  
# Copied from http://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
dat1 = dat %>% dplyr::filter(!extrema) %>% dplyr::group_by(dataset) %>%
  dplyr::group_modify(~ { .x %>% dplyr::mutate(density = get_density(ratio, gc))}) %>% dplyr::ungroup()
dat1 = dat1 %>% dplyr::group_by(dataset) %>%
  dplyr::group_modify(function(x,y){
    model <- mgcv::gam(ratio ~ s(gc, bs = "cs"), data=x)
    xrange <- range(x$gc)
    xseq <- seq(from=xrange[1], to=xrange[2], length=100)
    pred <- predict(model, newdata = data.frame(gc = x$gc))
    deviance = x$ratio - pred
    x$deviance = deviance;
    x$d_u = pred + 2 * sd(deviance);
    x$d_uu = pred + 3 * sd(deviance);
    x$d_l = pred - 2 * sd(deviance);
    x$d_ll = pred - 3 * sd(deviance);
    return(x)}) %>% dplyr::ungroup()


dat2 = dat %>% dplyr::filter(extrema) %>% dplyr::mutate(density = NA, deviance=NA, d_u=NA, d_uu=NA, d_l=NA, d_ll=NA)
dat = dplyr::bind_rows(dat1, dat2)
stopifnot(dim(dat)[[1]] == dim(udat)[[1]])

# First create plot with geom_density
gg <- ggplot(data = dat %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL"), !extrema), aes(x = gc, y = ratio)) +
  geom_point(cex = 0.5) + 
  facet_wrap(~dataset) +
  geom_density_2d(size = 0.8, n = 100, color="red", bins=150) + theme_cowplot() + 
  geom_smooth(aes(x=gc, y=ratio), color="blue", model="gam") +
  theme_pubclean(base_size=20)+
  ylab("Median of normalized reads per bin")
gg


gg <- ggplot(data = dat %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL"), !extrema),
             aes(x = gc, y = ratio, color=as.factor(ratio > d_u | ratio < d_l))) +
  geom_point(data = dat %>% dplyr::filter(dataset %in% c("10X", "DLP"), !extrema), aes(x = gc, y = ratio, color=as.factor(ratio > d_u | ratio < d_l)), cex = 0.5) + 
  geom_point(data = dat %>% dplyr::filter(dataset %in% c("JBL"), !extrema), aes(x = gc, y = ratio, color=as.factor(ratio > d_uu | ratio < d_ll)), cex = 0.5) + 
  facet_wrap(~dataset) +
  scale_colour_manual( values = c("TRUE" = "black","FALSE" = "lightblue")) + 
  geom_density_2d(size = 0.8, n = 100, color="red", bins=150) +
  geom_smooth(aes(x=gc, y=ratio), color="blue", method="gam") + 
  geom_point(data=dat %>% dplyr::filter(dataset %in% c("10X", "DLP"), !extrema) %>%
               dplyr::filter(gc > 32 & gc < 60), aes(x=gc, y=d_u), color="black", size=0.1) + 
  geom_point(data=dat %>% dplyr::filter(dataset %in% c("10X", "DLP"), !extrema) %>%
               dplyr::filter(gc > 32 & gc < 60), aes(x=gc, y=d_l), color="black", size=0.1) + 
  geom_point(data=dat %>% dplyr::filter(dataset %in% c("JBL"), !extrema) %>%
               dplyr::filter(gc > 32 & gc < 60), aes(x=gc, y=d_uu), color="black", size=0.1) + 
  geom_point(data=dat %>% dplyr::filter(dataset %in% c("JBL"), !extrema) %>%
               dplyr::filter(gc > 32 & gc < 60), aes(x=gc, y=d_ll), color="black", size=0.1) + 
  ylab("Median of normalized reads per bin") + 
  xlab("GC content") + theme_pubclean(base_size = 20) +
  theme(legend.position = "none") + theme(strip.background = element_blank(),
                                          strip.placement = "outside",
                                          strip.text.x = element_text(size = 14))

gg + coord_cartesian(ylim=c(0.5, 1.5)) 
gg

# ggplot(data=dat2) + geom_density(aes(x=deviance), bins=100) + facet_wrap(~dataset) + theme_pubclean()
#ggpubr::ggexport(gg + coord_cartesian(ylim=c(0.5, 1.5)), filename="~/scAbsolute/figures/Supplement-QC-density.pdf")


# Extract levels denoted by contours by going into the
#   ggplot build object. I found these coordinates by
#   examining the object in RStudio; Note, the coordinates
#   would change if the layer order were altered.
gb <- ggplot_build(gg)
# NOTE this value will change if you edit gg
contour_levels <- unique(gb[["data"]][[3]][["level"]])

# Add layer that relies on given contour level
gg2 <- gg +
  geom_point(data = dat %>% dplyr::filter(!extrema, density < contour_levels[1]),
             color = "blue", size = 0.3) + theme_cowplot()
gg2

decision = dat %>% mutate(outlier_density = density < contour_levels[1],
                          outlier_ratio = ratio < d_l | ratio > d_u,
                          outlier_ratio2 = ratio < d_ll | ratio > d_uu, )
outlier_10X = decision %>% dplyr::filter(dataset == "10X" & (outlier_density | outlier_ratio | extrema)) %>% dplyr::pull(position)
outlier_DLP = decision %>% dplyr::filter(dataset == "DLP" & (outlier_density | outlier_ratio | extrema)) %>% dplyr::pull(position)
outlier_JBL = decision %>% dplyr::filter(dataset == "JBL" & (outlier_density | outlier_ratio2 | extrema)) %>% dplyr::pull(position)
length(intersect(outlier_10X, outlier_DLP))
length(outlier_10X)
length(outlier_DLP)
length(outlier_JBL)
decision$valid = !(decision$position %in% union(union(outlier_10X, outlier_DLP), outlier_JBL))
# decision$valid = !(decision$position %in% union(outlier_10X, outlier_DLP))

## Excurs: plot this as VennDiagram
# Load library
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
v = venn.diagram(
  x = list(outlier_DLP, outlier_JBL, outlier_10X),
  category.names = c("DLP" , "JBL" , "10X"),
  filename = '~/scAbsolute/figures/Supplement-QC-Venn.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 920 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


tbl = table(decision %>% dplyr::filter(dataset == "10X") %>% dplyr::pull(valid))
sum(tbl)
tbl["FALSE"] / 30970

p = ggplot(data = decision %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL")))  +
  geom_point(data = decision %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL")) %>%
    dplyr::filter(gc > 32 & gc < 60), aes(x=gc, y=ratio, color=!valid, position=position), size=0.1) + 
  facet_wrap(~dataset) + 
  geom_point(data=decision %>% dplyr::filter(dataset %in% c("10X", "DLP")) %>% dplyr::filter(gc > 32 & gc < 60), 
             aes(x=gc, y=d_u), color="black", size=0.1) + 
  geom_point(data=decision %>% dplyr::filter(dataset %in% c("10X", "DLP")) %>% dplyr::filter(gc > 32 & gc < 60), 
             aes(x=gc, y=d_l), color="black", size=0.1) + 
  geom_point(data=decision %>% dplyr::filter(dataset %in% c("JBL")) %>% dplyr::filter(gc > 32 & gc < 60), 
             aes(x=gc, y=d_uu), color="black", size=0.1) + 
  geom_point(data=decision %>% dplyr::filter(dataset %in% c("JBL")) %>% dplyr::filter(gc > 32 & gc < 60), 
             aes(x=gc, y=d_ll), color="black", size=0.1) + 
  scale_colour_manual(values = c("FALSE" = "#00BFC4","TRUE" = "#F8766D")) + 
  coord_cartesian(ylim=c(0,2)) +
  theme_pubclean(base_size = 20) +
  guides(color = guide_legend("Outlier bin ", override.aes = list(size = 3))) +
  theme(
    legend.background=element_blank(),
    #legend.box.background = element_rect(color="white", size=2),
    legend.key = element_rect(colour = "white", fill=NA),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(face = "bold")
  ) + 
  ylab("Median of normalized reads per bin") + xlab("GC content") + 
  theme(strip.background = element_blank(),
                                          strip.placement = "outside",
                                          strip.text.x = element_text(size = 14))
# ggplotly(p)
p
#ggpubr::ggexport(p, filename="~/scAbsolute/figures/Supplement-QC-Filtering-decision.pdf")


p2 = ggplot(data = decision %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL")))  +
  geom_point(data = decision %>% dplyr::filter(dataset %in% c("10X", "DLP", "JBL")),
             aes(x=map, y=ratio, color=!valid, position=position), size=0.1) + 
  facet_wrap(~dataset) + 
  scale_colour_manual(values = c("FALSE" = "#00BFC4","TRUE" = "#F8766D")) + 
  #scale_colour_manual( values = c("TRUE" = "blue","FALSE" = "red")) + 
  coord_cartesian(ylim=c(0,2)) + 
  ylab("Median of normalized reads per bin") + xlab("Mappability") + theme_pubclean(base_size=20) + 
  theme(legend.position = "none") + theme(strip.background = element_blank(),
                                          strip.placement = "outside",
                                          strip.text.x = element_text(size = 14))
# ggplotly(p)
p2
#ggpubr::ggexport(p2, filename="~/scAbsolute/figures/Supplement-QC-Filtering-Mappability.pdf")

table(decision %>% dplyr::filter(dataset=="10X") %>% dplyr::select(valid))

# keep these positions in eye for JBL data
# a = decision %>% dplyr::filter(dataset == "JBL", ratio < m_l, valid) %>% dplyr::pull(position)
# warning("JBL checklist")
# warning(a)


## Export and merge with some other criteria
blacklist_regions = gtools::mixedsort(unique(decision %>% dplyr::filter(!valid, dataset %in% c("10X", "DLP")) %>% dplyr::pull(position)))
use = rep(TRUE, n_bins)
names(use) = rnames 
use[blacklist_regions] = FALSE
# only extending to X and Y areas - we don't have blacklist information here
use[!valid_bins] = FALSE

tbl = table(use)
stopifnot(sum(tbl) == 30970)
tbl["FALSE"] / 30970

# create GRanges object for blacklisting
gr = createGR(object[,1])
exclude_regions = gr[!use]
exclude_regions = GenomicRanges::reduce(exclude_regions)

## add centromere and telomere regions
centromers = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/centromere.tsv"),
                             col_names = c("chromosome", "start", "end", "arm", "acen"), col_types = "ciicc")
# custom filter extending regions around chromosome ends
conflicting_regions_telomere = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/extended_telomere_filter.tsv"),
                                      col_names = c("chromosome", "start", "end"),
                                      col_types = "cii")
conflicting_regions_manual = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/manual_filter.tsv"),
                                      col_names = c("chromosome", "start", "end"),
                                      col_types = "cii")
conflicting_regions = dplyr::bind_rows(conflicting_regions_telomere, conflicting_regions_manual)

subject1 = GRanges(seqnames=gsub("chr", "", centromers$chromosome),
                   ranges = IRanges(start=centromers$start+1, end=centromers$end), seqinfo = seqinfo(gr))
subject2 = GRanges(seqnames=gsub("chr", "", conflicting_regions$chromosome),
                   ranges = IRanges(start=conflicting_regions$start+1, end=conflicting_regions$end), seqinfo=seqinfo(gr))
subject = sort(c(exclude_regions, subject1, subject2))

# subject@seqnames
subject = GenomicRanges::reduce(subject[!(subject@seqnames %in% c("X", "Y"))])
sum(subject@ranges@width) / sum(gr@ranges@width)

library(png)
img_venn <- png::readPNG("~/scAbsolute/figures/Supplement-QC-Venn.png")
venn_plot <- ggplot() + background_image(img_venn) + 
  theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

p4 <- ggplot(data.frame(l = "Median of normalized read counts", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() + theme_pubclean() + 
  coord_cartesian(clip = "off")

p_legend = cowplot::get_legend(p)
S_fig_QC = ggpubr::ggarrange(ggpubr::ggarrange(gg + rremove("ylab"), 
                                               ggpubr::ggarrange(venn_plot + theme_nothing() + theme(plot.margin = unit(c(10,5,5,10), "pt")), p_legend, nrow=2, heights = c(2, 1)), widths = c(3, 1), nrow=1, labels=c("a", "b"),
                                               font.label = list(size=20)), 
                  ggpubr::ggarrange(p + rremove("ylab") + rremove("legend"),
                                    p2 + rremove("ylab"), widths = c(1, 1), nrow=2, labels=c("c", "d"),
                                    font.label = list(size=20)), 
                  ncol = 1, nrow = 2, heights = c(1, 2))

# Annotate the figure by adding a common labels
S_fig_QC = annotate_figure(S_fig_QC,
                left = text_grob("Median of normalized read counts", rot = 90, size=22))
S_fig_QC

ggpubr::ggexport(S_fig_QC,filename = "~/scAbsolute/figures/Sup_QC-Autosomes_fix.pdf", height = 9, width = 16)

stop("Manually uncomment file export to overwrite existing file!")
#rtracklayer::export.bed(subject, "~/scAbsolute/data/blacklisting/final_exclude_regions_hg19_1-22.bed")

reads_lost = apply(CN@assayData$probdloss, 2, function(x){
  return(sum(x[!use], na.rm=TRUE) / sum(x, na.rm=TRUE))})

# estimate number of reads removed

# 
# ## IDENTIFY BINS THAT WE WANT TO KEEP -> BINS with CANCER GENES
# gr = createGR(object[,1])
# BASEDIR="~/scAbsolute/"
# 
# # cancer_genes = readr::read_tsv(file=file.path(BASEDIR, "data/blacklisting/NCG6_cancergenes.tsv"), 
# #                                col_names=TRUE, col_types="dcdcccc")
# # cancer_genes = readr::read_csv(file=file.path(BASEDIR, "data/blacklisting/Cosmic-Census-Tier1-2021.csv"),
# #                                col_names = TRUE)
# cancer_genes = readr::read_tsv(file=file.path(BASEDIR, "data/blacklisting/NCG7-canconical-drivers.txt"),col_names="Gene Symbol") 
# gene_annotation = read_tsv(file.path(BASEDIR, "data/blacklisting/gencode.v19.tsv"), 
#                            col_names = c("chromosome", "annotation", "start", "end", "gene_id", "gene_name"),
#                            col_types = "ccddcc")
# 
# # do not remove regions containing cancer genes
# gene_annotation = gene_annotation[gene_annotation$gene_name %in% cancer_genes$`Gene Symbol`, ]
# gene_gr = GRanges(seqnames=gsub("chr", "", gene_annotation$chromosome), ranges = IRanges(start=gene_annotation$start, end=gene_annotation$end), 
#                   gene_id=gene_annotation$gene_id, gene_name=gene_annotation$gene_name)
# 
# hitsBad = GenomicRanges::findOverlaps(gene_gr, gr)
# 
# hits = GenomicRanges::findOverlaps(gene_gr, gr)
# keep_bins = rep(FALSE, times=length(gr))


# 
# a = sort(gene_gr[queryHits(hits[which(!use[keep_bins]),])]$gene_name)
# contentious_bins = gr[subjectHits(hits[which(!use[keep_bins])])]
# end=(contentious_bins@ranges@start + contentious_bins@ranges@width -1)
# contentious_bins = paste0(contentious_bins@seqnames, ":",
#                           contentious_bins@ranges@start, "-", end)
# 
# which(blacklist_regions %in% contentious_bins)
# which(rownames(CN) %in% contentious_bins)
# 
# table(use[keep_bins])
# cancer_gene_bins = rownames(CN)[keep_bins]
# # valid1 = !startsWith(rownames(object1), "Y:") & !startsWith(rownames(object1), "X:")
# # table(valid1)
# # 
# # table(keep_bins, valid)
# decision = decision %>% dplyr::mutate(is_cancer_gene = position %in% cancer_gene_bins,
#                                       is_contentious = position %in% contentious_bins)
# p = ggplot(data = decision) +
#   # geom_point(aes(x=gc, y=ratio, color=valid, shape=is_cancer_genes, position=position), size=0.1) + 
#   geom_point(data=decision %>% dplyr::filter(is_contentious), aes(x=map, y=ratio, position=position), size=0.5, color="black") + 
#   geom_hline(aes(yintercept=m_l), color="cyan") + 
#   geom_hline(aes(yintercept=m_u), color="cyan") + 
#   facet_wrap(~dataset) + coord_cartesian(ylim=c(0,2)) +
#   ylab("Median of normalized read counts") + theme_pubclean()
# # ggplotly(p)
# p
# 
# table(decision %>% dplyr::filter(dataset=="10X") %>% dplyr::select(valid))
# 









# ggplot(data = dat) +
#   geom_point(aes(x=map, y=ratio, color=valid), size=0.1) + facet_wrap(~dataset) + theme_pubclean()
# 
# ggplot(data = dat) +
#   geom_density_2d(aes(x=gc, y=ratio, color=map), color="red", adjust=4) +
#   geom_point(aes(x=gc, y=ratio), size=0.1, alpha=0.1) +
#   facet_wrap(~dataset) + theme_pubclean()


# ## NOT USED ====
# ## Extract outliers

# 
# 
# # plot(gc[valid_index], r[valid_index], col=r) + facet_wrap(slx)
# # plot(map[valid_index], r[valid_index]) #, use="complete.obs")
# 
# 
# # log transform seg level for better stability
# # dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
# #                  reads_normalized=as.vector(reads)/as.vector(segs),
# #                  gc=gc, map=map)
# 
# ## Build models
# formula.gam = formula(ratio ~ s(gc, map))  
# newdata=data.frame(ratio=dat$ratio, gc=dat$gc, map=dat$map, dataset=dat$dataset)
# dat_model = newdata %>% dplyr::filter(ratio > 1) %>% dplyr::filter(dataset == "10X")
# 
# model <- mgcv::gam(formula.gam, data=dat_model) #, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
# summary(model)
# hist(residuals(model), breaks=100)
# which(residuals(model) < 0.10)
# 
# terms = predict(model, newdata=newdata, type="terms")
# 
# 
# 
# # library(VennDiagramm)
# 
# 
# quantile(r, c(0.01, 0.05, 0.1, 0.9, 0.95, 0.99), na.rm=TRUE)
# table(r < mean(r, na.rm=T) - 2 * sd(r, na.rm=T))
# cor(gc[1:28823], r[1:28823],use = "complete.obs")
# 
# rownames(object)[r < 1.62 & !is.na(r)]
# 


# ratio = apply(object@assayData$calls / object@phenoData@data$rpc, 1, function(x) sum(x))

# ## MEDICC analysis - on full data ====
# df = getMediccTable(CN) # NOTE, by default accessing copynumber slot!
# print("Segments for MEDICC run")
# print(nrow(df %>% dplyr::filter(sample_id == unique(df$sample_id)[[1]])))
# print("Samples for MEDICC run")
# print(dim(CN)[[2]])
# 
# require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
# reticulate::py_discover_config()
# reticulate::source_python(file.path("~/scUnique/interface_medicc2.py"), convert=TRUE)
# 
# start_time_medicc <- Sys.time()
# results = interface_medicc(df)
# end_time_medicc <- Sys.time()
# print(paste0("Run MEDICC algorithm ",  difftime(end_time_medicc,start_time_medicc,units="mins")))
# 
# tree = results[[1]]; input_df = results[[2]]; output_df = results[[3]];
# summary = results[[4]]; pdms = results[[5]];
# uniqueEvents = results[[6]]; newick_tree = results[[7]]
# 
# medicc_tree <- read.dendrogram(text = newick_tree)
# # NOTE remove diploid
# tree <- phylogram::prune(medicc_tree, pattern = "diploid")
# 
# 


# 
# # unique(subjectHits(hits))
# # table(valid[queryHits(hits)])
# # table(valid[])
# # fO = findOverlaps(query, subject, type="within")
# 
# # Manually updating centromer regions ====
# # check centromer regions are okay -> and working, commit code
# # there are some other regions -> check with control data, identify regions and remove from valid
# ## analyse the variance
# CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
# valid = binsToUseInternal(CN)
# table(valid)
# 
# 
# dp = df %>% dplyr::filter(UID=="UID-JBL-NA12878") %>% dplyr::filter(cellcycle %in% c("G1"), abs(ploidy - round(ploidy)) < 0.05)
# # ggplot(data=dp) + geom_point(aes(x=cellcycle, y=hmm.alpha))
# 
# CN2 = CN[, colnames(CN) %in% dp$name]
# obj = CN2[, 5]
# plotCopynumberHeatmap(CN2, cluster_rows=FALSE)
# 
# gr = createGR(obj)
# customFilter = c(2488:2493, 4926, 6907, 8811:8815, 10631, 13932, 20855, 22002, 22003, 23069:23076, 25824, 25825, 27197, 28826:28850, 30369:30400)
# gr_export = reduce(gr[customFilter, ])
# # export.bed(gr_export, "~/region_filter.tsv")
# # 2493  4925  6906  8818 10628 12340 13932 15396 16809 18165 19516 20855 22007 23081 24107 25011 25823 26604 27196 27827
# # 28309 28823 30376 30970
# 
# plotCopynumberHeatmap(CN2[2480:2500,], cluster_rows=FALSE)
# CN2[2488:2493,]@assayData$copynumber
# apply(CN2[2485:2493,]@assayData$calls / CN2[2485:2493,]@phenoData@data$rpc, 1, median)
# # 2488-2493
# 
# plotCopynumberHeatmap(CN2[4924:4930,], cluster_rows=FALSE)
# CN2[4926:4930,]@assayData$copynumber
# apply(CN2[4924:4930,]@assayData$calls / CN2[4924:4930,]@phenoData@data$rpc, 1, median)
# # 4926
# 
# plotCopynumberHeatmap(CN2[6906:6920,], cluster_rows=FALSE)
# CN2[6907:6910,]@assayData$copynumber
# apply(CN2[6906:6910,]@assayData$calls / CN2[6906:6910,]@phenoData@data$rpc, 1, median)
# # 6907
# 
# plotCopynumberHeatmap(CN2[8811:8815,], cluster_rows=FALSE)
# CN2[8812:8818,]@assayData$copynumber
# apply(CN2[8810:8820,]@assayData$calls / CN2[8810:8820,]@phenoData@data$rpc, 1, median)
# # 8811 - 15
# 
# plotCopynumberHeatmap(CN2[10631:10640,], cluster_rows=FALSE)
# rownames(CN2)[10631]
# CN2[10631,]@assayData$copynumber
# # 10631
# 
# plotCopynumberHeatmap(CN2[12335:12350,], cluster_rows=FALSE)
# CN2[12335:12350,]@assayData$copynumber
# apply(CN2[12335:12350,]@assayData$calls / CN2[12335:12350,]@phenoData@data$rpc, 1, median)
# # fine
# 
# plotCopynumberHeatmap(CN2[13920:13950,], cluster_rows=FALSE)
# CN2[13930:13935,]@assayData$copynumber
# apply(CN2[13930:13935,]@assayData$calls / CN2[13930:13935,]@phenoData@data$rpc, 1, median)
# CN2[13932,]@assayData$copynumber
# # 13932
# 
# plotCopynumberHeatmap(CN2[15380:15400,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[16800:16820,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[18150:18180,], cluster_rows=FALSE)
# # fine 
# 
# plotCopynumberHeatmap(CN2[19500:19530,], cluster_rows=FALSE)
# CN2[19510:19520,]@assayData$copynumber
# apply(CN2[19510:19520,]@assayData$calls / CN2[19510:19520,]@phenoData@data$rpc, 1, median)
# CN2[13932,]@assayData$copynumber
# # 13932
# 
# plotCopynumberHeatmap(CN2[20840:20925,], cluster_rows=FALSE)
# CN2[20840:20870,]@assayData$copynumber
# apply(CN2[20855:20860,]@assayData$calls / CN2[20855:20860,]@phenoData@data$rpc, 1, median)
# CN2[20855,]@assayData$copynumber
# # 20855
# 
# plotCopynumberHeatmap(CN2[22000:22020,], cluster_rows=FALSE)
# CN2[22002:22020,]@assayData$copynumber
# apply(CN2[22001:22020,]@assayData$calls / CN2[22001:22020,]@phenoData@data$rpc, 1, median)
# # 22002, 22003
# 
# plotCopynumberHeatmap(CN2[23070:23090,], cluster_rows=FALSE)
# CN2[23065:23080,]@assayData$copynumber
# apply(CN2[23069:23077,]@assayData$calls / CN2[23069:23077,]@phenoData@data$rpc, 1, median)
# # 23069-23076
# 
# plotCopynumberHeatmap(CN2[24100:24120,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[25000:25030,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[25820:25840,], cluster_rows=FALSE)
# CN2[25824:25825,]@assayData$copynumber
# apply(CN2[25820:25840,]@assayData$calls / CN2[25820:25840,]@phenoData@data$rpc, 1, median)
# # 25824, 25825
# 
# plotCopynumberHeatmap(CN2[26590:26620,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[27180:27210,], cluster_rows=FALSE)
# CN2[27197:27200,]@assayData$copynumber
# apply(CN2[27197:27200,]@assayData$calls / CN2[27197:27200,]@phenoData@data$rpc, 1, median)
# # 27197
# 
# plotCopynumberHeatmap(CN2[28300:28330,], cluster_rows=FALSE)
# # fine
# 
# plotCopynumberHeatmap(CN2[28800:28850,], cluster_rows=FALSE)
# CN2[28800:28850,]@assayData$copynumber
# apply(CN2[28825:28850,]@assayData$calls / CN2[28825:28850,]@phenoData@data$rpc, 1, median)
# # 28826-28850
# 
# plotCopynumberHeatmap(CN2[30360:30490,], cluster_rows=FALSE)
# CN2[30369:30400,]@assayData$copynumber
# apply(CN2[30369:30400,]@assayData$calls / CN2[30369:30400,]@phenoData@data$rpc, 1, mean)
# # 30369-30400
# 
# 
# # table(g2 %in% g1)
# binsToUseInternal(obj)[8800:8820]
# # 1 high readout leads to jump in state, this is not fixed later...
# obj[30500:30550,] -> should not be 1 copy number state
# obj[30500:30550,]@assayData$calls
# obj[30500:30550,]@assayData$segmented
# 
# # second example
# obj[8800:8820, ]
# obj[8800:8820, ]@assayData$copynumber
# plotCopynumber(obj[8800:8820, ])
# obj[8800:8820, ]@assayData$calls
# obj[8800:8820, ]@assayData$probdloss
# obj[8800:8820, ]@assayData$segmented
# 
# a = CN2@assayData$calls / CN2@phenoData@data$rpc
# b = apply(a, 1, median, na.rm=TRUE)
# 
# 6:58700001-58800000  4:190500001-190600000   16:70900001-71000000   17:21200001-21300000             4:1-100000 
# 10.145935               7.955250               5.026639               4.990455               4.854051 
# 6:300001-400000   16:71000001-71100000   17:21300001-21400000   10:47000001-47100000   18:18500001-18600000 
# 4.729752               4.452376               4.392340               4.351512               4.214220 
# 5:134200001-134300000    6:57300001-57400000   10:80900001-81000000   17:44200001-44300000 10:127500001-127600000 
# 4.185474               4.138114               3.816044
# 
# 
# # CN3 = CN[, colnames(CN) %in% dp$name]
# 
# 
# 
# hist(CN@featureData@data$mappability, breaks=100)
# hist(CN@featureData@data$gc, breaks=100)
# table( & CN@featureData@data$mappability > 0)
# table(CN@featureData@data$gc > 56)
# 
# obj = CN[,1]
# 
# valid = binsToUseInternal(obj) 
# valid = binsToUseInternal(obj) & obj@featureData@data$mappability > 70 & obj@featureData@data$gc < 56 & obj@featureData@data$gc > 34
# 
# 
# reads = Biobase::assayDataElement(obj, "calls")[valid, 1, drop=TRUE]
# segs = Biobase::assayDataElement(obj, "segmented")[valid, 1, drop=TRUE]
# gc = Biobase::fData(obj[valid, 1])[["gc"]]
# map = Biobase::fData(obj[valid, 1])[["mappability"]]
# stopifnot(length(gc) == length(reads))
# 
# # log transform seg level for better stability
# dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
#                  reads_normalized=as.vector(reads)/as.vector(segs),
#                  gc=gc, map=map)
# 
# ## Build models
# if(all(map == 100)){
#   formula.gam = formula(reads ~ segs + s(gc))
#   newdata=data.frame(segs = dat$segs, gc=dat$gc)
# }else{
#   formula.gam = formula(reads ~ segs + s(gc, map))  
#   newdata=data.frame(segs = dat$segs, gc=dat$gc, map=dat$map)
# }
# dat_model = dat %>% dplyr::filter(segs > 0)
# model <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
# terms = predict(model, newdata=newdata, type="terms")
# 
# plot(model)
# summary(model)
# 
# 
# # check on diploid cell line
# # annotation = readRDS(file.path("~/scAbsolute/data/blacklisting/extendedBlacklisting-crosstechnologies-GRCh37.RDS"))
# 
# ggplot(data=df %>% dplyr::filter(UID=="UID-JBL-NA12878")) + geom_point(aes(x=cellcycle, y=ploidy))
# 
# dp = df %>% dplyr::filter(UID=="UID-JBL-NA12878") %>% dplyr::filter(cellcycle %in% c("G1"), abs(ploidy - round(ploidy)) < 0.05)
# ggplot(data=dp) + geom_point(aes(x=cellcycle, y=hmm.alpha))
# 
# CN2 = CN[, colnames(CN) %in% dp$name]
# plotCopynumberHeatmap(CN2, cluster_rows=FALSE)
# plotCopynumberHeatmap(CN2[8800:8820,], cluster_rows = FALSE, row_split = as.factor(dp$cellcycle))
# plotCopynumberHeatmap(CN2[30000:30970,], cluster_rows = FALSE, row_split = as.factor(dp$cellcycle))
# 
# table(apply(CN2@assayData$copynumber, 1, median) != 2)
# which(apply(CN2@assayData$copynumber, 1, median) != 2)
# 
# cor = estimate_gc_correction(CN2[, 1])
# 
# # length(cor)
# 
# apply(cor[8830:8840, ], 1, median)
# apply(cor[8820:8830, ], 1, median)
# plotCopynumber(CN2[8800:8850,8])
# plotCopynumber(CN2[30300:30570, 5], correction = TRUE, addSegmentation = TRUE)
# 
# table(CN2@assayData$copynumber[30300:30570, 5])
# 
# 
# plotCopynumberHeatmap(CN2[30400:30500,], cluster_rows = FALSE, row_split = as.factor(dp$cellcycle))
# a = apply(CN2@assayData$calls[30400:30500,], 1, mean, na.rm=TRUE)
# b = apply(CN2@assayData$calls[30500:30580,], 1, mean, na.rm=TRUE)
# 
# m1 = CN2@featureData@data$mappability[30400:30500]
# m2 = CN2@featureData@data$mappability[30500:30580]
# g1 = CN2@featureData@data$gc[30400:30500]
# g2 = CN2@featureData@data$gc[30500:30580]
# set.seed(42)
# p1 = hist(m1, breaks=100, col="blue")
# p2 = hist(m2, breaks=100, col="red")
# plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,100))  # first histogram
# plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,100), add=T)  # second
# dev.off()
# set.seed(42)
# p3 = hist(g1, breaks=100, col="blue")
# p4 = hist(g2, breaks=100, col="red")
# plot( p3, col=rgb(0,0,1,1/4), xlim=c(35,50))  # first histogram
# plot( p4, col=rgb(1,0,0,1/4), xlim=c(35,50), add=T)  # second
# 
# # dev.off()
# # a = apply(, 1, mean, na.rm=TRUE)
# # b = apply(CN2@featureData@data$mappability[30500:30580], 1, mean, na.rm=TRUE)
# # hist(a, breaks=100)
