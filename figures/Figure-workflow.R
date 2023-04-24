# Copyright 2022, Michael Schneider, All rights reserved.
## Fig. 1 - workflow


# cells selected:
# s-phase vs normal
# normal: "UID-JBL-CIOV1_SLX-18430_000036_SINCEL-31-Plate-95-D5-GGACTCCT-GCGTAAGA"
# s phase: "UID-JBL-CIOV1_SLX-18430_000053_SINCEL-33-Plate-100-E7-CTCTCTAC-CTATTAAG"
# G2 cell (DLP): "UID-DLP-SA1044_SLX-A96139A_000296_R10-C40"
#plotCopynumber(CN[, colnames(CN)[grepl("DLP-SA1044_SLX-A96139A_000296", colnames(CN))]], correction = TRUE)
#plotCopynumber(CN[, colnames(CN)[grepl("DLP-SA1044_SLX-A96139A_000604", colnames(CN))]], correction = TRUE)
# G1 cell (ACT, tetraploid): "UID-NNA-mb453_SLX-NAVINACT_000714_MDAMB453-Parental-3-D-C18-S2706-SRR13904692"
#plotCopynumber(CN[, colnames(CN)[grepl("NNA-mb453_SLX-NAVINACT_000714", colnames(CN))]], correction = TRUE)

require(QDNAseq)
library(ggplot2, quietly=TRUE, warn.conflicts = FALSE)
library(ggpubr, quietly=TRUE, warn.conflicts = FALSE)
library(png, quietly=TRUE, warn.conflicts = FALSE)
source("~/scAbsolute/R/core.R")
source("~/scAbsolute/R/visualization.R")
source("~/scAbsolute/R/scAbsolute.R")
source("~/scAbsolute/R/scSegment.R")
source("~/scAbsolute/R/mean-variance.R")
BASEDIR="~/scAbsolute/"

#cell_path = list.files("~/Data/Figures/scAbsolute/rawdata/", pattern = ".bam", full.names = TRUE)
#readCounts1 = readData(cell_path[1:4], binSize=500, genome="hg19", extendedBlacklisting = TRUE, filterChromosomes = c("MT"))
#readCounts2 = readData(cell_path[5:6], binSize=500, genome="hg19", extendedBlacklisting = TRUE, filterChromosomes = c("MT"))
#object = readRDS("~/Data/Figures/scAbsolute/plotExamples-scAbsolute.RDS")
#
#
### Show characteristic cell line with G1 / S phase:
#rpc1 = object[, 3]@phenoData@data$rpc
#p1 = plotCounts(object[, 3], main="", correction = FALSE, ylim = c(0, 2000),
#      #chromosome_break_label = c(as.character(1:22), "X", "Y"), 
#      addSegmentation=TRUE) +
#  scale_y_continuous(sec.axis = sec_axis( trans=~./rpc1, name="absolute copy number", breaks=seq(0, 5, by = 1))) +
#  ylab("read counts")
#p1
#
#p2 = plotCounts(object[, 4], main="", correction = FALSE, ylim = c(0, 2000),
#                #chromosome_break_label = c(as.character(1:22), "X", "Y"),
#                addSegmentation = TRUE)
#p2
#
#segmentIndex = startsWith(rownames(object), "7:") & binsToUseInternal(object) & object[, 3]@assayData$calls > 100
#correctionFactor = estimate_gc_correction(object[, 3])[segmentIndex]
#calls = object[, 3]@assayData$calls[segmentIndex]
#repTime = object[segmentIndex, 3]@featureData@data$replicationTiming
#p1b = ggplot(data = dplyr::tibble(y = calls, x=repTime, cor = correctionFactor)) + 
#  geom_point(aes(y=y, x = x, color=cor)) + theme_pubclean() +
#  geom_smooth(aes(x=x, y=y), method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
#  theme_pubclean() +
#  theme(legend.position = "right") +
#  scale_color_viridis_c("GC\ncorrection") +
#  xlab("replication time") + ylab("observed counts")
#p1b 
#idx = base::sort(repTime, index.return=TRUE)$ix
#cor.corrected1 = trend::partial.cor.trend.test(calls[idx], correctionFactor[idx], method = "spearman")[["estimate"]]
#cor.corrected1
#
#segmentIndex = startsWith(rownames(object), "7:") & binsToUseInternal(object) & object[, 4]@assayData$calls > 100
#correctionFactor = estimate_gc_correction(object[, 4])[segmentIndex]
#calls = object[, 4]@assayData$calls[segmentIndex]
#repTime = object[segmentIndex, 4]@featureData@data$replicationTiming
#idx = base::sort(repTime, index.return=TRUE)$ix
#cor.corrected2 = trend::partial.cor.trend.test(calls[idx], correctionFactor[idx], method = "spearman")[["estimate"]]
#cor.corrected2
#p2b = ggplot(data = dplyr::tibble(y = calls, x=repTime, cor = correctionFactor)) + 
#  geom_point(aes(y=y, x = x, color=cor)) + theme_pubclean() +
#  geom_smooth(aes(x=x, y=y), method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
#  theme_pubclean() +
#  theme(legend.position = "right") +
#  scale_color_viridis_c("GC\ncorrection") +
#  xlab("replication time") + ylab("observed counts")
#p2b 
object = readRDS("~/Data/Figures/scAbsolute/Example-SA1044-matched-rpc.RDS")
# cells with similar rpc
#[1] "UID-DLP-SA1044_SLX-A96139A_000007_R03-C17" - G1
#[2] "UID-DLP-SA1044_SLX-A96139A_000057_R04-C47" - G2

## Margin plots
index = startsWith(rownames(object), "20") | startsWith(rownames(object), "21") |
  startsWith(rownames(object), "22") |  startsWith(rownames(object), "23") |  
  startsWith(rownames(object), "X")
p1 = plotCopynumber(object[, 1], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3,
               showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange",
               chromosome_break_label = NULL, correction = TRUE) + 
  coord_cartesian(ylim=c(0, 8), clip="off") +
  theme_pubclean(base_size=20) +
  theme(plot.margin = margin(1.0, 0, 0, 0, "cm")) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 1.0, name = "correct absolute copy number", breaks = seq(0,8,2), labels=sprintf("%02s", seq(0, 8, 2))) #, color="black")
  )
  geom_text(x=1000, y=9.3, size=8, label="hypotriploid (G1)") #c(as.character(20:22), "X"), correction=TRUE)
p1

p2 = plotCopynumber(object[, 2], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3,
               showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange",
               chromosome_break_label = NULL, correction = TRUE) + #c(as.character(20:22), "X"), correction=TRUE)
  coord_cartesian(ylim=c(0, 8), clip="off") +
  theme_pubclean(base_size=20) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 2.0, name = "correct absolute copy number", breaks = seq(0,16,4), labels=sprintf("%02s", seq(0, 16, 4))) #, color="black")
  )
  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"), axis.ticks.x=element_blank()) +
  geom_text(x=1200, y=9.3, size=8, label="hypohexaploid (G2)") #c(as.character(20:22), "X"), correction=TRUE)
p2



object@phenoData@data$name
object@phenoData@data$rpc

# fit-and-scale
# Show only marginal plot for x axis
obj = object[index,1] 
n_bins = sum(index) #dim(object)[[1]]
x = seq(1, n_bins)
rpc = obj@phenoData@data$rpc
calls = Biobase::assayDataElement(obj, "probdloss")
segments = Biobase::assayDataElement(obj, "segmented")
copies = Biobase::assayDataElement(obj, "copynumber")
d = data.frame(x=x, y=as.numeric(calls), y2=as.numeric(segments)*rpc)
# x axis
chromosome_break_label = c("", "", "", "") #c(as.character(1:22), "X", "Y")
location = rownames(obj)
seqName = unlist(lapply(strsplit(location, ":"), `[`, 1))
chrOrder = base::rle(seqName)[["values"]]
chromosome <- factor(seqName, levels=chrOrder)
rl = S4Vectors::Rle(chromosome)
if(length(rl@values) == 1){
  chromosome_break_position = c(rl@lengths / 2)
  chromosome_boundaries = c()
}else{
  rl_cumsum = cumsum(rl@lengths)
  chromosome_boundaries = rl_cumsum[1:(length(rl_cumsum)-1)]
  chromosome_break_position = c(0, cumsum(rl@lengths)[1:(length(rl_cumsum)-1)]) + ceiling(rl@lengths/2)
}
po3 = ggplot(data=d, aes(x=x, y=y2)) + 
  geom_point(color="#E69F00", size=0.7, alpha=1.0) + ylab("read counts") +
  geom_point(aes(x=x, y=y), size=0.2, alpha=0.45) +
  scale_y_continuous(breaks=scales::pretty_breaks()) + 
  scale_x_continuous("chromosome", breaks=chromosome_break_position, labels=chromosome_break_label, expand=c(0,0)) +
  geom_vline(xintercept = chromosome_boundaries, linetype = "dashed", size=0.5, alpha=1.0) +
  ylim(c(0, 750)) +
  theme_cowplot() + theme_pubclean(base_size=20)

library(ggExtra)
# ggplot(data = data.frame(aes))
pp1 <- ggMarginal(po3 + theme(
                        text = element_text(size = 30),
                        axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.line.x=element_blank(),
                        axis.ticks.x=element_blank()),
                 y="y2", margins = 'y',type="density", color="orange", size=5, 
                 yparams = list(adjust = 1/4, fill = "orange"))
pp1


pp2 = ggplot(data = d) +  coord_cartesian(xlim=c(0, 600)) + geom_density(aes(x=y2), adjust=1/4) +
  expand_limits(x = c(0, 600)) + 
  geom_vline(xintercept = c(rpc, 2*rpc, 3*rpc, 4*rpc), color="red", linetype="dashed") +
  xlab("read counts") + 
  theme_pubclean()
pp2



## fit-and-scale
## Show only marginal plot for x axis
#obj = object[,5] 
#n_bins = dim(object)[[1]]
#x = seq(1, n_bins)
#rpc = obj@phenoData@data$rpc
#calls = Biobase::assayDataElement(obj, "probdloss")
#segments = Biobase::assayDataElement(obj, "segmented")
#copies = Biobase::assayDataElement(obj, "copynumber")
#d = data.frame(x=x, y=as.numeric(calls), y2=as.numeric(segments)*rpc)
## x axis annotation
#chromosome_break_label = c(as.character(1:22), "X", "Y")
#location = rownames(obj)
#seqName = unlist(lapply(strsplit(location, ":"), `[`, 1))
#chrOrder = base::rle(seqName)[["values"]]
#chromosome <- factor(seqName, levels=chrOrder)
#rl = S4Vectors::Rle(chromosome)
#if(length(rl@values) == 1){
#  chromosome_break_position = c(rl@lengths / 2)
#  chromosome_boundaries = c()
#}else{
#  rl_cumsum = cumsum(rl@lengths)
#  chromosome_boundaries = rl_cumsum[1:(length(rl_cumsum)-1)]
#  chromosome_break_position = c(0, cumsum(rl@lengths)[1:(length(rl_cumsum)-1)]) + ceiling(rl@lengths/2)
#}
#po4 = ggplot(data=d, aes(x=x, y=y2)) + 
#  geom_point(color="#E69F00", size=0.1, alpha=1.0) + ylab("read counts") +
#  geom_point(aes(x=x, y=y), size=0.1, alpha=0.05) +
#  scale_x_continuous("chromosome", breaks=chromosome_break_position, labels=chromosome_break_label, expand=c(0,0)) +
#  geom_vline(xintercept = chromosome_boundaries, linetype = "dashed", size=0.5, alpha=1.0) +
#  scale_y_continuous(breaks=scales::pretty_breaks()) + 
#  ylim(c(0, 350)) +
#  theme_cowplot()
#
#library(ggExtra)
## ggplot(data = data.frame(aes))
#pp4 <- ggMarginal(po4 + theme(axis.title.x=element_blank(),
#                           axis.text.x=element_blank(),
#                           axis.line.x=element_blank(),
#                           axis.ticks.x=element_blank()),
#                 y="y2", margins = 'y',type="density", color="orange", size=5, 
#                 yparams = list(binwidth = 1, fill = "orange"))
#pp4


filePath = "~/Data/Figures/scAbsolute/rawdata/UID-DLP-SA1044_SLX-A96139A_000007_R03-C17.bam"
test = readPosition(object[,1], filePath)
dens = test@protocolData@data

dty = dens$mean_overlap[dens$copy >= 1 & dens$n_elem >= 10]
#dty = dens$muf125[dens$copy >= 1  & dens$n_elem >= 10 & !is.na(dens$muf125)]
dtx = dens$copy[dens$copy >= 1 & dens$n_elem >= 10]
n_elem = dens$n_elem[dens$copy >= 1 & dens$n_elem >= 10]


filePath = "~/Data/Figures/scAbsolute/rawdata/UID-DLP-SA1044_SLX-A96139A_000057_R04-C47.bam"
test = readPosition(object[,2], filePath)
dens = test@protocolData@data

dty2 = dens$mean_overlap[dens$copy >= 1 & dens$n_elem >= 10]
#dty = dens$muf125[dens$copy >= 1  & dens$n_elem >= 10 & !is.na(dens$muf125)]
dtx2 = dens$copy[dens$copy >= 1 & dens$n_elem >= 10]
n_elem2 = dens$n_elem[dens$copy >= 1 & dens$n_elem >= 10]

# check prediction for copy=2
require(robustbase)
lmrob_control = robustbase::lmrob.control(setting="KS2014",
                                          max.it = 10000, maxit.scale = 10000,
                                          scale.tol=1e-5, solve.tol = 1e-5,
                                          subsampling = "nonsingular")
model = withCallingHandlers( robustbase::lmrob(prediction ~ I(copy), weights=weight,
                                                                       data = data.frame(copy=dtx,
                                                                                         prediction=dty,
                                                                                         weight=n_elem),
                                                                       control=lmrob_control))
                    prediction.value = predict(model, newdata=dplyr::tibble(copy=2.0))
prediction = predict(model, newdata=dplyr::tibble(copy = 2.0))

# check prediction for copy=2
require(robustbase)
lmrob_control = robustbase::lmrob.control(setting="KS2014",
                                          max.it = 10000, maxit.scale = 10000,
                                          scale.tol=1e-5, solve.tol = 1e-5,
                                          subsampling = "nonsingular")
model = withCallingHandlers( robustbase::lmrob(prediction ~ I(copy), weights=weight,
                                                                       data = data.frame(copy=dtx2,
                                                                                         prediction=dty2,
                                                                                         weight=n_elem2),
                                                                       control=lmrob_control))
                    prediction.value = predict(model, newdata=dplyr::tibble(copy=2.0))
prediction2 = predict(model, newdata=dplyr::tibble(copy = 2.0))

library(ggpubr)
pi1 = ggplot(data=dplyr::tibble(x=c(dtx, dtx2), y=c(dty, dty2), n=c(n_elem, n_elem2),
                                cell=c(rep("hypotriploid (G1)", length(dtx)), rep("hypohexaploid (G2)", length(dtx2))))) +
  geom_point(aes(x=x, y=y, size=n, color=as.factor(cell))) +
  geom_smooth(aes(x=x, y=y, color=as.factor(cell)), method=robustbase::lmrob, formula=y ~ poly(x,2)) + 
  theme_pubclean() + xlab("copy number") + ylab("read density") +
#  scale_y_continuous(breaks=c(0, 40, 80)) +
#  scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6)) +
  geom_point(aes(x=2.0, y=prediction), color="black", shape=5) +
  geom_point(aes(x=2.0, y=prediction2), color="black", shape=5) +
  #geom_text(x=2.5, y=1.57, angle=26, size=10, label="hypohexaploid (G2)") +
  #geom_text(x=4.2, y=1.37, angle=16, size=10, label="hypotriploid (G1)") +
  theme_pubclean(base_size=25) +
  guides(size="none", colour = guide_legend(title="", nrow=2, override.aes = list(alpha = 1, fill=NA))) +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 24, face = 2)) +
  theme(legend.position = "top")#c(0.2, 0.85))
pi1


## Model images
im1 <- ggdraw() +
  draw_image("~/scAbsolute/figures/Panel_graphical-model.pdf") +
  theme_nothing()

im2 <- ggdraw() +
  draw_image("~/scAbsolute/figures/Panel_read-density.pdf", scale=1000) + 
  theme_nothing()

# shared y axis with patchwork
#font_size=16
p1b_legend = get_legend(p2b)
result <- ((p1b + theme(legend.position = "none") + rremove("xlab") + rremove("ylab")) / (p2b + theme(legend.position = "none") + rremove("ylab")))
gt <- patchwork::patchworkGrob(result)

joint = ggpubr::ggarrange(p1+rremove('xlab')+rremove("ylab")+ggtitle("hypotriploid (G1)")+#rremove("y.title")
                          theme(axis.ticks.x=element_blank(),
                                                                        plot.margin = unit(c(5.0, 15, 5.0, 4.0), "pt")),
                          p2+ggtitle("hypohexaploid (G2)")+rremove("ylab")+ #rremove("y.title") +
                            theme(axis.ticks.x=element_blank(), plot.margin = unit(c(5.0, 10, 5.0, 4.0), "pt")), nrow=2, heights = c(1,1))
joint

joint = annotate_figure(joint, left = textGrob("Initial absolute copy number", rot = 90, vjust = 0.1, gp = gpar(cex = 2.0)),
                right = textGrob("Correct absolute copy number", rot = 270, vjust = 0.5, gp = gpar(cex = 2.0)))
joint

#joint = annotate_figure(joint, left = textGrob("absolute copy number", rot = 90, vjust = 1), fig.lab.size=20) #)font_size)
                #bottom = textGrob("Common x-axis", gp = gpar(cex = 1.3)))
#theme(axis.text.x=element_text(size=15))

Fig_workflow = ggpubr::ggarrange(
  ggpubr::ggarrange(pp1,#+theme(axis.title.y = element_text(size = 20)),
                    pp2 +
                      scale_x_continuous(sec.axis = sec_axis(~ . * (1/rpc), name = "absolute copy number")) +
                      geom_segment(aes(x=rpc,   y=0.011, xend=2*rpc, yend=0.011), arrow=arrow(length = unit(1, "mm"))) +
                      geom_segment(aes(x=2*rpc, y=0.011, xend=rpc, yend=0.011), arrow=arrow(length = unit(1, "mm"))) +
                      geom_segment(aes(x=2*rpc, y=0.010, xend=3*rpc, yend=0.010), arrow=arrow(length = unit(1, "mm"))) +
                      geom_segment(aes(x=3*rpc, y=0.010, xend=2*rpc, yend=0.010), arrow=arrow(length = unit(1, "mm"))) +
                      geom_segment(aes(x=3*rpc, y=0.009, xend=4*rpc, yend=0.009), arrow=arrow(length = unit(1, "mm"))) +
                      geom_segment(aes(x=4*rpc, y=0.009, xend=3*rpc, yend=0.009), arrow=arrow(length = unit(1, "mm"))) +
                      geom_text(aes(x = rpc+0.5*rpc, y = 0.009), label = expression(rho), size=9, parse=TRUE, colour = "red") +
                      geom_text(aes(x = 2*rpc+0.5*rpc, y = 0.008), label = expression(rho), size=10, parse=TRUE, colour = "red") +
                      geom_text(aes(x = 3*rpc+0.5*rpc, y = 0.007), label = expression(rho), size=11, parse=TRUE, colour = "red")
                    ,
                    pi1+xlab("absolute copy number"), ncol=3, widths = c(3, 2, 2)),
  ggpubr::ggarrange(joint, nrow=1, ncol=2, widths = c(2,1)), nrow=2, labels=c("", ""), heights = c(2,3)
)

Fig_workflow

Fig_workflow_A = pp1
# not a ggplot object, so manipulate directly above
Fig_workflow_A
ggexport(Fig_workflow_A, filename = "~/scAbsolute/figures/Fig_workflow_A.pdf", width = 8, height = 5)

Fig_workflow_B = pp2 +
                    scale_x_continuous(sec.axis = sec_axis(~ . * (1/rpc), name = "absolute copy number")) +
                    geom_segment(aes(x=rpc+0.05,   y=0.011, xend=2*rpc-0.05, yend=0.011), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_segment(aes(x=2*rpc+0.05, y=0.011, xend=rpc-0.05, yend=0.011), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_segment(aes(x=2*rpc+0.05, y=0.010, xend=3*rpc-0.05, yend=0.010), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_segment(aes(x=3*rpc+0.05, y=0.010, xend=2*rpc-0.05, yend=0.010), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_segment(aes(x=3*rpc+0.05, y=0.009, xend=4*rpc-0.05, yend=0.009), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_segment(aes(x=4*rpc+0.05, y=0.009, xend=3*rpc-0.05, yend=0.009), arrow=arrow(length = unit(6, "mm")), size=1) +
                    geom_text(aes(x = rpc+0.5*rpc, y = 0.010), label = expression(rho), size=10, parse=TRUE, colour = "black") +
                    geom_text(aes(x = 2*rpc+0.5*rpc, y = 0.009), label = expression(rho), size=10, parse=TRUE, colour = "black") +
                    geom_text(aes(x = 3*rpc+0.5*rpc, y = 0.008), label = expression(rho), size=10, parse=TRUE, colour = "black") +
  coord_flip() + theme_pubclean(base_size=30)
Fig_workflow_B

ggexport(Fig_workflow_B, filename = "~/scAbsolute/figures/Fig_workflow_B.pdf", width = 8, height = 5)

Fig_workflow_C = ggarrange(joint + theme_pubclean(base_size=30),
                           pi1+xlab("absolute copy number"),
  #ggpubr::ggarrange(joint, nrow=1, ncol=2, widths = c(2,1)),
  nrow=1, ncol=2, labels=c("", ""), widths = c(5,3))
Fig_workflow_C 
ggexport(Fig_workflow_C, filename = "~/scAbsolute/figures/Fig_workflow_C.pdf", width = 12, height = 5)

Fig_workflow_C_demo = ggpubr::ggarrange(p1+rremove('xlab')+rremove("ylab")+
                          theme(axis.ticks.x=element_blank(),
                                plot.margin = unit(c(5.0, 15, 5.0, 4.0), "pt")),
                          NULL,
                          p2+rremove("ylab")+ #rremove("y.title") +
                            theme(axis.ticks.x=element_blank(),
                                  plot.margin = unit(c(5.0, 10, 5.0, 4.0), "pt")),
                          nrow=3, heights = c(1, 0.15, 1))

Fig_workflow_C_demo
ggexport(Fig_workflow_C_demo, filename = "~/scAbsolute/figures/Fig_workflow_C_demo.pdf", width = 12, height = 5)

#  ggpubr::ggarrange(p2,
#                    gridExtra::grid.arrange(gt, left = "read counts", clip="on", padding = unit(0.0, "line")), p1b_legend, ncol=3, widths = c(30, 10, 3)),
#  ggpubr::ggarrange(im1, pp4, ncol=2, widths=c(1,3), labels=c("C", "D")),
#  ggpubr::ggarrange(pp3,
#                    ggpubr::ggarrange(im2, pi1, ncol=1, heights = c(1, 1.5)),
#                    ncol=2, widths = c(2, 1)),
# im2 + theme(plot.margin = margin(10, 5, 5, 10, "pt"))
  # ggpubr::ggarrange(pp3, im2 + theme(plot.margin = margin(10, 5, 5, 10, "pt")), pi1, ncol=3, widths = c(2, 1, 1)),
  # ggpubr::ggarrange(p5 + theme(axis.title.x=element_blank(),
  #                             axis.text.x=element_blank(),
  #                             axis.ticks.x=element_blank()), p6, nrow=2, labels=c("E", "F")), 
#  nrow=4, labels=c("A", "B", "", "E")#, vjust=-0.5
#)
#Fig_workflow

#ggexport(Fig_workflow, filename = "~/scAbsolute/figures/Fig_workflow.pdf", width = 8, height = 5)
