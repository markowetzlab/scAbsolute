## Analysis cellenone

rm(list=ls())
BASEDIR="~/scAbsolute/"
set.seed(2021)
library(tidyverse, quietly=TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts = FALSE)
library(gghalves, quietly=TRUE, warn.conflicts = FALSE)
library(cowplot, quietly=TRUE, warn.conflicts = FALSE)
library(ggpubr, quietly=TRUE, warn.conflicts = FALSE)
library(ggExtra, quietly=TRUE, warn.conflicts = FALSE)
library(ggbeeswarm, quietly=TRUE, warn.conflicts = FALSE)
library(rtracklayer, quietly=TRUE, warn.conflicts = FALSE)
library(phylogram, quietly=TRUE, warn.conflicts = FALSE)
library(dendextend, quietly=TRUE, warn.conflicts = FALSE)
source("~/scUnique/R/visualize.R")
source("~/scAbsolute/R/core.R")
source("~/scUnique/R/core.R")
source("~/scUnique/R/scUnique.R")
source("~/scAbsolute/R/visualization.R")

binSizeValue=500
files = c(
  paste0("~/Data/project1/30-scale/500/predict/",
         paste0(c(
           "UID-JBL-CIOV1",
           "UID-JBL-NA12878",
           "UID-JBL-PEO1",
           "UID-JBL-OVCAR3"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files,
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))


# load cellcycle info
df_cellcycle_DLP = readr::read_csv("~/mean-variance-model/data/cellcycle-info-DLP.csv", col_types = "ccccdd")
df_cellcycle_JBL = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv", col_types = "ccccdd")
df_dapi = readr::read_tsv("~/mean-variance-model/data/cellenone/DAPI.tsv") %>% setNames(paste0('dapi.', names(.)))
df_trans = readr::read_tsv("~/mean-variance-model/data/cellenone/Transmission.tsv") %>% setNames(paste0('transmission.', names(.)))

df_meta = dplyr::bind_rows(df_cellcycle_DLP, df_cellcycle_JBL) %>%
  dplyr::mutate(name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

df_cellenone = dplyr::bind_cols(df_dapi, df_trans)
df_cellenone$cellenone_ID = gsub("_", "-", df_cellenone$`dapi.Cell ID`)
stopifnot(all(df_cellenone$`dapi.Cell ID` == df_cellenone$`transmission.Cell ID`))

df1 = pd %>% dplyr::filter(SLX=="SLX-00000" | SLX %in% c("SLX-18430","SLX-18431", "SLX-18432", "SLX-20765") |
                          SLX %in% c("SLX-A96139A") | #DLP-SA1044
                          SLX %in% c("SLX-A73044A", "SLX-A90553C", "SLX-A90689B", "SLX-A90689C", "SLX-A90694B") | #DLP-SA928
                          SLX %in% c("SLX-A73046B", "SLX-A73056B", "SLX-A90554C", "SLX-A90600C", "SLX-A95660A", "SLX-A96146A", "SLX-A96172B", "SLX-A96211C") | #DLP-SA906
                          SLX %in% c("SLX-A73044B", "SLX-A73047D", "SLX-A96199B", "SLX-A96226B") | #DLP-SA039
                          SLX %in% c("SLX-NAVINACT")) %>%
  dplyr::left_join(df_meta, by="name", suffix=c("", ".y"))
df1$ID = gsub("\\-(\\w*)\\-(\\w*)$","",df1$celltag, perl=TRUE)

df = dplyr::inner_join(df1, df_cellenone, by=c("ID"="cellenone_ID"))

# load full data
CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
valid_cells = colnames(CN)[colnames(CN) %in% df$name]

## Plot OVCAR with full annotation ====
df_ovcar = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv") %>%
  dplyr::filter(UID == "UID-JBL-OVCAR3") %>%
  dplyr::mutate(status = case_when(min_ploidy == 1.3 ~ "diploid",
                                   min_ploidy == 2.3 ~ "triploid",
                                   min_ploidy == 3.3 ~ "diploid",
                                   min_ploidy == 5.1 ~ "triploid",
                                   TRUE ~ "unknown"),
                name = gsub("\\.bam$", "", filename)) %>%
  dplyr::filter(name %in% colnames(CN))

table(df_ovcar$status, df_ovcar$cellcycle)
df_ovcar$annotation = interaction(df_ovcar$status, df_ovcar$cellcycle)

p1 = plotCopynumberHeatmap(CN[, match(df_ovcar$name,colnames(CN))], row_split = factor(df_ovcar$annotation))
p_heatmap = grid.grabExpr(draw(p1))

p2 = ggplot(data = df1 %>% dplyr::filter(UID == "UID-JBL-OVCAR3") %>% dplyr::inner_join(df_ovcar, by="name")) + 
  geom_quasirandom(aes(x=annotation, y=ploidy)) + theme_pubclean()

p3 = ggplot(data = df1 %>% dplyr::filter(UID == "UID-JBL-OVCAR3") %>% dplyr::inner_join(df_ovcar, by="name")) + 
  geom_density(aes(x=ploidy, color=annotation)) + theme_pubclean() +
  geom_vline(xintercept=c(3.1, 3.3), color="red")

fig_ovcar = ggpubr::ggarrange(p_heatmap, 
                              ggpubr::ggarrange(p2, p3, ncol=2, widths = c(3,1)), ncol=1, heights = c(3, 1))
fig_ovcar

ggpubr::ggexport(fig_ovcar, filename="~/scAbsolute/figures/Check-OVCAR3.pdf")

newObject = readRDS("~/debug-newobject-OVCAR3.RDS")
initial_tree = readRDS("~/debug-tree-OVCAR3.RDS")
dt = phylogram::as.dendrogram(initial_tree)

ploidy = newObject[, labels(dt)]@phenoData@data$ploidy
ha = rowAnnotation(annotation = as.factor(df_ovcar$annotation[match(labels(dt), df_ovcar$name)]), 
                   facs_cellcycle = as.factor(df_ovcar$cellcycle[match(labels(dt), df_ovcar$name)]),
                   facs_ploidy = as.factor(df_ovcar$status[match(labels(dt), df_ovcar$name)]),
                   comp_ploidy = ploidy,
                   comp_ploidy.cens = ifelse(ploidy < 2.9 | ploidy > 3.4, NA, ploidy))
p4 = plotCopynumberHeatmap(newObject[, labels(dt)], cluster_rows = dt, har=ha)
p4_heatmap = grid.grabExpr(draw(p4))

fig_ovcar2 = ggpubr::ggarrange(p4_heatmap)
fig_ovcar2

ggpubr::ggexport(fig_ovcar2, filename="~/scAbsolute/figures/Check-OVCAR3-annotated.pdf", width = 10, height = 6)


## Check cellenone features

q1 = ggplot(data = df) + geom_quasirandom(aes(x=cellcycle, y=dapi.Intensity)) + theme_pubclean()
q2 = ggplot(data = df) + geom_quasirandom(aes(x=cellcycle, y=transmission.Diameter)) + theme_pubclean()
q3 = ggplot(data = df) + geom_quasirandom(aes(x=cellcycle, y=transmission.Elongation)) + theme_pubclean()
q4 = ggplot(data = df) + geom_quasirandom(aes(x=cellcycle, y=transmission.Circularity)) + theme_pubclean()
q5 = ggplot(data = df) + geom_quasirandom(aes(x=cellcycle, y=transmission.Intensity)) + theme_pubclean()

fig_cellenone = ggpubr::ggarrange(q1 + rremove("xlab"),
                  q2 + rremove("xlab"),
                  q3 + rremove("xlab"),
                  q4 + rremove("xlab"),
                  q5, ncol=1)
fig_cellenone

ggpubr::ggexport(fig_cellenone, filename="~/scAbsolute/figures/Check-Cellenone.pdf")

require(caret)
require(kernlab)
require(e1071)

df_predict = df %>% dplyr::filter(cellcycle %in% c("G1", "G2")) %>% dplyr::select(cellcycle, dapi.Intensity, transmission.Diameter, transmission.Elongation, transmission.Circularity, transmission.Intensity) %>%
  dplyr::filter(!is.na(dapi.Intensity))

control <- trainControl(method='cv',
                        number=10)
tunegrid <- expand.grid(.mtry = c(3, 4, 5, 10))

f1 <- caret::train(cellcycle ~ ., data = df_predict,
                   metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
f1$finalModel
