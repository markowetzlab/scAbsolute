# Copyright 2022, Michael Schneider, All rights reserved.
## Predict ploidy 2 - demonstrate ploidy problem with examples


## Data and library ====
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
source("~/scAbsolute/R/visualization.R")

binSizeValue=500
files = c(
  paste0("~/Data/project1/30-scale/", binSizeValue, "/predict/", 
         paste0(c(
           "UID-DLP-SA1044",
           "UID-DLP-SA928"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files, 
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", "")) %>%
  dplyr::mutate(cycling_activity = cellcycle.kendall.repTime.weighted.median.cor.corrected)

# load cellcycle info
df_cellcycle_DLP = readr::read_csv("~/mean-variance-model/data/shahlab/cellcycle-info-DLP.all.csv", col_types = "cccccd")  %>% dplyr::select(UID, SLX, filename, cellcycle, quality)
df_cellcycle_JBL = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv", col_types = "ccccdd") %>% dplyr::select(UID, SLX, filename, cellcycle)
df_meta = dplyr::bind_rows(df_cellcycle_DLP, df_cellcycle_JBL) %>%
  dplyr::mutate(name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

df1 = pd %>% dplyr::filter(SLX=="SLX-00000" | SLX %in% c("SLX-18430","SLX-18431", "SLX-18432", "SLX-20765") |
                             SLX %in% c("SLX-A73046B", "SLX-A73056B", "SLX-A90554C", "SLX-A90600C", "SLX-A95660A", "SLX-A96146A", "SLX-A96172B", "SLX-A96211C") |
                             SLX %in% c("SLX-A73044B", "SLX-A73047D", "SLX-A96199B", "SLX-A96226B") |
                             base::startsWith(SLX, "SLX-A") |
                             SLX %in% c("SLX-A73044A", "SLX-A90553C", "SLX-A96139A", "SLX-A96213A", "SLX-A90554A", "SLX-A90554B", "SLX-NAVINACT")) %>%
  dplyr::left_join(df_meta, by="name", suffix=c("", ".y"))  %>%
  dplyr::mutate(cellcycle = case_when(UID %in% c("UID-NNA-mb453", "UID-10X-Fibroblast-cell", "UID-10X-Fibroblast-nuclei", "UID-10X-Spikein-1pc") ~ "G1", 
                                      TRUE ~ cellcycle)) %>%
  dplyr::mutate(technology = dplyr::recode(as.factor(substr(UID, 5, 7)),
                                           "DLP" = "DLP",
                                           "NNA" = "ACT",
                                           "10X" = "10X",
                                           "JBL" = "JBL")) %>%
  dplyr::mutate(UID2 = dplyr::recode(UID,
                                     "UID-10X-Andor-2020-50941-SNU-638" = "SNU-638", 
                                     "UID-10X-Andor-2020-48959-MKN-45" = "MKN-45",
                                     "UID-10X-Andor-2020-49599-SNU-16" = "SNU-16",
                                     "UID-10X-Andor-2020-49600-SNU-668" = "SNU-668",
                                     "UID-10X-Andor-2020-49606-NCI-N87" = "NCI-N87" ,
                                     "UID-10X-Andor-2020-50941-SNU-638" = "SNU-638",
                                     "UID-10X-Andor-2020-59068-KATOIII" = "KATOIII" ,
                                     "UID-NNA-TN1" = "TN1",
                                     "UID-NNA-TN2" = "TN2",
                                     "UID-NNA-TN3" = "TN3",
                                     "UID-NNA-TN4" = "TN4",
                                     "UID-NNA-TN5" = "TN5",
                                     "UID-NNA-TN8" = "TN8",
                                     "UID-JBL-CIOV1" = "CIOV1",
                                     "UID-JBL-NA12878" = "NA12878",
                                     "UID-JBL-OVCAR3" = "OVCAR3",
                                     "UID-JBL-PEO1" = "PEO1",
                                     "UID-DLP-SA928" = "SA928",
                                     "UID-DLP-SA1044" = "SA1044",
                                     "UID-DLP-SA1089" = "SA1089",
                                     "UID-DLP-SA1090" = "SA1090",
                                     "UID-DLP-SA1135" = "SA1135",
                                     "UID-DLP-SA501X2XB00096" = "SA501",
                                     "UID-NNA-BT20" = "BT20",
                                     "UID-NNA-mb157" =  "mb157" ,
                                     "UID-NNA-mb453" =  "mb453",
                                     "UID-NNA-MDAMB231-popp31" =  "MDAMB231-popp31",
                                     "UID-NNA-MDAMB231c28" =  "MDAMB231c28",
                                     "UID-NNA-MDAMB231c8" =  "MDAMB231c8",
                                     "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))

# load full data
CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
valid_cells = colnames(CN)[colnames(CN) %in% df1$name]
pf = Biobase::protocolData(CN)@data %>% dplyr::filter(name %in% valid_cells)

## Analyse full data
colnames(pf) = paste0("pos.", colnames(pf))
pf1 = dplyr::inner_join(df1, pf, by=c("name"="pos.name"))

set.seed(2022)
v1 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle == "G1", rpc > 100) %>% dplyr::pull(name)
v2 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle == "G2", rpc > 100) %>% dplyr::pull(name)

# select example for Workflow figure
vv1 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle == "G1", rpc > 120, rpc < 125) %>% dplyr::pull(name)
vv2 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle == "G2", rpc > 120, rpc < 125) %>% dplyr::pull(name)
vv1[1]
vv2[1]


p1 = plotCopynumber(CN[, v1[1]], correction = TRUE, main="", readinfo = FALSE, showUnique = FALSE, showMarker = FALSE, ylim=c(0, 8))
p1
p2 = plotCopynumber(CN[, v2[1]], correction = TRUE, main="", readinfo = FALSE, showUnique = FALSE, showMarker = FALSE, ylim=c(0, 8))
p2

w1 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA928", cellcycle == "G1", rpc > 100) %>% dplyr::pull(name)
w2 = pf1 %>% dplyr::filter(UID == "UID-DLP-SA928", cellcycle == "G2", rpc > 100) %>% dplyr::pull(name)

p3 = plotCopynumber(CN[, w1[1]], correction = TRUE, main="", readinfo = FALSE, showUnique = FALSE, showMarker = FALSE, ylim=c(0, 8))
p3
p4 = plotCopynumber(CN[, w2[3]], correction = TRUE, main="", readinfo = FALSE, showUnique = FALSE, showMarker = FALSE, ylim=c(0, 8))
p4


Sup_ploidy_examples = ggpubr::ggarrange(p1 + rremove("x.title") + rremove("y.title"),
                                        p2 + rremove("y.title") + rremove("x.title"),
                                        p3 + rremove("x.title") + rremove("y.title"),
                                        p4 + rremove("y.title"),
                               labels=c("A", "", "B", ""), nrow=4, ncol=1, hjust=1.8)
Sup_ploidy_examples = annotate_figure(Sup_ploidy_examples, left = textGrob("absolute copy number", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
Sup_ploidy_examples

ggpubr::ggexport(Sup_ploidy_examples, filename = "~/scAbsolute/figures/Sup_ploidy_examples.pdf", width=12, height=16)

#Fig_cell_example = plotCopynumber(CN[, v1[1]], correction = TRUE, main="", readinfo = FALSE, showUnique = FALSE, showMarker = FALSE, ylim=c(0, 7)) + theme_pubclean() + theme(text = element_text(size=20))
#ggpubr::ggexport(Fig_cell_example, filename = "~/scAbsolute/figures/Fig_cell_example.pdf", width=24, height=9)
