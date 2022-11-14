# Copyright 2022, Michael Schneider, All rights reserved.
## Predict cell cycle correlation - on 10X cell lines
# Data in https://doi.org/10.1093/nargab/lqaa016, see supplementary info


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
           "UID-10X-Andor-2020-48959-MKN-45",
           "UID-10X-Andor-2020-49599-SNU-16",
           "UID-10X-Andor-2020-49600-SNU-668",
           "UID-10X-Andor-2020-49606-NCI-N87",
           "UID-10X-Andor-2020-49607-HGC-27",
           "UID-10X-Andor-2020-50941-SNU-638",
           "UID-10X-Andor-2020-59065-NUGC-4",
           "UID-10X-Andor-2020-59068-KATOIII",
           "UID-10X-Andor-2020-59069-SNU-601"
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
                                     "UID-10X-Andor-2020-49607-HGC-27" = "HGC-27",
                                     "UID-10X-Andor-2020-59065-NUGC-4" = "NUGC-4",
                                     "UID-10X-Andor-2020-59069-SNU-601" = "SNU-601",
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
source("~/scUnique/R/core.R")
df = predict_replicating(df1, cutoff_value = 2) %>% dplyr::group_by(UID2) %>% dplyr::summarise(n = n(), G1=sum(!replicating), S=sum(replicating), read_depth=median(rpc)) %>%
#df = predict_replicating(df1 %>% dplyr::mutate(UID="A", SLX="B"), cutoff_value = 2) %>% dplyr::group_by(UID2) %>% dplyr::summarise(n = n(), G1=sum(!replicating), S=sum(replicating)) %>%
  dplyr::mutate(ratio = S/(G1+S), scRNA_ratio = dplyr::case_when(
    UID2 == "SNU-16" ~ 0.35,
    UID2 == "KATOIII" ~ 0.33,
    UID2 == "HGC-27" ~ 0.51,
    UID2 == "SNU-668" ~ 0.19,
    UID2 == "NUGC-4" ~ 0.28,
    UID2 == "SNU-601" ~ 0.24,
    UID2 == "SNU-638" ~ 0.11,
    UID2 == "MKN-45" ~ 0.33,
    UID2 == "NCI-N87" ~ 0.17,
    TRUE ~ 0.0
  ), doubling_time=dplyr::case_when(
    UID2 == "SNU-16" ~ 27,
    UID2 == "KATOIII" ~ 32,
    UID2 == "HGC-27" ~ 17,
    UID2 == "SNU-668" ~ 74,
    UID2 == "NUGC-4" ~ 36,
    UID2 == "SNU-601" ~ 47,
    UID2 == "SNU-638" ~ 58,
    UID2 == "MKN-45" ~ 24,
    UID2 == "NCI-N87" ~ 47,
    TRUE ~ 0
  ))
df

# SNU-16 and SNU-668 have high passage difference, so are excluded here
p1 = ggplot(data = df %>% dplyr::filter(!(UID2 %in% c("SNU-668", "SNU-16"))) %>%
              dplyr::mutate(same_suspension = UID2 %in% c("KATOIII", "SNU-601", "SNU-638"))) +
  geom_point(aes(x=ratio, y=scRNA_ratio, color=same_suspension)) + #, size=read_depth)) + 
  geom_smooth(aes(x=ratio, y=scRNA_ratio), method='lm', formula= y~x) +
  xlab("% S-phase cells (scDNAseq)") +  ylab("% S-phase cells (scRNAseq)") +
  labs(colour="scDNA & scRNA sequenced from the same suspension") +
  theme_pubclean()

cor.test(df %>% dplyr::filter(!(UID2 %in% c("SNU-668", "SNU-16"))) %>% dplyr::pull(ratio),
         df %>% dplyr::filter(!(UID2 %in% c("SNU-668", "SNU-16"))) %>% dplyr::pull(scRNA_ratio))

p2 = ggplot(data = df) + geom_point(aes(x=ratio, y=doubling_time)) + #, size=read_depth)) + 
  geom_smooth(aes(x=ratio, y=doubling_time), method='lm', formula= y~x) +
  xlab("% S-phase cells (scDNAseq)") + ylab("Doubling time") +
  theme_pubclean()

cor.test(df$ratio, df$doubling_time)

x = predict_replicating(df1 %>% dplyr::filter(UID2 == "SNU-668"))
ggplot(data=x) + geom_quasirandom(aes(x="NA", y=cycling_activity, color=replicating))

#ggplot(data=df1) + geom_quasirandom(aes(x=UID, y=rpc))

#    https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nargab/2/2/10.1093_nargab_lqaa016/2/lqaa016_supplemental_files.pdf?Expires=1659216487&Signature=Owai0oMot2WUq9V1ahsLyYUHplr94qL-AGBBW7gmFGXOMWBNqVsx-HCTG0H3UcenWjDSrz2whPRt8iuqlcBFvY82ad1UCVLeOe0NLtiFAxGwJGlQpb7MI-XOZPiK-UmMXSsxgjusswl3Ywc6IioBOoaGEeqoNZ6YUQE1nLR6-O26H03ZKJQ2TmvuRFD7kZbPOFnwwNL06ViuCFuZRrQnf4U78TGW1-yAs951kQV-1ULg9qkFsRqcAUZgu8Jjle44appaf1eRW40heOwL5I-eMGZZG8SC1fhMftJb76wiqaby9CH7cTTXtoD8f6iTwbFhh-LYRr109HEAjeIhp4IexA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

leg = get_legend(p1)
Sup_cellcycle_cellines = ggpubr::ggarrange(p1, p2,
                               labels=c("A", "B"), nrow=1, ncol=2, common.legend = TRUE, legend.grob = leg) #, hjust=1.8)
#Fig_ploidy_examples = annotate_figure(Fig_ploidy_examples, left = textGrob("absolute copy number", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
Sup_cellcycle_cellines

ggpubr::ggexport(Sup_cellcycle_cellines, filename = "~/scAbsolute/figures/Sup_cellcycle_cellines.pdf", width=12, height=8)


p3 = ggplot(data=predict_replicating(df1)) +
  geom_quasirandom(aes(x=UID2, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=replicating)) +
  #facet_wrap(~UID, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme_pubclean() +
  theme(legend.direction='vertical', legend.position = "right") + 
  ylab("Cycling activity") + xlab("Cell line") + labs(color="Replicating") + 
  theme_pubclean()
  #scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p3

p4 = ggplot(data=predict_replicating(df1 %>% dplyr::mutate(UID = "A", SLX = "B"))) +
  geom_quasirandom(aes(x=UID2, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=replicating)) +
  #facet_wrap(~UID, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme_pubclean() +
  theme(legend.direction='vertical', legend.position = "right") + guides(color="none") +
  ylab("Cycling activity") + xlab("Cell line") + labs(color="Replicating") + 
  theme_pubclean()
  #scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p4


leg = get_legend(p3)
Sup_cellcycle_transfer = ggpubr::ggarrange(p3 + rremove("x.title") + rremove("y.title"), p4 + rremove("y.title"), legend.grob = leg, common.legend = TRUE,
                               labels=c("A", "B"), nrow=2, ncol=1, hjust=1.8, vjust=1.6)
Sup_cellcycle_transfer = annotate_figure(Sup_cellcycle_transfer, left = textGrob("Cycling activity", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
Sup_cellcycle_transfer

ggpubr::ggexport(Sup_cellcycle_transfer, filename = "~/scAbsolute/figures/Sup_cellcycle_transfer.pdf", width=12, height=16)

table(df1$UID)
predict_mix = predict_replicating(df1 %>% dplyr::mutate(UID = "A", SLX = "B"))
predict_unmix = predict_replicating(df1)
stopifnot(all(predict_mix$name == predict_unmix$name))

cor.test(df1$rpc, df1$cycling_activity)

table(predict_mix$replicating, predict_unmix$replicating)
sum(predict_mix$replicating)
sum(predict_unmix$replicating)
