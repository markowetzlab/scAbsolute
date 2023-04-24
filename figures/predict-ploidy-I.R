# Copyright 2022, Michael Schneider, All rights reserved.
## Predict ploidy - demonstrate ability to predict correct ploidy (without G2 details)


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
library(grid, quietly=TRUE, warn.conflicts = FALSE)
source("~/scUnique/R/visualize.R")
source("~/scAbsolute/R/core.R")
source("~/scAbsolute/R/visualization.R")

binSizeValue=500
UID_levels = c("HGC-27", "KATOIII", "MKN-45", "NCI-N87", "NUGC-4",
               "OV2295", "T-47D", "HeLa",
               "BT20","mb157","mb231", "mb453",
               "TN1", "TN2", "TN3", "TN4", "TN5", "TN6", "TN7", "TN8")
#mb231 here is mb231c28
files = c(
  paste0("~/Data/project1/30-scale/", binSizeValue, "/predict/", 
         paste0(c(
           #"UID-DLP-SA039",
           "UID-DLP-SA1044",
           #"UID-DLP-SA928",
           #"UID-DLP-SA906",
           "UID-DLP-SA1087",
           #"UID-DLP-SA1088",
           #"UID-DLP-SA1089",
           "UID-DLP-SA1090",
           #"UID-DLP-SA921",
           "UID-DLP-SA922",
           #"UID-DLP-SA1101",
           #"UID-DLP-SA1135",
#           "UID-DLP-SA501X2XB00096",
#           "UID-DLP-SA501X2XB00097",
           "UID-10X-BREAST-A",
           "UID-10X-BREAST-B",
           "UID-10X-BREAST-C",
           "UID-10X-BREAST-D",
           "UID-10X-BREAST-E",
           "UID-10X-Andor-2020-48959-MKN-45",
           "UID-10X-Andor-2020-49599-SNU-16", 
           #"UID-10X-Andor-2020-49600-SNU-668", 
           "UID-10X-Andor-2020-49606-NCI-N87", 
           #"UID-10X-Andor-2020-50941-SNU-638", 
           "UID-10X-Andor-2020-59068-KATOIII", 
           "UID-10X-Andor-2020-49607-HGC-27",
           "UID-10X-Andor-2020-59065-NUGC-4",
           #"UID-10X-Andor-2020-59069-SNU-601",
           "UID-NNA-BT20", 
           "UID-NNA-mb157", 
           "UID-NNA-mb453",
           "UID-NNA-MDAMB231c28",
           "UID-NNA-TN1",
           "UID-NNA-TN2",
           "UID-NNA-TN3",
           "UID-NNA-TN4",
           "UID-NNA-TN5",
           "UID-NNA-TN6",
           "UID-NNA-TN7",
           "UID-NNA-TN8"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files, 
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

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
                                     "UID-10X-BREAST-A" = "BREAST-A",
                                     "UID-10X-BREAST-B" = "BREAST-B",
                                     "UID-10X-BREAST-C" = "BREAST-C",
                                     "UID-10X-BREAST-D" = "BREAST-D",
                                     "UID-10X-BREAST-E" = "BREAST-E",
                                     "UID-NNA-TN1" = "TN1",
                                     "UID-NNA-TN2" = "TN2",
                                     "UID-NNA-TN3" = "TN3",
                                     "UID-NNA-TN4" = "TN4",
                                     "UID-NNA-TN5" = "TN5",
                                     "UID-NNA-TN6" = "TN6",
                                     "UID-NNA-TN7" = "TN7",
                                     "UID-NNA-TN6-paired" = "TN6-paired",
                                     "UID-NNA-TN7-paired" = "TN7-paired",
                                     "UID-NNA-TN8" = "TN8",
                                     "UID-JBL-CIOV1" = "CIOV1",
                                     "UID-JBL-NA12878" = "NA12878",
                                     "UID-JBL-OVCAR3" = "OVCAR3",
                                     "UID-JBL-PEO1" = "PEO1",
                                     "UID-DLP-SA928" = "SA928",
                                     "UID-DLP-SA921" = "SA921",
                                     #"UID-DLP-SA922" = "SA922",
                                     "UID-DLP-SA922" = "OV2295",
                                     "UID-DLP-SA609" = "SA609",
                                     "UID-DLP-SA906" = "SA906",
                                     #"UID-DLP-SA1044" = "SA1044",
                                     "UID-DLP-SA1044" = "T-47D",
                                     "UID-DLP-SA039" = "SA039",
                                     #"UID-DLP-SA1087" = "SA1087",
                                     "UID-DLP-SA1087" = "HeLa",
                                     "UID-DLP-SA1088" = "SA1088",
                                     "UID-DLP-SA1089" = "SA1089",
                                     "UID-DLP-SA1090" = "SA1090",
                                     "UID-DLP-SA1101" = "SA1101",
                                     "UID-DLP-SA1135" = "SA1135",
                                     "UID-DLP-SA501X2XB00096" = "SA501",
                                     "UID-DLP-SA501X2XB00097" = "SA501",
                                     "UID-NNA-BT20" = "BT20",
                                     "UID-NNA-mb157" =  "mb157" ,
                                     "UID-NNA-mb453" =  "mb453",
                                     "UID-NNA-MDAMB231-popp31" =  "MDAMB231-popp31",
                                     "UID-NNA-MDAMB231c28" =  "mb231",
                                     "UID-NNA-MDAMB231c8" =  "mb231c8",
                                     "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))

# load full data
CN = combineQDNASets(lapply(files, function(x){
  return(selectChromosomes(readRDS(x), exclude=c("X", "Y")))}), dropProtocol=TRUE, reduceMetadata=TRUE)
valid_cells = colnames(CN)[colnames(CN) %in% df1$name]
pf = Biobase::protocolData(CN)@data %>% dplyr::filter(name %in% valid_cells)



## Analyse full data
#colnames(pf) = paste0("pos.", colnames(pf))
#pf1 = dplyr::inner_join(df1, pf, by=c("name"="pos.name"))

read_depth_table = df1 %>% dplyr::group_by(UID, SLX) %>% dplyr::summarise(median_rpc = median(rpc)) %>%
  dplyr::filter(median_rpc < 25)

## Ploidy prediction ====
df_ploidy2 = df1 %>% dplyr::filter(startsWith(UID, "UID-10X") | startsWith(UID, "UID-NNA") | startsWith(UID, "UID-DLP")) %>%
  dplyr::filter(!(UID %in% read_depth_table$UID & SLX %in% read_depth_table$SLX))
                                    #UID %in% c("UID-DLP-SA1089", "UID-DLP-SA1090", "UID-DLP-SA1135", "UID-DLP-SA501X2XB00096")) %>% 
  #dplyr::filter(is.na(cellcycle) | cellcycle %in% c("G1") | UID == "UID-DLP-SA1090")
# DLP-SA1090 has wrong cellcycle annotation in DLP metadata (it's all s phase)

# remove cycling cells
q1a = ggplot(data = predict_replicating(df_ploidy2, cutoff_value = 2)) +
  geom_quasirandom(aes(x=UID2, y=cycling_activity, color=replicating)) +
  facet_wrap(~technology, scales = "free_x") +
  geom_boxplot(aes(x=UID2, y=cycling_activity), alpha=0.0, color="black") + 
  theme_pubclean()
q1b = ggplot(data = predict_replicating(df1 %>% dplyr::filter(startsWith(UID, "UID-10X")), cutoff_value = 2)) +
  geom_quasirandom(aes(x=UID2, y=cycling_activity, color=replicating)) +
  geom_point(aes(x=UID2, y=cycling_median), size=1.0, color="black") + 
  theme_pubclean()
q1 = ggpubr::ggarrange(q1a, q1b)
q1


#df_test = qc_gini_norm(qc_mapd(qc_gini(predict_replicating(qc_rpc(df_ploidy2, cutoff_value = 2))  %>% dplyr::filter(!replicating, !rpc.outlier), cutoff_value = 2.0), cutoff_value = 2.0), cutoff_value = 2.0)
#df_test$outlier = df_test$gini.outlier | df_test$dmapd.outlier | df_test$dgini.outlier
#ggplot(data = df_test) +
#  geom_quasirandom(aes(x=UID, y=gini, color=outlier)) + 
#  theme_pubclean()
#prop.table(table(df_test$outlier, df_test$UID), margin = 2)

ggplot(data = df_ploidy2 %>% dplyr::filter(!startsWith(UID, "UID-NNA-T"))) +
  geom_density(aes(x=gini, color=interaction(technology, UID)))

ggplot(data = df_ploidy2 %>% dplyr::filter(startsWith(UID, "UID-NNA-T"))) +
  geom_point(aes(x=gini, y=gini_normalized, color=interaction(technology, UID)))

p_check = ggplot(data = df_ploidy2 %>% dplyr::filter(!startsWith(UID, "UID-NNA-T"))) +
  geom_quasirandom(aes(x=UID2, y = gini_normalized, color=technology)) + #, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=UID2, y=gini_normalized), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ technology, scales = "free_x", strip.position="bottom", nrow=1) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + 
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside")
p_check



q2 = ggplot(data = df_ploidy2 %>% dplyr::filter(technology %in% c("10X", "ACT"))) +
  geom_quasirandom(aes(x=UID2, y=ploidy, color=technology)) +
  theme_pubclean()
q2
# no error removal
#df_ploidy = qc_gini_norm(qc_error(predict_replicating(df_ploidy2, cutoff_value = 2.0))) %>% dplyr::filter(!replicating, !error.outlier, !dgini.outlier)
df_ploidy = predict_replicating(df_ploidy2, cutoff_value = 1.5) %>% dplyr::filter(!replicating)
#df_ploidy_cleaned = qc_error(qc_mapd(qc_gini(qc_gini_norm(df_ploidy)))) %>%
#        dplyr::filter(!error.outlier, !dgini.outlier, !dmapd.outlier, !gini.outlier)



ggplot(data = df_ploidy) + geom_quasirandom(aes(x=UID2, y = ploidy, color=technology)) +
  theme_pubclean() + theme(axis.text.x = element_text(angle = 90)) 
df_ploidy = qc_ploidy(df_ploidy, range=0.5)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$technology),margin = 2)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$UID2),margin = 2)

## Add meta data information about 10X and ACT cell lines
#source: 
#  10X:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7079336/
#  10X supplementary information: Supplementary Figure 1: Karyotyping and SNP-array analysis confirm scDNA-Seq d
#  manual reading of y-axis on figure F
#  
# ACT:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7079336/#sup1
# Supplementary Table 1 ??? Clinical Tumor Information and Sequencing Metrics
ploidy_table = data.frame(UID2=c("HGC-27", "KATOIII", "MKN-45", "NCI-N87", "NUGC-4", "OV2295", "T-47D", "HeLa", "BT20", "mb157", "mb231", "mb453",
                          "TN1", "TN2", "TN3", "TN4", "TN5", "TN6", "TN7", "TN8"),
                          ploidy=c(3.2, 3.4, 1.9, 1.75, 2.3, 4.2, 2.7, 3.5, 2.70, 2.55, 2.41, 4.17, 3.45, 3.03, 3.44, 3.76, 2.65, 3.17, 3.15, 3.95),
                          #ploidy=c(3.2, 3.4, 1.9, 1.75, 2.3, NA, NA, NA, 2.70, 2.55, 2.41, 4.17, 3.45, 3.03, 3.44, 3.76, 2.65, 3.17, 3.15, 3.95),
                          dtype=c(rep("cellline", 12), rep("patient", 8)),
                          technology=c(rep("10X", 5), rep("DLP", 3), rep("ACT - cell lines", 4), rep("ACT - tumours", 8))) %>%
  dplyr::mutate(ymin = ploidy - 0.2, ymax=ploidy+0.2)


# clean outlier cells
p_ploidy_prediction = ggplot(data = df_ploidy %>% dplyr::filter(rpc > 25, technology != "ACT") %>% 
                               dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"), !startsWith(UID, "UID-DLP-SA1090"), !startsWith(UID, "UID-DLP-SA1101"), !startsWith(UID, "UID-DLP-SA1088"),
                                             !startsWith(UID, "UID-DLP-SA1135"), !startsWith(UID, "UID-DLP-SA928"), !startsWith(UID, "UID-DLP-SA039"), !startsWith(UID, "UID-DLP-SA906"))) +
  geom_quasirandom(aes(y=UID2, x = ploidy, color=technology, shape=ploidy.outlier), groupOnX = FALSE) +
  geom_boxplot(aes(y=UID2, x=ploidy), color="black", outlier.shape = NA, alpha=0.2 ) +
  geom_tile(data=ploidy_table %>% dplyr::filter(technology!="ACT"), aes(y=UID2,x=ploidy, width=1.0,height=0.8),alpha=0.6,fill="grey") +
  geom_point(data=ploidy_table %>% dplyr::filter(technology!="ACT"), aes(y=UID2, x=ploidy), size=2, shape=3, color="red") +
  #geom_tile(data=ploidy_table, aes(x=UID2,y=ploidy, width=0.8,height=0.2),alpha=0.8,fill="grey") +
  facet_wrap(~ technology, scales = "free_y", strip.position="top", nrow=3) +
  theme_pubclean(base_size = 20, flip=TRUE) + theme(legend.position="right") + ylab("sample") + xlab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  guides(shape="none", colour = guide_legend(title="Sequenced with ", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 14, face = 2)) +
  theme(strip.text.y=element_blank(), panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8), labels=c("1", "2", "3", "4", "5", "6", "7", "8")) +
  scale_color_manual(values = c("#56B4E9", "#56B4E9"), labels=c("10X", "DLP"))
#  theme(axis.text.x = element_text(angle = 90))
p_ploidy_prediction

#ggplot(data = qc_error(qc_mapd(qc_gini(qc_gini_norm(df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA928"))))) %>%
#        dplyr::filter(!error.outlier, !dgini.outlier, !dmapd.outlier, !gini.outlier)) +
#  geom_point(aes(x=rpc, y=ploidy))

#table(df_ploidy %>% dplyr::filter(rpc > 25) %>% dplyr::filter(startsWith(UID, "UID-DLP-SA906") | startsWith(UID, "UID-DLP-SA039")) %>% dplyr::pull(UID))
#p_ploidy_prediction_p53 = ggplot(data = df_ploidy %>% dplyr::filter(rpc > 25) %>% dplyr::filter(startsWith(UID, "UID-DLP-SA906") | startsWith(UID, "UID-DLP-SA039"))) +
#  geom_quasirandom(aes(x=SLX, y = ploidy, color=technology, shape=ploidy.outlier)) +
#  geom_boxplot(aes(x=SLX, y=ploidy), color="black", outlier.shape = NA, alpha=0.2 ) +
#  facet_wrap(~ UID2, scales = "free_x", strip.position="bottom", nrow=1) +
#  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
#  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("ploidy") +
#  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
#        strip.background = element_rect(colour="white", fill="white"),
#        strip.placement = "outside") +
#  theme(axis.text.x = element_text(angle = 90))
#p_ploidy_prediction_p53

outlier_SA928 = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA928", ploidy.outlier) %>% dplyr::pull(name)
p_outlier_SA928 = plotCopynumberHeatmap(CN[, outlier_SA928], cluster_rows = "ploidy")

# Relationship rpc and deviation
df_ploidy_outliers = df_ploidy %>% dplyr::group_by(UID, SLX) %>% dplyr::mutate(median_ploidy = median(ploidy)) %>% dplyr::ungroup() %>%
  dplyr::mutate(diff_ploidy = ploidy - median_ploidy)

p_ploidy_rpc = ggplot(data = df_ploidy_outliers %>% dplyr::filter(rpc > 25) %>% dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"))) + 
  geom_point(aes(x=used.reads, y=diff_ploidy, color=technology)) +
  facet_wrap(~technology+UID2) +
  theme_pubclean()
p_ploidy_rpc

p_cellline_reads = ggplot(data = df_ploidy %>% 
                               dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"), !startsWith(UID, "UID-DLP-SA1090"), !startsWith(UID, "UID-DLP-SA1101"), !startsWith(UID, "UID-DLP-SA1088"),
                                             !startsWith(UID, "UID-DLP-SA1135"), !startsWith(UID, "UID-DLP-SA928"), !startsWith(UID, "UID-DLP-SA039"), !startsWith(UID, "UID-DLP-SA906"))) +
  geom_quasirandom(aes(x=UID2, y = rpc, color=technology, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=UID2, y=rpc), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ technology, scales = "free_x", strip.position="bottom", nrow=1) +
  geom_hline(yintercept = 25, color="red", linetype="dashed") +
#  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("rpc") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside")
p_cellline_reads


p_ploidy_prediction_tumor = ggplot(data = df_ploidy %>% dplyr::filter(rpc > 25) %>%
    dplyr::filter(startsWith(UID, "UID-NNA-T") | startsWith(UID, "UID-10X-BR") | startsWith(UID, "UID-DLP-SA1135"))) +
  geom_quasirandom(aes(x=UID2, y = ploidy, color=technology, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=UID2, y=ploidy), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ technology, scales = "free_x", strip.position="bottom", nrow=1) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 90))
p_ploidy_prediction_tumor

p_tumor_reads = ggplot(data = df_ploidy %>% dplyr::filter(startsWith(UID, "UID-NNA-T") | startsWith(UID, "UID-10X-BR") | startsWith(UID, "UID-DLP-SA1135"))) + #!(UID %in% c("UID-NNA-TN6-paired", "UID-NNA-TN7-paired", "UID-NNA-TN8")))) +
  geom_quasirandom(aes(x=UID2, y = rpc, color=technology, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=UID2, y=rpc), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ technology, scales = "free_x", strip.position="bottom", nrow=1) +
  geom_hline(yintercept = 25, color="red", linetype="dashed") +
#  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("rpc") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside")
p_tumor_reads


## Add HMMCopy derived copy number values for comparison
ploidy_table_ACT_2 = data.frame(UID2=c("BT20", "mb157", "mb231", "mb453", "TN1", "TN2", "TN3", "TN4", "TN5", "TN6", "TN7", "TN8"),
                          ploidy=c(2.70, 2.55, 2.41, 4.17, 3.45, 3.03, 3.44, 3.76, 2.65, 3.17, 3.15, 3.95),
                          technology=c(rep("ACT", 12))) %>%
  dplyr::mutate(ymin = ploidy - 0.2, ymax=ploidy+0.2) %>%
  dplyr::mutate(UID3 = factor(UID2, ordered = TRUE, levels=UID_levels),
                dtype=ifelse(base::startsWith(as.character(UID2), "TN"), "patient", "cellline"))

ploidy_table_ACT_2$method = "groundtruth"
ploidy_table_ACT = dplyr::bind_rows(ploidy_table_ACT_2, ploidy_table_ACT_2)#, ploidy_table_ACT_2)
ploidy_table_ACT$method = c(rep("scAbsolute", times=nrow(ploidy_table_ACT_2)),
                            rep("hmmcopy", times=nrow(ploidy_table_ACT_2)))
                            #rep("groundtruth", times=nrow(ploidy_table_ACT_2)))
# groundtruth data for ACT data (from supplement, based on FACS)

## hmmcopy data for comparison
sample_mb453 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-mb453/NAVINACT_segments.csv.gz")
hmmcopy_d1 = sample_mb453 %>% dplyr::mutate(sample="mb453", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_BT20 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-BT20/NAVINACT_segments.csv.gz")
hmmcopy_d2 = sample_BT20 %>% dplyr::mutate(sample="BT20", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_mb157= readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-mb157/NAVINACT_segments.csv.gz")
hmmcopy_d3 = sample_mb157 %>% dplyr::mutate(sample="mb157", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_MDAMB231c28 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-MDAMB231c28/NAVINACT_segments.csv.gz")
hmmcopy_d4 = sample_MDAMB231c28 %>% dplyr::mutate(sample="mb231", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

# tumour samples
sample_TN1 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN1/NAVINACT_segments.csv.gz")
hmmcopy_d5 = sample_TN1 %>% dplyr::mutate(sample="TN1", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN2 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN2/NAVINACT_segments.csv.gz")
hmmcopy_d6 = sample_TN2 %>% dplyr::mutate(sample="TN2", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN3 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN3/NAVINACT_segments.csv.gz")
hmmcopy_d7 = sample_TN3 %>% dplyr::mutate(sample="TN3", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN4 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN4/NAVINACT_segments.csv.gz")
hmmcopy_d8 = sample_TN4 %>% dplyr::mutate(sample="TN4", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN5 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN5/NAVINACT_segments.csv.gz")
hmmcopy_d9 = sample_TN5 %>% dplyr::mutate(sample="TN5", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN6 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN6/NAVINACT_segments.csv.gz")
hmmcopy_d10 = sample_TN6 %>% dplyr::mutate(sample="TN6", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN7 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN7/NAVINACT_segments.csv.gz")
hmmcopy_d11 = sample_TN7 %>% dplyr::mutate(sample="TN7", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

sample_TN8 = readr::read_csv(file="~/Data/project1/HMMCopy/output/UID-NNA-TN8/NAVINACT_segments.csv.gz")
hmmcopy_d12 = sample_TN8 %>% dplyr::mutate(sample="TN8", size=end-start, algorithm="HMMCopy") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(hmmcopy= stats::weighted.mean(state, size))

hmmcopy_df = dplyr::bind_rows(hmmcopy_d1, hmmcopy_d2, hmmcopy_d3, hmmcopy_d4, hmmcopy_d5, hmmcopy_d6,
                              hmmcopy_d7, hmmcopy_d8, hmmcopy_d9, hmmcopy_d10, hmmcopy_d11, hmmcopy_d12)

# Ginkgo results
ginkgo_d1 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-BT20_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="BT20", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d2 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-mb157_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="mb157", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d3 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-MDAMB231c28_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="mb231", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d4 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-mb453_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="mb453", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d5 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN1_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN1", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d6 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN2_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN2", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d7 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN3_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN3", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d8 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN4_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN4", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d9 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN5_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN5", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d10 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN6_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN6", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d11 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN7_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN7", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_d12 = readr::read_tsv(file="~/Data/project1/Ginkgo/UID-NNA-TN8_SegCopy.txt") %>%
  tidyr::pivot_longer(cols=dplyr::starts_with("UID.NNA"), names_to="cell_id", values_to="state") %>%
  dplyr::mutate(sample="TN8", size=END-START, algorithm="ginkgo") %>% dplyr::group_by(cell_id, sample, algorithm) %>%
  dplyr::summarise(ginkgo=stats::weighted.mean(state, size))

ginkgo_df = dplyr::bind_rows(ginkgo_d1, ginkgo_d2, ginkgo_d3, ginkgo_d4, ginkgo_d5, ginkgo_d6,
                             ginkgo_d7, ginkgo_d8, ginkgo_d9, ginkgo_d10, ginkgo_d11, ginkgo_d12) %>%
  dplyr::mutate(cell_id = gsub("\\.", "-", cell_id))


# CHISEL
CHISEL_PATH="~/Data/project1/CHISEL/NONORMAL/"
# note: impute data have slightly better performance overall than just phasing, so we use those
# note: we add pseudolables to the cells, we do not care about individual cells
chisel_d1 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN1/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN1", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d1$cell_id = gsub("\\.", "-", ginkgo_d5$cell_id[1:nrow(chisel_d1)])

chisel_d2 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN2/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN2", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d2$cell_id = gsub("\\.", "-", ginkgo_d6$cell_id[1:nrow(chisel_d2)])

chisel_d3 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN3/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN3", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d3$cell_id = gsub("\\.", "-", ginkgo_d7$cell_id[1:nrow(chisel_d3)])

chisel_d4 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN4/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN4", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d4$cell_id = gsub("\\.", "-", ginkgo_d8$cell_id[1:nrow(chisel_d4)])

chisel_d5 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN5/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN5", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d5$cell_id = gsub("\\.", "-", ginkgo_d9$cell_id[1:nrow(chisel_d5)])

chisel_d6 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN6/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN6", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d6$cell_id = gsub("\\.", "-", ginkgo_d10$cell_id[1:nrow(chisel_d6)])

chisel_d7 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN7/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN7", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d7$cell_id = gsub("\\.", "-", ginkgo_d11$cell_id[1:nrow(chisel_d7)])

chisel_d8 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-TN8/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="TN8", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d8$cell_id = gsub("\\.", "-", ginkgo_d12$cell_id[1:nrow(chisel_d8)])
  
chisel_d9 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-BT20/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="BT20", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d9$cell_id = gsub("\\.", "-", ginkgo_d1$cell_id[1:nrow(chisel_d9)])

chisel_d10 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-mb157/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="mb157", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d10$cell_id = gsub("\\.", "-", ginkgo_d2$cell_id[1:nrow(chisel_d10)])

chisel_d11 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-MDAMB231c28/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="mb231", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d11$cell_id = gsub("\\.", "-", ginkgo_d3$cell_id[1:nrow(chisel_d11)])

chisel_d12 = readr::read_tsv(file = paste0(CHISEL_PATH, "UID-NNA-mb453/calls/calls.tsv")) %>%
  dplyr::rename("CHR"="#CHR") %>% 
  tidyr::separate(col="CORRECTED_HAP_CN", into=c("cn_A", "cn_B"), sep="\\|", convert=TRUE, remove=FALSE) %>%
  dplyr::mutate(sample="mb453", total_copynumber = cn_A + cn_B, size=END-START, algorithm="CHISEL") %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(total_copynumber, size))
chisel_d12$cell_id = gsub("\\.", "-", ginkgo_d4$cell_id[1:nrow(chisel_d12)])

chisel_df = dplyr::bind_rows(chisel_d1, chisel_d2, chisel_d3, chisel_d4,chisel_d5, chisel_d6, chisel_d7, chisel_d8,
                             chisel_d9, chisel_d10, chisel_d11, chisel_d12)

## Compare HMMCopy with scAbsolute calls
compare_algorithms = dplyr::left_join(dplyr::inner_join(dplyr::inner_join(dplyr::select(df1 %>% #predict_replicating(df1) %>% dplyr::filter(!replicating) %>%
                            dplyr::mutate(algorithm="scAbsolute", scAbsolute=ploidy), name, scAbsolute, UID2, SLX, technology, used.reads),
                            hmmcopy_df, by=c("name"="cell_id")), ginkgo_df, by=c("name"="cell_id")), chisel_df, by=c("name"="cell_id")) %>% 
  tidyr::pivot_longer(c("hmmcopy", "scAbsolute", "ginkgo", "chisel"), names_to="method", values_to="ploidy") %>%
  dplyr::mutate(UID3 = factor(UID2, ordered = TRUE, levels=UID_levels),
                dtype=ifelse(base::startsWith(as.character(UID3), "TN"), "patient", "cellline"))
  
##ploidy_table_ACT_3 = ploidy_table_ACT %>% #%>% dplyr::filter(UID2 %in% c("BT20", "mb157", "mb231", "mb453", "TN1", "TN3"))
#baseline_comparison_cellline = ggplot(data=compare_algorithms %>% dplyr::filter(dtype=="cellline")) +
#  geom_quasirandom(aes(y=method, color=method, x=ploidy), alpha=1.0, size=1.0, groupOnX = FALSE) +
#  ##geom_density(aes(x=ploidy, y=UID2, color=method), stat = "density", adjust=2, color="red") +
#  geom_tile(data=ploidy_table_ACT %>% dplyr::filter(dtype=="cellline") , aes(y=method,x=ploidy, height=1.0,width=1.0),alpha=0.5,fill="grey") +
#  geom_point(data=ploidy_table_ACT %>% dplyr::filter(dtype=="cellline") %>% dplyr::filter(method=="scAbsolute"), aes(y=method, x=ploidy), position=position_nudge(x = 0, y = -0.5), color="red", shape=3, size=2) +
#  ##geom_tile(data=ploidy_table_ACT_3, aes(y=,x=ploidy, height=0.8,width=0.2),alpha=0.95,fill="grey") +
#  facet_wrap(.~UID3, scales = "free_y", drop=TRUE, ncol=1, strip.position = "left") +
#  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
#  theme_pubclean(base_size = 20, flip = TRUE) +
#  scale_color_manual(values = c("#E69F00", "#56B4E9"), labels=c("HMMCopy ", "scAbsolute ")) +
#  theme(legend.position="top",
#        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
#        strip.background = element_blank(), 
#        strip.text.y.left = element_text(angle = 0),
#        panel.spacing = unit(1, "lines"))  +
#  guides(shape="none",
#         colour = guide_legend(title.position="left", #title.hjust = 0.5, #title.align=0.0, 
#                               keywidth = 1.5, title="Method", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
#  theme(legend.key = element_rect(color = "white", fill="white"),
#        legend.title = element_text(#family = "Playfair",
#                                    #color = "chocolate",
#                                    size = 14, face = 2))
# 
#baseline_comparison_cellline
#
#
#baseline_comparison_tumour = ggplot(data=compare_algorithms %>% dplyr::filter(dtype == "patient")) +
#  geom_quasirandom(aes(y=method, color=method, x=ploidy), alpha=1.0, size=1.0, groupOnX = FALSE) +
#  ##geom_density(aes(x=ploidy, y=UID2, color=method), stat = "density", adjust=2, color="red") +
#  geom_tile(data=ploidy_table_ACT %>% dplyr::filter(dtype=="patient"), aes(y=method,x=ploidy, height=1.0,width=1.0),alpha=0.5,fill="grey") +
#  geom_point(data=ploidy_table_ACT %>% dplyr::filter(dtype=="patient") %>% dplyr::filter(method=="scAbsolute"), aes(y=method, x=ploidy), position=position_nudge(x = 0, y = -0.5), color="red", shape=3, size=2) +
#  ##geom_tile(data=ploidy_table_ACT_3, aes(y=,x=ploidy, height=0.8,width=0.2),alpha=0.95,fill="grey") +
#  facet_wrap(.~UID3, scales = "free_y", drop=TRUE, ncol=1, strip.position = "left") +
#  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
#  theme_pubclean(base_size = 20, flip = TRUE) +
#  scale_color_manual(values = c("#E69F00", "#56B4E9"), labels=c("HMMCopy ", "scAbsolute ")) +
#  theme(legend.position="top",
#        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
#        strip.background = element_blank(), 
#        strip.text.y.left = element_text(angle = 0),
#        panel.spacing = unit(1, "lines"))  +
#  guides(shape="none",
#         colour = guide_legend(title.position="left", #title.hjust = 0.5, #title.align=0.0, 
#                               keywidth = 1.5, title="Method", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
#  theme(legend.key = element_rect(color = "white", fill="white"),
#        legend.title = element_text(#family = "Playfair",
#                                    #color = "chocolate",
#                                    size = 14, face = 2))
#
#baseline_comparison_tumour

## results in Table 
stats_table = dplyr::left_join(compare_algorithms, ploidy_table_ACT_2, by="UID3", suffix=c("", ".groundtruth")) %>%
  dplyr::mutate(rpc_unbiased = used.reads / (ploidy * 6206))
# not possible to filter for read depth (problem is this introduces bias in which methods are used)
  #dplyr::filter(rpc_unbiased > 25)
table(stats_table$UID3, stats_table$method)
tb = table(stats_table$ploidy > stats_table$ploidy.groundtruth + 0.5 | 
      stats_table$ploidy < stats_table$ploidy.groundtruth - 0.5,
      stats_table$method, stats_table$UID3)
tb
results = stats_table %>% dplyr::group_by(UID3, method) %>% dplyr::summarise(n=n(), md = mean(abs(ploidy-ploidy.groundtruth)), sd=sd(abs(ploidy-ploidy.groundtruth)))
print(results, n=100)
prop.table(tb, c(3,2))

#results %>% dplyr::filter(method == "hmmcopy")


## double check TN5
#table(stats_table %>% dplyr::filter(UID3 == "TN5") %>% dplyr::filter(ploidy > 3) %>% dplyr::pull(method))


## Create a single figure with everything

# clean outlier cells
dat = dplyr::bind_rows(df_ploidy %>% dplyr::mutate(method="scAbsolute") %>% dplyr::filter(technology != "ACT") %>% 
                               dplyr::mutate(UID3=factor(UID2, ordered=TRUE, levels=UID_levels), technology = factor(technology, ordered=TRUE, levels=c("10X", "DLP", "ACT - cell lines",  "ACT - tumours"))) %>%
                               dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"), !startsWith(UID, "UID-DLP-SA1090"), !startsWith(UID, "UID-DLP-SA1101"), !startsWith(UID, "UID-DLP-SA1088"),
                                             !startsWith(UID, "UID-DLP-SA1135"), !startsWith(UID, "UID-DLP-SA928"), !startsWith(UID, "UID-DLP-SA039"), !startsWith(UID, "UID-DLP-SA906")),
                       compare_algorithms %>% dplyr::mutate(technology2=dplyr::case_when(technology == "ACT" & dtype == "patient" ~ "ACT - tumours",
                                                                                  technology == "ACT" & dtype!="patient" ~ "ACT - cell lines",
                                                                                  TRUE ~ "NA")) %>%
  dplyr::mutate(technology = factor(technology2, ordered=TRUE, levels=c("10X", "DLP", "ACT - cell lines", "ACT - tumours")))) %>%
  dplyr::mutate(method2 = factor(method, ordered=TRUE, levels=c("exp", "hmmcopy", "scAbsolute", "ginkgo", "chisel"))) %>%
  dplyr::mutate(rpc_unbiased = used.reads / (ploidy * 6206)) %>%
  dplyr::filter(rpc_unbiased > 25) %>%
  dplyr::select(-method) %>% dplyr::rename(method=method2) %>%
  dplyr::mutate(method2 = base::interaction(UID3, method)) %>%
  dplyr::mutate(grid1 = ifelse(technology %in% c("10X", "DLP"), "noACT", "yesACT"), 
                grid2 = dplyr::case_when(technology == "ACT - cell lines" ~ "left",
                                         technology == "ACT - tumours" ~ "right",
                                         technology == "10X" ~ "left",
                                         technology == "DLP" ~ "right",
                                         TRUE ~ "other"
                                         ))

table(dat$method2)
levels(dat$method)
levels(dat$method2)

#custom_labeller = 


Fig_ploidy_I_a = ggplot(data = dat %>% dplyr::filter(technology %in% c("10X", "DLP"))) +
  geom_point(aes(x=1, y=NA, color=factor(method, levels=c("exp", "hmmcopy", "scAbsolute", "ginkgo", "chisel"))), size=0) +
  geom_quasirandom(aes(y=base::interaction(UID2, method), color=method, x=ploidy), alpha=1.0, size=1.0, groupOnX = FALSE) +
  geom_tile(data=ploidy_table %>% dplyr::filter(technology=="10X") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=0.8),alpha=0.5,fill="grey") +
  geom_tile(data=ploidy_table %>% dplyr::filter(technology=="DLP") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=0.8),alpha=0.0,fill="white") +
  geom_point(data=ploidy_table %>% dplyr::filter(technology=="10X") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy), size=2, shape=3, color="red") +
  geom_point(data=ploidy_table %>% dplyr::filter(technology=="10X") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*0.5), size=2, shape=8, color="blue", position=position_nudge(x=0,y=0.0)) +
  geom_point(data=ploidy_table %>% dplyr::filter(technology=="10X") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*2.0), size=2, shape=8, color="blue", position=position_nudge(x=0,y=0.0)) +
  #geom_tile(data=ploidy_table %>% dplyr::filter(startsWith(as.character(technology), "ACT")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=2.4),alpha=0.5,fill="grey", position=position_nudge(x=0,y=-1.0)) +
  #geom_point(data=ploidy_table %>% dplyr::filter(startsWith(as.character(technology), "ACT")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy), size=2, shape=3, color="red", position=position_nudge(x=0,y=-1.0)) +
  #geom_point(data=ploidy_table %>% dplyr::filter(startsWith(as.character(technology), "ACT")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*0.5), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.0)) +
  #geom_point(data=ploidy_table %>% dplyr::filter(startsWith(as.character(technology), "ACT")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*2.0), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.0)) +
  facet_wrap(~technology, scales = "free_y", drop=TRUE, strip.position = "top", ncol=2) +
  theme_pubclean(base_size = 20, flip=TRUE) + theme(legend.position="right") + ylab("sample") + xlab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  guides(shape="none", colour = guide_legend(title="Method  ", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.box.just = "right",
        legend.margin = margin(c(5, 0, 5, 0)),
        legend.text = element_text(margin = margin(r = 10, unit = "pt")),
        legend.position = "none",
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 14, face = 2)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), labels=c("1", "2", "3", "4", "5", "6", "7")) +
  coord_cartesian(xlim=c(1, 7)) +
  scale_y_discrete(labels=function(x){
    #print(x);
    #return(x)
    if(startsWith(x, "BT") || startsWith(x, "mb") || startsWith(x, "TN")){
      return(ifelse(endsWith(x, ".scAbsolute") | endsWith(x, ".ginkgo"), "", gsub(".hmmcopy", "", x)))
    }else{
      return(ifelse(endsWith(x, ".hmmcopy") | endsWith(x, ".ginkgo"), "", gsub(".scAbsolute", "", x)))
    }
    }) +
  scale_color_manual(values = c("exp"="#DCDCDC", "scAbsolute"="#E69F00", "hmmcopy"="#56B4E9", "ginkgo"="#009E73", "chisel"="#F0E442"),
                     labels=c("exp"="Experimental ", "hmmcopy"="HMMCopy ", "ginkgo"="Ginkgo ", "chisel"="CHISEL", "scAbsolute"="scAbsolute "))

Fig_ploidy_I_a


Fig_ploidy_I_b = ggplot(data = dat %>% dplyr::filter(technology %in% c("ACT - cell lines", "ACT - tumours"))) +
  geom_point(aes(x=1, y=NA, color=factor(method, levels=c("exp", "hmmcopy", "scAbsolute", "ginkgo", "chisel"))), size=0) +
  geom_quasirandom(aes(y=base::interaction(UID2, method), color=method, x=ploidy), alpha=1.0, size=1.0, groupOnX = FALSE) +
  #geom_tile(data=ploidy_table %>% dplyr::filter(technology=="10X") %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=0.8),alpha=0.5,fill="grey") +
  geom_tile(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "lines")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=3.6),alpha=0.5,fill="grey", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "lines")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy), size=2, shape=3, color="red", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "lines")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*0.5), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "lines")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*2.0), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.5)) +
  geom_tile(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "tumours")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method),x=ploidy, width=1.0,height=3.6),alpha=0.5,fill="grey", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "tumours")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy), size=2, shape=3, color="red", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "tumours")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*0.5), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.5)) +
  geom_point(data=ploidy_table %>% dplyr::filter(endsWith(as.character(technology), "tumours")) %>% dplyr::mutate(method="scAbsolute"), aes(y=interaction(UID2, method), x=ploidy*2.0), size=2, shape=8, color="blue", position=position_nudge(x=0,y=-1.5)) +
  facet_wrap(~technology, scales = "free_y", drop=TRUE, strip.position = "top", ncol=2) +
  theme_pubclean(base_size = 20, flip=TRUE) + theme(legend.position="right") + ylab("sample") + xlab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  guides(shape="none", colour = guide_legend(title="Method ", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
  theme(axis.text.y = element_text(vjust=1.0),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(color = "white", fill="white"),
        legend.position = "top",
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 14, face = 2)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), labels=c("1", "2", "3", "4", "5", "6", "7")) +
  coord_cartesian(xlim=c(1, 7)) +
  scale_y_discrete(labels=function(x){
    #print(x);
    #return(x)
    if(startsWith(x, "BT") || startsWith(x, "mb") || startsWith(x, "TN")){
      return(ifelse(endsWith(x, ".scAbsolute") | endsWith(x, ".ginkgo") | endsWith(x, ".chisel"), "", gsub(".hmmcopy", "", x)))
    }else{
      return(ifelse(endsWith(x, ".hmmcopy") | endsWith(x, ".ginkgo") | endsWith(x, ".chisel"), "", gsub(".scAbsolute", "", x)))
    }
    }) +
  scale_color_manual(values = c("exp"="#DCDCDC", "scAbsolute"="#E69F00", "hmmcopy"="#56B4E9", "ginkgo"="#009E73", "chisel"="#F0E442"),
                     labels=c("exp"="Experimental ", "hmmcopy"="HMMCopy ", "ginkgo"="Ginkgo ", "chisel"="CHISEL", "scAbsolute"="scAbsolute "))

Fig_ploidy_I_b

# todo combine top and bottom figure for different size scales

leg = get_legend(Fig_ploidy_I_b)
Fig_ploidy_I = ggpubr::ggarrange(leg, 
                                 Fig_ploidy_I_a + rremove("xlab") + rremove("ylab"),
                                 Fig_ploidy_I_b + rremove("ylab") + theme(legend.position = "none"),
                                 nrow=3, labels=c("", "", ""), heights = c(1, 3, 5))
Fig_ploidy_I = annotate_figure(Fig_ploidy_I, left = textGrob("sample", rot = 90, vjust = 0.5, gp = gpar(cex = 1.5)))
Fig_ploidy_I

ggpubr::ggexport(Fig_ploidy_I, filename = "~/scAbsolute/figures/Fig_ploidy_I.pdf", width=12, height=12)

#baseline_comparison_tumour = ggplot(data=compare_algorithms %>% dplyr::filter(dtype == "patient")) +
#  geom_quasirandom(aes(y=method, color=method, x=ploidy), alpha=1.0, size=1.0, groupOnX = FALSE) +
#  ##geom_density(aes(x=ploidy, y=UID2, color=method), stat = "density", adjust=2, color="red") +
#  geom_tile(data=ploidy_table_ACT %>% dplyr::filter(dtype=="patient"), aes(y=method,x=ploidy, height=1.0,width=1.0),alpha=0.5,fill="grey") +
#  geom_point(data=ploidy_table_ACT %>% dplyr::filter(dtype=="patient") %>% dplyr::filter(method=="scAbsolute"), aes(y=method, x=ploidy), position=position_nudge(x = 0, y = -0.5), color="red", shape=3, size=2) +
#  ##geom_tile(data=ploidy_table_ACT_3, aes(y=,x=ploidy, height=0.8,width=0.2),alpha=0.95,fill="grey") +
#  facet_wrap(.~UID3, scales = "free_y", drop=TRUE, ncol=1, strip.position = "left") +
#  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
#  theme_pubclean(base_size = 20, flip = TRUE) +
#  scale_color_manual(values = c("#E69F00", "#56B4E9"), labels=c("HMMCopy ", "scAbsolute ")) +
#  theme(legend.position="top",
#        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
#        strip.background = element_blank(), 
#        strip.text.y.left = element_text(angle = 0),
#        panel.spacing = unit(1, "lines"))  +
#  guides(shape="none",
#         colour = guide_legend(title.position="left", #title.hjust = 0.5, #title.align=0.0, 
#                               keywidth = 1.5, title="Method", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
#  theme(legend.key = element_rect(color = "white", fill="white"),
#        legend.title = element_text(#family = "Playfair",
#                                    #color = "chocolate",
#                                    size = 14, face = 2))
#


## Create panel figure
#leg = get_legend(p_ploidy_prediction)
#leg_method = get_legend(baseline_comparison)
#Fig_ploidy_I = ggpubr::ggarrange(p_ploidy_prediction + theme(legend.position = "none"),
#                                 baseline_comparison_cellline + ggtitle("ACT"),
#                                 baseline_comparison_tumour + theme(legend.position = "none"),
#                                 nrow=3, ncol=1, heights = c(1,1,1), labels=c("A", "B", ""))
#                                 #ggpubr::ggarrange(
#                                 #  ggpubr::ggarrange(leg_method),
#                                 #  baseline_comparison + theme(legend.position = "none"), heights = c(1, 12),
#                                 #                  ncol=1, labels=c("", "B"), vjust=8),
#Fig_ploidy_I
#
#ggpubr::ggexport(Fig_ploidy_I, filename = "~/scAbsolute/figures/Fig_ploidy_I.pdf", width=9, height=18)
#leg = get_legend(p3_nice)
#leg = get_legend(p_ploidy_prediction2)
#Fig_ploidy_I = ggpubr::ggarrange(p_ploidy_prediction,
#           ggpubr::ggarrange(
#             ggpubr::ggarrange(p3_nice + theme(legend.position="none") + rremove("x.title"),
#                               p_predict_nice + theme(legend.position="none") + theme(strip.background = element_blank(),
#                                                      strip.text.x = element_blank()), nrow=2, labels=c("", "C")),
#                               p_ploidy_prediction2 + theme(legend.position="left"), ncol=2,
#             labels=c("", "D"), hjust=1.1, vjust=1.5, widths = c(3,2)),
#                               labels=c("A", "B"), nrow=2, ncol=1,heights = c(1, 2))
#Fig_ploidy_I

## BEGIN CHISEL supplementary figure ====
## Supplementary figure for CHISEL

# use 1000kb BREAST data as very shallow sequencing
files2 = c(
  paste0("~/Data/project1/30-scale/", 1000, "/predict/", 
         paste0(c(
           "UID-10X-BREAST-A",
           "UID-10X-BREAST-B",
           "UID-10X-BREAST-C",
           "UID-10X-BREAST-D",
           "UID-10X-BREAST-E"
         ), "_", 1000, ".rds")))
pdB = dplyr::bind_rows(lapply(files2, 
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", "")) %>%
  dplyr::mutate(cycling_activity = cellcycle.kendall.repTime.weighted.median.cor.corrected)
CNB = combineQDNASets(lapply(files2, function(x) readRDS(x))) 

chisel_df_A_full = readr::read_tsv("~/Data/project1/CHISEL/chisel-data/patientS0/calls/sectionA/calls.tsv.gz", col_types = "ciiciidiidic")
chisel_df_A = chisel_df_A_full %>%
  tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
  dplyr::mutate(sample="UID-10X-BREAST-A", size=END-START, algorithm="CHISEL", cn=cn_a+cn_b) %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(cn, size))

chisel_df_B_full = readr::read_tsv("~/Data/project1/CHISEL/chisel-data/patientS0/calls/sectionB/calls.tsv.gz", col_types = "ciiciidiidic")
chisel_df_B = chisel_df_B_full %>%
  tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
  dplyr::mutate(sample="UID-10X-BREAST-B", size=END-START, algorithm="CHISEL", cn=cn_a+cn_b) %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(cn, size))

chisel_df_C_full = readr::read_tsv("~/Data/project1/CHISEL/chisel-data/patientS0/calls/sectionC/calls.tsv.gz", col_types = "ciiciidiidic")
chisel_df_C = chisel_df_C_full %>%
  tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
  dplyr::mutate(sample="UID-10X-BREAST-C", size=END-START, algorithm="CHISEL", cn=cn_a+cn_b) %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(cn, size))

chisel_df_D_full = readr::read_tsv("~/Data/project1/CHISEL/chisel-data/patientS0/calls/sectionD/calls.tsv.gz", col_types = "ciiciidiidic")
chisel_df_D = chisel_df_D_full %>%
  tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
  dplyr::mutate(sample="UID-10X-BREAST-D", size=END-START, algorithm="CHISEL", cn=cn_a+cn_b) %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(cn, size))

chisel_df_E_full = readr::read_tsv("~/Data/project1/CHISEL/chisel-data/patientS0/calls/sectionE/calls.tsv.gz", col_types = "ciiciidiidic")
chisel_df_E = chisel_df_E_full %>%
  tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
  dplyr::mutate(sample="UID-10X-BREAST-E", size=END-START, algorithm="CHISEL", cn=cn_a+cn_b) %>% dplyr::group_by(CELL, sample, algorithm) %>%
  dplyr::summarise(chisel=stats::weighted.mean(cn, size))

df_CHISEL = dplyr::bind_rows(chisel_df_A, chisel_df_B, chisel_df_C, chisel_df_D, chisel_df_E) %>% 
  dplyr::group_by(sample) %>%
  tidyr::pivot_longer(chisel, names_to="method", values_to="ploidy") %>%
  dplyr::rename(name="CELL") %>% dplyr::select(name, sample, method, ploidy)
df_scAbsolute = pdB %>% dplyr::filter(startsWith(UID, "UID-10X-BREAST")) %>% dplyr::select(UID, ploidy, celltag, used.reads) %>%
  dplyr::rename(scAbsolute="ploidy", sample="UID") %>% dplyr::mutate(name = gsub("-1$", "", celltag)) %>% dplyr::select(-celltag) %>%
  dplyr::group_by(sample) %>%
  tidyr::pivot_longer(scAbsolute, names_to="method", values_to="ploidy") %>% dplyr::select(name, sample, method, ploidy, used.reads)
df_sup_CHISEL = dplyr::bind_rows(df_CHISEL, df_scAbsolute)

breast_replicating = gsub("-1", "", predict_replicating(pdB %>% dplyr::filter(startsWith(UID, "UID-10X-BREAST"))) %>%
  dplyr::filter(replicating | used.reads<250000) %>% dplyr::pull(celltag))

ggplot(data=pdB %>% dplyr::filter(startsWith(UID, "UID-10X-BR"))) +
  geom_quasirandom(aes(x=UID, y=used.reads))

ggplot(data=df_sup_CHISEL) +
  geom_quasirandom(aes(x=method, y=ploidy, color=as.factor(method))) +
  facet_wrap(~sample) +
  theme_pubclean() + coord_cartesian(ylim=c(0, 8))

df_sup_CHISEL_joined = dplyr::inner_join(df_CHISEL, df_scAbsolute, by=c("name"="name","sample"="sample"), suffix=c(".chisel", ".scAbsolute")) %>%
         dplyr::filter(ploidy.chisel <= 8) %>% dplyr::mutate(sample = as.factor(sample)) %>%
  dplyr::group_by(sample) %>% dplyr::slice_max(used.reads, prop=0.75) %>% dplyr::ungroup()


ggplot(data=df_sup_CHISEL_joined) +
  geom_quasirandom(aes(x=sample, y=used.reads))

gc_correction = TRUE
# look into some of the outliers (false ploidy outliers)
set.seed(2022)
examples = df_sup_CHISEL_joined %>% dplyr::filter(!(name %in% breast_replicating), ploidy.scAbsolute < 2.30, ploidy.scAbsolute > 1.70, ploidy.chisel > 3.0, ploidy.chisel < 6) %>%
 dplyr::ungroup() %>% dplyr::slice_sample(n=12) %>% dplyr::mutate(fulltag=paste0(name, "-1")) %>% dplyr::select(fulltag, sample)
examples1 = dplyr::inner_join(examples, df1, by=c("fulltag"="celltag", "sample"="UID")) %>%
  dplyr::mutate(shorttag = gsub("-1", "", fulltag)) %>%
  dplyr::select(name, ploidy, fulltag, shorttag, sample, rpc)
names1 = examples1 %>% dplyr::filter(rpc > 25) %>% dplyr::pull(name)

# helper function to plot chisel calls on top of other data
addCHISEL_segmentation <- function(CN1, name){
  object = CN1[, name]
  l = str_split(name, "_") 
  shortt = gsub("-1", "", l[[1]][[4]])
  samplet = l[[1]][[1]]
  dt = switch(samplet, "UID-10X-BREAST-A"=chisel_df_A_full,
                       "UID-10X-BREAST-B"=chisel_df_B_full,
                       "UID-10X-BREAST-C"=chisel_df_C_full,
                       "UID-10X-BREAST-D"=chisel_df_D_full,
                       "UID-10X-BREAST-E"=chisel_df_E_full)
  
  a = dt %>% dplyr::filter(CELL == shortt) %>% dplyr::rename("chr2"="#CHR") %>% dplyr::mutate(chr = as.integer(gsub("chr", "", chr2))) %>%
    dplyr::arrange(chr) %>% 
    tidyr::separate(CN_STATE, into=c("cn_a", "cn_b"), sep="\\|") %>% dplyr::mutate(cn_a = as.integer(cn_a), cn_b = as.integer(cn_b)) %>%
    dplyr::mutate(TOTAL=cn_a+cn_b)
  
  sub = object[, 1]
  call = GenomicRanges::GRanges(seqnames = a$chr, ranges = IRanges(a$START+1, a$END), total=a$TOTAL)
  ref = createGR(sub)
  hits = findOverlaps(call, ref)
  
  b = sub@assayData$segmented
  valid = binsToUseInternal(sub)
  b[valid] = NA
  # needed to fix plotting of X and Y chromosomes (not called in CHISEL)
  b[startsWith(rownames(sub), "X") | startsWith(rownames(sub), "Y"),] = 100
  b[subjectHits(hits)] = call$total[queryHits(hits)]  + 0.05
  object = Biobase::assayDataElementReplace(object, "segmented", b)
  
  return(object)
}

test = addCHISEL_segmentation(CNB, names1[[1]])

chromosome_break_label = c("1", "", "3", "", "5", "", "7", "", "9", "", "11", "", "", "", "15", "", "", "", "19", "", "", "", "", "")
pe1 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[1]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe1 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[1]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe2 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[2]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe3 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[3]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe4 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[4]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe5 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[5]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe6 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[6]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe7 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[7]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pe8 = plotCopynumber(addCHISEL_segmentation(CNB, names1[[8]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))

#pe9 = plotCopynumber(CN[, names[[9]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe10 = plotCopynumber(CN[, names[[10]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe11 = plotCopynumber(CN[, names[[11]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe12 = plotCopynumber(CN[, names[[12]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))


joint_sup_1 = ggpubr::ggarrange(pe1 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe2 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe3 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe4 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe5 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe6 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pe7 +rremove("ylab"),
                              pe8 +rremove("ylab"),
                              nrow=4, ncol=2)
joint_sup_1

set.seed(2022)
examples = df_sup_CHISEL_joined %>% dplyr::filter(!(name %in% breast_replicating), ploidy.chisel < 2.30, ploidy.chisel > 1.70, ploidy.scAbsolute > 3.0, ploidy.scAbsolute < 6.0) %>%
 dplyr::ungroup() %>% dplyr::slice_sample(n=12) %>% dplyr::mutate(fulltag=paste0(name, "-1")) %>% dplyr::select(fulltag, sample)
examples2 = dplyr::inner_join(examples, pdB, by=c("fulltag"="celltag", "sample"="UID")) %>%
  dplyr::select(name, ploidy, fulltag, sample, rpc)
names2 = examples2 %>% dplyr::filter(rpc > 25) %>% dplyr::pull(name)

test = CNB[, names2[[1]]]
pee1 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[1]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee2 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[2]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee3 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[3]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee4 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[4]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))

pee5 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[5]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee6 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[6]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee7 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[7]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
pee8 = plotCopynumber(addCHISEL_segmentation(CNB, names2[[8]]), addSegmentation = TRUE, addCopynumber = TRUE, ylim = c(0, 7), copyColor = "orange", segmentationColor="blue", alphaLevel = 0.5, alphaLevelSeg=0.1,  showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = gc_correction, chromosome_break_label = chromosome_break_label) + coord_cartesian(ylim=c(0, 7), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))

#pe9 = plotCopynumber(CN[, names[[9]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe10 = plotCopynumber(CN[, names[[10]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe11 = plotCopynumber(CN[, names[[11]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))
#pe12 = plotCopynumber(CN[, names[[12]]], addSegmentation = FALSE, addCopynumber = TRUE, ylim = c(0, 8), copyColor = "orange", alphaLevel = 0.3, showUnique = FALSE,showMarker = FALSE, main="",readinfo = FALSE, uniqueColor="orange", correction = TRUE) + coord_cartesian(ylim=c(0, 8), clip="off") +  theme_pubclean(base_size=20) +  theme(plot.margin = margin(1.0, 0, 0, 0, "cm"))


joint_sup_2 = ggpubr::ggarrange(pee1 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee2 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee3 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee4 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee5 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee6 + rremove('xlab')+rremove("x.text")+rremove("ylab") + theme(axis.ticks.x=element_blank()),
                              pee7 +rremove("ylab"),
                              pee8 +rremove("ylab"),
                              nrow=4, ncol=2, heights = c(1,1,1,1.5,1,1,1,1.5))
joint_sup_2


joint_sup_1_f = annotate_figure(joint_sup_1, left = textGrob("Absolute copy number", rot = 90, vjust = 0.1, gp = gpar(cex = 1.5)))
joint_sup_1_f
joint_sup_2_f = annotate_figure(joint_sup_2, left = textGrob("Absolute copy number", rot = 90, vjust = 0.1, gp = gpar(cex = 1.5)))
joint_sup_2_f


sup_chisel_overview = ggplot(data=df_sup_CHISEL_joined %>% dplyr::filter(!(name %in% breast_replicating)) %>% dplyr::mutate(UIDX = gsub("UID-10X-", "", sample))) +
  geom_point(aes(x=ploidy.chisel, y=ploidy.scAbsolute, color=UIDX), size=0.8, alpha=0.6) +
  geom_point(data = dplyr::inner_join(df_sup_CHISEL_joined %>% dplyr::mutate(name=paste0(name, "-1")), dplyr::bind_rows(examples1, examples2),
                                      by=c("name"="fulltag", "sample"="sample")),
             aes(x=ploidy.chisel, y=ploidy.scAbsolute), color="black", size=3.0, shape=1, alpha=1.0, position="jitter") +
  #geom_point(data = dplyr::inner_join(df_sup_CHISEL_joined %>% dplyr::mutate(name=paste0(name, "-1")), dplyr::bind_rows(examples2),
  #                                    by=c("name"="fulltag", "sample"="sample")),
  #           aes(x=ploidy.chisel, y=ploidy.scAbsolute), color="blue", size=3.0, shape=1, alpha=1.0, position="jitter") +
  xlim(1,8) + ylim(1,8) +
  geom_density2d(aes(x=ploidy.chisel, y=ploidy.scAbsolute), color="black", alpha=0.3) +
  xlab("Ploidy estimate - CHISEL") + ylab("Ploidy estimate - scAbsolute") +
  theme_pubclean(base_size = 20) +
  guides(shape="none", colour = guide_legend(title="Section", nrow=1, override.aes = list(alpha = 1, size=3, fill=NA))) +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "top",
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 14, face = 2))
sup_chisel_overview 

# examples are not stratified by sequencing depth
ggplot(data=df_sup_CHISEL_joined) +
  geom_quasirandom(aes(x=sample, y=used.reads), alpha=0.3) +
  geom_point(data=df_sup_CHISEL_joined %>% dplyr::filter(name %in% gsub("-1", "", c(examples1$fulltag, examples2$fulltag))), aes(x=sample, y=used.reads), color="red")


ggpubr::ggexport(sup_chisel_overview, filename="~/scAbsolute/figures/Sup_chisel_10X-BREAST.pdf", width=12, height=8)

#ploidy.scAbsolute < 2.05, ploidy.scAbsolute > 1.95, ploidy.chisel > 2.5, ploidy.chisel < 6
ggpubr::ggexport(joint_sup_1_f, filename="~/scAbsolute/figures/Sup_chisel_examples_1.pdf", width=6, height=9)


# ploidy.chisel < 2.05, ploidy.chisel > 1.95, ploidy.scAbsolute > 2.5, ploidy.scAbsolute < 6 %>%
ggpubr::ggexport(joint_sup_2_f, filename="~/scAbsolute/figures/Sup_chisel_examples_2.pdf", width=6, height=9)


section_E = predict_replicating(pdB%>% dplyr::filter(UID == "UID-10X-BREAST-E")) %>%
  dplyr::filter(!replicating) %>% dplyr::pull(name)
  #dplyr::filter(!(name %in% breast_replicating)) %>% dplyr::pull(name)
p_heatmap_E = plotCopynumberHeatmap(CNB[, section_E], cluster_rows = "ploidy", show_cell_names = FALSE)
p_heatmap_E
ggpubr::ggexport(p_heatmap_E, filename="~/scAbsolute/figures/Sup_chisel_sectionE.pdf")


## END CHISEL supplementary figure

prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$technology),margin = 2)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$UID2),margin = 2)

#Sup_ploidy = ggpubr::ggarrange(p_cellline_reads)
#Sup_ploidy
#ggpubr::ggexport(Sup_ploidy, filename = "~/scAbsolute/figures/Sup_ploidy.pdf", width=16, height=8)
#
#ggpubr::ggexport(p_sup_heatmap_SA1044, filename = "~/scAbsolute/figures/Sup_ploidy_profiles_SA1044.pdf", width=16, height=8)

#ggpubr::ggexport(p_sup_heatmap_SA1135, filename = "~/scAbsolute/figures/Sup_ploidy_profiles_SA1135.pdf", width=16, height=8)

df_select_TN4 = df_ploidy2 %>% dplyr::filter(UID == "UID-NNA-TN4") %>% dplyr::select(name, ploidy) %>% dplyr::mutate(pout = ploidy < 2.5) %>%
  dplyr::arrange(desc(ploidy))
Sup_ploidy_TN4 = plotCopynumberHeatmap(CN[, df_select_TN4$name], row_split = as.factor(df_select_TN4$pout))
Sup_ploidy_TN4

ggpubr::ggexport(Sup_ploidy_TN4, filename = "~/scAbsolute/figures/Sup_ploidy_TN4.pdf", width=12, height=16)

# predict ploidy for cycling cells
df_ploidy_replicating = predict_replicating(df_ploidy2, cutoff_value = 1.5) %>% dplyr::filter(replicating)

df_ploidy_replicating = qc_ploidy(df_ploidy_replicating, range=0.5)
prop.table(table(df_ploidy_replicating$ploidy.outlier, df_ploidy_replicating$technology),margin = 2)
prop.table(table(df_ploidy_replicating$ploidy.outlier, df_ploidy_replicating$UID2),margin = 2)

# clean outlier cells
Sup_ploidy_prediction_replicating = ggplot(data = df_ploidy_replicating %>% dplyr::filter(rpc > 25) %>%
                               dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"), !startsWith(UID, "UID-DLP-SA1090"), !startsWith(UID, "UID-DLP-SA1101"), !startsWith(UID, "UID-DLP-SA1088"),
                                             !startsWith(UID, "UID-DLP-SA1135"), !startsWith(UID, "UID-DLP-SA928"), !startsWith(UID, "UID-DLP-SA039"), !startsWith(UID, "UID-DLP-SA906"))) +
  geom_quasirandom(aes(x=UID2, y = ploidy, color=technology, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=UID2, y=ploidy), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ technology, scales = "free_x", strip.position="bottom", nrow=1) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside")
#  theme(axis.text.x = element_text(angle = 90))
Sup_ploidy_prediction_replicating

ggpubr::ggexport(Sup_ploidy_prediction_replicating, filename = "~/scAbsolute/figures/Sup_ploidy_prediction_replicating.pdf", width=12, height=16)

#ggpubr::ggexport(p_predict_nice_woS, filename = "~/scAbsolute/figures/Sup_ploidy_prediction_no_replicating.pdf", width=12, height=16)


subset_ploidy_replicating_conserved = predict_replicating(df1) %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle %in% c("G1", "S")) %>%
  dplyr::filter(!(cellcycle == "G1" & replicating)) %>% 
  dplyr::mutate(cellcycle_interact = interaction(cellcycle, replicating))

Sup_ploidy_replicating_conserved = ggplot(data = subset_ploidy_replicating_conserved %>% dplyr::filter(!(cellcycle == "G1" & replicating))) +
  geom_quasirandom(aes(x=cellcycle, y=ploidy)) +
  theme_pubclean(base_size = 20) +
  #scale_x_discrete(labels=c("G1", "S-phase (predicted G1)", "S-phase")) +
  scale_x_discrete(labels=c("G1", "S-phase")) +
  xlab("Cell cycle annotation")
Sup_ploidy_replicating_conserved

mean(subset_ploidy_replicating_conserved$ploidy[subset_ploidy_replicating_conserved$cellcycle == "G1"])
prop.table(table(subset_ploidy_replicating_conserved$cellcycle_interact, subset_ploidy_replicating_conserved$ploidy < 2.53 | subset_ploidy_replicating_conserved$ploidy > 3.13),margin = 1)
prop.table(table(subset_ploidy_replicating_conserved$cellcycle, subset_ploidy_replicating_conserved$ploidy < 2.53 | subset_ploidy_replicating_conserved$ploidy > 3.13),margin = 1)

ggpubr::ggexport(Sup_ploidy_replicating_conserved, filename = "~/scAbsolute/figures/Sup_ploidy_replicating_conserved.pdf", width=6, height=4)

#cn_profiles_TN4 = CN[, df_select_TN4$name]
#valid = binsToUseInternal(CN)
#dmat = cn_profiles_TN4@assayData$copynumber[valid, ]
#dmat[, df_select_TN4$pout] = dmat[, df_select_TN4$pout] * 2
## normalize to max state
#dmat_norm = ifelse(dmat > 8, 8, dmat)
#dim(dmat_norm)
#
#dist = copyNumberDistance(t(dmat_norm), BASEDIR="~/")

### Supp figure for ACT examples -> ACT doesn't reach saturation, so not useful in this case
#CN_scaled = combineQDNASets(list(readRDS("~/Data/project1/30-scale/500/scaling/UID-NNA-TN4_500.rds"),
#                                 readRDS("~/Data/project1/30-scale/500/scaling/UID-NNA-MDAMB231c28_500.rds")))
#df_scaled = Biobase::phenoData(CN_scaled)@data %>%
#  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
#  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", "")) %>%
#  dplyr::mutate(cycling_activity = cellcycle.kendall.repTime.weighted.median.cor.corrected)
#  
#df_scaled$scaling = "scaled"
#
#df_division = dplyr::bind_rows(df_scaled, 
#                               df1 %>% dplyr::filter(UID %in% c("UID-NNA-TN4", "UID-NNA-MDAMB231c28")) %>%
#                                 dplyr::mutate(scaling = "predict")) %>%
#  dplyr::mutate(correct_ploidy = dplyr::case_when((UID == "UID-NNA-MDAMB231c28" & ploidy > 3.0) ~ FALSE,
#                                                  (UID == "UID-NNA-TN4" & ploidy < 3.0) ~ FALSE,
#                                                  TRUE ~ TRUE))
#
#ggplot(data = df_division) + geom_point(aes(x=rpc, y=prediction, color=ploidy)) +
#  facet_wrap(~scaling+UID) + theme_pubclean(base_size = 20) +
#  geom_abline(slope=0.00066, intercept=1.0, color="red")
#

