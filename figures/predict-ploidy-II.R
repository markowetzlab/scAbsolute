# Copyright 2022, Michael Schneider, All rights reserved.
## Predict ploidy - demonstrate ability to predict correct ploidy for G2 cells (identify WGD)


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
           "UID-DLP-SA039",
           "UID-DLP-SA1044",
           "UID-DLP-SA928",
           "UID-DLP-SA906",
           "UID-DLP-SA1087",
           "UID-DLP-SA1088",
           "UID-DLP-SA1089",
           "UID-DLP-SA1090",
           "UID-DLP-SA921",
           "UID-DLP-SA922",
           "UID-DLP-SA1101",
           "UID-DLP-SA1135",
           "UID-DLP-SA501X2XB00096",
           "UID-DLP-SA501X2XB00097",
           "UID-10X-BREAST-A",
           "UID-10X-BREAST-B",
           "UID-10X-BREAST-C",
           "UID-10X-BREAST-D",
           "UID-10X-BREAST-E",
           "UID-10X-Andor-2020-48959-MKN-45",
           "UID-10X-Andor-2020-49599-SNU-16", 
           "UID-10X-Andor-2020-49600-SNU-668", 
           "UID-10X-Andor-2020-49606-NCI-N87", 
           "UID-10X-Andor-2020-50941-SNU-638", 
           "UID-10X-Andor-2020-59068-KATOIII", 
           "UID-10X-Andor-2020-49607-HGC-27",
           "UID-10X-Andor-2020-59065-NUGC-4",
           "UID-10X-Andor-2020-59069-SNU-601",
           "UID-NNA-BT20", 
           "UID-NNA-mb157", 
           "UID-NNA-mb453", 
           "UID-NNA-MDAMB231c8",
           "UID-NNA-TN1",
           "UID-NNA-TN2",
           "UID-NNA-TN3",
           "UID-NNA-TN4",
           "UID-NNA-TN5",
           "UID-NNA-TN6",
           "UID-NNA-TN7",
           "UID-NNA-TN6-paired",
           "UID-NNA-TN7-paired",
           "UID-NNA-TN8"
#           "UID-JBL-NA12878",
#           "UID-JBL-CIOV1",
#           "UID-JBL-PEO1",
#           "UID-JBL-OVCAR3"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files, 
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", "")) %>%
  dplyr::mutate(cycling_activity = cellcycle.kendall.repTime.weighted.median.cor.corrected) %>%
  dplyr::mutate(UID3 = case_when(UID2 == "UID-DLP-SA1044" ~ "T-47D",
                                 UID2 == "UID-DLP-SA928" ~ "SA928",
                                 TRUE ~ UID2))

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
                                     "UID-DLP-SA922" = "SA922",
                                     "UID-DLP-SA609" = "SA609",
                                     "UID-DLP-SA906" = "SA906",
                                     "UID-DLP-SA1044" = "SA1044",
                                     "UID-DLP-SA039" = "SA039",
                                     "UID-DLP-SA1087" = "SA1087",
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
                                     "UID-NNA-MDAMB231c28" =  "MDAMB231c28",
                                     "UID-NNA-MDAMB231c8" =  "mb231c8",
                                     "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))

# load full data
CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
valid_cells = colnames(CN)[colnames(CN) %in% df1$name]
pf = Biobase::protocolData(CN)@data %>% dplyr::filter(name %in% valid_cells)

## Analyse full data
colnames(pf) = paste0("pos.", colnames(pf))
pf1 = dplyr::inner_join(df1, pf, by=c("name"="pos.name"))

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
  geom_point(aes(x=UID2, y=cycling_mode), size=1.0, color="black") + 
  theme_pubclean()
q1 = ggpubr::ggarrange(q1a, q1b)
q1


df_test = qc_gini_norm(qc_mapd(qc_gini(predict_replicating(qc_rpc(df_ploidy2, cutoff_value = 2))  %>% dplyr::filter(!replicating, !rpc.outlier), cutoff_value = 2.0), cutoff_value = 2.0), cutoff_value = 2.0)
df_test$outlier = df_test$gini.outlier | df_test$dmapd.outlier | df_test$dgini.outlier
ggplot(data = df_test) +
  geom_quasirandom(aes(x=UID, y=gini, color=outlier)) + 
  theme_pubclean()
prop.table(table(df_test$outlier, df_test$UID), margin = 2)

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
df_ploidy = predict_replicating(df_ploidy2, cutoff_value = 2.0) %>% dplyr::filter(!replicating)
df_ploidy_cleaned = qc_error(qc_mapd(qc_gini(qc_gini_norm(df_ploidy)))) %>%
        dplyr::filter(!error.outlier, !dgini.outlier, !dmapd.outlier, !gini.outlier)



ggplot(data = df_ploidy) + geom_quasirandom(aes(x=UID2, y = ploidy, color=technology)) +
  theme_pubclean() + theme(axis.text.x = element_text(angle = 90)) 
df_ploidy = qc_ploidy(df_ploidy, range=0.5)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$technology),margin = 2)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$UID2),margin = 2)

# clean outlier cells
p_ploidy_prediction = ggplot(data = df_ploidy %>% dplyr::filter(rpc > 25) %>%
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
p_ploidy_prediction


ggplot(data = qc_error(qc_mapd(qc_gini(qc_gini_norm(df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA928"))))) %>%
        dplyr::filter(!error.outlier, !dgini.outlier, !dmapd.outlier, !gini.outlier)) +
  geom_point(aes(x=rpc, y=ploidy))

table(df_ploidy %>% dplyr::filter(rpc > 25) %>% dplyr::filter(startsWith(UID, "UID-DLP-SA906") | startsWith(UID, "UID-DLP-SA039")) %>% dplyr::pull(UID))
p_ploidy_prediction_p53 = ggplot(data = df_ploidy %>% dplyr::filter(rpc > 25) %>% dplyr::filter(startsWith(UID, "UID-DLP-SA906") | startsWith(UID, "UID-DLP-SA039"))) +
  geom_quasirandom(aes(x=SLX, y = ploidy, color=technology, shape=ploidy.outlier)) +
  geom_boxplot(aes(x=SLX, y=ploidy), color="black", outlier.shape = NA, alpha=0.2 ) +
  facet_wrap(~ UID2, scales = "free_x", strip.position="bottom", nrow=1) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_pubclean() + theme(legend.position="none") + xlab("cell line") + ylab("ploidy") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(axis.text.x = element_text(angle = 90))
p_ploidy_prediction_p53

outlier_SA928 = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA928", ploidy.outlier) %>% dplyr::pull(name)
p_outlier_SA928 = plotCopynumberHeatmap(CN[, outlier_SA928], cluster_rows = "ploidy")

# Relationship rpc and deviation
df_ploidy_outliers = df_ploidy %>% dplyr::group_by(UID, SLX) %>% dplyr::mutate(median_ploidy = median(ploidy)) %>% dplyr::ungroup() %>%
  dplyr::mutate(diff_ploidy = ploidy - median_ploidy)

p_ploidy_rpc = ggplot(data = df_ploidy_outliers %>% dplyr::filter(rpc > 25) %>% dplyr::filter(!startsWith(UID, "UID-NNA-T"), !startsWith(UID, "UID-10X-BR"))) + 
  geom_point(aes(x=used.reads, y=diff_ploidy, color=technology)) +
  facet_wrap(~technology+UID2)
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

#a = df_ploidy %>% dplyr::filter(startsWith(UID, "UID-NNA-TN5"), ploidy > 3.0) %>% dplyr::pull(name)
#plotCopynumber(CN[, a[5]], correction=TRUE)
#ggplot(data = df_ploidy %>% dplyr::filter(startsWith(UID, "UID-NNA-TN5"))) +
#  geom_quasirandom(aes(x=NA, y=error/rpc, color=ploidy)) + theme_pubclean()

#X_ploidy = apply(CN[, df_ploidy$name[df_ploidy$UID == "UID-10X-Andor-2020-49599-SNU-16"]]@assayData$copynumber, 2,
#  function(x){
#    y = x[startsWith(names(x), "X:")]
#    return(mean(y, na.rm=TRUE));
#  })
#hist(X_ploidy, breaks=100)


## G1/G2 prediction ====
threshold_rpc = 75
rpc_cutoff = 25 # minimum rpc to be included

df_ppredict = df1 %>% dplyr::filter(rpc > rpc_cutoff) %>%
  dplyr::filter(startsWith(UID, "UID-DLP") | startsWith(UID, "UID-NNA"), is.na(cellcycle) | cellcycle %in% c("G1", "G2", "S")) # | UID %in% c("UID-DLP-SA1044", "UID-DLP-SA928", "UID-DLP-SA1089", "UID-DLP-SA1090", "UID-DLP-SA1135", "UID-DLP-SA501X2XB00096"), is.na(cellcycle) | cellcycle %in% c("G1", "G2", "S"))
# NOTE: USE 125 for DLP data - different read size TODO no longer necessary -> use fragment size instead
# df_ppredict$prediction = df_ppredict$prediction.125
p1 = ggplot(data = df_ppredict) + geom_quasirandom(aes(x=UID, y=ploidy, color=UID)) +
  theme_pubclean() + theme(axis.text.x = element_text(angle = 90))
p1

p1 = ggplot(data = df_ppredict) + geom_quasirandom(aes(x=UID, y=prediction.quality, color=rpc)) +
  theme_pubclean() + theme(axis.text.x = element_text(angle = 90))
p1

df_predict = predict_replicating(df_ppredict) %>% dplyr::filter((SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")) | 
  (!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")) & !replicating)) %>%
  dplyr::filter(rpc > rpc_cutoff)
print(dim(df_predict))

#df_predict = qc_ploidy(qc_predict(df_ppredict, cutoff_percentile = 0.10), range=0.3) %>% dplyr::filter(startsWith(UID, "UID-DLP"), !ploidy.outlier, !predict.outlier_percentile)
#p2 = ggplot(data = df_predict %>% dplyr::filter(rpc > rpc_cutoff)) + geom_quasirandom(aes(x=SLX, y=prediction.quality, color=UID)) +
#  theme_pubclean()
#p2 
index_autosomes = !(startsWith(rownames(CN), "X:") | startsWith(rownames(CN), "Y:"))
valid = binsToUseInternal(CN)
normalness_ploidy = apply(CN[index_autosomes, df_predict$name]@assayData$copynumber, 2, function(x)
  sum(x == 2, na.rm=TRUE))
df_predict$normalness_ploidy = normalness_ploidy
df_predict$normal_cell = df_predict$normalness_ploidy > 0.9 * sum(index_autosomes & valid)
#save.image(file="~/Data/debug-ploidy.RData")

#ggplot(data = df_predict) + 
#  geom_histogram(aes(x=normalness_ploidy), bins=100) + 
#  geom_vline(xintercept=4255, color="red") + theme_pubclean()

#diploid_datasets = df_predict %>% dplyr::mutate(diploidness = (normalness_ploidy > 0.9 * sum(index_autosomes & valid))) %>%
#  dplyr::group_by(UID, SLX) %>% dplyr::summarize(mean_diploidness = mean(diploidness)) %>% dplyr::ungroup()

#p1 = ggplot(data = df_predict %>% dplyr::filter(rpc>rpc_cutoff) %>%
#              dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")))) +
#  geom_quasirandom(aes(x=UID, y=prediction.quality, color=UID)) +
#  facet_wrap(~ normal_cell)
#  theme_pubclean()
#  p1

df_predict$prediction.value = df_predict$prediction#ifelse(df_predict$normal_cell, df_predict$prediction.diploid, df_predict$prediction)

ggplot(data = df_predict %>% dplyr::filter(!startsWith(UID, "UID-NNA"),!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")))) +
  geom_quasirandom(aes(x=SLX, y=ploidy, color=UID2))
# Model 
p3 = ggplot(data=df_predict %>% dplyr::filter(rpc > rpc_cutoff) %>% 
             dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")))) + #, UID != "UID-DLP-SA928")) +
  geom_point(aes(x=rpc, y=prediction.value, color=startsWith(UID, "UID-NNA"), ploidy=ploidy), alpha=0.6, size=0.9) +
  geom_point(aes(x=200, y=1.25), shape=4, color="red") +
  geom_vline(xintercept = 80, color="black", linetype="dashed") +
  geom_smooth(aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,1)) +
  geom_density_2d(aes(x=rpc, y=prediction)) +
  geom_abline(aes(slope=slope, intercept=intercept), data=dplyr::tibble(slope=c(0.0015, 0.0015), intercept=c(1.10, 1.15), normal_cell=c(FALSE, TRUE)), color="black", linetype="dashed") +
  facet_wrap(~normal_cell) +
  theme_pubclean() + ylab("read density") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
  ylim(c(0, 2))
p3

# cleanup step
df_clean = df_predict %>% dplyr::filter(!startsWith(UID, "UID-NNA")) %>%
  dplyr::filter(rpc > rpc_cutoff) %>% #, ploidy < 2.25, ploidy > 1.75) %>%
  dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C"))) %>%
  dplyr::filter((((prediction.value <= (0.0015 * rpc) + 1.10)) & !normal_cell) | 
                (((prediction.value <= (0.0015 * rpc) + 1.15)) & normal_cell))

p3_post = ggplot(data=df_predict %>% dplyr::filter(!startsWith(UID, "UID-NNA")) %>%
             dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")))) + #, UID != "UID-DLP-SA928")) +
  geom_point(aes(x=rpc, y=prediction.value, color=UID, ploidy=ploidy), alpha=0.6, size=0.9) +
  geom_point(aes(x=200, y=1.25), shape=4, color="red") +
  geom_point(data=df_predict %>% dplyr::filter(!startsWith(UID, "UID-NNA")) %>% dplyr::filter(!(name %in% df_clean$name), !(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C"))),
             aes(x=rpc, y=prediction.value, ploidy=ploidy), color="grey", alpha=0.6, size=0.9) +
  geom_vline(xintercept = 80, color="black", linetype="dashed") +
  geom_smooth(aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
  #geom_smooth(aes(x=rpc, y=prediction), method=loess) + #, formula=y ~ poly(x,1)) +
  geom_density_2d(aes(x=rpc, y=prediction)) +
  geom_abline(aes(slope=slope, intercept=intercept), data=dplyr::tibble(slope=c(0.0015, 0.0015), intercept=c(1.10, 1.15), normal_cell=c(FALSE, TRUE)), color="black", linetype="dashed") +
  #geom_point(data = df_predict %>% dplyr::filter(cellcycle == "G1") %>%
  #           dplyr::filter((SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C"))), aes(x=rpc, y=prediction.value), color="red") +
  facet_wrap(~normal_cell) +
  theme_pubclean() + ylab("read density") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0)))
p3_post


ggplot(data = df_clean) +
  geom_point(aes(x=rpc, y=prediction.value)) +
  geom_density_2d(aes(x=rpc, y=prediction.value)) +
  geom_point(aes(x=100, y=1.2), color="red", shape=4) +
  geom_point(aes(x=200, y=1.3), color="red", shape=4) +
  facet_wrap(~normal_cell) + theme_pubclean()

#ggplot(data = df_predict %>% dplyr::filter((SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")), cellcycle=="G1")) +
#  geom_point(aes(x=rpc, y=ploidy))

model_all = robustbase::lmrob(prediction ~ poly(rpc,2), data=df_clean)
#model_all = loess(prediction ~ rpc, data=df_clean) #robustbase::lmrob(prediction ~ poly(rpc,1), data=df_clean)
summary(model_all)

df_model_diploid = df_clean %>% dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")), normal_cell)
#model_diploid = loess(prediction ~ rpc, data=df_model_diploid)
model_diploid = robustbase::lmrob(prediction ~ poly(rpc,2), data=df_model_diploid)
summary(model_diploid)

df_model = df_clean %>% dplyr::filter(!(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")), !(UID %in% c("UID-DLP-SA928", "UID-DLP-SA1089")), !normal_cell, ploidy < 2.5)
#model = loess(prediction ~ rpc, data=df_model) 
model = robustbase::lmrob(prediction ~ poly(rpc,2), data=df_model)
summary(model)

ggplot(data = df_model %>% dplyr::filter(UID %in% c("UID-DLP-SA1135", "UID-DLP-SA1090"))) + 
  geom_point(aes(x=rpc, y=prediction, color=UID))

df_prediction_diploid = df_predict %>% dplyr::filter(SLX %in% c("SLX-A73044A", "SLX-A90553C"))
df_prediction_diploid$predictor = predict(model_all, newdata = df_prediction_diploid)
#df_prediction_diploid$predictor = predict(model_diploid, newdata = df_prediction_diploid)

df_prediction = df_predict %>% dplyr::filter(SLX %in% c("SLX-A96139A"))
df_prediction$predictor = predict(model_all, newdata = df_prediction)
#df_prediction$predictor = predict(model, newdata = df_prediction)

df_prediction_all = df_predict %>% dplyr::filter(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C"))
df_prediction_all$predictor = predict(model_all, newdata = df_prediction_all)

p_diploid = ggplot(data=df_prediction_diploid) +
  geom_point(aes(x=rpc, y=prediction, color=cellcycle, shape=SLX), alpha=0.8, size=1.1) +
  geom_smooth(data=df_clean, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
  geom_point(data=df_clean, aes(x=rpc, y=prediction, ploidy=ploidy), color="black", alpha=0.1, size=1.1) +
  geom_vline(xintercept = 75, color="black", linetype="dashed") +
#  coord_cartesian(ylim=c(0, 80)) +
#  facet_wrap(~UID+SLX) +
  theme_pubclean() + ylab("read density") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 8, face = 2)) +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0)),
         shape = guide_legend(title = "Library", override.aes = list(size=3)))
p_diploid


p_other = ggplot(data=df_prediction) +
  geom_point(aes(x=rpc, y=prediction, ploidy=ploidy, color=cellcycle),  alpha=0.6, size=0.9) +
  geom_smooth(data=df_clean, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
  geom_point(data=df_clean, aes(x=rpc, y=prediction, ploidy=ploidy), color="black", alpha=0.1, size=1.1) +
#  geom_point(aes(x=rpc, y=predictor, ploidy=ploidy), color="cyan", alpha=0.6, size=0.9) +
  geom_vline(xintercept = 75, color="black", linetype="dashed") +
#  geom_smooth(aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x, 2)) +
#  facet_wrap(~UID+SLX) +
  theme_pubclean() + ylab("read density") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0)))
p_other


## prediction

df_prediction_all$delta = df_prediction_all$predictor - df_prediction_all$prediction
p_predict_nice = ggplot(data = df_prediction_all) + # %>% dplyr::filter(!replicating)) +
  geom_point(aes(x=rpc, y=delta, color=cellcycle), size=1.1) +
  geom_vline(xintercept = 75, linetype="dashed") +
  #geom_abline(slope = 0, intercept = -0.05, color="red", linetype="dotted") +
  geom_abline(data = dplyr::tibble(slope=c(0, 0, 0), intercept=c(-0.05, -0.04, -0.02),
                                   UID3=c("T-47D", "SA928", "SA928"), SLX=c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")), aes(slope=slope, intercept=intercept), color="red", linetype="dotted") +
  facet_wrap(~UID3+SLX) + #, scales="free_x") +
  ylab(expression(Delta ~ " prediction")) +
  theme_pubclean(base_size=20) +
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 8, face = 2)) +
    theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside") +
#      legend.position = "right") +
  guides(color = guide_legend(title="Ploidy (Cell cycle)    ", override.aes = list(fill="white", size = 3, alpha=1.0))) +
  scale_color_manual(values=c("G1"="#E69F00", "G2"="#56B4E9", "S"="#99999940"),
                     labels=c("P=2.7 (G1)", "P=5.4 (G2)", "P=2.7-5.4 (S)")) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank())
p_predict_nice


p_predict_nice_woS = ggplot(data = df_prediction_all %>% dplyr::filter(!replicating)) +
  geom_point(aes(x=rpc, y=delta, color=cellcycle), size=1.1) +
  geom_vline(xintercept = 75, linetype="dashed") +
  #geom_abline(slope = 0, intercept = -0.05, color="red", linetype="dotted") +
  geom_abline(data = dplyr::tibble(slope=c(0, 0, 0), intercept=c(-0.05, -0.04, -0.02), UID2=c("SA1044", "SA928", "SA928"), SLX=c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C")), aes(slope=slope, intercept=intercept), color="red", linetype="dotted") +
  facet_wrap(~UID2+SLX, scales="free_x") +
  ylab(expression(Delta ~ " prediction")) +
  theme_pubclean() +
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 8, face = 2)) +
    theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(fill="white", size = 3, alpha=1.0))) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))
p_predict_nice_woS



df_prediction$delta = df_prediction$predictor - df_prediction$prediction
p_other_predict = ggplot(data = df_prediction) + # %>% dplyr::filter(!replicating)) +
  geom_point(aes(x=rpc, y=delta, color=cellcycle), size=1.1) +
  geom_vline(xintercept = 75, linetype="dashed") +
  geom_abline(slope = 0, intercept = -0.05, color="red", linetype="dotted") +
  ylab(expression(Delta ~ " prediction")) +
  theme_pubclean() +
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 8, face = 2)) +
    theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))
p_other_predict

df_prediction_diploid$delta = df_prediction_diploid$predictor - df_prediction_diploid$prediction
p_diploid_predict = ggplot(data = df_prediction_diploid) + #%>% dplyr::filter(!replicating)) +
  geom_point(aes(x=rpc, y=delta, color=cellcycle)) +
  geom_abline(data = dplyr::tibble(slope=c(0, 0), intercept=c(-0.04, -0.02), SLX=c("SLX-A73044A", "SLX-A90553C")), aes(slope=slope, intercept=intercept), color="red", linetype="dotted") +
  facet_wrap(~SLX) + theme_pubclean() +
  geom_vline(xintercept = 75, linetype="dashed") + 
  ylab(expression(Delta ~ " prediction")) +
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(#family = "Playfair",
                                    #color = "chocolate",
                                    size = 8, face = 2)) +
    theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Library", breaks = NULL, labels = NULL)) +
  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999"))
p_diploid_predict

df_prediction$cellcycle_prediction = ifelse(df_prediction$rpc > 75, ifelse(df_prediction$delta < -0.05, "G2", "G1"), "NA")
df_prediction_diploid = df_prediction_diploid %>%
  dplyr::mutate(cellcycle_prediction = case_when(SLX == "SLX-A73044A" & delta < -0.04 & rpc > threshold_rpc ~ "G2", 
                                                 SLX == "SLX-A73044A" & delta >= -0.04 & rpc > threshold_rpc ~ "G1", 
                                                 SLX == "SLX-A90553C" & delta < -0.02 & rpc > threshold_rpc ~ "G2",
                                                 SLX == "SLX-A90553C" & delta >= -0.02 & rpc > threshold_rpc ~ "G1",
                                                 TRUE ~ "NA"
                                                 ))

# SA1044
table(df_prediction$cellcycle, df_prediction$cellcycle_prediction)
(191 + 128) / (191 + 4 + 20 + 128)

# SA928
table(df_prediction_diploid$cellcycle, df_prediction_diploid$cellcycle_prediction, df_prediction_diploid$SLX)
# SLX-A73044A
(266 + 119) / (266 + 65 + 7 + 119)
# SLX-A90553C
(144 + 115) / (144 + 2 + 47 + 115)

p3_nice = ggplot(data=predict_replicating(df_predict) %>% dplyr::filter(!replicating) %>% dplyr::filter(SLX %in% c("SLX-A96139A", "SLX-A73044A", "SLX-A90553C"))) +
  geom_point(aes(x=rpc, y=prediction, color=cellcycle), alpha=0.4, size=0.8) +
  geom_vline(xintercept = 80, color="black", linetype="dashed") +
  facet_wrap(~UID3+SLX, scales="free_x") +
  theme_pubclean() + ylab("read density") +
  guides(color = guide_legend(title="Ploidy (Cell cycle)  ", override.aes = list(fill="white", size = 3, alpha=1.0))) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank()) +
  scale_color_manual(values=c("G1"="#E69F00", "G2"="#56B4E9", "S"="#999999"),
                     labels=c("P=2.7 (G1)", "P=5.4 (G2)", "P=2.7-5.4 (S)"))
p3_nice


  guides(color = guide_legend(ncol=1,title="Data set", override.aes = list(size = 3, alpha=1.0)),
         fill = guide_legend(ncol=1, title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
  scale_fill_manual(values=c("G1"="#E69F00", "G2"="#56B4E9", "S"="#999999")) +

# Identify ploidy in DLP-SA1135
ggplot(data = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1135")) +
  geom_quasirandom(aes(x=UID, y=ploidy))

cells = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1135") %>% dplyr::pull(name)
length(cells)
nrow(df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1135", ploidy < 2.1, ploidy > 1.9))
77 / 360
p_sup_heatmap_SA1135 = plotCopynumberHeatmap(CN[, colnames(CN) %in% cells], cluster_rows = "ploidy", show_cell_names = FALSE)

cells = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1090") %>% dplyr::pull(name)
length(cells)
nrow(df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1090", ploidy > 2.5))
72 / 415
p_sup_heatmap_SA1090 = plotCopynumberHeatmap(CN[, colnames(CN) %in% cells], cluster_rows = "ploidy")

# This shows issue with identifiability (but rpc to low to fix it)
cells = df_ploidy %>% dplyr::filter(UID == "UID-NNA-TN4") %>% dplyr::pull(name)
length(cells)
nrow(df_ploidy %>% dplyr::filter(UID == "UID-NNA-TN4", ploidy < 2.1))
200 / 1204
p_sup_heatmap_TN4 = plotCopynumberHeatmap(CN[, colnames(CN) %in% cells], cluster_rows = "ploidy")

# Demonstrate SA1044 is not separable based on copy number profiles
cells = predict_replicating(df_ppredict) %>% dplyr::filter(UID == "UID-DLP-SA1044", cellcycle %in% c("G1", "G2"), !replicating) %>% dplyr::select(name, cellcycle)
ha = rowAnnotation(cellcycle = cells$cellcycle, col = list(QC = c("G1" = "#63ACBE", "G2" = "#EE442F")))
p_sup_heatmap_SA1044 = plotCopynumberHeatmap(CN[, colnames(CN) %in% cells$name], cluster_rows = "ploidy", har=ha, show_cell_names = FALSE)
p_sup_heatmap_SA1044

#ggplot(data = df_ploidy %>% dplyr::filter(UID == "UID-DLP-SA1135")) +
#  geom_point(aes(x=rpc, y=prediction.125, color=ploidy < 2.5))

df_ploidy_prediction_DLP = predict_replicating(df_predict) %>% dplyr::filter(technology %in% c("DLP")) %>% 
  dplyr::filter(cellcycle %in% c("G1") | is.na(cellcycle), !replicating, rpc > rpc_cutoff)
table(df_ploidy_prediction_DLP$UID)

## ACT excurs - problem we don't reach saturation of read density at all
#df_ploidy_prediction_ACT = predict_replicating(df1) %>% dplyr::filter(technology %in% c("ACT"), !(UID %in% c("UID-NNA-TN6-paired", "UID-NNA-TN7-paired"))) %>%
#  dplyr::filter(!replicating) %>% 
#  dplyr::mutate(correct_ploidy = dplyr::case_when((UID == "UID-NNA-MDAMB231c28" & ploidy > 3.0) ~ FALSE,
#                                                  (UID == "UID-NNA-TN4" & ploidy < 3.0) ~ FALSE,
#                                                  TRUE ~ TRUE))
#table(df_ploidy_prediction_ACT$UID)
#
#p_ploidy_ACT = ggplot() +
#  geom_point(data=df_ploidy_prediction_ACT, aes(x=rpc, y=prediction, color=UID3, shape=correct_ploidy, ploidy=ploidy), alpha=0.6, size=0.9) +
#  geom_smooth(data=df_ploidy_prediction_ACT, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
#  geom_vline(data=dplyr::tibble(x= 75, technology="DLP"), aes(xintercept=x), color="black", linetype="dashed") +
#  theme_pubclean() + ylab("read density") +
#  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
#      strip.background = element_rect(colour="white", fill="white"),
#      strip.placement = "outside",
#      legend.position = "right") +
#  guides(color = guide_legend(ncol=1,title="Data set", override.aes = list(size = 3, alpha=1.0)),
#         fill = guide_legend(ncol=1, title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
#  scale_fill_manual(values=c("G1"="#E69F00", "G2"="#56B4E9", "S"="#999999")) +
#  theme(legend.background=element_blank()) + theme(legend.key=element_blank())
#p_ploidy_ACT
### end ACT excurs

p_ploidy_prediction2 = ggplot() +
  geom_point(data=df_ploidy_prediction_DLP, aes(x=rpc, y=prediction, color=UID3, ploidy=ploidy), alpha=0.7, size=0.5) +
  geom_smooth(data=df_ploidy_prediction_DLP, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
  #geom_point(data=dplyr::tibble(x=c(1,1,1), y=c(1,1,1), state=as.factor(c("G1", "S", "G2"))), aes(x=x, y=y, fill=state), shape=21) +
#  geom_smooth(data=df_clean, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2), color="cyan") +
#  geom_point(data=df_ploidy_prediction_ACT, aes(x=rpc, y=prediction, color=interaction(technology, UID2, sep=" - "), ploidy=ploidy), alpha=0.6, size=0.9) +
#  geom_smooth(data=df_ploidy_prediction_ACT, aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x,2)) +
#  geom_point(data=qc_ploidy(df_ploidy_prediction_ACT) %>% dplyr::filter(UID == "UID-NNA-TN8", ploidy.outlier), aes(x=rpc, y=prediction, ploidy=ploidy), color="black", alpha=0.6, size=0.9) +
#  geom_point(aes(x=rpc, y=predictor, ploidy=ploidy), color="cyan", alpha=0.6, size=0.9) +
  geom_density2d(data=df_ploidy_prediction_DLP, aes(x=rpc, y=prediction), color="black") +
  geom_vline(data=dplyr::tibble(x= 75, technology="DLP"), aes(xintercept=x), color="black", linetype="dashed") +
#  geom_smooth(aes(x=rpc, y=prediction), method=robustbase::lmrob, formula=y ~ poly(x, 2)) +
  #facet_wrap(~technology, scales = "free_y") +
  theme_pubclean(base_size=20) + ylab("read density") +
  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
      strip.background = element_rect(colour="white", fill="white"),
      strip.placement = "outside",
      legend.position = "right") +
  guides(color = guide_legend(ncol=1,title="Data set", override.aes = list(size = 3, alpha=1.0)),
         fill = guide_legend(ncol=1, title="Cell cycle", override.aes = list(size = 3, alpha=1.0))) +
  scale_fill_manual(values=c("G1"="#E69F00", "G2"="#56B4E9", "S"="#999999")) +
  theme(legend.background=element_blank()) + theme(legend.key=element_blank())
p_ploidy_prediction2

model_DLP = robustbase::lmrob(prediction ~ poly(rpc, 2), data=df_ploidy_prediction_DLP)
summary(model_DLP)
saveRDS(model_DLP, file = "~/scAbsolute/data/ploidyPrediction/model_DLP.RDS")

#df_ploidy_prediction_DLP$predictor = predict(model_DLP, newdata = df_ploidy_prediction_DLP)
#df_ploidy_prediction_DLP$delta = df_ploidy_prediction_DLP$predictor - df_ploidy_prediction_DLP$prediction
#p_other_predict = ggplot(data = df_ploidy_prediction_DLP) +
#  geom_point(aes(x=rpc, y=delta, color=UID), size=1.1) +
#  geom_vline(xintercept = 75, linetype="dashed") +
#  geom_abline(slope = 0, intercept = -0.05, color="red", linetype="dotted") +
#  ylab(expression(Delta ~ " prediction")) +
#  theme_pubclean() +
#  theme(legend.key = element_rect(fill = "white"),
#        legend.title = element_text(#family = "Playfair",
#                                    #color = "chocolate",
#                                    size = 8, face = 2)) +
#    theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
#      strip.background = element_rect(colour="white", fill="white"),
#      strip.placement = "outside",
#      legend.position = "right") +
#  guides(color = guide_legend(title="Cell cycle", override.aes = list(size = 3, alpha=1.0)))
#p_other_predict



#model_ACT = robustbase::lmrob(prediction ~ poly(rpc, 2), data=df_ploidy_prediction_ACT)
#summary(model_ACT)
#saveRDS(model_ACT, file = "~/scAbsolute/data/ploidyPrediction/model_ACT.RDS")

# show it doesn't work for other technologies
#p5 = ggplot(data = df1 %>% dplyr::filter(technology %in% c("ACT"))) + #c("10X", "JBL"))) +
#  geom_point(aes(x=rpc, y=prediction, color=cellcycle), alpha=0.8, size=0.9) +
#  facet_wrap(~technology, scales = "free", ncol=1) +
#  theme_pubclean() + ylab("read density") +
#  theme(#axis.text.x = element_text(angle = 45, vjust = 0.0, hjust=0.0),
#      strip.background = element_rect(colour="white", fill="white"),
#      strip.placement = "outside",
#      legend.position = "none")
#p5


## Create panel figure
leg1 = get_legend(p3_nice)
leg2 = get_legend(p_ploidy_prediction2)
Fig_ploidy_II = ggpubr::ggarrange(
  p_ploidy_prediction2 + theme(legend.position="right"),
  p_predict_nice + theme(legend.position="top"),
  #p3_nice + theme(legend.position = "none", 
  #                strip.background = element_blank(), #element_rect(colour="white", fill="white"),
  #                strip.text.x = element_blank()),
  nrow=2, labels=c("A", "B"), heights = c(1.2,2)) #, hjust = 0.5, vjust=-0.5)
Fig_ploidy_II


ggpubr::ggexport(Fig_ploidy_II, filename = "~/scAbsolute/figures/Fig_ploidy_II.pdf", width=12, height=14)

prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$technology),margin = 2)
prop.table(table(df_ploidy$ploidy.outlier, df_ploidy$UID2),margin = 2)

Sup_ploidy = ggpubr::ggarrange(p_cellline_reads)
Sup_ploidy
ggpubr::ggexport(Sup_ploidy, filename = "~/scAbsolute/figures/Sup_ploidy.pdf", width=16, height=8)

ggpubr::ggexport(p_sup_heatmap_SA1044, filename = "~/scAbsolute/figures/Sup_ploidy_profiles_SA1044.pdf", width=16, height=8)

#ggpubr::ggexport(p_sup_heatmap_SA1135, filename = "~/scAbsolute/figures/Sup_ploidy_profiles_SA1135.pdf", width=16, height=8)

df_select_TN4 = df_predict %>% dplyr::filter(UID == "UID-NNA-TN4") %>% dplyr::select(name, ploidy) %>% dplyr::mutate(pout = ploidy < 2.5) %>%
  dplyr::arrange(desc(ploidy))
Sup_ploidy_TN4 = plotCopynumberHeatmap(CN[, df_select_TN4$name], row_split = as.factor(df_select_TN4$pout))
Sup_ploidy_TN4

ggpubr::ggexport(Sup_ploidy_TN4, filename = "~/scAbsolute/figures/Sup_ploidy_TN4.pdf", width=12, height=16)

# predict ploidy for cycling cells
df_ploidy_replicating = predict_replicating(df_ploidy2, cutoff_value = 2.0) %>% dplyr::filter(replicating)

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

ggpubr::ggexport(p_predict_nice_woS, filename = "~/scAbsolute/figures/Sup_ploidy_prediction_no_replicating.pdf", width=12, height=16)

#cn_profiles_TN4 = CN[, df_select_TN4$name]
#valid = binsToUseInternal(CN)
#dmat = cn_profiles_TN4@assayData$copynumber[valid, ]
#dmat[, df_select_TN4$pout] = dmat[, df_select_TN4$pout] * 2
## normalize to max state
#dmat_norm = ifelse(dmat > 8, 8, dmat)
#dim(dmat_norm)
#
#dist = copyNumberDistance(t(dmat_norm), BASEDIR="~/")
