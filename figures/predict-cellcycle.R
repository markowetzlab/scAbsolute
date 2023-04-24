# Copyright 2022, Michael Schneider, All rights reserved.
## Predict cellcycle on cellines - on DLP + JBL data

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
           "UID-10X-Fibroblast-cell", 
           "UID-10X-Fibroblast-nuclei",
           "UID-DLP-SA928",
           "UID-DLP-SA928-cellcycle",
           "UID-DLP-SA1044",
           "UID-DLP-SA922",
           "UID-DLP-SA1090",
           "UID-NNA-TN2",
           "UID-NNA-mb453",
           "UID-10X-Andor-2020-48959-MKN-45",
           "UID-10X-Andor-2020-50941-SNU-638", 
           "UID-10X-Andor-2020-49599-SNU-16", 
           "UID-10X-Andor-2020-49600-SNU-668", 
           "UID-10X-Andor-2020-49606-NCI-N87", 
           "UID-10X-Andor-2020-59068-KATOIII", 
           "UID-10X-Andor-2020-49607-HGC-27",
           "UID-10X-Andor-2020-59065-NUGC-4",
           "UID-10X-Andor-2020-59069-SNU-601",
           "UID-FML-PEO1-PLOIDY-2N",
           "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%",
           "UID-FML-PEO1-FUCCI"
         ), "_", binSizeValue, ".rds")))
pd = dplyr::bind_rows(lapply(files,
                             function(x) Biobase::pData(readRDS(x)) )) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

# modify UID names for facet plots
custom_label_function <- function(x, y){
  if(!is.data.frame(x)){
    z = dplyr::tibble(x)
  }else{
    z=x
  }
  
  if("UID" %in% colnames(x)){
    y = z %>% dplyr::transmute(dplyr::recode(UID, 
                                             "UID-10X-Andor-2020-50941-SNU-638" = "10X SNU-638", 
                                             "UID-10X-Fibroblast" = "10X Fibroblast", 
                                             "UID-DLP-SA1044" = "DLP T-47D", 
                                             "UID-DLP-SA928" = "DLP NA12878", 
                                             "UID-JBL-NA12878" = "JBL NA12878", 
                                             "UID-JBL-PEO1" = "JBL PEO1", 
                                             "UID-NNA-mb453" = "ACT mb453", 
                                             "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))
  }else{
    if("SLX" %in% colnames(x)){
      y = z %>% dplyr::transmute(dplyr::recode(SLX, 
                                               "SLX-00000" = "10X FB cell", 
                                               "SLX-00001" = "10X FB nuclei", 
                                               "SLX-18431" = "JBL NA12878",
                                               "SLX-20749" = "JBL NA12878",
                                               "SLX-A73044A" = "DLP A73044A",
                                               "SLX-A90553C" = "DLP A90553C",
                                               "SLX-A90689B" = "DLP A90689B",
                                               "SLX-A90689C" = "DLP A90689C",
                                               "SLX-A96139A" = "DLP T-47D",
                                               "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))
    }else{
      y = z
    }
  }
  
  if(!("UID" %in% colnames(x)) && ("SLX" %in% colnames(x))){
    y = z %>% dplyr::transmute(dplyr::recode(SLX, 
                                             "SLX-00000" = "10X FB cell", 
                                             "SLX-00001" = "10X FB nuclei", 
                                             "SLX-18431" = "JBL NA12878",
                                             "SLX-A73044A" = "DLP A73044A",
                                             "SLX-A90553C" = "DLP A90553C",
                                             "SLX-A90689B" = "DLP A90689B",
                                             "SLX-A90689C" = "DLP A90689C",
                                             "SLX-A96139A" = "DLP T-47D",
                                             "UID-10X-Fibroblast-cell" = "UID-10X-Fibroblast"))
  }
  
  if(!is.data.frame(x)){
    y = y[,1]
  }
  return(y)
}

# load cellcycle info
df_cellcycle_DLP = readr::read_csv("~/mean-variance-model/data/cellcycle-info-DLP.csv", col_types = "ccccdd")
df_cellcycle_JBL = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv", col_types = "ccccdd")
df_meta1 = dplyr::bind_rows(df_cellcycle_DLP, df_cellcycle_JBL) %>%
  dplyr::mutate(name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

df_meta2 = pd %>% dplyr::filter(UID %in% c("UID-FML-PEO1-FUCCI", "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%")) %>% 
  tidyr::separate(celltag, sep="-", into=c("A", "B", "cellcycle"), remove=FALSE, extra="drop") %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(cellcycle = case_when(UID == "UID-FML-PEO1-PLOIDY-2N" ~ "G1", 
                                      TRUE ~ as.character(cellcycle)))
df_meta = dplyr::bind_rows(df_meta1, df_meta2)

# SLX-18431 is early stage NA12878
df = pd %>% dplyr::filter(SLX=="SLX-00000" | SLX %in% c("SLX-18430", "SLX-18432", "SLX-20765", "SLX-20749") |
                            SLX %in% c("SLX-A96139A") | #DLP-SA1044
                            SLX %in% c("SLX-23300", "SLX-23301", "SLX-23302", "SLX-23303") | # FML-PEO1-FUCCI and PLOIDY
                            SLX %in% c("SLX-A90554B", "SLX-A96213A") | 
                            SLX %in% c("SLX-A73044A", "SLX-A90553C", "SLX-A90689B", "SLX-A90689C", "SLX-A90694B") | #DLP-SA928
                            SLX %in% c("SLX-A73046B", "SLX-A73056B", "SLX-A90554C", "SLX-A90600C", "SLX-A95660A", "SLX-A96146A", "SLX-A96172B", "SLX-A96211C") | #DLP-SA906
                            SLX %in% c("SLX-A73044B", "SLX-A73047D", "SLX-A96199B", "SLX-A96226B") | #DLP-SA039
                            SLX %in% c("SLX-NAVINACT")) %>%
  dplyr::left_join(df_meta, by="name", suffix=c("", ".y"))

table(df$cellcycle, df$UID)

# NOTE
## FIX READ COVERAGE
# we fix uneven read coverage in PEO1-FUCCI samples by combining original G1 with DOWNSAMPLED S, G2 cells
# it works for uneven read coverage, and below we demonstrate it works for even, too
ggplot(data = df %>% dplyr::filter(UID %in% c("UID-FML-PEO1-FUCCI", "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%"))) + 
  facet_wrap(~SLX+UID) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=replicating))

ggplot(data = df %>% dplyr::filter(UID %in% c("UID-FML-PEO1-FUCCI", "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%"))) + 
  facet_wrap(~SLX+UID) +
  geom_quasirandom(aes(x=cellcycle, y=used.reads, color=replicating))

df = df %>% dplyr::mutate(valid = dplyr::case_when(UID == "UID-FML-PEO1-FUCCI" & cellcycle == "G1" ~ TRUE,
                                       UID == "UID-FML-PEO1-FUCCI" & cellcycle %in% c("G2", "S") ~ FALSE,
                                       UID == "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%" & cellcycle == "G1" ~ FALSE,
                                       UID == "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%" & cellcycle %in% c("S", "G2") ~ TRUE,
                                       TRUE ~ TRUE)) %>%
  dplyr::filter(valid) %>%
  dplyr::mutate(UID = ifelse(UID == "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%", "UID-FML-PEO1-FUCCI", UID))

## END FIX


# load full data
CN = combineQDNASets(lapply(files,function(x) readRDS(x)), dropProtocol = TRUE, reduceMetadata = TRUE)
valid_cells = colnames(CN)[colnames(CN) %in% df$name]


## 1. Predict cycling cells

## diploid cells should have fixed number of CN2
diploid_baseline = CN[!(startsWith(rownames(CN), "X:") | startsWith(rownames(CN), "Y:")), colnames(CN)[colnames(CN) %in% (df %>% dplyr::filter(UID %in% c("UID-DLP-SA928", "UID-10X-Fibroblast-cell", "UID-10X-Fibroblast-nuclei", "UID-JBL-NA12878")) %>% dplyr::pull(name))]]
ppl = apply(diploid_baseline@assayData$copynumber, 2, function(x) return(sum(x==2, na.rm=TRUE)))
#max(table(ppl)) -> 2369 is 4727 mode of distribution
CN_quality = dplyr::tibble(filename = paste0(names(ppl), ".bam"), nc2=ppl) %>%
  dplyr::mutate(high_quality = nc2 >= 4711 & nc2 <= 4729, quality_gradient = abs(nc2-4727))

df_A = df %>% dplyr::mutate(filename=ifelse(is.na(filename), paste0(name, ".bam"), filename)) %>%
  dplyr::left_join(CN_quality, by=c("filename"="filename"), suffix=c("", ".q2"))

## Show off cellcycle working
measure_name = "Cycling activity"


## NOTE: For unbalanced data, we compute the cutoff value on D or G1 data, otherwise bias in highly unbalanced data
df_cutoff = predict_replicating(df_A %>%
  dplyr::filter((startsWith(UID, "UID-FML") & cellcycle %in% c("G1")) | 
                (startsWith(UID, "UID-DLP") & cellcycle %in% c("G1", "D")) ),
  cutoff_value=1.5, iqr_value = 1.5) %>%
  dplyr::group_by(UID, SLX) %>%
  dplyr::summarise(cycling_cutoff_subset = unique(cycling_cutoff_sd),
                   cycling_cutoff_iqr_subset = unique(cycling_cutoff_iqr),
                   cycling_median_subset = unique(cycling_median))

df = predict_replicating(dplyr::left_join(df_A, df_cutoff, by=c("SLX", "UID"))) %>%
  dplyr::mutate(replicating = ifelse(!is.na(cellcycle), cycling_activity > cycling_cutoff_subset | cellcycle.cmi_yrT > cycling_cutoff_iqr_subset, replicating)) %>%
  dplyr::mutate(cycling_median_override = ifelse(!is.na(cycling_median_subset), cycling_median_subset, cycling_median),
                cycling_cutoff_override = ifelse(!is.na(cycling_cutoff_subset), cycling_cutoff_subset, cycling_cutoff_sd))

# filtered out -> only valid if we use both
#ggplot(data = df %>% dplyr::filter(UID == "UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%")) +
#  geom_hline(aes(yintercept=cycling_cutoff_override)) +
#  geom_vline(aes(xintercept=cycling_cutoff_iqr_subset)) +
#  facet_wrap(~SLX) +
#  geom_point(aes(x=cellcycle.cmi_yrT, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=replicating))


## compare mean vs median
p0 = ggplot(data=df %>% 
  #dplyr::filter(nc2 > 1000, SLX != "SLX-A90694B", UID %in% c("UID-10X-Fibroblast-cell", "UID-DLP-SA928", "UID-JBL-NA12878"))) +
  dplyr::filter(UID %in% c("UID-DLP-SA1044", "UID-DLP-SA928", "UID-FML-PEO1-FUCCI"))) +
  #geom_point(aes(x=cellcycle.cmi_yrT, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=as.factor(high_quality))) +
  geom_point(aes(x=cellcycle.cmi_yrT, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=cellcycle)) +
  facet_wrap(~UID+SLX, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") +
  ylab("Kendall") + xlab("MUTINFO") + labs(color="CN state")
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p0


#  normal (diploid) cells without CNAs
p1 = ggplot(data=df %>% dplyr::mutate(SLX = case_when(UID == "UID-10X-Fibroblast-nuclei" ~ "SLX-00001", TRUE ~ SLX)) %>%
              dplyr::filter(nc2 > 1000, !(SLX %in% c("SLX-A96146A","SLX-A96172B","SLX-A96199B","SLX-A96211C","SLX-A96226B")),
                            SLX != "SLX-A90694B", UID %in% c("UID-10X-Fibroblast-cell", "UID-10X-Fibroblast-nuclei", "UID-DLP-SA928", "UID-JBL-NA12878"))) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=as.factor(high_quality), name=name, alpha=0.9)) +
  geom_hline(aes(x=cellcycle, yintercept=cycling_cutoff_override)) +
  facet_wrap(~SLX, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") + guides(alpha="none") +
  ylab(measure_name) + xlab("Cellcycle status") + labs(color="CN state")
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p1

# G1 and S phase only
data_cutoff = df %>% dplyr::mutate(SLX3 = ifelse(UID=="UID-FML-PEO1-FUCCI", paste0("DLPi ", SLX), as.character(SLX))) %>%
                 dplyr::mutate(SLX2 = factor(SLX, levels=c("SLX-A73044A", "SLX-A90553C", "SLX-A96139A",
                                                           "SLX-23300", "SLX-23301", "SLX-23302"), 
                                             labels = c("DLP A73044A", "DLP A90553C", "DLP T47D", "DLPi 23300", "DLPi 23301", "DLPi 23302"), ordered=TRUE)) %>%
                 dplyr::filter((nc2 > 1000 &
                                 UID %in% c("UID-DLP-SA928")) |
                                 UID == "UID-DLP-SA1044" | UID == "UID-FML-PEO1-FUCCI",
                               SLX != "SLX-A90694B", cellcycle %in% c("G1", "S"))

p2 = ggplot(data=data_cutoff) +
  geom_density(data=data_cutoff %>% dplyr::filter(cellcycle == "G1"),
               aes(x=cycling_activity, group=cellcycle, linetype="G1", color="G1"), bins=100) +
  geom_density(data=data_cutoff %>% dplyr::filter(cellcycle == "S"),
               aes(x=cycling_activity, group=cellcycle, linetype="S", color="S"), bins=100) +
  geom_vline(aes(type="Median", xintercept=cycling_median_override, SLX=SLX, color="Median", linetype="Median"), show.legend = TRUE) + 
  geom_vline(aes(type="Cutoff", xintercept=cycling_cutoff_override, SLX=SLX, color="Cutoff", linetype="Cutoff"), show.legend=TRUE) + 
  facet_wrap(~SLX2, labeller = custom_label_function, nrow=1, scales="free_x") +
  theme_pubclean(base_size=20) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_colour_manual(name="Cellcycle annotation", values = c("G1" = "black", "S" = "black", "Median" = "blue", "Cutoff" = "red"),
                      labels=c("Threshold", "G1", "Median", "S")) +
  scale_linetype_manual(name="Cellcycle annotation", values = c("G1" = "solid", "S" = "dotted", "Median" = "dashed", "Cutoff" = "dotdash"),
                      labels=c("Threshold", "G1", "Median", "S")) +
  xlab("Cycling activity") + ylab("Density") +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "right",
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 16, face = 1)) +
theme(#legend.text=element_text(size=12), legend.title=element_text(size=14),
      legend.key.height  = grid::unit(0.15, "npc"),
      legend.key.size=unit(0.8, "cm"), 
      legend.key = element_rect(fill = NA, color = NA),
      legend.background=element_rect(colour="white")) +
  guides(color=guide_legend(ncol=2))
  #guides(color = guide_legend(override.aes = list(linetype = c(1, 2, 2, 1) ))) 
p2


# Extend to cell lines
p3 = ggplot(data=df %>% dplyr::filter(!(SLX %in% c("SLX-A90694B", "SLX-A90689B", "SLX-A90689C", "SLX-A73044A")), UID %in% c("UID-DLP-SA1044", "UID-DLP-SA928", "UID-JBL-PEO1", "UID-JBL-NA12878"))) +
  #geom_hline(yintercept=0.5, color="red", linetype="dashed") + 
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.cmi_yrT, color=as.factor(high_quality))) +
  #geom_boxplot(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected)) +
  facet_wrap(~UID, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") +
  ylab(measure_name) + xlab("Cellcycle status") + labs(color="Uniform\nCN state") + 
  theme_pubclean()
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p3

p3b = ggplot(data=df %>% dplyr::filter(!(SLX %in% c("SLX-A90694B", "SLX-A90689B", "SLX-A90689C", "SLX-A73044A")), UID %in% c("UID-DLP-SA1044", "UID-DLP-SA928", "UID-JBL-PEO1", "UID-JBL-NA12878"))) +
  #geom_hline(yintercept=0.5, color="red", linetype="dashed") + 
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.cmi_yrT, color=replicating)) +
  geom_boxplot(aes(x=cellcycle, y=cellcycle.cmi_yrT),  alpha=0.0) +
  geom_hline(aes(yintercept=cycling_cutoff_iqr_subset), color="black", linetype="dashed") +
  facet_wrap(~UID, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") +
  ylab(measure_name) + xlab("Cellcycle status") + labs(color="Uniform\nCN state") + 
  theme_pubclean()
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p3b


p4 = ggplot(data=df %>% dplyr::filter(!(SLX %in% c("SLX-A90694B", "SLX-A90689B", "SLX-A90689C", "SLX-A73044A")), UID %in% c("UID-10X-Andor-2020-50941-SNU-638", "UID-NNA-mb453"))) +
  #geom_hline(yintercept=0.5, color="red", linetype="dashed") + 
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.cmi_yrT, color=replicating)) +
  facet_wrap(~UID, nrow=1, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") +
  ylab(measure_name) + xlab("Cellcycle status") + labs(color="Uniform\nCN state") + 
  theme_pubclean()# + coord_cartesian(ylim=c(-0.2, 0.3))
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p4

require(caret)
require(kernlab)
require(e1071)

#df_A_predict = df_A[rowSums(is.na(df_A)) > 0,] %>% dplyr::filter(startsWith(UID, "UID-DLP") | startsWith(UID, "UID-JBL")) %>% dplyr::filter(cellcycle %in% c("G1",  "S")) %>% dplyr::select(!starts_with(c("hmm2"))) %>%
#  dplyr::filter(!is.na(cellcycle.kendall.repTime.mean)) %>% dplyr::filter(!is.na(cellcycle.kendall.gcContent.sd.tau)) %>%
#  mutate(across(where(is_character),as_factor)) %>% dplyr::select(-min_ploidy, -max_ploidy, -runtime, -scaling.rpc_robust, -scaling.rpc_95, -cellcycle.kendall.gcContent.weighted.mean.tau, -gc_correction.alpha, -n_reads, -scaling.rpc_p95) %>%
#  select_if(~ !any(is.na(.))) %>% select_if(~ (!any(is.logical(.)) & !any(is.character(.))) & !(any(is.factor(.)) & (length(unique(.)) == 1))) %>%
#  dplyr::filter(UID != "UID-JBL-OVCAR3")
## dplyr::select(starts_with("cellcycle")) %>% 

#threshold_SA1044 = predict_replicating(df_A_predict %>% dplyr::filter(cellcycle == "D", UID=="UID-DLP-SA1044")) %>% dplyr::mutate(threshold = cycling_mode + 1 * cycling_sd_est) %>% dplyr::pull(threshold)
#threshold_SA928 = predict_replicating(df_A_predict %>% dplyr::filter(cellcycle == "D", UID=="UID-DLP-SA928")) %>% dplyr::mutate(threshold = cycling_mode + 1 * cycling_sd_est) %>% dplyr::pull(threshold)
#df_A_predict = df_A_predict %>% dplyr::filter(cellcycle != "D")
df_A_predict = df %>% dplyr::select(starts_with("cellcycle."), cellcycle, rpc, ploidy, hmm.alpha, UID, SLX, replicating, name) %>%
  dplyr::filter(!is.na(hmm.alpha)) %>%
  dplyr::filter(UID %in% c("UID-DLP-SA1044") |
                UID == "UID-DLP-SA928" & SLX %in% c("SLX-A73044A", "SLX-A90553C"))  %>%
  dplyr::select_if(~ !any(is.na(.)))



table(df_A_predict$cellcycle, df_A_predict$UID)
#UID-DLP-SA1044 UID-DLP-SA928 UID-FML-PEO1-FUCCI
#D             277           461                  0
#G1            240           818                286
#G2            386           532                288
#S             402           575                285

df_A_predict_cycling = df_A_predict %>% dplyr::filter(cellcycle %in% c("G1", "S")) %>%
  dplyr::select(-replicating, -name)

control <- trainControl(method='cv',
                        number=10)
tunegrid <- expand.grid(.mtry = c(3, 4, 5, 10))

f1 <- caret::train(cellcycle ~ ., data = df_A_predict_cycling,
                   metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
f1$finalModel
#randomForest::importance(f1$finalModel)
#a = data.frame(importance=randomForest::importance(f1$finalModel)) %>% dplyr::arrange(desc(MeanDecreaseGini))
#b = base::sort(a)

f2 <- caret::train(cellcycle ~ rpc + hmm.alpha + cellcycle.cmi_yrT + cellcycle.kendall.repTime.weighted.median.cor.corrected, data = df_A_predict_cycling,
                   metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
f2$finalModel

# NOTE: derive threshold from G1 cells, as other dataset is massively unbalanced
# okay to do, as we only look at lower part of distribution anyways and 
# mode should always be on G1
prediction_simple = df %>% 
  dplyr::filter(UID %in% c("UID-DLP-SA1044") |
                UID == "UID-DLP-SA928" & SLX %in% c("SLX-A73044A", "SLX-A90553C")) %>%
  dplyr::filter(cellcycle %in% c("G1", "S")) %>%
  dplyr::select(replicating, cellcycle, SLX, UID)


prediction_simple$cellcycle_prediction = ifelse(prediction_simple$replicating, "S", "G1")
f3 = caret::confusionMatrix(table(prediction_simple$cellcycle, prediction_simple$cellcycle_prediction))
f3

# simple model - diploid
prediction_simple_diploid = prediction_simple %>% dplyr::filter(UID == "UID-DLP-SA928")
f3_diploid = caret::confusionMatrix(table(prediction_simple_diploid$cellcycle, prediction_simple_diploid$cellcycle_prediction))

prediction_simple_1044 = prediction_simple %>% dplyr::filter(UID == "UID-DLP-SA1044")
f3_1044 = caret::confusionMatrix(table(prediction_simple_1044$cellcycle, prediction_simple_1044$cellcycle_prediction))

f3_diploid

f3_1044

# prediction on gold standard ground truth:
prediction_groundtruth = df %>% dplyr::filter(UID == "UID-FML-PEO1-FUCCI") %>%
  dplyr::filter(cellcycle %in% c("G1", "S")) %>%
  dplyr::mutate(cellcycle_prediction = ifelse(replicating, "S", "G1"))
f4_FUCCI = caret::confusionMatrix(table(prediction_groundtruth$cellcycle, prediction_groundtruth$cellcycle_prediction))

f4_FUCCI

#caret::confusionMatrix(table(prediction_simple, df_A_predict$cellcycle_prediction))

#f3 <- caret::train(cellcycle_prediction ~ rpc + hmm.alpha_estimate + cellcycle.kendall.repTime.weighted.median.cor.corrected + cellcycle.kendall.repTime.3q.cor.corrected, data = df_A_predict %>% dplyr::select(-cellcycle),
#                   metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
#f3$finalModel
#
#
## Underestimate of performance -> predict only on G1/S
#f4 <- caret::train(cellcycle ~ rpc + hmm.alpha_estimate + cellcycle.kendall.repTime.weighted.mean.tau.corrected + cellcycle.kendall.repTime.mean.tau.corrected + cellcycle.kendall.repTime.weighted.median.tau.corrected,
#                   data = df_A_predict %>% dplyr::filter(cellcycle %in% c("G1", "S")) %>% dplyr::select(-cellcycle_prediction),
#                   metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
#f4$finalModel

# Visualize issue with S/G2 in case of SA1044
## Skip diploid cells
#ggplot(data = df1) + geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.mean.tau.corrected, color=dmapd.outlier), alpha=0.3) + 
#  facet_wrap(~UID+SLX) + geom_hline(yintercept=0.5, linetype="dashed", color="red") + theme_pubclean()
#threshold = 0.5
#c1 = df2$name[df2$UID == "UID-DLP-SA928" & df2$cellcycle == "S" & df2$cellcycle.kendall.repTime.weighted.mean.tau.corrected < threshold]
#c2 = df2$name[df2$UID == "UID-DLP-SA928" & df2$cellcycle == "G1" & df2$cellcycle.kendall.repTime.weighted.mean.tau.corrected < threshold]
#c3 = df2$name[df2$UID == "UID-DLP-SA928" & df2$cellcycle == "G2" & df2$cellcycle.kendall.repTime.weighted.mean.tau.corrected < threshold]
#c4 = df2$name[df2$UID == "UID-DLP-SA928" & df2$cellcycle == "S" & df2$cellcycle.kendall.repTime.weighted.mean.tau.corrected > threshold]
#c5 = df2$name[df2$UID == "UID-DLP-SA928" & df2$cellcycle == "G2" & df2$cellcycle.kendall.repTime.weighted.mean.tau.corrected > threshold]
#ordered_samples = c(c1, c2, c3, c4, c5)
#ordered_groups = rep(c("S (b)", "G1 (b)", "G2 (b)", "S (a)", "G2 (a)"), times=c(length(c1), length(c2), length(c3), length(c4), length(c5)))
#subset = sample(1:length(ordered_samples), size=300, replace = FALSE)
#plotCopynumberHeatmap(CN[, ordered_samples[subset]], row_split = as.factor(ordered_groups[subset]))


# G2 and S phase overlap (example DLP-SA1044)
df_B = df_A_predict
df_B$cellcycle_interaction = interaction(df_B$cellcycle, df_B$replicating)
df_B$cellcycle_interaction = ifelse(df_B$cellcycle == "G1", "G1", df_B$cellcycle_interaction)
df_B$cellcycle_interaction = ifelse(df_B$cellcycle == "D", "D", df_B$cellcycle_interaction)
df_B$cellcycle_interaction = dplyr::recode_factor(df_B$cellcycle_interaction, `3`="G2", `4`="S (predicted G1)", `7`="G2 (predicted S)", `8`="S", `G1`="G1", `D`="D")
df_B$cellcycle_interaction = factor(ordered=TRUE, df_B$cellcycle_interaction, levels=c("G2", "S (predicted G1)", "G2 (predicted S)", "D", "S", "G1"))

df2 = df_B
c1 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "G1" & !df2$replicating]
c2 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "S" & df2$replicating]
c3 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "G2" & !df2$replicating]
c4 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "G2" & df2$replicating]
c5 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "S" & !df2$replicating]
c6 = df2$name[df2$UID == "UID-DLP-SA1044" & df2$cellcycle == "D"]

ordered_samples = c(c1, c2, c3, c4, c5)
ordered_groups = rep(c("G1", "S", "G2", "G2 (predicted S)", "S (predicted G1)"), times=c(length(c1), length(c2), length(c3), length(c4), length(c5)))
#ordered_samples = c(c1, c2, c5)
#ordered_groups = rep(c("G1", "S", "S (predicted G1)"), times=c(length(c1), length(c2), length(c5)))
index = startsWith(rownames(CN), "2:") | startsWith(rownames(CN), "3:")
p6 = plotCopynumberHeatmap(CN[index, ordered_samples], row_split = as.factor(ordered_groups),show_cell_names = FALSE)
p6

#ggpubr::ggexport(p6, path="~/test-cases.pdf")

valid = binsToUseInternal(CN)
set.seed(2023)
model = uwot::umap(t(Biobase::assayDataElement(CN[, ordered_samples], "copynumber")[valid,]), n_components = 2, ret_model = TRUE)
D_embedding = uwot::umap_transform(t(Biobase::assayDataElement(CN[, c6], "copynumber")[valid,]), model)
embedding = model$embedding

ordered_groups_recoded = factor(c(ordered_groups, rep("D", times=length(c6))), levels=c("G2", "S (predicted G1)", "G2 (predicted S)", "D", "S", "G1"))#levels(df_A$cellcycle_interaction))
p7 = ggplot(data = dplyr::tibble(x1=c(embedding[,1], D_embedding[,1]), x2=c(embedding[, 2], D_embedding[,2]), groups = ordered_groups_recoded)) +
  geom_point(aes(x=x1, y=x2, color=groups), size=1.0, alpha=0.3) + 
  geom_point(aes(x=-3, y=0), color="red", shape=3) +
  scale_color_brewer("Cellcycle ", palette="Set2") +
  facet_wrap(~groups) +
  #geom_density2d(aes(x=x1, y=x2, groups=groups), color="black", bins=100) +
  xlab("UMAP 1") + ylab("UMAP 2") + 
  theme_pubclean() +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(legend.key = element_rect(size=8, color = "white", fill = "white"),
        #legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(margin = margin(r = 0.3, unit = 'cm')),
        #legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 14, face = 1)) +
  guides(colour = guide_legend(nrow=1, override.aes = list(alpha = 1, size=2), title.hjust = 2))
p7


p5 = ggplot(data=df_B %>% dplyr::filter(UID == "UID-DLP-SA1044")) +
  #geom_density(aes(x=cellcycle.kendall.repTime.weighted.mean.tau.corrected, color=as.factor(high_quality))) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=cellcycle_interaction), bins=100) +
  #  geom_hline(yintercept=0.01, color="red", linetype="dashed") + 
  ylab(measure_name) + xlab("Cellcycle annotation") +
  scale_colour_brewer(palette = "Set2") +
  theme_pubclean()
p5


## Collate figures ====

# rework what goes together (figure for data with cell cycle annotation)

p_annotation = ggplot(data=df %>% dplyr::mutate(SLX = case_when(UID == "UID-10X-Fibroblast-nuclei" ~ "SLX-00001", TRUE ~ SLX)) %>%
                        dplyr::filter(SLX %in% c("SLX-A90553C","SLX-A73044A","SLX-A96139A"))) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=as.factor(high_quality), name=name, alpha=0.9)) +
  #geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.mean.cor.corrected, color=as.factor(high_quality))) +
  facet_wrap(~SLX, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") + guides(alpha="none") +
  ylab("Cycling\nactivity") + xlab("Cellcycle annotation (FACS)") + labs(color="CN state") +
  theme_pubclean(base_size=20) +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "right",
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 14, face = 1))
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p_annotation

df_seq = dplyr::bind_rows(df_A %>% dplyr::filter(startsWith(UID, "UID-10X-F") | startsWith(UID, "UID-DLP-SA928")), df %>% dplyr::filter(!(startsWith(UID, "UID-10X-F") | startsWith(UID, "UID-DLP-SA928")))) %>%
  dplyr::mutate(technology = dplyr::case_when(startsWith(UID, "UID-10X") ~ "10X",
                                              startsWith(UID, "UID-DLP") ~ "DLP",
                                              startsWith(UID, "UID-NNA") ~ "ACT",
                                              TRUE ~ "other"),
                UID3 = dplyr::case_when(UID=="UID-10X-Fibroblast-nuclei" ~ "Fibroblast",
                                        UID=="UID-10X-Andor-2020-50941-SNU-638" ~ "SNU-638",
                                        UID=="UID-NNA-mb453" ~ "mb453",
                                        UID=="UID-NNA-TN2" ~ "Patient TN2",
                                        UID=="UID-DLP-SA1090" ~ "OV2295",
                                        UID=="UID-DLP-SA928" ~ "NA12878",
                                        TRUE ~ UID))

p_sequencing = ggplot(data=df_seq %>% dplyr::filter( #!(SLX %in% c("SLX-A90694B", "SLX-A90689B", "SLX-A90689C", "SLX-A73044A")),
  UID %in% c("UID-10X-Fibroblast-nuclei", "UID-10X-Andor-2020-50941-SNU-638", "UID-NNA-mb453", "UID-NNA-TN2", "UID-DLP-SA1090") | SLX == "SLX-A90689C")) +
  #geom_hline(yintercept=0.5, color="red", linetype="dashed") + 
  geom_quasirandom(aes(x=UID3, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=as.factor(high_quality))) +
  theme_pubclean() + 
  facet_wrap(~technology, nrow=1, scales="free_x") + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") +
  ylab("Cycling\nactivity") + xlab("Sample") + labs(color="CN states") + 
  theme_pubclean(base_size=20) + #coord_cartesian(ylim=c(-0.2, 0.3)) +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  scale_color_discrete(labels=c("Varying", "Uniform", "not\navailable")) +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "right",
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 16, face = 1)) +
  guides(color = guide_legend(nrow=2, override.aes = list(size=3)))
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p_sequencing

p_groundtruth = ggplot(data=df %>% dplyr::filter(UID=="UID-FML-PEO1-FUCCI") %>%
                         dplyr::mutate(SLX2 = gsub("SLX-", "DLPi ", SLX))) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, name=name), alpha=0.9) +
  #geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.mean.cor.corrected, color=as.factor(high_quality))) +
  facet_wrap(~SLX2, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") + guides(alpha="none") +
  ylab("Cycling\nactivity") + xlab("Cellcycle annotation (FUCCI)") + labs(color="CN state") +
  theme_pubclean(base_size=20) +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "right",
        axis.title.y = element_blank(),
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 14, face = 1))
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
p_groundtruth


leg = get_legend(p_sequencing)
leg_lines = get_legend(p2)
Fig_cellcycle = ggpubr::ggarrange(ggpubr::ggarrange(
  ggpubr::ggarrange(leg, leg_lines, nrow=2, labels=c("", ""), heights=c(1,2)), p_sequencing + theme(legend.position = "none"), nrow =1, widths=c(2, 5), labels=c("", "A")),
                                  p2 + labs(x="Cycling activity") + theme(legend.position="none"),
                                  ggpubr::ggarrange(p_annotation + theme(legend.position="none"), p_groundtruth, nrow=1, widths=c(1,1), labels=c("C", "D")),
                                  nrow=3, labels=c("", "B", ""), heights = c(1.05, 1, 1))
Fig_cellcycle

Fig_cellcycle_sup = ggpubr::ggarrange(p5 + theme(legend.position = "none"),
                                      p7,
                                      ncol=1, nrow=2, labels=c("", ""), heights = c(2, 5))
Fig_cellcycle_sup


ggplot(data=df %>% dplyr::filter(UID=="UID-FML-PEO1-FUCCI")) +
  geom_quasirandom(aes(x=cellcycle, y=total.reads, name=name), alpha=0.9)

sup_groundtruth = ggplot(data=df %>% dplyr::filter(UID=="UID-FML-PEO1-FUCCI") %>%
                         dplyr::mutate(SLX2 = gsub("SLX-", "DLPi ", SLX))) +
  geom_quasirandom(aes(x=cellcycle, y=cellcycle.cmi_yrT, name=name), alpha=0.9) +
  #geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.mean.cor.corrected, color=as.factor(high_quality))) +
  facet_wrap(~SLX2, labeller=custom_label_function) + theme_pubclean() + #scale_color_continuous(name="# of CNAs", trans="log") +
  theme(legend.direction='vertical', legend.position = "right") + guides(alpha="none") +
  ylab("Cycling\nactivity") + xlab("Cellcycle annotation (FUCCI)") + labs(color="CN state") +
  theme_pubclean(base_size=20) +
  theme(panel.spacing = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white"),
        strip.placement = "outside") +
  theme(legend.key = element_rect(color = "white", fill="white"),
        legend.position = "right",
        axis.title.y = element_blank(),
        legend.title = element_text(#family = "Playfair",
          #color = "chocolate",
          size = 14, face = 1))
#scale_color_viridis_c(direction = -1, alpha=0.9, limits = c(1, 4000), na.value = 0, guide = guide_legend("title"))
sup_groundtruth

## Main figure - cell cycle - S phase
#leg = get_legend(p1)
#Fig_cellcycle = ggpubr::ggarrange(p1 + theme(legend.position="none"),
#                  ggpubr::ggarrange(ggpubr::ggarrange(leg, p2 + labs(x="Cycling activity"), nrow=1, widths = c(1, 5), labels=c("", "C")), 
#                                    p5 + theme(legend.position = "none"), ncol=1, nrow=2, labels=c("", "D"), heights = c(3, 2), vjust = -1),
#                  cowplot::plot_grid(p3 + theme(plot.margin = margin(t=6, r=0, b=9, l=0)), 
#                                     p4 + theme(axis.ticks.y=element_blank(), axis.title.y = element_blank(), axis.text.y=element_blank(),
#                                                axis.ticks.x=element_blank()) + 
#                                       theme(#plot.background = element_rect(color = 1,
#                                             #                                size = 1),
#                                               panel.spacing = unit(5, "pt"),
#                                               plot.margin = margin(t = 6,  # Top margin
#                                                                    r = 0,  # Right margin
#                                                                    b = 9,  # Bottom margin
#                                                                    l = 0)) + labs(y=NULL, x=""),
#                  nrow=1, rel_widths = c(2, 1)),
#                  p7,
#                  labels=c("A", "", "B",G "E"), nrow=2, ncol=2)
#Fig_cellcycle


ggpubr::ggexport(Fig_cellcycle, filename = "~/scAbsolute/figures/Fig_cellcycle.pdf", width=16, height=12)

# Supplementary figure - DLP-SA1044
ggpubr::ggexport(p6, filename="~/scAbsolute/figures/Sup_Cellcycle-SA1044_A.pdf")
ggpubr::ggexport(Fig_cellcycle_sup, filename="~/scAbsolute/figures/Sup_Cellcycle-SA1044_B.pdf")


## Compare direct to PERT model - sample 100 G1 cells, build model and evaluate on remaining subsets
dfcomp = Biobase::pData(readRDS("~/Data/project1/30-scale/500/predict/UID-FML-PEO1-FUCCI-DOWNSAMPLED-EQUAL-FORDLP_500.rds"))
membership = readr::read_tsv("~/scAbsolute/data/cellcycle/results-PERT/membership_200.tsv", col_types = c("cc"))

train_set = dfcomp %>% dplyr::filter(name %in% membership$Name[membership$Set == "Train"])
test_set_g1 = dfcomp %>% dplyr::filter(name %in% membership$Name[membership$Set == "Test_G1"])
test_set_s = dfcomp %>% dplyr::filter(name %in% membership$Name[membership$Set == "Test_S"])
train_set$UID = "dummy"
train_set$SLX = "dummy"
test_set_g1$true_cellcycle = "G1"
test_set_s$true_cellcycle = "S"
test_set = dplyr::bind_rows(test_set_g1, test_set_s)
test_set$cycling_activity = test_set$cellcycle.kendall.repTime.weighted.median.cor.corrected

train_set = predict_replicating(train_set, cutoff_value = 1.5, iqr_value = 1.5)
train_set$cycling_activity


test_set$predict_cellcycle = ifelse(test_set$cycling_activity > unique(train_set$cycling_cutoff_sd) | test_set$cellcycle.cmi_yrT > unique(train_set$cycling_cutoff_iqr),
                              "S", "G1")

#test_set$predict_cellcycle = ifelse(test_set$cycling_activity > unique(train_set$cycling_cutoff_sd),
#                              "S", "G1")

ff = caret::confusionMatrix(table(test_set$true_cellcycle, test_set$predict_cellcycle))

ff

ggplot(data = test_set) + 
  geom_quasirandom(aes(x=true_cellcycle, y=cycling_activity, color=predict_cellcycle)) +
  geom_hline(yintercept = unique(train_set$cycling_cutoff_sd), color="red")

ggplot(data = test_set) + 
  geom_quasirandom(aes(x=true_cellcycle, y=cellcycle.cmi_yrT, color=predict_cellcycle)) +
  geom_hline(yintercept = unique(train_set$cycling_cutoff_iqr), color="red")

#df %>% dplyr::filter(UID == "UID-FML-PEO1-FUCCI")


###  Gemini stained data ====
#df_gemini = df %>% dplyr::filter(UID %in% c("UID-FML-PEO1-FUCCI")) %>% dplyr::mutate(UID="UID-FML-PEO1", SLX="SLX-UNIFY")
## not enough cells per library and condition (imbalanced data sets)
#df2_gem = predict_replicating(df_gemini %>% dplyr::filter(cellcycle %in% c("G1", "G2")), cutoff_value = 1.0, hard_threshold = 0.1)
#df_gemini$cycling_activity = df_gemini$cellcycle.kendall.repTime.weighted.median.cor.corrected
#
#df_gemini$replicating = df_gemini$cycling_activity > (unique(df2_gem$cycling_mode) + (1.0 * unique(df2_gem$cycling_sd_est)))
##df_gemini$replicating = df_gemini$cycling_activity > 0.1 #(df_gemini$cycling_mode + (1.0 * df_gemini$cycling_sd_est))
#ggplot(data = df_gemini) +
#  geom_quasirandom(aes(x=cellcycle, y=cellcycle.cmi_yrT, color=replicating))
#
#ggplot(data = df_gemini) +
#  geom_quasirandom(aes(x=cellcycle, y=cellcycle.kendall.repTime.weighted.median.cor.corrected, color=replicating))
#
#table(df_gemini$cellcycle, df_gemini$replicating)
#
#d = table(df_gemini$cellcycle[df_gemini$cellcycle %in% c("G1", "S")], df_gemini$replicating[df_gemini$cellcycle %in% c("G1", "S")])
#sum(diag(d))/sum(d) #overall accuracy
#
### plot CN profiles
#g1_name = df_gemini[!df_gemini$replicating & df_gemini$cellcycle == "G1",] %>%
#  dplyr::select(name, ploidy) %>% dplyr::arrange(ploidy) %>% dplyr::pull(name)
#s_correct = df_gemini[df_gemini$replicating & df_gemini$cellcycle == "S",] %>%
#  dplyr::select(name, ploidy) %>% dplyr::arrange(ploidy) %>% dplyr::pull(name)
#s_incorrect = df_gemini[!df_gemini$replicating & df_gemini$cellcycle == "S",] %>%
#  dplyr::select(name, ploidy) %>% dplyr::arrange(ploidy) %>% dplyr::pull(name)
#g2_cycling = df_gemini[df_gemini$replicating & df_gemini$cellcycle == "G2",] %>%
#  dplyr::select(name, ploidy) %>% dplyr::arrange(ploidy) %>% dplyr::pull(name)
#
##plotCopynumberHeatmap(CN[, g1_name])
#
#ordered_samples = c(g1_name, s_correct, s_incorrect, g2_cycling)
#ordered_groups = rep(c("G1", "S", "S (G1?)", "G2 (S?)"), times=c(length(g1_name), length(s_correct), length(s_incorrect), length(g2_cycling)))
##index = startsWith(rownames(CN), "2:") | startsWith(rownames(CN), "3:")
#p_gemini_copynumber = plotCopynumberHeatmap(CN[, ordered_samples], row_split = as.factor(ordered_groups),show_cell_names = FALSE)
#p_gemini_copynumber
#
#
### Random forest - feature selection
#require(caret)
#require(kernlab)
#require(e1071)
#
#df_predict_gemini = df_gemini %>% dplyr::select(starts_with("cellcycle."), cellcycle) %>%
#  dplyr::filter(cellcycle %in% c("G1", "S"))
#
#control <- trainControl(method='cv',
#                        number=10)
#tunegrid <- expand.grid(.mtry = c(3, 4, 5, 10))
#
#f_gemini <- caret::train(cellcycle ~ ., data = df_predict_gemini,
#                         metric = 'Accuracy',  method = "rf", tunegrid = tunegrid)
#f_gemini$finalModel
#randomForest::importance(f_gemini$finalModel)
#a = data.frame(importance=randomForest::importance(f_gemini$finalModel)) %>% dplyr::arrange(desc(MeanDecreaseGini))
#b = base::sort(a)


#cellcycle.cmi_yrT
