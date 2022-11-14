## Example for copy number fits in cells undergoing S phase (SA1044, SA928)


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
           "UID-DLP-SA928",
           "UID-DLP-SA1044"
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
                        "UID-DLP-SA1044" = "DLP SA1044", 
                        "UID-DLP-SA928" = "DLP SA928", 
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
                        "SLX-A96139A" = "DLP SA1044",
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
df_meta = dplyr::bind_rows(df_cellcycle_DLP, df_cellcycle_JBL) %>%
  dplyr::mutate(name = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename))) %>%
  tidyr::separate(name, sep="_", into=c("UID2", "SLX", "cellid", "celltag"), remove=FALSE) %>%
  dplyr::mutate(UID = str_replace(UID2, "-D11-1-B4", ""))

# SLX-18431 is early stage NA12878
df = pd %>% dplyr::filter(SLX=="SLX-00000" | SLX %in% c("SLX-18430", "SLX-18432", "SLX-20765", "SLX-20749") |
                          SLX %in% c("SLX-A96139A") | #DLP-SA1044
                          SLX %in% c("SLX-A73044A", "SLX-A90553C", "SLX-A90689B", "SLX-A90689C", "SLX-A90694B") | #DLP-SA928
                          SLX %in% c("SLX-A73046B", "SLX-A73056B", "SLX-A90554C", "SLX-A90600C", "SLX-A95660A", "SLX-A96146A", "SLX-A96172B", "SLX-A96211C") | #DLP-SA906
                          SLX %in% c("SLX-A73044B", "SLX-A73047D", "SLX-A96199B", "SLX-A96226B") | #DLP-SA039
                          SLX %in% c("SLX-NAVINACT")) %>%
  dplyr::left_join(df_meta, by="name", suffix=c("", ".y"))


# load full data
CN = combineQDNASets(lapply(files,function(x) readRDS(x)))
valid_cells = colnames(CN)[colnames(CN) %in% df$name]

predict_replicating(df) %>% dplyr::filter(!replicating, cellcycle == "G1", UID=="UID-DLP-SA1044") %>% dplyr::pull(name)
predict_replicating(df) %>% dplyr::filter(replicating, cellcycle == "S", UID=="UID-DLP-SA928") %>% dplyr::pull(name)

predict_replicating(df) %>%
  dplyr::filter(replicating, cellcycle == "S", UID=="UID-DLP-SA1044", cycling_activity < 0.1, cycling_activity > 0.0) %>%
  dplyr::pull(name)

p = ggplot(data = predict_replicating(df)) +
  geom_quasirandom(aes(x=UID, y=cycling_activity, name=name))

p0 = plotCopynumber(CN[, "UID-DLP-SA928_SLX-A90553C_001571_R66-C09"],main = "",correction = FALSE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE)
p01 = plotCopynumber(CN[, "UID-DLP-SA928_SLX-A90553C_000037_R43-C39"],main = "",correction = FALSE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE)
p02 = plotCopynumber(CN[, "UID-DLP-SA928_SLX-A90553C_000038_R43-C40"],main = "",correction = FALSE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE)

p1 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000003_R03-C07"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE)
p2 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000014_R03-C31"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # > 0.2
p3 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000015_R03-C34"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # > 0.2

p4 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000048_R04-C25"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # <0.2 & > 0.1
p5 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000055_R04-C37"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # <0.2 & > 0.1

p6 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000047_R04-C24"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # <0.1 & > 0.0
p7 = plotCopynumber(CN[, "UID-DLP-SA1044_SLX-A96139A_000049_R04-C26"],main = "",correction = TRUE,verbose = FALSE,readinfo = FALSE, showUnique = FALSE, showMarker = FALSE) # <0.1 & > 0.0

Sup_cellcycle_example_replicating = ggpubr::ggarrange(
  ggpubr::ggarrange(p0 + rremove("x.title"), p01 + rremove("x.title"), p02, ncol=1, nrow=3),
  ggpubr::ggarrange(p1 + rremove("y.title") + rremove("x.title"),
                    p2 + rremove("y.title") + rremove("x.title"),
                    p4 + rremove("y.title") + rremove("x.title"),
                    p5 + rremove("y.title") + rremove("x.title") ,
                    p6 + rremove("y.title"), #+ rremove("x.title"),
                    p7 + rremove("y.title"), #+ rremove("x.title"),
                    ncol=2, nrow=3), ncol=2, widths = c(1, 2), labels=c("A", "B"))
Sup_cellcycle_example_replicating

ggpubr::ggexport(Sup_cellcycle_example_replicating, filename = "~/scAbsolute/figures/Sup_cellcycle_example_replicating.pdf", width=16, height=16)

