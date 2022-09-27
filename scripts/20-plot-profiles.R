# plot cells individually

## Variables ####
args = commandArgs(trailingOnly=TRUE)

if (interactive()){
  rm(list=ls()) 
  
  cn = readRDS("~/UID-JBD-PROTOCOL_500.rds") 
  RESULTPATH = "~/ExampleImages"
  
  reticulate::use_condaenv(condaenv = "rstudio-server", conda = "~/.anaconda3/bin/conda")
  BASEDIR="~/scAbsolute"
}else{
  
  cn = args[1]
  RESULTPATH = args[2]
  
  BASEDIR="/opt/scAbsolute"
}

# PACKAGE DEPENDENCIES
library(reticulate, quietly=TRUE, warn.conflicts = FALSE)
library(future.apply, quietly=TRUE, warn.conflicts = FALSE)
suppressWarnings(plan(multiprocess))
library(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)

library(tidyverse, quietly=TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly=TRUE, warn.conflicts = FALSE)

library(Biobase, quietly=TRUE, warn.conflicts = FALSE)
library(BiocGenerics, quietly=TRUE, warn.conflicts = FALSE)
library(devtools, quietly=TRUE, warn.conflicts = FALSE)
library(digest, quietly=TRUE, warn.conflicts = FALSE)
library(IRanges, quietly=TRUE, warn.conflicts = FALSE)
library(MASS, quietly=TRUE, warn.conflicts = FALSE)
library(robustbase, quietly=TRUE, warn.conflicts = FALSE)
library(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
library(matrixStats, quietly=TRUE, warn.conflicts = FALSE)
print(paste0("Base directory: ", BASEDIR))

source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
source(file.path(BASEDIR, "R/scSegment.R"))
source(file.path(BASEDIR, "R/scAbsolute.R"))
source(file.path(BASEDIR, "R/core.R"))
source(file.path(BASEDIR, "R/mean-variance.R"))
source(file.path(BASEDIR, "R/visualization.R"))
source(file.path(BASEDIR, "R/cellcycle.R"))
# END DEPENDENCIES

library(ggpubr)
ifelse(!dir.exists(file.path(RESULTPATH)), dir.create(file.path(RESULTPATH)), FALSE)
for(i in 1:ncol(cn)){
  print(colnames(cn)[i])
  p1 = plotCopynumber(cn[, i])
  p2 = plotCopynumber(cn[, i], correction = TRUE, showMarker = FALSE)
  wrong=prop.table(table(cn[,i]@assayData$copynumber != 2))["TRUE"]
  print(table(cn[,i]@assayData$copynumber))
  
  p = ggpubr::ggarrange(p1 + ggpubr::rremove('xlab') + ggpubr::rremove("y.title"),
                        p2 + ggtitle(label=paste0("Percentage != 2: ", round(wrong, digits=4))) + ggpubr::rremove("y.title"), ncol = 1)
  ggpubr::ggexport(p, filename=file.path(RESULTPATH, paste0(colnames(cn)[i], ".png")))
  print("-----")
}
