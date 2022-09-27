## Output images for wetlab

## Variables ####
args = commandArgs(trailingOnly=TRUE)
set.seed(2020)

if (interactive()){
  rm(list=ls())
  
  sampleFile = c(
    "~/Data/project1/30-scale/500/predict/UID-JBL-NA12878_500.rds"
  )
  RESULTPATH = "~/Data/Ania/SLX-20749"
  
  # sample subset of data for debugging purposes
  sampleFile = do.call("c", (base::strsplit(sampleFile, split=",")))
  
}else{
  
  sampleFile = args[1]
  RESULTPATH = args[2]
  
  # sample subset of data for debugging purposes
  sampleFile = do.call("c", (base::strsplit(sampleFile, split=",")))
}

## Setup ====
BASEDIR="~/"
BASEDIR=normalizePath(BASEDIR)
require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE)
require(ggbeeswarm, quietly = TRUE, warn.conflicts = FALSE)
require(ggpubr, quietly = TRUE, warn.conflicts = FALSE)
require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
require(ComplexHeatmap, quietly = TRUE, warn.conflicts = FALSE)
source(file.path(BASEDIR, "scUnique/R/core.R"))
source(file.path(BASEDIR, "scUnique/R/visualize.R"))
source(file.path(BASEDIR, "scAbsolute/R/core.R"))
source(file.path(BASEDIR, "scAbsolute/R/visualization.R"))

if(length(sampleFile) == 1){
  object = readRDS(sampleFile) 
}else{
  object = combineQDNASets(future.apply::future_lapply(sampleFile, readRDS))
}

cutoff_replicating = 2.0
df = predict_replicating(Biobase::pData(object) %>%
 tidyr::separate(name, sep="_", into=c("UID", "SLX", "cellid", "celltag"), remove=FALSE), cutoff_value=cutoff_replicating)
df = df %>% dplyr::filter(SLX %in% c("SLX-20749"))
object = object[, df$name]
#annotation = readr::read_csv("~/mean-variance-model/data/cellcycle-info-PROTOCOL.csv") %>% dplyr::mutate(name = gsub("\\.bam", "", filename))
#df = df %>% dplyr::left_join(annotation, by="name", suffix=c("", ".y"))


FILENAME = paste0(lapply(sampleFile, function(x) sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(x))), collapse = "+")
if(!dir.exists(dirname(file.path(RESULTPATH, FILENAME, "a")))) dir.create(dirname(file.path(RESULTPATH, FILENAME, "a")), recursive = TRUE)
file.path(RESULTPATH, FILENAME)
for (celln in colnames(object)){
  p = plotCopynumber(object[, celln], correction = TRUE, showMarker = FALSE)
  ggpubr::ggexport(filename=file.path(RESULTPATH, FILENAME, paste0(celln, "_profile.png")), p, width = 480 * 2)
}

pa = plotCopynumberHeatmap(object, cluster_rows = FALSE, abbreviate_cell_names = TRUE, show_cell_names = TRUE)
ggpubr::ggexport(pa, filename=file.path(RESULTPATH, FILENAME, paste0(FILENAME, "_plot-heatmap.pdf")))
