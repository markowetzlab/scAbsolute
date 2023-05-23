# meta function to load dependencies for interactive analysis/vignette of scAbsolute package
load_dependencies <- function(BASEDIR="~/scAbsolute/") {
  
  # Load libraries
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
  # Load source files
  source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
  source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
  source(file.path(BASEDIR, "R/scSegment.R"))
  source(file.path(BASEDIR, "R/scAbsolute.R"))
  source(file.path(BASEDIR, "R/core.R"))
  source(file.path(BASEDIR, "R/mean-variance.R"))
  source(file.path(BASEDIR, "R/visualization.R"))
  source(file.path(BASEDIR, "R/cellcycle.R"))
}
