
args = commandArgs(trailingOnly=TRUE)
set.seed(2020)
require(QDNAseq)
require(Biobase)

path = args[1]
# path = "~/Data/project1/30-scale/30/scaling/UID-10X-Fibroblast-nuclei_30.rds"
cn = readRDS(path)

if(is.null(Biobase::assayDataElement(cn, "probdloss"))){
  corrected_counts = Biobase::assayDataElement(cn, "calls")
}else{
  corrected_counts = Biobase::assayDataElement(cn, "probdloss")  
}

valid = !is.na(corrected_counts[,1])

basepath = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(path))
write.table(t(corrected_counts[valid,]), file = file.path(dirname(path), paste0(basepath, ".csv")),
            sep = ",", row.names=FALSE, col.names=FALSE)

