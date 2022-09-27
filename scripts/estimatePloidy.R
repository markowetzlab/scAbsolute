## estimate ploidy for CN object

args = commandArgs(trailingOnly=TRUE)
path = args[1]
set.seed(2020)
require(QDNAseq)
require(Biobase)

# path = "~/Data/project1/30-scale/30/scaling/UID-10X-Fibroblast-nuclei_30.rds"
cn = readRDS(path)
ploidies = apply(Biobase::assayDataElement(cn, "copynumber"), 2, mean, na.rm=TRUE)
ploidy = max(2, round(median(ploidies, na.rm=TRUE), digits=0))
cat(ploidy)
