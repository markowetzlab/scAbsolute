source(file.path("R", "core.R"))

stopifnot(identical(
  .normalizeChromosomesToSeqlevels(c("chr1", "chrX", "chrY_random"), c("1", "X", "Y")),
  c("1", "X", "Y_random")))
stopifnot(identical(
  .normalizeChromosomesToSeqlevels(c("1", "X", "Y"), c("chr1", "chrX", "chrY")),
  c("chr1", "chrX", "chrY")))
stopifnot(identical(
  .normalizeChromosomesToSeqlevels(c("1", "X", "Y"), c("1", "X", "Y")),
  c("1", "X", "Y")))

valid = c(as.character(1:19), "X", "Y")
bed_chromosomes = c("chr1", "chrX", "chrY", "chrUn_GL456239")
normalized = .normalizeChromosomesToSeqlevels(bed_chromosomes, valid)
stopifnot(identical(normalized[normalized %in% valid], c("1", "X", "Y")))

cat("mouse seqlevel normalization tests passed\n")
