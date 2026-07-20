source(file.path("R", "scSegment.R"))

stopifnot(is.na(.shortChromosomeChangepoints(numeric())))
stopifnot(identical(.shortChromosomeChangepoints(5), 1L))
stopifnot(is.null(.shortChromosomeChangepoints(c(5, 6))))

stopifnot(identical(formals(.selectHMMSummaryChromosomes)$sex, "auto"))
stopifnot(identical(formals(copynumberSegmentation)$sex, "auto"))

select_chromosomes <- function(chrom_names=c("1", "X", "Y"),
                               alpha=c(0.1, 0.2, 0.3),
                               rpc=c(10, 20, 30),
                               logprob=c(-10, -20, -30),
                               epoch1=c(100, 100, 100),
                               sex="auto"){
  n = length(chrom_names)
  .selectHMMSummaryChromosomes(
    chrom_names=chrom_names, alpha=alpha,
    alpha_estimate=rep(0.1, n), rpc=rpc,
    rpc_zero=rep(1, n), alpha_zero=rep(0.1, n),
    logprob=logprob, epoch1=epoch1,
    n_oscillations=rep(0, n), sd_oscillations=rep(0, n),
    magnitude_oscillations=rep(0, n), sex=sex, cell_name="test-cell")
}

female = select_chromosomes(sex="female")
stopifnot(identical(female$keep, c(TRUE, TRUE, FALSE)))
stopifnot(identical(female$dropped_chromosomes, "Y"))

male = select_chromosomes(sex="male")
stopifnot(all(male$keep))

auto = suppressWarnings(select_chromosomes(alpha=c(0.1, 0.2, NaN),
                                           logprob=c(-10, -20, -Inf), sex="auto"))
stopifnot(identical(auto$keep, c(TRUE, TRUE, FALSE)))

failed_autosome = suppressWarnings(select_chromosomes(
  alpha=c(NaN, 0.2, 0.3), logprob=c(-Inf, -20, -30), sex="female"))
stopifnot(identical(failed_autosome$keep, c(FALSE, TRUE, FALSE)))
stopifnot(identical(failed_autosome$dropped_chromosomes, c("1", "Y")))

failed_epoch = suppressWarnings(select_chromosomes(epoch1=c(0, 100, 100)))
stopifnot(identical(failed_epoch$keep, c(FALSE, TRUE, TRUE)))

all_invalid_error = tryCatch({
  suppressWarnings(select_chromosomes(alpha=rep(NaN, 3),
                                      logprob=rep(-Inf, 3),
                                      epoch1=rep(0, 3)))
  FALSE
}, error=function(e) grepl("no chromosomes have a valid HMM summary", conditionMessage(e)))
stopifnot(all_invalid_error)

invalid_sex_error = tryCatch({
  select_chromosomes(sex="unknown")
  FALSE
}, error=function(e) TRUE)
stopifnot(invalid_sex_error)

length_error = tryCatch({
  select_chromosomes(rpc=c(10, 20))
  FALSE
}, error=function(e) grepl("do not match", conditionMessage(e)))
stopifnot(length_error)

sex = c("female", "male")
global_sex_does_not_leak = select_chromosomes(sex="auto")
stopifnot(all(global_sex_does_not_leak$keep))

cat("sex-aware chrY tests passed\n")
