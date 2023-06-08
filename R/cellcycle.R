# Copyright 2022, Michael Schneider, All rights reserved.
#' cellcycleMetadata
#'
#' \code{cellcycleMetadata}
#'
#' @param segmentedCounts QDNAseq object with segmented slot
#' @param segment_size_cutoff minimum segment size to compute kendall test for
#'
#' @return a QDNAseq object with a set of quantitative measures related to a cell's cell cycle status 
#' (as part of the phenoData in segmentedCounts object)
#'
#' @export
cellcycleMetadata <- function(segmentedCounts,segment_size_cutoff=20){

  require(Kendall, quietly = TRUE, warn.conflicts = FALSE)
  require(trend, quietly = TRUE, warn.conflicts = FALSE)
  require(Biobase)
  require(S4Vectors)
  # NOTE not possibly to use limitPloidy as cutoff,
  # as we do not have a way to determine relative ploidy at this stage
  limitPloidy = Inf
  valid = binsToUseInternal(segmentedCounts)

  iterate_processing <- function(i){

    gr_raw = createGR(segmentedCounts[,i,drop=FALSE])

    if (species == "Human"){
      rep_rle = readRDS(file.path(BASEDIR, "data/replicationTiming/replicationInfo.RDS"))
    } else if (species == "Mouse"){
      rep_rle = readRDS(file.path(BASEDIR, "data/replicationTiming/replicationInfo_mouse.RDS"))
    } else {
      stop("Species is not supported")
    }
    # accommodate restricted genome regions, i.e. only selected chromosomes
    rep_rle = rep_rle[names(rep_rle) %in% seqlevels(gr_raw)]
    replicationTiming = binnedAverage(gr_raw, rep_rle, "replicationTime", na.rm=TRUE)

    segmentation = Biobase::assayDataElement(segmentedCounts[valid, i], "segmented")
    counts = Biobase::assayDataElement(segmentedCounts[valid, i], "calls")
    gc_map_covariate = estimate_gc_correction(segmentedCounts[valid,i,drop=FALSE])

    repTime = replicationTiming$replicationTime[valid]

    ## Approach based on different replication time for different subsets of the same segment ->
    # no longer optimal segmentation as can be observed ====
    rl = S4Vectors::Rle(segmentation)
    n_segs = length(rl@lengths)
    index = rep(1:n_segs, times = rl@lengths)
    results = lapply(1:n_segs, function(x){
      # print(x)
      if(sum(index == x) < segment_size_cutoff){
        return(list("pval"=NA, "tau"=NA, "tau_corrected"=NA, "cor_corrected"=NA))
      }
      co = counts[index == x]
      if(all(co == 0)){
        return(list("pval"=NA, "tau"=NA, "tau_corrected"=NA, "cor_corrected"=NA))
      }
      time = repTime[index == x]
      covar = gc_map_covariate[index == x]
      idx = base::sort(time, index.return=TRUE)$ix
      # ggplot(data=dplyr::tibble(co=co[idx], covar=covar[idx], time=time[idx])) + geom_point(aes(x=time, y=co, color=covar))
      trend.kendall_test = Kendall::MannKendall(co[idx])
      trend.kendall_test.corrrected = trend::partial.mk.test(co[idx], covar[idx])
      pv = trend.kendall_test[["sl"]][[1]]
      tau = trend.kendall_test[["tau"]][[1]]
      tau.corrected = trend.kendall_test.corrrected[["statistic"]]
      cor.corrected = trend::partial.cor.trend.test(co[idx], covar[idx], method = "spearman")[["estimate"]]
      return(list("pval"=pv, "tau"=tau, "tau_corrected"=tau.corrected, "cor_corrected"=cor.corrected))
    })

    pvals = unlist(sapply(results, "[", 1))
    tau = unlist(sapply(results, "[", 2))
    tau.corrected = unlist(sapply(results, "[", 3))
    cor.corrected = unlist(sapply(results, "[", 4))
    
    rll = rl@lengths[!is.na(pvals)]
    rll2 = rl@lengths[!is.na(cor.corrected)]
    cor.corrected = cor.corrected[!is.na(cor.corrected)]
    tau = tau[!is.na(pvals)]
    tau.corrected = tau.corrected[!is.na(pvals)]
    pvals = pvals[!is.na(pvals)]

    weighted.mean = stats::weighted.mean(pvals, rll)
    weighted.median = stats::median(rep(pvals, times=rll))
    weighted.sd = stats::sd(rep(pvals, times=rll))

    weighted.mean.tau = stats::weighted.mean(tau, rll)
    weighted.median.tau = stats::median(rep(tau, times=rll))
    weighted.sd.tau = stats::sd(rep(tau, times=rll))

    weighted.3q.tau.corrected = quantile(rep(tau.corrected, times=rll), c(0.75))
    weighted.1q.tau.corrected = quantile(rep(tau.corrected, times=rll), c(0.25))
    weighted.mean.tau.corrected = stats::weighted.mean(tau.corrected, rll)
    weighted.median.tau.corrected = stats::median(rep(tau.corrected, times=rll))
    weighted.sd.tau.corrected = stats::sd(rep(tau.corrected, times=rll))
    weighted.skewness.tau.corrected = spssSkew(rep(tau.corrected, times=rll))

    weighted.3q.cor.corrected = quantile(rep(cor.corrected, times=rll), c(0.75))
    weighted.1q.cor.corrected = quantile(rep(cor.corrected, times=rll), c(0.25))
    weighted.mean.cor.corrected = stats::weighted.mean(cor.corrected, rll2)
    weighted.median.cor.corrected = stats::median(rep(cor.corrected, times=rll2))
    weighted.sd.cor.corrected = stats::sd(rep(cor.corrected, times=rll2))
    weighted.skewness.cor.corrected = spssSkew(rep(cor.corrected, times=rll2))
    
    df_result = dplyr::tibble(
      kendall.repTime.mean = mean(pvals),
      kendall.repTime.median = stats::median(pvals),
      kendall.repTime.wmean = weighted.mean,
      kendall.repTime.wmedian = weighted.median,
      kendall.repTime.wsd = weighted.sd,
      kendall.repTime.sd = sd(pvals),
      kendall.repTime.weighted.mean.tau = weighted.mean.tau,
      kendall.repTime.weighted.median.tau = weighted.median.tau,
      kendall.repTime.weighted.sd.tau = weighted.sd.tau,
      kendall.repTime.3q.tau.corrected = weighted.3q.tau.corrected,
      kendall.repTime.1q.tau.corrected = weighted.1q.tau.corrected,
      kendall.repTime.weighted.mean.tau.corrected = weighted.mean.tau.corrected,
      kendall.repTime.weighted.median.tau.corrected = weighted.median.tau.corrected,
      kendall.repTime.weighted.sd.tau.corrected = weighted.sd.tau.corrected,
      kendall.repTime.weighted.skewness.tau.corrected = weighted.skewness.tau.corrected,
      kendall.repTime.mean.tau = mean(tau),
      kendall.repTime.median.tau = median(tau),
      kendall.repTime.mean.tau.corrected = mean(tau.corrected),
      kendall.repTime.median.tau.corrected = median(tau.corrected),
      kendall.repTime.sd.tau = sd(tau),
      kendall.repTime.weighted.mean.cor.corrected = weighted.mean.cor.corrected,
      kendall.repTime.weighted.median.cor.corrected = weighted.median.cor.corrected,
      kendall.repTime.weighted.sd.cor.corrected = weighted.sd.cor.corrected,
      kendall.repTime.weighted.skewness.cor.corrected = weighted.skewness.cor.corrected,
      kendall.repTime.3q.cor.corrected = weighted.3q.cor.corrected,
      kendall.repTime.1q.cor.corrected = weighted.1q.cor.corrected
    )

    return(df_result)
  }

  dfs = lapply(1:dim(segmentedCounts)[[2]], iterate_processing)

  df = do.call("rbind", dfs)
  colnames(df) = paste0("cellcycle.", colnames(df))
  description <- pData(segmentedCounts)

  # necessary in order not to loose row names
  cnames = rownames(description)
  description = dplyr::bind_cols(description, df)
  rownames(description) = cnames
  pData(segmentedCounts) = description

  return(segmentedCounts)
}


cellcycle_gctest <- function(segmentedCounts, segment_size_cutoff=20){
  
  require(Kendall, quietly = TRUE, warn.conflicts = FALSE)
  stopifnot("segmented" %in% Biobase::assayDataElementNames(segmentedCounts))
  
  iterate_processing <- function(i){
  
    object = segmentedCounts[,i]
    valid = (!(base::startsWith(rownames(object), "X:") | base::startsWith(rownames(object), "Y:"))) &
      binsToUseInternal(object) & (Biobase::assayDataElement(object, "calls") >= 1) &
      (Biobase::assayDataElement(object, "segmented") >= 1e-6)

    counts = Biobase::assayDataElement(object[valid, 1], "calls")
    
    # This should eliminate any effect
    # counts = counts * (1/estimate_gc_correction(object[valid, 1]))
    
    segmentation = Biobase::assayDataElement(object[valid, 1], "segmented")
    gc = Biobase::fData(object[valid, 1])[["gc"]]
    rl = S4Vectors::Rle(segmentation)
    n_segs = length(rl@lengths)
    index = rep(1:n_segs, times = rl@lengths)
    results = lapply(1:n_segs, function(x){
      # print(x)
      if(sum(index == x) < segment_size_cutoff){
        return(list("pval"=NA, "tau"=NA))
      }
      co = counts[index == x]
      time = gc[index == x]
      idx = base::sort(time, index.return=TRUE)$ix
      trend.kendall_test = Kendall::MannKendall(co[idx])
      pv = trend.kendall_test[["sl"]][[1]]
      tau = trend.kendall_test[["tau"]][[1]]
      return(list("pval"=pv, "tau"=tau))
    })
    
    pvals = unlist(sapply(results, "[", 1))
    tau = unlist(sapply(results, "[", 2))
    rll = rl@lengths[!is.na(pvals)]
    tau = tau[!is.na(pvals)]
    pvals = pvals[!is.na(pvals)]
    
    weighted.mean = stats::weighted.mean(pvals, rll)
    weighted.median = stats::median(rep(pvals, times=rll))
    weighted.sd = stats::sd(rep(pvals, times=rll))
    
    weighted.mean.tau = stats::weighted.mean(tau, rll)
    weighted.median.tau = stats::median(rep(tau, times=rll))
    weighted.sd.tau = stats::sd(rep(tau, times=rll))
    
    df_result = dplyr::tibble(
      kendall.gcContent.mean = mean(pvals),
      kendall.gcContent.median = stats::median(pvals),
      kendall.gcContent.wmean = weighted.mean,
      kendall.gcContent.wmedian = weighted.median,
      kendall.gcContent.wsd = weighted.sd,
      kendall.gcContent.sd = sd(pvals),
      kendall.gcContent.weighted.mean.tau = weighted.mean.tau,
      kendall.gcContent.weighted.median.tau = weighted.median.tau,
      kendall.gcContent.weighted.sd.tau = weighted.sd.tau,
      kendall.gcContent.mean.tau = mean(tau),
      kendall.gcContent.median.tau = median(tau),
      kendall.gcContent.sd.tau = sd(tau)
    )
    
    return(df_result)
  }
  
  dfs = future.apply::future_lapply(1:dim(segmentedCounts)[[2]], iterate_processing)

  df = do.call("rbind", dfs)

  colnames(df) = paste0("cellcycle.", colnames(df))
  description <- pData(segmentedCounts)
  
  # necessary in order not to loose row names
  cnames = rownames(description)
  description = dplyr::bind_cols(description, df)
  rownames(description) = cnames
  pData(segmentedCounts) = description
  
  return(segmentedCounts)
}


#' #'
#' #' NOT IMPLEMENTED - option to extract genomic information outside of QDNAseq
#' #'
#' per_sample_gc <- function(bamfile, resolution=10){
#'
#'   require(Rsamtools, quietly=TRUE)
#'   require(GenomicAlignments, quietly = TRUE)
#'   stopifnot(file.exists(bamfile))
#'
#'   # load precomputed bins with 1kb for GRCh37
#'   # gr = readRDS(file=paste0("~/scAbsolute/data/replicationTiming/GenomicRanges-GRCh37-", resolution, "kb.RDS"))
#'   gc_content_bias = values(gr)[["mcols.gc_baseline"]]
#'
#'   ## process set of regions at a time to speed things up
#'   calculateSubCoverage <- function(range, bam){
#'     ## read bam file from given ranges
#'     ## filter out duplicated reads, secondary reads and unmapped reads
#'
#'     ## exclude reads with mapQ==0
#'     param = Rsamtools::ScanBamParam(flag=scanBamFlag(isProperPair=TRUE,
#'                                                      isUnmappedQuery=FALSE,
#'                                                      isSecondaryAlignment=FALSE,
#'                                                      isDuplicate=FALSE),
#'                                     which=range, what="seq",
#'                                     mapqFilter=1)
#'
#'     ## read alignment
#'     alignment = GenomicAlignments::readGAlignments(bam,param=param)
#'
#'     # split alignments up in per range groups
#'     overlap = countOverlaps(alignment, range, ignore.strand=TRUE, type="within")
#'     alignment[overlap == 0] = NULL
#'     hits = findOverlaps(alignment, range, ignore.strand=TRUE, type="within")
#'     elements = S4Vectors::split(alignment, subjectHits(hits))
#'     indices = as.numeric(names(elements))
#'
#'     ## calculate gc content
#'     gc_cont = NULL
#'     nre = NULL
#'     for(j in 1:length(elements)){
#'       # print(j)
#'       al = elements[[j]]
#'       alignedSeqs <- mcols(al)$seq
#'       alf0 <- alphabetFrequency(alignedSeqs, as.prob=TRUE, collapse=TRUE)
#'       gc_co <- sum(alf0[c("G", "C")])
#'       reads = length(alignedSeqs)
#'       gc_cont = c(gc_cont, gc_co)
#'       nre = c(nre, reads)
#'     }
#'
#'     # create per range group values
#'     gc_content = rep(NA, length(range))
#'     reads = rep(NA, length(range))
#'     gc_content[indices] = gc_cont
#'     reads[indices] = nre
#'
#'     ## calculate coverage
#'     cov = GenomicAlignments::coverage(alignment)
#'     cov = cov[range]
#'     coverage = as.numeric(mean(cov))
#'
#'     return(list(gc_content=gc_content, reads=reads, coverage=coverage))
#'   }
#'
#'
#'   gc_content = NULL
#'   coverage = NULL
#'   reads = NULL
#'   nRegion = length(gr)
#'   # print(nRegion)
#'   for (i in seq(1,nRegion,10000)){
#'     # print(i)
#'     ## report progress
#'     if(i%%100000==1&i>1) cat(paste(i-1,"regions processed\n"))
#'     end = ifelse(i+9999>nRegion,nRegion,i+9999)
#'     sub_depth = calculateSubCoverage(gr[i:end],bamfile)
#'
#'     gc_content = c(gc_content, sub_depth$gc_content)
#'     reads = c(reads, sub_depth$reads)
#'     coverage = c(coverage, sub_depth$coverage)
#'   }
#'
#'   values(gr) <- S4Vectors::DataFrame(gc_content = gc_content, coverage=coverage, reads=reads, gc_content_bias=gc_content_bias)
#'
#'   return(gr)
#' }
