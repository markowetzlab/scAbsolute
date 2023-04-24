# Copyright 2022, Michael Schneider, All rights reserved.
#' copynumberSegmentation
#'
#' \code{copynumberSegmentation} use HMM for copy number inference 
#' (requires single cell input, that has been scaled and preprocessed (alpha, rpc, and name metadata))
#'
#' @param countsObject a QDNAseq object with segmented slot
#' @param change_prob numeric total probability of changing state (is divided by num_states-1)
#' @param max_iterations maximum number of iterations for optimization
#' @param max_states maximum number of copy number states
#' @param hmm_path path to store hmm object (for post-processing)
#' @param marginal_path path to store state posterior marginals
#' @param learning_rates list of learning rates to use in optimization
#' @param verbose
#'
#' @return a list with the updated object and inference results (result)
#'
#' @export
copynumberSegmentation <- function(countsObject, change_prob=1e-1,
                                   max_iterations=101, max_states=8,
                                   hmm_path=NULL,
                                   learning_rates = list(0.1, 0.01),
                                   verbose=FALSE, gc_correction=TRUE,
                                   splitPerChromosome=FALSE){
  
  
  stopifnot(dim(countsObject)[[2]] == 1)
  stopifnot(all(c("alpha", "rpc", "name") %in% colnames(Biobase::pData(countsObject))))
  reticulate::source_python(file.path(BASEDIR, "R/segmentation.py"), convert=TRUE)
  
  valid = binsToUseInternal(countsObject)
  data = Biobase::assayDataElement(countsObject, "calls")[valid,1]
  
  rpc = Biobase::pData(countsObject)[["rpc"]]
  name = Biobase::pData(countsObject)[["name"]]
  
  # estimate gc correction is part of scAbsolute/core
  if(gc_correction){
    covariate_info = estimate_gc_correction(countsObject[valid,1])[,1]  
  }else{
    covariate_info = rep(1.0, sum(valid))
  }
  
  # estimate alpha - use segmentation to estimate copy number state
  cn_estimate = round(Biobase::assayDataElement(countsObject, "segmented")[,1], digits=0)
  cn_estimate = ifelse(cn_estimate >= max_states, max_states-1, cn_estimate)
  countsObject = Biobase::assayDataElementReplace(countsObject, "copynumber", matrix(cn_estimate, ncol=1))
  alpha_estimate = estimate_overdispersion(countsObject, robust=c(0, max_states-1))
  ifelse(is.na(alpha_estimate) || is.nan(alpha_estimate), 0.01, alpha_estimate)
  
  if(splitPerChromosome){
    
    # identify chromosome boundaries (prior for breakpoints)
    # chromosomal breakpoints (breakpoints should be more likely at chromosome boundaries)
    n_bins = sum(valid)
    info = rownames(countsObject)[valid]
    seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)
    rl = Rle(chromosome)
    chromosomeIndices = cumsum(rl@lengths)
    chromosomeIndices = chromosomeIndices[1:(length(chromosomeIndices)-1)]
    
    chroms = unique(chromosome)
    df_result = list()
    state_posterior = c();
    state_marginals = matrix(data=NA, nrow=0, ncol=max_states)
    c_alpha=c(); c_rpc=c(); c_rpc_zero=c(); c_alpha_zero=c(); c_logprob = c(); c_alpha_estimate = c();
    c_oscillations = c(); c_sd_oscillations = c(); c_magnitude_oscillations = c();
    counter = 0
    
    for(c in chroms){
      counter = counter + 1
      
      if(is.null(hmm_path)){
        hmm_path = ""
      }
      out_path = gsub(pattern="//", replacement="/", paste0(base::path.expand(hmm_path),"/"))
      hmm_path2 = paste0(out_path, name, "/", name, ".", c, ".hmm")
      marginal_path = paste0(out_path, name, "/", name, ".", c, ".marginals.rds")
      
      # inference part - for a set of learning rates
      result <- interface_hmm(matrix(data[chromosome == c], ncol = 1), as.double(rpc), as.double(alpha_estimate), matrix(covariate_info[chromosome == c], nrow=1),
                              as.double(change_prob),
                              as.character(name),
                              chromosome_breakpoints=list(),
                              max_iterations=as.integer(max_iterations),
                              num_states=as.integer(max_states),
                              learning_rates=learning_rates, verbose=verbose,
                              hmm_path=as.character(hmm_path2))
      
      # save parameters, and export marginal posterior
      df_result_tmp = dplyr::tibble(alpha=result$opt$alpha, rpc=result$opt$rpc, rpc_zero=result$opt$rpc_zero, alpha_zero=result$opt$alpha_zero,
                                logprob=result$opt$log_prob, learning_rate=result$opt$learning_rate, epoch1=result$opt$epoch1,
                                marginal_path=marginal_path, hmm_path=hmm_path2, out_path=out_path, alpha_estimate=alpha_estimate)
      colnames(df_result_tmp) = paste0("hmm2.", c, ".", colnames(df_result_tmp))
      df_result[[counter]] = df_result_tmp
      state_posterior = c(state_posterior, result[["state_posterior"]])
      state_marginals = base::rbind(state_marginals, result$state_marginals)
      
      c_alpha=c(c_alpha, result$opt$alpha)
      c_alpha_estimate = c(c_alpha_estimate, alpha_estimate)
      c_rpc=c(c_rpc, result$opt$rpc)
      c_rpc_zero=c(c_rpc_zero, result$opt$rpc_zero)
      c_alpha_zero=c(c_alpha_zero, result$opt$alpha_zero)
      c_logprob = c(c_logprob, result$opt$log_prob)
      c_oscillations = c(c_oscillations, result$opt$n_oscillations)
      c_sd_oscillations = c(c_sd_oscillations, result$opt$sd_oscillations)
      c_magnitude_oscillations = c(c_magnitude_oscillations, result$opt$magnitude_oscillations)
    }
    
    df_result = dplyr::bind_cols(df_result)
    
    df_result[["hmm.alpha"]] = median(c_alpha)
    df_result[["hmm.alpha_estimate"]] = median(c_alpha_estimate)
    df_result[["hmm.rpc"]] = median(c_rpc)
    df_result[["hmm.rpc_zero"]] = median(c_rpc_zero)
    df_result[["hmm.alpha_zero"]] = median(c_alpha_zero)
    df_result[["hmm.n_nan"]] = sum(is.na(c_logprob))
    df_result[["hmm.n_oscillations"]] = sum(c_oscillations)
    df_result[["hmm.med_n_oscillations"]] = median(c_oscillations)
    df_result[["hmm.magnitude_oscillations"]] = sum(c_magnitude_oscillations)
    df_result[["hmm.sd_oscillations"]] = median(c_sd_oscillations)
    
    
    
  }else{
    

    if(is.null(hmm_path)){
      hmm_path = ""
    }
    out_path = gsub(pattern="//", replacement="/", paste0(base::path.expand(hmm_path),"/"))
    hmm_path = paste0(out_path, name, ".hmm")
    marginal_path = paste0(out_path, name, ".marginals.rds")
  
    # identify chromosome boundaries (prior for breakpoints)
    # chromosomal breakpoints (breakpoints should be more likely at chromosome boundaries)
    n_bins = sum(valid)
    info = rownames(countsObject)[valid]
    seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)
    rl = Rle(chromosome)
    chromosomeIndices = cumsum(rl@lengths)
    chromosomeIndices = chromosomeIndices[1:(length(chromosomeIndices)-1)]
    
    # inference part - for a set of learning rates
    result <- interface_hmm(matrix(data, ncol = 1), as.double(rpc), as.double(alpha_estimate), matrix(covariate_info, nrow=1),
                            as.double(change_prob),
                            as.character(name),
                            chromosome_breakpoints=as.list(chromosomeIndices),
                            max_iterations=as.integer(max_iterations),
                            num_states=as.integer(max_states),
                            learning_rates=learning_rates, verbose=verbose,
                            hmm_path=as.character(hmm_path))
    
    # save parameters, and export marginal posterior
    df_result = dplyr::tibble(alpha=result$opt$alpha, rpc=result$opt$rpc, rpc_zero=result$opt$rpc_zero, alpha_zero=result$opt$alpha_zero,
                              logprob=result$opt$log_prob, learning_rate=result$opt$learning_rate, epoch1=result$opt$epoch1,
                              marginal_path=marginal_path, hmm_path=hmm_path, out_path=out_path, alpha_estimate=alpha_estimate,
                              n_oscillations=result$opt$n_oscillations, sd_oscillations=result$opt$sd_oscillations, magnitude_oscillations=result$opt$magnitude_oscillations)
    colnames(df_result) = paste0("hmm.", colnames(df_result))
  }
    
  pd = Biobase::pData(countsObject)
  pd = dplyr::bind_cols(pd, df_result)
  Biobase::pData(countsObject) = pd
  copynumber_estimate = rep(NA, length(data))
  if(splitPerChromosome){
    copynumber_estimate[valid] = as.integer(state_posterior)
  }else{
    copynumber_estimate[valid] = as.integer(result[["state_posterior"]])  
  }
  
  countsObject = Biobase::assayDataElementReplace(countsObject, "copynumber", matrix(copynumber_estimate, ncol=1))
  
  # save marginal posterior to file, state_marginals are already logits
  if(splitPerChromosome){
    marginal_posterior = state_marginals
  }else{
    marginal_posterior = result[["state_marginals"]]  
  }
  
  if(marginal_path != ""){
    saveRDS(marginal_posterior, marginal_path)  
  }
  
  return(list("object"=countsObject))
}
                    

#' segment
#'
#' \code{segment} segment raw read counts using the PELT (Pruned Exact Linear Time) algorithm with a suitable choice of likelihood function
#'
#' @param countsObject QDNAseq object with absolute read counts in copynumber slot
#' @param penalty (see changepoint package) to use with PELT method.
#' Choice of "None", "SIC", "BIC", "MBIC", AIC", "Hannan-Quinn", "Manual"
#' If Manual is specified, the manual penalty is contained in the pen.value parameter. The predefined penalties listed DO count the
#' changepoint as a parameter, postfix a 0 e.g."SIC0" to NOT count the changepoint as a parameter.
#' Default is "BIC".
#' @param pen.value (see changepoint package) The value of the penalty when
#' using the Manual penalty option - this can be a numeric value or text giving the formula to use.
#' Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood,
#' tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
#' @param testStatistic model to use to determine optimal segmentation
#' By default we are using a Negative Binomial model, but it is also possible to specify a Normal or Poisson model here
#' (mostly for reasons of baseline models). Options are "NegBinML", "NegBinMM", "NegBinBCMM", "NegBinCML", "Poisson", "Normal", "non-parametric", "FPOP-NegBin".
#' @param gap numeric optional minimum distance between means in FPOP NegBin model (only in conjucntion with FPOP-NegBin-GAP)
#' @param alpha numeric global alpha value for NegBin
#' @param rpc numeric global reads per copy per bin value
#' 
#' @return A QDNAseq object with segmented slot
#'
#' @export
segment <- function(countsObject, penalty="MBIC", pen.value=NULL, testStatistic="NegBinMM", 
                    splitPerChromosome=TRUE, debug=FALSE, gap=NULL, alpha=NA, rpc=NA, rescale=FALSE){
  stopifnot(class(countsObject)[[1]] == "QDNAseqCopyNumbers")
  stopifnot(!(testStatistic == "NegBin" && is.na(alpha)))
  stopifnot(length(alpha) == 1)
  stopifnot(length(alpha) == length(rpc))

  stopifnot(testStatistic %in% c("NegBinML", "NegBinMM", "NegBinBCMM", 
                                 "NegBinCML", "NegBin", "NegBinAlpha", "Poisson", "Normal",
                                 "non-parametric", "np", 
                                 "FPOP-NegBin", "FPOP-NegBin-GAP",
                                 "FPOP-Normal-Posthoc",
                                 "OptSeg-SI"))

  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)

  n_cells = ncol(countsObject)
  n_bins = dim(countsObject)[[1]]
  location = rownames(countsObject)

  segmentCell <-  function(xs) {

      n = length(xs)
      idx_na = which(is.na(xs))
      idx_ok = which(!is.na(xs))
      tmp = xs[idx_ok]
      n_good = length(tmp)

      # NOTE for compilation of c dependencies
      # R CMD SHLIB cost_general_functions.c PELT_one_func_minseglen.c -o PELT.so
      source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
      if (testStatistic != "non-parametric" && testStatistic != "np"){
        dyn.load(file.path(BASEDIR, "data/changepoint/PELT.so"))
      }

      if(testStatistic == "NegBinML"){
        costfunc = "meanvar.negbin_ml"
      }else if (testStatistic == "NegBinCML"){
        costfunc = "meanvar.negbin_cml"
      }else if (testStatistic == "NegBinMM"){
        costfunc = "meanvar.negbin_mm"
      }else if (testStatistic == "NegBinBCMM"){
        costfunc = "meanvar.negbin_bcmm"
      }else if (testStatistic == "NegBinAlpha"){
        costfunc = "meanvar.negbin_alpha"
      }else if (testStatistic == "NegBin"){
        costfunc = "meanvar.negbin"
      }else if (testStatistic == "Normal"){
        costfunc = "meanvar.norm"
      } else if (testStatistic == "Poisson"){
        costfunc = "meanvar.poisson"
      } else if (testStatistic == "non-parametric" || testStatistic == "np"){
        costfunc = "non_parametric"
      } else if (testStatistic == "FPOP-NegBin"){
        costfunc = "fpop.negbin"
      } else if (testStatistic == "FPOP-NegBin-GAP"){
        costfunc = "fpop.negbin.gap"
      } else if (testStatistic == "FPOP-Normal-Posthoc"){
        costfunc = "fpop.normal.posthoc"
      }
      
      if(splitPerChromosome){
        ## segment per chromosome
        seqName = unlist(lapply(strsplit(location, ":"), `[`, 1))
        chrOrder = base::rle(seqName)[["values"]]
        # NOTE important to subset for items used in cpt analysis
        chromosome <- factor(seqName, levels=chrOrder)[idx_ok]
        stopifnot(length(chromosome) == length(tmp))
        rl = Rle(chromosome)
        if(length(rl@values) == 1){
          offsets = c(0)
        }else{
          offsets = c(0, cumsum(rl@lengths)[1:(length(rl@lengths)-1)])
        }
        perChromosomeCounts = split(tmp, chromosome)
        cpt_split = lapply(perChromosomeCounts, function(x){
      
            if(length(x) == 0){
              # Fix for Y chromosome
              return(NA)
            }
            
            if(costfunc == "non_parametric"){
              require(changepoint.np)
              s = changepoint.np::cpt.np(data = as.integer(x), pen.value=pen.value, penalty = penalty, class=FALSE)
            }else{
              dyn.load(file.path(BASEDIR, "data/changepoint/PELT.so")) # necessary for future.apply to work
              if(costfunc == "meanvar.negbin_alpha"){
              alpha <- tryCatch(
                      {
                          alpha = 1.0 / suppressWarnings(MASS::fitdistr(as.integer(x), densfun = "negative binomial")[["estimate"]][["size"]])
                      },
                      error=function(cond) {
                          return(0.01)
                      }
              )
              }
              s = wrap_segmentation(as.integer(x), pen.value=pen.value,
                        penalty=penalty, costfunc=costfunc, debug=debug, gap=gap, alpha=alpha, rpc=rpc)
              debug = FALSE
            }

            return(s)})
        # Remove chromsome with only NA values (Y chromosome)
        cpt_split[sapply(cpt_split,function(x) (length(x) != 0) && all(is.na(x)))] <- NULL
        cpt = unlist(lapply(1:length(cpt_split), function(x) return(cpt_split[[x]] + offsets[x])))
      }else{
        if(costfunc == "non_parametric"){
          require(changepoint.np)
          cpt = changepoint.np::cpt.np(data = as.integer(tmp), pen.value=pen.value, penalty = penalty, class=FALSE)
        }else{
              if(costfunc == "meanvar.negbin_alpha"){
                alpha <- tryCatch(
                        {
                            alpha = 1.0 / suppressWarnings(MASS::fitdistr(as.integer(tmp), densfun = "negative binomial")[["estimate"]][["size"]])
                        },
                        error=function(cond) {
                            return(0.01)
                        }
                )
              }
              cpt = wrap_segmentation(as.integer(tmp), pen.value=pen.value,
                                      penalty=penalty, costfunc = costfunc, debug=debug, gap=gap, alpha=alpha, rpc=rpc)
          }
        }
      
      K = length(cpt)

      sol = cpt
      if(sol[length(sol)] != n_good){
        sol = c(sol, n_good)
      }
      
      segs = rep(c(1:length(sol)), times=c(sol[1], diff(sol)))
      reg = S4Vectors::Rle(segs)
      ir = IRanges::ranges(reg)
      
      segmentedValues = smoothOperator(tmp, ir, function(xts){return(stats::median(xts, na.rm=TRUE))}, minLength=0, trimLength=0)
      lengthSegmentedValues = smoothOperator(tmp, ir, function(xts){return(length(xts))}, minLength=0, trimLength=0)
      
      segmentedGenome = rep(segmentedValues, times=lengthSegmentedValues)
      stopifnot(length(segmentedGenome) == length(idx_ok))
      
      if(rescale){
        segmentedGenome = round(segmentedGenome / rpc, digits = 0)
      }
      
      xs[idx_ok] = segmentedGenome
      return(xs)

  }
  z <- apply(countsObject@assayData$copynumber, 2, FUN=segmentCell)

  if ("segmented" %in% names(countsObject@assayData)) {
    countsObject = assayDataElementReplace(countsObject, 'segmented', z)
  } else {
    assayDataElement(countsObject, 'segmented') = z
  }

  return(countsObject)
}


#' segment_CROPS
#'
#' \code{segment_CROPS} segment raw read counts using the CROPS algorithm. Selection is based on mean-variance fit.
#'
#' @param countsObject QDNAseq object with absolute read counts in copynumber slot
#' @param pen.value (see changepoint package) The value of the penalty when
#' using the Manual penalty option - this can be a numeric value or text giving the formula to use.
#' Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood,
#' tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
#' @param trimLength bins at both ends of segments that are not included in the median read counts per segment computation
#' Default is set to 1, in order to remove potential breakpoint bins which might not be representative of the median segment level
#' @param testStatistic model to use to determine optimal segmentation
#' By default we are using a Negative Binomial model, but it is also possible to specify a Normal or Poisson model here
#' (mostly for reasons of baseline models). Options are "NegBinML", "NegBinMM", "NegBinBCMM", "NegBinCML", "Poisson", "Normal", "non-parametric".
#'
#' @return A QDNAseq object with segmented slot
#'
#' @export
segment_CROPS <- function(countsObject, pen.value=NULL, testStatistic="NegBinMM", debug=TRUE){
  stopifnot(class(countsObject)[[1]] == "QDNAseqCopyNumbers")
  
  stopifnot(testStatistic %in% c("NegBinML", "NegBinMM", "NegBinBCMM", "NegBinCML", "Poisson", "Normal", "np", "non-parametric"))
  
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  
  n_cells = ncol(countsObject)
  n_bins = dim(countsObject)[[1]]
  location = rownames(countsObject)
  
  segmentCell <-  function(xs) {
    
    n = length(xs)
    idx_na = which(is.na(xs))
    idx_ok = which(!is.na(xs))
    tmp = xs[idx_ok]
    n_good = length(tmp)
    
    # NOTE update for package deployment
    # R CMD SHLIB cost_general_functions.c PELT_one_func_minseglen.c -o PELT.so
    if(testStatistic != "non-parametric" && testStatistic != "np"){
      source(file.path(BASEDIR, "data/changepoint/wrap_PELT.R"))
      dyn.load(file.path(BASEDIR, "data/changepoint/PELT.so"))
    }
    
    if(testStatistic == "NegBinML"){
      costfunc = "meanvar.negbin_ml"
    }else if (testStatistic == "NegBinCML"){
      costfunc = "meanvar.negbin_cml"
    }else if (testStatistic == "NegBinMM"){
      costfunc = "meanvar.negbin_mm"
    }else if (testStatistic == "NegBinBCMM"){
      costfunc = "meanvar.negbin_bcmm"
    }else if (testStatistic == "Normal"){
      costfunc = "meanvar.norm"
    } else if (testStatistic == "Poisson"){
      costfunc = "meanvar.poisson"
    } else if (testStatistic == "non-parametric" || testStatistic == "np"){
      costfunc = "non_parametric"
    }

    
    cpt_list = wrap_CROPS_segmentation(as.integer(tmp), costfunc = costfunc, pen.value=pen.value)
    
    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
    stopifnot(length(cpt_list$changepoints) >= 1)
    
    rsq = function(xt){
      splits = splitAt(as.integer(tmp), xt)
      ms = unlist(lapply(splits, mean, na.rm=TRUE))
      vars = unlist(lapply(splits, var, na.rm=TRUE))
      weights = unlist(lapply(splits, length))
      valid = !(is.na(ms)  | is.na(vars) | is.na(weights))
      stopifnot(length(ms) == length(vars))
      dat = data.frame(l_mean=ms[valid], l_var=vars[valid], l_weight=weights[valid])
      
      beta_rsquared <- tryCatch(
        {
          # message("This is the 'try' part")
          fit_beta = lm(l_var - I(l_mean) ~ I(l_mean^2) + 0, weights=l_weight, data=dat)
          rsquared = summary(fit_beta)[["adj.r.squared"]]
          rsquared
        },
        error=function(cond) {
          message("Error when fitting beta regression in scSegment - segment_CROPS")
          print(cond)
          return(NA)
        },
        warning=function(cond) {
          message("Warning when fitting beta regression in scSegment - segment_CROPS")
          print(cond)
          return(NA)
          # return(NULL)
        },
        finally={
          # NOTE:
          # Here goes everything that should be executed at the end, regardless of success or error.
        }
      )
      
      # ggplot(data = dat) + geom_point(aes(x=l_mean, y=l_var, size=l_weight)) +
      #   stat_smooth(method=rlm, formulat = y ~ x, aes(x=l_mean,y=l_var, weight=l_weight))

      return(beta_rsquared)
    }
    
    rsq_list = unlist(lapply(cpt_list$changepoints, rsq))
    rsq_list[is.na(rsq_list)] = 0
    if(debug){
      print(paste0("rsquared list: ", rsq_list))  
    }
    cpt = cpt_list$changepoints[[which.max(rsq_list)]]
    
    K = length(cpt)
    
    sol = cpt
    if(sol[length(sol)] != n_good){
      sol = c(sol, n_good)
    }
    
    segs = rep(c(1:length(sol)), times=c(sol[1], diff(sol)))
    reg = S4Vectors::Rle(segs)
    ir = IRanges::ranges(reg)
    
    segmentedValues = smoothOperator(tmp, ir, function(xts){return(stats::median(xts, na.rm=TRUE))}, minLength=0, trimLength=0)
    lengthSegmentedValues = smoothOperator(tmp, ir, function(xts){return(length(xts))}, minLength=0, trimLength=0)
    
    segmentedGenome = rep(segmentedValues, times=lengthSegmentedValues)
    stopifnot(length(segmentedGenome) == length(idx_ok))
    
    xs[idx_ok] = segmentedGenome
    return(xs)
    
  }
  z <- apply(countsObject@assayData$copynumber, 2, FUN=segmentCell)
  
  if ("segmented" %in% names(countsObject@assayData)) {
    countsObject = assayDataElementReplace(countsObject, 'segmented', z)
  } else {
    assayDataElement(countsObject, 'segmented') = z
  }
  
  return(countsObject)
}


#' @DEPRECATED
#' #' segmentSingleBreakpoint
#' #'
#' #' \code{segmentSingleBreakpoint} find single optimal breakpoint in sequence based on Negative Binomial count model
#' #'
#' #' @param data count data vector
#' 
#' #' @return out vector with changepoint position and likelihood for null model and alternative model
#' #'
#' #' @export
#' segmentSingleBreakpoint <- function(data){
#'   minseglen=5
#'   require(MASS, warn.conflicts = FALSE)
#'   stopifnot(length(data) >= 10)
#'   stopifnot(is.null(dim(data)))
#'   n = length(data)
#' 
#'   ## robust computation of negative binomial - fall back on poisson
#'   compute_ll <- function(x){
#'     suppressWarnings(
#'       tmp <- tryCatch(
#'         {
#'           null = MASS::fitdistr(x, densfun="negative binomial")[["loglik"]]
#'           tmp2 = rep(0, times=n)
#'           for(i in (minseglen):(n-minseglen+1)){
#'             tmp2[i] = tmp2[i] + MASS::fitdistr(x[1:i], densfun="negative binomial")[["loglik"]]
#'             tmp2[i] = tmp2[i] + MASS::fitdistr(x[(i+1):n], densfun="negative binomial")[["loglik"]]
#'           }
#'           tmp2[1:(minseglen-1)] = NA
#'           tmp2[(n-minseglen+2):n] = NA
#'           return(list(tmp=tmp2, null=null))
#'         },
#'         error=function(cond) {
#'           null = MASS::fitdistr(x, densfun="Poisson")[["loglik"]]
#'           tmp2 = rep(0, times=n)
#'           for(i in (minseglen):(n-minseglen+1)){
#'             tmp2[i] = tmp2[i] + MASS::fitdistr(x[1:i], densfun="Poisson")[["loglik"]]
#'             tmp2[i] = tmp2[i] + MASS::fitdistr(x[(i+1):n], densfun="Poisson")[["loglik"]]
#'           }
#'           tmp2[1:(minseglen-1)] = NA
#'           tmp2[(n-minseglen+2):n] = NA
#'           return(list(tmp=tmp2, null=null))
#'         }))
#'   }
#'   
#'   model = compute_ll(data)
#'   tmp = model[["tmp"]]
#'   tau=which(tmp==max(tmp,na.rm=T))[1]
#'   taulike=tmp[tau]
#'   out=list(cpt=tau,null=model[["null"]],alt=taulike)
#'   # names(out)=c('cpt','null','alt')
#'   return(out)
#' }
