# Copyright 2022, Michael Schneider, All rights reserved.
## Mean-Variance relationship

#' computeModel
#'
#' \code{computeModel} compute alpha and beta for a given segmentation
#'
#' @param counts read counts (as taken from calls or copynumber slot)
#' @param gr GRanges object, representing a segmentation of the genome
#' @param limitPloidy Numeric cutoff values for ploidy values to include in fitting
#' @param scale Numeric scale used to transform read counts
#' @param minLength Numeric minimum segment size to include in fitting
#' @param trimLength Numeric number of bins to remove from start and end of segments for estimates
#' @param debug Boolean return additional output such as plot, excluded elements and fit
#' @param scaffold Boolean return dummy output for elements that cannot be fitted reasonably
#' @param skip_alpha Boolean skip computation of alpha altogether
#' 
#' NOTE:
#' A) counts is raw count data
#' B) scale is transform, s.t. raw * scale = (scaled copy number counts)
#' C) segmentation can be a) rounded scaled data or b) not rounded scaled data (scaled vs not scaled only matters for rounding)
#'    Rounding is only relevant for single cell DNAseq data, where expected levels can only occur at integer states.
#' 
#' @return data frame containing the model parameters
computeModel <- function(counts, gr, limitPloidy, scale=1.0, minLength=10, trimLength=1, roundInteger=FALSE, debug=FALSE, scaffold=FALSE){
  
  if(!scaffold){
  
      ## Initialize variables ####
      suppressPackageStartupMessages({
      require(MASS, quietly = TRUE, warn.conflicts = FALSE)
      require(gamlss, quietly = TRUE, warn.conflicts = FALSE)
      require(mgcv, quietly = TRUE, warn.conflicts = FALSE)
      require(robustbase, quietly = TRUE, warn.conflicts = FALSE)
      })
      lmrob_control = robustbase::lmrob.control(setting="KS2014",
                                                max.it = 10000, maxit.scale = 10000,
                                                scale.tol=1e-5, solve.tol = 1e-5,
                                                subsampling = "nonsingular")
      warning_handler <- function(w){ 
        
        if(length(w) > 1 && any(grepl("Error", w[[1]]))){
          stop()
        }
        
        if( any( grepl( "Detected possible local breakdown of SMDM-estimate", w) | 
                                               grepl("find_scale\\(\\) did not converge in", w) | 
                                               grepl("Convergence failure: singular convergence", w)) ){
        invokeRestart( "muffleWarning" )}
      }
      
      # # remove extreme outliers, 99% quantile -> ISSUE WITH OVERDISPERSION VALUE
      # xy = lapply(ccounts, function(x) return(decode(x)))
      # xy_cutoff = quantile(unlist(xy), 0.99, na.rm=TRUE)
      # counts = RleList(lapply(ccounts, function(x){y=decode(x); y[!is.na(y) & y >= xy_cutoff] = NA; return(Rle(y))}))

      # compute per segments estimates
      l_mean = smoothOperator(counts * scale, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
      ## Single cell property -> means at integer state
      if(roundInteger){
        l_mean = base::round(l_mean, digits = 0)
      }
      l_var = smoothOperator(counts * scale, gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
      l_weight = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = trimLength, minLength=minLength)
      l_weight_nseg = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = 0, minLength=1)
      
      l_var_unscaled = smoothOperator(counts, gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
      l_mean_unscaled = smoothOperator(counts, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
      l_mean_squared = smoothOperator(counts^2, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
      
      # required at gam model stage
      filter_lw_counts = smoothOperator(counts, gr, function(x){return(x[!is.na(x)])}, trimLength = trimLength, minLength=minLength)
      filter_l_weight = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = trimLength, minLength=minLength)
      filter_l_mean_unscaled = smoothOperator(counts, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
        
      lw_counts = smoothOperator(counts, gr, function(x){return(x[!is.na(x)])}, trimLength = 0, minLength=1)
      lw_mean = smoothOperator(counts, gr, function(x){return(mean(x, na.rm=TRUE))}, trimLength = 0, minLength=1)
      lw_mean = lw_mean[!is.na(lw_mean)]
      lw_mean = rep(lw_mean, times=l_weight_nseg[l_weight_nseg != 0]) * scale
      stopifnot(length(lw_counts) == length(lw_mean))
      
      # NOTE compute values for rpc and ploidy, more versions are computed below
      zero_elements = lw_mean < 0.5
      mu_zero = mean(lw_counts[zero_elements])
      mu_zero_ratio = sum(lw_counts[zero_elements] > 0) / sum(zero_elements)
      n_zero = sum(zero_elements)
      rpc.meanvar = mean(lw_counts[!zero_elements] / round(lw_mean[!zero_elements], digits=0), na.rm=TRUE)
      ploidy.meanvar = mean(round(lw_mean, digits=0), na.rm=TRUE)

      infer_alpha <- function(x){
        if(length(x) == 0 | length(x[!is.na(x)]) == 0){
          return(NA)
        }  
        
        fit <- tryCatch(
          {
            # message("This is the 'try' part")
            x = x[!is.na(x)]
            MASS::fitdistr(x, densfun = "negative binomial")
          },
          error=function(cond) {
            # message("Error in mean-variance - infer_alpha")
            # print(cond)
            # Choose a return value in case of error
            return(NA)
          },
          warning=function(cond) {
            # message("Here's the original warning message:")
            # message(cond)
            # Choose a return value in case of warning
            return(NA)
          },
          finally={
            # NOTE:
            # Here goes everything that should be executed at the end,
            # regardless of success or error.
          }
        )
        
        if(all(is.na(fit))){
          return(NA)
        }else{
          # we return inverse of size parameter for compatibility
          return(1. / as.numeric(fit$estimate["size"]))
        }
        # end function
      }

      infer_omega <- function(x){
        if(length(x) == 0 | length(x[!is.na(x)]) == 0){
          return(NA)
        }  
        
          fit <- tryCatch(
            {
              # message("This is the 'try' part")
              x = x[!is.na(x)]
              fit <- stats::glm(x ~ 1, quasipoisson)
              omega = summary(fit)[["dispersion"]]
              return(omega)
            },
            error=function(cond) {
              # message("Error in mean-variance - infer_omega")
              # print(cond)
              # Choose a return value in case of error
              return(NA)
            },
            warning=function(cond) {
              # message("Here's the original warning message:")
              # message(cond)
              # Choose a return value in case of warning
              return(NA)
            },
            finally={
              # NOTE:
              # Here goes everything that should be executed at the end,
              # regardless of success or error.
            }
          )
          
          if(all(is.na(fit))){
            return(NA)
          }else{
            # we return inverse of size parameter for compatibility
            return(fit)
          }
          # end function
        }

      l_alpha = smoothOperator(counts, gr, infer_alpha, trimLength = trimLength, minLength=minLength)
      l_omega = smoothOperator(counts, gr, infer_omega, trimLength = trimLength, minLength=minLength)
      lw_alpha = rep(l_alpha, times=l_weight)

      #############################################
      # estimates that large are very unreliable and would unduly bias the overall beta estimate
      # -> we can use limitPloidy here as l_mean should have been scaled!
      stopifnot(length(l_alpha) == length(l_mean))
      n_seg_effective = length(l_mean)
      keep_ix_1 = (l_mean <= limitPloidy & l_mean >= 0.5 & !is.na(l_var) & l_var > 0 & l_weight >= (minLength - 2*trimLength))
      keep_ix_2 = !is.na(l_alpha) & !is.na(l_omega) & keep_ix_1
      
      l_alpha = l_alpha[keep_ix_2]
      l_omega = l_omega[keep_ix_2]
      l_mean_alpha = l_mean[keep_ix_2]
      l_mean_unscaled_alpha = l_mean_unscaled[keep_ix_2]
      l_var_alpha = l_var[keep_ix_2]
      l_var_unscaled_alpha = l_var_unscaled[keep_ix_2]
      l_weight_alpha = l_weight[keep_ix_2]
      l_mean_squared_alpha = l_mean_squared[keep_ix_2]
      if(length(l_alpha) == 0){
        return(NULL)        
      }

      l_mean = l_mean[keep_ix_1]
      l_var = l_var[keep_ix_1]
      l_var_unscaled = l_var_unscaled[keep_ix_1]
      l_mean_unscaled = l_mean_unscaled[keep_ix_1]
      l_weight = l_weight[keep_ix_1]
      l_mean_squared = l_mean_squared[keep_ix_1]
      
      if(length(l_mean) == 0){
        return(NULL)
      }
      
      rpc_meanvar = sum(lw_counts, na.rm=TRUE) / sum(round(lw_mean, digits = 0))
      rpc_meanvar2 = 1/scale
      ploidy_meanvar = mean(round(lw_mean, digits=0), na.rm=TRUE)
      ploidy_meanvar2 = mean(decode(scale*unlist(counts)), na.rm=TRUE)
      ploidy_meanvar3 = sum(gr@ranges@width * gr$segmentation, na.rm=TRUE) / sum(gr@ranges@width, na.rm=TRUE)
      n_segs = sum(l_weight_nseg != 0)
      n_segs_na = sum(l_weight_nseg == 0)

      #############################################
      stopifnot(all(l_weight >= minLength - 2*trimLength))
      
      ## covariates
      beta.n_seg = length(l_var)
      beta.percent_used = sum(l_weight) / sum(gr@ranges@width)
      
      ## Beta model ####
 
      ## beta raw model - predict variance given mean ----
      fit_unscaled_beta <- tryCatch(
        {
          # message("This is the 'try' part")
          fit_unscaled_beta = withCallingHandlers( robustbase::lmrob(var_unscaled ~ I(mean_unscaled) + I(mean_unscaled^2) + 0, weights=weight,
                                                data = data.frame(var_unscaled=l_var_unscaled,
                                                                  mean_unscaled=l_mean_unscaled,
                                                                  weight=l_weight),
                                                control=lmrob_control), warning = warning_handler )
          fit_unscaled_beta
        },
        error=function(cond) {
          message("Error in fit_unscaled_beta - mean-variance")
          message(cond)
          fit_unscaled_beta = lm(var_unscaled ~ I(mean_unscaled) + I(mean_unscaled^2) + 0, weights=weight,
                                 data = data.frame(var_unscaled=l_var_unscaled, mean_unscaled=l_mean_unscaled, weight=l_weight))
          
          fit_unscaled_beta
        },
        warning=function(cond) {
          message("Warning in fit_unscaled_beta - mean-variance")
          message(cond)
          fit_unscaled_beta = lm(var_unscaled ~ I(mean_unscaled) + I(mean_unscaled^2) + 0, weights=weight,
                                 data = data.frame(var_unscaled=l_var_unscaled, mean_unscaled=l_mean_unscaled, weight=l_weight))
          
          fit_unscaled_beta
        },
        finally={
          # NOTE:
          # Here goes everything that should be executed at the end, regardless of success or error.
        }
      )
      
      # information about variance given (unscaled) mean
      var.beta_rsquared_unscaled = summary(fit_unscaled_beta)[["adj.r.squared"]]
      var.beta_cn4 = predict(fit_unscaled_beta, newdata=data.frame(mean_unscaled = rpc_meanvar * 4))
      var.beta_cn2 = predict(fit_unscaled_beta, newdata=data.frame(mean_unscaled = rpc_meanvar * 2))
      var.beta_cn1 = predict(fit_unscaled_beta, newdata=data.frame(mean_unscaled = rpc_meanvar * 1))
      var.meanvar_ratio_4 = var.beta_cn4 / (rpc_meanvar * 4)
      var.meanvar_ratio_2 = var.beta_cn2 / (rpc_meanvar * 2)
      var.meanvar_ratio_1 = var.beta_cn1 / (rpc_meanvar * 1)
      var.coefficient_beta_unscaled_1 = fit_unscaled_beta$coefficients["I(mean_unscaled)"]
      var.coefficient_beta_unscaled_2 = fit_unscaled_beta$coefficients["I(mean_unscaled^2)"]
      if(class(fit_unscaled_beta) == "lmrob"){
        beta.converged_unscaled_beta = TRUE
        beta.rweights_unscaled_beta = fit_unscaled_beta[["rweights"]]
      }else{
        beta.converged_unscaled_beta = FALSE
        beta.rweights_unscaled_beta = rep(1.0, length(l_mean_unscaled))
      }
      beta.quality_unscaled_beta = sum(beta.rweights_unscaled_beta * l_weight) / sum(rep(1.0, length(l_mean_unscaled)) * l_weight)
      
      
      
      ## beta model for overdispersion in NB ----
      fit_negbin_beta <- tryCatch(
        {
          fit_negbin_beta = withCallingHandlers(robustbase::lmrob(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                                                                  data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight)), warning=warning_handler)
          fit_negbin_beta
        },
        error=function(cond) {
          message("Error in fit_negbin_beta - mean-variance")
          message(cond)
          fit_negbin_beta = lm(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                               data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight))
          
          fit_negbin_beta
        },
        warning=function(cond) {
          message("Warning in fit_negbin_beta - mean-variance")
          message(cond)
          fit_negbin_beta = withCallingHandlers(robustbase::lmrob(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                                                                  data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight)), warning=warning_handler)
          fit_negbin_beta
        },
        finally={
          # NOTE:
          # Here goes everything that should be executed at the end, regardless of success or error.
        }
      )
      beta.rsquared_negbin = summary(fit_negbin_beta)[["adj.r.squared"]]
      beta.coefficient_negbin = fit_negbin_beta$coefficients[["I(mean^2)"]]
      # NOTE this is alpha in the strict NB sense
      if(class(fit_negbin_beta) == "lmrob"){
        beta.converged_negbin = TRUE
        beta.rweights_negbin = fit_negbin_beta[["rweights"]]
      }else{
        beta.converged_negbin = FALSE
        beta.rweights_negbin = rep(1.0, length(l_mean_unscaled))
      }
      beta.quality_negbin = sum(beta.rweights_negbin * l_weight) / sum(rep(1.0, length(l_mean_unscaled)) * l_weight)
 
      fit_combined_beta <- tryCatch(
          {
            fit_combined_beta = withCallingHandlers(robustbase::nlrob(var ~ coefficient_omega * I(mean) + coefficient_alpha * I(mean^2) + 0, weights=weight,
                                                data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight),
                                                start = list(coefficient_omega = 1.0, coefficient_alpha = 1e-3),
                                                lower=list(coefficient_omega=0.9, coefficient_alpha=1e-6),
                                                upper=list(coefficient_omega=Inf, coefficient_alpha=Inf),
                                                maxit = 500, tol=1e-06,
                                                algorithm = "port", control=stats::nls.control(maxiter = 1000, tol = 1e-07, minFactor=1/4096, warnOnly = TRUE)), warning=warning_handler)

            fit_combined_beta
          },
          error=function(cond) {
            message("Error in fit_combined_beta - mean-variance")
            message(cond)
            fit_combined_beta = lm(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                                 data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight))
            
            fit_combined_beta
          },
          warning=function(cond) {
            message("Warning in fit_combined_beta - mean-variance")
            message(cond)
            fit_combined_beta = withCallingHandlers(robustbase::nlrob(var ~ coefficient_omega * I(mean) + coefficient_alpha * I(mean^2) + 0, weights=weight,
                                               data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight),
                                               start = list(coefficient_omega = 1.0, coefficient_alpha = 1e-6),
                                               lower=list(coefficient_omega=0.0, coefficient_alpha=0.0),
                                               upper=list(coefficient_omega=Inf, coefficient_alpha=Inf),
                                               maxit = 1000, tol=1e-06,
                                               algorithm = "port", control=stats::nls.control(maxiter = 1000, tol = 1e-07, minFactor=1/4096, warnOnly = TRUE)), warning=warning_handler)
            
            fit_combined_beta
          },
          finally={
            # NOTE:
            # Here goes everything that should be executed at the end, regardless of success or error.
          }
        )

      if(class(fit_combined_beta)[1] == "nlrob"){
        if(fit_combined_beta$status == "converged"){
          beta.converged_combined = TRUE
        }else{
          beta.converged_combined = FALSE
        }
        beta.coefficient_combined_1 = fit_combined_beta$coefficients[["coefficient_omega"]]
        beta.coefficient_combined_2 = fit_combined_beta$coefficients[["coefficient_alpha"]]
        beta.rweights_combined = fit_combined_beta[["rweights"]]
      }else{
        beta.converged_combined = FALSE
        beta.coefficient_combined_1 = 1
        beta.coefficient_combined_2 = fit_combined_beta$coefficients[["I(mean^2)"]]
        beta.rweights_combined = rep(1.0, length(l_mean_unscaled))
      }
      beta.quality_combined = sum(beta.rweights_combined * l_weight) / sum(rep(1.0, length(l_mean_unscaled)) * l_weight)
      # NOTE this is alpha and omega in the combined sense
      
      ## Classical statistical models
      keep_index_3 = !(is.na(lw_counts) | is.na(lw_mean) | lw_mean < 0.5)
      lw_counts = lw_counts[keep_index_3]
      lw_mean = lw_mean[keep_index_3]
      
      fit_0 = gamlss::gamlss(lw_counts ~ lw_mean, sigma.formula=~1, family=NBI,
                             control=gamlss.control(c.crit = 0.0001, n.cyc = 100,trace = FALSE))
      fit_1 = gamlss::gamlss(lw_counts ~ lw_mean, sigma.formula=~lw_mean, family=NBI,
                             control=gamlss.control(c.crit = 0.0001, n.cyc = 100,trace = FALSE))
      lr.test <- tryCatch(
        {
          lr.test = LR.test(fit_0, fit_1, print=FALSE)
          lr.test
        },
        error=function(cond) {
          message("Error in LR-test - mean-variance")
          message(cond)
          lr.test = list()
          lr.test[["chi"]]=NA
          lr.test[["p.val"]]=NA
          lr.test
        },
        warning=function(cond) {
          message("Warning in fit_combined_beta - mean-variance")
          message(cond)
          lr.test = LR.test(fit_0, fit_1, print=FALSE)
          lr.test
        },
        finally={
          # NOTE:
          # Here goes everything that should be executed at the end, regardless of success or error.
        }
      )
      gamlss.chi = lr.test[["chi"]]
      gamlss.pvalue = lr.test[["p.val"]]
      gamlss.mu_coefficient_0 = fit_1$mu.coefficients[["(Intercept)"]]
      gamlss.mu_coefficient_1 = fit_1$mu.coefficients[["lw_mean"]]
      gamlss.sigma_coefficient_0 = fit_1$sigma.coefficients[["(Intercept)"]]
      gamlss.sigma_coefficient_1 = fit_1$sigma.coefficients[["lw_mean"]]
      
      ## Compute classical ML model for negative binomial and quasi-poisson
      formula.gam = formula(counts ~ mean)
      model.nb <- mgcv::gam(formula.gam, data=data.frame(counts=lw_counts, mean=lw_mean), 
                         family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
      gam.alpha.devexpl = summary(model.nb)[["dev.expl"]]
      gam.alpha = 1/model.nb$family$getTheta(TRUE)
      gam.alpha.coefficient_0 = model.nb$coefficients[["(Intercept)"]]
      gam.alpha.coefficient_1 = model.nb$coefficients[["mean"]]
      
      model.qp <- mgcv::gam(formula.gam, data=data.frame(counts=lw_counts, mean=lw_mean), 
                            family=quasipoisson(link = "log"), method = "REML")
      gam.omega.devexpl = summary(model.qp)[["dev.expl"]]
      gam.omega = summary(model.qp)[["dispersion"]]
      gam.omega.coefficient_0 = model.qp$coefficients[["(Intercept)"]]
      gam.omega.coefficient_1 = model.qp$coefficients[["mean"]]
      
      ## Compute (more robust?) classical ML model for negative binomial, 
      # 1. filtering based on beta estimation unscaled
      keep_ix_3 = keep_ix_1
      condition_1 = beta.rweights_unscaled_beta >= 0.95
      stopifnot(length(condition_1) == length(l_mean))
      stopifnot(length(condition_1) == length(l_weight))

      keep_ix_3[which(keep_ix_1)[!condition_1]] = FALSE
      filter_valid = rep(keep_ix_3, times=filter_l_weight)

      filter_lw_mean_1 = rep(l_mean[condition_1], l_weight[condition_1])
      stopifnot(length(filter_lw_mean_1) == length(filter_lw_counts[filter_valid]))
      
      if(sum(condition_1) >= 3){
        formula.gam = formula(counts ~ mean)
        model.nb.filter <- mgcv::gam(formula.gam, 
                                     data=data.frame(counts=filter_lw_counts[filter_valid], 
                                                     mean=filter_lw_mean_1), 
                                     family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
        gam.filter_beta_alpha.devexpl = summary(model.nb.filter)[["dev.expl"]]
        gam.filter_beta_alpha = 1/model.nb.filter$family$getTheta(TRUE)
        gam.filter_beta_alpha.coefficient_0 = model.nb.filter$coefficients[["(Intercept)"]]
        gam.filter_beta_alpha.coefficient_1 = model.nb.filter$coefficients[["mean"]]
      }else{
        gam.filter_beta_alpha.devexpl = NA
        gam.filter_beta_alpha = NA
        gam.filter_beta_alpha.coefficient_0 = NA
        gam.filter_beta_alpha.coefficient_1 = NA
      }
 
      # 2. filtering based on neg bin via mean variance
      keep_ix_3 = keep_ix_1
      condition_2 = beta.rweights_negbin >= 0.95
      stopifnot(length(condition_2) == length(l_mean))
      stopifnot(length(condition_2) == length(l_weight))
      
      keep_ix_3[which(keep_ix_1)[!condition_2]] = FALSE
      filter_valid = rep(keep_ix_3, times=filter_l_weight)
      
      filter_lw_mean_1 = rep(l_mean[condition_2], l_weight[condition_2])
      stopifnot(length(filter_lw_mean_1) == length(filter_lw_counts[filter_valid]))
      
      if(sum(condition_2) >= 3){
        formula.gam = formula(counts ~ mean)
        model.nb.filter <- mgcv::gam(formula.gam, 
                                     data=data.frame(counts=filter_lw_counts[filter_valid], 
                                                     mean=filter_lw_mean_1), 
                                     family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
        gam.filter_negbin_alpha.devexpl = summary(model.nb.filter)[["dev.expl"]]
        gam.filter_negbin_alpha = 1/model.nb.filter$family$getTheta(TRUE)
        gam.filter_negbin_alpha.coefficient_0 = model.nb.filter$coefficients[["(Intercept)"]]
        gam.filter_negbin_alpha.coefficient_1 = model.nb.filter$coefficients[["mean"]]
      }else{
        gam.filter_negbin_alpha.devexpl = NA
        gam.filter_negbin_alpha = NA
        gam.filter_negbin_alpha.coefficient_0 = NA
        gam.filter_negbin_alpha.coefficient_1 = NA
      }
      
      
      
      ## alpha model ####
      alpha.n_seg = NA; alpha.percent_used=NA;
      alpha.rsquared_meanvar = NA; alpha.coefficient_meanvar = NA;
      alpha.rsquared_invmean = NA; alpha.coefficient_invmean = NA;
      alpha.coefficient_combined_1 = NA; alpha.coefficient_combined_2 = NA; alpha.rsquared_combined = NA;
      alpha.converged_meanvar_alpha = FALSE; alpha.converged_invmean_alpha = FALSE; alpha.converged_combined_alpha = FALSE;
      
      if(!is.null(l_alpha) && length(l_alpha) >= 2){
 
        # NOTE: l_mean_unscaled_alpha is already reduced to number of items in keep_ix_2
        alpha.n_seg = length(l_alpha)
        alpha.percent_used = sum(l_weight_alpha) / sum(gr@ranges@width)
        
        # fit_meanvar_alpha = withCallingHandlers(robustbase::lmrob(l_alpha ~ I(l_var/l_mean^2) + 0, weights = l_weight_alpha,
        fit_meanvar_alpha <- tryCatch(
              {
                fit_meanvar_alpha = withCallingHandlers(robustbase::lmrob(alpha ~ I(var/mean^2) + 0, weights = weight_alpha,
                                                                          data = data.frame(alpha = l_alpha,
                                                                                            var=l_var_unscaled_alpha,
                                                                                            mean=l_mean_unscaled_alpha,
                                                                                            weight_alpha=l_weight_alpha),
                                                                          control=lmrob_control), warning=warning_handler)
                fit_meanvar_alpha
              },
              error=function(cond) {
                message("Error in fit_meanvar_alpha -mean-variance")
                message(cond)

                fit_meanvar_alpha = lm(alpha ~ I(var/mean^2) + 0, weights = weight_alpha,
                                                          data = data.frame(alpha = l_alpha,
                                                                            var=l_var_unscaled_alpha,
                                                                            mean=l_mean_unscaled_alpha,
                                                                            weight_alpha=l_weight_alpha))
                fit_meanvar_alpha
              },
              warning=function(cond) {
                message("Warning in fit_meanvar_alpha -mean-variance")
                message(cond)
                fit_meanvar_alpha = withCallingHandlers(robustbase::lmrob(alpha ~ I(var/mean^2) + 0, weights = weight_alpha,
                                                                          data = data.frame(alpha = l_alpha,
                                                                                            var=l_var_unscaled_alpha,
                                                                                            mean=l_mean_unscaled_alpha,
                                                                                            weight_alpha=l_weight_alpha),
                                                                          control=lmrob_control), warning=warning_handler)
                fit_meanvar_alpha
              },
              finally={
                # NOTE:
                # Here goes everything that should be executed at the end, regardless of success or error.
              }
            )

        
        # fit_invmean_alpha <- withCallingHandlers(robustbase::lmrob(l_alpha ~ I(1/l_mean) + 0, weights = l_weight_alpha,                               
        fit_invmean_alpha <- tryCatch(
          {
            fit_invmean_alpha <- withCallingHandlers(robustbase::lmrob(alpha ~ I(1/mean) + 0, weights = weight_alpha,                               
                                                                       data = data.frame(alpha = l_alpha,
                                                                                         var=l_var_unscaled_alpha,
                                                                                         mean=l_mean_unscaled_alpha,
                                                                                         weight_alpha=l_weight_alpha),
                                                                       control=lmrob_control), warning=warning_handler)
            fit_invmean_alpha
          },
          error=function(cond) {
            message("Error in fit_invmean_alpha - mean-variance")
            message(cond)
            
            fit_invmean_alpha <- lm(alpha ~ I(1/mean) + 0, weights = weight_alpha,                               
                                                                       data = data.frame(alpha = l_alpha,
                                                                                         var=l_var_unscaled_alpha,
                                                                                         mean=l_mean_unscaled_alpha,
                                                                                         weight_alpha=l_weight_alpha))
            fit_invmean_alpha
          },
          warning=function(cond) {
            message("Warning in fit_invmean_alpha - mean-variance")
            message(cond)
            
            fit_invmean_alpha <- withCallingHandlers(robustbase::lmrob(alpha ~ I(1/mean) + 0, weights = weight_alpha,                               
                                                                       data = data.frame(alpha = l_alpha,
                                                                                         var=l_var_unscaled_alpha,
                                                                                         mean=l_mean_unscaled_alpha,
                                                                                         weight_alpha=l_weight_alpha),
                                                                       control=lmrob_control), warning=warning_handler)
            fit_invmean_alpha
          },
          finally={
            # NOTE:
            # Here goes everything that should be executed at the end, regardless of success or error.
          }
        )

        # fit_combined_alpha = withCallingHandlers(robustbase::lmrob(l_alpha ~ I(l_var/l_mean^2) + I(1/l_mean) + 0, weights = l_weight_alpha,
        fit_combined_alpha <- tryCatch(
          {
            fit_combined_alpha = withCallingHandlers(robustbase::lmrob(alpha ~ I(var/mean^2) + I(1/mean) + 0, weights = weight_alpha,
                                                                       data = data.frame(alpha = l_alpha,
                                                                                         var=l_var_unscaled_alpha,
                                                                                         mean=l_mean_unscaled_alpha,
                                                                                         weight_alpha=l_weight_alpha),
                                                                       control=lmrob_control), warning=warning_handler)
            fit_combined_alpha
          },
          error=function(cond) {
            message("Error in fit_combined_alpha - mean-variance")
            message(cond)
            
            fit_combined_alpha = lm(alpha ~ I(var/mean^2) + I(1/mean) + 0, weights = weight_alpha,
                                       data = data.frame(alpha = l_alpha,
                                                         var=l_var_unscaled_alpha,
                                                         mean=l_mean_unscaled_alpha,
                                                         weight_alpha=l_weight_alpha))
            fit_combined_alpha
          },
          warning=function(cond) {
            message("Warning in fit_combined_alpha - mean-variance")
            message(cond)
            
            fit_combined_alpha = lm(alpha ~ I(var/mean^2) + I(1/mean) + 0, weights = weight_alpha,
                                    data = data.frame(alpha = l_alpha,
                                                      var=l_var_unscaled_alpha,
                                                      mean=l_mean_unscaled_alpha,
                                                      weight_alpha=l_weight_alpha))
            fit_combined_alpha
            
            # fit_combined_alpha = withCallingHandlers(robustbase::lmrob(alpha ~ I(var/mean^2) + I(1/mean) + 0, weights = weight_alpha,
            #                                                            data = data.frame(alpha = l_alpha,
            #                                                                              var=l_var_unscaled_alpha,
            #                                                                              mean=l_mean_unscaled_alpha,
            #                                                                              weight_alpha=l_weight_alpha),
            #                                                            control=lmrob_control), warning=warning_handler)
            # fit_combined_alpha
          },
          finally={
            # NOTE:
            # Here goes everything that should be executed at the end, regardless of success or error.
          }
        )
  
        alpha.rsquared_meanvar = summary(fit_meanvar_alpha)[["adj.r.squared"]]
        alpha.coefficient_meanvar = fit_meanvar_alpha$coefficients[["I(var/mean^2)"]]
        alpha.rsquared_invmean = summary(fit_invmean_alpha)[["adj.r.squared"]]
        alpha.coefficient_invmean = fit_invmean_alpha$coefficients[["I(1/mean)"]]
        
        alpha.coefficient_combined_1 = fit_combined_alpha$coefficients[["I(var/mean^2)"]]
        alpha.coefficient_combined_2 = fit_combined_alpha$coefficients[["I(1/mean)"]]
        alpha.rsquared_combined = summary(fit_combined_alpha)[["adj.r.squared"]]
        
        if(class(fit_meanvar_alpha) == "lmrob"){alpha.converged_meanvar_alpha=TRUE}
        if(class(fit_invmean_alpha) == "lmrob"){alpha.converged_invmean_alpha=TRUE}
        if(class(fit_combined_alpha) == "lmrob"){alpha.converged_combined_alpha=TRUE}
      }
      ## end alpha model ##

      # NOTE, * rpc_meanvar^2 leads to prediction in reads space, these values are on copy number scale
      var_cn2 = var.beta_cn2 / rpc_meanvar^2
      var_cn1 = var.beta_cn1 / rpc_meanvar^2
      delta_var = var_cn2 - var_cn1

      ## output ####
      df1 = dplyr::tibble(n_segs=n_segs, rpc.meanvar = rpc.meanvar, ploidy.meanvar=ploidy.meanvar,
                          omega=beta.coefficient_combined_1, alpha=beta.coefficient_negbin,
                          alpha_nb=beta.coefficient_negbin, alpha_gam=gam.alpha, zero_offset=mu_zero)
      
      df2 = dplyr::tibble(
        mu_zero=mu_zero, mu_zero_ratio=mu_zero_ratio, n_zero=n_zero,
        rpc.meanvar=rpc.meanvar, ploidy.meanvar=ploidy.meanvar, n_seg_effective=n_seg_effective,
        rpc_meanvar=rpc_meanvar, rpc_meanvar2=rpc_meanvar2, 
        
        ploidy_meanvar=ploidy_meanvar, ploidy_meanvar2=ploidy_meanvar2, ploidy_meanvar3=ploidy_meanvar3,
        n_segs=n_segs, n_segs_na=n_segs_na,

        var.beta_rsquared_unscaled=var.beta_rsquared_unscaled,
        var.beta_cn4=var.beta_cn4, var.beta_cn2=var.beta_cn2, var.beta_cn1=var.beta_cn1,
        var.meanvar_ratio_4=var.meanvar_ratio_4, var.meanvar_ratio_2=var.meanvar_ratio_2, var.meanvar_ratio_1=var.meanvar_ratio_1,
        var.coefficient_beta_unscaled_1=var.coefficient_beta_unscaled_1, var.coefficient_beta_unscaled_2=var.coefficient_beta_unscaled_2,
        beta.rweights_unscaled_beta=beta.rweights_unscaled_beta,

        beta.rsquared_negbin=beta.rsquared_negbin, beta.coefficient_negbin=beta.coefficient_negbin,
        beta.converged_negbin=beta.converged_negbin, beta.quality_negbin=beta.quality_negbin,

        beta.converged_combined=beta.converged_combined, beta.coefficient_combined_1=beta.coefficient_combined_1,
        beta.coefficient_combined_2=beta.coefficient_combined_2, beta.quality_combined=beta.quality_combined,

        gamlss.chi=gamlss.chi, gamlss.pvalue=gamlss.pvalue,
        gamlss.mu_coefficient_0=gamlss.mu_coefficient_0, gamlss.mu_coefficient_1=gamlss.mu_coefficient_1,
        gamlss.sigma_coefficient_0=gamlss.sigma_coefficient_0, gamlss.sigma_coefficient_1=gamlss.sigma_coefficient_1,

        gam.alpha.devexpl=gam.alpha.devexpl, gam.alpha=gam.alpha,
        gam.alpha.coefficient_0=gam.alpha.coefficient_0, gam.alpha.coefficient_1=gam.alpha.coefficient_1,

        gam.omega.devexpl=gam.omega.devexpl, gam.omega=gam.omega,
        gam.omega.coefficient_0=gam.omega.coefficient_0, gam.omega.coefficient_1=gam.omega.coefficient_1,

        gam.filter_beta_alpha.devexpl = gam.filter_beta_alpha.devexpl, gam.filter_beta_alpha = gam.filter_beta_alpha,
        gam.filter_beta_alpha.coefficient_0 = gam.filter_beta_alpha.coefficient_0,
        gam.filter_beta_alpha.coefficient_1 = gam.filter_beta_alpha.coefficient_1,
        gam.filter_negbin_alpha.devexpl = gam.filter_negbin_alpha.devexpl, gam.filter_negbin_alpha = gam.filter_negbin_alpha,
        gam.filter_negbin_alpha.coefficient_0 = gam.filter_negbin_alpha.coefficient_0,
        gam.filter_negbin_alpha.coefficient_1 = gam.filter_negbin_alpha.coefficient_1,

        alpha.n_seg=alpha.n_seg, alpha.percent_used=alpha.percent_used,
        alpha.rsquared_meanvar=alpha.rsquared_meanvar, alpha.coefficient_meanvar=alpha.coefficient_meanvar,
        alpha.rsquared_invmean=alpha.rsquared_invmean, alpha.coefficient_invmean=alpha.coefficient_invmean,
        alpha.coefficient_combined_1=alpha.coefficient_combined_1,
        alpha.coefficient_combined_2=alpha.coefficient_combined_2,
        alpha.rsquared_combined=alpha.rsquared_combined,

        alpha.converged_meanvar_alpha=alpha.converged_meanvar_alpha,
        alpha.converged_invmean_alpha=alpha.converged_invmean_alpha,
        alpha.converged_combined_alpha=alpha.converged_combined_alpha,

        var_cn2=var_cn2, var_cn1=var_cn1,
        delta_var=delta_var)

      colnames(df2) = paste0("meanvar.", colnames(df2))
      df = dplyr::bind_cols(df1, df2)
      
      # # # # ## DEBUG ####
      if(debug){
        require(ggplot2, quietly=TRUE, warn.conflicts = FALSE)
        # filename = paste0("~/debug", ".", slope, ".pdf")
        p = ggplot(data=data.frame(x=l_mean, y=l_var, z = l_weight)) + geom_point(aes(x=x,y=y,size=z), alpha=0.3) +
          theme_minimal() + xlab("mean") + ylab("variance") +
          geom_smooth(aes(x=x, y=y, weight=z), method="lm", formula=y ~ I(x) + I(x^2) + 0) 
          # geom_abline(intercept = 0, slope = beta, color="#E69F00") +
          # geom_abline(intercept = 0, slope = beta_linear, color="#A3C1AD") +
          # ggtitle(paste0("Two: ", mean(counts), "-", var(counts))) #, "-", beta_linear, "-", beta_robust))

        # ggsave(filename, plot=p, width = 9, height = 4)
        # keep indeces are required to save to segtable (plus used in cellcycle and in visualization)

        # identify chromosome on which a segment is located
        chrom_index = names(counts)
        rl = RleList(lapply(names(counts), function(x){a = counts[[x]]; return(Rle(rep(match(x, chrom_index), sum(a@lengths))))}))
        names(rl) = chrom_index
        l_chromosome = smoothOperator(rl, gr, function(x){x = unique(x); return(x)}, trimLength = trimLength, minLength=minLength)
        l_chromosome = chrom_index[l_chromosome]
        l_chromosome_alpha = l_chromosome[keep_ix_2]
        l_chromosome = l_chromosome[keep_ix_1]
        
        return(list(df=df, fit=fit_unscaled_beta, gr=gr, plot=p, keep_index_alpha=keep_ix_2, keep_index_beta=keep_ix_1,
                    l_mean=l_mean, l_var=l_var, l_weight=l_weight, l_chromosome=l_chromosome,
                    l_mean_unscaled=l_mean_unscaled, l_var_unscaled=l_var_unscaled,
                    lw_mean=lw_mean, lw_counts=lw_counts,
                    l_alpha=l_alpha, l_omega=l_omega, l_mean_alpha=l_mean_alpha, l_var_alpha=l_var_alpha,
                    l_weight_alpha=l_weight_alpha, l_chromosome_alpha=l_chromosome_alpha,
                    l_mean_unscaled_alpha=l_mean_unscaled_alpha, l_var_unscaled_alpha=l_var_unscaled_alpha,
                    l_mean_squared=l_mean_squared, l_mean_squared_alpha=l_mean_squared_alpha
                    ))
      }
      
      return(list(df=df, fit=fit_unscaled_beta))
    
  } else {
    # SCAFFOLD OUTPUT - required for failed instances
    df1 = dplyr::tibble(n_segs=NA, rpc.meanvar = NA, ploidy.meanvar=NA,
                        omega=NA, alpha=NA,
                        alpha_nb=NA, alpha_gam=NA, zero_offset=NA)
    
    df2 = dplyr::tibble(
      mu_zero=NA,mu_zero_ratio=NA,n_zero=NA,
      rpc.meanvar=NA, ploidy.meanvar=NA, n_seg_effective=NA,
      
      rpc_meanvar=NA, rpc_meanvar2=NA,
      ploidy_meanvar=NA, ploidy_meanvar2=NA, ploidy_meanvar3=NA,
      n_segs=NA, n_segs_na=NA,
      
      var.beta_rsquared_unscaled=NA,
      var.beta_cn4=NA, var.beta_cn2=NA, var.beta_cn1=NA,
      var.meanvar_ratio_4=NA, var.meanvar_ratio_2=NA, var.meanvar_ratio_1=NA,
      var.coefficient_beta_unscaled_1=NA, var.coefficient_beta_unscaled_2=NA,
      beta.rweights_unscaled_beta=NA,
      
      beta.rsquared_negbin=NA, beta.coefficient_negbin=NA, beta.converged_negbin=NA, beta.quality_negbin=NA,
      
      beta.converged_combined=NA, beta.coefficient_combined_1=NA, 
      beta.coefficient_combined_2=NA, beta.quality_combined=NA,
      
      gamlss.chi=NA, gamlss.pvalue=NA, gamlss.mu_coefficient_0=NA, gamlss.mu_coefficient_1=NA,
      gamlss.sigma_coefficient_0=NA, gamlss.sigma_coefficient_1=NA,

      gam.alpha.devexpl=NA, gam.alpha=NA, gam.alpha.coefficient_0=NA, gam.alpha.coefficient_1=NA,

      gam.omega.devexpl=NA, gam.omega=NA, gam.omega.coefficient_0=NA, gam.omega.coefficient_1=NA,

      gam.filter_beta_alpha.devexpl = NA, gam.filter_beta_alpha = NA,
      gam.filter_beta_alpha.coefficient_0 = NA, gam.filter_beta_alpha.coefficient_1 = NA,
      gam.filter_negbin_alpha.devexpl = NA, gam.filter_negbin_alpha = NA,
      gam.filter_negbin_alpha.coefficient_0 = NA, gam.filter_negbin_alpha.coefficient_1 = NA,

      alpha.n_seg=NA, alpha.percent_used=NA, alpha.rsquared_meanvar=NA, alpha.coefficient_meanvar=NA,
      alpha.rsquared_invmean=NA, alpha.coefficient_invmean=NA, alpha.coefficient_combined_1=NA,
      alpha.coefficient_combined_2=NA, alpha.rsquared_combined=NA,
      alpha.converged_meanvar_alpha=NA, alpha.converged_invmean_alpha=NA, alpha.converged_combined_alpha=NA,
      
      var_cn2=NA, var_cn1=NA, delta_var=NA)
    
    colnames(df2) = paste0("meanvar.", colnames(df2))
    df = dplyr::bind_cols(df1, df2)

    return(df)
  }
}


#' computeBasicModel
#'
#' \code{computeBasicModel} compute alpha and rpc for a given segmentation, basic version of computeModel method
#'
#' @param counts read counts (as taken from calls or copynumber slot)
#' @param gr GRanges object, representing a segmentation of the genome
#' @param limitPloidy Numeric cutoff values for ploidy values to include in fitting
#' @param scale Numeric scale used to transform read counts
#' @param minLength Numeric minimum segment size to include in fitting
#' @param trimLength Numeric number of bins to remove from start and end of segments for estimates

#' 
#' NOTE:
#' A) counts is raw count data
#' B) scale is transform, s.t. raw * scale = (scaled copy number counts)
#' C) segmentation can be a) rounded scaled data or b) not rounded scaled data (scaled vs not scaled only matters for rounding)
#'    Rounding is only relevant for single cell DNAseq data, where expected levels can only occur at integer states.
#' 
#' @return data frame containing the model parameters
computeBasicModel <- function(counts, gr, limitPloidy, scale=1.0, minLength=10, trimLength=1, roundInteger=FALSE){
  
    ## Initialize variables ####
    suppressPackageStartupMessages({
      require(MASS, quietly = TRUE, warn.conflicts = FALSE)
      require(gamlss, quietly = TRUE, warn.conflicts = FALSE)
      require(mgcv, quietly = TRUE, warn.conflicts = FALSE)
      require(robustbase, quietly = TRUE, warn.conflicts = FALSE)
    })
    lmrob_control = robustbase::lmrob.control(setting="KS2014",
                                              max.it = 10000, maxit.scale = 10000,
                                              scale.tol=1e-5, solve.tol = 1e-5,
                                              subsampling = "nonsingular")
    warning_handler <- function(w){ 
      
      if(length(w) > 1 && any(grepl("Error", w[[1]]))){
        stop()
      }
      
      if( any( grepl( "Detected possible local breakdown of SMDM-estimate", w) | 
               grepl("find_scale\\(\\) did not converge in", w) | 
               grepl("Convergence failure: singular convergence", w)) ){
        invokeRestart( "muffleWarning" )}
    }
    
    # # remove extreme outliers, 99% quantile -> ISSUE WITH OVERDISPERSION VALUE
    # xy = lapply(ccounts, function(x) return(decode(x)))
    # xy_cutoff = quantile(unlist(xy), 0.99, na.rm=TRUE)
    # counts = RleList(lapply(ccounts, function(x){y=decode(x); y[!is.na(y) & y >= xy_cutoff] = NA; return(Rle(y))}))
    
    # compute per segments estimates
    l_mean = smoothOperator(counts * scale, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
    ## Single cell property -> means at integer state
    if(roundInteger){
      l_mean = base::round(l_mean, digits = 0)
    }
    l_var = smoothOperator(counts * scale, gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
    l_weight = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = trimLength, minLength=minLength)
    l_weight_nseg = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = 0, minLength=1)
    
    l_var_unscaled = smoothOperator(counts, gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
    l_mean_unscaled = smoothOperator(counts, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
    l_mean_squared = smoothOperator(counts^2, gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=minLength)
    
    lw_counts = smoothOperator(counts, gr, function(x){return(x[!is.na(x)])}, trimLength = 0, minLength=1)
    lw_mean = smoothOperator(counts, gr, function(x){return(mean(x, na.rm=TRUE))}, trimLength = 0, minLength=1)
    lw_mean = lw_mean[!is.na(lw_mean)]
    lw_mean = rep(lw_mean, times=l_weight_nseg[l_weight_nseg != 0]) * scale
    stopifnot(length(lw_counts) == length(lw_mean))
    
    # NOTE compute values for rpc and ploidy, more versions are computed below
    zero_elements = lw_mean < 0.5
    rpc.meanvar = mean(lw_counts[!zero_elements] / round(lw_mean[!zero_elements], digits=0), na.rm=TRUE)

    
    #############################################
    # estimates that large are very unreliable and would unduly bias the overall beta estimate
    # -> we can use limitPloidy here as l_mean should have been scaled!
    # stopifnot(length(l_alpha) == length(l_mean))
    n_seg_effective = length(l_mean)
    keep_ix_1 = (l_mean <= limitPloidy & l_mean >= 0.5 & !is.na(l_var) & l_var > 0 & l_weight >= (minLength - 2*trimLength))
    # keep_ix_2 = !is.na(l_alpha) & !is.na(l_omega) & keep_ix_1
    
    l_mean = l_mean[keep_ix_1]
    l_var = l_var[keep_ix_1]
    l_var_unscaled = l_var_unscaled[keep_ix_1]
    l_mean_unscaled = l_mean_unscaled[keep_ix_1]
    l_weight = l_weight[keep_ix_1]
    l_mean_squared = l_mean_squared[keep_ix_1]
    
    if(length(l_mean) == 0){
      return(NULL)
    }
    
    rpc_meanvar = sum(lw_counts, na.rm=TRUE) / sum(round(lw_mean, digits = 0))
    n_segs = sum(l_weight_nseg != 0)

    #############################################
    stopifnot(all(l_weight >= minLength - 2*trimLength))
    
    ## Beta model ####
    
    ## beta model for overdispersion in NB ----
    fit_negbin_beta <- tryCatch(
      {
        fit_negbin_beta = withCallingHandlers(robustbase::lmrob(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                                                                data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight)), warning=warning_handler)
        fit_negbin_beta
      },
      error=function(cond) {
        message("Error in fit_negbin_beta - mean-variance")
        message(cond)
        fit_negbin_beta = lm(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                             data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight))
        
        fit_negbin_beta
      },
      warning=function(cond) {
        message("Warning in fit_negbin_beta - mean-variance")
        message(cond)
        fit_negbin_beta = withCallingHandlers(robustbase::lmrob(var - I(mean) ~ I(mean^2) + 0, weights=weight,
                                                                data = data.frame(var=l_var_unscaled, mean=l_mean_unscaled, weight=l_weight)), warning=warning_handler)
        fit_negbin_beta
      },
      finally={
        # NOTE:
        # Here goes everything that should be executed at the end, regardless of success or error.
      }
    )
    beta.rsquared_negbin = summary(fit_negbin_beta)[["adj.r.squared"]]
    beta.coefficient_negbin = fit_negbin_beta$coefficients[["I(mean^2)"]]
    # NOTE this is alpha in the strict NB sense
    if(class(fit_negbin_beta) == "lmrob"){
      beta.converged_negbin = TRUE
      beta.rweights_negbin = fit_negbin_beta[["rweights"]]
    }else{
      beta.converged_negbin = FALSE
      beta.rweights_negbin = rep(1.0, length(l_mean_unscaled))
    }
    beta.quality_negbin = sum(beta.rweights_negbin * l_weight) / sum(rep(1.0, length(l_mean_unscaled)) * l_weight)
    
    
    
    ## output ####
    df = dplyr::tibble(n_segs=n_segs, rpc.meanvar = rpc.meanvar,
                        alpha=beta.coefficient_negbin)
    
    return(list(df=df))
}


## Doesn't work very well unfortunately
#' #' updateSegmentation
#' #'
#' #' \code{updateSegmentation} update a given segmentation gr, by replacing gr_select (in gr) with a new segmentation gr_replace.
#' #' gr_replace introduces an additional breakpoint in gr_select, computed using update_test.stat
#' #'
#' #' @param gr GRanges object, full segmentation of genome represented by GRanges object
#' #' @param gr_select GRanges object representing a single segment that is be split in two
#' #' @param counts_to_resegment vector of read values to resegment (can contain NaNs)
#' #' @param update_test.stat Character likelihood to use for finding single new breakpoint in sequence (supports Negative Binomial, Poisson and Normal)
#' #' 
#' #' @returna list with a) a new full GRanges object containing the updated segmentation of the full genome (gr)
#' #' b) the two new segments replacing the gr_select segment as a GRanges object (gr_replacement) and 
#' #' c) the old full genome segmentation (gr_old)
#' updateSegmentation <- function(gr, gr_select, counts_to_resegment, update_test.stat="NegBin"){
#'   require(changepoint, quietly = TRUE, warn.conflicts = FALSE)
#'   stopifnot(update_test.stat %in% c("NegBin", "Poisson", "Normal"))
#'   
#'   na_index = is.na(counts_to_resegment)
#'   if(update_test.stat == "NegBin"){
#'     stopifnot(length(counts_to_resegment[!na_index]) >= 10)
#'     out = segmentSingleBreakpoint(as.integer(counts_to_resegment[!na_index]))
#'     if(out[["null"]] >= out[["alt"]]){
#'       # warning("NULL hypothesis rejected")
#'       return(list(gr_new=gr, gr_replacement=gr_select, gr_old=gr))
#'     }
#'     cpt = out[["cpt"]]
#'   }else{
#'     cpt = changepoint::cpt.meanvar(as.integer(counts_to_resegment[!na_index]), penalty="BIC", pen.value=NULL, method="AMOC",
#'                                    test.stat=update_test.stat, class=FALSE, param.estimates=FALSE, minseglen=2)    
#'   }
#' 
#'   # NOTE necessary to add nan position to cpt value (have to be excluded to compute changepoint)
#'   if(update_test.stat == "Normal"){
#'     cpt = as.numeric(cpt["cpt"]) + cumsum(na_index)[as.numeric(cpt["cpt"])]
#'   }else{
#'     cpt = as.numeric(cpt) + cumsum(na_index)[as.numeric(cpt)]
#'   }
#'   
#'   stopifnot(length(cpt) == 1)
#'   
#'   ## Add new ranges
#'   iranges = gr_select@ranges
#'   start1 = iranges@start
#'   end2 = iranges@start + iranges@width - 1
#'   end1 = start1 + cpt - 1
#'   start2 = start1 + cpt
#'   width1 = end1 - start1 + 1
#'   width2 = end2 - start2 + 1
#' 
#'   median.reads = smoothOperator(as.numeric(counts_to_resegment), IRanges(start=c(1, cpt), end=c(cpt-1, length(counts_to_resegment))), 
#'                                 function(x){stats::median(x, na.rm=TRUE)}, trimLength = 0, minLength=0)
#'   gr_replace = GRanges(seqnames(gr_select),
#'                    IRanges(start=c(start1, start2), end=c(end1, end2), width=c(width1, width2)))
#'   gr_replace$segmentation = median.reads
#'   
#'   # some sanity checks
#'   # print("DEBUG")
#'   # print(gr_select)
#'   # print("----")
#'   # print(gr_replace)
#'   # print("-------------------")
#'   stopifnot(gr_select@ranges@start[1] == gr_replace@ranges@start[1])
#'   stopifnot((gr_select@ranges@start[length(gr_select)] + gr_select@ranges@width[length(gr_select)] - 1) == 
#'               (gr_replace@ranges@start[length(gr_replace)] + gr_replace@ranges@width[length(gr_replace)] - 1))
#'   stopifnot(sum(gr_replace@ranges@width) == sum(gr_select@ranges@width))
#'   stopifnot(runValue(seqnames(gr_select)) == runValue(seqnames(gr_replace)))
#'   
#'   
#'   ## Manipulate gr object, remove gr_select and insert gr_replace
#'   grs = length(gr@ranges@width)
#'   gr_old = gr
#' 
#'   gr[gr == gr_select] = NULL
#'   gr_new = c(gr, gr_replace)
#'   gr_new <- sortSeqlevels(gr_new)
#'   gr_new <- sort(gr_new)
#'   
#'   stopifnot(grs - length(gr_select@ranges@width) + length(gr_replace@ranges@width) == length(gr_new@ranges@width))
#'   
#'   return(list(gr_new=gr_new, gr_replacement=gr_replace, gr_old=gr_old))
#' }
#' 
#' #' optimizeSegmentation
#' #'
#' #' \code{optimizeSegmentation} iterative algorithm to optimize segmentation
#' #'
#' #' @param object QDNAseq object with segmented counts (segmented slot)
#' #' @param trimLength Numeric number of bins to remove from start and end of segments for estimates
#' #' @param minLength Numeric minimum segment size to include in fitting
#' #' @param n_iterations Numeric number of iterations of optimisation procedure
#' #' @param delta_rsquared Numeric Citerion to keep a novel breakpoint: Only if the overall rsquared is improved by at least delta_rsquared
#' #' # in regard to the robustbase::lmrob fit the new breakpoint is accepted
#' #' @param rweight_criteria Numeric stopping criteria for single round of optimisation procedure, 
#' #' In one round, rweight < rweight_criteria elements are broken up (rweight is determined through robustbase::lmrob method)
#' #' @param update_test.stat Character test statistic to use to split individual segments (supported: Poisson (default) and Normal)
#' #'
#' #' @return QDNAseq object with updated segmentation (and a number of metadata detailing model changes (prefix optim))
#' optimizeSegmentation <- function(object,
#'                                  trimLength=1, minLength=30,
#'                                  n_iterations=10, delta_rsquared=1e-4,
#'                                  rweight_criteria=0.9, update_test.stat="NegBin"){
#' 
#'   require(robustbase, quietly = TRUE, warn.conflicts = FALSE)
#'   require(IRanges, quietly = TRUE, warn.conflicts = FALSE)
#'   require(changepoint, quietly = TRUE, warn.conflicts = FALSE)
#'   
#'   # valid = binsToUseInternal(object)
#'   n_samples = dim(object)[[2]]
#'   n_bins = dim(object)[[1]]
#'   
#'   optimize_segmentation_call <- function(cell){
#'     
#'     ## initial segmentation (per chromosome split)
#'     gr = retrieveSegmentation(object, cell)
#'     counts = retrieveAssayData(object, cell, "copynumber")
#'     
#'     # this segmentation needs to be run before scaling -> do not scale, i.e. scale is at default of 1.0
#'     initial_model = computeModel(counts, gr, limitPloidy=Inf, minLength=minLength, trimLength=trimLength, debug=TRUE, skip_alpha = TRUE)
#'     initial_gr = gr
#'     
#'     model_fit = initial_model
#'     rsquared = model_fit$df$beta_rsquared_linear
#'     stopifnot(class(model_fit$fit) == "lmrob")
#'     rweights = model_fit$fit$rweights
#'     previous_rsquared = 0.0
#'     counter = 0
#'     
#'     # # DEBUG
#'     # debug_list = list()
#'     
#'     for(j in 1:n_iterations){
#'       # Needs to be reset in each iteration loop, when iterating over all segments
#'       broken_gr = GRanges()
#'       n_segs = sum(gr@ranges@width > minLength)
#'       
#'       # print(paste0("j: ", j))
#'       # flush(stdout())
#'       # Sys.sleep(1)
#'       
#'       for(i in 1:sum(rweights < rweight_criteria)){
#'         # print(i)
#'         
#'         stopifnot(sum(gr@ranges@width) == n_bins)
#'         
#'         residuals = residuals(model_fit$fit)
#'         rweights = model_fit$fit$rweights
#'    
#'         ## Select gr based on 3 criteria
#'         # 1. minLength
#'         # 2. not in previously broken ir region
#'         # 3. highest residual (potentially weighted, see above)
#'         # gr_sel = gr[gr@ranges@width >= minLength & !is.na(gr$segmentation)]
#'         gr_sel = gr[gr@ranges@width >= minLength][model_fit$keep_index_beta]
#'         stopifnot(all(!is.na(gr_sel)))
#'         stopifnot(length(rweights) == length(gr_sel))
#'         rws = rweights[(gr_sel %outside% broken_gr)]
#'         rws_weights = model_fit$l_weight[(gr_sel %outside% broken_gr)]
#'         stopifnot(length(rws) == length(rws_weights))
#'         gr_select = gr_sel[gr_sel %outside% broken_gr]
#'         stopifnot(length(rws) == length(gr_select@ranges@width))
#'         
#'         # Prioritze segments that do not fit (rweights small or 0)
#'         
#'         # select single item with lowest rws
#'         # sorted = sort(rws, decreasing = FALSE, index.return=TRUE)
#'         # segment_to_break_index = sorted$ix[1]
#'         
#'         # select random item with probability iverse to rweight
#'         prob_vector = ifelse(rws < rweight_criteria, (1-rws)*rws_weights, 0)
#'         if(sum(prob_vector) < 1e-6){
#'           # case of perfect fit, i.e. all rws are 1
#'           break;
#'         }
#'         segment_to_break_index = sample.int(length(rws), size=1, prob=prob_vector)
#'         
#'         # alternatively possible to use the residuals
#'         # sorted = sort(abs(residuals)[gr@ranges@width > minLength], decreasing = TRUE, index.return=TRUE)
#'         
#'         # Choose next segment to break up
#'         gr_select = gr_select[segment_to_break_index]
#'         broken_gr = c(broken_gr, gr_select)
#' 
#'         # Raw data in segments that needs to be further split up -> should always be on single chromosome
#'         counts_to_resegment = decode(counts[[seqnames(gr_select)]])[gr_select@ranges@start:(gr_select@ranges@start + gr_select@ranges@width - 1)]
#' 
#'         if((length(unique(counts_to_resegment)) < 3) || length(counts_to_resegment[!is.na(counts_to_resegment)]) <= minLength){
#'           # fix case of idential values (e.g. zeros) and require minimum of 4 observations to find new breakpoint
#'           break;
#'         }
#'         
#'         # Update segmentation -> break gr_select region into two pieces
#'         # print(sample)
#'         grlist = updateSegmentation(gr, gr_select, counts_to_resegment, update_test.stat=update_test.stat)
#'         if(all(grlist$gr_replacement@ranges@width == gr_select@ranges@width)){
#'           next
#'         }
#'         gr = grlist$gr_new
#'         gr_old = grlist$gr_old
#'         gr_replacement = grlist$gr_replacement
#'         
#'         old_rsquared_linear = model_fit$df$beta_rsquared_linear
#'         old_rsquared_robust = model_fit$df$beta_rsquared_robust
#'         model_fit = computeModel(counts, gr, limitPloidy=Inf, minLength=minLength, trimLength=trimLength, debug=TRUE, skip_alpha=TRUE)
#'         
#'         target_residual = residuals[segment_to_break_index]
#'         residuals = residuals(model_fit$fit)
#'         stopifnot(length(gr_replacement@ranges@width) == 2)
#'         novel_residuals = residuals[segment_to_break_index:(segment_to_break_index+1)]
#'         size = sum(gr_replacement@ranges@width)
#'         w = gr_replacement@ranges@width / size
#'         relative_size = size / sum(gr@ranges@width)
#'         
#'         novel_rsquared_linear = model_fit$df$beta_rsquared_linear
#'         novel_rsquared_robust = model_fit$df$beta_rsquared_robust
#'         
#'         # CRITERION for accepting a further breakup of a segment
#'         CRITERION1 = !any(is.na(novel_residuals)) && all(gr_replacement@ranges@width >= minLength)
#'         CRITERION2 = (novel_rsquared_robust - old_rsquared_robust) > delta_rsquared
#'         
#'         if(is.na(CRITERION2)){CRITERION2=FALSE}
#'         
#'         # # DEBUG
#'         # print(paste("CRITERIA: ", CRITERION1, CRITERION2, CRITERION3, sep=" | "))
#'         # debug_list[["criteria1"]] = c(debug_list[["criteria1"]], CRITERION1)
#'         # debug_list[["criteria2"]] = c(debug_list[["criteria2"]], CRITERION2)
#'         # debug_list[["criteria3"]] = c(debug_list[["criteria3"]], CRITERION3)
#'         
#'         CRITERION = CRITERION1 && CRITERION2
#'         if(!CRITERION){
#'           # if no improvement, do not update the segmentation and keep old segmentation, i.e. gr_old
#'           gr = gr_old
#'           model_fit = computeModel(counts, gr, limitPloidy=Inf, minLength=minLength, trimLength=trimLength, debug=TRUE, skip_alpha=TRUE)
#'           previous_rsquared = model_fit$df$beta_rsquared_robust
#'         } else {
#'           counter = counter + 1
#'           # counter is a direct measure of changes to model (some segments might disappear)
#'           previous_rsquared = model_fit$df$beta_rsquared_robust
#'         }
#'         
#'         # end:    for(i in 1:n_segs){
#'       }
#'       
#'       # end:    for(j in 1:n_iterations){
#'     }
#'     
#'     # # DEBUG
#'     # a = debug_list[["criteria1"]]
#'     # b = debug_list[["criteria2"]]
#'     # c = debug_list[["criteria3"]]
#'     # print(paste0("Crosstable: ", sample))
#'     # print(table(b, c, a))
#'     
#'     initial_rsquared_robust = initial_model$df$beta_rsquared_robust
#'     final_rsquared_robust = model_fit$df$beta_rsquared_robust
#'     delta_beta_rsquared_robust = final_rsquared_robust - initial_rsquared_robust
#'     
#'     initial_rsquared_linear = initial_model$df$beta_rsquared_linear
#'     final_rsquared_linear = model_fit$df$beta_rsquared_linear
#'     delta_beta_rsquared_linear = final_rsquared_linear - initial_rsquared_linear
#'     
#'     initial_n_seg = initial_model$df$beta_n_seg
#'     final_n_seg = model_fit$df$beta_n_seg
#'     delta_n_seg = final_n_seg - initial_n_seg
#'     
#'     initial_slope_linear = initial_model$df$beta_linear
#'     final_slope_linear = model_fit$df$beta_linear
#'     delta_beta_linear = final_slope_linear - initial_slope_linear
#'     
#'     initial_slope_robust = initial_model$df$beta_robust
#'     final_slope_robust = model_fit$df$beta_robust
#'     delta_beta_robust = final_slope_robust - initial_slope_robust
#'     
#'     # obtain new per sample segmentations
#'     medians = smoothOperator(counts, gr, function(x){stats::median(x, na.rm=TRUE)}, trimLength = 0, minLength=0)
#'     new_segments = rep(medians, times=gr@ranges@width)
#'     
#'     beta = model_fit$df$beta
#'     beta_linear = model_fit$df$beta_linear
#'     beta_robust = model_fit$df$beta_robust
#'     beta_linear_initial = initial_model$df$beta_linear
#'     beta_robust_initial = initial_model$df$beta_robust
#' 
#'     beta_rsquared_robust = model_fit$df$beta_rsquared_robust
#'     beta_rsquared_linear = model_fit$df$beta_rsquared_linear
#'     beta_n_seg = model_fit$df$beta_n_seg
#'     beta_percent_used = model_fit$df$beta_percent_used
#'     
#'     return(list(new_segments=new_segments,
#'                 df=dplyr::tibble(delta_beta_rsquared_linear=delta_beta_rsquared_linear, delta_beta_rsquared_robust=delta_beta_rsquared_robust, 
#'                               delta_nseg=delta_n_seg, delta_changes = counter,
#'                               delta_beta_linear=delta_beta_linear, delta_beta_robust=delta_beta_robust,
#'                               beta=beta, beta_linear=beta_linear, beta_robust=beta_robust,
#'                               beta_linear_initial=beta_linear_initial, beta_robust_initial=beta_robust_initial,
#'                               beta_rsquared_robust=beta_rsquared_robust, beta_rsquared_linear=beta_rsquared_linear,
#'                               beta_n_seg=beta_n_seg, beta_percent_used=beta_percent_used)))
#'   }
#'   # end per sample part
#'   
#'   results = future.apply::future_lapply(1:n_samples, optimize_segmentation_call)
#'   
#'   # update object with new information
#'   new_segments = do.call(cbind, lapply(results, `[[`, 1))
#'   dfs = do.call("rbind", lapply(results, `[[`, 2))
#'   
#'   # safe metadata
#'   colnames(dfs) = paste0("optim.", colnames(dfs))
#'   description = Biobase::pData(object)
#'   cnames = rownames(description)
#'   description = dplyr::bind_cols(description, dfs)
#'   rownames(description) = cnames
#'   Biobase::pData(object) = description
#'   
#'   segments = Biobase::assayDataElement(object, "segmented")
#'   stopifnot(dim(segments) == dim(new_segments))
#'   
#'   object = Biobase::assayDataElementReplace(object, "segmented", new_segments)
#'   return(object)
#' }
