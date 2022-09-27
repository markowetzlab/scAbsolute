# Copyright 2022, Michael Schneider, All rights reserved.
#' @details
#' Please refer to the vignette to see how to use the package.
#' @keywords internal
"_PACKAGE"

#' scAbsolute: A package for estimating absolute copy number state in single cell DNA sequencing data.
#'
#' @docType package
#' @name scAbsolute
#' @author Michael P Schneider, \email{michael.schneider@cruk.cam.ac.uk}
#' @keywords single cell, absolute scaling, DNA sequencing
#' @export
#' @section The scAbsolute package provides the following functions:
#' - scAbsolute
#' - computeScale
#' - selectSolution
#' 
#' - segment
#' 
#' - binsToUseInternal
#' - createGR
#' - retrieveAssayData
#' - retrieveSegmentation
#' - selectChromosomes
#' - applyScale
#' - computeRPC
#' - computeGini
#' - computeHigherOrderMoments
#' - addObservedVariance
#' - readData
#' - readFilter
#' - readFlagstat
#' - combineQDNASets
#' - getSegTable
#' - writeSegTable
#' 
#' - cellcycleMetadata
#' 
#' - computeModel
#' - updateSegmentation
#' - optimizeSegmentation
#' 
#' - plotCounts
#' - plotQDNAseq
#' - plotCopynumber
#' - plotFit
#' - plotDataset


#' computeScale
#' 
#' identify a scaling factor transforming a relative read per bin scale to an absolute copy number scale
#'
#' \code{computeScale} infers an optimal scaling factor per input sample and returns the infered fit using a python implementation of 
#' stochastic variational inference with a Dirichlet Process Gaussian Mixture Model.
#'
#' @param segmentedCounts QDNAseq object with segmented slot
#' @param truncation Numeric Truncation value for Variational Inference of DPGMM
#' @param batchSize Numeric How many cells to process in parallel (only relevant for local computation - we advise to parallelise on cluster)
#' @param n_init Numeric Number of random initialisations for variational inference algorithm (3 is generally a good default)
#' @param n_steps Numeric Number of epochs across data set per initialisation. We find 10 epochs sufficient for 5 million reads and 30kb resolution. Increase if the fits appear imprecise.
#' @param randomSeed Numeric Random seed to use for random initialisation at start of SVI algorithm.
#' @param verbosity Debugging output from python. Default 2
#' 
#' @return fit object containing best fit among the iterations from the SVI algorithm
#'   
#' @export
computeScale <- function(segmentedCounts, truncation = 64, batchSize=1, n_init = 3, n_steps = 10, randomSeed=2018, verbosity=2){

  ## Access the python code that does the heavy lifting
  require(reticulate, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  stopifnot(class(segmentedCounts)[[1]] == "QDNAseqCopyNumbers")

  data = t(Biobase::assayDataElement(segmentedCounts, "segmented"))
  reticulate::source_python(file.path(BASEDIR, "data/scAbsolute.py"), convert=TRUE)
  
  # remove NaN columns
  complete_data = data[,colSums(is.na(data))==0]
  if (!is.null(dim(complete_data))){
    input = as.matrix(complete_data)
  }else{
    input = as.matrix(complete_data)
    input = t(unname(input))
  }

  result <- run_scAbsolute(data=input, truncation = as.integer(truncation), n_steps = as.integer(n_steps), 
                           n_init = as.integer(n_init), batch_size = as.integer(batchSize),
                           random_state=as.integer(randomSeed), verbose=as.integer(verbosity))
  
  solution <- result[[1]]
  fit <- result[[2:length(result)]]
  
  # unpack data frame to more convenient format
  unpack <- function(fit){
    n = length(fit$xi)
    liste = list()
    for(i in 1:n){
      item = list(xi=fit$xi[[i]], prescale=fit$prescale[[i]], scale=fit$scale[[i]], 
               weights=fit$weights[i, ], means=fit$means[i, ], covariance=fit$covariance[i, ])
      liste[[i]] = item
    }
    return(liste)
  }
  fit <- unpack(fit)

  description <- Biobase::pData(segmentedCounts)
  description$scale = as.numeric(unname(solution$scale))
  Biobase::pData(segmentedCounts) <- description

  return(fit)
}

#' evaluateScalings
#'
#' \code{evaluateScalings} helper function to evaluate a given discrete scaling solution in terms of error or beta
#'
#' The parameters are directly passed through the selectSolution function. See this function for more details on the parameters.
#' The evaluateScalings method work on invidividual cells and makes it possible to parallelize across multiple cells.
evaluateScalings <- function(segmCN, fiti, cellname, 
                             component_cutoff, weight_cutoff,
                             minPloidy, maxPloidy, ploidyRegion,
                             trimLength, minLength, limitPloidy, maxStates,
                             quick=FALSE, readPositionModel=FALSE){
  
  require(S4Vectors, quietly = TRUE, warn.conflicts = FALSE)
  require(IRanges, quietly = TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(robustbase, quietly = TRUE, warn.conflicts = FALSE)
  lmrob_control = robustbase::lmrob.control(setting="KS2014",
                                            max.it = 10000, maxit.scale = 10000,
                                            scale.tol=1e-5, solve.tol = 1e-5,
                                            subsampling = "nonsingular")
  
  
  stopifnot(dim(segmCN)[[2]] == 1)
  valid_bins=binsToUseInternal(segmCN)
  chr_value = unlist(lapply(strsplit(rownames(segmCN), ":"), `[`, 1))
  stopifnot("xi" %in% names(fiti)) # make sure fiti is unpacked, i.e. references a single cell fiti
  n_bins = dim(segmCN)[[1]]
  
  if(mean(segmCN[,cellname]@assayData$copynumber, na.rm=TRUE) < 0.5){
    base::warning(paste0("Cell has on average less than 0.5 reads per bin! Check quality (reads: ", segmCN[,cellname]@phenoData@data$used.reads, ")"))

    df = dplyr::tibble()
  } else {
    
    # select peaks
    means = fiti$means * (1/fiti$prescale)
    w = fiti$weights
    idx = base::sort(w, index.return=TRUE, decreasing=TRUE, method="quick")$ix
    select_idx = idx[w[idx] > component_cutoff]
    
    peaks = means[select_idx]
    weights = fiti$weights[select_idx]
    cov = fiti$covariance[select_idx] / (fiti$prescale^2)
    
    # use peak differences and distance from 0
    gaps = c(abs(diff(means[select_idx])), means[select_idx])
    
    gaps_weights_1 = c()
    for (i in 1:(length(fiti$weights[select_idx])-1)){
      if(length(fiti$weights[select_idx]) == 1){
        break
      }
      j = i+1
      gaps_weights_1[i] = fiti$weights[select_idx][i] * fiti$weights[select_idx][j]
    }
    gaps_weights = c(gaps_weights_1, (fiti$weights[select_idx])^2)
    stopifnot(length(gaps_weights) == length(gaps))

    # X needs to be correctly filtered for NaN values
    tmp = Biobase::assayDataElement(segmCN, "segmented")[,cellname,drop=FALSE]
    is_nan = is.na(as.numeric(tmp))
        
    # subset data based on ploidyRegion
    stopifnot(is.null(ploidyRegion) || all(startsWith(ploidyRegion, "chr")))
    if(!is.null(ploidyRegion)){
      valid = valid_bins & (!is_nan) & chr_value %in% base::sub("^chr", "", ploidyRegion)
    } else {
      valid = valid_bins & (!is_nan)
    }
    
    X = Biobase::assayDataElement(segmCN, "segmented")[valid,cellname,drop=FALSE]
    Yraw = Biobase::assayDataElement(segmCN, "copynumber")[valid,cellname,drop=FALSE]
    stopifnot(all(!is.na(X)))
    stopifnot(all(!is.na(Yraw)))
    
    # NOTE raw is used for mean-variance relationship
    raw = retrieveAssayData(segmCN, cellname, value="copynumber", valid=valid)
    for(chr in names(raw)){
      if(all(is.na(S4Vectors::runValue(raw[[chr]])))){
        raw[[chr]] = NULL
      }
    }
    
    # n_reads needs to be adapted in case of active ploidyRegion restriction
    n_reads = sum(unlist(lapply(raw, sum, na.rm=TRUE)))
    if("expected.variance" %in% colnames(Biobase::pData(segmCN)[cellname,])){
      expected.variance = Biobase::pData(segmCN)[cellname,][["expected.variance"]]
    }else{
      expected.variance = NA
    }
    
    ## Compute position information - only available if read distance file has been created
    if("bamfile" %in% colnames(Biobase::pData(segmCN)[cellname,]) & "bampath" %in% colnames(Biobase::pData(segmCN)[cellname,])){
      position_path = sub("\\/\\/", "\\/", file.path(Biobase::pData(segmCN)[cellname,"bampath"], (paste0(sub('\\.bam$', '', Biobase::pData(segmCN)[cellname,"bamfile"]), ".position.tsv"))))
      if(readPositionModel && file.exists(position_path) && !(file.info(position_path)$size == 0)){
        
        # extract diff read information
        high_quality_regions = createGR(segmCN[,cellname])
        high_quality_regions = high_quality_regions[valid]
        chromosome = decode(seqnames(high_quality_regions))
        
        read_position_table = readr::read_tsv(position_path, col_types = "cdd", col_names=TRUE) %>%
          dplyr::filter(!is.na(start))
        
        read_table = GRanges(read_position_table$chromosome,
                             IRanges(read_position_table$start, read_position_table$start+abs(read_position_table$length)-1),
                             length=abs(read_position_table$length),
                             dist=c(diff(read_position_table$start), Inf))
        
        hits = findOverlaps(high_quality_regions, read_table)
        hq = read_table[subjectHits(hits)]
        n_overlap = countOverlaps(hq, hq)
        
        filter_dists <- function(x){
          return(x[x < quantile(x, 0.99, na.rm=TRUE)])
        }
        
        validate_fitdistr <- function(x){
          x = x[!is.na(x)]
          if(length(x) == 0 | all(x == 0)){
            return(0)
          }else{
            return(x)
          }
        }
        
        temp = dplyr::tibble(element = queryHits(hits), dists = read_table[subjectHits(hits)]$dist,
                             chrom = chromosome[queryHits(hits)], n_overlap=n_overlap)
      }
    }
    # end compute position information

    mylist <- list()
    counter = 1
    
    for (j in 1:length(gaps)){
      delta = gaps[[j]]
      weightsquared = gaps_weights[[j]]
      # this is the scale

      for (delta_map in 1:9){
        ##
        
        if (weightsquared < weight_cutoff){
          # skip implausible fits
          next
        }
        
        int_scale = delta_map / delta
        
        # save the variables of transform
        Y = (X * int_scale)
        
        # NOTE rpc and ploidy is computed on rounded, and segmented values of Y_prime -> only use for single cell data
        Y_prime = round(Y, digits=0)

        ploidy = mean(Y_prime)
        ploidy.mod = as.numeric(names(sort(table(Y_prime),decreasing=TRUE)[1])[1])
        ploidy.continuous = mean(Y)


        ## ploidy criteria for solutions
        if( ((!is.null(minPloidy)) && (ploidy < minPloidy)) || ((!is.null(maxPloidy)) && (ploidy > maxPloidy)) ){
          next
        }

        ## different types of error ##
        
        # rounding error, i.e. how strong is the segmented data different from a integer solution
        diff = Y-Y_prime
        error_seg_l1 = sum(abs(diff))
        error_seg_l2 = sqrt(sum(diff * diff))
        error_seg_sd = sd(abs(diff), na.rm = TRUE)
        error_seg_median = median(abs(diff), na.rm = TRUE)
        
        # data error, i.e. how different is all the raw data from the final integer solution
        Y_raw = (Yraw * int_scale)
        diff = Y_raw - Y_prime
        error_all_l1 = sum(abs(diff))
        error_all_l2 = sqrt(sum(diff * diff))
        error_all_sd = sd(abs(diff), na.rm = TRUE)
        error_all_median = median(abs(diff), na.rm = TRUE)
        
        ## compute rpc (reads per copy)
        # NOTE rpc and ploidy is computed on rounded, and segmented values!
        rpc = n_reads / sum(Y_prime)
        stopifnot(length(Yraw) == length(Y_prime))
        
        ## some alternative rpc metrics
        cut = floor(0.99 * length(Yraw))
        rpc.99 = sum(sort(Yraw)[1:cut]) / sum(sort(Y_prime[1:cut]))
        
        cut = floor(0.95 * length(Yraw))
        rpc.95 = sum(sort(Yraw)[1:cut]) / sum(sort(Y_prime[1:cut]))
 
        rvec = Yraw / Y_prime
        rpc.var = var(rvec, na.rm=TRUE)
        rpc.median = median(rvec, na.rm=TRUE)
        rpc.robust = base::mean(rvec, trim=0.05, na.rm=TRUE)
        p95 <- stats::quantile(rvec, 0.95, na.rm=TRUE)
        rpc.p95 <- mean(rvec[which(rvec <= p95)], na.rm=TRUE)
        
        ## compute alpha and beta for given scaling ##
        gr = retrieveSegmentation(Y_prime)
        if(quick){
          # ddpl = list()
          ddpl$df = dplyr::tibble()
        }else{
          # ddpl = computeModel(raw, gr, scale=int_scale, roundInteger=FALSE, limitPloidy=limitPloidy, minLength=minLength, trimLength=trimLength, debug=FALSE, scaffold=FALSE)  
          ddpl = computeBasicModel(raw, gr, scale=int_scale, roundInteger=FALSE, limitPloidy=limitPloidy, minLength=minLength, trimLength=trimLength)
        }
        
        
        if("bamfile" %in% colnames(Biobase::pData(segmCN)[cellname,]) & "bampath" %in% colnames(Biobase::pData(segmCN)[cellname,])){
          position_path = sub("\\/\\/", "\\/", file.path(Biobase::pData(segmCN)[cellname,"bampath"], (paste0(sub('\\.bam$', '', Biobase::pData(segmCN)[cellname,"bamfile"]), ".position.tsv"))))
          if(readPositionModel && file.exists(position_path) && !(file.info(position_path)$size == 0)){

            temp$copy = Y_prime[queryHits(hits)]

            dens = temp %>% dplyr::filter(chrom %in% sub("chr", "", selectRegion), copy > 0, dists >= 0) %>% dplyr::group_by(chrom, element, copy) %>%
              dplyr::summarise(n=n(),
                               elem_overlap = mean(n_overlap, na.rm=TRUE), 
                               p_geom = MASS::fitdistr(validate_fitdistr(dists[1:(length(dists)-1)]), densfun = "geometric")[["estimate"]],
                               p_geom_censored = MASS::fitdistr(validate_fitdistr(filter_dists(dists[1:(length(dists)-1)])), densfun = "geometric")[["estimate"]],
                               .groups = 'keep') %>%
              dplyr::group_by(chrom, copy) %>%
              dplyr::summarise(n_elem = n(),
                               mean_overlap = mean(elem_overlap, na.rm=TRUE), median_overlap = median(elem_overlap, na.rm=TRUE),
                               mgp = mean(p_geom, na.rm=TRUE), mgpc = mean(p_geom_censored, na.rm=TRUE), medgp = median(p_geom, na.rm=TRUE), medgpc = median(p_geom_censored, na.rm=TRUE),
                               .groups = 'keep') %>% dplyr::ungroup()
            
            dty = dens$mean_overlap[dens$copy >= 1 & dens$n_elem >= 10 & dens$copy <= (maxStates-2)]
            dtz.mgeom = dens$mgp[dens$copy >= 1  & dens$n_elem >= 10 & dens$copy <= (maxStates-2)] 
            dtz.mgeom_censored = dens$mgpc[dens$copy >= 1  & dens$n_elem >= 10 & dens$copy <= (maxStates-2)] 
            dtx = dens$copy[dens$copy >= 1 & dens$n_elem >= 10 & dens$copy <= (maxStates-2)]
            n_elem = dens$n_elem[dens$copy >= 1 & dens$n_elem >= 10 & dens$copy <= (maxStates-2)]
            
            if(sum(n_elem) <= 10 || length(dtx) == 0 || length(dty) == 0 || all(is.na(dtx) | dtx == 0) || all(is.na(dty) | dty == 0)){
              #warning(paste0("Issue in readPosition - ", length(dens$copy[dens$copy > 0]), "/", length(dens$muf36[dens$copy > 0]), "/", sum(is.na(dens$copy[dens$copy > 0])), "/", sum(is.na(dens$muf36[dens$copy > 0]))))
              prediction.quality = 0
              prediction.diploid = NA
              prediction.value = NA
              prediction.value.pgeom = NA
              prediction.value.pgeom_censored = NA
            }else{
              
              if(length(unique(dens$copy)) == 1){
                # single copy number state
                prediction.diploid = median(rep(dty, times=n_elem), na.rm=TRUE)
                prediction.value = median(rep(dty, times=n_elem), na.rm=TRUE)
                prediction.value.pgeom = median(rep(dtz.mgeom, times=n_elem), na.rm=TRUE)
                prediction.value.pgeom_censored = median(rep(dtz.mgeom_censored, times=n_elem), na.rm=TRUE)
                prediction.quality = NA
              }else{
                prediction.diploid = median(rep(dty, times=n_elem), na.rm=TRUE)
                
                ## predict robust regression
                prediction.output.overlap <- tryCatch(
                  {
                    model = withCallingHandlers( robustbase::lmrob(prediction ~ I(copy), weights=weight,
                                                                   data = data.frame(copy=dtx,
                                                                                     prediction=dty,
                                                                                     weight=n_elem),
                                                                   control=lmrob_control))
                    pred = predict(model, newdata=dplyr::tibble(copy=2.0))
                    prediction.quality = summary(model)[["adj.r.squared"]]
                    list("prediction"=pred[[1]], "prediction.quality"=prediction.quality)
                  },
                  error=function(cond) {
                    # Choose a return value in case of error
                    return(list("prediction"=NA, "prediction.quality"=NA))
                  }
                )
                prediction.value = prediction.output.overlap[["prediction"]]
                prediction.quality = prediction.output.overlap[["prediction.quality"]]
                
                ## predict robust regression
                prediction.value.pgeom <- tryCatch(
                  {
                    model.pgeom = withCallingHandlers( robustbase::lmrob(prediction.pgeom ~ copy + I(copy^2), weights=weight,
                                                                       data = data.frame(copy=dtx,
                                                                                         prediction.pgeom=dtz.mgeom,
                                                                                         weight=n_elem),
                                                                       control=lmrob_control))
                    prediction.value.pgeom= predict(model.pgeom, newdata=dplyr::tibble(copy=2.0))
                    prediction.value.pgeom[[1]]
                    # prediction.quality = summary(model)[["adj.r.squared"]]
                    # list("prediction"=pred[[1]], "prediction.quality"=prediction.quality)
                  },
                  error=function(cond) {
                    # Choose a return value in case of error
                    return(NA)
                  }
                )
                
                ## predict robust regression
                prediction.value.pgeom_censored <- tryCatch(
                  {
                    model.pgeom = withCallingHandlers( robustbase::lmrob(prediction.pgeom_censored ~ copy + I(copy^2) + 0, weights=weight,
                                                                       data = data.frame(copy=dtx,
                                                                                         prediction.pgeom_censored=dtz.mgeom_censored,
                                                                                         weight=n_elem),
                                                                       control=lmrob_control))
                    prediction.value.pgeom_censored= predict(model.pgeom, newdata=dplyr::tibble(copy=2.0))
                    prediction.value.pgeom_censored[[1]]
                    # prediction.quality = summary(model)[["adj.r.squared"]]
                    # list("prediction"=pred[[1]], "prediction.quality"=prediction.quality)
                  },
                  error=function(cond) {
                    # Choose a return value in case of error
                    return(NA)
                  }
                )
                
              }
              
            }
            
            }else{
              prediction.value=NA
              prediction.diploid=NA
              prediction.quality=NA
              prediction.value.pgeom = NA
              prediction.value.pgeom_censored = NA
            }
          }else{
            prediction.value=NA
            prediction.diploid=NA
            prediction.quality=NA
            prediction.value.pgeom = NA
            prediction.value.pgeom_censored = NA
          }
        
        if(is.null(ddpl)){
          ## computeModel returns NULL, if model cannot be fitted in cases where not enough segments are valid
          next
        }

        ## store results
        # NOTE: alpha, beta, scale and rpc are required in other parts and should be part of this output
        mylist[[counter]] = dplyr::bind_cols(dplyr::tibble(
          name=cellname, n_reads=n_reads, scale=int_scale, rpc=rpc,
          error=error_seg_l2,  ploidy=ploidy, n_bins = n_bins,
          prediction=prediction.value, prediction.diploid=prediction.diploid, prediction.quality=prediction.quality,
          prediction.pgeom = prediction.value.pgeom, prediction.pgeom_censored = prediction.value.pgeom_censored,
          scaling.rpc_var=rpc.var, scaling.rpc_median=rpc.median, scaling.rpc_robust=rpc.robust,
          scaling.rpc_p95=rpc.p95, scaling.rpc_95=rpc.95, scaling.rpc_99=rpc.99), ddpl$df,
          dplyr::tibble(fit_flag=TRUE,
          ploidy.continuous=ploidy.continuous, ploidy.mod=ploidy.mod,
          expected.variance=expected.variance,
          delta=delta, weight=weightsquared, delta_map=delta_map,
          error_seg_l1 = error_seg_l1, error_seg_l2 = error_seg_l2, error_seg_sd=error_seg_sd, error_seg_median=error_seg_median,
          error_all_l1 = error_all_l1, error_all_l2 = error_all_l2, error_all_sd=error_all_sd, error_all_median=error_all_median))
        
        counter = counter + 1
      }
      
    }
    
    df <- dplyr::bind_rows(mylist)

  }
  
  if(!(nrow(df) > 0)){

    ddpl_scaffold = computeModel(NULL, NULL, NULL, NULL, NULL, NULL, debug=FALSE, scaffold=TRUE)
    df = dplyr::bind_cols(dplyr::tibble(
      name=cellname, n_reads=segmCN[,cellname]@phenoData@data$used.reads, scale=1.0, rpc=0.0,
      error=Inf, ploidy=0.0, n_bins=n_bins,
      prediction=NA, prediction.36 =NA, prediction.125=NA, prediction.0=NA, prediction.diploid=NA, prediction.quality=NA,
      prediction.pgeom = NA, prediction.pgeom_censored = NA,
      scaling.rpc_var=NA, scaling.rpc_median=NA, scaling.rpc_robust=NA,
      scaling.rpc_p95=NA, scaling.rpc_95=NA, scaling.rpc_99=NA), ddpl_scaffold,
      dplyr::tibble(fit_flag=FALSE,
      ploidy.continuous=NA, ploidy.mod=NA,
      expected.variance=NA,
      delta=NA, weight=0.0, delta_map=NA,
      error_seg_l1 = NA, error_seg_l2 = NA, error_seg_sd=NA, error_seg_median=NA,
      error_all_l1 = NA, error_all_l2 = NA, error_all_sd=NA, error_all_median=NA))
    
    return(df)
  }

  return(df)
}

#' selectSolution
#'
#' \code{selectSolution} given a fit obtained by the computeScale function, we select the final scaling among the remaining discrete solutions. 
#' 
#' This function supports multiple methods for selecting the final solution. If there is no previous data available, you can use the 'error' method. 
#' Ideally, you have some prior information on the ploidy limits of your sample and can use these to further restrict the solutions obtained by the 'error' method.
#' In case that you know the approximate ploidy level of your data, you can use the error function in combination with a tight limit on the possible ploidy levels via minPloidy and maxPloidy.
#' 
#' Alternatively, one can use a model based approach. Make sure to choose the same bin size between the model and your data! 
#' Here, a model needs to be specified in the model argument. You can either specify the full path to the model files in the source directory or use the available pretrained models.
#' Alternatively, you can build your own model using a large set of data with known ploidy.
#' 
#' @param segCN a segmented QDNAseq object
#' @param fit fit object obtained from running the computeScale function
#' @param method Choose among 'error' and 'model', i.e a model-free approach (ideally combined with prior knowledge of the ploidy) and 'model' as a model-based approach. 
#' 
#' @param globalModel character or path. Can either specify the path to a custom model file stored as an RDS file
#' or can specify one of the pre-trained models. We provide pre-built models for a range of bin sizes, please see the vignette for details or the data folder on github.
#' A model needs to support the predict function and predict a value as a function of rpc. 
#' Specifying a model is required when the method chosen is model. Specifying a model with the error method will throw an error.
#' The predict functionality can be overloaded with a custom predictFunction to be specified as described below. 
#' @param predictFunction option to supply custom function to predict alpha value.
#' @param minPloidy Numeric Hard limit on possible solution. Only solutions with higher or equal ploidy will be considered. 
#' @param maxPloidy Numeric Hard limit on possible solution. Only solutions with lower or equal ploidy will be considered. 
#' @param ploidyWindow Numeric (default: 0.1) Given a model based fit, we allow all solutions that lie within a ploidy window of 0.1 of the optimal solution. 
#' In other words, we select (ploidy-)equivalent solutions based on this threshold. Among these, the solution with minimal error is selected.
#' @param ploidyRegion Character (default: NULL) Chromosomes to restrict analyis to. This refers only to the selection of the ploidy region, i.e. it would be possible to select ploidy of a given chromosome.
#' @param trimLength Numeric Start and end of segments will be removed, to reduce bias in the estimate of the beta (\eqn{\beta_{1}}).
#' @param minLength Numeric Minimum length of segments to be included in \eqn{\alpha} and \eqn{\beta} prediction. This is used to minimize variance of the beta estimate (\eqn{\beta_{1}}).
#' @param limitPloidy (segments to discard from fitting of linear model - high values tend to be outlier and unreliable)
#' 
#' @param component_cutoff Numeric Used to select discrete subset of solutions from the fit output. Value is component weight as inferred by DPGMM. The lower, the more discrete solutions will be considered. Default 0.01
#' @param weight_cutoff Numeric Used to select discrete subset of solutions from the fit output. Value is product of two component weights. The lower, the more discrete solutions will be considered. Default 0.01
#' @param debug Add additional information to output. Useful for debugging purposes or deep dive into the algorithm.
#'
#' @return A QDNAseq object rescaled with the estimate of the 'best' solution. Some additional information is added to the phenoData slot.
#'
#' @export
selectSolution <- function(segCN, fit, method, globalModel, predictFunction=NULL,
                           minPloidy=1.2, maxPloidy=10.0,
                           ploidyWindow=0.1, ploidyRegion=NULL,
                           trimLength=1, minLength=30, limitPloidy=16,
                           component_cutoff=0.01, weight_cutoff = 0.001, maxStates=8,
                           debug=FALSE, quick=FALSE, readPositionModel=FALSE){

  method = tolower(method)
  if(!(method %in% c("error", "model"))) stop("Please choose method as error or model based!")
  stopifnot((minLength - 2*trimLength) >= 1)

  if(method %in% c("model")){
    stopifnot(!is.null(globalModel))
  }else{
    stopifnot(is.null(globalModel))
  }

  # select appropriate ploidy level
  # e.g. list("cella" = 1.2,"cellb" = 1.2,"cellc" = 1.3)
  assignValue <- function(argument, celle){
    if(is.null(argument)){
      return(NULL)
    }

    if(typeof(argument) == "list"){
      return(argument[[celle]])
    } else {
      return(argument)
    }
  }

  cellnames = colnames(segCN)
  dd = lapply(1:ncol(segCN),
          function(ij) evaluateScalings(segCN[,ij], fit[[ij]], cellnames[[ij]],
                                       component_cutoff=component_cutoff, weight_cutoff=weight_cutoff,
                                       minPloidy=assignValue(minPloidy, cellnames[ij]), maxPloidy=assignValue(maxPloidy, cellnames[ij]),
                                       ploidyRegion=assignValue(ploidyRegion, cellnames[ij]),
                                       trimLength=trimLength, minLength=minLength, 
                                       limitPloidy=limitPloidy, maxStates=maxStates, quick=FALSE, readPositionModel=readPositionModel))
  df = dplyr::bind_rows(dd)

  df = df %>% dplyr::mutate(norm_error = error / scale, norm_error_ploidy = error / ploidy)

  ##################################################################################
  ### Crucial part - select final scaling solution                               ###
  # because of numeric instability, beta should not be included in the distinct here
  df = df %>% dplyr::distinct(name, error, ploidy, .keep_all = TRUE)
  
  if(dim(df)[[1]] == 1 && df$rpc == 0){
    transform = df
  }else{
    
    # BEGIN select transform
    if (method == "error") {
      # simply select minimum error solution
      transform = (df %>% dplyr::group_by(name) %>% dplyr::slice(which.min(error)))
    } else {
        
      # either load globalModel from file or use globalModel directly from R object
      if (is.character(globalModel)) {
        
        stopifnot(file.exists(globalModel) || file.exists(paste0("~/scAbsolute/data/models/", globalModel)) || 
                    file.exists(normalizePath(paste0(.libPaths(), "/scAbsolute/data/models/", globalModel))))
        if(file.exists(globalModel)){
          prediction_fit = readRDS(globalModel)  
        }else{
          if(file.exists(paste0("~/scAbsolute/data/models/", globalModel))){
            prediction_fit = readRDS(paste0("~/scAbsolute/data/models/", globalModel))  
          }else{
            prediction_fit = readRDS(normalizePath(paste0(.libPaths(), "/scAbsolute/data/models/", globalModel))) 
          }
        }
      }
        
        if(is.null(predictFunction)){
          df$predicted = predict(prediction_fit, newdata = dplyr::tibble(rpc = df[["rpc"]]))
        }else{
          df$predicted = predictFunction(prediction_fit, df)
        }
  
        df$name = as.character(df$name)
        # crucial computation of residual, i.e. difference between predicted and actual regression coefficients # 
        df$residual = df$prediction - df$predicted
  
        stopifnot(all(df$error >= 0))
        # create variables for different criteria to select best one
        df = df %>% dplyr::group_by(name) %>% dplyr::mutate(
          tmp_min_res = min(abs(residual)),
          tmp_min_error = min(error)
        ) %>%
          dplyr::mutate(
            sel_error = case_when(error == tmp_min_error ~ TRUE,
                                  TRUE ~ FALSE),
            sel_residual = case_when(abs(residual) == tmp_min_res ~ TRUE,
                                     TRUE ~ FALSE)
          ) %>% # below, we select min residual, then min error within the min residual ploidy equivalence class
          dplyr::group_by(name) %>% dplyr::mutate(tmp_resmin_ploidy = ploidy[which.min(abs(residual))]) %>%
          dplyr::group_by(name) %>% dplyr::mutate(mod_ploidy = abs(ploidy - tmp_resmin_ploidy)) %>%
          dplyr::group_by(name) %>% dplyr::mutate(sel_ploidy_equivalence = case_when(mod_ploidy <= ploidyWindow ~ TRUE,
                                                                                     TRUE ~ FALSE)) %>%
          dplyr::group_by(name, sel_ploidy_equivalence) %>% dplyr::mutate(equi_min_error = min(error)) %>%
          dplyr::ungroup() %>% dplyr::group_by(name) %>%
          dplyr::mutate(sel_equivalence = case_when(((sel_ploidy_equivalence == TRUE) &
                                                       (error == equi_min_error)) ~ TRUE,
          TRUE ~ FALSE)) %>%
          dplyr::select(-c(
            tmp_min_error,
            tmp_min_res,
            tmp_resmin_ploidy
          ))
        
        # debug only
        #   transform = df %>% dplyr::group_by(name) %>% dplyr::slice(which(sel_residual))
        # we don't select minimum residual (above code), but we choose an equivalence class of similar ploidy (based on size of ploidyWindow)
        # and we select minimum error within this equivalence class
        transform = df %>% dplyr::group_by(name) %>% dplyr::slice(which(sel_equivalence))
  
    }
  }
  # END select transform
  stopifnot(length(transform$name) == ncol(segCN))
  transform = transform[match(colnames(segCN), transform$name), ]

  ### END ###
  ##################################################################################
  if(debug){
    meta = base::data.frame(labelDescription = rep(NA, times=dim(df)[[2]]))
    rownames(meta) = colnames(df)
    ProtData <- new("AnnotatedDataFrame",
                       data=as.data.frame(df), 
                       varMetadata = meta)
    
    protocolData(segCN) = ProtData
  }

  stopifnot(all(transform$name == colnames(segCN)))
  scale = transform %>% dplyr::pull(scale)

  ## Join selection data
  sample_names = rownames(Biobase::pData(segCN))
  description <- Biobase::pData(segCN)
  transform$name = as.character(transform$name)
  transform = as.data.frame(transform)
  stopifnot(all(sample_names == transform$name))
  
  # NOTE remove duplicate data frame names
  transform_subset = transform[c("name", base::setdiff(names(transform), names(description)))]
  transform_ext = transform[c("name", base::intersect(names(transform), names(description)))]
  colnames(transform_ext) = paste0(colnames(transform_ext), ".seg")
  transform_ext["name"] = transform_ext["name.seg"]; transform_ext["name.seg"] = NULL;
  description = dplyr::left_join(description, transform_subset, by="name") %>% dplyr::left_join(transform_ext, by="name")
  
  rownames(description) = sample_names
  stopifnot(all(rownames(description) == transform$name))
  Biobase::pData(segCN) = description
  
  
  # Transform data with final fit
  scaledCN = applyScale(segCN, scale=Biobase::pData(segCN)[["scale"]])
  
  return(scaledCN)
}


#' scAbsolute
#'
#' \code{scAbsolute} wrapper function to run scAbsolute for absolute copy number calling in single cell DNA sequencing data
#'
#' @param input Character file paths to bam files OR alternatively QDNAseq object (output of readData) OR pre-segmented QDNAseq object (skip read/segmentation steps)
#' @param method Choose among 'error' and 'model', i.e a model-free approach (ideally combined with prior knowledge about the ploidy) and 'model' as a model-based approach. 
#' @param globalModel Character specify path to model (model needs to implement a predict function of rpc)
#' Specifically: df$pred_alpha = predict(prediction_fit, newdata = dplyr::tibble(rpc = df[["rpc"]]))
#' More generally, we also support a signature of the type: df$pred_alpha = predictFunction(prediction_fit, df)
#' @param binSize Numeric Genomic bin size in kbp (supported values from QDNAseq: 1, 5, 10, 15, 30, 50, 100, 500, 1000)
#' @param species Character We support "Human" and "Mouse" at the moment.
#' @param genome Choice only for Human between hg19 and hg38.
#' @param filterChromosomes List of chromosomes to exclude when reading the data. Default: c("MT")
#' 
#' @param minLength Numeric Minimum segment length, NOTE this should be chosen in conjunction with binSize
#' @param trimLength Numeric Number of bins at start and end of segment to remove when calculating sufficient statistics (default is 1)
#' @param limitPloidy Numeric Cutoff for segments to consider in estimation of alpha, segments with a scaled ploidy higher than limitPloidy are not included in the model building process
#' 
#' @param minPloidy Numeric Minimum ploidy of solutions (can be specified for each sample as in:
#' Example: cellnames should be file name of bam file without ".bam" ending
#' cellnames = c("UID-10X-COLO829_SLX-00000_001313_TGAAAGACACCGTCGA-1", "UID-10X-COLO829_SLX-00000_000022_AACCATGAGGCTATCT-1"))
#' minPloidy = list(1.7, 3.5); names(minPloidy) = cellnames
#' maxPloidy = list(2.3, 5.5); names(maxPloidy) = cellnames
#' @param maxPloidy Numeric Maximum ploidy of solution. See minPloidy for details.
#' @param ploidyWindow Numeric (default: 0.1) Given a model based fit, we allow all solutions that lie within a ploidy window of 0.1 of the optimal solution. 
#' In other words, we select (ploidy-)equivalent solutions based on this range Among these, the solution with minimal error is selected.
#' @param ploidyRegion Character vector Region to use for ploidy calculation (default is all chromosomes except Y) and alpha/beta calculation
#' @param selectRegion Character vector Region to select for processing (default is all chromosomes (1:22, X, Y))
#' 
#' @param batchSize Numeric Number of cells to process in one batch (relevant for tensorflow implementation of SVI) (default: 1)
#' @param n_init Numeric Number of initialisations of SVI algorithm (default: 3)
#' @param n_steps Numeric Number of epochs for each SVI initialisation before stopping (default: 12)
#' @param truncation Numeric Truncation value for copy number levels in SVI algorithm. Default: 128
#' @param randomSeed Numeric Random seed used with SVI algorithm (for reproducibility)
#' 
#' @param testStatistic Character statistical model to use to determine optimal segmentation
#' By default we are using a Negative Binomial model with Methods of Moments estimator, but it is also possible to specify a Normal or Poisson model here
#' Options are: "NegBinML", "NegBinMM", "NegBinBCMM", "NegBinCML", "Poisson", "Normal"
#' @param penalty penalty to use for segmentation, default is MBIC, alternative see changepoint package
#' @param pen.value penalty to use when penalty is set to "Manual", see changepoint package for details
#' @param splitPerChromosome Boolean flag to speed up segmentation by segmentation each chromosome separately (default: false)
#' @param initialTestStatistics testStatistic to use for first pass through data (by default same as testStatistic)
#' 
#' @param outputPath Character write segmentation information to this path, (default NULL, works only when outputSegmentation is TRUE)
#' @param outputSegmentation Boolean flag to confirm writing segmentation table to outputPath (default: False)
#' @param optimizeSegmentation Boolean flag to run optimization procedure instead of single point segmentation (default: False)
#' @param gcCorrection Boolean flag to run gc and mappability correction procedure (default: True)
#' @param doCellcycle Boolean flag to run cellcycle tests
#' @param debug Boolean flag, include debug information from scaling procedure in output object (default: False)
#' @param quick Boolean flag, skip meanvar analysis, useful when running error model and runtime is important (default: False)
#' @param readPositionModel Boolean run read position analysis part of pipeline
#'
#' @return A QDNAseq object scaled to absolute copy number state
#'
#' @export
scAbsolute <- function(input, method="error", globalModel=NULL,
                       binSize=binSize, species="Human", genome="hg19", filterChromosomes=c("MT"),
                       minLength=30, trimLength=1, limitPloidy=16,
                       minPloidy=1.2, maxPloidy=10.0, ploidyWindow=0.1, ploidyRegion=paste0("chr", as.character(seq(1:22))),
                       selectRegion=paste0("chr", c(as.character(seq(1:22)), "X", "Y")),
                       batchSize=1, n_init=3, n_steps=12, truncation = 128, randomSeed=2020,
                       testStatistic="NegBinMM", penalty="MBIC", pen.value=NULL, splitPerChromosome=FALSE,  initialTestStatistic=NULL,
                       max_iterations=201, hmm_path=NULL, change_prob=1e-1, max_states=10,
                       outputPath=NULL, outputSegmentation=FALSE, optimizeSegmentation=FALSE,
                       gcCorrection=TRUE, doCellcycle=TRUE, debug=FALSE, quick=FALSE, skipForEvaluation=FALSE, readPositionModel=TRUE){
  
  # Make sure either error or model approach is specified!
  if(method %in% c("model")){
    stopifnot(!is.null(globalModel))
  }else{
    stopifnot(is.null(globalModel))
  }
  
  program_start_time <- Sys.time()

  if(class(input)[[1]] == "QDNAseqCopyNumbers" && ("segmented" %in% Biobase::assayDataElementNames(input))){
    ## Skip reading and segmenting when input is segmented QDNAseq object
    warning("Input is segmented - skipping default steps")
    segmentedCounts = input
    
  }else{
    ## 1. read in data or use data
    start_time <- Sys.time()
    if(class(input)[[1]] == "QDNAseqCopyNumbers"){
      readCounts = input
    } else {
      readCounts = readData(input, binSize=binSize,
                            species=species, genome=genome,
                            filterChromosomes = filterChromosomes)
    }
    end_time <- Sys.time()
    stopifnot(is.null(selectRegion) || all(ploidyRegion %in% selectRegion))
    if(base::class(minPloidy) == "list"){
      stopifnot(all(Biobase::pData(readCounts)[["name"]] %in% names(minPloidy)))
    }
    if(base::class(maxPloidy) == "list"){
      stopifnot(all(Biobase::pData(readCounts)[["name"]] %in% names(maxPloidy)))
    }
    if(!is.null(selectRegion)){
      readCounts = selectChromosomes(readCounts, include=selectRegion)
    }
    print(paste0("AA readData runtime ",  difftime(end_time,start_time,units="mins")))
  
    # initial segmentation to provide segments for gc correction
    if(gcCorrection){
      start_time <- Sys.time()
      # initial segmentation is best done with AIC - better to slightly oversegment at this stage
      if(is.null(initialTestStatistic)){
        initialTestStatistic = testStatistic
      }
      initialSegmentedCounts = segment(readCounts, penalty = "MBIC", pen.value = NULL, testStatistic = initialTestStatistic,
                                       splitPerChromosome=TRUE)
      initialSegmentedCountsCorrected = initial_gc_correction(initialSegmentedCounts)
      end_time <- Sys.time()
      print(paste0("AB initialSegment runtime ",  difftime(end_time,start_time,units="mins")))
    }else{
      initialSegmentedCountsCorrected = readCounts
    }
  
    ## 2. Segmentation
    if(!optimizeSegmentation){
      start_time <- Sys.time()
      segmentedCounts = segment(initialSegmentedCountsCorrected, penalty = penalty, pen.value = pen.value, testStatistic = testStatistic,
                                       splitPerChromosome=splitPerChromosome, debug=debug)
      end_time <- Sys.time()
      print(paste0("B segment runtime ",  difftime(end_time,start_time,units="mins")))
     
    }else{
      start_time <- Sys.time()
      if(splitPerChromosome){
        stop("splitPerChromosome TRUE is not supported with segment_CROPS")
      }
      segmentedCounts = segment_CROPS(initialSegmentedCountsCorrected, pen.value = pen.value, testStatistic = testStatistic)
      end_time <- Sys.time()
      print(paste0("B segmentCROPS runtime ",  difftime(end_time,start_time,units="mins")))
    }
  }
  
  if(!is.null(outputPath) & outputSegmentation){
    ifelse(!dir.exists(dirname(outputPath)), dir.create(dirname(outputPath), recursive = TRUE), FALSE)
    outputName = paste0(dirname(outputPath), "/", sub('[.][^.]+$', '', basename(outputPath)))
    writeSegTable(segmentedCounts, filename=paste0(outputName, paste0(".", testStatistic, "-raw.bed")), trimLength=trimLength)
  }
  
  ## 4. cellcycle metadata (S phase) 
  # NOTE: classify before scaling
  if(doCellcycle){
    start_time <- Sys.time()
    segmentedCounts = cellcycleMetadata(segmentedCounts)
    segmentedCounts = cellcycle_gctest(segmentedCounts)
    segmentedCounts = computeHigherOrderMoments(segmentedCounts, minLength = minLength)
    end_time <- Sys.time()
    print(paste0("D cellcycleMetadata runtime ",  difftime(end_time,start_time,units="mins")))
  }

  
  ## 5. compute best fit using absolute graphical model and stochastic variational inference
  # (NOTE that this computation is done in python and requires reticulate and other python dependencies)
  # we use randomSeed to set the python random seed
  start_time <- Sys.time()
  fit = computeScale(segmentedCounts, truncation = truncation, batchSize=batchSize, 
                     n_init = n_init, n_steps = n_steps,
                     randomSeed=randomSeed, verbosity=2)
  end_time <- Sys.time()
  print(paste0("E computeScale runtime ",  difftime(end_time,start_time,units="mins")))

  ## 6. select optimal copy number fit based either on error or on existing model
  # model approach is based on a global model of the mean-variance relationship in these kind of count data
  start_time <- Sys.time()
  scaledCN = selectSolution(segmentedCounts, fit, method=method, globalModel=globalModel, debug=debug, limitPloidy=limitPloidy,
                            maxStates=max_states,
                            ploidyWindow=ploidyWindow, minPloidy=minPloidy, maxPloidy=maxPloidy, ploidyRegion=ploidyRegion,
                            trimLength = trimLength, minLength = minLength, quick=quick, readPositionModel=readPositionModel)
  end_time <- Sys.time()
  if(any(Biobase::pData(scaledCN)[["rpc"]] <= 0.0)){
         warning(paste0("Sample failed - ploidy constraint unsatisfiable\n", Biobase::pData(scaledCN)[["name"]]))
  }
  stopifnot(all(Biobase::pData(scaledCN)[["rpc"]] > 0.0))
  print(paste0("F selectSolution runtime ",  difftime(end_time,start_time,units="mins")))
  
  ## 7. finalize segmentation by using HMM
  protData = Biobase::protocolData(scaledCN)
  start_time <- Sys.time()
  scaledsegmentedCN = combineQDNASets(lapply(1:ncol(scaledCN), function(li){cs=copynumberSegmentation(scaledCN[,li], change_prob=change_prob,
                                                                       max_iterations=max_iterations, max_states=max_states,
                                                                       hmm_path=hmm_path, verbose=debug, gc_correction=gcCorrection, splitPerChromosome = splitPerChromosome); return(cs[["object"]])}))
  end_time <- Sys.time()
  Biobase::protocolData(scaledsegmentedCN) = protData
  print(paste0("G finalizeSegmentation runtime ",  difftime(end_time,start_time,units="mins")))
  
  program_end_time <- Sys.time()

  if(!skipForEvaluation){
      countdata = selectChromosomes(scaledsegmentedCN, exclude=c("X", "Y"))
      valid = binsToUseInternal(countdata)
      countdat = Biobase::assayDataElement(countdata, "calls")
      fitdistrs = suppressWarnings(lapply(1:ncol(scaledsegmentedCN), function(li) return(MASS::fitdistr(countdat[valid, li], 'negative binomial'))))
      description <- Biobase::pData(scaledsegmentedCN)
      description$runtime = difftime(program_end_time,program_start_time,units="mins")
      description$fitdistr.mu = unlist(lapply(fitdistrs, function(x) return(x[["estimate"]][["mu"]])))
      description$fitdistr.alpha = unlist(lapply(fitdistrs, function(x) return(1/x[["estimate"]][["size"]])))
      Biobase::pData(scaledsegmentedCN) = description
  }
  
  if(!is.null(outputPath) & outputSegmentation){
    ifelse(!dir.exists(dirname(outputPath)), dir.create(dirname(outputPath), recursive = TRUE), FALSE)
    outputName = paste0(dirname(outputPath), "/", sub('[.][^.]+$', '', basename(outputPath)))
    writeSegTable(scaledsegmentedCN, filename=paste0(outputName, paste0(".", testStatistic, "-scaled.bed")), trimLength=trimLength, minLength=minLength)
  }
  
  return(scaledsegmentedCN)
}
