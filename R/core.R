## core methods for copynumber analysis
# Copyright 2022, Michael Schneider, All rights reserved.

#' binsToUseInternal
#'
#' \code{binsToUseInternal} modified method from QDNAseq package to select bins to use
binsToUseInternal <- function(object){
  require(Biobase, warn.conflicts = FALSE)
  if ("use" %in% colnames(Biobase::fData(object))) {
    return(Biobase::fData(object)$use)
  }else{
    warning("No filter available for this object.")
  }

  return(rep(TRUE, times=nrow(object)))
}

#' getGenomeInformation
#'
#' \code{getGenomeInformation} access QDNAseq object genome information
getGenomeInformation <- function(object){
  return(Biobase::experimentData(object)@other)
}

#' createGR 
#'
#' \code{createGR} create a GRanges representation 
#'
#' @param object QDNAseq object
#' @param value character add QNDAseq assayData element of object to the GRanges object as a metadata column [copynumber, calls, segmented]
#' @param valid index include only elements at valid positions
#' @param asRleList convert GRanges output to RleList
#' 
#' @return GRanges object (or RleList, if asRleList is set)
createGR <- function(object, value=NULL, valid=NULL, asRleList=FALSE, variable=NULL){
  require(GenomicRanges, warn.conflicts = FALSE)
  stopifnot(dim(object)[[2]] == 1)

  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)

  tmp = strsplit(unlist(lapply(strsplit(info, ":"), `[`, 2)), "-")
  if(asRleList){
    chunkStart = do.call("c", lapply(rle(seqName)$lengths, function(x) 1:x))
    chunkEnd = do.call("c", lapply(rle(seqName)$lengths, function(x) 1:x))
  }else{
    chunkStart = as.numeric(unlist(lapply(tmp, `[`, 1)))
    chunkEnd = as.numeric(unlist(lapply(tmp, `[`, 2)))
  }
  if(!is.null(variable)){
    gr = GRanges(chromosome, IRanges(chunkStart, chunkEnd), annotation=variable)
  }else{
    gr = GRanges(chromosome, IRanges(chunkStart, chunkEnd))
  }
  

  if(!is.null(value)){
    val = Biobase::assayDataElement(object, value)
    if(!is.null(valid)){
      val[!valid] = NA
    }

    df = data.frame(val=as.numeric(val))
    colnames(df) = value
    mcols(gr) = df

    if(asRleList){
      gr = GenomicRanges::mcolAsRleList(gr, value)
    }
  }else{
    if(asRleList){
      stop("asRleList not supported if value is not set")
    }
    if(!is.null(valid)){
      gr = gr[valid]
    }
  }

  return(gr)
}

#' retrieveAssayData 
#'
#' \code{retrieveAssayData} extract assayData from QDNAseq object
#'
#' @param object QDNAseq object
#' @param index sample index into QDNAseq object
#' @param value character QNDAseq assayData element of object to the GRanges object as a metadata column [copynumber, calls, segmented]
#' @param valid index include only elements at valid positions
#' @param asRleList convert GRanges output to RleList
#' 
#' @return RleList with value as metadata
#' 
#' @export
retrieveAssayData <- function(object, index, value, valid=NULL){
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  stopifnot(length(index) == 1)
  return(createGR(object[,index,drop=FALSE], value=value, valid=valid, asRleList=TRUE))
}

#' retrieveSegmentation 
#'
#' \code{retrieveSegmentation} create a GRanges representation 
#'
#' @param object QDNAseq object
#' @param index sample index into QDNAseq object
#' @param valid index include only elements at valid positions
#' @param formatIR return IRanges instead of GRanges
#' 
#' @return return GRanges object (or IRanges if flag is set) containing segmentation
#' 
#' NOTE
#' segment encoding is processed as follows:
#   1. Remove NaNs that are covered by two equal segments (30 NA 30) -> 30
#   2. Keep NaNs that border two unequal segments (30 NA 31) -> 30 NA 31
#   3. Remove NaNs at start or end (NA 30) -> 30, (30 NA) -> 30
#   4. Do this on a per chromosome level (chromosome encoding in chromosome variable)
#   5. Keep NaN elements in segmented
retrieveSegmentation <- function(object, index, valid=NULL, formatIR=FALSE, slot="segmented"){

  if(class(object)[[1]] != "QDNAseqCopyNumbers"){
    segments = object
    stopifnot(!is.null(rownames(object)))
    stopifnot(is.null(valid))
  }else{
    stopifnot(!is.null(Biobase::assayDataElement(object, slot)))
    if(is.null(valid)){
      segments = Biobase::assayDataElement(object, slot)[,index,drop=FALSE]
    }else{
      segments = Biobase::assayDataElement(object, slot)[valid,index,drop=FALSE]
    }
  }
  
  replacena <- function(rl) {
    lv = rl@values
    ll = rl@lengths
    stopifnot(is.numeric(lv))
    indx <- is.na(lv)
    mat = matrix(c(lv, ll), nrow=2, byrow = TRUE)
    
    if(all(is.na(lv))){
      return(list(lv=lv, ll=ll))
    }
    
    for(x in 1:length(lv)){
      if(is.na(lv[x])){
        
        # left border
        if(x == 1){
          mat[1,x+1] = mat[1,2]
          mat[2,x+1] = mat[2,1] + mat[2,2]
          # mat[1,x] = NA
          mat[2,x] = 0
          next
        }
        
        # right border
        if(x == length(lv)){
          mat[1,x-1] = mat[1,x-1]
          mat[2,x-1] = mat[2,x-1] + mat[2,x]
          # mat[1,x] = NA
          mat[2,x] = 0
          next
        }
        
        ## normal case
        left = x-1
        right = x+1
        
        if(mat[1,left] == mat[1,right]){
          mat[1,x] = mat[1,left]
          mat[2,x] = mat[2,left] + mat[2,x] + mat[2,right]
          # mat[1,left] = NA
          mat[2,left] = 0
          # mat[1,right] = NA
          mat[2,right] = 0
        }
        
      }
    }
    
    lv = mat[1,]
    ll = mat[2,]
    return(list(lv=lv, ll=ll))
  }
  
  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  
  segments.cleaned = do.call("RleList", lapply(split(segments, chromosome), function(x){
    tmp = replacena(Rle(x))
    lv = tmp$lv
    ll = tmp$ll
    remove_idx = ll == 0
    lv = lv[!remove_idx]
    ll = ll[!remove_idx]
    vec = as.vector(Rle(values=lv, lengths = ll))
    stopifnot(length(vec) == length(x))
    ## NOTE: dummy encode NA as -1, to avoid filtering segments in bindAsGRanges -> replace values again after transform
    rl = Rle(as.numeric(vec))
    rl@values[is.na(rl@values)] = -1
    return(rl)
  }))
  
  # transform back to GRanges representation
  grs_new = bindAsGRanges(segmentation = segments.cleaned)
  grs_new$segmentation[grs_new$segmentation == -1] = NA
  # see note in function documentation
  
  stopifnot(sum(grs_new@ranges@width) == dim(object)[[1]])
  
  if(formatIR){
    n_bins = dim(object)[[1]]
    ir = grs_new@ranges
    ir_start = c(0, cumsum(ir@width))+1
    ir_end = cumsum(ir@width)
    ir = IRanges(start=ir_start[1:(length(ir_start)-1)], end=ir_end)
    return(ir)    
  }
  
  return(grs_new)
}

#' selectChromosomes 
#'
#' \code{selectChromosomes} select subregion of QDNAseq based on chromosome names either to be included or excluded
#'
#' @param object QDNAseq object
#' @param include character chromosome names (either in format "1" or "chr1") to include
#' @param exclude character chromosome names (either in format "1" or "chr1") to exclude
#' 
#' @return object containing only selected chromosomes
selectChromosomes <- function(object, include=NULL, exclude=NULL){
  require(S4Vectors, warn.conflicts = FALSE)
  
  if(!is.null(include) && !is.null(exclude)){
    stop("Choose either include or exclude regions.")
  }
  stopifnot(is.null(include) || all(is.character(include)))
  stopifnot(is.null(exclude) || all(is.character(exclude)))
  
  
  ploidyRegion=c(as.character(seq(1,22)), "X", "Y")
  if(!is.null(include)){
    if(all(base::startsWith(include, "chr"))){
      include = base::sub("^chr", "", include)
    }
    ploidyRegion = include
  }
  
  
  if(!is.null(exclude)){
    if(all(base::startsWith(exclude, "chr"))){
      exclude = base::sub("^chr", "", exclude)
    }
    ploidyRegion = ploidyRegion[!(ploidyRegion %in% exclude)]
  }
  
  gr = createGR(object[,1,drop=FALSE])
  indices = S4Vectors::decode(seqnames(gr) %in% ploidyRegion)
  
  return(object[indices, ,drop=FALSE])
}

#' smoothOperator
#' 
#' \code{smoothOperator} apply the summary function func to the data z. The summary function is computed 
#' for each segmented defined in the GRanges or IRanges object gr. Trim length reads are removed from the start and end of each segment
#'
#' @param z RleList or numeric vector (this will usually be the raw read count vector)
#' @param gr GRanges object with segmentation (z is RleList) or alternatively IRanges object (z is vector)
#' @param func function returning summary statistics for raw data (for a given segment)
#' @param minLength Numeric minimum number of bins to form a segment
#' @param trimLength Numeric number of bins to trim from start and end of segments
#' 
#' @return fit object containing best fit among the iterations from the SVI algorithm
smoothOperator <- function(z, gr, func, minLength=30, trimLength=0){
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  require(IRanges, quietly=TRUE, warn.conflicts = FALSE)
  
  if(class(gr)[[1]] == "IRanges"){
    index = gr@width >= minLength
  }else{
    if(class(gr)[[1]] == "GRanges"){
      index = gr@ranges@width >= minLength
    }else{
      stop("RangeObject not supported.")
    }
  }
  
  vw = viewApply(Views(z, gr[index]), function(x) {
    
    # necessary to avoid XDouble errors
    x = as.numeric(x)
    
    y = x[(1+trimLength):(length(x)-trimLength)]
    y = func(y)
    
    return(y)
  }, simplify = TRUE)
  
  if(class(gr)[[1]] == "GRanges"){
    vw = do.call("c", lapply(names(vw), function(x) return(unlist(vw[[x]]))))
  }
  
  return(vw)
}

#' applyScale
#'
#' \code{applyScale} helper function to apply a given scaling factor to the counts in the form of a linear transformation
#'
#' @param CN QDNAseq object
#' @param scale Numeric Linear scaling factor, i.e. a in y = a x + b
#' @param shift Numeric Offest factor in linear transformation, i.e. b in y = a x + b
#'
#' @return scaled copy number object
applyScale <- function(CN, scale=1, shift=0){
  
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  result = new("QDNAseqCopyNumbers",
               bins=Biobase::featureData(CN),
               copynumber=Biobase::assayDataElement(CN, "copynumber"), #t((diag(as.numeric(unname(scale)), nrow=dim(CN)[[2]]) %*% t(Biobase::assayDataElement(CN, "copynumber"))) + shift),
               phenodata=Biobase::pData(CN))
  Biobase::assayDataElement(result, "segmented") <- t(diag(as.numeric(unname(scale)), nrow=dim(CN)[[2]]) %*% t(Biobase::assayDataElement(CN, "segmented"))) + shift
  if(!is.null(Biobase::assayDataElement(CN, "calls"))){
    Biobase::assayDataElement(result, "calls") <- Biobase::assayDataElement(CN, "calls")
  }
  if(!is.null(Biobase::assayDataElement(CN, "probdloss"))){
    Biobase::assayDataElement(result, "probdloss") <- Biobase::assayDataElement(CN, "probdloss")
  }
  Biobase::protocolData(result) <- Biobase::protocolData(CN)
  return(result)
}

#' computeRPC
#'
#' \code{computeRPC} compute average reads per copy per bin (rpc) value for a given QDNAseq object (based on calls and segments).
#' segments need to be a proxy for copy number state, i.e. a scaled data
#' 
#' @param object a QDNAseq object
#' @param rpc_value 
#' @param ploidy_value 
#' 
#' @return average reads per copy per bin
computeRPC <- function(object, rpc_value="rpc", ploidy_value="ploidy"){
  
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  
  valid = binsToUseInternal(object)
  rpc_raw = Biobase::assayDataElement(object, "calls")[valid,,drop=FALSE] / Biobase::assayDataElement(object, "segmented")[valid,,drop=FALSE]
  rpc_raw[is.infinite(rpc_raw)] <- NA
  rpc = apply(rpc_raw, 2, sum, na.rm=TRUE) / sum(valid, na.rm=TRUE)
  ploidy = mean(Biobase::assayDataElement(object, "segmented")[valid,,drop=FALSE], na.rm=TRUE)
  
  valid = binsToUseInternal(object) & !(base::startsWith(rownames(object), "X:") | base::startsWith(rownames(object), "Y:"))
  rpc_raw = Biobase::assayDataElement(object, "calls")[valid,,drop=FALSE] / Biobase::assayDataElement(object, "copynumber")[valid,,drop=FALSE]
  rpc_raw[is.infinite(rpc_raw)] <- NA
  rpc_final = apply(rpc_raw, 2, sum, na.rm=TRUE) / apply(rpc_raw, 2, function(x) sum(!is.na(x)))
  
  pd = Biobase::pData(object)
  pd[[rpc_value]] = rpc
  pd[[ploidy_value]] = ploidy
  if(!is.null(pd[["rpc"]])){
    pd[["rpc.old"]] = pd[["rpc"]]  
  }
  pd[["rpc"]] = rpc_final
  Biobase::pData(object) = pd
  
  return(object)
}

#' computeGini
#'
#' \code{computeGini} compute Gini coefficient for a given QDNAseq object (based on calls)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return Gini coefficient
computeGini <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  Y = Biobase::assayDataElement(scaledCN, "calls")
  Y = Y[which(!is.na(Y[, 1])), ,drop=FALSE]
  
  igini <- function(x, weights=rep(1,length=length(x))){
    # source: https://rdrr.io/cran/reldist/src/R/gini.r
    ox <- order(x)
    x <- x[ox]
    weights <- weights[ox]/sum(weights)
    p <- cumsum(weights)
    nu <- cumsum(weights*x)
    n <- length(nu)
    nu <- nu / nu[n]
    sum(nu[-1]*p[-n]) - sum(nu[-n]*p[-1])
  }
  
  ngini <- function(i){
    obj = scaledCN[, i]
    stopifnot(dim(obj)[[2]] == 1)
    valid = binsToUseInternal(obj)
    Y = obj@assayData$calls[valid, ]
    C = obj@assayData$copynumber[valid, ]

    info = rownames(obj)
    seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)

    dft = dplyr::tibble(copy=C, count=Y, chromosome = chromosome[valid])
    dftp = dft %>% dplyr::group_by(chromosome, copy) %>% dplyr::summarise(n = n(), normalized_gini = igini(count), .groups="keep") %>%
      dplyr::ungroup() %>% dplyr::filter(copy > 0)
    
    va = dftp$n > 10
    normalized_gini = mean(rep(dftp$normalized_gini[va], times=dftp$n[va]), na.rm=TRUE)
    return(normalized_gini)
  }
  
  gini = apply(Y, 2, igini)
  gini_normalized = unlist(lapply(1:ncol(scaledCN), ngini))
  # scaledCN@phenoData@data$gini = gini
  pd = Biobase::pData(scaledCN)
  pd$gini = gini
  pd$gini_normalized = gini_normalized
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

#' computeMAPD
#'
#' \code{computeMAPD} compute MAPD coefficient for a given QDNAseq object (based on calls)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return MAPD coefficient
computeMAPD <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  Y = Biobase::assayDataElement(scaledCN, "calls")
  Y = Y[which(!is.na(Y[, 1])), ,drop=FALSE]
  
  mapd <- function(x){
    dd = (-1 * diff(x)) / mean(x)
    return(median( abs(dd - median(dd)) ))
  }
  
  mapd = apply(Y, 2, mapd)
  pd = Biobase::pData(scaledCN)
  pd$mapd = mapd
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

#' computeL2
#'
#' \code{computeL2} compute l2 error for a given QDNAseq object (based on calls and scaling)
#' 
#' @param scaledCN a QDNAseq object
#' 
#' @return l2 error (normalized)
computeL2 <- function(scaledCN){
  
  stopifnot(class(scaledCN)[[1]] == "QDNAseqCopyNumbers")
  stopifnot("rpc" %in% colnames(Biobase::pData(scaledCN)))
  valid = binsToUseInternal(scaledCN)
  cnv = Biobase::assayDataElement(scaledCN, "copynumber")[valid, ,drop=FALSE]
  rpcs = Biobase::pData(scaledCN)[["rpc"]]
  error = (t(t(Biobase::assayDataElement(scaledCN, "calls")[valid, ,drop=FALSE]) / rpcs) - cnv)
  l2 = apply(error, 2, function(x) sqrt(sum(x^2)))
  
  pd = Biobase::pData(scaledCN)
  pd$l2 = l2 
  Biobase::pData(scaledCN) = pd
  return(scaledCN)
}

# Compute information theoretic measures for scaled and segmented QDNAseq object
computeInfotheo <- function(scaledCN){
  
  require(infotheo)
  
  computeInfotheoSingle <- function(i){
    object = scaledCN[, i]
    valid = binsToUseInternal(object)
    rT = featureData(object)[["replicationTiming"]][valid]
    gc = featureData(object)[["gc"]][valid]
    X = Biobase::assayDataElement(object, "copynumber")[valid]
    Y = Biobase::assayDataElement(object, "calls")[valid]
    covariate_info = estimate_gc_correction(object)[valid]
    
    # replication time
    cellcycle.cmi_xy = condinformation(X, Y, S=infotheo::discretize(rT))
    cellcycle.cmi_yrT = condinformation(Y, infotheo::discretize(rT), S=X)
    
    # cell quality
    mi.cmi_yGC = condinformation(Y, infotheo::discretize(gc), S=X)
    mi.yGC = mutinformation(Y, infotheo::discretize(gc))
    
    mi.cmi_xy = condinformation(X, Y, infotheo::discretize(covariate_info))
    mi.xy = mutinformation(X, Y)
    mi.xm = mutinformation(X, infotheo::discretize(covariate_info))
    
    return(c(cellcycle.cmi_xy, cellcycle.cmi_yrT, mi.cmi_yGC, mi.yGC, mi.cmi_xy, mi.xy, mi.xm))
  }
  
  ve = do.call("rbind", lapply(1:ncol(scaledCN), computeInfotheoSingle))
  colnames(ve) = c("cellcycle.cmi_xy", "cellcycle.cmi_yrT", "mi.cmi_yGC", "mi.yGC", "mi.cmi_xy", "mi.xy", "mi.xm")
  
  pd = Biobase::pData(scaledCN)
  ve = as.data.frame(ve)
  rownames(ve) = rownames(pd)
  pd2 = dplyr::bind_cols(pd, ve)
  Biobase::pData(scaledCN) = pd2
  return(scaledCN)
}

# Skewness and kurtosis and their standard errors as implement by SPSS
#
# Reference: pp 451-452 of
# http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
# 
# See also: Suggestion for Using Powerful and Informative Tests of Normality,
# Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
# The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321
spssSkew=function(x){
  x = x[!is.na(x)]
  if(length(x) < 10){
    return( NA)
  }
  w=length(x)
  m1=mean(x)
  # m2=sum((x-m1)^2)
  m3=sum((x-m1)^3)
  # m4=sum((x-m1)^4)
  s1=sd(x)
  skew=w*m3/(w-1)/(w-2)/s1^3
  # sdskew=sqrt( 6*w*(w-1) / ((w-2)*(w+1)*(w+3)) )
  # kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
  # sdkurtosis=sqrt( 4*(w^2-1) * sdskew^2 / ((w-3)*(w+5)) )
  # se only valid for normality assumption
  # mat=matrix(c(skew,kurtosis, sdskew,sdkurtosis), 2,
  #            dimnames=list(c("skew","kurtosis"), c("estimate","se")))
  # return(mat)
  return(skew)
}

# Skewness and kurtosis and their standard errors as implement by SPSS
#
# Reference: pp 451-452 of
# http://support.spss.com/ProductsExt/SPSS/Documentation/Manuals/16.0/SPSS 16.0 Algorithms.pdf
# 
# See also: Suggestion for Using Powerful and Informative Tests of Normality,
# Ralph B. D'Agostino, Albert Belanger, Ralph B. D'Agostino, Jr.,
# The American Statistician, Vol. 44, No. 4 (Nov., 1990), pp. 316-321
spssKurtosis=function(x){
  x = x[!is.na(x)]
  if(length(x) < 10){
    return(NA)
  }
  w=length(x)
  m1=mean(x)
  m2=sum((x-m1)^2)
  m4=sum((x-m1)^4)
  s1=sd(x)
  kurtosis=(w*(w+1)*m4 - 3*m2^2*(w-1)) / ((w-1)*(w-2)*(w-3)*s1^4)
  return(kurtosis)
}


#' computeHigherOrderMoments
#'
#' \code{computeHigherOrderMoments} compute median skewness across segments for a given QDNAseq object (based on calls)
#' 
#' @param object a (segmented) QDNAseq object
#' 
#' @return higher order moments estimate
computeHigherOrderMoments <- function(object, minLength=20, trimLength=1){
  
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  stopifnot(all(c("calls", "segmented") %in% Biobase::assayDataElementNames(object)))
  require(future.apply, warn.conflicts = FALSE)
  require(matrixStats, warn.conflicts = FALSE)

  processSkewness <- function(obj){
    stopifnot(dim(obj)[[2]] == 1)
    
    valid_bins=binsToUseInternal(obj)
    counts = retrieveAssayData(obj, 1, value="calls", valid=valid_bins)
    gr = retrieveSegmentation(obj, 1)
  
    l_skewness = smoothOperator(counts, gr, function(x){spssSkew(x)}, trimLength = trimLength, minLength=minLength)
    l_kurtosis = smoothOperator(counts, gr, function(x){spssKurtosis(x)}, trimLength = trimLength, minLength=minLength)
    l_weight = smoothOperator(counts, gr, function(x){x = x[!is.na(x)]; return(length(x))}, trimLength = trimLength, minLength=minLength)
    
    return(list(skewness=matrixStats::weightedMedian(l_skewness, w=l_weight, na.rm=TRUE), kurtosis=matrixStats::weightedMedian(l_kurtosis, w=l_weight, na.rm=TRUE)))
  }
  
  results = future.apply::future_lapply(lapply(colnames(object), function(x) object[, x]), processSkewness)
  skewness = sapply(results, "[[", 1)
  kurtosis = sapply(results, "[[", 2)
  Biobase::pData(object)[["skewness"]] = skewness
  Biobase::pData(object)[["kurtosis"]] = kurtosis
  return(object)
}

#' addObservedVariance
#'
#' \code{addObservedVariance} compute observed variance (as defined by QDNAseq and displayed on the QDNAseq plots as standard deviation). 
#' See sdDiffTrim computation from QDNAseq.
#' 
#' @param object QDNAseq object
#' @param trim 
#' @param scale
#' 
#' @return object with additional slot for "observed.variance"
addObservedVariance <- function(object, trim=0.001, scale=TRUE) {
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  require(Biobase, warn.conflicts = FALSE)
  require(matrixStats, warn.conflicts = FALSE)

  x = Biobase::assayDataElement(object, "copynumber")
  if (scale)
    x <- t(t(x) / base::colMeans(x, na.rm=TRUE))

  vars = apply(x, 2L, function(y) matrixStats::varDiff(y, trim=trim, na.rm=TRUE))
  Biobase::pData(object)[["observed.variance"]] = vars

  return(object)
}

#' readData
#' 
#' Read bam files into a QDNAseq object.
#'
#' \code{readData} returns a QDNAseq copynumber object containing all the cells specified in the bamfiles argument.
#' This is a wrapper function around QDNAseq functionality to read in bam files.
#'
#' @param bamfiles A character vector of full (BAM) file paths.
#' @param binSize Numeric See the QDNAseq package for valid values for binSize. 
#' @param species Character We support Human and Mouse.
#' @param filterChromosomes A vector of chromosome names. By default, the mitochondrial DNA are filtered. 
#' @param extendedBlacklisting filter additional bins based on low read counts (centromer regions)
#'
#' @return A QDNAseq object containing the reads from cells specified through the bamfiles argument
#'
#' @export
readData <- function(bamfiles, binSize, species = "Human", filterChromosomes=c("MT"), genome="hg19",
                     extendedBlacklisting=TRUE,...) {
  
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)
  require(methods, quietly=TRUE, warn.conflicts = FALSE) # Rscript bugfix for QDNAseq (https://github.com/ccagc/QDNAseq/issues/49)
  require(GenomicRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(Rsamtools, quietly=TRUE, warn.conflicts = FALSE)
  require(dplyr, quietly=TRUE, warn.conflicts = FALSE)
  require(readr, quietly=TRUE, warn.conflicts = FALSE)
  
  transformFun = "none" # non-linear transformations are not supported
  # QDNAseq only supports these bin sizes out of the box
  bins = c(1, 5, 10, 15, 30, 50, 100, 500, 1000)
  stopifnot(binSize %in% bins)
  

  # required for smaller binsizes
  options(future.globals.maxSize= 4096*1024^2)
  
  if(species == "Human"){
    if(!(genome %in% c("GRCh37", "GRCh38", "hg19", "hg38"))){
      stop("Genome not supported")
    }
    
    if(genome == "hg19"){
      genome = "GRCh37"
    }
    
    if(genome == "hg38"){
      genome = "GRCh38"
    }
    
    if(genome == "GRCh37"){
      require(QDNAseq.hg19, quietly=TRUE, warn.conflicts = FALSE)
      binannotations = paste0("hg19.", binSize, "kbp.SR50")
      data(list = c(binannotations), package="QDNAseq.hg19")
    }else{
      if(genome == "GRCh38"){
        require(QDNAseq.hg38, quietly=TRUE, warn.conflicts = FALSE)
        binannotations = paste0("hg38.", binSize, "kbp.SR50")
        data(list = c(binannotations), package="QDNAseq.hg38")
      }
    }
    assign("bins", get(binannotations))
    # bins <- QDNAseq::getBinAnnotations(binSize=binSize)
  } else {
    if(species == "Mouse"){
      require(QDNAseq.mm10, quietly=TRUE, warn.conflicts = FALSE)
      binannotations = paste0("mm10.", binSize, "kbp.SR50")
      data(list = c(binannotations), package="QDNAseq.mm10")
      assign("bins", get(binannotations))
      # bins <- QDNAseq::getBinAnnotations(binSize=binSize, genome = "mm10")
    }else {
      stop("Species is not supported")
    }
  }
  
  isPairedReads = sapply(Rsamtools::BamFileList(bamfiles), Rsamtools::testPairedEndBam)
  if(!(all(isPairedReads) || all(!isPairedReads))){
    stop("Mixing single-end and paired-end bam files is not supported!")
  }
  isPaired = all(isPairedReads)
  

  
  if (all(endsWith(bamfiles, ".bam"))){
    readCounts = QDNAseq::binReadCounts(bins, bamfiles=bamfiles, isDuplicate = FALSE,
                                        isPaired = isPaired, isProperPair=NA,
                                        isFirstMateRead=ifelse(isPaired, TRUE, NA),
                                        isSecondMateRead=ifelse(isPaired, FALSE, NA))
    rawCounts = Biobase::assayDataElement(readCounts, "counts")
    Biobase::pData(readCounts)["bampath"] = dirname(bamfiles)
    Biobase::pData(readCounts)["bamfile"] = basename(bamfiles)
    Biobase::pData(readCounts)["isPaired"] = isPaired
  } else {
    stop("Check input file names, only bam files supported")
  }
  Biobase::pData(readCounts)["binSize"] = binSize
  Biobase::pData(readCounts)["species"] = species
  
  # Compute coverage (as an alternative to reads per cell metric)
  compute_read_size <- function(bamFile){
    yieldSize(bamFile) <- 10000
    # open(bamFile)
    reads = scanBam(bamFile)[[1]]$seq
    
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    
    return(getmode(reads@ranges@width))
  }
  read_size = sapply(Rsamtools::BamFileList(bamfiles), compute_read_size)
  Biobase::pData(readCounts)["read_size"] = as.numeric(read_size)
  
  if(species == "Human"){
    Biobase::pData(readCounts)["genome"] = genome
    coverage = (readCounts@phenoData@data$used.reads * (as.numeric(isPaired) + 1) * read_size) / 3.2e9
    Biobase::pData(readCounts)["coverage"] = coverage
  }
  
  if(species == "Human" & extendedBlacklisting){
    if(binSize < 10) stop("blacklisting doesn't support binsizes smaller than 10kb")
    # this is because our blacklisted regions have a resolution of 10kb
    readCounts = readFilter(readCounts, genome=genome)
    extendedFilter = Biobase::fData(readCounts)[["use"]]
  }else{
    extendedFilter = Biobase::fData(readCounts)[["use"]] | base::startsWith(rownames(readCounts), "X:") | base::startsWith(rownames(readCounts), "Y:")
  }
  
  # create QNDAseq object ####
  readCounts <- suppressWarnings(QDNAseq::estimateCorrection(readCounts, verbose=FALSE))
  readCountsFiltered <- QDNAseq::applyFilters(readCounts, residual=FALSE, blacklist=TRUE, mappability = NA, chromosomes = filterChromosomes, verbose=FALSE)
  Biobase::fData(readCountsFiltered)[["use"]] = extendedFilter & !is.na(Biobase::fData(readCounts)[["gc"]])
  fdat = Biobase::fData(readCountsFiltered)
  if(binSize %in% c(30, 50, 100, 500, 1000)){
    repTime = readRDS(file.path(BASEDIR, "data/replicationTiming/replicationTiming_per_binSize.RDS"))[[as.character(binSize)]]
    fdat$replicationTiming = repTime$replicationTime
  }

  copyNumbers <- new("QDNAseqCopyNumbers", bins=fdat,
                     copynumber=Biobase::assayDataElement(readCountsFiltered, "counts"),
                     phenodata=Biobase::pData(readCountsFiltered))
  
  cn = Biobase::assayDataElement(copyNumbers, "copynumber")
  keep_idx = binsToUseInternal(copyNumbers)
  cn[!keep_idx] = NA

  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "copynumber", cn)

  # add value for observed variance (sqrt of value in plots)
  copyNumbers = addObservedVariance(copyNumbers)

  # remove dummy variables - we do not apply GC correction at this stage
  Biobase::pData(copyNumbers)[["loess.span"]] = NULL
  Biobase::pData(copyNumbers)[["loess.family"]] = NULL

  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "calls", assayDataElement(copyNumbers, "copynumber"))
  copyNumbers = Biobase::assayDataElementReplace(copyNumbers, "probdloss", rawCounts)
  
  if(species == "Human"){
    Biobase::experimentData(readCounts)@other = list(species="Human", genome=genome)
  }else{
    Biobase::experimentData(readCounts)@other = list(species=species)
  }
  
  return(copyNumbers)
}

#' \code{readFilter} apply an updated bin filter to exclude low quality and unreliable bins from the estimation
#'
#' @param readCounts QDNAseq object
#' 
#' @return QDNAseq object with updated bin usability
readFilter <- function(readCounts, genome="GRCh37"){

  stopifnot(genome %in% c("GRCh37", "GRCh38"))
  if(genome == "GRCh37"){
    gr = createGR(readCounts[,1,drop=FALSE])
    excluded_regions_bed = readr::read_tsv(file.path(BASEDIR, "data/blacklisting/final_exclude_regions_hg19.bed"),
                                           col_names = c("chromosome", "start", "end", "A" ,"B", "C"), col_types="ciiccc")
    excluded_regions = GRanges(seqnames=excluded_regions_bed$chromosome,
                               ranges = IRanges(start=excluded_regions_bed$start+1, end=excluded_regions_bed$end),
                               seqinfo = seqinfo(gr))
  }
  if(genome == "GRCh38"){
    warning("blacklist has been optimized for GRCh37 only")
    stop("Not supported")
  }

  hits <- findOverlaps(excluded_regions, gr)

  # Compute percent overlap and filter the hits:
  overlaps <- pintersect(excluded_regions[queryHits(hits)], gr[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr[subjectHits(hits)])
  hits <- hits[percentOverlap > 0.25]

  pass_qc = Biobase::fData(readCounts)[["use"]] | base::startsWith(rownames(readCounts), "X:") | base::startsWith(rownames(readCounts), "Y:")
  pass_qc[subjectHits(hits)] = FALSE
   
  Biobase::fData(readCounts)[["use"]] = pass_qc
  return(readCounts)
}


#' readFlagstat
#' 
#' Read flagstat file to get additional information on read quality
#'
#' @param bamfiles A character vector of (BAM) file names.
#'
#' @return a data frame with the flagstat information.
#'
#' @export
readFlagstat <- function(bamfiles){
  
  read_flagstat_sub <- function(fh){
    # Read flagstat results
    require(stringr, quietly=TRUE, warn.conflicts = FALSE)
    f <- readLines(fh)
    a <- data.frame( as.character(sub("\\.flagstat$", "", base::basename(fh))),
                     as.numeric(stringr::word(f[1]) )
                     , as.numeric(stringr::word(f[2]) )
                     , as.numeric(stringr::word(f[4])  )
                     , as.numeric(stringr::word(f[5])  )
                     , as.numeric(stringr::word(f[6])  )
                     , as.numeric(stringr::word(f[9])  )
                     , as.numeric(stringr::word(f[10]) )
                     , as.numeric(stringr::word(f[11]) )
                     , as.numeric(stringr::word(f[12]) )
                     , as.numeric(stringr::word(f[13]) )
    )
    names(a) <- NULL
    colnames(a) <- c("name", "qcTotal", "qcSecondary", "qcDuplicates", "qcMapped", "qcPairedSequencing", "qcProperlyPaired", "qcPairedWithItself",
                     "qcSingletons", "qcMateMapped", "qcMateMappedmapQ5")
    return(a)
  }
  
  fh = bamfiles
  fh_list = c()
  for (i in 1:length(fh)){
    if(endsWith(fh[i], ".bam")){
      fh_list <- c(fh_list, fh[i])
    }else{
      if(dir.exists(fh[i])){
        files <- list.files(path = fh[i], pattern = "\\.bam$")
        files = paste0(fh[i], files)
        fh_list <- c(fh_list, files)
      }
    }
  }
  
  flagstat_description = list()
  counter = 1
  flaglist = sub("\\.bam$", ".flagstat", fh_list)
  for(i in 1:length(flaglist)){
    if(file.exists(flaglist[i]) && !dir.exists(flaglist[i])){
      flagstat_description[[counter]] = read_flagstat_sub(flaglist[i])
      counter = counter + 1
    }
  }
  
  result <- do.call("rbind", flagstat_description)
  return(result)
}




readPosition <- function(object, filePaths, selectRegion=c(as.character(c(1:22)), "X"), debug=FALSE){
  
  require(Biobase)
  require(QDNAseq)
  require(MASS)
  stopifnot("name" %in% colnames(Biobase::pData(object)))
  
  readPositionSingle <- function(obj, fP){
    valid = binsToUseInternal(obj)
    high_quality_regions = createGR(obj[,1])
    high_quality_regions = high_quality_regions[valid]
    chromosome = decode(seqnames(high_quality_regions))
    copynumber = Biobase::assayDataElement(obj, "copynumber")[valid]

    read_position_table = readr::read_tsv(paste0(sub('\\.bam$', '', fP), ".position.tsv"), col_types = "cdd")
    
    read_table = GRanges(read_position_table$chromosome,
                     IRanges(read_position_table$start, read_position_table$start+abs(read_position_table$length-1)),
                     length=abs(read_position_table$length),
                     dist=c(diff(read_position_table$start), Inf))
    
    if(!isSorted(read_table)){
      warning(paste0("Read table is not sorted - ", fP))
    }
    
    hits = findOverlaps(high_quality_regions, read_table)
    hq = read_table[subjectHits(hits)]
    n_overlap = countOverlaps(hq, hq)
    
    temp = dplyr::tibble(element = queryHits(hits), dists = read_table[subjectHits(hits)]$dist,
                         n_overlap = n_overlap,
                         chrom = chromosome[queryHits(hits)], copy=copynumber[queryHits(hits)]) %>%
      dplyr::filter(dists >= 0)
    
    
    compute_overlaps <- function(diffs, threshold){
      pos = diffinv(diffs)
      dis = as.matrix(dist(pos))
      n = length(pos)-1
      overlaps = sapply(1:n, function(index) sum(dis[index,(index+1):(n+1)] <= threshold))
      return(sum(overlaps))
    }
    
    filter_dists <- function(x){
      return(x[x < quantile(x, 0.99)])
    }
    
    dens = temp %>% dplyr::filter(chrom %in% sub("chr", "", selectRegion), copy > 0) %>% dplyr::group_by(chrom, element, copy) %>%
      dplyr::summarise(n=n(),
                       p_geom = MASS::fitdistr(dists[1:(length(dists)-1)], densfun = "geometric")[["estimate"]],
                       p_geom_censored = MASS::fitdistr(filter_dists(dists[1:(length(dists)-1)]), densfun = "geometric")[["estimate"]],
                       elem_overlap = mean(n_overlap, na.rm=TRUE), 
                       o0 = compute_overlaps(dists[1:(length(dists)-1)], 0),
                       o36 = compute_overlaps(dists[1:(length(dists)-1)], 36),
                       o125 = compute_overlaps(dists[1:(length(dists)-1)], 125),
                       .groups = 'keep') %>%
      dplyr::group_by(chrom, copy) %>%
      dplyr::summarise(n_elem = n(),
                       mean_overlap = mean(elem_overlap, na.rm=TRUE), median_overlap = median(elem_overlap, na.rm=TRUE),
                       nz_36 = sum(o36 != 0), nz_125 = sum(o125 != 0), nz_0 = sum(o0 != 0),
                       mgp = mean(p_geom, na.rm=TRUE), mgpc = mean(p_geom_censored, na.rm=TRUE), medgp = median(p_geom, na.rm=TRUE), medgpc = median(p_geom_censored, na.rm=TRUE),
                       mu0 = mean(o0), sum0 = sum(o0), sumf0 = sum(o0[o0 < quantile(o0, c(0.99))]), muf0 = mean(o0[o0 < quantile(o0, 0.99)]),
                       mu36 = mean(o36), sum36 = sum(o36), sumf36 = sum(o36[o36 < quantile(o36, c(0.99))]), muf36 = mean(o36[o36 < quantile(o36, 0.99)]),
                       mu125 = mean(o125), sum125 = sum(o125), sumf125 = sum(o125[o125 < quantile(o125, c(0.99))]), muf125 = mean(o125[o125 < quantile(o125, 0.99)]),
                       .groups = 'keep') %>% dplyr::ungroup()

    dens$name = Biobase::pData(obj)[["name"]]
    dens$rpc = Biobase::pData(obj)[["rpc"]]
    
    return(dens)
  } 
  
  densTable = lapply(1:ncol(object), function(index) readPositionSingle(object[, index], filePaths[index]))
  densTable = do.call("rbind", densTable) %>% dplyr::ungroup()
  
  meta = base::data.frame(labelDescription = rep(NA, times=dim(densTable)[[2]]))
  rownames(meta) = colnames(densTable)
  outProtData <- new("AnnotatedDataFrame",
                      data=densTable, 
                      varMetadata = meta)
  Biobase::protocolData(object) = outProtData

  return(object)
}


#' initial_gc_correction
#'
#' \code{initial_gc_correction} corrects segmentedCounts object for gc and map variation
#'
#' @param segmentedCounts QDNAseq object with segmented slot
#' @param debug Boolean flag to add debug information to output
#' 
#'
#' @return A QDNAseq object with raw reads from call slot corrected and added to copynumber slot
#'
#' @export
initial_gc_correction <- function(segmentedCounts, debug=FALSE){
  suppressPackageStartupMessages(require(mgcv, quietly=TRUE, warn.conflicts = FALSE))
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  
  iterate_initial_gc_correction <- function(object){
    
    stopifnot(dim(object)[[2]] == 1)
    n_bins = dim(object)[[1]]
    valid = binsToUseInternal(object)

    reads = Biobase::assayDataElement(object[valid, 1], "calls")
    segs = Biobase::assayDataElement(object[valid, 1], "segmented")
    gc = Biobase::fData(object[valid, 1])[["gc"]]
    map = Biobase::fData(object[valid, 1])[["mappability"]]
    stopifnot(length(gc) == length(reads))

    # log transform seg level for better stability
    dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
                     reads_normalized=as.vector(reads)/as.vector(segs),
                     gc=gc, map=map)
    
    ## Build models
    if(all(map == 100)){
      formula.gam = formula(reads ~ segs + s(gc))
      newdata=data.frame(segs = dat$segs, gc=dat$gc)
    }else{
      formula.gam = formula(reads ~ segs + s(gc, map))
      newdata=data.frame(segs = dat$segs, gc=dat$gc, map=dat$map)
    }
    dat_model = dat %>% dplyr::filter(segs > 0)
    model <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    terms = predict(model, newdata=newdata, type="terms")
    resp = predict(model, newdata=newdata, type="response")
    # quantile(resp)
    # max_reads = sort(resp, decreasing = TRUE)[1:10]
    # dat[resp %in% max_reads,]
    # sum(dat$reads) - sum(round(resp))
    
    
    # baseline without gc vs gc-map model
    # formula.gam.baseline = formula(reads ~ segs + s(gc))
    formula.gam.baseline = formula(reads ~ segs)
    model.baseline <- mgcv::gam(formula.gam.baseline, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    
    # formula.gam.both = formula(reads ~ segs + s(gc) + s(map))
    # formula.gam.both = formula(reads ~ segs + s(gc, map))
    model.baseline.both <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    
    av = anova(model.baseline, model.baseline.both, test="Chisq")
    ## end baseline anova
    
    ## Remove covariate effect
    
    # response = exp(terms[,1] + terms[,2] + coef(model)[1])
    # stopifnot(all(abs(response - resp) < 1e-6))
    # mean(terms[,2])
    response_gc_and_map = exp(terms[,2])
    residuum = dat$reads * (1/response_gc_and_map)
    
    dat$pred_seg = terms[,1]
    dat$pred_cov = terms[,2]
    dat$residuum = residuum

    stopifnot(length(residuum) == length(reads))
    corrected = round(residuum, digits=0)
    dat$corrected_reads = corrected
    negative_elements = sum(corrected < 0, na.rm=TRUE)
    corrected[corrected < 0] = 0
    stopifnot(all(!is.na(corrected)))
    
    corrected_reads = Biobase::assayDataElement(object, "calls")
    stopifnot(length(corrected_reads[valid]) == length(corrected))
    corrected_reads[valid,1] = corrected
    
    corrected_reads_residuum = Biobase::assayDataElement(object, "calls")
    corrected_reads_residuum[valid,1] = residuum
    
    # sum(dat$reads) - sum(dat$corrected)
    # table(dat$corrected != dat$reads)
    # pl1 = ggplot(data=dat) + geom_point(aes(x=gc, y=reads_normalized)) +
    #   # geom_point(aes(x=gc, y=reads_normalized, color=as.factor(state)), size=0.1) +
    #   geom_point(aes(x=gc, y=reads_normalized), size=0.1) +
    #   geom_smooth(aes(x=gc, y=reads_normalized), formula="y ~ x", method="loess", span=0.7, color="red", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot()
    # pl1
    
    # pl2 = ggplot(data=dat1) + geom_point(aes(x=gc, y=reads)) +
    #   geom_point(aes(x=gc, y=reads, color=as.factor(state)), size=0.1) +
    #   geom_smooth(aes(x=gc, y=reads), formula="y ~ x", method="loess", span=0.5, color="red", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot()
    # pl2

    # ggplot(data=dat1) + geom_histogram(aes(x=ratio1)) + theme_cowplot() + coord_cartesian(xlim=c(-0.3, 0.3))
    # dat1$corrected_reads1 = dat1$reads - round(dat1$ratio1 * dat1$reads, digits=0)
    # dat1$corrected_reads2 = dat1$reads - round(dat1$ratio2 * dat1$reads, digits=0)
    # ggplot(data = dat1 %>%
    #     dplyr::select(reads, corrected_reads2) %>% tidyr::gather(condition, r, c(corrected_reads2, reads))) +
    #     geom_histogram(aes(x=r, fill=condition), binwidth = 1, position="identity", alpha=0.75) +
    #     theme_cowplot()
    # 
    # pl3 = ggplot(data=dat1) + geom_point(aes(x=gc, y=corrected_reads2)) +
    #   geom_smooth(aes(x=gc, y=corrected_reads2), formula="y ~ x", method="loess", span=0.75, color="red", se=FALSE) +
    #   geom_smooth(aes(x=gc, y=reads), formula="y ~ x", method="loess", span=0.75, color="blue", se=FALSE) +
    #   theme_cowplot()
    # pl3
    
    # # #Visualize mgcv model
    # options(rgl.useNULL=TRUE)
    # library(mgcViz) # conda dependencies: rgl, qgam, make surn KernSmooth is correct version (conda-forge)
    # 
    # model.viz <- getViz(model)
    # o <- plot( sm(model.viz, 1) )
    # o + theme_pubr()
    
    ## Read correction
    # dat1$corrected_reads = dat1$reads - round(ratio * dat1$reads, digits=0)
    
    # # smoothing per cn state
    # pl3 = ggplot(data=dat1) + geom_point(aes(x=gc, y=corrected_reads, color=as.factor(state))) + 
    #   geom_smooth(aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="red", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==1), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==2), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_smooth(data=dat1 %>% dplyr::filter(state==3), aes(x=gc, y=corrected_reads), formula="y ~ x", method="loess", span=span, color="blue", se=FALSE) +
    #   geom_vline(xintercept = c(37, 43), color="red") + theme_cowplot() #+ coord_cartesian(ylim=c(200, 300)
    # pl3

    # ggplot(data = dat %>%
    #     dplyr::select(reads, corrected_reads) %>% tidyr::gather(condition, r, c(corrected_reads, reads))) +
    #     geom_histogram(aes(x=r, fill=condition), binwidth = 1, position="identity", alpha=0.75) +
    #     coord_cartesian(xlim=c(200, 500)) +
    #     theme_cowplot()
    
    object = Biobase::assayDataElementReplace(object, "copynumber", matrix(data=corrected_reads, ncol=1))
    desc = Biobase::pData(object)
    
    desc$gc_correction = TRUE
    desc$gc_correction.negative_elements = negative_elements
    desc$gc_correction.av.deviance = av[["Deviance"]][2]
    desc$gc_correction.av.p_value = av[["Pr(>Chi)"]][[2]]
    desc$gc_correction.dev.expl = summary(model)[["dev.expl"]]
    desc$gc_correction.alpha = 1/model$family$getTheta(TRUE)
    desc$gc_correction.delta_n_reads = sum(Biobase::assayDataElement(object, "calls"), na.rm = TRUE) - sum(Biobase::assayDataElement(object, "copynumber"), na.rm=TRUE)
    Biobase::pData(object) = desc
    return(list("obj"=object, list("name"=colnames(object), "dat"=dat, "model"=model)))
  }
  
  results = future.apply::future_lapply(lapply(colnames(segmentedCounts), function(x) segmentedCounts[, x]), iterate_initial_gc_correction)
  new_object = combineQDNASets(sapply(results, "[[", 1))
  debugInformation = sapply(results, "[[", 2)
  
  if(debug){
    return(list("new_object"=new_object, "debugInformation"=debugInformation) )
  }else{
    return(new_object)  
  }
}


#' estimate_gc_correction
#'
#' \code{estimate_gc_correction} estimate gc correction per cell
#'
#' @param object a QDNAseq object with segmented slot
#'
#' @return a n_bins x n_cells matrix with offset for gc/map correction
#'
#' @export
estimate_gc_correction <- function(segmentedCounts){
  
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)  
  n_cells = dim(segmentedCounts)[[2]]
  valid = binsToUseInternal(segmentedCounts)
  
  iterate_processing <- function(x){
    suppressPackageStartupMessages(require(QDNAseq, quietly = TRUE, warn.conflicts = FALSE))
    obj = segmentedCounts[, x]
    
    stopifnot(dim(obj)[[2]] == 1)
    n_bins = dim(obj)[[1]]
    valid = binsToUseInternal(obj) & !is.na(Biobase::fData(obj)[["gc"]])
    
    reads = Biobase::assayDataElement(obj, "calls")[valid, 1, drop=TRUE]
    segs = Biobase::assayDataElement(obj, "segmented")[valid, 1, drop=TRUE]
    gc = Biobase::fData(obj[valid, 1])[["gc"]]
    map = Biobase::fData(obj[valid, 1])[["mappability"]]
    stopifnot(length(gc) == length(reads))
    
    # log transform seg level for better stability
    dat = data.frame(reads=as.vector(reads), segs=log(as.vector(segs) + 1),
                     reads_normalized=as.vector(reads)/as.vector(segs),
                     gc=gc, map=map)
    
    ## Build models
    if(all(map == 100)){
      formula.gam = formula(reads ~ segs + s(gc))
      newdata=data.frame(segs = dat$segs, gc=dat$gc)
    }else{
      formula.gam = formula(reads ~ segs + s(gc, map))  
      newdata=data.frame(segs = dat$segs, gc=dat$gc, map=dat$map)
    }
    dat_model = dat %>% dplyr::filter(segs > 0)
    model <- mgcv::gam(formula.gam, data=dat_model, family=mgcv::nb(theta = NULL, link = "log"), method = "REML")
    terms = predict(model, newdata=newdata, type="terms")
    
    gc_terms = terms[,2]
    cor = rep(0, n_bins)
    cor[valid] = exp(gc_terms)
    return(cor)
  }
  
  xs <- seq(1, dim(segmentedCounts)[[2]])
  results = future.apply::future_lapply(as.list(xs), iterate_processing)
  mat = t(do.call("rbind", results))
  return(mat)
}


#' combineQDNASets
#'
#' \code{combineQDNASets} combines multiple QDNAseq objects to a single object, useful for parallelisation across cells.
#'
#' @param ... list of QDNAseq objects
#'
#' @return A QDNAseq object containing all the input objects
#'
#' @export
combineQDNASets <- function(..., reduceMetadata=FALSE){
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  require(QDNAseq, quietly=TRUE, warn.conflicts = FALSE)


  x = ...[[1]]

  if(length(...) == 1){
    # idempotence
    return(x)
  }

  # keep information from x, append assayData and other data from other elements
  if("segmented" %in% names(Biobase::assayData(x))){
    outSegmented <- do.call(cbind, lapply(..., Biobase::assayDataElement, "segmented"))
  }
  outCopynumber <- do.call(cbind, lapply(..., Biobase::assayDataElement, "copynumber"))
  if("calls" %in% names(Biobase::assayData(x))){
    outCalls <- do.call(cbind, lapply(..., Biobase::assayDataElement, "calls"))
  }
  if("probdloss" %in% names(Biobase::assayData(x))){
    outProbdloss <- do.call(cbind, lapply(..., Biobase::assayDataElement, "probdloss"))
  }

  ## check pheno data is compatible and indicate source of error
  test = lapply(..., pData)
  n_names = unlist(lapply(..., BiocGenerics::colnames))
  n_items = unlist(lapply(test, length))
  n_max = max(n_items)
  if(any(n_items < n_max)){
    print(n_names[n_items == n_max])
    print("---")
    print(n_names[n_items < n_max])
    warning(paste0("INCONSISTENT METADATA -> CHECK ", sum(n_items < n_max), " OFFENDING SAMPLES"))
    if(reduceMetadata){
      commonNames = colnames(test[[1]])
      for(item in 1:length(test)){
        commonNames = intersect(commonNames, colnames(test[[item]]))
        if(length(setdiff(colnames(test[[item]]), commonNames)) > 0){
          print(paste0("Items in question: ", setdiff(colnames(test[[item]]), commonNames)))
        }
      }
    }else{
      stop("stop to avoid inconsistent metadata")
    }
  }else{
    commonNames = colnames(test[[1]])
  }

  meta = Biobase::varMetadata(x)
  meta = meta[rownames(meta) %in% commonNames,,drop=FALSE]
  outPhenoData <- new("AnnotatedDataFrame",
                      data=base::do.call(base::rbind, lapply(..., function(x){y = pData(x); return(y[, commonNames, drop=FALSE])})),
                      varMetadata = meta)
  
  outFeatureData <- new("AnnotatedDataFrame",
                        data=fData(x),
                        varMetadata = Biobase::fvarMetadata(x))

  outProtocolData = protocolData(x)
  if(dim(outProtocolData)[[2]] > 0){
    protDat = data.frame(do.call(base::rbind, lapply(..., function(x) Biobase::protocolData(x)@data)))
    rownames(protDat) = as.character(seq(1, dim(protDat)[[1]]))
    outProtocolData@data = protDat
  }

  result = new("QDNAseqCopyNumbers",
               bins=outFeatureData,
               copynumber=outCopynumber,
               phenodata=outPhenoData)
  if("segmented" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "segmented") <- outSegmented
  }
  if("calls" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "calls") <- outCalls
  }
  if("probdloss" %in% names(Biobase::assayData(x))){
    Biobase::assayDataElement(result, "probdloss") <- outProbdloss
  }
  Biobase::protocolData(result) <- outProtocolData
  
  return(result)
}

#' getSegTable
#'
#' \code{getSegTable} return a bed style representation of a QDNAseq object's segmentation
#'
#' @param object QDNAseq object with segmentation slot
#' @param trimLength Numeric number of bins to trim from start and end of segments
#'
#' @return a GRanges List object with n_samples elements, representing the individual samples segmentation profiles
#'
#' @export
getSegTable <- function(object, trimLength = 1, slot="segmented"){
  require(GenomicRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(IRanges, quietly=TRUE, warn.conflicts = FALSE)
  require(S4Vectors, quietly=TRUE, warn.conflicts = FALSE)
  require(Biobase, quietly=TRUE, warn.conflicts = FALSE)
  
  export_genomic_ranges_call <- function(sample){
    
    gr = retrieveSegmentation(object, sample, slot=slot)
    counts = retrieveAssayData(object, sample, "copynumber")
    
    n_segments = length(gr@ranges)

    if(!is.null(Biobase::assayDataElement(object, "calls"))){
      gr$reads = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){stats::median(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$var = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){stats::var(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$mean = smoothOperator(retrieveAssayData(object, sample, "calls"), gr, function(x){base::mean(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
      gr$sample_mean = base::mean(assayDataElement(object[, sample], "calls"), na.rm=TRUE)
      gr$sample_var = stats::var(assayDataElement(object[, sample], "calls"), na.rm=TRUE)
    }
    gr$copynumber = smoothOperator(retrieveAssayData(object, sample, "copynumber"), gr, function(x){stats::median(x, na.rm=TRUE)}, trimLength = trimLength, minLength=1)
    gr$support = smoothOperator(retrieveAssayData(object, sample, "copynumber"), gr, function(x){length(na.omit(x))}, trimLength = 0, minLength=1)
    gr$length = gr@ranges@width
    gr$sample = sample
    
    ## convert to genomic coordinates (bp coordinates)
    gra = createGR(object[,sample])
    ir_new = do.call("c", lapply(seqlevels(gr), function(x){
      grat = gra[seqnames(gra) == x]
      grt = gr[seqnames(gr) == x]
      
      if(length(grt) == 0){
        return(IRanges())
      }
      
      stopifnot(all(length(grat@ranges) >= grt@ranges@start))
      binsize = grat@ranges@width[[1]]
      irStart = c();irEnd = c(); irWidth = c();
      for(idx in 1:length(grt)){
        start_grt = grt@ranges@start[idx]
        end_grt = grt@ranges@start[idx] + grt@ranges@width[idx] - 1
        start = grat@ranges@start[start_grt]
        end = grat@ranges@start[start_grt] + sum(grat@ranges@width[start_grt:end_grt]) - 1
        width = end - start + 1
        irStart = c(irStart, start); irEnd = c(irEnd, end);
        irWidth = c(irWidth, width)
      }
      
      irnew = IRanges::IRanges(start = irStart,
                               end = irEnd, 
                               width = irWidth)
      
      return(irnew)
    }))
    
    gr@ranges = ir_new
    # end coordinate transformation
    
    return(gr)
  }
  
  result = lapply(colnames(object), export_genomic_ranges_call)
  return(GRangesList(result))
}

#' writeSegTable
#'
#' \code{writeSegTable} combines multiple QDNAseq objects to a single object, useful for parallelisation across cells.
#'
#' @param object QDNAseq object
#' @param filename filename to write bed file to
#' @param trimLength Numeric number of bins to trim from start and end of segments
#'
#' @export
writeSegTable <- function(object, filename, trimLength=1, minLength=NULL){
  # minLength is ignored, keep for backward compatibility
  
  gr = getSegTable(object, trimLength=trimLength)
  
  df = do.call(rbind, lapply(gr, function(gr1){return(
    dplyr::bind_cols(
      data.frame(seqnames=seqnames(gr1),
                 start=start(gr1),
                 end=end(gr1),
                 width=width(gr1)),
      data.frame(gr1@elementMetadata)))}))
  
  write.table(df, file=filename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
}
