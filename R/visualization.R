# Copyright 2022, Michael Schneider, All rights reserved.
## Visualization methods

#' plotCounts
#'
#' \code{plotCounts} wrapper around QDNAseq::plot to plot raw count data (calls) and segmentation
#'
#' @param object a QDNAseq object for a single sample
#' @param ylim tuple of minimum and maximum count value to show. Default: c(0, 90)
#' @param ybreaks format of y axis breaks
#' @param readinfo Boolean add information about sample (as title)
#' @param ylab y axis label
#' @param chromosome_break_label x axis names to show (default is medium dense, with not all chromosomes shown)
#' @param main title of plot
#'
#' @return A ggplot plot object
#'
#' @export
plotCounts <- function(object, ylim=c(0, 90),  ybreaks=NULL, ylab="read counts", main=NULL,
                       copynumber=FALSE, correction=FALSE, addSegmentation=FALSE, segmentationColor="#E69F00",
  chromosome_break_label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                             "11", "", "13", "", "15", "", "17", "", "19", "",
                             "21", "", "X", "Y")){
  
  
  require(S4Vectors, warn.conflicts = FALSE)
  require(ggplot2, warn.conflicts = FALSE)
  require(cowplot, warn.conflicts = FALSE)
  require(scales, warn.conflicts = FALSE)

  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  if(length(unique(seqName)) != 24){
    chromosome_break_label = base::rle(seqName)[["values"]]
  }
    
  if(copynumber){
    dat = Biobase::assayDataElement(object, "copynumber")
  }else{
    dat = Biobase::assayDataElement(object, "calls")
  }
  
  if(addSegmentation){
    dat2 = Biobase::assayDataElement(object, "segmented") * Biobase::pData(object)[["rpc"]]
  }else{
    dat2 = NA
  }

  if(correction){
    stopifnot("probdloss" %in% names(object@assayData))
    stopifnot("segmented" %in% names(object@assayData))
    covariate_info = estimate_gc_correction(object)
    dat = Biobase::assayDataElement(object, "calls") * (1/covariate_info)
  }

  obj = new("QDNAseqCopyNumbers",
            bins=Biobase::featureData(object),
            copynumber=dat,
            phenodata=Biobase::pData(object))

  cellname = Biobase::pData(object)[["name"]]
  title = unlist(strsplit(cellname, "_"))
  if(length(title) >= 3){
    title = paste0(title[[1]], " - ", title[[3]])
  } else {
    title = cellname
  }
  if(is.null(main)){
    main = title
  }
  
  if(is.null(ybreaks)){
    ybreaks=scales::pretty_breaks()
  }
  
  n_bins = dim(object)[[1]]
  x = seq(1, n_bins)
  
  location = rownames(object)
  seqName = unlist(lapply(strsplit(location, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  rl = S4Vectors::Rle(chromosome)
  if(length(rl@values) == 1){
    chromosome_break_position = c(rl@lengths / 2)
    chromosome_boundaries = c()
  }else{
    rl_cumsum = cumsum(rl@lengths)
    chromosome_boundaries = rl_cumsum[1:(length(rl_cumsum)-1)]
    chromosome_break_position = c(0, cumsum(rl@lengths)[1:(length(rl_cumsum)-1)]) + ceiling(rl@lengths/2)
  }
  
  pointsize=0.1
  d = data.frame(x=x, y=as.numeric(dat), y2=as.numeric(dat2))
  p = ggplot(data=d) + geom_point(aes(x=x, y=y), color="black", size=pointsize) +
    {if(main!="")ggtitle(main)}+
    ylab("read counts") + # geom_point(aes(x=x, y=y),size=pointsize, alpha=0.3) + 
    geom_vline(xintercept = chromosome_boundaries, linetype = "dashed", size=0.5, alpha=1.0) +
    {if(addSegmentation)geom_point(aes(x=x, y=y2), color=segmentationColor, size=pointsize)} +
    scale_x_continuous("chromosome", breaks=chromosome_break_position, labels=chromosome_break_label, expand=c(0,0)) +
    scale_y_continuous(breaks=ybreaks, limits=ylim) +
    theme_cowplot()
  
  return(p)
}

#' plotQDNAseqCounts
#'
#' \code{plotQDNAseqCounts} wrapper around QDNAseq::plot to plot raw count data (calls) and segmentation
#'
#' @param object a QDNAseq object for a single sample
#' @param rescale rescale segmentation to fit the original raw read counts (for visualization purposes)
#' @param ylim tuple of minimum and maximum count value to show. Default: c(0, 90)
#' @param yaxp y axis tick frequency. Default: c(0, 90, 9)
#' @param ylab y axis title
#' @param ... you can give other parameters to QDNAseq::plot when plotting
#'
#' @return A QDNAseq plot object
#'
#' @export
plotQDNAseqCounts <- function(object, rescale=FALSE, copynumber=TRUE, ylim=c(0, 90), yaxp=c(0, 90, 9), ylab="read counts", ...){
  stopifnot(dim(object)[[2]] == 1)

  if(!is.null(Biobase::assayDataElement(object, "calls"))){
    dat = Biobase::assayDataElement(object, "calls")
  }else{
    dat = Biobase::assayDataElement(object, "copynumber")
  }
  obj = new("QDNAseqCopyNumbers",
            bins=Biobase::featureData(object),
            copynumber=dat,
            phenodata=Biobase::pData(object))
  if(!is.null(Biobase::assayDataElement(object, "segmented"))){
    if(rescale){
      stopifnot(!is.null(object@phenoData@data[["scale"]]))
      scale = object@phenoData@data[["scale"]]
    }else{
      scale = rep(1, dim(object)[[2]])
    }
    
    if(copynumber){
      dat = (diag(1/as.numeric(unname(scale)), nrow=dim(object)[[2]]) %*% t(object@assayData$copynumber))
      Biobase::assayDataElement(obj, "segmented") <- t(ifelse(is.na(dat), -1, dat))
    }else{
      dat = (diag(1/as.numeric(unname(scale)), nrow=dim(object)[[2]]) %*% round(t(object@assayData$segmented), digits = 0))
      Biobase::assayDataElement(obj, "segmented") <- t(ifelse(is.na(dat), -1, dat))
    }
  }
  cellname = Biobase::pData(object)[["name"]]
  title = unlist(strsplit(cellname, "_"))
  if(length(title) >= 3){
    title = paste0(title[[1]], " - ", title[[3]])
  } else {
    title = cellname
  }

  QDNAseq::plot(obj, logTransform=FALSE, ylim=ylim, yaxp=yaxp, ylab=ylab, main=title, ...)
}

#' plotQDNAseqCopynumber
#'
#' \code{plotQDNAseqCopynumber} wrapper around QDNAseq::plot to plot copy number profile for a single sample
#'
#' @param object a scaled QDNAseq object for a single sample
#' @param round apply rounding to segmentation, i.e. fix all segmented values at integer values (might lead to merging of some segments)
#' @param copynumber use copynumber slot for plotting data (TRUE) or alternative use (scaled) calls slot (FALSE)
#' @param ylim tuple of minimum and maximum count value to show. Default: c(0, 9)
#' @param yaxp y axis tick frequency. Default: c(0, 9, 9)
#' @param ylab y axis title
#' @param main custom title
#' @param ... you can give other parameters to QDNAseq::plot when plotting
#'
#' @return A QDNAseq plot object
#'
#' @export
plotQDNAseqCopynumber <- function(object, round=FALSE, copynumber=TRUE, ylim=c(0, 9), yaxp=c(0, 9, 9), ylab="absolute copy number", main=NULL, ...){
  stopifnot(dim(object)[[2]] == 1)
  stopifnot(!is.null(Biobase::pData(object)[["scale"]]))

  scale = Biobase::pData(object)[["scale"]]
  stopifnot(!is.null(Biobase::assayDataElement(object, "calls")))
  dat = scale * Biobase::assayDataElement(object, "calls")

  obj = new("QDNAseqCopyNumbers",
            bins=Biobase::featureData(object),
            copynumber=dat,
            phenodata=Biobase::pData(object))

  if(copynumber){
    if(round){
      dat = round(Biobase::assayDataElement(object, "segmented"), digits = 0)
      Biobase::assayDataElement(obj, "segmented") <- ifelse(is.na(dat), -1, dat)
    }else{
      dat = Biobase::assayDataElement(object, "copynumber")
      Biobase::assayDataElement(obj, "segmented") <- ifelse(is.na(dat), -1, dat)
    }
  }

  cellname = Biobase::pData(obj)[["name"]]
  title = unlist(strsplit(cellname, "_"))
  if(length(title) >= 3){
    title = paste0(title[[1]], " - ", title[[3]])
  } else {
    title = cellname
  }
  if(is.null(main)){
    main = title
  }

  QDNAseq::plot(obj, logTransform=FALSE, ylim=ylim, yaxp=yaxp, ylab=ylab, main=main, ...)
}


#' plotCopynumber
#'
#' \code{plotCopynumber} wrapper around ggplot to plot copy number profile for a single sample
#'
#' @param object a scaled QDNAseq object for a single sample
#' @param round apply rounding to segmentation, i.e. fix all segmented values at integer values (might lead to merging of some segments)
#' @param copynumber use copynumber slot for plotting data (TRUE) or alternative use (scaled) calls slot (FALSE)
#' @param ylim tuple of minimum and maximum count value to show. Default: c(0, 9)
#' @param ybreaks format of y axis breaks
#' @param readinfo Boolean add information about sample (as title)
#' @param xlab x axis label
#' @param ylab y axis label
#' @param chromosome_break_label x axis names to show (default is medium dense, with not all chromosomes shown)
#' @param main title of plot
#' @param copyColor color to use for copynumber slot
#' @param showUnique show separate color for unique events 
#' @param showMarker show symbol to indicate unique events
#' @param addCopynumber add separate lines for copynumber slot
#' @param addSegmentation add separate lines for segmentation slot segments
#' @param correction apply correction to data points plotted
#' @param plotSegmentation use segmentation slot to plot segments
#' @param alphaLevel alpha level for raw data points display
#'
#' @return A ggplot plot object
#'
#' @export
plotCopynumber <- function(object, showUnique=TRUE, round=FALSE, correction=FALSE, ylim=c(0, 10),
                           ybreaks=NULL, readinfo=TRUE, ylab="absolute copy number", xlab=NULL,
                           chromosome_break_label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                      "11", "", "13", "", "15", "", "17", "", "19", "",
                                                      "21", "", "X", "Y"),
                           main=NULL, copyColor="#E69F00", addSegmentation=FALSE, addCopynumber=TRUE,
                           showMarker=TRUE, alphaLevel=0.3, uniqueColor="#FFFFFF", segmentationColor="#E69F00", alphaLevelSeg=1.0, verbose=TRUE){

  require(S4Vectors, warn.conflicts = FALSE)
  require(ggplot2, warn.conflicts = FALSE)
  require(cowplot, warn.conflicts = FALSE)
  require(scales, warn.conflicts = FALSE)
  require(Biobase, warn.conflicts = FALSE)
  require(QDNAseq, warn.conflicts = FALSE)
  
  stopifnot(dim(object)[[2]] == 1)
  stopifnot(!is.null(assayDataElement(object, "segmented")))
  # stopifnot(!is.null(Biobase::pData(object)[["scale"]]))
  info = rownames(object)
  seqName = unlist(lapply(strsplit(info, ":"), `[`, 1))
  valid=binsToUseInternal(object)
  if(length(unique(seqName)) != 24){
    chromosome_break_label = base::rle(seqName)[["values"]]
  }
  
  if(is.null(ybreaks)){
    ybreaks=scales::pretty_breaks()
  }

  if(is.null(Biobase::pData(object)[["scale"]])){
    scale = 1
  }else{
    scale = Biobase::pData(object)[["scale"]]
  }

  # if(copynumber){
  #   calls = scale * Biobase::assayDataElement(object, "copynumber")
  # }else{
    stopifnot(!is.null(Biobase::assayDataElement(object, "calls")))
    calls = scale * Biobase::assayDataElement(object, "calls")
  # }
  
  if(correction){
    stopifnot("probdloss" %in% names(object@assayData))
    covariate_info = estimate_gc_correction(object)
    calls = Biobase::assayDataElement(object, "calls") * (1/covariate_info) * scale
  }

  n_bins = dim(object)[[1]]
  x = seq(1, n_bins)

  if(round){
    segments = round(Biobase::assayDataElement(object, "segmented"), digits = 0)
    copies = round(Biobase::assayDataElement(object, "copynumber"), digits = 0)
  }else{
    segments = Biobase::assayDataElement(object, "segmented")
    copies = Biobase::assayDataElement(object, "copynumber")
  }

  cellname = Biobase::pData(object)[["name"]]
  title = unlist(strsplit(cellname, "_"))
  if(length(title) > 3){
    if(base::startsWith(title[[4]], "SINCEL")){
      title = paste0(title[[1]], " - ", title[[3]], " - ", paste0(unlist(strsplit(title[[4]], "-"))[2:5], collapse="-"))
    }else{
      title = paste0(title[[1]], " - ", title[[3]])  
    }
  }else {
    title = cellname
  }
  
  if(is.null(main)){
    main = title
  }

  ## add read information to plot title
  if(is.logical(readinfo) && readinfo){
    if("rpc" %in% colnames(Biobase::pData(object))){
      rpc = base::format(round(Biobase::pData(object)[["rpc"]], digits=1), nsmall=0L, big.mark=",", digits=10)
    }else{
      rpc = "-"
    }
    if("used.reads" %in% colnames(Biobase::pData(object))){
      used_reads = base::format(Biobase::pData(object)[["used.reads"]], nsmall=0L, big.mark=",")
    }else{
      used_reads = "-"
    }
    main = paste0(main, " (", used_reads, " / ", rpc , ")")
  }
  if(!is.logical(readinfo) && readinfo == "FULL"){
    if("name" %in% colnames(Biobase::pData(object))){
      main = Biobase::pData(object)[["name"]]
    }else{
      warning("name missing in metadata")
    } 
  }

  location = rownames(object)
  seqName = unlist(lapply(strsplit(location, ":"), `[`, 1))
  chrOrder = base::rle(seqName)[["values"]]
  chromosome <- factor(seqName, levels=chrOrder)
  rl = S4Vectors::Rle(chromosome)
  rl_cumsum = cumsum(rl@lengths)
  chromosome_boundaries = rl_cumsum[1:(length(rl_cumsum)-1)]
  chromosome_break_position = c(0, cumsum(rl@lengths)[1:(length(rl_cumsum)-1)]) + ceiling(rl@lengths/2)
  isSingleChromosome = length(chromosome_break_label) == 1

  pointsize=0.1
  
  # Look at unique events / events that differ between segmentation and copynumber slot
  unique_index = round(as.numeric(segments), digits=0) != round(as.numeric(copies), digits=0)
  colorUnique = rep(copyColor, times=n_bins)
  #colorUnique[unique_index] = uniqueColor
  sizeUnique = rep(pointsize, times=n_bins)
  sizeUnique[unique_index] = 10 * pointsize
  
  # add symbols for unique events/deviation between copynumber and segmentation to top of plot
  unique_marker = rep(0, sum(valid))
  segment_vector = paste0(as.numeric(segments[valid]), "-", round(as.numeric(segments)[valid], digits=0) != round(as.numeric(copies)[valid], digits=0))
  rle_l = Rle(segment_vector)@lengths
  if(length(rle_l) > 1){
    unique_marker[cumsum(c(0, rle_l[1:(length(rle_l)-1)])) + floor(rle_l/2)] = 1
    unique_marker_position = rep(FALSE, n_bins)
    unique_marker_position[valid][unique_marker == 1 & unique_index[valid]] = TRUE
  }else{
    unique_marker_position = rep(FALSE, n_bins)
  }
  
  
  if(verbose){
    print("Unique items:")
    print(table(unique_index))
    print(sum(unique_marker_position))
  }
  
  if(is.na(sum(unique_marker_position)) | sum(unique_marker_position)==0){
    showMarker=FALSE
  }
  
  d = data.frame(x=x, y2=as.numeric(calls), y=as.numeric(segments), z=as.numeric(copies))
  p = ggplot(data=d) + 
    {if(showUnique)scale_colour_manual("KEYS", values = colorUnique)} +
    {if(main!="")ggtitle(main)} +
    geom_point(aes(x=x, y=y2),size=pointsize, alpha=alphaLevel) + ylab(ylab) +
    {if(addSegmentation)geom_point(aes(x=x, y=y), color=segmentationColor, size=pointsize, alpha=alphaLevelSeg)} +
    {if(addCopynumber)geom_point(aes(x=x, y=z), color=colorUnique, size=pointsize)} +
    {if(showMarker)geom_point(data=data.frame(x2=x[unique_marker_position], y2=(ylim[[2]]-1e-6)),
      aes(x=x2, y=y2), shape=25, size=3, fill="#D55E00")} +
    {if(!(is.logical(ybreaks) && ybreaks==FALSE))scale_y_continuous(breaks=ybreaks)} +
    coord_cartesian(ylim=ylim) +
    theme_cowplot()

  if(isSingleChromosome){
    if(is.null(xlab)){
      p = p + scale_x_continuous(paste0("chromosome ", chromosome_break_label), breaks=NULL, labels=chromosome_break_label, expand=c(0,0))
    }else{
      p = p + scale_x_continuous(xlab, breaks=NULL, labels=chromosome_break_label, expand=c(0,0))
    }
  }else{
    if(is.null(xlab)){
      p = p + geom_vline(xintercept = chromosome_boundaries, linetype = "dashed", size=0.5, alpha=1.0) +
        scale_x_continuous("chromosome", breaks=chromosome_break_position, labels=chromosome_break_label, expand=c(0,0))
    }else{
      p = p + geom_vline(xintercept = chromosome_boundaries, linetype = "dashed", size=0.5, alpha=1.0) +
        scale_x_continuous(xlab, breaks=chromosome_break_position, labels=chromosome_break_label, expand=c(0,0))
    }
  }

  return(p)
}

#' plotFit
#'
#' \code{plotFit} plot mean-variance relationship for single sample (or alternatively mean-alpha relationship)
#'
#' @param object a scaled QDNAseq object for a single sample
#' @param scale Boolean scale mean and variance with fit
#' @param alpha Boolean plot mean-alpha instead of mean-variance relationship
#' @param trimLength Numeric number of bins to remove from start and end of segments for estimates
#' @param minLength Numeric minimum segment size to include in fitting
#' @param limitPloidy only applied if scale is TRUE -> remove l_mean values larger than limitPloidy
#'
#' @return a plot object detailing the mean-variance relationship
#'
#' @export
plotFit <- function(object, scale=TRUE, alpha=FALSE, trimLength=1, minLength=10, limitPloidy=16){
  stopifnot(dim(object)[[2]] == 1)
  stopifnot(!is.null(Biobase::assayDataElement(object, "segmented")))

  require(ggplot2, quietly = TRUE, warn.conflicts = FALSE)

  valid = binsToUseInternal(object)
  counts = retrieveAssayData(object, 1, value="calls", valid=valid)
  segs = Biobase::assayDataElement(object, "segmented")
  gr = retrieveSegmentation(base::round(segs, digits = 0))

  if(scale){
    stopifnot(!is.null(Biobase::pData(object)[["scale"]]))
    scalar = Biobase::pData(object)[["scale"]]
    model = computeModel(counts, gr, scale=scalar, roundInteger=FALSE, limitPloidy=limitPloidy, minLength=minLength,trimLength=trimLength,debug=TRUE)
  }else{
    model = computeModel(counts, gr, scale=1.0, limitPloidy=Inf, minLength=minLength,trimLength=trimLength,debug=TRUE)
  }

  if(alpha){
    dat2 = data.frame(l_mean=model$l_mean_unscaled_alpha, l_alpha=model$l_alpha, l_weight=model$l_weight_alpha, l_chromosome=model$l_chromosome_alpha)
    sample_points = seq(0, max(1/dat2$l_mean), 0.0001)
    points_robust = model$df$meanvar.alpha.coefficient_combined_1 + (sample_points) * model$df$meanvar.alpha.coefficient_combined_2
    # points_linear = model$df$meanvar.alpha_coefficient_linear_1 + (sample_points) * model$df$meanvar.alpha_coefficient_linear_2
    
    p = ggplot(data=dat2) + geom_point(aes(x=1/l_mean, y=l_alpha, size=l_weight, l_chromosome=l_chromosome), alpha=0.5) +
      theme_minimal() + xlab("mean") + ylab("alpha") +
      scale_size_continuous() + 
      # scale_y_continuous(limits=c(0, max(model$l_alpha))) + 
      # scale_x_continuous(limits=c(0, max(1/dat2$l_mean))) + 
      geom_point(data=data.frame(x=sample_points, y=points_robust), aes(x=x, y=y), size=0.1, alpha=0.3, color="red")
      # geom_point(data=data.frame(x=sample_points, y=points_linear), aes(x=x, y=y), size=0.1, alpha=0.3, color="orange")
      # 
      # geom_abline(intercept = 0, slope = model$df$alpha*scalar, color="#E69F00") +
      # geom_abline(intercept = 0, slope = model$df$alpha_linear*scalar, color="#A3C1AD", lty=3)
    
    # print(paste0("Observed alpha: ", round(model$df$alpha, digits = 3), "(", round(model$df$alpha_rsquared_robust, digits = 3), "/",
                 # round(model$df$alpha_rsquared_linear, digits = 3), ")"))
    
  }else{
    dat = data.frame(l_mean=model$l_mean, l_mean_unscaled=model$l_mean_unscaled, l_var=model$l_var, l_weight=model$l_weight, l_chromosome=model$l_chromosome)
    sample_points = seq(0, max(dat$l_mean) + 3, 0.01)
    points_robust = sample_points * model$df$meanvar.beta.coefficient_combined_1 + sample_points^2 * model$df$meanvar.beta.coefficient_combined_2 + 0
    # points_linear = sample_points * model$df$meanvar.beta_coefficient_linear_1 + sample_points^2 * model$df$meanvar.beta_coefficient_linear_2 + 0
    p = ggplot(data=dat) + geom_point(aes(x=l_mean, y=l_var, size=l_weight, chromosome=l_chromosome), alpha=0.5) +
      theme_minimal() + xlab("mean") + ylab("variance") +
      scale_size_continuous() + scale_color_continuous() + scale_shape_manual(values=c(18,19)) +
      geom_point(data=data.frame(x=sample_points, y=points_robust), aes(x=x, y=y), size=0.1, alpha=0.3, color="red")
      # geom_point(data=data.frame(x=sample_points, y=points_linear), aes(x=x, y=y), size=0.1, alpha=0.3, color="orange")
      # geom_abline(intercept = 0, slope = model$df$beta, color="#E69F00") +
      # geom_abline(intercept = 0, slope = model$df$beta_linear, color="#A3C1AD", lty=3)
  }

  # print(paste0("Observed beta: ", round(model$df$beta, digits = 3), "(", round(model$df$beta_rsquared_robust, digits = 3), "/",
               # round(model$df$beta_rsquared_linear, digits = 3), ")"))
  return(p)
}

#' plotDataset
#'
#' \code{plotDataset} wrapper to plot a whole dataset using plottting function f to a single file
#'
#' @param object QDNAseq object for a single sample
#' @param file file to save plot to
#' @param chunkSize split data set into chunks in order to keep size ok, a corresponding number of files will be
#' generated
#' @param f function to apply to all samples
#'
#' @return A QDNAseq plot
#'
#' @export
plotDataset <- function(object, file=NULL,chunkSize=500, f=plotQDNAseq,...){
  stopifnot(grepl("\\.pdf$", file))
  stopifnot(class(object)[[1]] == "QDNAseqCopyNumbers")
  require(QDNAseq, warn.conflicts = FALSE)
  require(Biobase, warn.conflicts = FALSE)


  if(dim(object)[[2]] > chunkSize){
    m = ceiling(dim(object)[[2]] / chunkSize)
    j = 1
    while(j <= m){
      new_file = paste0(sub('\\.pdf$', '', file), "-", j, ".pdf")
      pdf(new_file, width = 9, height = 4)
      start = 1+(j-1)*chunkSize
      end = min(j * chunkSize, dim(object)[[2]])
      for (i in start:end){
        tmp = object[, i]
        a = f(tmp, ...)
        if(all(c("layers", "theme") %in% names(a))){
          print(f(tmp, ...))
        }else{
          f(tmp, ...)
        }
      }
      dev.off()
      j = j + 1
   }
  }else{
    pdf(file, width = 9, height = 4)
    for (i in 1:dim(object)[[2]]){
      tmp = object[, i]
      a = f(tmp, ...)
      if(all(c("layers", "theme") %in% names(a))){
        print(f(tmp, ...))
      }else{
        f(tmp, ...)
      }
    }
    dev.off()
  }
}
