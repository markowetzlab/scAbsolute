# This code is adapted from the original changepoint package.
# Original author: Rebecca Killick <r.killick@lancs.ac.uk>
# https://github.com/rkillick/changepoint/
#
# Changes and updates:
# Michael P Schneider
# 2020/2021

penalty_decision_custom = function(penalty, pen.value, n, diffparam, asymcheck){

  if((penalty=="SIC0") || (penalty=="BIC0")){
    pen.return=diffparam*log(n)
  }
  else if((penalty=="SIC") || (penalty=="BIC")){
    pen.return=(diffparam+1)*log(n)
  }
  else if(penalty=="MBIC"){
    pen.return=(diffparam+2)*log(n)
  }
  else if((penalty=="SIC1") || (penalty=="BIC1")){
    stop("SIC1 and BIC1 have been depreciated, use SIC or BIC for the same result.")
  }
  else if(penalty=="AIC0"){
    pen.return=2*diffparam
  }
  else if(penalty=="AIC"){
    pen.return=2*(diffparam+1)
  }
  else if(penalty=="AIC1"){
    stop("AIC1 has been depreciated, use AIC for the same result.")
  }
  else if(penalty=="Hannan-Quinn0"){
    pen.return=2*diffparam*log(log(n))
  }
  else if(penalty=="Hannan-Quinn"){
    pen.return=2*(diffparam+1)*log(log(n))
  }
  else if(penalty=="Hannan-Quinn1"){
    stop("Hannan-Quinn1 has been depreciated, use Hannan-Quinn for the same result.")
  }
  else if(penalty=="None"){
    pen.return=0
  }
  else if((penalty!="Manual")&&(penalty!="Asymptotic")){
    stop('Unknown Penalty')
  }
  if((penalty=="Manual")&&(is.numeric(pen.value)==FALSE)){
    pen.value=try(eval(parse(text=paste(pen.value))),silent=TRUE)
    if(class(pen.value)=='try-error'){
      stop('Your manual penalty cannot be evaluated')
    }else{
      pen.return=pen.value
    }
  }

  if((penalty=="Manual")&&(is.numeric(pen.value)==TRUE)){
    pen.return=pen.value
  }

  if(pen.return < 0){
    stop('pen.value cannot be negative, please change your penalty value')
  }else{
    return(pen.return)
  }
}

wrap_segmentation <- function(data, costfunc, pen.value=0, penalty="BIC", debug=FALSE, gap=NULL, alpha=NA, rpc=NA){

  stopifnot(is.null(dim(data))==TRUE)

  stopifnot(costfunc %in% c("meanvar.poisson", "meanvar.norm", "meanvar.negbin", "meanvar.negbin_alpha",
                            "meanvar.negbin_ml", "meanvar.negbin_mm", "meanvar.negbin_bcmm", "meanvar.negbin_cml",
                            "fpop.negbin", "fpop.negbin.gap", "fpop.normal.posthoc"))

  n = length(data)
  # if(n<(2*minseglen)){stop('Minimum segment length is too large to include a change in this data')}
  minseglen=2
  # costfunc = "meanvar.poisson"
  if(penalty=="MBIC"){
    MBIC = 1
  }else{
    MBIC = 0
  }
  #   stop("we do not support MBIC")
  #   costfunc = paste0(costfunc, ".mbic")
  # }
  # NOTE asymcheck functionality is disabled here
  if(costfunc %in% c("meanvar.poisson", "meanvar.norm")){
    diffparam = 1
  }else{
    diffparam = 2
  }
  pen.value = penalty_decision_custom(penalty, pen.value, n, diffparam=diffparam, asymcheck=NULL)
  # pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck=costfunc, method=mul.method)
  # print(pen.value)

  # NOTE the addition of the raw data in a third vector -> required for negbin computations
  sumstat=cbind(c(0,cumsum(data)),c(0,cumsum(data^2)),c(0, data))
  # original value of sumstat:
  # sumstat=cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))
  if(penalty == "MBIC" && costfunc %in% c("fpop.negbin", "fpop.negbin.gap", "fpop.normal.posthoc")){
    warning("FPOP doesn't support MBIC penalty! - changing to BIC")
    pen.value = penalty_decision_custom("BIC", pen.value, n, diffparam=diffparam, asymcheck=NULL)
  }

  if(debug){
    print(paste0("Penalty choice: ", penalty, " - ", pen.value))
  }

  if(costfunc == "fpop.negbin"){
    require(gfpop, quietly = TRUE, warn.conflicts = FALSE)
    myGraph <- graph(
      Edge(0, 0,"abs", penalty = pen.value),
      Edge(0, 0,"null"))
    fit = gfpop(data, myGraph, type="negbin")
    cpts = fit$changepoints
    return(cpts)
  }

  if(costfunc == "fpop.negbin.gap"){
    require(gfpop, quietly = TRUE, warn.conflicts = FALSE)
    stopifnot(!is.null(gap))
    myGraph <- graph(
      Edge(0, 0,"abs", penalty = pen.value, gap=gap),
      Edge(0, 0,"null"))
    fit = gfpop(data, myGraph, type="negbin")
    cpts = fit$changepoints
    return(cpts)
  }

  if(costfunc == "fpop.normal.posthoc"){
    require(ChangepointInference, quietly=TRUE, warn.conflicts = FALSE)
    dat = log2(data + 1 + .Machine$double.xmin)

    h = 100
    fit_inference <- changepoint_inference(dat, 'L0-fixed', pen.value, window_size = h, sig = NULL)

    cpts_all = fit_inference$change_pts
    pvalues = fit_inference$pvals

    cpts = cpts_all[p.adjust(pvalues, method="BH") < 0.05]
    cpts = cpts[cpts > 0]
    return(cpts)
  }

  # print(dim(sumstat))
  out=PELTE(sumstat, pen=pen.value, cost_func = costfunc, minseglen=minseglen, MBIC=MBIC, alpha=alpha, rpc=rpc)
  #cpts=out[[2]]
  return(out[[2]])
}


wrap_state_space_segmentation <- function(data, pen.value=0, penalty="NONE", method="PELT", n_cells=1, debug=FALSE, fullOutput=FALSE, costfunc = "state_space"){
  # NOTE: requires data shape (n_states, n_bins)
  # NOTE: data are logprobabilities -> cumsum with logsumexp

  stopifnot(dim(data)[[1]] > dim(data)[[2]])

  if(!is.null(dim(data))){
    n = dim(data)[[1]]
  }else{
    n = length(data)
  }

  if(penalty=="MBIC"){
    MBIC = 1
  }else{
    MBIC = 0
  }

  # NOTE asymcheck functionality is disabled here
  if(penalty == "NONE"){
    pen.value = 0
    MBIC = 0
  }else{
    pen.value = penalty_decision_custom(penalty, pen.value, n, diffparam=1, asymcheck=NULL)
  }

  # pen.value = penalty_decision(penalty, pen.value, n, diffparam=1, asymcheck=costfunc, method=mul.method)

  # NOTE we need to use logsumexp as data are logprobabilities!
  require(matrixStats, quietly = TRUE, warn.conflicts = FALSE)
  cumsum.logsumexp <- function(x){

    cums = -Inf
    cs = rep(NA, length(x))
    for(i in 1:length(x)){
      lx = c(cums, x[i])
      cums = logSumExp(lx)
      cs[i] = cums
    }
    return(cs)
  }

  sumstat = apply(data, 2, function(x) cumsum.logsumexp(x))

  # NOTE, bug if data is logprobs
  # sumstat = apply(data, 2, cumsum)

  # NOTE to add leading zeros to increase dimension of problem
  # so that output is of correct dimensionality
  sumstat = rbind(rep(0, dim(sumstat)[[2]]), sumstat)

  if(debug){
    print(paste0("Penalty choice: ", penalty, " - ", pen.value, " with ", n_cells, " cells (", pen.value/n_cells, ")"))
  }

  # costfunc = "state_space_penalized"
  minseglen = 2
  if(method == "PELT"){
    out=PELTE(sumstat,pen=pen.value/n_cells, cost_func = costfunc, minseglen=minseglen)
    #cpts=out[[2]]
    if(fullOutput){
      return(out)
    }
    return(out[[2]])
  }

  if(method == "CROPS"){
    pen.min=0 #penalty_decision_custom("AIC", 0.0, n, diffparam=1, asymcheck=NULL)
    stopifnot(pen.value > 0)
    pen.max=pen.value
    out = CROPSE(sumstat, cost=costfunc, min_pen=pen.min, max_pen=pen.max, minseglen=minseglen)
    if(fullOutput){
      return(out)
    }
    return(out[[2]])
  }

  stop("Method not supported!")

  out=PELTE(sumstat,pen=pen.value/n_cells, cost_func = costfunc, minseglen=minseglen)
  #cpts=out[[2]]
  return(out[[2]])
  # return(out)
}

PELTE = function(sumstat, pen=0, cost_func = "mean.norm", minseglen = 1, MBIC=0, alpha=0.0, rpc=0.0){
  n = length(sumstat[,1])-1
  m = length(sumstat[1,])
  tol = 0
  shape = 1
  if(is.na(alpha)){
    alpha = 0.0;
  }
  if(is.na(rpc)){
    rpc = 0.0;
  }
  
  # if(cost_func == "mean.norm" || cost_func == "var.norm" || cost_func == "meanvar.norm" || cost_func == "meanvar.exp" || cost_func == "meanvar.gamma" || cost_func == "meanvar.poisson" || cost_func == "meanvar.negbin_ml" || cost_func == "meanvar.negbin_mm"){
  #   MBIC = 0
  # }
  # MBIC = 0
  
  if(n<2){stop('Data must have at least 2 observations to fit a changepoint model.')}
  
  storage.mode(sumstat) = 'double'
  error=0
  
  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)
  
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
  answer=list()
  answer[[7]]=1
  on.exit(.C("FreePELT",answer[[7]]))
  
  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'
  
  min = 0
  optimal = 0
  max = 0
  
  answer=.C('PELT', cost_func = cost_func, sumstat = sumstat, n = as.integer(n), m = as.integer(m), pen = as.double(pen), cptsout = cptsout, error = as.integer(error), shape = as.double(shape), minorder = as.integer(min), optimalorder = as.integer(optimal), maxorder = as.integer(max), minseglen = as.integer(minseglen), tol = as.double(tol), lastchangelike = lastchangelike, lastchangecpts = lastchangecpts, numchangecpts = numchangecpts, MBIC = as.integer(MBIC), as.double(alpha), as.double(rpc))
  
  if(answer$error>0){
    stop("C code error:",answer$error,call.=F)
  }
  
  return(list(lastchangecpts=answer$lastchangecpts,cpts=sort(answer$cptsout[answer$cptsout>0]), lastchangelike=answer$lastchangelike, ncpts=answer$numchangecpts))
  
}

######## The whole thing could also possibly be extended to CROPS functionality. ############

wrap_CROPS_segmentation <- function(data, costfunc, pen.value=0){
  method = "PELT"
  penalty="CROPS"
  minseglen = 2
  ni = length(data)
# class=TRUE, param.est=TRUE, minseglen, func
  # mu <- mean(data)
  
  # NOTE the addition of the raw data in a third vector -> required for negbin computations
  sumstat=cbind(c(0,cumsum(data)),c(0,cumsum(data^2)),c(0, data))
  # sumstat=cbind(c(0,cumsum(data)),c(0,cumsum(data^2)),cumsum(c(0,(data-mu)^2)))
  # needs fix
  # costfunc="meanvar.negbin"
  # costfunc="meanvar.normal"
  if(costfunc %in% c("meanvar.poisson", "meanvar.norm")){
    diffparam = 1
  }else{
    diffparam = 2
  }
  pen.max.mbic = penalty_decision_custom("MBIC", pen.value, ni, diffparam=diffparam, asymcheck=NULL)
  pen.min.bic = penalty_decision_custom("BIC", pen.value, ni, diffparam=diffparam, asymcheck=NULL)
  pen.min.aic = penalty_decision_custom("AIC", pen.value, ni, diffparam=diffparam, asymcheck=NULL)
  pen.min = max(min(pen.min.aic, pen.min.bic), 0)
  dist = pen.max.mbic - pen.min
  pen.max = pen.max.mbic + 1.0
  print(paste0("Penalty range: ", pen.min, "(", pen.min.aic , "|", pen.min.bic, ") - ", pen.max, "(", pen.max.mbic, ")"))

  if(costfunc == "non_parametric"){
    require(changepoint.np)
    out = changepoint.np::cpt.np(data = data, pen.value=c(pen.min, pen.max), penalty = "CROPS", class=FALSE)
  }else{
    out = CROPSE(sumstat, cost=costfunc, min_pen=pen.min, max_pen=pen.max, minseglen=minseglen)
  }
  
  return(out)
  # if(class==TRUE){
      # ans = class_input(data=data,cpttype=cpttype, method="PELT", test.stat=test.stat, penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.est, out=out,shape=shape)
      # if(func=="var"){
        # param.est(ans)=c(param.est(ans),mean=mu)
      # }
    # return(ans)
  # }else{return(out)}
}

CROPSE <- function(sumstat, cost ,min_pen=log(length(sumstat)/3-1),max_pen=10*log(length(sumstat)/3-1), minseglen) {
  # THis function is taken from the R package
  # https://github.com/rkillick/changepoint
  PELT=TRUE
  shape = 1

  NCALC=0
  pen_interval <- c(min_pen,max_pen)
  # NOTE change this for reasons of compatibility
  # nj = length(sumstat)/3 - 1
  nj = dim(sumstat)[[1]] - 1

  test_penalties <- NULL
  numberofchangepoints <- NULL
  penal <- NULL
  overall_cost <- array()
  segmentations <- NULL
  b_between <- array()

  ##### Want to store and use Func, M and CP in PELT

  count <- 0

  while (length(pen_interval) > 0){

    new_numcpts <- array()
    new_penalty <- array()
    new_cpts <- array()

    for (b in 1:length(pen_interval)) {

      ans<- PELTE(sumstat,pen=pen_interval[b], cost_func = cost, minseglen = minseglen)
      resultingcpts <- ans[[2]]
      new_numcpts[b] <- length(resultingcpts)
      new_cpts[b] <- list(resultingcpts[-length(resultingcpts)])
      new_penalty[b] <- ans[[3]][nj+1]-(ans[[4]][nj+1]-1)*pen_interval[b]
    }

    if (count == 0){
      print(paste("Maximum number of runs of algorithm = ", new_numcpts[1] - new_numcpts[2] + 2, sep = ""))
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }

    else{
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }


    ## Add the values calculated to the already stored values
    test_penalties <- unique((sort(c(test_penalties,pen_interval))))
    new_numcpts <- c(numberofchangepoints,new_numcpts)
    new_penalty <- c(penal,new_penalty)

    new_cpts <- c(segmentations,new_cpts)
    numberofchangepoints <- -sort(-new_numcpts) ##can use sort to re-order
    penal <- sort(new_penalty)

    ls <- array()

    for (l in 1:length(new_cpts)){
      ls[l] <- length(new_cpts[[l]])
    }


    ls1 <- sort(ls,index.return = T, decreasing = T)
    ls1 <- ls1$ix


    segmentations <- new_cpts[c(ls1)]

    pen_interval <- NULL
    tmppen_interval <- NULL

    for (i in 1:(length(test_penalties)-1)){
      if(abs(numberofchangepoints[i]-numberofchangepoints[i+1])>1){ ##only need to add a beta if difference in cpts>1
        j <- i+1
        tmppen_interval <- (penal[j] - penal[i]) * (((numberofchangepoints[i]) - (numberofchangepoints[j]))^-1)
        pen_interval <- c(pen_interval, tmppen_interval )
      }
    }

    if(length(pen_interval)>0){
      for(k in length(pen_interval):1){
        index <- which.min(abs(pen_interval[k]-test_penalties))
        if (isTRUE(all.equal(pen_interval[k], test_penalties[index]))){
          pen_interval=pen_interval[-k]
        }
      }
    }
  }


  ##PRUNE VALUES WITH SAME num_cp
  for(j in length(test_penalties):2){
    if(numberofchangepoints[j]==numberofchangepoints[j-1]){
      numberofchangepoints=numberofchangepoints[-j]
      test_penalties=test_penalties[-j]
      penal=penal[-j]
      segmentations = segmentations[-j]
    }
  }

  ###calculate beta intervals
  nb=length(test_penalties)
  beta.int=rep(0,nb)
  beta.e=rep(0,nb)
  for(k in 1:nb){
    if(k==1){
      beta.int[1]=test_penalties[1]
    }else{
      beta.int[k]=beta.e[k-1]
    }
    if(k==nb){
      beta.e[k]=test_penalties[k]
    }else{
      beta.e[k]=(penal[k]-penal[k+1])/(numberofchangepoints[k+1]-numberofchangepoints[k])
    }

  }

  return(list(cpt.out = rbind(beta_interval = beta.int,numberofchangepoints,penalised_cost = penal),changepoints = segmentations))
  #segmentations is output matrix
  #beta.int
}
