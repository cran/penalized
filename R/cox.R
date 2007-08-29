.coxfit <- function(response) {

  n <- nrow(response)
  type <- attr(response, "type")
  if (!type %in% c("right", "counting"))
    stop("Cox model doesn't support \"", type, "\" survival data")
 
  if (ncol(response) == 2) {
    time <- response[,1]
    status <- response[,2]
    Riskset <- outer(time, time, "<=")
  } else {
    time <- response[,2]
    start <- response[,1]
    status <- response[,3]
    Riskset <- outer(time, time, "<=") & outer(time, start, ">")
  }    

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {

    if (!missing(leftout)) {
      status <- status[!leftout]
      time <- time[!leftout]
      Riskset <- Riskset[!leftout, !leftout]
    }
    ws <- as.vector(exp(lp))
    if (any(ws == Inf | ws == 0)) { 
      ws <- 1e-10 + 1e10 * status
      exploded <- TRUE
    } else {
      exploded <- FALSE
    }

    breslows <- as.vector(status / Riskset %*% ws)
    breslow <- as.vector(breslows[status==1] %*% Riskset[status==1,,drop=FALSE])
    
    # The martingale residuals
    residuals <- status - breslow * ws
    
    # The loglikelihood
    if (!exploded)
      loglik <- -sum(ws * breslow) + sum(log(breslows[status==1])) + sum(lp[status==1])
    else
      loglik <- NA

    # The weights matrix
    Pij <- outer(ws, breslows) * t(Riskset)
    W <- - crossprod(t(Pij[,status==1]))
    diag(W) <- diag(W) + breslow * ws
    
    # The fitted baseline
    dtimes <- time[status==1]
    basesurv <- cumprod(1-breslows[status==1][sort.list(dtimes)])
    if (max(dtimes) < max(time)) {
      basetimes <- c(0, sort(dtimes), max(time))
      basesurv <- c(1, basesurv, basesurv[length(basesurv)])
    } else {
      basetimes <- c(0, sort(dtimes))
      basesurv <- c(1, basesurv)
    }
    
    baseline <- new("breslow")
    baseline@time <- basetimes
    baseline@curves <- matrix(basesurv,1,byrow=TRUE)

    return(list(residuals = residuals, loglik = loglik, W = W, lp = lp, fitted = exp(lp), nuisance = list(baseline = baseline)))
  }
  
  #cross-validated likelihood
  cvl <- function(lp, leftout)
  { 
    ws <- exp(lp)
    somw <- apply(Riskset, 1, function(rr) sum(ws[rr]))
    cvls <- numeric(length(leftout))
    for (k in which(leftout)) {
      pij <- ws[k] / somw
      cvls[k] <- sum(log(1 - pij[(status ==1) & !Riskset[k,]])) + status[k] * log(pij[k])
    }
    return(sum(cvls[leftout]))
  }
 
  # mapping from the linear predictor lp to an actual prediction
  prediction <- function(lp, nuisance) {
    out <- nuisance$baseline
    out@curves <- nuisance$baseline@curves ^ exp(lp)
    out
  }
 
  return(list(fit = fit, cvl = cvl, prediction = prediction))
}



.coxgamma <- function(response, unpenalized, data) {

  if (is.matrix(unpenalized)) {
    if (ncol(unpenalized) > 0) 
      startgamma <- coefficients(coxph(response ~ ., data = as.data.frame(unpenalized), method = "breslow"))
    else
      startgamma <- numeric(0)
  } else {
    .response <- response
    form <- as.formula(paste(".response~", paste(c("1", attr(terms(unpenalized),"term.labels")), collapse="+"))) 
    startgamma <- coefficients(coxph(form, data = data, method = "breslow"))
  }
}

                                        
# merges predicted survival curves with different time points
# input: a list of breslow objects
.coxmerge <- function(predictions) {

  times <- sort(unique(unlist(lapply(predictions, time))))
  curves <- sapply(predictions, function(pred) {
    res <- rep(NA, length(pred@time))
    res[times %in% time(pred)] <- pred@curves[1,]
    res <- sapply(1:length(res), function(i) {
      if (is.na(res[i]) && any(!is.na(res[-(1:i)])))
        res[min(which(!is.na(res) & (1:length(res)>i)))]
      else
        res[i]
    })
  })
  out <- new("breslow")
  out@time <- times
  out@curves <-  t(curves)
  out
}





