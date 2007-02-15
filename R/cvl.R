######################################
# Finds the cross-validated loglikelihood for a given penalty
######################################
cvl <- function(response, penalized, unpenalized, lambda1 = 0, lambda2= 0, positive = FALSE, 
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,  
  epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE) {

  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0 && !positive) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # check for the presence of penalty parameters
  if (ncol(prep$X) >= nrow(prep$X) && lambda1 == 0 && lambda2 == 0)
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response),
    logistic = .logitfit(prep$response),
    linear = .lmfit(prep$response)
  )

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  
  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  res <- .cvl(prep$X, lambda1 * prep$baselambda1, lambda2 * prep$baselambda2, 
    positive = prep$positive, beta = prep$beta, fit=fit$fit, cvl=fit$cvl, 
    prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, 
    trace = trace, quit.if.failed = FALSE)
  res$predictions <- switch(prep$model, 
    cox = .coxmerge(res$predictions),
    logistic = .logitmerge(res$predictions),
    linear = .lmmerge(res$predictions)
  )

  return(list(
    cvl = res$cvl, 
    predictions = res$predictions, 
    fold = groups, 
    fullfit = .makepenfit(res$fit, pu, prep$model, lambda1, 
      lambda2, prep$orthogonalizer, prep$weights)
  ))
}


######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0, 
  positive = FALSE, data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE,
  steps = 100, minsteps = steps/4, log = FALSE) {

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response),
    logistic = .logitfit(prep$response),
    linear = .lmfit(prep$response)
  )

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)

  # find the maxlambda1 and minlambda1
  if (missing(maxlambda1)) {
    if (pu > 0) 
      lp <- drop(prep$X[,1:pu,drop=FALSE] %*% prep$nullgamma)
    else 
      lp <- numeric(n)
    gradient <- drop(crossprod(prep$X[,pu+1:pp,drop=FALSE], fit$fit(lp)$residuals))
    rel <- gradient / prep$baselambda1[pu+1:pp]
    maxlambda1 <- max(ifelse(prep$positive[pu+1:pp],  rel, abs(rel)))
  }
  if (missing(minlambda1)) {
    if (log) 
      stop("argument \"minlambda1\" is missing. please specify \"minlambda1\" or set log = FALSE", call. = FALSE)
    else
      minlambda1 <- 0
  }  
  
  # find the sequence from maxlambda1 to minlambda1
  if (steps < 2) stop("please set \"steps\" >= 2", call. = FALSE)
  if (log) {
    lambda1s <- exp(seq(log(maxlambda1), log(minlambda1), length.out = steps))
  } else {
    lambda1s <- seq(maxlambda1, minlambda1, length.out = steps)
  }

  # benchmark: cvl at infinite penalty
  if (pu > 0) {
    nullfit <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), 
      lambda2 = rep(0,pu), positive = FALSE, beta = prep$nullgamma, fit=fit$fit, 
      cvl=fit$cvl, prediction = fit$prediction, groups=groups, epsilon=epsilon, 
      maxiter=maxiter, trace = FALSE)
    nullcvl <- nullfit$cvl
    nullfit <- nullfit$fit
  } else {
    nullcvl <- fit$cvl(numeric(n), !logical(n))
    nullfit <- list()
    nullfit$fit <- fit$fit(numeric(n))
    nullfit$iterations <- 1
    nullfit$converged <- TRUE
  }

  # the actual repeated cvl-calculation
  betas <- NULL
  beta <- prep$beta
  cvls <- rep(NA,length(lambda1s))
  finished <- FALSE
  iter <- 0
  fits <- vector("list", length = length(lambda1s))
  predictions <- vector("list", length = length(lambda1s))
  while (!finished) {
    iter <- iter + 1
    rellambda <- lambda1s[iter]
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, rellambda*prep$baselambda1, lambda2*prep$baselambda2, 
      positive = prep$positive, beta = beta, fit=fit$fit, cvl=fit$cvl, 
      prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, 
      trace = trace, betas = betas, quit.if.failed=FALSE)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && ((cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) && (iter >= minsteps))) || (iter == length(lambda1s))
  }

  # remove the tail of the output
  if (fold > 1) {
    lambda1s <- lambda1s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]                                           
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }

  # merge the cross-validated predictions
  predictions <- lapply(predictions, function(preds) {
    switch(prep$model, 
      cox = .coxmerge(preds),
      logistic = .logitmerge(preds),
      linear = .lmmerge(preds)
    )
  })
  
  # create all the penfit objects
  makethisfit <- function(iter)
    .makepenfit(fits[[iter]], pu, prep$model, lambda1s[[iter]], lambda2, 
    prep$orthogonalizer, prep$weights)

  return(list(
    lambda = lambda1s, 
    fold = groups, 
    cvl = cvls, 
    predictions = predictions, 
    fullfit = lapply(1:length(cvls), makethisfit)
  ))
} 

######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, 
  positive = FALSE, data, model = c("cox", "logistic", "linear"), startbeta, startgamma, 
  fold, epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
  steps = 100, minsteps = steps/4, log = TRUE) {

  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0 && !positive) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response),
    logistic = .logitfit(prep$response),
    linear = .lmfit(prep$response)
  )

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)

  # Find the sequence from maxlambda2 to minlambda2
  if (!log && missing(minlambda2)) minlambda2 <- 0
  if (steps < 2) stop("please set \"steps\" >= 2", call. = FALSE)
  if (log) 
    lambda2s <- exp(seq(log(maxlambda2), log(minlambda2), length.out = steps))
  else
    lambda2s <- seq(maxlambda2, minlambda2, length.out = steps)

  # benchmark: cvl at infinite penalty
  if (pu > 0) {
    nullfit <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), lambda2 = rep(0,pu),
      positive = FALSE, beta = prep$nullgamma, fit=fit$fit, cvl=fit$cvl, 
      prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, trace = FALSE)
    nullcvl <- nullfit$cvl
    nullfit <- nullfit$fit
  } else {
    nullcvl <- fit$cvl(numeric(n), !logical(n))
    nullfit <- list()
    nullfit$fit <- fit$fit(numeric(n))
    nullfit$iterations <- 1
    nullfit$converged <- TRUE
  }

  # the actual repeated cvl-calculation
  betas <- NULL
  beta <- prep$beta
  cvls <- rep(NA,length(lambda2s))
  finished <- FALSE
  iter <- 0
  fits <- vector("list", length = length(lambda2s))
  predictions <- vector("list", length = length(lambda2s))
  while (!finished) {
    iter <- iter + 1
    rellambda <- lambda2s[iter]
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, lambda1*prep$baselambda1, rellambda*prep$baselambda2, 
      positive = prep$positive, beta = beta, fit=fit$fit, cvl=fit$cvl, 
      prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, 
      trace = trace, betas = betas, quit.if.failed=FALSE)
    if (trace) if (fold > 1) cat("cvl=", out$cvl, "\n") else cat("\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && (cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) &&
      (iter >= minsteps)) || (iter == length(lambda2s))
  }

  # remove the tail of the output
  if (fold > 1) {
    lambda2s <- lambda2s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }
         
  # merge the cross-validated predictions               
  predictions <- lapply(predictions, function(preds) {
    switch(prep$model, 
      cox = .coxmerge(preds),
      logistic = .logitmerge(preds),
      linear = .lmmerge(preds)
    )
  })
    
  # create the penfit objects                          
  makethisfit <- function(iter)
    .makepenfit(fits[[iter]], pu, prep$model, lambda1, lambda2s[[iter]], 
      prep$orthogonalizer, prep$weights)

  return(list(
    lambda = lambda2s, 
    fold = groups, 
    cvl = cvls, 
    predictions = predictions, 
    fullfit = lapply(1:length(cvls), makethisfit)
  ))
} 

######################################
# Finds the optimal cross-validated L1-penalty for a given L2-penalty
######################################
optL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0, 
  positive = FALSE, data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE, 
  tol = .Machine$double.eps^0.25) {

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response),
    logistic = .logitfit(prep$response),
    linear = .lmfit(prep$response)
  )

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  # Find the maxlambda1 and minlambda1
  if (missing(maxlambda1)) {
    if (pu > 0) 
      lp <- drop(prep$X[,1:pu,drop=FALSE] %*% prep$nullgamma)
    else 
      lp <- numeric(n)
    gradient <- drop(crossprod(prep$X[,pu+1:pp,drop=FALSE], fit$fit(lp)$residuals))
    rel <- gradient / prep$baselambda1[pu+1:pp]
    maxlambda1 <- max(ifelse(prep$positive[pu+1:pp],  rel, abs(rel)))
  }
  if (missing(minlambda1)) minlambda1 <- 0

  # The function to be optimized
  # Note: in passing it keeps track of the fit with the best cvl so far
  # It does this to be able to return more than just the optimal cvl-value
  betas <- NULL
  beta <- prep$beta
  maxcvl <- -Inf
  best <- NULL
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, rellambda*prep$baselambda1, lambda2 * prep$baselambda2, 
      positive = prep$positive, beta = beta, fit=fit$fit, cvl=fit$cvl, 
      prediction = fit$prediction, groups=groups, 
      epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <<- out$fit$beta
    if (is.null(betas)) {
      between <- mean(abs(beta - out$fit$beta))
      within <- mean(abs(out$betas - matrix(out$fit$beta, pp+pu, fold)))
      if (between < within) {
        betas <<- out$betas
      } 
    } else {
      betas <<- out$betas
    }
    if (out$cvl > maxcvl) {
      maxcvl <<- out$cvl
      best <<- out
    }
    out$cvl
  }
  
  #optimize it
  opt <- opt.brent(thiscvl, c(minlambda1, maxlambda1), maximum = TRUE, tol = tol)
  
  # merge the cross-validated predictions of the optimal model
  best$predictions <- switch(prep$model, 
    cox = .coxmerge(best$predictions),
    logistic = .logitmerge(best$predictions),
    linear = .lmmerge(best$predictions)
  )

  return(list(
    lambda = opt$argmax, 
    cvl = opt$max, 
    predictions = best$predictions, 
    fold = groups, 
    fullfit = .makepenfit(best$fit, pu, prep$model, opt$argmax, lambda2, 
      prep$orthogonalizer, prep$weights)
  ))
}

######################################
# Finds the optimal cross-validated L2-penalty for a given L1-penalty
######################################
optL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, 
  positive = FALSE, data, model = c("cox", "logistic", "linear"), startbeta, startgamma, 
  fold, epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE, 
  tol = .Machine$double.eps^0.25) {

  # maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0 && !positive) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response),
    logistic = .logitfit(prep$response),
    linear = .lmfit(prep$response)
  )

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  
  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)

  # benchmark: cvl at infinite penalty
  if (pu > 0) {
    null <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), lambda2 = rep(0,pu),
      positive = FALSE, beta = prep$nullgamma, fit=fit$fit, cvl=fit$cvl, 
      prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, 
      trace = FALSE)
    null$fit$beta <- c(null$fit$beta, numeric(pp))
  } else {
    null <- list()
    null.lp <- numeric(n)
    names(null.lp) <- rownames(prep$X)
    null$cvl <- fit$cvl(numeric(n), !logical(n))
    null$fit <- list()
    null$fit$beta <- c(numeric(pu+pp))
    names(null$fit$beta) <- colnames(prep$X)
    null$fit$fit <- fit$fit(null.lp)
    null$predictions <- lapply(as.list(null.lp), fit$prediction, 
      nuisance= null$fit$fit$nuisance)
    null$fit$iterations <- 1
    null$fit$converged <- TRUE
  }
  
  # The function to be optimized
  # Note: in passing it keeps track of the fit with the best cvl so far
  # It does this to be able to return more than just the optimal cvl-value
  betas <- NULL
  beta <- prep$beta
  best <- null
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, lambda1 * prep$baselambda1, rellambda*prep$baselambda2, 
      positive = prep$positive, beta = beta, fit=fit$fit, cvl=fit$cvl, groups=groups, 
      prediction = fit$prediction, epsilon=epsilon, maxiter=maxiter, 
      trace = trace, betas = betas)
    if (trace) cat("cvl=", out$cvl, "\n")
    if (out$cvl > - Inf) {
      beta <<- out$fit$beta
      if (is.null(betas)) {
        between <- mean(abs(beta - out$fit$beta))
        within <- mean(abs(out$betas - matrix(out$fit$beta, pp+pu, fold)))
        if (between < within) {
          betas <<- out$betas
        } 
      } else {
        betas <<- out$betas
      }
      if (out$cvl >= best$cvl) {
        best <<- out
      }
    }
    out$cvl
  }
  
  # phase 1: find the order of magnitude of lambda if not given by the user
  if (missing(minlambda2) || missing(maxlambda2)) {
    if (missing(minlambda2) && missing(maxlambda2)) {
      left <- 1
    } else if (missing(maxlambda2)) {
      left <- minlambda2
    } else {
      left <- maxlambda2 / 10
    }
    leftcvl <- thiscvl(left)
    right <- 10*left
    rightcvl <- thiscvl(right)
    if (leftcvl < rightcvl || rightcvl == -Inf) {
      high <- right
      highcvl <- rightcvl
      low <- left 
      lowcvl <- leftcvl
      fac <- 10
    } else {
      high <- left
      highcvl <- leftcvl
      low <- right
      lowcvl <- rightcvl
      fac <- 0.1
    }
    ready <- FALSE
    # infmax: the maximum is (numerically) at infinite penalty
    # infmin: the maximum is (numerically) at zero penalty
    infmax <- ((abs(lowcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon) && 
      (abs(highcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon))
    infmin <- FALSE
    while (!ready && !infmax) {
      nxt <- high*fac
      nxtcvl <- thiscvl(nxt)
      ready <- nxtcvl < highcvl
      if (!ready) {
        low <- high
        lowcvl <- highcvl
        high <- nxt
        highcvl <- nxtcvl
      }
      infmax <- ((abs(lowcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon) && 
        (abs(highcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon))
      infmin <- (fac < 1) && (abs(lowcvl - nxtcvl) / abs(nxtcvl + 0.1) < epsilon) 
    }
    minlambda2 <- min(low, nxt)
    maxlambda2 <- max(low, nxt)
  } else {
    infmax <- infmin <- FALSE
  }
  
  # phase 2: optimize lambda within the order of magnitude found
  if (!infmax && !infmin) {
    opt <- opt.brent(thiscvl, sort(c(minlambda2,maxlambda2)), maximum = TRUE, 
      tol = tol)
  } else {
    if (infmax) {
      opt <- list(argmax = Inf, max = null$cvl)
      names(best$fit$beta) <- colnames(prep$X)
      best$fit$penalty <- c(L1 = 0, L2 = Inf)
    } else {
      best$cvl <- -Inf    # clears bestfit
      best$cvl <- thiscvl(0)
      opt <- list(argmax = 0, max = best$cvl) 
    } 
  } 

  # merge the cross-validated predictions at the optimal fit
  best$predictions <- switch(prep$model, 
    cox = .coxmerge(best$predictions),
    logistic = .logitmerge(best$predictions),
    linear = .lmmerge(best$predictions)
  )

  return(list(
    lambda = opt$argmax, 
    cvl = opt$max, 
    predictions = best$predictions, 
    fold = groups, 
    fullfit = .makepenfit(best$fit, pu, prep$model, lambda1, opt$argmax, 
      prep$orthogonalizer, prep$weights)
  ))
}
      