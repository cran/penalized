####################################
# Fits the penalized regression model
####################################                                    
penalized <- function(response, penalized, unpenalized, lambda1=0, lambda2=0, data, model = c("cox", "logistic", "linear"),
  startbeta, startgamma, steps =1, epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data)
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data)[[attr(terms(response), "response")]]
  }

  
  # determine the model if missing
  if (missing(model)) {
    if (is(response, "Surv")) model <- "cox"
    else if (all(response %in% 0:1)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  }
  model <- match.arg(model)
  
  # make sure the response has he correctformat and determine the sample size
  if (model == "cox" && !is(response, "Surv")) response <- Surv(response)
  n <- if (model == "cox") length(response)/2 else length(response)

  # fill in miscelaneous missing input
  if (missing(unpenalized)) unpenalized <- matrix(,n,0)
  if (missing(data)) data <-  as.data.frame(matrix(,n,0))
  if (missing(maxiter)) maxiter <- if (all(lambda1 == 0)) 25 else Inf

  # turn penalized into a matrix if it is not a matrix
  if (is.data.frame(penalized) || is.vector(penalized)) penalized <- as.matrix(penalized)
  if (is(penalized, "formula")) {
    oldcontrasts <- options("contrasts")[[1]]
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- model.matrix(penalized, data)
    options(contrasts = oldcontrasts)
    penalized <- penalized[,colnames(penalized) != "(Intercept)", drop = FALSE]
  }
  m <- ncol(penalized)
  if (missing(startbeta)) 
    startbeta <- rep(0, m)

  if (nrow(penalized) != n) {
    stop("The length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  
  # startgamma is set at the model of infinite penalty
  if (missing(startgamma)) {
    startgamma <- switch(model,
      cox = .coxgamma(response, unpenalized, data),
      logistic = .logitgamma(response, unpenalized, data),
      linear = .lmgamma(response, unpenalized, data))
  }

  intercept <- (model != "cox")
  if (is(unpenalized, "formula")) intercept <- intercept && (attr(terms(unpenalized), "intercept") == 1)

  # Join unpenalized and penalized covariates and standardize to the same scale
  inputlambda1 <- lambda1
  inputlambda2 <- lambda2
  prepared <- .prepare(penalized, unpenalized, lambda1, lambda2, data, startbeta, 
    startgamma, intercept = intercept, standardize = standardize)
  X <- prepared$X 
  lambda1 <- prepared$lambda1 
  lambda2 <- prepared$lambda2
  weights <- prepared$weights
  beta <- prepared$beta
  orthogonalizer <- prepared$orthogonalizer
  rm(prepared)
  p <- ncol(X)

  # Prepare the model
  fit <- switch(model,
    cox = .coxfit(response)$fit,
    logistic = .logitfit(response)$fit,
    linear = .lmfit(response)$fit
  )

  # If a steps argument is given, determine where to start
  if (steps > 1) {
    nullgamma <- switch(model, 
      cox = .coxgamma(response, unpenalized, data),
      logistic = .logitgamma(response, unpenalized, data),
      linear = .lmgamma(response, unpenalized, data))
    null <- c(nullgamma, numeric(m)) * weights
    nzn <- (null != 0)
    lp <- X[,nzn,drop=FALSE] %*% null[nzn]
    gradient <- as.vector(crossprod(X[,!nzn,drop=FALSE], fit(lp)$residuals))
    rel <- gradient / lambda1[!nzn]
    from <- max(abs(rel))
  } else {
    from <- 1
  }
  lambda1s <- as.list(seq(from, 1, length.out=steps+1))
  if (steps == 1) lambda1s <- lambda1s[-1]

  # fit the model for all lambdas
  outs <- lapply(lambda1s, function(rellambda) {
    if (!all(lambda1 == 0)) {
      if (all(lambda2 == 0)) {
        out <- .steplasso(beta = beta, lambda = rellambda * lambda1, lambda2 = 0, X = X,
          fit = fit, trace = trace, epsilon = epsilon, maxiter = maxiter)
      } else {
        out <- .lasso(beta = beta, lambda = rellambda * lambda1, lambda2 = lambda2, X = X,
          fit = fit, trace = trace, epsilon = epsilon, maxiter = maxiter)
      }
    } else {
      if (m > n) {
        P <- .makeP(X, lambda2)
        gams <- solve(crossprod(t(P)), P %*% beta)
        PX <- P %*% t(X)
        Pl <- P * matrix(sqrt(lambda2), nrow(P), ncol(P), byrow = TRUE)
        PlP <- crossprod(t(Pl))
        # Do ridge regression on gamma
        out <- .ridge(beta = gams, Lambda = PlP, X = t(PX), fit = fit, trace = trace,
          epsilon = epsilon, maxiter = maxiter)
        # and transform back
        out$beta <- as.vector(crossprod(P, out$beta)) 
      } else {
        out <- .ridge(beta = beta, Lambda = lambda2, X = X, fit = fit, trace = trace,
          epsilon = epsilon, maxiter = maxiter)
      }
    }
    if (trace) cat("\n")
    beta <<- out$beta

    out$beta <- out$beta / weights
    out
  })
  
  # put the output in a penfit object
  outs <- sapply(1:length(outs), function(nr) {
    thislambda1 <- inputlambda1 * lambda1s[[nr]]
    .makepenfit(outs[[nr]], length(startgamma), model, thislambda1, inputlambda2, orthogonalizer)
  })

  if(length(outs)==1) 
    outs <- outs[[1]]
  else
    names(outs) <- paste("lambda1 =", unlist(lambda1s) * inputlambda1)

  outs
}

######################################
# Finds the cross-validated loglikelihood for a given penalty
######################################
cvl <- function(response, penalized, unpenalized, lambda1 = 0, lambda2= 0, data, 
  model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, groups, 
  epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data)
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data)[[attr(terms(response), "response")]]
  }

  # determine the model if missing
  if (missing(model)) {
    if (is(response, "Surv")) model <- "cox"
    else if (all(response %in% 0:1)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  }
  model <- match.arg(model)
  
  # make sure the response has he correctformat and determine the sample size
  if (model == "cox" && !is(response, "Surv")) response <- Surv(response)
  n <- if (model == "cox") length(response)/2 else length(response)

  # fill in miscelaneous missing input
  if (missing(unpenalized)) unpenalized <- matrix(,n,0)
  if (missing(data)) data <-  as.data.frame(matrix(,n,0))
  if (missing(maxiter)) maxiter <- if (all(lambda1 == 0)) 25 else Inf

  # turn penalized into a matrix if it is not a matrix
  if (is.data.frame(penalized) || is.vector(penalized)) penalized <- as.matrix(penalized)
  if (is(penalized, "formula")) {
    oldcontrasts <- options("contrasts")[[1]]
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- model.matrix(penalized, data)
    options(contrasts = oldcontrasts)
    penalized <- penalized[,colnames(penalized) != "(Intercept)", drop = FALSE]
  }
  m <- ncol(penalized)
  if (missing(startbeta)) 
    startbeta <- rep(0, m)

  if (nrow(penalized) != n) {
    stop("The length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  
  # startgamma is set at the model of infinite penalty
  if (missing(startgamma)) {
    startgamma <- switch(model,
      cox = .coxgamma(response, unpenalized, data),
      logistic = .logitgamma(response, unpenalized, data),
      linear = .lmgamma(response, unpenalized, data))
  }
  
  intercept <- (model != "cox")
  if (is(unpenalized, "formula")) intercept <- intercept && (attr(terms(unpenalized), "intercept") == 1)

  # Join unpenalized and penalized covariates and standardize to the same scale
  inputlambda1 <- lambda1
  inputlambda2 <- lambda2
  prepared <- .prepare(penalized, unpenalized, lambda1, lambda2, data, startbeta, 
    startgamma, intercept = intercept, standardize = standardize)
  X <- prepared$X 
  lambda1 <- prepared$lambda1 
  lambda2 <- prepared$lambda2
  weights <- prepared$weights
  beta <- prepared$beta
  orthogonalizer <- prepared$orthogonalizer
  rm(prepared)
  p <- ncol(X)

  # Prepare the model
  fit <- switch(model,
    cox = .coxfit(response),
    logistic = .logitfit(response),
    linear = .lmfit(response)
  )
  
  if (missing(groups)) {
    if (missing(fold) || fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  }
  
  res <- .cvl(X, lambda1, lambda2, beta, fit=fit$fit, cvl=fit$cvl, groups=groups, epsilon=epsilon, maxiter=maxiter, trace = trace)
  out <- res["cvl"]
  return(list(cvl = res$cvl, fullfit = .makepenfit(res$fit, length(startgamma), model, inputlambda1, inputlambda2, orthogonalizer)))
}

######################################
# Finds the optimal cross-validated L1-penalty for a given L2-penalty
######################################
optL1 <- function(response, penalized, unpenalized, lambda2 = 0, data, 
  model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, groups, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE, 
  accuracy = 1e-3) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data)
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data)[[attr(terms(response), "response")]]
  }

  # determine the model if missing
  if (missing(model)) {
    if (is(response, "Surv")) model <- "cox"
    else if (all(response %in% 0:1)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  }
  model <- match.arg(model)
  
  # make sure the response has he correctformat and determine the sample size
  if (model == "cox" && !is(response, "Surv")) response <- Surv(response)
  n <- if (model == "cox") length(response)/2 else length(response)

  # fill in miscelaneous missing input
  if (missing(unpenalized)) unpenalized <- matrix(,n,0)
  if (missing(data)) data <-  as.data.frame(matrix(,n,0))

  # turn penalized into a matrix if it is not a matrix
  if (is.data.frame(penalized) || is.vector(penalized)) penalized <- as.matrix(penalized)
  if (is(penalized, "formula")) {
    oldcontrasts <- options("contrasts")[[1]]
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- model.matrix(penalized, data)
    options(contrasts = oldcontrasts)
    penalized <- penalized[,colnames(penalized) != "(Intercept)", drop = FALSE]
  }
  m <- ncol(penalized)
  if (missing(startbeta)) 
    startbeta <- rep(0, m)

  if (nrow(penalized) != n) {
    stop("The length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  
  # startgamma is set at the model of infinite penalty
  if (missing(startgamma)) {
    startgamma <- switch(model,
      cox = .coxgamma(response, unpenalized, data),
      logistic = .logitgamma(response, unpenalized, data),
      linear = .lmgamma(response, unpenalized, data))
  }
  
  intercept <- (model != "cox")
  if (is(unpenalized, "formula")) intercept <- intercept && (attr(terms(unpenalized), "intercept") == 1)

  # Join unpenalized and penalized covariates and standardize to the same scale
  inputlambda2 <- lambda2
  prepared <- .prepare(penalized, unpenalized, lambda1 = 1, lambda2, data, startbeta, 
    startgamma, intercept = intercept, standardize = standardize)
  X <- prepared$X 
  lambda1 <- prepared$lambda1 
  lambda2 <- prepared$lambda2
  weights <- prepared$weights
  beta <- prepared$beta
  orthogonalizer <- prepared$orthogonalizer
  rm(prepared)
  p <- ncol(X)

  # Prepare the model
  fit <- switch(model,
    cox = .coxfit(response),
    logistic = .logitfit(response),
    linear = .lmfit(response)
  )
  
  if (missing(groups)) {
    if (missing(fold) || fold == n) {
      groups <- 1:n
      fold <- n
    } else {
      groups <- sample(n) %% fold + 1
    }
  }
  
  # Find the highest interesting value of lambda
  nullgamma <- switch(model, 
    cox = .coxgamma(response, unpenalized, data),
    logistic = .logitgamma(response, unpenalized, data),
    linear = .lmgamma(response, unpenalized, data)
  )
  null <- c(nullgamma, numeric(m)) * weights
  nzn <- (null != 0)
  lp <- X[,nzn,drop=FALSE] %*% null[nzn]
  gradient <- as.vector(crossprod(X[,!nzn,drop=FALSE], fit$fit(lp)$residuals))
  rel <- gradient / lambda1[!nzn]
  maxlambda <- max(abs(rel))

  # The function to be optimized
  betas <- NULL
  maxcvl <- -Inf
  bestfit <- NULL
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(X, rellambda*lambda1, lambda2, beta, fit=fit$fit, cvl=fit$cvl, groups=groups, 
      epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <<- out$fit$beta
    if (is.null(betas)) {
      between <- mean(abs(beta - out$fit$beta))
      within <- mean(abs(out$betas - matrix(out$fit$beta, p, fold)))
      if (between < within) {
        betas <<- out$betas
      } 
    } else {
      betas <<- out$betas
    }
    if (out$cvl > maxcvl) {
      maxcvl <<- out$cvl
      bestfit <<- out$fit
    }
    out$cvl
  }
  
  #optimize it
  opt <- opt.brent(thiscvl, c(0, 1.1*maxlambda), maximum = TRUE, tol = accuracy)

  return(list(lambda = opt$argmax, cvl = opt$max, fullfit = .makepenfit(bestfit, length(startgamma), model, opt$argmax, inputlambda2, orthogonalizer)))
}

######################################
# Finds the optimal cross-validated L2-penalty for a given L1-penalty
######################################
optL2 <- function(response, penalized, unpenalized, lambda1 = 0, startlambda2 = 1, data, 
  model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, groups, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE, 
  accuracy = 1e-3) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data)
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data)[[attr(terms(response), "response")]]
  }

  # determine the model if missing
  if (missing(model)) {
    if (is(response, "Surv")) model <- "cox"
    else if (all(response %in% 0:1)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  }
  model <- match.arg(model)
  
  # make sure the response has he correctformat and determine the sample size
  if (model == "cox" && !is(response, "Surv")) response <- Surv(response)
  n <- if (model == "cox") length(response)/2 else length(response)

  # fill in miscelaneous missing input
  if (missing(unpenalized)) unpenalized <- matrix(,n,0)
  if (missing(data)) data <-  as.data.frame(matrix(,n,0))
  if (missing(maxiter)) maxiter <- if (all(lambda1 == 0)) 25 else Inf

  # turn penalized into a matrix if it is not a matrix
  if (is.data.frame(penalized) || is.vector(penalized)) penalized <- as.matrix(penalized)
  if (is(penalized, "formula")) {
    oldcontrasts <- options("contrasts")[[1]]
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- model.matrix(penalized, data)
    options(contrasts = oldcontrasts)
    penalized <- penalized[,colnames(penalized) != "(Intercept)", drop = FALSE]
  }
  m <- ncol(penalized)
  if (missing(startbeta)) 
    startbeta <- rep(0, m)

  if (nrow(penalized) != n) {
    stop("The length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  
  # startgamma is set at the model of infinite penalty
  if (missing(startgamma)) {
    startgamma <- switch(model,
      cox = .coxgamma(response, unpenalized, data),
      logistic = .logitgamma(response, unpenalized, data),
      linear = .lmgamma(response, unpenalized, data))
  }
  
  intercept <- (model != "cox")
  if (is(unpenalized, "formula")) intercept <- intercept && (attr(terms(unpenalized), "intercept") == 1)

  # Join unpenalized and penalized covariates and standardize to the same scale
  inputlambda1 <- lambda1
  prepared <- .prepare(penalized, unpenalized, lambda1, lambda2= 1, data, startbeta, 
    startgamma, intercept = intercept, standardize = standardize)
  X <- prepared$X 
  lambda1 <- prepared$lambda1 
  lambda2 <- prepared$lambda2
  weights <- prepared$weights
  beta <- prepared$beta
  orthogonalizer <- prepared$orthogonalizer
  rm(prepared)
  p <- ncol(X)

  # Prepare the model
  fit <- switch(model,
    cox = .coxfit(response),
    logistic = .logitfit(response),
    linear = .lmfit(response)
  )
  
  if (missing(groups)) {
    if (missing(fold) || fold == n) {
      groups <- 1:n
      fold <- n
    } else {
      groups <- sample(n) %% fold + 1
    }
  }

  # benchmark: cvl at infinite penalty
  g <- length(startgamma)
  if (g > 0) {
    nullfit <- .cvl(X[,1:g, drop=FALSE], lambda1 = rep(0,g), lambda2 = rep(0,g), startgamma, fit=fit$fit, cvl=fit$cvl, groups=groups, 
      epsilon=epsilon, maxiter=maxiter, trace = FALSE)
    nullcvl <- nullfit$cvl
    nullfit <- nullfit$fit
    nullgamma <- nullfit$fit$beta
  } else {
    nullcvl <- fit$cvl(numeric(n), !logical(n))
    nullfit <- list()
    nullgamma <- numeric(0)
    nullfit$fit <- fit$fit(numeric(n))
    nullfit$iterations <- 1
    nullfit$converged <- TRUE
  }
  
  # The function to be optimized
  betas <- NULL
  maxcvl <- nullcvl
  bestfit <- nullfit
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(X, lambda1, rellambda*lambda2, beta, fit=fit$fit, cvl=fit$cvl, groups=groups, 
      epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas)
    if (trace) cat("cvl=", out$cvl, "\n")
    if (out$cvl > - Inf) {
      beta <<- out$fit$beta
      if (is.null(betas)) {
        between <- mean(abs(beta - out$fit$beta))
        within <- mean(abs(out$betas - matrix(out$fit$beta, p, fold)))
        if (between < within) {
          betas <<- out$betas
        } 
      } else {
        betas <<- out$betas
      }
      if (out$cvl >= maxcvl) {
        maxcvl <<- out$cvl
        bestfit <<- out$fit
      }
    }
    out$cvl
  }
  
  # phase 1: find the order of magnitude of lambda
  left <- startlambda2
  leftcvl <- thiscvl(left)
  right <- 10*left
  rightcvl <- thiscvl(right)
  if (leftcvl < rightcvl || rightcvl == -Inf) {
    high <- right; highcvl <- rightcvl; low <- left; lowcvl <- leftcvl; fac <- 10
  } else {
    high <- left; highcvl <- leftcvl; low <- right; lowcvl <- rightcvl; fac <- 0.1
  }
  ready <- FALSE
  infmax <- ((abs(lowcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon) && (abs(highcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon))
  infmin <- FALSE
  while (!ready && !infmax) {
    nxt <- high*fac
    nxtcvl <- thiscvl(nxt)
    ready <- nxtcvl < highcvl
    if (!ready) {
      high <- nxt; highcvl <- nxtcvl; low <- high; lowcvl <- highcvl
    }
    infmax <- ((abs(lowcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon) && (abs(highcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon))
    infmin <- (fac < 1) && (abs(lowcvl - nxtcvl) / abs(nxtcvl + 0.1) < epsilon) 
  }
  
  # phase 2: optimize lambda within the order of magnitude found
  if (!infmax && ! infmin) {
    opt <- opt.brent(thiscvl, sort(c(low,nxt)), maximum = TRUE, tol = accuracy)
  } else {
    if (infmax) {
      opt <- list(argmax = Inf, max = nullcvl)
      bestfit$beta <- c(nullgamma, numeric(m))
      names(bestfit$beta) <- colnames(X)
      bestfit$penalty <- c(L1 = 0, L2 = Inf)
    } else {
      maxcvl <- -Inf    # clears bestfit
      bestcvl <- thiscvl(0)
      opt <- list(argmax = 0, max = bestcvl) 
    } 
  }

  return(list(lambda = opt$argmax, cvl = opt$max, fullfit = .makepenfit(bestfit, length(startgamma), model, inputlambda1, opt$argmax, orthogonalizer)))
}


