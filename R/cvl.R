######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0, 
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE,
  steps = 100, minsteps = steps/4, log = FALSE) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data, globalenv())
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data))
      response <- eval(attr(terms(response), "variables"), globalenv())[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data, globalenv())[[attr(terms(response), "response")]]
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

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  names(groups) <- rownames(X)


  # Find the highest interesting value of lambda
  nullgamma <- switch(model,
    cox = .coxgamma(response, unpenalized, data),
    logistic = .logitgamma(response, unpenalized, data),
    linear = .lmgamma(response, unpenalized, data)
  )
  null <- c(nullgamma, numeric(m)) * weights
  nzn <- (null != 0)
  lp <- X[,nzn,drop=FALSE] %*% null[nzn]
  gradient <- drop(crossprod(X[,!nzn,drop=FALSE], fit$fit(lp)$residuals))  
  rel <- gradient / lambda1[!nzn]
  
  # which lambda-values?
  if (missing(maxlambda1)) maxlambda1 <- max(abs(rel))
  if (missing(minlambda1)) {
    if (log) 
      stop("argument \"minlambda1\" is missing. please specify \"minlambda1\" or set log = FALSE", call. = FALSE)
    else
      minlambda1 <- 0
  }  
  if (log) {
    lambda1s <- exp(seq(log(maxlambda1), log(minlambda1), length.out = steps+1))
  } else {
    lambda1s <- seq(maxlambda1, minlambda1, length.out = steps+1)
  }

  # benchmark: cvl at infinite penalty
  g <- length(nullgamma)
  if (g > 0) {
    nullfit <- .cvl(X[,1:g, drop=FALSE], lambda1 = rep(0,g), lambda2 = rep(0,g), nullgamma, fit=fit$fit, cvl=fit$cvl, 
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


  # The function to be optimized
  betas <- NULL
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
    out <- .cvl(X, rellambda*lambda1, lambda2, beta, fit=fit$fit, cvl=fit$cvl, prediction = fit$prediction, groups=groups,
      epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas, quit.if.failed=FALSE)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && ((cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) && (iter >= minsteps))) || (iter == length(lambda1s))
  }

  if (fold > 1) {
    lambda1s <- lambda1s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]                                           
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }

  predictions <- lapply(predictions, function(preds) {
    switch(model, 
      cox = .coxmerge(preds),
      logistic = .logitmerge(preds),
      linear = .lmmerge(preds)
    )
  })
  
  makethisfit <- function(iter)
    .makepenfit(fits[[iter]], length(startgamma), model, lambda1s[[iter]], inputlambda2, orthogonalizer, weights)

  return(list(lambda = lambda1s, fold = groups, cvl = cvls, predictions = predictions, fullfit = lapply(1:length(cvls), makethisfit)))
} 

######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, 
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
  steps = 100, minsteps = steps/4, log = TRUE) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data, globalenv())
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data))
      response <- eval(attr(terms(response), "variables"), globalenv())[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data, globalenv())[[attr(terms(response), "response")]]
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
  prepared <- .prepare(penalized, unpenalized, lambda1, lambda2 = 1, data, startbeta,
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

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  names(groups) <- rownames(X)

  # Find the highest interesting value of lambda
  nullgamma <- switch(model,
    cox = .coxgamma(response, unpenalized, data),
    logistic = .logitgamma(response, unpenalized, data),
    linear = .lmgamma(response, unpenalized, data)
  )
  g <- length(nullgamma)
  nullgamma <- nullgamma * weights[1:g]
  
  # which lambda-values?
  if (!log && missing(minlambda2)) {
      minlambda2 <- 0
  }  
  if (log) {
    lambda2s <- exp(seq(log(maxlambda2), log(minlambda2), length.out = steps+1))
  } else
    lambda2s <- seq(maxlambda2, minlambda2, length.out = steps+1)

  # benchmark: cvl at infinite penalty
  if (g > 0) {
    nullfit <- .cvl(X[,1:g, drop=FALSE], lambda1 = rep(0,g), lambda2 = rep(0,g), nullgamma, fit=fit$fit, cvl=fit$cvl, 
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

  # The function to be optimized
  betas <- NULL
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
    out <- .cvl(X, lambda1, rellambda*lambda2, beta, fit=fit$fit, cvl=fit$cvl, prediction = fit$prediction, groups=groups,
      epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas, quit.if.failed=FALSE)
    if (trace) if (fold > 1) cat("cvl=", out$cvl, "\n") else cat("\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && (cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) && (iter >= minsteps)) || (iter == length(lambda2s))
  }

  if (fold > 1) {
    lambda2s <- lambda2s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }
                        
  predictions <- lapply(predictions, function(preds) {
    switch(model, 
      cox = .coxmerge(preds),
      logistic = .logitmerge(preds),
      linear = .lmmerge(preds)
    )
  })
                              
  makethisfit <- function(iter)
    .makepenfit(fits[[iter]], length(startgamma), model, inputlambda1, lambda2s[[iter]], orthogonalizer, weights)

  return(list(lambda = lambda2s, fold = groups, cvl = cvls, predictions = predictions, fullfit = lapply(1:length(cvls), makethisfit)))
} 

######################################
# Finds the optimal cross-validated L1-penalty for a given L2-penalty
######################################
optL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0, 
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE, 
  tol = .Machine$double.eps^0.25) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data, globalenv())
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"), globalenv())[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data, globalenv())[[attr(terms(response), "response")]]
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
  
  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  names(groups) <- rownames(X)
  
  # Find the highest interesting value of lambda
  nullgamma <- switch(model, 
    cox = .coxgamma(response, unpenalized, data),
    logistic = .logitgamma(response, unpenalized, data),
    linear = .lmgamma(response, unpenalized, data)
  )
  null <- c(nullgamma, numeric(m)) * weights
  nzn <- (null != 0)
  lp <- X[,nzn,drop=FALSE] %*% null[nzn]
  gradient <- drop(crossprod(X[,!nzn,drop=FALSE], fit$fit(lp)$residuals))
  rel <- gradient / lambda1[!nzn]
  if (missing(maxlambda1)) maxlambda1 <- max(abs(rel))
  if (missing(minlambda1)) minlambda1 <- 0

  # The function to be optimized
  betas <- NULL
  maxcvl <- -Inf
  best <- NULL
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(X, rellambda*lambda1, lambda2, beta, fit=fit$fit, cvl=fit$cvl, prediction = fit$prediction, 
      groups=groups, epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas)
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
      best <<- out
    }
    out$cvl
  }
  
  #optimize it
  opt <- opt.brent(thiscvl, c(minlambda1, maxlambda1), maximum = TRUE, tol = tol)
  
  best$predictions <- switch(model, 
    cox = .coxmerge(best$predictions),
    logistic = .logitmerge(best$predictions),
    linear = .lmmerge(best$predictions)
  )


  return(list(lambda = opt$argmax, cvl = opt$max, predictions = best$predictions, fold = groups, 
    fullfit = .makepenfit(best$fit, length(startgamma), model, opt$argmax, inputlambda2, orthogonalizer, weights)))
}

######################################
# Finds the optimal cross-validated L2-penalty for a given L1-penalty
######################################
optL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, 
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
  epsilon = 1e-10, maxiter = Inf, standardize = FALSE, trace = TRUE, 
  tol = .Machine$double.eps^0.25) {

  # determine the response
  if (!missing(data)) response <- eval(as.list(match.call())$response, data, globalenv())
  if (is(response, "formula")) {
    if (!missing(penalized)) warning("Ignoring \"penalized\" argument because \"response\" is a formula.", call.=FALSE)
    penalized <- response
    if (missing(data)) 
      response <- eval(attr(terms(response), "variables"), globalenv())[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data, globalenv())[[attr(terms(response), "response")]]
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
  
  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  if(length(fold) == 1) {
    if (fold == n) {
      groups <- 1:n
    } else {
      groups <- sample(n) %% fold + 1
    }
  } else {
    if (length(fold) == n) {
      groups <- fold 
      fold <- max(groups)
    } else {
      stop("incorrect input of \"fold\"", call.=FALSE)
    }
  }
  if (!all(1:fold %in% groups)) stop("incorrect input of \"fold\"", call.=FALSE)
  names(groups) <- rownames(X)

  # benchmark: cvl at infinite penalty
  g <- length(startgamma)
  if (g > 0) {
    null <- .cvl(X[,1:g, drop=FALSE], lambda1 = rep(0,g), lambda2 = rep(0,g), beta[1:g], fit=fit$fit, cvl=fit$cvl,
      prediction = fit$prediction, groups=groups, epsilon=epsilon, maxiter=maxiter, trace = FALSE)
    nullgamma <- null$fit$fit$beta
  } else {
    null <- list()
    null.lp <- numeric(n)
    names(null.lp) <- rownames(X)
    null$cvl <- fit$cvl(numeric(n), !logical(n))
    null$fit <- list()
    nullgamma <- numeric(0)
    null$fit$beta <- c(numeric(m))
    names(null$fit$beta) <- colnames(X)
    null$fit$fit <- fit$fit(null.lp)
    null$predictions <- lapply(as.list(null.lp), fit$prediction, nuisance= null$fit$fit$nuisance)
    null$fit$iterations <- 1
    null$fit$converged <- TRUE
  }
  
  # The function to be optimized
  betas <- NULL
  best <- null
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(X, lambda1, rellambda*lambda2, beta, fit=fit$fit, cvl=fit$cvl, groups=groups, 
      prediction = fit$prediction, epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas)
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
      high <- right; highcvl <- rightcvl; low <- left; lowcvl <- leftcvl; fac <- 10
    } else {
      high <- left; highcvl <- leftcvl; low <- right; lowcvl <- rightcvl; fac <- 0.1
    }
    ready <- FALSE
    # infmax: the maximum is (numerically) at infinite penalty; infmin: the maximum is (numerically) at zero penalty
    infmax <- ((abs(lowcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon) && (abs(highcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon))
    infmin <- FALSE
    while (!ready && !infmax) {
      nxt <- high*fac
      nxtcvl <- thiscvl(nxt)
      ready <- nxtcvl < highcvl
      if (!ready) {
        low <- high; lowcvl <- highcvl; high <- nxt; highcvl <- nxtcvl
      }
      infmax <- ((abs(lowcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon) && (abs(highcvl - null$cvl) / abs(null$cvl + 0.1) < epsilon))
      infmin <- (fac < 1) && (abs(lowcvl - nxtcvl) / abs(nxtcvl + 0.1) < epsilon) 
    }
    minlambda2 <- min(low, nxt)
    maxlambda2 <- max(low, nxt)
  } else {
    infmax <- infmin <- FALSE
  }
  
  # phase 2: optimize lambda within the order of magnitude found
  if (!infmax && !infmin) {
    opt <- opt.brent(thiscvl, sort(c(minlambda2,maxlambda2)), maximum = TRUE, tol = tol)
  } else {
    if (infmax) {
      opt <- list(argmax = Inf, max = null$cvl)
      names(best$fit$beta) <- colnames(X)
      best$fit$penalty <- c(L1 = 0, L2 = Inf)
    } else {
      best$cvl <- -Inf    # clears bestfit
      best$cvl <- thiscvl(0)
      opt <- list(argmax = 0, max = best$cvl) 
    } 
  } 

  best$predictions <- switch(model, 
    cox = .coxmerge(best$predictions),
    logistic = .logitmerge(best$predictions),
    linear = .lmmerge(best$predictions)
  )

  return(list(lambda = opt$argmax, cvl = opt$max, predictions = best$predictions, fold = groups, 
    fullfit = .makepenfit(best$fit, length(startgamma), model, inputlambda1, opt$argmax, orthogonalizer, weights)))
}
      