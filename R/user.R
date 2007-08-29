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
  if ((p >= n) && all(inputlambda1 == 0) && all(inputlambda2 == 0)) 
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

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
        gams <- .solve(crossprod(t(P)), P %*% beta)
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

  outs
}

######################################
# Finds the cross-validated loglikelihood for a given penalty
######################################
cvl <- function(response, penalized, unpenalized, lambda1 = 0, lambda2= 0, data, 
  model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,  
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
  if ((p >= n) && all(inputlambda1 == 0) && all(inputlambda2 == 0)) 
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

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
  
  res <- .cvl(X, lambda1, lambda2, beta, fit=fit$fit, cvl=fit$cvl, prediction = fit$prediction, 
    groups=groups, epsilon=epsilon, maxiter=maxiter, trace = trace, quit.if.failed = FALSE)
  res$predictions <- switch(model, 
    cox = .coxmerge(res$predictions),
    logistic = .logitmerge(res$predictions),
    linear = .lmmerge(res$predictions)
  )

  return(list(cvl = res$cvl, predictions = res$predictions, fold = groups, fullfit = .makepenfit(res$fit, length(startgamma), model, inputlambda1, inputlambda2, orthogonalizer)))
}



