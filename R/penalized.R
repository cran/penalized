####################################           
# Fits the penalized regression model
####################################                                    
penalized <- function(response, penalized, unpenalized, lambda1=0, lambda2=0, data, model = c("cox", "logistic", "linear"),
  startbeta, startgamma, steps =1, epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE) {
                                    
  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # check for the presence of penalty parameters
  if (ncol(prep$X) >= nrow(prep$X) && lambda1 == 0 && lambda2 == 0)
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

  # prepare the model
  fit <- switch(prep$model,
    cox = .coxfit(prep$response)$fit,
    logistic = .logitfit(prep$response)$fit,
    linear = .lmfit(prep$response)$fit
  )
  
  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)

  # If a steps argument is given, determine where to start
  if (steps > 1) { 
    if (pu > 0) 
      lp <- drop(prep$X[,1:pu,drop=FALSE] %*% prep$nullgamma)
    else 
      lp <- numeric(n)
    gradient <- drop(crossprod(prep$X[,pu+1:pp,drop=FALSE], fit(lp)$residuals))
    rel <- gradient / prep$baselambda1[pu+1:pp]
    from <- max(abs(rel))
  } else {
    from <- lambda1
  }
  lambda1s <- as.list(seq(from, lambda1, length.out=steps+1))
  if (steps == 1) lambda1s <- lambda1s[-1]

  # fit the model for all lambdas
  beta <- prep$beta
  outs <- lapply(lambda1s, function(rellambda) {
    if (!(rellambda == 0)) {
      if (lambda2 == 0) {
        out <- .steplasso(beta = beta, lambda = rellambda * prep$baselambda1, 
          lambda2 = 0, X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
          maxiter = maxiter)
      } else {
        out <- .lasso(beta = beta, lambda = rellambda * prep$baselambda1, 
          lambda2 = lambda2 * prep$baselambda2, X = prep$X, fit = fit, 
          trace = trace, epsilon = epsilon, maxiter = maxiter)
      }
    } else {
      if (pp > n) {
        P <- .makeP(prep$X, lambda2 * prep$baselambda2)
        gams <- .solve(crossprod(t(P)), P %*% beta)
        PX <- P %*% t(prep$X)
        Pl <- P * matrix(sqrt(lambda2 * prep$baselambda2), nrow(P), ncol(P), 
          byrow = TRUE)
        PlP <- crossprod(t(Pl))
        # Do ridge regression on gamma
        out <- .ridge(beta = gams, Lambda = PlP, X = t(PX), fit = fit, 
          trace = trace, epsilon = epsilon, maxiter = maxiter)
        # and transform back
        out$beta <- drop(crossprod(P, out$beta)) 
      } else {
        out <- .ridge(beta = beta, Lambda = lambda2 * prep$baselambda2, 
          X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
          maxiter = maxiter)
      }
    }
    if (trace) cat("\n")
    beta <<- out$beta

    out$beta <- out$beta / prep$weights
    out
  })
  
  # put the output in a penfit object
  outs <- sapply(1:length(outs), function(nr) {
    thislambda1 <- lambda1s[[nr]]
    .makepenfit(outs[[nr]], pu, prep$model, thislambda1, lambda2, 
      prep$orthogonalizer, prep$weights)
  })

  if(length(outs)==1) 
    outs <- outs[[1]]

  outs
}





