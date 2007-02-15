# "penfit" object to store the result of a penalized regression
# Theseobjects are meant to be accessed by the user, but not made directly by the user
setClass("penfit", 
  representation(
    penalized = "vector", 
    unpenalized = "vector",
    residuals = "vector",
    fitted = "vector",
    lin.pred = "vector",
    loglik = "numeric",
    penalty = "vector",
    iterations = "numeric",
    converged = "logical",
    model = "character",
    lambda1 = "vector",
    lambda2 = "vector",
    nuisance = "list",
    weights = "vector" 
  )
)

# creation method for a penfit object 
.makepenfit <- function(object, unpenalized, model, lambda1, lambda2, orthogonalizer, weights) {
  out <- new("penfit")
  
  beta <- object$beta[unpenalized + seq_len(length(object$beta) - unpenalized)]
  gamma <- object$beta[seq_len(unpenalized)] - as.vector(orthogonalizer %*% beta)
   
  out@unpenalized <- gamma
  out@penalized <- beta
  
  out@residuals <- object$fit$residuals
  out@fitted <- object$fit$fitted
  out@lin.pred <- object$fit$lp
  
  out@loglik <- if (is.na(object$fit$loglik)) -Inf else object$fit$loglik
  out@penalty <- object$penalty
  
  out@iterations <- object$iter
  
  out@converged <- object$converged
  
  out@model <- model
  
  out@nuisance <- object$fit$nuisance
  
  out@lambda1 <- lambda1
  out@lambda2 <- lambda2
  
  out@weights <- weights
  
  out
}

# show method
setMethod("show", "penfit", function(object) {
  cat("Penalized", object@model, "regression object\n")
  if (object@converged) {
    coefs <- unlist(c(object@penalized, object@unpenalized))
    cat(length(coefs), "regression coefficients")
    if (any(coefs == 0)) cat(" of which", sum(coefs != 0), "are non-zero")
    cat("\n\n")
    cat("Loglikelihood =\t", object@loglik, "\n")
    if (any(object@lambda1 > 0))
      if (length(object@lambda1) > 3) 
        cat("L1 penalty =\t", object@penalty[1], "\tat lambda1 = ", object@lambda1[1:3], "...\n")
      else
        cat("L1 penalty =\t", object@penalty[1], "\tat lambda1 = ", object@lambda1, "\n") 
    if (any(object@lambda2 > 0))
      if (length(object@lambda2) > 3) 
        cat("L1 penalty =\t", object@penalty[2], "\tat lambda2 = ", object@lambda2[1:3], "...\n")
      else
        cat("L2 penalty =\t", object@penalty[2], "\tat lambda2 = ", object@lambda2, "\n") 
  } else {
    cat("Model failed to converge\n")
  }
})

# extracts the coefficients
setMethod("coefficients", "penfit", function(object, which = c("nonzero", "all", "penalized", "unpenalized"), standardize = FALSE) {
  which <- match.arg(which)
  nunp <- length(object@unpenalized)
  np <- length(object@penalized)
  whichunp <- switch(which, 
    all =, unpenalized =, nonzero = rep(TRUE,nunp),
    penalized = rep(FALSE, nunp))
  whichp <- switch(which,
    all =, penalized = rep(TRUE, np),
    unpenalized = rep(FALSE, np),
    nonzero = (object@penalized != 0))
  out <- c(object@unpenalized[whichunp], object@penalized[whichp])
  if (standardize) out <- out * object@weights[c(whichunp, whichp)]
  out
})

# extracts the residuals
setMethod("residuals", "penfit", function(object, ...) {
  object@residuals
})

# extracts the linear predictors
setGeneric("linear.predictors", function(object, ...) standardGeneric("linear.predictors"))
setMethod("linear.predictors", "penfit", function(object, ...) {
  object@lin.pred
})

# extracts the fitted values
setMethod("fitted.values", "penfit", function(object, ...) {
  object@fitted
})

# extracts the weights
setMethod("weights", "penfit", function(object, ...) {
  object@weights
})

# extracts the baseline hazard (survival models only)
setGeneric("basesurv", function(fit, centered = TRUE, ...) standardGeneric("basesurv"))
setMethod("basesurv", "penfit", function(fit, centered = TRUE) {
  if (fit@model == "cox") 
    if (centered) {
      meanlp <- mean(linear.predictors(fit))  
      out <- fit@nuisance$baseline
      out@curves <- out@curves^exp(meanlp)
      return(out)
    } else {
      return(fit@nuisance$baseline)
    }
  else
    return(NULL)
})

setGeneric("basehaz")
setMethod("basehaz", "penfit", function(fit, centered = TRUE) {
  if (fit@model == "cox") {
    bs <- basesurv(fit, centered)
    out <- data.frame(hazard = -log(drop(bs@curves)), time = time(bs))
    return(out)
  } else
    return(NULL)
})

# extracts the penalty
setGeneric("penalty", function(object, ...) standardGeneric("penalty"))
setMethod("penalty", "penfit", function(object, ...) {
  object@penalty
})

# extracts the likelihood
setGeneric("loglik", function(object, ...) standardGeneric("loglik"))
setMethod("loglik", "penfit", function(object, ...) {
  object@loglik
})
