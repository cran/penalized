# "penfit" object to store the result of a penalized regression
# Theseobjects are meant to be accessed by the user, but not made directly by the user
setClass("penfit", 
  representation(
    penalized = "vector", 
    unpenalized = "vector",
    residuals = "vector",
    fitted = "vector",
    loglik = "numeric",
    penalty = "vector",
    iterations = "numeric",
    converged = "logical",
    model = "character",
    lambda1 = "vector",
    lambda2 = "vector",
    nuisance = "list" 
  )
)

# creation method for a penfit object 
.makepenfit <- function(object, unpenalized, model, lambda1, lambda2, orthogonalizer) {
  out <- new("penfit")
  
  beta <- object$beta[unpenalized + seq_len(length(object$beta) - unpenalized)]
  gamma <- object$beta[seq_len(unpenalized)] - as.vector(orthogonalizer %*% beta)
   
  out@unpenalized <- gamma
  out@penalized <- beta
  
  out@residuals <- object$fit$residuals
  out@fitted <- object$fit$fitted
  
  out@loglik <- if (is.na(object$fit$loglik)) -Inf else object$fit$loglik
  out@penalty <- object$penalty
  
  out@iterations <- object$iter
  
  out@converged <- object$converged
  
  out@model <- model
  
  out@nuisance <- object$fit$nuisance
  
  out@lambda1 <- lambda1
  out@lambda2 <- lambda2
  
  out
}

# show method
setMethod("show", "penfit", function(object) {
  cat("Penalized", object@model, "regression object\n")
  if (object@converged) {
    coefs <- unlist(c(object@penalized, object@unpenalized))
    cat(length(coefs), "regression coefficients")
    if (any(object@lambda1>0) || any(object@lambda2 == Inf)) cat(" of which", sum(coefs != 0), "are non-zero")
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
setMethod("coefficients", "penfit", function(object, which = c("nonzero", "all", "penalized", "unpenalized")) {
  which <- match.arg(which)
  switch(which, 
    all = c(object@unpenalized, object@penalized),
    penalized = object@penalized,
    unpenalized = object@unpenalized,
    nonzero = c(object@unpenalized, object@penalized[object@penalized != 0]))
})

# extracts the residuals
setMethod("residuals", "penfit", function(object, ...) {
  object@residuals
})

# extracts the fitted values
setMethod("fitted.values", "penfit", function(object, ...) {
  object@fitted
})


# extracts the baseline hazard (survival models only)
setGeneric("basehaz")
setMethod("basehaz", "penfit", function(fit, centered) {
  if (fit@model == "cox") 
    return(fit@nuisance$baseline)
  else
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
