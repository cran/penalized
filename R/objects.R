setClass("penfit", 
  representation(
    penalized = "vector", 
    unpenalized = "vector",
    residuals = "vector",
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

.makepenfit <- function(object, unpenalized, model, lambda1, lambda2, orthogonalizer) {
  out <- new("penfit")
  
  beta <- object$beta[unpenalized + seq_len(length(object$beta) - unpenalized)]
  gamma <- object$beta[seq_len(unpenalized)] - as.vector(orthogonalizer %*% beta)
   
  out@unpenalized <- gamma
  out@penalized <- beta
  
  out@residuals <- object$fit$residuals
  
  out@loglik <- object$fit$loglik
  out@penalty <- object$penalty
  
  out@iterations <- object$iter
  
  out@converged <- object$converged
  
  out@model <- model
  
  if (model == "cox") 
    out@nuisance <- list(baseline = object$fit$baseline())
  if (model == "linear")
    out@nuisance <- list(sigma2 = object$fit$sigma2)
  
  out@lambda1 <- lambda1
  out@lambda2 <- lambda2
  
  out
}


setMethod("show", "penfit", function(object) {
  cat("Penalized", object@model, "regression object\n")
  if (object@converged) {
    coefs <- unlist(c(object@penalized, object@unpenalized))
    cat(length(coefs), "regression coefficients")
    if (all(object@lambda1>0)) cat(" of which", sum(coefs != 0), "are non-zero")
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

setMethod("coefficients", "penfit", function(object, which = c("all", "penalized", "unpenalized")) {
  which <- match.arg(which)
  switch(which, 
    all = c(object@unpenalized, object@penalized),
    penalized = object@penalized,
    unpenalized = object@unpenalized)
})

setMethod("residuals", "penfit", function(object, ...) {
  object@residuals
})

setGeneric("basehaz")
setMethod("basehaz", "penfit", function(fit, centered) {
  if (fit@model == "cox") 
    return(fit@nuisance$baseline)
  else
    return(NULL)
})

setGeneric("penalty", function(object, ...) standardGeneric("penalty"))
setMethod("penalty", "penfit", function(object, ...) {
  object@penalty
})

setGeneric("loglik", function(object, ...) standardGeneric("loglik"))
setMethod("loglik", "penfit", function(object, ...) {
  object@loglik
})