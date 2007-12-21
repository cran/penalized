.logitfit <- function(response) {

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {

    if (!missing(leftout))
      response <- response[!leftout]
    probs <- exp(lp) / (1+exp(lp))
    ws <- probs * (1-probs)

    # The residuals
    residuals <- response - probs

    # The loglikelihood
    loglik <- sum(log(probs[response == 1])) + sum(log(1-probs[response == 0]))
    if (!is.na(loglik) && (loglik == - Inf)) loglik <- NA

    return(list(residuals = residuals, loglik = loglik, W = ws, lp = lp, fitted = probs, nuisance = list()))
  }

  cvl <- function(lp, leftout) {
    probs <- exp(lp) / (1+exp(lp))
    return(sum(log(probs[response == 1 & leftout])) + sum(log(1-probs[response == 0 & leftout])))
  }

  # mapping from the linear predictor lp to an actual prediction
  prediction <- function(lp, nuisance) {
    out <- exp(lp) / (1+exp(lp))
    out
  }


  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


# merges predicted probalities
.logitmerge <- function(predictions) {
  out <- unlist(predictions)
  out
}
