.lmfit <- function(response) {

  # Finds local gradient and subject weights
  fit <- function(lp, leftout) {
    if (!missing(leftout))
      response <- response[!leftout]

    # The residuals
    residuals <- response - lp

    # The loglikelihood
    ss <- sum(residuals * residuals)
    if (missing(leftout)) n <- length(lp) else n <- sum(!leftout)
    loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin))

    return(list(residuals = residuals, loglik = loglik, W = 1, lp = lp, fitted = lp, nuisance = list(sigma2 = ss/n), lp = lp))
  }

  cvl <- function(lp, leftout) {

    residuals <- response - lp
    sigma2 <- sum(residuals[!leftout] * residuals[!leftout]) / sum(!leftout)
    ss <- sum(residuals[leftout] * residuals[leftout])

    return(-(sum(leftout)/2) * log(2*pi*sigma2) - ss / (2*sigma2))
  }
  
  # mapping from the linear predictor lp to an actual prediction
  prediction <- function(lp, nuisance) {
    out <- c(mu = lp, sigma2 = nuisance$sigma2)
    out
  }


  return(list(fit = fit, cvl = cvl, prediction = prediction))
}


.lmgamma <- function(response, unpenalized, data) {

  if (is.matrix(unpenalized)) {
    if (ncol(unpenalized) > 0) {
      .response <- response
      startgamma <- coefficients(lm(.response ~ ., data = as.data.frame(unpenalized)))
    } else
      startgamma <- mean(response)
  } else {
    .response <- response
    terms <- c(attr(terms(unpenalized),"intercept"), attr(terms(unpenalized),"term.labels"))
    form <- as.formula(paste(".response~", paste(terms, collapse="+")))
    startgamma <- coefficients(lm(form, data = data))
  }
}

# merges predicted means and variances
.lmmerge <- function(predictions) {
  out <- matrix(unlist(predictions), length(predictions), 2, byrow=TRUE)
  colnames(out) <- c("mu", "sigma2")
  rownames(out) <- names(predictions)
  out
}
