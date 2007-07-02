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

    return(list(residuals = residuals, loglik = loglik, W = ws, lp = lp))
  }

  cvl <- function(lp, leftout) {
    probs <- exp(lp) / (1+exp(lp))
    return(sum(log(probs[response == 1 & leftout])) + sum(log(1-probs[response == 0 & leftout])))
  }

  return(list(fit = fit, cvl = cvl))
}


.logitgamma <- function(response, unpenalized, data) {

  if (is.matrix(unpenalized)) {
    if (ncol(unpenalized) > 0) {
      .response <- response
      startgamma <- coefficients(glm(.response ~ ., data = as.data.frame(unpenalized), family = binomial))
    } else                                                   
      startgamma <- mean(response)
  } else {
    .response <- response
    terms <- c(attr(terms(unpenalized),"intercept"), attr(terms(unpenalized),"term.labels"))
    form <- as.formula(paste(".response~", paste(terms, collapse="+")))
    startgamma <- coefficients(glm(form, data = data, family = binomial))
  }
  return(startgamma)
}