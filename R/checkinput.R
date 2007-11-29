########################################################
# A function that does the common input checking of
# the functions penalized, cvl, profL1, profL2, optL1 and optL2.
########################################################
.checkinput <- function(call, env) {
                      
  # Functions to extract the original input variables
  call <- as.list(call)
  input <- function(str) eval(call[[str]], env)
  input.data <- function(str) eval(call[[str]], data, env)
  missing <- function(str) !str %in% names(call)

  # determine the response
  if (!missing("data")) data <- input("data") else data <- NULL
  response <- input.data("response")
  if (is(response, "formula")) {
    formula.response <- response
    if (missing("data"))
      response <- eval(attr(terms(response), "variables"), environment(response))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response), "variables"), data, environment(response))[[attr(terms(response), "response")]]
  } else {
    formula.response <- NULL
  }

  # determine the model if missing
  if (missing("model")) {
    if (is(response, "Surv")) model <- "cox"
    else if (all(response %in% 0:1)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  } else {
    model <- match.arg(input("model"), c("cox", "logistic", "linear"))
  }

  # determine penalized and unpenalized
  if (!missing("penalized")) {
    penalized <- input.data("penalized")
  } else {
    if (!is.null(formula.response)) {
      penalized <- formula.response
      formula.response <- NULL
    } else
      stop("argument \"penalized\" is missing, with no default")
  }
  if (!missing("unpenalized")) {
    unpenalized <- input.data("unpenalized")
  } else {
    if (!is.null(formula.response)) {
      unpenalized <- formula.response
      formula.response <- NULL
    } else {
      unpenalized <- response ~ 1
    }
  }

  # Has the response formula been used?
  if (!is.null(formula.response))
    warning("right hand side of response formula ignored")

  # coerce unpenalized into a matrix
  if (is.data.frame(unpenalized) || is.vector(unpenalized))
    unpenalized <- as.matrix(unpenalized)
  if (is(unpenalized, "formula")) {
    unpenalized <- terms(unpenalized)
    # prevent problems for input ~1 or ~0:
    if ((attr(unpenalized, "response") == 0) && (length(attr(unpenalized, "term.labels")) == 0)) {
      if (attr(unpenalized, "intercept") == 1)
        unpenalized <- terms(response ~ 1)
      else
        unpenalized <- terms(response ~ 0)
    }
    unpenalized <- model.matrix(unpenalized, data)
    # suppress intercept if necessary
    if (model == "cox" ) unpenalized <- unpenalized[,-1,drop=FALSE]
  }
  
  # coerce penalized into a matrix
  if (is.data.frame(penalized) || is.vector(penalized))
    penalized <- as.matrix(penalized)
  if (is(penalized, "formula")) {
    oldcontrasts <- unlist(options("contrasts"))
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    penalized <- terms(penalized)
    # prevent problems for input ~1 or ~0:
    if (length(attr(penalized, "term.labels")) == 0) 
      penalized <- terms(response ~ 1)
    attr(penalized, "intercept") <- 1
    penalized <- model.matrix(penalized, data)
    options(contrasts = oldcontrasts)
    penalized <- penalized[,-1,drop=FALSE]
  }

  # check dimensions of response, penalized and unpenalized
  n <- if (model == "cox") length(response)/2 else length(response)
  if (nrow(penalized) != n) {
    stop("the length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  if (nrow(unpenalized) != n) {
    stop("the length of \"response\" (",n, ") does not match the row count of \"unpenalized\" (", nrow(penalized), ")")
  }

  # get the value of startbeta
  if (missing("startbeta"))
    startbeta <- numeric(ncol(penalized))
  else {
    startbeta <- input("startbeta")
    if (length(startbeta) != ncol(penalized))
      stop("The length of \"startbeta\" (", length(startbeta), ") does not match the column count of \"penalized\" (", ncol(penalized), ")")
  }
  if (is.null(names(startbeta)))
    names(startbeta) <- colnames(penalized)

  # get the value of startgamma
  if (ncol(unpenalized) > 0) {
    nullgamma <- switch(model,
      cox = coefficients(coxph(response ~ unpenalized)),
      logistic = coefficients(glm(response ~ 0 + unpenalized, family = binomial)),
      linear = coefficients(lm(response ~ 0 + unpenalized))
    )
    names(nullgamma) <- colnames(unpenalized)
  } else nullgamma <- numeric(0)
  if (missing("startgamma")) {
    startgamma <- nullgamma
  } else {
    startgamma <- input("startgamma")
    if (length(startgamma) != ncol(unpenalized))
      stop("The length of \"startgamma\" (", length(startgamma), ") does not match the column count of \"unpenalized\" (", ncol(unpenalized), ")")
  }
  if (is.null(names(startgamma)))
    names(startgamma) <- colnames(unpenalized)


  # orthogonalize penalized with respect to unpenalized
  if (ncol(unpenalized) > 0) {
    orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
    penalized <- penalized - unpenalized %*% orthogonalizer
  } else {
    orthogonalizer <- matrix(,0,ncol(penalized))
  }

  # Join penalized and unpenalized together
  X <- cbind(unpenalized, penalized)
  beta <- c(startgamma, startbeta)

  # stabilize/standardize
  vars <- apply(X,2,var) * (n-1)/n
  vars[vars == 0] <- 1
  sds <- sqrt(vars)
  X <- X / matrix(sds, nrow(X), ncol(X), byrow=T)
  standardize <- if (missing("standardize")) FALSE else input("standardize")
  beta[beta != 0] <- beta[beta != 0] * sds[beta != 0]
  nullgamma <- nullgamma * sds[1:length(nullgamma)]
  
  # find baselambda1 and baselambda2 
  # This lambda1 and lambda2 for unit input lambda1=1 and lambda2=1
  if (standardize) {
    baselambda1 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
    baselambda2 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
  } else {
    sel <- ncol(unpenalized) + 1:ncol(penalized)
    baselambda1 <- c(numeric(ncol(unpenalized)), 1/sds[sel])
    baselambda2 <- c(numeric(ncol(unpenalized)), 1/vars[sel])
  }

  return(list(
    response = response,
    X = X, 
    beta = beta, 
    weights = sds, 
    baselambda1 = baselambda1, 
    baselambda2 = baselambda2,
    orthogonalizer = orthogonalizer, 
    model = model, 
    nullgamma = nullgamma
  ))
}