\name{penfit object}
\alias{penfit-class}
\alias{penfit}
\docType{class}

\alias{show,penfit-method}
\alias{coefficients,penfit-method}
\alias{coef,penfit-method}
\alias{residuals,penfit-method}
\alias{fitted.values,penfit-method}
\alias{fitted,penfit-method}
\alias{penalty,penfit-method}
\alias{weights,penfit-method}
\alias{loglik,penfit-method}
\alias{basehaz,penfit-method}
\alias{basesurv,penfit-method}
\alias{basesurv}
\alias{linear.predictors}
\alias{linear.predictors,penfit-method}
\alias{as.data.frame,penfit-method}
\alias{penalty}
\alias{loglik}


\title{Penalized regression object}

\description{Stores the result of a call to \code{\link{penalized}}.}

\section{Slots}{
  \describe{
    \item{\code{penalized}:}{Object of class "vector". Regression coefficients for the penalized covariates.}
    \item{\code{unpenalized}:}{Object of class "vector". Regression coefficients for the unpenalized covariates.}
    \item{\code{residuals}:}{Object of class "vector". Unstandardized residuals for each subject in the fitted model. Martingale residuals are given for the cox model.}
    \item{\code{fitted}:}{Object of class "vector". Fitted values (means) for each subject in the fitted model. In the cox model, this slot holds the relative risks.}
    \item{\code{lin.pred}:}{Object of class "vector". Linear predictors for each subject in the fitted model.}
    \item{\code{loglik}:}{Object of class "numeric". Log likelihood of the fitted model. For the Cox model, reports the full likelihood rather than the partial likelihood.}
    \item{\code{penalty}:}{Object of class "vector". L1 and L2 penalties of the fitted model.}
    \item{\code{iterations}:}{Object of class "numeric". Number of iterations used in the fitting process.}
    \item{\code{converged}:}{Object of class "logical". Whether the fitting process was judged to be converged.}
    \item{\code{model}:}{Object of class "character". The name of the generalized linear model used.}
    \item{\code{lambda1}:}{Object of class "vector". The lambda1 parameter(s) used.}
    \item{\code{lambda2}:}{Object of class "vector". The lambda2 parameter(s) used.}
    \item{\code{nuisance}:}{Object of class "list". The maximum likelihood estimates of any nuisance parameters in the model.}
    \item{\code{weights}:}{Object of class "vector". The weights of the covariates used for standardization.}
    \item{\code{formula}:}{Object of class "list". A named list containing the unpenalized and penalized formula objects, if present.}
  }
}

\section{Methods}{
  \describe{
    \item{basehaz}{(penfit): Returns the baseline hazard (a \code{data.frame}) if a cox model was fitted, \code{NULL} otherwise. An additional argument \code{center} (default (\code{TRUE}) can be used to give the survival curve at the covariate mean (\code{center = TRUE}) rather than at zero.}
    \item{basesurv}{(penfit): Returns the baseline survival curve (a \code{\link{breslow}} object) if a cox model was fitted, \code{NULL} otherwise. An additional argument \code{center} (default (\code{TRUE}) can be used to give the survival curve at the covariate mean (\code{center = TRUE}) rather than at zero.}
    \item{coef}{(penfit): Returns the regression coefficients. Accepts a second argument "which", that takes values "nonzero" (the default), "all", "penalized" or "unpenalized" for extracting only the non-zero, the penalized or the unpenalized regression coefficients. A third argument "standardize" (default FALSE) can be used to let the method return the regression coefficients for the standardized covariates.}
    \item{coefficients}{(penfit): synonym for \code{coef} above.}
    \item{fitted}{(penfit): Returns the fitted values for each subject (i.e. the predicted means). In the Cox model, this method returns the relative risks for each individual.}
    \item{fitted.values}{(penfit): synonym for \code{fitted} above.}
    \item{linear.predictors}{(penfit): Returns the linear predictors for each subject.}
    \item{loglik}{(penfit): Returns the log likelihood of the fitted model.}
    \item{penalty}{(penfit): Returns the L1 and L2 penalties of the fitted model.}
    \item{residuals}{(penfit): Returns the residuals.}
    \item{show}{(penfit): Summarizes the fitted model.}
    \item{weights}{(penfit): Returns the weights used for standardization.}
    \item{predict}{(penfit): Calculates predictions for new subjects. See \code{\link[penalized:predict]{predict}}.}
  }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{penalized}}.}

\keyword{multivariate}
\keyword{regression} 