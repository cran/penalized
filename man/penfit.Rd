\name{penfit object} 
\alias{penfit-class} 
\alias{penfit} 
\docType{class}

\alias{show,penfit-method} 
\alias{coefficients,penfit-method} 
\alias{residuals,penfit-method} 
\alias{fitted.values,penfit-method} 
\alias{penalty,penfit-method} 
\alias{loglik,penfit-method} 
\alias{basehaz,penfit-method} 
\alias{basehaz} 
\alias{penalty} 
\alias{loglik}


\title{Penalized regression object}

\description{Stores the result of a call to \code{\link{penalized}}.}

\section{Slots}{
  \describe{
    \item{\code{res}:}{Object of class "matrix". Test results summary.}
    \item{\code{penalized}:}{Object of class "vector". Regression coefficients for the penalized covariates.}
    \item{\code{unpenalized}:}{Object of class "vector". Regression coefficients for the unpenalized covariates.}
    \item{\code{residuals}:}{Object of class "vector". Unstandardized residuals of the fitted model. Martingale residuals are given for the cox model.}
    \item{\code{loglik}:}{Object of class "numeric". Log likelihood of the fitted model.}
    \item{\code{penalty}:}{Object of class "vector". L1 and L2 penalties of the fitted model.}
    \item{\code{iterations}:}{Object of class "numeric". Number of iterations used in the fitting process.}
    \item{\code{converged}:}{Object of class "logical". Whether the fitting process was judged to be converged.}
    \item{\code{model}:}{Object of class "character". The name of the generalized linear model used.}
    \item{\code{lambda1}:}{Object of class "vector". The lambda1 parameter(s) used.}
    \item{\code{lambda2}:}{Object of class "vector". The lambda2 parameter(s) used.}
    \item{\code{nuisance}:}{Object of class "list". The maximum likelihood estimates of any nuisance parameters in the model.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{(penfit): Summarizes the test result.}
    \item{residuals}{(penfit): Returns the residuals.}
    \item{coefficients}{(penfit): Returns the regression coefficients. Accepts a second argument "which", that takes
    values "nonzero" (the default), "all", "penalized" or "unpenalized" for extracting only the non-zero, the penalized
    or the unpenalized regression coeffiicents.}
    \item{basehaz}{(penfit): Returns the baseline hazard (a \code{\link{breslow}} object) if a cox model was fitted, \code{NULL} otherwise.}
    \item{penalty}{(penfit): Returns the L1 and L2 penalties of the fitted model.}
    \item{loglik}{(penfit): Returns the log likelihood of the fitted model.}
  }
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}} 

\seealso{\code{\link{penalized}}.}

\keyword{multivariate} 
\keyword{regression} 