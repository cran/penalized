\name{Prediction from penalized models}
\alias{predict,penfit-method}

\title{Prediction based on penfit objects}

\description{Predicting a response for new subjects based on a fitted penalized regression model.}

\usage{

\S4method{predict}{penfit}(object, penalized, unpenalized, data)

}

\arguments{
\item{object}{The fitted model (a \code{\link{penfit}} object).}

\item{penalized}{The penalized covariates for the new subjects. These may be specified either as a matrix or as a (one-sided) \code{\link[stats:formula]{formula}} object. In the latter case, the formula is interpreted in terms of the \code{data} argument.}

\item{unpenalized}{The unpenalized covariates for the new subjects. These may be specified either as a matrix or as a (one-sided) \code{\link[stats:formula]{formula}} object. In the latter case, the formula is interpreted in terms of the \code{data} argument.}

\item{data}{A \code{data.frame} used to evaluate the terms of  \code{penalized} or \code{unpenalized} when these have been specified as a  \code{\link[stats:formula]{formula}} object.}
}

\details{The terms or columns of the penalized and unpenalized arguments must be exactly the same as in the original call that produced the \code{\link{penfit}} object. Any factors in \code{data} must have the same levels.}

\value{The predictions, either as a \code{vector} (logistic and Poisson models), a \code{matrix} (linear model), or a \code{\link{breslow}} object (Cox model).}

\examples{

data(nki70)

pen <- penalized(Surv(time, event), penalized = nki70[1:50,8:77],
    unpenalized = ~ER+Age+Diam+N+Grade, data = nki70[1:50,], lambda1 = 10)

predict(pen, nki70[51:52,8:77], ~ER+Age+Diam+N+Grade, data = nki70[51:52,])
}

\keyword{multivariate}
\keyword{regression}
