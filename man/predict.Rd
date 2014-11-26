\name{Prediction from penalized models}
\alias{predict,penfit-method}
\alias{predict}

\title{Prediction based on penfit objects}

\description{Predicting a response for new subjects based on a fitted penalized regression model.}

\usage{

\S4method{predict}{penfit}(object, penalized, unpenalized, data)

}

\arguments{
\item{object}{The fitted model (a \code{\link{penfit}} object).}

\item{penalized}{The matrix of penalized covariates for the new subjects.}

\item{unpenalized}{The unpenalized covariates for the new subjects.}

\item{data}{A \code{data.frame} used to evaluate the terms of  \code{penalized} or \code{unpenalized} when these have been specified as a  \code{\link[stats:formula]{formula}} object.}
}

\details{The user need only supply those terms from the original call that are different relative to the original call that produced the \code{\link{penfit}} object. In particular, if penalized and/or unpenalized was specified in matrix form, a matrix must be given with the new subjects' data. The columns of these matrices must be exactly the same as in the matrices supplied in the original call that produced the \code{\link{penfit}} object. If either penalized or unpenalized was given as a \code{formula} in the original call, the user of \code{predict} must supply a new \code{data} argument. As with matrices, the new \code{data} argument must have a similar make-up as the \code{data} argument in the original call that produced the \code{\link{penfit}} object. In particular, any factors in \code{data} must have the same levels.}

\value{The predictions, either as a \code{vector} (logistic and Poisson models), a \code{matrix} (linear model), or a \code{\link{breslow}} object (Cox model).}

\examples{

data(nki70)

pen <- penalized(Surv(time, event), penalized = nki70[1:50,8:77],
    unpenalized = ~ER+Age+Diam+N+Grade, data = nki70[1:50,], lambda1 = 10)

predict(pen, nki70[51:52,8:77], data = nki70[51:52,])
}

\keyword{multivariate}
\keyword{regression}
