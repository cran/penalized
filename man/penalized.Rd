\name{Penalized generalized linear models}
\alias{penalized}

\title{Penalized regression}

\description{Fitting generalized linear models with L1 (lasso and fused lasso)
and/or L2 (ridge) penalties, or a combination of the two.}

\usage{

penalized (response, penalized, unpenalized, lambda1=0, 
  lambda2=0, positive = FALSE, data, fusedl=FALSE,
  model = c("cox", "logistic", "linear", "poisson"), 
  startbeta, startgamma, steps =1, epsilon = 1e-10, 
  maxiter, standardize = FALSE, trace = TRUE)

}

\arguments{
\item{response}{The response variable (vector). This should be a numeric vector for linear regression, a \code{\link[survival:Surv]{Surv}} object for Cox regression and \code{factor} or a vector of 0/1 values for logistic regression.}

\item{penalized}{The penalized covariates. These may be specified either as a matrix or as a (one-sided) \code{\link[stats:formula]{formula}} object. See also under \code{data}.}

\item{unpenalized}{Additional unpenalized covariates. Specified as under \code{penalized}. Note that an unpenalized intercept is included in the model by default (except in the Cox model). This can be suppressed by specifying \code{unpenalized = ~0}.}

\item{lambda1, lambda2}{The tuning parameters for L1 and L2 penalization. Each must be either a single positive numbers or a vector with length equal to the number of covariates in \code{penalized} argument. In the latter case, each covariate is given its own penalty weight.}

\item{positive}{If \code{TRUE}, constrains the estimated regression coefficients of all penalized covariates to be non-negative. If a logical vector with the length of the number of covariates in \code{penalized}, constrains the estimated regression coefficients of a subset of the penalized covariates to be non-negative.}

\item{data}{A \code{data.frame} used to evaluate \code{response}, and the terms of  \code{penalized} or \code{unpenalized} when these have been specified as a  \code{\link[stats:formula]{formula}} object.}

\item{fusedl}{If \code{TRUE} or a vector, the penalization method used is fused lasso. The value for \code{lambda1} is used as the tuning parameter for L1 penalization on the coefficients and the value for \code{lambda2} is used as the tuning parameter for L1 penalization on the differences of the coefficients. Default value is \code{FALSE}.}

\item{model}{The model to be used. If missing, the model will be guessed from the \code{response} input.}

\item{startbeta}{Starting values for the regression coefficients of the penalized covariates.}

\item{startgamma}{Starting values for the regression coefficients of the unpenalized covariates.}

\item{steps}{If greater than 1, the algorithm will fit the model for a range of \code{steps} \code{lambda1}-values, starting from the maximal value down to the value of \code{lambda1} specified. This is useful for making plots as in \code{\link{plotpath}}. With \code{steps = "Park"} it is possible to choose the steps in such a way that they are at the approximate value at which the active set changes, following Park and Haste (2007).}

\item{epsilon}{The convergence criterion. As in \code{\link[stats:glm.control]{glm}}. Convergence is judged separately on the likelihood and on the penalty.}

\item{maxiter}{The maximum number of iterations allowed. Set by default at 25 when only an L2 penalty is present, infinite otherwise.}

\item{standardize}{If \code{TRUE}, standardizes all penalized covariates to unit central L2-norm before applying penalization.}

\item{trace}{If \code{TRUE}, prints progress information. Note that setting \code{trace=TRUE} may slow down the algorithm up to 30 percent (but it often feels quicker)}
}

\details{The \code{penalized} function fits  regression models for a given combination of L1 and L2 penalty parameters.}

\value{\code{penalized} returns a \code{\link{penfit}} object when \code{steps = 1} or a list of such objects if \code{steps > 1}.}

\note{The \code{response} argument of the function also accepts formula input as in  \code{\link{lm}} and related functions. In that case, the right hand side of the \code{response} formula is used as the \code{penalized} argument or, if that is already given, as the \code{unpenalized} argument. For example, the input \code{penalized(y~x)} is equivalent to \code{penalized(y, ~x)} and \code{penalized(y~x, ~z)} to \code{penalized(y, ~z, ~x)}.

In case of tied survival times, the function uses Breslow's version of the partial likelihood.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\references{Goeman J.J. (2010). L-1 Penalized Estimation in the Cox Proportional Hazards Model. Biometrical Journal 52 (1) 70-84.}

\seealso{\code{\link{penfit}} for the \code{penfit} object returned, \code{\link{plotpath}}
for plotting the solution path, and \code{\link{cvl}} for cross-validation and
optimizing the tuning parameters.}

\examples{
# More examples in the package vignette:
#  type vignette("penalized")

data(nki70)

# A single lasso fit predicting survival
pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
    unpenalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1 = 10)
show(pen)
coefficients(pen)
coefficients(pen, "penalized")
basehaz(pen)

# A single lasso fit using the clinical risk factors
pen <- penalized(Surv(time, event), penalized = ~ER+Age+Diam+N+Grade,
    data = nki70, lambda1=10, standardize=TRUE)

# using steps
pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
    data = nki70, lambda1 = 1,steps = 20)
plotpath(pen)


# A fused lasso fit predicting survival
pen <- penalized(Surv(time, event), penalized = nki70[,8:77], data = nki70, 
     lambda1 = 1, lambda2 = 2, fusedl = TRUE)
plot(coefficients(pen, "all"),type="l",xlab = "probes",ylab = "coefficient value")
plot(predict(pen,penalized=nki70[,8:77]))
 }

\keyword{multivariate}
\keyword{regression} 