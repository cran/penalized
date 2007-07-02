\name{Penalized generalized linear models} 
\alias{penalized} 
\alias{cvl} 
\alias{optL1} 
\alias{optL2} 

\title{Penalized regression} 

\description{Fitting and cross-validating generalized linear models with L1 (lasso) 
and/or L2 (ridge) penalties, or a combination of the two.} 

\usage{ 

penalized (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0, data, 
    model = c("cox", "logistic", "linear"), startbeta, startgamma, steps = 1, 
    epsilon = 1e-10, maxiter, standardize = TRUE, trace = TRUE) 

cvl (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0, data, 
    model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
    groups, epsilon = 1e-8, maxiter, standardize = TRUE, trace = TRUE) 

optL1 (response, penalized, unpenalized, lambda2 = 0, data, 
    model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
    groups, epsilon = 1e-8, maxiter, standardize = TRUE, trace = TRUE, accuracy = 1e-3) 

optL2 (response, penalized, unpenalized, lambda1 = 0, startlambda2 = 1, data, 
    model = c("cox", "logistic", "linear"), startbeta, startgamma, fold, 
    groups, epsilon = 1e-8, maxiter, standardize = TRUE, trace = TRUE, accuracy = 1e-3) 

} 

\arguments{ 
\item{response}{The response variable (vector). This should be a numeric vector for 
linear regression, a \code{\link[survival:Surv]{Surv}} object for Cox regression and 
a vector of 0/1 values for logistic regression.} 

\item{penalized}{The penalized covariates. These may be specified 
either as a matrix or as a (one-sided) \code{\link[stats:formula]{formula}} object. 
See also under \code{data}.} 

\item{unpenalized}{Additional unpenalized covariates. 
Specified as under \code{penalized}. 
Note that an unpenalized intercept is included in the model by default (except in the cox model). 
This can be suppressed by specifying \code{unpenalized = ~0}.}

\item{lambda1, lambda2}{The tuning parameters for L1 and L2 penalization. May be 
vectors if different cavariates are to be penalized differently.} 

\item{data}{A \code{data.frame} used to evaluate \code{response}, and the terms of
 \code{penalized} or \code{unpenalized} when these have been specified as a 
 \code{\link[stats:formula]{formula}} object.} 

\item{model}{The model to be used. If missing, the model will be guessed from the \code{response} input.} 

\item{startbeta}{Starting values for the regression coefficients of the penalized covariates.} 

\item{startgamma}{Starting values for the regression coefficients of the unpenalized covariates.} 

\item{steps}{If greater than 1, the algorithm will fit the model for a range of \code{steps} 
\code{lambda1}-values, starting from the maximal value down to the value of \code{lambda1} specified.
This is useful for making plots as in \code{\link{plotpath}}.} 

\item{epsilon}{The convergence criterion. As in \code{\link[stats:glm.control]{glm}}. 
Convergence is judged separately on the likelihood and on the penalty.} 

\item{maxiter}{The maximum number of iterations allowed. Set by default at 25 
when lambda1 = 0, infinite otherwise.} 

\item{standardize}{If \code{TRUE}, standardizes all penalized covariates to 
unit central L2-norm before applying penalization.} 

\item{trace}{If \code{TRUE}, prints progress information. Note that setting 
\code{trace=TRUE} may slow down the algorithm up to 30 percent (but it often feels quicker)} 

\item{fold}{The fold for cross-validation. Leave-one-out (n-fold) is the default.} 

\item{groups}{Forces the allocation of samples into cross-validation folds. 
Specified as a length \code{n} vector with values in \code{1:fold}.} 

\item{accuracy}{Target relative precision of the optimal L1 or L2 penalty parameter.} 

\item{startlambda2}{The starting value for optimizing lambda2.}
} 

\details{The \code{penalized} function fits  regression models for a given 
combination of L1 and L2 penalty parameters. \code{cvl} calculates the corresponding 
cross-validated likelihood. the twin functions \code{optL1} and \code{optL2} find optimal 
cross-validated values of the L1 and L2 penalty parameter, respectively, when the value 
of the other penalty parameter is given. } 

\value{\code{penalized} returns a \code{link{penfit}} object when \code{steps = 1} 
or a list of such objects if \code{steps > 1}. 

\code{cvl} returns a list with elements
\code{cvl}, the cross-validated likelihood and, \code{fullfit}, a \code{link{penfit}} 
object with the full data fit. 

\code{optL1} and \code{optL2} each return a list with elements \code{lambda}, the 
optimal penalty parameter, \code{cvl}, the associated cross-validated likelihood, 
and \code{fullfit}, a \code{link{penfit}} object with the 
full data fit at the optimal lambda value.}

\note{The functions also accept formula input as in \code{\link{lm}}. The input 
\code{penalized(y~x)} is equivalent to \code{penalized(y, ~x)}.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}} 

\seealso{\code{\link{penfit}}, \code{\link{plotpath}}.} 

\examples{ 

data(nki70)

# A single lasso fit predicting survival
pen <- penalized(Surv(time, event), penalized = nki70[,8:77], unpenalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1 = 10)
show(pen)
coefficients(pen)
coefficients(pen, "penalized")
basehaz(pen)

# A single lasso fit using using the clinical risk factors 
pen <- penalized(Surv(time, event), penalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1=10, standardize=TRUE)

# Find the associated leave-one-out-cross-validated likelihood
cvl(Surv(time, event), penalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1=10, standardize=TRUE)

# Finding an optimal crossvalidated likelihood
attach(nki70)
pen <- optL1(Surv(time, event), penalized = nki70[,8:77], fold = 10)
coefficients(pen$fullfit)
} 

\keyword{multivariate} \keyword{regression} 