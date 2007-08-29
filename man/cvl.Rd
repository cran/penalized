\name{Cross-validation in penalized generalized linear models}
\alias{cvl}
\alias{profL1}
\alias{profL2}
\alias{optL1}
\alias{optL2}

\title{Penalized regression}

\description{Cross-validating generalized linear models with L1 (lasso)
and/or L2 (ridge) penalties, using likelihood cross-validation.}

\usage{

cvl (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0, data,
    model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,
    epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE)

optL1 (response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0,
    data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,
    epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
    tol = .Machine$double.eps^0.25)

optL2 (response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2,
    data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,
    epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
    tol = .Machine$double.eps^0.25)

profL1 (response, penalized, unpenalized, minlambda1, maxlambda1, lambda2 = 0,
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,
  epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
  steps = 100, autominsteps = steps/5, log = FALSE)
  
profL2 (response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2,
  data, model = c("cox", "logistic", "linear"), startbeta, startgamma, fold,
  epsilon = 1e-10, maxiter, standardize = FALSE, trace = TRUE,
  steps = 100, autominsteps = steps/5, log = TRUE)
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

\item{lambda1, lambda2}{The fixed values of the tuning parameters for L1 and L2
penalization. Both may be vectors if different cavariates are to be penalized differently.}

\item{minlambda1, minlambda2, maxlambda1, maxlambda2}{The values of the tuning
parameters for L1 or L2 penalization between which the cross-validated likelihood
is to be profiled or optimized.}

\item{data}{A \code{data.frame} used to evaluate \code{response}, and the terms of
 \code{penalized} or \code{unpenalized} when these have been specified as a
 \code{\link[stats:formula]{formula}} object.}

\item{model}{The model to be used. If missing, the model will be guessed from the \code{response} input.}

\item{startbeta}{Starting values for the regression coefficients of the penalized covariates.
These starting values will be used only for the first values of \code{lambda1} and \code{lambda2}.}

\item{startgamma}{Starting values for the regression coefficients of the unpenalized covariates.
These starting values will be used only for the first values of \code{lambda1} and \code{lambda2}.}

\item{fold}{The fold for cross-validation. May be supplied as a single number
(between 2 and n) giving the number of folds, or, alternatively, as a length \code{n}
vector with values in \code{1:fold}, specifying exactly which subjects are assigned
to which fold. The default is \code{fold = 1:n}, resulting in leave-one-out (n-fold)
cross-validation.}

\item{epsilon}{The convergence criterion. As in \code{\link[stats:glm.control]{glm}}.
Convergence is judged separately on the likelihood and on the penalty.}

\item{maxiter}{The maximum number of iterations allowed. Set by default at 25
when lambda1 = 0, infinite otherwise.}

\item{standardize}{If \code{TRUE}, standardizes all penalized covariates to
unit central L2-norm before applying penalization.}

\item{trace}{If \code{TRUE}, prints progress information. Note that setting
\code{trace=TRUE} may slow down the algorithm (but it often feels quicker)}

\item{steps}{The maximum number of steps between \code{minlambda1}/\code{minlambda1} and
\code{maxlambda1}/\code{maxlambda1} at which the cross-validated likelihood is to
be calculated.}

\item{autominsteps}{The minimum number of steps between \code{minlambda1}/\code{minlambda1} and
\code{maxlambda1}/\code{maxlambda1} at which the cross-validated likelihood is to
be calculated. If \code{autominsteps} is smaller than \code{steps}, the algorithm will
automatically stop when the cross-validated likelihood drops below the cross-validated
likelihood of the null model, provided it has done at least \code{autominsteps} steps.}

\item{log}{If \code{FALSE}, the steps between \code{minlambda1}/\code{minlambda1} and
\code{maxlambda1}/\code{maxlambda1} are equidistant on a linear scale, if \code{TRUE}
on a logaritmic scale. Please note the different default between \code{optL1} (\code{FALSE}) 
and \code{optL2} (\code{TRUE}).}

\item{tol}{The tolerance of the Brent algorithm used for minimization.
See also \code{\link[stats:optimize]{optimize}}.}
}


\details{All five functions return a list with the following named elements:
  \describe{
    \item{\code{lambda}:}{For \code{optL1} and
      \code{optL2} \code{lambda} gives the optimal value of the tuning parameters found. For
      \code{profL1} and \code{profL2} \code{lambda} is the vector of values of the
      tuning parameter for which the cross-validated likelihood has been calculated.
      Absent in the output of \code{cvl}.}
    \item{\code{cvl}:}{The value(s) of the cross-validated likelihood. For \code{optL1},
      \code{optL2} this is the cross-validated likelihood at the optimal value of the
      tuning parameter.}
    \item{\code{fold}:}{Returns the precise allocation of the subjects into
      the cross-validation folds. Note that the same allocation is used for all
      cross-validated likelihood calculations in each call to \code{optL1}, \code{optL2},
      \code{profL1}, \code{profL2}.}
    \item{\code{predictions}:}{The cross-validated predictions for the left-out samples.
      The precise format of the cross-validated predictions depends on the type of generalized
      linear model (see \code{\link{breslow}} for survival models. The functions \code{profL1}
      and \code{profL2} return a list here, whereas \code{optL1}, \code{optL2} return the
      predictions for the optimal value of the tuning parameter only.}
    \item{\code{fullfit}:}{The fitted model on the full data. The functions \code{profL1}
      and \code{profL2} return a list here, whereas \code{optL1}, \code{optL2} return the
      full data fit for the optimal value of the tuning parameter only.}
  }
}

\value{A named list. See details.}

\note{The \code{optL1} and \code{optL2} functions use Brent's algorithm for 
  minimization without derivatives (see also \code{\link[stats:optimize]{optimize}}).
  There is a risk that these functions converge to a local instead of to a global 
  optimum. This is especially the case for \code{optL1}, as the cross-validated 
  likelihood as a function of \code{lambda1} sometimes has local optima. It is 
  recommended to use \code{optL1} in combination with \code{profL1} to check whether
  \code{optL1} has converged to the right optimum.
  
  See also the notes under \code{\link{penalized}}.
} 

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{penalized}}, \code{\link{penfit}}, \code{\link{plotpath}}.}

\examples{

data(nki70)

# Finding an optimal crossvalidated likelihood
attach(nki70)
opt <- optL1(Surv(time, event), penalized = nki70[,8:77], fold = 10)
coefficients(opt$fullfit)
plot(opt$predictions)

# Plotting the profile of the crossvalidated likelihood
prof <- profL1(Surv(time, event), penalized = nki70[,8:77], fold = opt$fold, steps=20)
plot(prof$lambda, prof$cvl, type="l")
plotpath(prof$fullfit)
}

\keyword{multivariate} \keyword{regression}