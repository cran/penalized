\name{Plotting the LASSO path}
\alias{plotpath}

\title{Plotting the LASSO path}

\description{Plotting a LASSO path fitted using \code{\link{penalized}} with \code{steps > 1}.}

\usage{plotpath(object, labelsize = 0.6,...) }

\arguments{
\item{object}{A \code{list} of \code{link{penfit}} objects calculated for the same data but different \code{lambda1} values. This object can be created
using the \code{link{penalized}} function using the \code{steps} argument.}

\item{labelsize}{Sets the size of the variable labels in the plot. Set to zero for no variable labels.}

\item{...}{Any other arguments will be forwarded to the plot function.}
}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{penalized}}, \code{\link{penfit}}.}

\examples{

data(nki70)

# Fit the model with the steps argument and plot
pen <- penalized(Surv(time, event), penalized = nki70[,8:77], data = nki70, lambda1=1, steps = 100)

plotpath(pen)

}

\keyword{multivariate} \keyword{regression}