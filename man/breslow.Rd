\name{breslow object}
\alias{breslow-class}
\alias{breslow}
\docType{class}

\alias{show,breslow-method}
\alias{plot,breslow-method}
\alias{[,breslow-method}
\alias{[[,breslow-method}
\alias{as.list,breslow-method}
\alias{as.matrix,breslow-method}
\alias{as.data.frame,breslow-method}
\alias{time,breslow-method}
\alias{survival}
\alias{survival,breslow-method}


\title{Breslow estimator object}

\description{Stores one or more survival stepfunctions of the Kaplan-Meier and Breslow type.}

\section{Slots}{
  \describe{
    \item{\code{time}:}{Object of class "vector". The time points at which the function starts, jumps or ends.}
    \item{\code{curves}:}{Object of class "matrix". The values of the curve at each time point.
      Note that the function is left-continuous. Curves are stored as a matrix of dimension (\# curves) x (\# timepoints).}
  }
}

\section{Methods}{
  \describe{
    \item{"["}{(breslow): Coerces the object to a matrix (as as.matrix) and selects rows and/or columns.}
    \item{"[["}{(breslow): Selects a subset of the curves.}
    \item{as.data.frame}{(breslow): Returns the "curves" slot together with the "time" slot in a \code{data.frame}.}
    \item{as.list}{(breslow): Coerces the object to a list.}
    \item{as.matrix}{(breslow): Returns the "curves" slot with the "time" slot as column names.}
    \item{plot}{(breslow): Straightforward plotting (all curves in one plot).}
    \item{show}{(breslow): Summarizes the dimensions of the object.}
    \item{survival}{(breslow): Takes a second argument (\code{time}). Extracts the value of the survival curve(s) at any time point.}
  }
}

\details{Breslow objects are generally created by the \code{\link{penalized}} function and related functions.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{penalized}}.}

\keyword{survival}
