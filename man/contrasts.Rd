\name{Penalized regression contrasts}
\alias{contr.none}
\alias{contr.diff}

\title{Contrast functions for penalized regression}

\description{Contrast functions for factors that are appropriate for penalized regression.}

\usage{

contr.none(n, contrasts)

contr.diff(n, contrasts = TRUE)
}

\arguments{
\item{n}{A vector of levels for a factor, or the number of levels.}

\item{contrasts}{A logical indicating whether contrasts should be computed. This
  argument is ignored in \code{contr.none}.}
}

\details{These functions are used for creating contrast matrices for use in fitting penalized
analysis of variance and regression models. The columns of the resulting matrices
contain contrasts which can be used for coding a factor with n levels. The returned
value contains the computed contrasts.

\code{contr.none} returns an identity matrix except when the number of levels is 2,
in which case it returns a single contrast. \code{contr.none} ensures that all levels
of an unordered factor are treated symmetrically in a penalized regression model.

\code{contr.diff} returns a lower triangular matrix of ones if \code{contrasts=FALSE}
and the same matrix without its first column if \code{contrasts=TRUE}. This makes sure
that penalization is done on the difference between successive levels of an ordered factor.
It is not appropriate for unordered factors.
}

\value{\code{contr.diff} returns a matrix with \code{n} rows and \code{k} columns, with
\code{k=n-1} if contrasts is \code{TRUE} and \code{k=n} if contrasts is \code{FALSE}.

\code{contr.none} returns a matrix with \code{n} rows and \code{n} columns, except when \code{n=2}
when it returns a matrix with 2 rows and one column.}

\author{Jelle Goeman: \email{j.j.goeman@lumc.nl}}

\seealso{\code{\link{penalized}}, \code{\link{contr.treatment}}, \code{\link{contr.poly}}, \code{\link{contr.helmert}},
\code{\link{contr.SAS}}, \code{\link{contr.sum}}.}

\examples{

# Three levels
levels <- LETTERS[1:3]
contr.none(levels)
contr.diff(levels)

# Two levels
levels <- LETTERS[1:2]
contr.none(levels)
}

\keyword{regression}