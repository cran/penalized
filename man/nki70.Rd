\name{Netherlands Cancer Institute 70 gene signature}

\docType{data}

\alias{nki70}

\title{The 70 gene signature of the Netherlands Cancer Institute for prediction of metastasis-free survival, measured on 144 independent lymph node positive patients.}

\description{A \code{\link{data.frame}}
of 144 lymph node positive breast cancer patients on metastasis-free survival, 5 clinical risk factors, and gene expression measurements of 70 genes found to be prognostic for metastasis-free survival in an earlier study. The included variables are
  \describe{
  \item{time}{Metastasis-free follow-up time (months).}
  \item{event}{Event indicator. 1 = metastasis or death; 0 = censored.}
  \item{Diam}{Diameter of the tumor (two levels).}
  \item{N}{Number of affected lymph nodes (two levels).}
  \item{ER}{Estrogen receptor status (two levels).}
  \item{Grade}{Grade of the tumor (three ordered levels).}
  \item{Age}{Patient age at diagnosis (years).}
  \item{TSPYL5...C20orf46}{Gene expression measurements of 70 prognostic genes.}
  }
}

\usage{data(nki70)}

\format{A \code{\link{data.frame}}.}

\note{The 70 gene signature was trained on lymph node negative patients only.}

\references{M.J. van de Vijver, Y.D. He, L.J. van 't Veer, H. Dai, A.A.M. Hart, D.W. Voskuil, G.J. Schreiber, J.L. Peterse, C. Roberts, M.J. Marton, M. Parrish, D. Atsma, A. Witteveen, A. Glas, L. Delahaye, T. van der Velde, H. Bartelink, S. Rodenhuis, E.T. Rutgers, S.H. Friend, and R. Bernards (2002). A gene-expression signature as a predictor of survival in breast cancer. New England Journal of Medicine 347 (25), 1999--2009.}

\keyword{datasets}
