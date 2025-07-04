% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jointmetaSE.object.R
\docType{data}
\name{jointmeta1SE.object}
\alias{jointmeta1SE.object}
\title{A \code{jointmeta1SE} object}
\format{
An object of class \code{NULL} of length 0.
}
\usage{
jointmeta1SE.object
}
\value{
A list with the following components. \describe{

  \item{\code{results}}{a data frame containing the estimates, standard
  errors and 95\% confidence intervals for the parameters from the model and
  any overall effects requested.}

  \item{\code{covmat}}{the covariance matrix for the model parameters}

  \item{\code{bootstraps}}{a data frame containing the results of each
  bootstrap}

  }
}
\description{
An object returned by the \code{jointmetaSE} function, inheriting from class
\code{jointmeta1SE} representing the results of bootstrapping a fit from the
\code{jointmeta1} function.  Objects of this class have methods for the
\code{\link{print.jointmeta1SE}} and \code{\link{vcov.jointmeta1SE}} functions.
}
\seealso{
\code{\link{jointmeta1}}, \code{\link{jointmetaSE}}.
}
\author{
Maria Sudell (\email{mesudell@liverpool.ac.uk})
}
\keyword{datasets}
