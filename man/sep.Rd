% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sep.R
\name{sep}
\alias{sep}
\title{Function to determine if estimates from separate models were requested}
\usage{
sep(ests, logical)
}
\arguments{
\item{ests}{estimates from initial longitudinal or survival analyses}

\item{logical}{a logical value indicating whether or not results from
separate longitudinal and survival analyses were requested.}
}
\value{
A list of results from the separate longitudinal and survival fits.
  The components of this list are: \describe{

  \item{\code{longests}}{a list containing estimates from the initial
  longitudinal fit.  The components of this list are: \describe{

  \item{\code{beta1}}{a data frame of the estimates of the fixed effects from
  the longitudinal sub-model}

  \item{\code{sigma.e}}{the value of the variance of the measurement error
  from the longitudinal sub-model}

  \item{\code{D}}{the estimate of the covariance matrix for the individual
  level random effects.  Individual level random effects are always included
  in the joint model}

  \item{\code{A}}{the estimate of the covariance matrix for the study level
  random effects.  This is only present if study level random effects are
  specified in the \code{jointmeta1} function call.}

  \item{\code{log.like.long}}{the numeric value of the log likelihood for the
  initial longitudinal model.}

  \item{\code{randstart.ind}}{a list of the conditional modes of the
  individual level random effects in each study given the data and the
  estimates of the separate longitudinal model parameters}

  \item{\code{randstart.ind.cov}}{a list of the conditional covariance
  matrices for each individual for the individual level random effects given
  the data and the estimates of the separate longitudinal model parameters}

  \item{\code{randstart.stud}}{a data frame containing the conditional modes
  of the study level random effects given the data and the estimates of the
  separate longitudinal model parameters.  This is only present if study
  level random effects were specified in the \code{jointmeta1} function call.
  }

  \item{\code{randstart.stud.cov}}{a list of conditional covariance matrices
  for each study for the study level random effects given the data and the
  estimates of the separate longitudinal model parameters. This is only
  present if study level random effects were specified in the
  \code{jointmeta1} function call.}

  \item{\code{modelfit}}{the initial longitudinal model fit.  The model has
  the same specification as the longitudinal sub-model for the joint model,
  fitted using the \code{\link[lme4]{lmer}} function from package
  \code{lme4}}

  }}

  \item{\code{survests}}{a list containing estimates from the initial
  survival fit.  The components of this list are: \describe{

  \item{\code{beta2}}{vector of the estimates of the fixed effects included
  in the survival model.}

  \item{\code{haz}}{if \code{strat = TRUE} then this is a list of numeric
  vectors of length equal to the number of studies in the dataset, giving the
  study specific baseline hazard.  If \code{strat = FALSE} then the baseline
  is not stratified by study, and this is one numeric vector giving the
  common baseline across studies.}

  \item{\code{rs}}{a counter to indicate the last how many unique event times
  had occured by the individual's survival time - this is for use during
  further calculation in the joint model EM algorithm.  If a stratified
  baseline this is a list of numerical vectors, whereas if the baseline is
  not stratified this is a single numeric vector.}

  \item{\code{sf}}{the unique event times observed in the dataset. If a
  stratified baseline this is a list of numerical vectors, whereas if the
  baseline is not stratified this is a single numeric vector. }

  \item{\code{nev}}{a counter of the number of events that occur at each
  event time.If a stratified baseline this is a list of numerical vectors,
  whereas if the baseline is not stratified this is a single numeric vector.}

  \item{\code{log.like.surv}}{a numeric containing two values, the
  log-likelihood with the initial values and the log-likelihood with the
  final values, see \code{\link[survival]{coxph.object}}}

  \item{\code{modelfit}}{the initial survival model fit.  The model has the
  same specification as the survival sub-model for the joint model, fitted
  using the \code{\link[survival]{coxph}} function from package
  \code{survival}}

  }}

  }
}
\description{
This function assesses during the \code{\link{jointmeta1}} fit whether
results from separate longitudinal and time-to-event models were requested,
and supplies their results if they were.
}
\keyword{internal}
