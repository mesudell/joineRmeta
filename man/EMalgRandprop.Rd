% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMalgRandprop.R
\name{EMalgRandprop}
\alias{EMalgRandprop}
\title{EM algorithm function used in jointmeta1}
\usage{
EMalgRandprop(
  data,
  longdat,
  survdat,
  long.rand.ind,
  long.rand.stud = NULL,
  id.name,
  study.name,
  gpt,
  max.it,
  tol,
  time.long,
  surv.formula,
  long.formula,
  long.formula.orig,
  paraests,
  studies,
  p1,
  p2,
  strat,
  print.detail,
  bootrun = FALSE,
  q,
  r = NULL
)
}
\arguments{
\item{data}{the original \code{jointdata} as supplied to the
\code{jointmeta1}] call}

\item{longdat}{the longitudinal data with factors and interaction terms
expanded, ordered by increasing survival time}

\item{survdat}{the survival data with factors and interaction terms expanded,
ordered by increasing survival time}

\item{long.rand.ind}{a vector of character strings to indicate what variables
to assign individual level random effects to.  A maximum of three
individual level random effects can be assigned.  To assign a random
intercept include 'int' in the vector.  To not include an individual level
random intercept include 'noint' in the vector.  For example to fit a model
with individual level random intercept and random slope set
\code{long.rand.ind = c('int', 'time')}, where \code{'time'} is the
longitudinal time variable in the \code{data}.}

\item{long.rand.stud}{a vector of character strings to indicate what
variables to assign study level random effects to.  If no study level
random effects then this either not specified in function call or set to
\code{NULL}.  If a study level random intercept is required, include the
name of the study membership variable for example \code{long.rand.stud =
'study'}.}

\item{id.name}{character string specifying the id variable in the dataset}

\item{study.name}{a character string denoting the name of the variable in the
baseline dataset in \code{data} holding study membership, for example
\code{study.name = 'study'}.}

\item{gpt}{the number of quadrature points across which the integration with
respect to the random effects will be performed.  If random effects are
specified at both the individual and the study level, the same number of
quadrature points is used in both cases.  Defaults to \code{gpt = 5}.}

\item{max.it}{the maximum number of iterations of the EM algorithm that the
function will perform.  Defaults to \code{max.it = 350} although more
iterations could be required for large complex datasets.}

\item{tol}{the tolerance level used to determine convergence in the EM
algorithm.  Defaults to \code{tol = 0.001}.}

\item{time.long}{the name of the variable holding the longitudinal time
covariate}

\item{surv.formula}{a formula object with the survival time, censoring
indicator and the covariates to include in the survival sub-model.  The
response must be a survival object as returned by the
\code{\link[survival]{Surv}} function.}

\item{long.formula}{a formula object with the response varaible, and the
covariates to include in the longitudinal sub-model}

\item{long.formula.orig}{the original longitudinal formula as supplied to the
function call}

\item{paraests}{a list of the estimates present from the separate
longitudinal and survival fits.  Same structure as \code{sepests} if
requested in a \code{\link{jointmeta1.object}}}

\item{studies}{the names of the studies present in the supplied data}

\item{p1}{the number of fixed effects included in the longitudinal sub-model}

\item{p2}{the number of fixed effects included in the survival sub-model}

\item{strat}{logical value: if \code{TRUE} then the survival sub-model is
calculated with a baseline stratified by study.  Otherwise baseline is
unstratified}

\item{print.detail}{logical value: if \code{TRUE} then details of the
parameter estimates at each iteration of the EM algorithm are printed to
the console.}

\item{bootrun}{logical value: if \code{TRUE} then the log-likelihood for the
model is not calculated.  This option is available so that when
bootstrapping to obtain standard errors, as the log-likelihood is not
needed, it is not calculated, thus speeding up the bootstrapping process.}

\item{q}{the number of individual level random effects}

\item{r}{the number of study level random effects, set to \code{NULL} if no
study level random effects included in the model}
}
\value{
This function returns a list of the estimates of parameters and other
  information from the run of the EM algorithm.  The list has the following
  components: \describe{

  \item{\code{beta1}}{a data frame containing the estimates of the fixed
  effect parameters from the longitudinal sub-model.}

  \item{\code{beta2}}{a data frame containing the estimates of the fixed
  effect parameters from the survival sub-model.}

  \item{\code{sigma.e}}{the estimate of the variance of the measurement
  errors.}

  \item{\code{haz}}{the estimated baseline hazard.  If \code{strat = TRUE} in
  the function call to \code{jointmeta1} then this is a list of length equal
  to the number of studies in the supplied dataset, each element of the list
  being the baseline hazard for the corresponding study. Otherwise there is a
  common baseline across all studies in the dataset and this is one vector.}

  \item{\code{D}}{the estimated covariance matrix for the individual level
  random effects}

  \item{\code{A}}{the estimated covariance matrix for the study level random
  effects.  This is only present in the output if study level random effects
  were specified in the function call to \code{jointmeta1}.}

  \item{\code{random2}}{a list of matrices containing the conditional modes
  of the individual level random effects given the supplied data and the
  estimated parameters of the joint model. The list is of length equal to the
  number of studies in the dataset, and each element of the list has number
  of rows equal to the number of individuals in the study, and number of
  columns equal to the number of specified individual level random effects.}

  \item{\code{random3}}{a matrix containing the conditional modes of the
  study level random effects given the supplied data and the estimated
  parameters of the joint model.  The matrix has number of rows equal to the
  number of studies, and number of columns equal to the number of specified
  study level random effects.}

  \item{\code{long.rand.ind.form}}{a character string giving the formulation
  of the individual level random effects.}

  \item{\code{long.rand.stud.form}}{a character string giving the formulation
  of the study level random effects if included in the model.}

  \item{\code{conv}}{a logical value indicating whether convergence of the EM
  algorithm was achieved or not.}

  \item{\code{iters}}{the number of iterations completed by the EM algorithm}

  \item{\code{n.bystudy}}{the number of individuals present in each study in
  the data supplied to the function.} }
}
\description{
Function to run EM algorithm during one stage model fit.  Used when the
jointmeta1 function is called.
}
\seealso{
\code{\link{jointmeta1}}, \code{\link{tojointdata}},
  \code{\link[joineR]{jointdata}},\code{\link[lme4]{lmer}},
  \code{\link[survival]{coxph}}
}
\keyword{internal}
