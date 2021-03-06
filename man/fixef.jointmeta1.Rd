% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixef.jointmeta1.R
\name{fixef.jointmeta1}
\alias{fixef.jointmeta1}
\title{Extract fixed effects}
\usage{
\method{fixef}{jointmeta1}(object, type = c("Longitudinal", "Survival", "Latent"), ...)
}
\arguments{
\item{object}{A joint model fit from the function jointmeta1}

\item{type}{Type of fixed effects to extract.  To extract fixed effects from
longitudinal sub-model set \code{type = 'Longitudinal'}, from the
time-to-event sub-model set \code{type = 'Survival'}, or extract the latent
association parameters set to \code{type = 'Latent'}.}

\item{...}{additional arguments; currently none are used.}
}
\value{
The function returns a vector of the fixed effects from the specified
  part of the supplied jointmeta1 model fit.
}
\description{
Function to extract the estimated fixed effects from a jointmeta1 model fit
}
\examples{
   #change example data to jointdata object
   jointdat2<-tojointdata(longitudinal = simdat2$longitudinal,
   survival = simdat2$survival, id = 'id',longoutcome = 'Y',
   timevarying = c('time','ltime'),
   survtime = 'survtime', cens = 'cens',time = 'time')

   #set variables to factors
   jointdat2$baseline$study <- as.factor(jointdat2$baseline$study)
   jointdat2$baseline$treat <- as.factor(jointdat2$baseline$treat)

   #fit multi-study joint model
   #note: for demonstration purposes only - max.it restricted to 2
   #model would need more iterations to truely converge
   onestagefit<-jointmeta1(data = jointdat2, long.formula = Y ~ 1 + time +
                           + treat + study, long.rand.ind = c('int'),
                           long.rand.stud = c('treat'),
                           sharingstrct = 'randprop',
                           surv.formula = Surv(survtime, cens) ~ treat,
                           study.name = 'study', strat = TRUE, max.it = 2)

    #extract longitudinal fixed effects
    fixef(onestagefit, type = 'Longitudinal')
}
\seealso{
\code{\link{jointmeta1}}
}
