#' Simulated joint longitudinal and survival dataset containing 5 studies
#'
#' A simulated dataset containing a single continuous longitudinal outcome and a
#' single survival outcome, with data available from 5 studies.
#'
#' @format A list of three objects: \describe{ \item{\code{longitudinal}}{A list
#'   of long format longitudinal datasets one for each of the 5 studies included
#'   in the dataset.  Each of these datasets contains the following variables:
#'   \describe{ \item{\code{id}}{long version of the id variable for the data.
#'   Identical ids between the longitudinal and the survival datasets identify
#'   the same individual} \item{\code{Y}}{a continuous longitudinal outcome}
#'   \item{\code{time}}{the longitudinal time variable} \item{\code{study}}{a
#'   long version of the study membership indicator} \item{\code{intercept}}{a
#'   long version of the intercept, always takes a value of 1}
#'   \item{\code{treat}}{a long version of the binary treatment group indicator}
#'   \item{\code{ltime}}{a duplicate of the longitudinal time variable,
#'   duplicated as part of the longitudinal data simulation process} }}
#'
#'   \item{\code{survival}}{A list of survival datasets, one for each of the 5
#'   studies included in the dataset.  Each of these datasets contains the
#'   following variables: \describe{ \item{\code{id}}{the id variable for the
#'   data. Identical ids between the longitudinal and the survival datasets
#'   identify the same individual} \item{\code{survtime}}{the survival time for
#'   each individual at which they experienced the event or were censored.  This
#'   is on the same scale as the longitudinal time measurements.}
#'   \item{\code{cens}}{censoring indicator for the survival data where 1
#'   indicates an event and 0 indicates censoring} \item{\code{study}}{study
#'   membership indicator} \item{\code{treat}}{binary treatment group indicator}
#'   }}
#'
#'   \item{\code{percentevent}}{A list of the percentage of events experienced
#'   in each datasets.  The first element contains the percentage of events
#'   observed for the first simulated study and so on.} }
#'
#' @details This is a simulated dataset generated using the
#'   \code{\link{simjointmeta}} function using the following function call:
#'
#'   \code{ simdat<-simjointmeta(k = 5, n = rep(500, 5), sepassoc = FALSE, ntms
#'   = 5, longmeasuretimes = c(0, 1, 2, 3, 4), beta1 = c(1, 2, 3), beta2 = 3,
#'   rand_ind = 'intslope', rand_stud = 'inttreat', gamma_ind = 1, gamma_stud =
#'   1, sigb_ind = matrix(c(1,0.5,0.5,1.5),nrow=2), sigb_stud =
#'   matrix(c(1,0.5,0.5,1.5),nrow=2), vare = 0.01, theta0 = -3, theta1 = 1,
#'   censoring = TRUE, censlam = exp(-3), truncation = TRUE, trunctime = 6) }
#'
#'   Note that this will not give you identical data to that held in
#'   \code{simdat} due to the differences in starting seed.
#'
#' @seealso \code{\link{simjointmeta}}
#'
"simdat"

#' Simulated joint longitudinal and survival dataset containing 15 studies
#'
#' A simulated dataset containing a single continuous longitudinal outcome and a
#' single survival outcome, with data available from 15 studies.
#'
#' @format A list of three objects: \describe{ \item{\code{longitudinal}}{A list
#'   of long format longitudinal datasets one for each of the 15 studies
#'   included in the dataset.  Each of these datasets contains the following
#'   variables: \describe{ \item{\code{id}}{long version of the id variable for
#'   the data. Identical ids between the longitudinal and the survival datasets
#'   identify the same individual} \item{\code{Y}}{a continuous longitudinal
#'   outcome} \item{\code{time}}{the longitudinal time variable}
#'   \item{\code{study}}{a long version of the study membership indicator}
#'   \item{\code{intercept}}{a long version of the intercept, always takes a
#'   value of 1} \item{\code{treat}}{a long version of the binary treatment
#'   group indicator} \item{\code{ltime}}{a duplicate of the longitudinal time
#'   variable, duplicated as part of the longitudinal data simulation process}
#'   }}
#'
#'   \item{\code{survival}}{A list of survival datasets, one for each of the 15
#'   studies included in the dataset.  Each of these datasets contains the
#'   following variables: \describe{ \item{\code{id}}{the id variable for the
#'   data. Identical ids between the longitudinal and the survival datasets
#'   identify the same individual} \item{\code{survtime}}{the survival time for
#'   each individual at which they experienced the event or were censored.  This
#'   is on the same scale as the longitudinal time measurements.}
#'   \item{\code{cens}}{censoring indicator for the survival data where 1
#'   indicates an event and 0 indicates censoring} \item{\code{study}}{study
#'   membership indicator} \item{\code{treat}}{binary treatment group indicator}
#'   }}
#'
#'   \item{\code{percentevent}}{A list of the percentage of events experienced
#'   in each datasets.  The first element contains the percentage of events
#'   observed for the first simulated study and so on.} }
#'
#' @details This is a simulated dataset generated using the
#'   \code{\link{simjointmeta}} function using the following function call:
#'
#'   \code{ simdat2<-simjointmeta(k = 15, n = rep(500, 15), sepassoc = FALSE,
#'   ntms = 5, longmeasuretimes = c(0, 1, 2, 3, 4), beta1 = c(1, 2, 3), beta2 =
#'   3, rand_ind = 'intslope', rand_stud = 'inttreat', gamma_ind = 1, gamma_stud
#'   = 1, sigb_ind = matrix(c(1,0.5,0.5,1.5),nrow=2), sigb_stud =
#'   matrix(c(1,0.5,0.5,1.5),nrow=2), vare = 0.01, theta0 = -3, theta1 = 1,
#'   censoring = TRUE, censlam = exp(-3), truncation = TRUE, trunctime = 6) }
#'
#'   Note that this will not give you identical data to that held in
#'   \code{simdat2} due to the differences in starting seed.
#'
#' @seealso \code{\link{simjointmeta}}
#'
"simdat2"

#' Simulated joint longitudinal and survival dataset containing 5 studies
#'
#' A simulated dataset containing a single continuous longitudinal outcome and a
#' single survival outcome, with data available from 5 studies.  This dataset
#' does not have longitudinal measurements capped at each individual's survival
#' time.
#'
#' @format A list of three objects: \describe{ \item{\code{longitudinal}}{A list
#'   of long format longitudinal datasets one for each of the 5 studies included
#'   in the dataset.  Each of these datasets contains the following variables:
#'   \describe{ \item{\code{id}}{long version of the id variable for the data.
#'   Identical ids between the longitudinal and the survival datasets identify
#'   the same individual} \item{\code{Y}}{a continuous longitudinal outcome}
#'   \item{\code{time}}{the longitudinal time variable} \item{\code{study}}{a
#'   long version of the study membership indicator} \item{\code{intercept}}{a
#'   long version of the intercept, always takes a value of 1}
#'   \item{\code{treat}}{a long version of the binary treatment group indicator}
#'   \item{\code{ltime}}{a duplicate of the longitudinal time variable,
#'   duplicated as part of the longitudinal data simulation process} }}
#'
#'   \item{\code{survival}}{A list of survival datasets, one for each of the 5
#'   studies included in the dataset.  Each of these datasets contains the
#'   following variables: \describe{ \item{\code{id}}{the id variable for the
#'   data. Identical ids between the longitudinal and the survival datasets
#'   identify the same individual} \item{\code{survtime}}{the survival time for
#'   each individual at which they experienced the event or were censored.  This
#'   is on the same scale as the longitudinal time measurements.}
#'   \item{\code{cens}}{censoring indicator for the survival data where 1
#'   indicates an event and 0 indicates censoring} \item{\code{study}}{study
#'   membership indicator} \item{\code{treat}}{binary treatment group indicator}
#'   }}
#'
#'   \item{\code{percentevent}}{A list of the percentage of events experienced
#'   in each datasets.  The first element contains the percentage of events
#'   observed for the first simulated study and so on.} }
#'
#' @details This is a simulated dataset generated in order to show how the
#'   function \code{\link{removeafter}} works, as it does not cap longitudinal
#'   outcomes at each individual's survival time.  This data was generated by
#'   manually stepping through the code available in simjointmeta and retaining
#'   instead of discarding any simulated longitudinal measurements recorded
#'   after the individual in question's survival time.
#'
#' @seealso \code{\link{simjointmeta}}, \code{\link{removeafter}}
#'
"simdat3"


#' Study specific joint model fits using the joineR package
#'
#' A dataset containing a list of the model fits for joint models fitted to the
#' data for each study in the \code{simdat} dataset using the joineR package.
#' Further details of model fits supplied below.
#'
#' @format A list of 10 objects: \describe{ \item{\code{joineRfit1}}{an object
#'   of class \code{joint}, the result of using the \code{joint} function to fit
#'   a joint model to the data from the first study in the \code{simdat}
#'   dataset.} \item{\code{joineRfit1SE}}{an object of class \code{data.frame},
#'   the result of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit1}.} \item{\code{joineRfit2}}{ an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the second study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit2SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit2}.} \item{\code{joineRfit3}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the third study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit3SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit3}.} \item{\code{joineRfit4}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the fourth study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit4SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit4}.} \item{\code{joineRfit5}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the fifth study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit5SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit5}.} }
#'
#' @details These are the results of fitting a joint model using the
#'   \code{joineR} package separately to the data from each study present in the
#'   \code{simdat} dataset.  This data has three levels, namely the longitudinal
#'   measurements at level 1, nested within individuals (level 2) who are
#'   themselves nested within studies (level 3). The joint models fitted to each
#'   study's data had the same format.  The longitudinal sub-model contained a
#'   fixed intercept, time and treatment assignment term, and random intercept
#'   and slope.  The survival sub-model contained a fixed treatment assignment
#'   term.  A proportional association structure was used to link the
#'   sub-models.  More formally, the longitudinal sub-model had the following
#'   format:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   b^{(2)}_{0ki} + b^{(2)}_{1ki}time + \epsilon_{kij}}
#'
#'   Where \eqn{Y} represents the continuous longitudinal outcome, fixed effect
#'   coefficients are represented by \eqn{\beta}, random effects coefficients by
#'   \eqn{b} and the measurement error by \eqn{\epsilon}.  For the  random
#'   effects the superscript of 2 indicates that these are individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. The longitudinal time variable is
#'   represented by \eqn{time}, and the treatment assignment variable (a binary
#'   factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha(b^{(2)}_{0ki} + b^{(2)}_{1ki}time)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the baseline hazard.  The fixed effect coefficient is
#'   represented by \eqn{\beta_{21}}.  The association parameter quantifying the
#'   link between the sub-models is represented by \eqn{\alpha}.  Again
#'   \eqn{treat} represents the binary factor treatment assignment variable, and
#'   \eqn{b^{(2)}_{0ki}} and \eqn{b^{(2)}_{1ki}} are the zero mean random
#'   effects shared from the longitudinal sub-model.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#' @seealso \code{\link[joineR]{jointdata}}, \code{\link[joineR]{joint}},
#'   \code{\link[joineR]{jointSE}}
#'
"joineRfits"

#' Study specific joint model fits using the joineR package
#'
#' A dataset containing a list of the model fits for joint models fitted to the
#' data for each study in the \code{simdat} dataset using the joineR package.
#' Further details of model fits supplied below.
#'
#' @format A list of 10 objects: \describe{ \item{\code{joineRfit6}}{an object
#'   of class \code{joint}, the result of using the \code{joint} function to fit
#'   a joint model to the data from the first study in the \code{simdat}
#'   dataset.} \item{\code{joineRfit6SE}}{an object of class \code{data.frame},
#'   the result of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit6}.} \item{\code{joineRfit7}}{ an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the second study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit7SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit7}.} \item{\code{joineRfit8}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the third study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit8SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit8}.} \item{\code{joineRfit9}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the fourth study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit9SE}}{an object of class \code{data.frame}, the result
#'   of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit9}.} \item{\code{joineRfit10}}{an object of class
#'   \code{joint}, the result of using the \code{joint} function to fit a joint
#'   model to the data from the fifth study in the \code{simdat} dataset.}
#'   \item{\code{joineRfit10SE}}{an object of class \code{data.frame}, the
#'   result of applying the function \code{jointSE} to the joint model fit
#'   \code{joineRfit10}.} }
#'
#' @details These are the results of fitting a joint model using the
#'   \code{joineR} package separately to the data from each study present in the
#'   \code{simdat} dataset.  This data has three levels, namely the longitudinal
#'   measurements at level 1, nested within individuals (level 2) who are
#'   themselves nested within studies (level 3). The joint models fitted to each
#'   study's data had the same format.  The longitudinal sub-model contained a
#'   fixed intercept, time and treatment assignment term, and random intercept.
#'   The survival sub-model contained a fixed treatment assignment term.  A
#'   proportional association structure was used to link the sub-models.  More
#'   formally, the longitudinal sub-model had the following format:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   b^{(2)}_{0ki} + \epsilon_{kij}}
#'
#'   Where \eqn{Y} represents the continuous longitudinal outcome, fixed effect
#'   coefficients are represented by \eqn{\beta}, random effects coefficients by
#'   \eqn{b} and the measurement error by \eqn{\epsilon}.  For the  random
#'   effects the superscript of 2 indicates that these are individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. The longitudinal time variable is
#'   represented by \eqn{time}, and the treatment assignment variable (a binary
#'   factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha(b^{(2)}_{0ki})) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the baseline hazard.  The fixed effect coefficient is
#'   represented by \eqn{\beta_{21}}.  The association parameter quantifying the
#'   link between the sub-models is represented by \eqn{\alpha}.  Again
#'   \eqn{treat} represents the binary factor treatment assignment variable, and
#'   \eqn{b^{(2)}_{0ki}} are the zero mean random effects shared from the
#'   longitudinal sub-model.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#' @seealso \code{\link[joineR]{jointdata}}, \code{\link[joineR]{joint}},
#'   \code{\link[joineR]{jointSE}}
#'
"joineRfits2"

#' Study specific joint model fits using the JM package
#'
#' A dataset containing a list of the model fits for joint models fitted to the
#' data for each study in the \code{simdat} dataset using the JM package.
#' Further details of model fits supplied below.
#'
#' @format A list of 5 \code{jointModel} objects, the result of fitting a joint
#'   model using the JM package to the data from each study in the \code{simDat}
#'   dataset in turn.
#'
#' @details These are the results of fitting a joint model using the \code{JM}
#'   package separately to the data from each study present in the \code{simdat}
#'   dataset.  This data has three levels, namely the longitudinal measurements
#'   at level 1, nested within individuals (level 2) who are themselves nested
#'   within studies (level 3). The joint models fitted to each study's data had
#'   the same format.  The longitudinal sub-model contained a fixed intercept,
#'   time and treatment assignment term, and random intercept and slope.  The
#'   survival sub-model contained a fixed treatment assignment term.  A current
#'   value association structure was used to link the sub-models.  More
#'   formally, the longitudinal sub-model had the following format:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   b^{(2)}_{0ki} + b^{(2)}_{1ki}time + \epsilon_{kij}}
#'
#'   Where \eqn{Y} represents the continuous longitudinal outcome, fixed effect
#'   coefficients are represented by \eqn{\beta}, random effects coefficients by
#'   \eqn{b} and the measurement error by \eqn{\epsilon}.  For the  random
#'   effects the superscript of 2 indicates that these are individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. The longitudinal time variable is
#'   represented by \eqn{time}, and the treatment assignment variable (a binary
#'   factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha(\beta_{10} + \beta_{11}time + \beta_{12}treat + b^{(2)}_{0ki} +
#'   b^{(2)}_{1ki}time)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the baseline hazard, which was modelled using splines. The fixed
#'   effect coefficient is represented by \eqn{\beta_{21}}.  The association
#'   parameter quantifying the link between the sub-models is represented by
#'   \eqn{\alpha}.  Again \eqn{treat} represents the binary factor treatment
#'   assignment variable, and \eqn{b^{(2)}_{0ki}} and \eqn{b^{(2)}_{1ki}} are
#'   the zero mean random effects from the longitudinal sub-model.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#' @seealso \code{\link[JM]{jointModel}}, \code{\link[JM]{jointModelObject}}
#'
"JMfits"

#' Study specific joint model fits using the JM package
#'
#' A dataset containing a list of the model fits for joint models fitted to the
#' data for each study in the \code{simdat} dataset using the JM package.
#' Further details of model fits supplied below.
#'
#' @format A list of 5 \code{jointModel} objects, the result of fitting a joint
#'   model using the JM package to the data from each study in the \code{simDat}
#'   dataset in turn.
#'
#' @details These are the results of fitting a joint model using the \code{JM}
#'   package separately to the data from each study present in the \code{simdat}
#'   dataset.  This data has three levels, namely the longitudinal measurements
#'   at level 1, nested within individuals (level 2) who are themselves nested
#'   within studies (level 3). The joint models fitted to each study's data had
#'   the same format.  The longitudinal sub-model contained a fixed intercept,
#'   time and treatment assignment term, as well as a fixed time by treatment
#'   assignment interaction term, and random intercept and slope.  The survival
#'   sub-model contained a fixed treatment assignment term. The sub-models were
#'   linked by inserting both the current value of the longitudinal trajectory
#'   and its first derivative with respect to time into the survival sub-model.
#'   More formally, the longitudinal sub-model had the following format:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   \beta_{13}time*treat+ b^{(2)}_{0ki} + b^{(2)}_{1ki}time + \epsilon_{kij}}
#'
#'   Where \eqn{Y} represents the continuous longitudinal outcome, fixed effect
#'   coefficients are represented by \eqn{\beta}, random effects coefficients by
#'   \eqn{b} and the measurement error by \eqn{\epsilon}.  For the  random
#'   effects the superscript of 2 indicates that these are individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. The longitudinal time variable is
#'   represented by \eqn{time}, and the treatment assignment variable (a binary
#'   factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha_{1}(\beta_{10} + \beta_{11}time + \beta_{12}treat + +
#'   \beta_{13}time*treat b^{(2)}_{0ki} + b^{(2)}_{1ki}time)+
#'   \alpha_{2}(\beta_{11} + \beta_{13}treat + b^{(2)}_{1ki})) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the baseline hazard, which was modelled using splines. The fixed
#'   effect coefficient is represented by \eqn{\beta_{21}}.  Association
#'   parameters representing the link between the sub-models are represented by
#'   \eqn{\alpha} terms, where \eqn{\alpha_{1}} represents the effect of the
#'   current value of the longitudinal outcome on the risk of an event, whilst
#'   \eqn{\alpha_{2}} represents the effect of the slope, or rate of change of
#'   the longitudinal trajectory (the value of the first derivative of the
#'   longitudinal trajectory with respect to time) on the risk of an event.
#'   Again \eqn{treat} represents the binary facator treatment assignment
#'   variable, and \eqn{b^{(2)}_{0ki}} and \eqn{b^{(2)}_{1ki}} are the zero mean
#'   random effects from the longitudinal sub-model.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#' @seealso \code{\link[JM]{jointModel}}, \code{\link[JM]{jointModelObject}}
#'
"JMfits2"

#' One stage jointmeta1 fit and bootstrapped standard errors
#'
#' A list of length two containing a one stage jointmeta1 fit and corresponding
#' bootstrapped standard errors.
#'
#' @format A list of 2 objects: \describe{ \item{\code{onestagefit0}}{an object
#'   of class \code{jointmeta1}} \item{\code{onestagefit0SE}}{an object of class
#'   \code{jointmeta1SE}} }
#'
#' @details These are the results of using the \code{jointmeta1} function to fit
#'   a one stage joint meta model for multi-study data, and also the bootstrap
#'   results of applying the \code{jointmetaSE} function to the resulting model
#'   fit.  The data used is the \code{simdat} data available in the
#'   \code{joineRmeta} package.  This data has three levels, namely the
#'   longitudinal measurements at level 1, nested within individuals (level 2)
#'   who are themselves nested within studies (level 3).
#'
#'   The format of this model is as follows.  The structure of the longitudinal
#'   sub-model is:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   b^{(2)}_{0ki} + b^{(2)}_{1ki}time + \epsilon_{kij}}
#'
#'   \eqn{Y_{kij}} represents the continuous longitudinal outcome for the
#'   \eqn{i}th individual in the \eqn{k}th study at the \eqn{j}th time point,
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects
#'   coefficients by \eqn{b} and the measurement error by \eqn{\epsilon}.  For
#'   the  random effects the superscript of 2 indicates that these are
#'   individual level, or level 2 random effects.  This means they take can take
#'   a unique value for each individual in the dataset. The longitudinal time
#'   variable is represented by \eqn{time}, and the treatment assignment
#'   variable (a binary factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha^{(2)}(b^{(2)}_{0ki} + b^{(2)}_{1ki}time)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the unspecified baseline hazard.  This baseline was not
#'   stratified by study. The fixed effect coefficient is represented by
#'   \eqn{\beta_{21}}.  A proportional random effects only association structure
#'   links the sub-models, with \eqn{\alpha^{(2)}} representing the association
#'   between the longitudinal and survival outcomes attributable to the
#'   deviation of the individual in question from the population mean
#'   longitudinal trajectory.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   This is a naive model as it analyses data from all the studies in the
#'   dataset but does not account for between study heterogeneity (differences
#'   between the studies included in the dataset) in any way.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#'   The code used to fit this one stage model was:
#'
#'   \code{ onestagefit0<-jointmeta1(data = jointdat, long.formula = Y ~ 1 +
#'   time + treat, long.rand.ind = c('int', 'time'), sharingstrct = 'randprop',
#'   surv.formula = Surv(survtime, cens) ~ treat, study.name = 'study', strat =
#'   F) }
#'
#'   And the code used to bootstrap the model was:
#'
#'   \code{ onestagefit0SE<-jointmetaSE(fitted = onestagefit0, n.boot = 200) }
#'
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}
#'
"onestage0"

#' One stage jointmeta1 fit and bootstrapped standard errors
#'
#' A list of length two containing a one stage jointmeta1 fit and corresponding
#' bootstrapped standard errors.
#'
#' @format A list of 2 objects: \describe{ \item{\code{onestagefit1}}{an object
#'   of class \code{jointmeta1}} \item{\code{onestagefit1SE}}{an object of class
#'   \code{jointmeta1SE}} }
#'
#' @details These are the results of using the \code{jointmeta1} function to fit
#'   a one stage joint meta model for multi-study data, and also the bootstrap
#'   results of applying the \code{jointmetaSE} function to the resulting model
#'   fit. The data used is the \code{simdat} data available in the
#'   \code{joineRmeta} package.  This data has three levels, namely the
#'   longitudinal measurements at level 1, nested within individuals (level 2)
#'   who are themselves nested within studies (level 3).
#'
#'   The format of this model is as follows.  The structure of the longitudinal
#'   sub-model is:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   \beta_{13}study + \beta_{14}treat*study + b^{(2)}_{0ki} + b^{(2)}_{1ki}time
#'   + \epsilon_{kij}}
#'
#'   \eqn{Y_{kij}} represents the continuous longitudinal outcome for the
#'   \eqn{i}th individual in the \eqn{k}th study at the \eqn{j}th time point,
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects
#'   coefficients by \eqn{b} and the measurement error by \eqn{\epsilon}.  For
#'   the  random effects the superscript of 2 indicates that these are
#'   individual level, or level 2 random effects.  This means they take can take
#'   a unique value for each individual in the dataset. The longitudinal time
#'   variable is represented by \eqn{time}, and the treatment assignment
#'   variable (a binary factor) is represented by \eqn{treat}.  The study
#'   membership variable, a factor variable, is represented by \eqn{study}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat + \beta_{22}study
#'   + \beta_{23}treat*study + \alpha^{(2)}(b^{(2)}_{0ki} + b^{(2)}_{1ki}time))
#'   }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the unspecified baseline hazard.  This baseline was not
#'   stratified by study. The fixed effect coefficients are represented by
#'   \eqn{\beta} terms.  A proportional random effects only association
#'   structure links the sub-models, with \eqn{\alpha^{(2)}} representing the
#'   association between the longitudinal and survival outcomes attributable to
#'   the deviation of the individual in question from the population mean
#'   longitudinal trajectory.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   This model accounts for between study heterogeneity by including study
#'   membership and interactions between the study membership and the treatment
#'   assignment in the fixed effects of both sub-models.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#'   The code used to fit this one stage model was:
#'
#'   \code{ onestagefit1<-jointmeta1(data = jointdat, long.formula = Y ~ 1 +
#'   time + treat*study, long.rand.ind = c('int', 'time'), sharingstrct =
#'   'randprop', surv.formula = Surv(survtime, cens) ~ treat*study, study.name =
#'   'study', strat = F) }
#'
#'   And the code used to bootstrap the model was:
#'
#'   \code{ onestagefit1SE<-jointmetaSE(fitted = onestagefit1, n.boot = 200,
#'   overalleffects = list(long = list(c('treat1', 'treat1:study2'), c('treat1',
#'   'treat1:study3'), c('treat1', 'treat1:study4'), c('treat1',
#'   'treat1:study5')), surv = list(c('treat1', 'treat1:study2'), c('treat1',
#'   'treat1:study3'), c('treat1', 'treat1:study4'), c('treat1',
#'   'treat1:study5')))) }
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}
#'
"onestage1"

#' One stage jointmeta1 fit and bootstrapped standard errors
#'
#' A list of length two containing a one stage jointmeta1 fit and corresponding
#' bootstrapped standard errors.
#'
#' @format A list of 2 objects: \describe{ \item{\code{onestagefit2}}{an object
#'   of class \code{jointmeta1}} \item{\code{onestagefit2SE}}{an object of class
#'   \code{jointmeta1SE}} }
#'
#' @details These are the results of using the \code{jointmeta1} function to fit
#'   a one stage joint meta model for multi-study data, and also the bootstrap
#'   results of applying the \code{jointmetaSE} function to the resulting model
#'   fit. The data used is the \code{simdat} data available in the
#'   \code{joineRmeta} package.  This data has three levels, namely the
#'   longitudinal measurements at level 1, nested within individuals (level 2)
#'   who are themselves nested within studies (level 3).
#'
#'   The format of this model is as follows.  The structure of the longitudinal
#'   sub-model is:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   b^{(2)}_{0ki} + b^{(2)}_{1ki}time + b^{(3)}_{0k} + b^{(3)}_{1k}treat +
#'   \epsilon_{kij}}
#'
#'   Where \eqn{Y_{kij}} represents the continuous longitudinal outcome for the
#'   \eqn{i}th individual in the \eqn{k}th study at the \eqn{j}th time point,
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects
#'   coefficients by \eqn{b} and the measurement error by \eqn{\epsilon}.  For
#'   the  random effects the superscript of 2 indicates individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. A superscript of 3 indicates study level
#'   random effects, or level 3 random effects.  This means that they can take a
#'   unique value for each study in the dataset. The longitudinal time variable
#'   is represented by \eqn{time}, and the treatment assignment variable (a
#'   binary factor) is represented by \eqn{treat}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(\beta_{21}treat +
#'   \alpha^{(2)}(b^{(2)}_{0ki} + b^{(2)}_{1ki}time) + \alpha^{(3)}(b^{(2)}_{0k}
#'   + b^{(3)}_{1k}treat)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the unspecified baseline hazard.  This baseline was not
#'   stratified by study. The fixed effect coefficients are represented by
#'   \eqn{\beta} terms.  A proportional random effects only association
#'   structure links the sub-models, with \eqn{\alpha^{(2)}} representing the
#'   association between the longitudinal and survival outcomes attributable to
#'   the deviation of the individual in question from the population mean
#'   longitudinal trajectory, and \eqn{\alpha^{(3)}} representing the
#'   association between the longitudinal and survival outcomes attributable to
#'   the deviation
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   This model accounts for between study heterogeneity using study level
#'   random effects.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#'   The code used to fit this one stage model was:
#'
#'   \code{ onestagefit2<-jointmeta1(data = jointdat, long.formula = Y ~ 1 +
#'   time + treat, long.rand.ind = c('int', 'time'), long.rand.stud = c('study',
#'   'treat'), sharingstrct = 'randprop', surv.formula = Surv(survtime, cens) ~
#'   treat, study.name = 'study', strat = F) }
#'
#'   And the code used to bootstrap the model was:
#'
#'   \code{ onestagefit2SE<-jointmetaSE(fitted = onestagefit2, n.boot = 200) }
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}
#'
"onestage2"

#' One stage jointmeta1 fit and bootstrapped standard errors
#'
#' A list of length two containing a one stage jointmeta1 fit and corresponding
#' bootstrapped standard errors.
#'
#' @format A list of 2 objects: \describe{ \item{\code{onestagefit3}}{an object
#'   of class \code{jointmeta1}} \item{\code{onestagefit3SE}}{an object of class
#'   \code{jointmeta1SE}} }
#'
#' @details These are the results of using the \code{jointmeta1} function to fit
#'   a one stage joint meta model for multi-study data, and also the bootstrap
#'   results of applying the \code{jointmetaSE} function to the resulting model
#'   fit. The data used is the \code{simdat} data available in the
#'   \code{joineRmeta} package.  This data has three levels, namely the
#'   longitudinal measurements at level 1, nested within individuals (level 2)
#'   who are themselves nested within studies (level 3).
#'
#'   The format of this model is as follows.  The structure of the longitudinal
#'   sub-model is:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   \beta_{13}study + \beta_{14}treat*study + b^{(2)}_{0ki} + b^{(2)}_{1ki}time
#'   + \epsilon_{kij}}
#'
#'   \eqn{Y_{kij}} represents the continuous longitudinal outcome for the
#'   \eqn{i}th individual in the \eqn{k}th study at the \eqn{j}th time point,
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects
#'   coefficients by \eqn{b} and the measurement error by \eqn{\epsilon}.  For
#'   the  random effects the superscript of 2 indicates that these are
#'   individual level, or level 2 random effects.  This means they take can take
#'   a unique value for each individual in the dataset. The longitudinal time
#'   variable is represented by \eqn{time}, and the treatment assignment
#'   variable (a binary factor) is represented by \eqn{treat}.  The study
#'   membership variable, a factor variable, is represented by \eqn{study}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0k}(t)exp(\beta_{21}treat +
#'   \alpha^{(2)}(b^{(2)}_{0ki} + b^{(2)}_{1ki}time)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0k}(t)}
#'   represents the unspecified baseline hazard which is stratified by study.
#'   The fixed effect coefficients are represented by \eqn{\beta} terms.  A
#'   proportional random effects only association structure links the
#'   sub-models, with \eqn{\alpha^{(2)}} representing the association between
#'   the longitudinal and survival outcomes attributable to the deviation of the
#'   individual in question from the population mean longitudinal trajectory.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   This model accounts for between study heterogeneity by including study
#'   membership and interactions between the study membership and the treatment
#'   assignment in the fixed effects of the longitudinal sub-model, and by
#'   stratifying the baseline hazard by study in the survival sub-model.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#'   The code used to fit this one stage model was:
#'
#'   \code{ onestagefit3<-jointmeta1(data = jointdat, long.formula = Y ~ 1 +
#'   time + treat*study, long.rand.ind = c('int', 'time'), sharingstrct =
#'   'randprop', surv.formula = Surv(survtime, cens) ~ treat, study.name =
#'   'study', strat = T) }
#'
#'   And the code used to bootstrap the model was:
#'
#'   \code{ onestagefit3SE<-jointmetaSE(fitted = onestagefit3, n.boot = 200,
#'   overalleffects = list(long = list(c('treat1', 'treat1:study2'), c('treat1',
#'   'treat1:study3'), c('treat1', 'treat1:study4'), c('treat1',
#'   'treat1:study5')))) }
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}
#'
"onestage3"

#' One stage jointmeta1 fit and bootstrapped standard errors
#'
#' A list of length two containing a one stage jointmeta1 fit and corresponding
#' bootstrapped standard errors.
#'
#' @format A list of 2 objects: \describe{ \item{\code{onestagefit4}}{an object
#'   of class \code{jointmeta1}} \item{\code{onestagefit4SE}}{an object of class
#'   \code{jointmeta1SE}} }
#'
#' @details These are the results of using the \code{jointmeta1} function to fit
#'   a one stage joint meta model for multi-study data, and also the bootstrap
#'   results of applying the \code{jointmetaSE} function to the resulting model
#'   fit. The data used is the \code{simdat} data available in the
#'   \code{joineRmeta} package.  This data has three levels, namely the
#'   longitudinal measurements at level 1, nested within individuals (level 2)
#'   who are themselves nested within studies (level 3).
#'
#'   The format of this model is as follows.  The structure of the longitudinal
#'   sub-model is:
#'
#'   \deqn{Y_{kij} = \beta_{10} + \beta_{11}time + \beta_{12}treat +
#'   \beta_{13}study + b^{(2)}_{0ki} + b^{(2)}_{1ki}time + b^{(3)}_{1k}treat +
#'   \epsilon_{kij}}
#'
#'   Where \eqn{Y_{kij}} represents the continuous longitudinal outcome for the
#'   \eqn{i}th individual in the \eqn{k}th study at the \eqn{j}th time point,
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects
#'   coefficients by \eqn{b} and the measurement error by \eqn{\epsilon}.  For
#'   the  random effects the superscript of 2 indicates individual level, or
#'   level 2 random effects.  This means they take can take a unique value for
#'   each individual in the dataset. A superscript of 3 indicates study level
#'   random effects, or level 3 random effects.  This means that they can take a
#'   unique value for each study in the dataset. The longitudinal time variable
#'   is represented by \eqn{time}, and the treatment assignment variable (a
#'   binary factor) is represented by \eqn{treat}.  The study membership, a
#'   factor variable, is represented by \eqn{study}.
#'
#'   The survival sub-model had format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0k}(t)exp(\beta_{21}treat +
#'   \alpha^{(2)}(b^{(2)}_{0ki} + b^{(2)}_{1ki}time) + \alpha^{(3)}(b^{(2)}_{0k}
#'   + b^{(3)}_{1k}treat)) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0k}(t)}
#'   represents the unspecified baseline hazard which is stratified by study.
#'   The fixed effect coefficients are represented by \eqn{\beta} terms. A
#'   proportional random effects only association structure links the
#'   sub-models, with \eqn{\alpha^{(2)}} representing the association between
#'   the longitudinal and survival outcomes attributable to the deviation of the
#'   individual in question from the population mean longitudinal trajectory,
#'   and \eqn{\alpha^{(3)}} representing the association between the
#'   longitudinal and survival outcomes attributable to the deviation.
#'
#'   We differentiate between the fixed effect coefficients in the longitudinal
#'   and the survival sub-models by varying the first number present in the
#'   subscript of the fixed effect, which takes a 1 for coefficients from the
#'   longitudinal sub-model and a 2 for coefficients from the survival
#'   sub-model.
#'
#'   This model accounts for between study heterogeneity using study level
#'   random effects.
#'
#'   These fits have been provided in this package for use with the package
#'   vignette, see the vignette for more information.
#'
#'   The code used to fit this one stage model was:
#'
#'   \code{ onestagefit4<-jointmeta1(data = jointdat, long.formula = Y ~ 1 +
#'   time + treat + study, long.rand.ind = c('int', 'time'), long.rand.stud =
#'   c('treat'), sharingstrct = 'randprop', surv.formula = Surv(survtime, cens)
#'   ~ treat, study.name = 'study', strat = T) }
#'
#'   And the code used to bootstrap the model was:
#'
#'   \code{ onestagefit4SE<-jointmetaSE(fitted = onestagefit4, n.boot = 200) }
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}
#'
"onestage4"
