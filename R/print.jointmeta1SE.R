#' Print function for \code{jointmeta1SE} objects
#'
#' This function extracts the results of the bootstrapping function
#' \code{jointmetaSE} applied to the results of a one stage joint meta model
#' fitted using \code{jointmeta1}, namely \code{jointmeta1SE} objects. The
#' bootstrapping function returns not only results but also the covariance
#' matrix for the estimated parameters and the results of each bootstrap.  This
#' print function allows just the results of the bootstraps to easily be
#' displayed.
#'
#' @param x a \code{jointmeta1SE} object, the result of applying
#'   \code{jointmetaSE} to the joint model fit obtained by applying
#'   \code{jointmeta1} to a multi-study joint data set.
#' @param ... additional arguments; currently none are used.
#'
#' @return a data frame containing the estimates, standard errors and 95\%
#'   confidence intervals for the parameters from the model and any overall
#'   effects requested in the \code{jointmetaSE} function call.
#'
#' @export
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmetaSE}},
#'   \code{\link{jointmeta1SE.object}}, \code{\link{jointmeta1.object}}
#'
#' @examples
#'    #change example data to jointdata object
#'    jointdat2<-tojointdata(longitudinal = simdat2$longitudinal,
#'    survival = simdat2$survival, id = 'id',longoutcome = 'Y',
#'    timevarying = c('time','ltime'),
#'    survtime = 'survtime', cens = 'cens',time = 'time')
#'
#'    #set variables to factors
#'    jointdat2$baseline$study <- as.factor(jointdat2$baseline$study)
#'    jointdat2$baseline$treat <- as.factor(jointdat2$baseline$treat)
#'
#'    #fit multi-study joint model
#'    #note: for demonstration purposes only - max.it restricted to 5
#'    #model would need more iterations to truely converge
#'    onestagefit<-jointmeta1(data = jointdat2, long.formula = Y ~ 1 + time +
#'                            + treat + study, long.rand.ind = c('int', 'time'),
#'                            long.rand.stud = c('treat'),
#'                            sharingstrct = 'randprop',
#'                            surv.formula = Surv(survtime, cens) ~ treat,
#'                            study.name = 'study', strat = TRUE, max.it=5)
#'
#'     \dontrun{
#'         #calculate the SE
#'         onestagefitSE <- jointmetaSE(fitted = onestagefit, n.boot = 200)
#'
#'         #print the results of the bootstrap function
#'         print(onestagefitSE)
#'     }
#'
print.jointmeta1SE <- function(x, ...) {
  print(x$results)
}
