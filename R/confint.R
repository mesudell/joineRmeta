#' Extract confidence intervals
#'
#' \code{confint} returns the bootstrapped confidence intervals from a
#'     \code{jointmeta1SE} object, which holds the results of bootstrapping the
#'     fit from a \code{jointmeta1} function fit.
#'
#' @param fittedSE A jointmeta1SE object - the result of running the
#'     bootstrapping function \code{jointmetaSE} on a jointmeta1 object
#'
#' @return Returns the name of variables, the part of the joint model they
#'     relate to (sub-model, variance parameter...), their estimate and
#'     their 95\% confidence interval
#'
#' @seealso \code{\link{jointmeta1}}
#'
#' @export
#'
#' @examples
#'     #change data to jointdata format
#'     jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                           survival = simdat$survival, id = "id",
#'                           longoutcome = "Y", timevarying = c("time","ltime"),
#'                           survtime = "survtime", cens = "cens",
#'                           time = "time")
#'
#'     #ensure variables are correctly formatted
#'     jointdat$baseline$study <- as.factor(jointdat$baseline$study)
#'     jointdat$baseline$treat <- as.factor(jointdat$baseline$treat)
#'
#'     #fit multi-study joint model and calculate the standard errors
#'     onestagefit<-jointmeta1(data = jointdat, long.formula = Y ~ 1 + time +
#'                             treat + study, long.rand.ind = c("int", "time"),
#'                             long.rand.stud = c("treat"),
#'                             sharingstrct = "randprop",
#'                             surv.formula = Surv(survtime, cens) ~ treat,
#'                             study.name = "study", strat = T)
#'
#'     onestagefitSE<-jointmetaSE(fitted = onestagefit, n.boot = 200)
#'
#'     #extract confidence intervals
#'     confint(onestagefitSE)
#'

confint.jointmeta1SE <- function(fittedSE) {
  if(class(fittedSE) != "jointmeta1SE") {
    stop("fittedSE should be of class jointmeta1SE")
  }
  fittedSE$results[, c(1:3, 5, 6)]
}
