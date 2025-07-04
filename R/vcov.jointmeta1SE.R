#' Extract the variance covariance matrix from the bootstrapped results
#'
#' Function applied to a \code{jointmeta1SE} object, the result of the
#' \code{jointmetaSE} function to extract the variance covariance matrix for the
#' estimated model parameters
#'
#' @param object an object of class \code{jointmeta1SE}
#' @param ... additional arguments; currently none are used.
#'
#' @return a variance covariance matrix for the fixed effects from the
#'   longitudinal sub-model, the time-to-event sub-model, the association
#'   parameters, the random effects and the error term.
#'
#' @export
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmetaSE}},
#'   \code{\link{jointmeta1SE.object}}
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
#'         #extract the variance covariance matrix
#'         vcov(onestagefitSE)
#'     }
#'
vcov.jointmeta1SE <- function(object, ...) {
  if (!inherits(object,"jointmeta1SE")) {
    stop("object should be of class jointmeta1SE")
  }
  object$covmat
}

