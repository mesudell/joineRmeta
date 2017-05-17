#' Extract fixed effects
#'
#' Function to extract the estimated fixed effects from a jointmeta1 model fit
#'
#' @param object A joint model fit from the function jointmeta1
#' @param type Type of fixed effects to extract.  To extract fixed effects from
#'   longitudinal sub-model set \code{type = 'Longitudinal'}, from the
#'   time-to-event sub-model set \code{type = 'Survival'}, or extract the latent
#'   association parameters set to \code{type = 'Latent'}.
#' @param ... additional arguments; currently none are used.
#'
#' @return The function returns a vector of the fixed effects from the specified
#'   part of the supplied jointmeta1 model fit.
#'
#' @seealso \code{\link{jointmeta1}}
#'
#' @export
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
#'    #note: for demonstration purposes only - max.it restricted to 2
#'    #model would need more iterations to truely converge
#'    onestagefit<-jointmeta1(data = jointdat2, long.formula = Y ~ 1 + time +
#'                            + treat + study, long.rand.ind = c('int'),
#'                            long.rand.stud = c('treat'),
#'                            sharingstrct = 'randprop',
#'                            surv.formula = Surv(survtime, cens) ~ treat,
#'                            study.name = 'study', strat = TRUE, max.it = 2)
#'
#'     #extract longitudinal fixed effects
#'     fixef(onestagefit, type = 'Longitudinal')

fixef.jointmeta1 <- function(object, type = c("Longitudinal", "Survival",
                                              "Latent"), ...) {
  if (class(object) != "jointmeta1") {
    stop("Variable object should be of class jointmeta1")
  }
  if (missing(type) || !(type %in% c("Longitudinal", "Survival", "Latent"))) {
    stop("type should be one of \"Longitudinal\", \"Survival\", \"Latent\"")
  }
  if (type == "Longitudinal") {
    fixef.long <- object$coefficients$fixed$longitudinal
    fixef.long <- as.vector(fixef.long)[, 1]
    names(fixef.long) <- rownames(object$coefficients$fixed$longitudinal)
    return(fixef.long)
  } else if (type == "Survival") {
    if (is.null(object$coefficients$fixed$survival)) {
      stop("No fixed effects for survival sub-model")
    } else {
      fixef.surv <- object$coefficients$fixed$survival
      return(fixef.surv)
    }
  } else if (type == "Latent") {
    fixef.latent <- object$coefficients$latent
    return(fixef.latent)
  }
}
