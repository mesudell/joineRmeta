#' Extract fixed effects
#'
#' Function to extract the estimated fixed effects from a jointmeta1 model fit
#'
#' @param fitted A joint model fit from the function jointmeta1
#' @param type Type of fixed effects to extract.  To extract fixed effects from
#'     longitudinal sub-model set \code{type = "Longitudinal"}, from the
#'     time-to-event sub-model set \code{type = "Survival"}, or extract
#'     the latent association parameters set to \code{type = "Latent"}.
#'
#' @return The function returns a vector of the fixed effects from the
#'     specified part of the supplied jointmeta1 model fit.
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
#'     #fit multi-study joint model
#'     onestagefit<-jointmeta1(data = jointdat, long.formula = Y ~ 1 + time +
#'                             treat + study, long.rand.ind = c("int", "time"),
#'                             long.rand.stud = c("treat"),
#'                             sharingstrct = "randprop",
#'                             surv.formula = Surv(survtime, cens) ~ treat,
#'                             study.name = "study", strat = T)
#'
#'     #extract longitudinal fixed effects
#'     fixef(onestagefit, type = "Longitudinal")
#'
#'     #extract time-to-event fixed effects
#'     fixef(onestagefit, type = "Survival")
#'
#'     #extract latent (association parameter) coefficients
#'     fixef(onestagefit, type = "Latent")


fixef.jointmeta1 <- function(fitted,
                             type=c("Longitudinal", "Survival", "Latent")) {
  if(class(fitted) != "jointmeta1") {
    stop("Variable fitted should be of class jointmeta1")
  }
  if(missing(type) || !(type %in% c("Longitudinal", "Survival", "Latent"))) {
    stop("type should be one of \"Longitudinal\", \"Survival\", \"Latent\"")
  }
  if(type == "Longitudinal") {
    fixef.long <- fitted$coefficients$fixed$longitudinal
    fixef.long <- as.vector(fixef.long)[, 1]
    names(fixef.long) <- rownames(fitted$coefficients$fixed$longitudinal)
    return(fixef.long)
  }else if(type == "Survival") {
    if(is.null(fitted$coefficients$fixed$survival)) {
      stop("No fixed effects for survival sub-model")
    }else {
      fixef.surv <- fitted$coefficients$fixed$survival
      return(fixef.surv)
    }
  }else if(type == "Latent") {
    fixef.latent <- fitted$coefficients$latent
    return(fixef.latent)
  }
}
