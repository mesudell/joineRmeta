#' Extract formulae from joint model fit
#'
#' Extract the formula of various parts of the joint model fit
#'
#' @param fitted A \code{jointmeta1} object, the result of applying the
#'   \code{\link{jointmeta1}} function to joint longitudinal and survival data.
#' @param type A character string indicating what part of the joint model the
#'   formula should be returned for.  Specifying \code{'Longitudinal'} will
#'   result in the function returning the formula for the fixed effect portion
#'   of the longitudinal sub-model, \code{'Survival'} will result in the formula
#'   for the fixed effect portion of the survival sub-model being returned.  To
#'   return the formula for the individual level random effect, specify
#'   \code{type = 'Rand_ind'}, or to return the formula for the study level
#'   random effects (if included in the joint model), specify \code{type =
#'   'Rand_stud'}.
#'
#' @return This function returns a formula for the specified portion of the
#'   joint model fitted in the supplied \code{jointmeta1} object.
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1.object}}
#'
#' @export
#'
#' @examples
#'     #change data to jointdata format
#'     jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                           survival = simdat$survival, id = 'id',
#'                           longoutcome = 'Y', timevarying = c('time','ltime'),
#'                           survtime = 'survtime', cens = 'cens',
#'                           time = 'time')
#'
#'     #ensure variables are correctly formatted
#'     jointdat$baseline$study <- as.factor(jointdat$baseline$study)
#'     jointdat$baseline$treat <- as.factor(jointdat$baseline$treat)
#'
#'     #fit multi-study joint model
#'     onestagefit<-jointmeta1(data = jointdat, long.formula = Y ~ 1 + time +
#'                             treat + study, long.rand.ind = c('int', 'time'),
#'                             long.rand.stud = c('treat'),
#'                             sharingstrct = 'randprop',
#'                             surv.formula = Surv(survtime, cens) ~ treat,
#'                             study.name = 'study', strat = T)
#'
#'     #return the formula for the longitudinal fixed effects
#'     formula(onestagefit, type = 'Longitudinal')
#'
#'     #return the formula for the time-to-event fixed effects
#'     formula(onestagefit, type = 'Survival')
#'
#'     #return the formula for the individual level random effects
#'     formula(onestagefit, type = 'Rand_ind')
#'
#'     #return the formula for the study level random effects
#'     formula(onestagefit, type = 'Rand_study')
#'
#'     #fit multi-study joint model
#'     onestagefit2<-jointmeta1(data = jointdat, long.formula = Y ~ 1 + time +
#'                              treat + study, long.rand.ind = c('int', 'time'),
#'                              sharingstrct = 'randprop',
#'                              surv.formula = Surv(survtime, cens) ~ treat,
#'                              study.name = 'study', strat = T)
#'
#'     #return the formula for the longitudinal fixed effects
#'     formula(onestagefit2, type = 'Longitudinal')
#'
#'     #return the formula for the time-to-event fixed effects
#'     formula(onestagefit2, type = 'Survival')
#'
#'     #return the formula for the individual level random effects
#'     formula(onestagefit2, type = 'Rand_ind')
#'
#'     #return the formula for the study level random effects
#'     #note this should produce an error as no study level random effects
#'     #were specified for fit onestagefit2
#'     formula(onestagefit2, type = 'Rand_study')
#'
formula.jointmeta1 <- function(fitted, type = c("Longitudinal", "Survival",
                                                "Rand_ind", "Rand_stud")) {
  if (class(fitted) != "jointmeta1") {
    stop("Variable fitted should be of class jointmeta1")
  }
  if (missing(type) || !(type %in% c("Longitudinal", "Survival", "Rand_ind",
                                     "Rand_stud"))) {
    stop("type should be one of \"Longitudinal\", \"Survival\",
         \"Rand_ind\", \"Rand_stud\"")
  }
  if (type == "Longitudinal") {
    out <- as.formula(fitted$Call$long.formula, env = globalenv())
  } else if (type == "Survival") {
    out <- as.formula(fitted$Call$surv.formula, env = globalenv())
  } else if (type == "Rand_ind") {
    randnames <- rownames(fitted$rand_cov$D)
    randnames[grepl("(Intercept)", randnames)] <- "1"
    out <- as.formula(paste("~", paste(randnames, collapse = " +"),
                            sep = ""), env = globalenv())
  } else if (type == "Rand_stud") {
    if (is.null(fitted$rand_cov$A)) {
      stop("No study level random effects in supplied model fit")
    }
    randnames <- rownames(fitted$rand_cov$A)
    randnames[grepl("(Intercept)", randnames)] <- "1"
    out <- as.formula(paste("~", paste(randnames, collapse = " +"),
                            sep = ""), env = globalenv())
  }
  return(out)
}
