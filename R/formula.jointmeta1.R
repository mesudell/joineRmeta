#' Extract formulae from joint model fit
#'
#' Extract the formula of various parts of the joint model fit
#'
#' @param x A \code{jointmeta1} object, the result of applying the
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
#' @param ... additional arguments; currently none are used.
#'
#' @return This function returns a formula for the specified portion of the
#'   joint model fitted in the supplied \code{jointmeta1} object.
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1.object}}
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
#'    #note: for demonstration purposes only - max.it restricted to 5
#'    #model would need more iterations to truely converge
#'    onestagefit<-jointmeta1(data = jointdat2, long.formula = Y ~ 1 + time +
#'                            + treat + study, long.rand.ind = c('int', 'time'),
#'                            long.rand.stud = c('treat'),
#'                            sharingstrct = 'randprop',
#'                            surv.formula = Surv(survtime, cens) ~ treat,
#'                            study.name = 'study', strat = TRUE, max.it=5)
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
#'     formula(onestagefit, type = 'Rand_stud')
#'
#'
formula.jointmeta1 <- function(x, type = c("Longitudinal", "Survival",
                                           "Rand_ind", "Rand_stud"), ...) {
  if (!inherits(x, "jointmeta1")) {
    stop("Variable x should be of class jointmeta1")
  }
  if (missing(type) || !(type %in% c("Longitudinal", "Survival", "Rand_ind",
                                     "Rand_stud"))) {
    stop("type should be one of \"Longitudinal\", \"Survival\",
         \"Rand_ind\", \"Rand_stud\"")
  }
  if (type == "Longitudinal") {
    out <- as.formula(x$Call$long.formula, env = globalenv())
  } else if (type == "Survival") {
    out <- as.formula(x$Call$surv.formula, env = globalenv())
  } else if (type == "Rand_ind") {
    randnames <- rownames(x$rand_cov$D)
    randnames[grepl("(Intercept)", randnames)] <- "1"
    out <- as.formula(paste("~", paste(randnames, collapse = " +"),
                            sep = ""), env = globalenv())
  } else if (type == "Rand_stud") {
    if (is.null(x$rand_cov$A)) {
      stop("No study level random effects in supplied model fit")
    }
    randnames <- rownames(x$rand_cov$A)
    randnames[grepl("(Intercept)", randnames)] <- "1"
    out <- as.formula(paste("~", paste(randnames, collapse = " +"),
                            sep = ""), env = globalenv())
  }
  return(out)
}
