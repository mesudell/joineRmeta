#' Function to extract estimated random effects
#'
#' This function extracts the estimated values of the random effects from a
#' supplied \code{jointmeta1} fit.
#'
#' @param object a \code{jointmeta1} object (the result of fitting a model using
#'   \code{\link{jointmeta1}}, see \code{\link{jointmeta1.object}})
#' @param type the type of random effects to return.  Set \code{type =
#'   'individual'} to return the estimates of the individual level random
#'   effects.  Set \code{type = 'study'} to return the estimates of the level
#'   random effects if included in the model.
#' @param ... additional arguments; currently none are used.
#'
#' @return If \code{type = 'individual'} then a list of matrices containing the
#'   individual level random effects is returned.  This list is of length equal
#'   to the number of studies in the dataset.  Each matrix has number of rows
#'   equal to the number of individuals in the corresponding study, and number
#'   of columns equal to the number of individual level random effects.
#'
#'   If \code{type = 'study'} then if study level random effects are present in
#'   the supplied model fit, a matrix of the estimated study level random
#'   effects is returned, with number of rows equal to the number of studies in
#'   the dataset, and number of columns equal to the number of study level
#'   random effects.  If study level random effects are requested but are not
#'   present in the supplied model fit, an error message is returned.
#'
#' @export
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1.object}},
#'   \code{\link{fixef.jointmeta1}}
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
#'     #extract the individual level random effects covariance matrix
#'     ranef(onestagefit, type = 'individual')
#'
#'     #extract the study level random effects covariance matrix
#'     ranef(onestagefit, type = 'study')
ranef.jointmeta1 <- function(object, type = c("individual", "study"), ...) {
  if (!inherits(object,"jointmeta1")) {
    stop("Variable object should be of class jointmeta1")
  }
  if (missing(type) || !(type %in% c("individual", "study"))) {
    stop("type should be one of \"individual\", \"study\"")
  }
  if (type == "individual") {
    ranef.ind <- object$coefficients$random$random_ind
    return(ranef.ind)
  } else if (type == "study") {
    if (is.null(object$coefficients$random$random_stud)) {
      stop("No study level random effects in supplied model fit")
    }
    ranef.stud <- object$coefficients$random$random_stud
    return(ranef.stud)
  }
}

