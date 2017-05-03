#' Function to extract estimated random effects
#'
#' This function extracts the estimated values of the random effects from a
#' supplied \code{jointmeta1} fit.
#'
#' @param fitted a \code{jointmeta1} object (the result of fitting a model
#'     using \code{\link{jointmeta1}}, see \code{\link{jointmeta1.object}})
#' @param type the type of random effects to return.  Set
#'     \code{type = "individual"} to return the estimates of the individual
#'     level random effects.  Set \code{type = "study"} to return the estimates
#'     of the level random effects if included in the model.
#'
#' @return If \code{type = "individual"} then a list of matrices containing the
#'     individual level random effects is returned.  This list is of length
#'     equal to the number of studies in the dataset.  Each matrix has number of
#'     rows equal to the number of individuals in the corresponding study, and
#'     number of columns equal to the number of individual level random effects.
#'
#'     If \code{type = "study"} then if study level random effects are present
#'     in the supplied model fit, a matrix of the estimated study level random
#'     effects is returned, with number of rows equal to the number of
#'     studies in the dataset, and number of columns equal to the number of
#'     study level random effects.  If study level random effects are requested
#'     but are not present in the supplied model fit, an error message is
#'     returned.
#'
#' @export
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1.object}},
#'     \code{\link{fixef}}
#'
#' @examples
#'     #change data to jointdata format
#'     jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                           survival = simdat$survival, id = "id",
#'                           longoutcome = "Y",
#'                           timevarying = c("time","ltime"),
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
#'     #extract the individual level random effects covariance matrix
#'     ranef(onestagefit, type = "individual")
#'
#'     #extract the study level random effects covariance matrix
#'     ranef(onestagefit, type = "study")
ranef.jointmeta1 <- function(fitted, type=c("individual", "study")) {
  if(class(fitted) != "jointmeta1") {
    stop("Variable fitted should be of class jointmeta1")
  }
  if(missing(type) || !(type %in% c("individual", "study"))) {
    stop("type should be one of \"individual\", \"study\"")
  }
  if(type == "individual") {
    ranef.ind <- fitted$random_cond$random_ind
    return(ranef.ind)
  }else if(type == "study"){
    if(is.null(fitted$random_cond$random_stud)) {
      stop("No study level random effects in supplied model fit")
    }
    ranef.stud <- fitted$random_cond$random_stud
    return(ranef.stud)
  }
}
