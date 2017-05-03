#' Function to extract the estimated covariance matrices for the random effects
#'     specified in the model
#'
#' A function to allow the random effects covariance matrix for a particular
#'     level of random effects specified in the sub-model to be extracted from
#'     the \code{jointmeta1} model fit.
#'
#' @param fitted a \code{jointmeta1.object}
#' @param type a character string indicating what level the random effects
#'     covariance matrix should be returned for.  If the individual level random
#'     effects covariance matrix is required then \code{type = "individual"}.
#'     If the study level random effects covariance matrix is required then
#'     \code{type = "study"}.  Note that if study level random effects are not
#'     included in the model, then attempting to extract them will result in an
#'     error message.
#'
#' @return a matrix of dimensions equal to the number of random effects at the
#'     level specified by the \code{type} parameter.
#'
#' @export
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1.object}}
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
#'     #extract the individual level random effects covariance matrix
#'     rancov(onestagefit, type = "individual")
#'
#'     #extract the study level random effects covariance matrix
#'     rancov(onestagefit, type = "study")
#'
#'
rancov<-function(fitted, type=c("individual", "study")) {
  if(class(fitted) != "jointmeta1") {
    stop("Variable fitted should be of class jointmeta1")
  }
  if(missing(type) || !(type %in% c("individual", "study"))) {
    stop("type should be one of \"individual\", \"study\"")
  }
  if(type == "individual") {
    ranef.ind <- fitted$rand_cov$D
    return(ranef.ind)
  }else if(type == "study"){
    if(is.null(fitted$rand_cov$A)) {
      stop("No study level random effects specified in this model")
    }else {
      ranef.stud <- fitted$rand_cov$A
    }
    return(ranef.stud)
  }
}
