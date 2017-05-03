#' Print function for \code{jointmeta1} objects
#'
#' A function to print a
#'     \code{\link{jointmeta1.object}}.
#'
#' @param fitted a jointmeta1.object, the result of fitting a jointmeta1 model
#'
#' @return An object inheriting from class \code{print.jointmeta1} with all
#'     components included in \code{fitted} (see
#'     \code{\link{jointmeta1.object}}).
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
#'     #print the fitted multi-study joint model
#'     print(onestagefit)
#'
print.jointmeta1 <- function(fitted) {
  if(class(fitted) != "jointmeta1") {
    stop("Variable fitted should be of class jointmeta1")
  }
  cat("\nCall:\n", paste(deparse(fitted$Call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Random effects joint meta model\n")
  cat(" Data:", deparse(fitted$Call$data), "\n")
  cat("\n")
  cat("Longitudinal sub-model fixed effects:",
      deparse(fitted$Call$long.formula))
  longfixed <- fitted$coefficients$fixed$longitudinal
  colnames(longfixed) <- ""
  print(longfixed)
  cat("\n")
  cat("Time-to-event sub-model fixed effects:",
      deparse(fitted$Call$surv.formula))
  cat("\n")
  if(as.logical(as.character(fitted$Call$strat)) == TRUE){
    cat("Strat: TRUE\n")
  }else{
    cat("Strat: FALSE\n")
  }
  survfixed <- fitted$coefficients$fixed$survival
  if(is.null(survfixed)) {
    cat("\n", "No survival baseline covariates specified",
        "\n")
  }else {
    print(survfixed)
    cat("\n")
  }
  cat("\n")
  cat("Latent association:")
  lat <- data.frame(fitted$coefficients$latent)
  names(lat) <- ""
  print(lat)
  cat("\n")
  cat("Variance components:\n")
  v1 <- 1
  D.diag <- diag(fitted$rand_cov$D)^v1
  names(D.diag) <- NULL
  sige <- as.numeric((fitted$sigma.e)^v1)
  names(sige) <- ""
  if(is.null(fitted$rand_cov$A) == FALSE) {
    A.diag <- diag(fitted$rand_cov$A)^v1
    names(A.diag) <- NULL
    type <- c("Individual level", rep("", length(D.diag)-1),
              "Study level", rep("", length(A.diag)-1),
              "Residual")
    paranames <- c(rownames(fitted$rand_cov$D),
                   rownames(fitted$rand_cov$A), "")
    values <- round(c(D.diag, A.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type,
                                  Name = paranames, Value = values))
    print(varoutput)
  }else {
    type <- c("Individual level", rep("", length(D.diag)-1),
              "Residual")
    paranames <- c(rownames(fitted$rand_cov$D), "")
    values <- round(c(D.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type,
                                  Name = paranames, Value = values))
    print(varoutput)
  }
  cat("\n")
  cat("Number of studies:", fitted$numstudies)
  cat("\n\n")
  cat("Number of individuals per study:\n")
  print(fitted$n.bystudy)
  cat("\n")
  cat("Number of longitudinal observations:\n")
  print(fitted$nobs)
  class(fitted) <- c("summary.jointmeta1", class(fitted))
  invisible(fitted)
}
