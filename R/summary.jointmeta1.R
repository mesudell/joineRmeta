#' Summary function for jointmeta1
#'
#' A function to provide a summary of a \code{\link{jointmeta1.object}}.
#'
#' @param object a jointmeta1.object, the result of fitting a jointmeta1 model
#' @param variance a logical if set to \code{TRUE} the variances of the
#'   measurement errors and the random effects are returned, else if
#'   \code{FALSE} then the standard errors are returned.
#' @param ... additional arguments; currently none are used.
#'
#' @return An object inheriting from class \code{summary.jointmeta1} with all
#'   components included in \code{object} (see \code{\link{jointmeta1.object}}).
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
#'     #request a summary of the fitted model, with variances printed
#'     summary(onestagefit, variance = TRUE)
#'
#'     #request a summary of the fitted model, with standard errors printed
#'     summary(onestagefit, variance = FALSE)
summary.jointmeta1 <- function(object, variance = TRUE, ...) {
  if (class(object) != "jointmeta1") {
    stop("Variable object should be of class jointmeta1")
  }
  cat("\nCall:\n", paste(deparse(object$Call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Random effects joint meta model\n")
  cat(" Data:", deparse(object$Call$data), "\n")
  cat(" Log-likelihood:", object$loglik$jointlhood, "\n")
  cat(" AIC:", object$AIC, "\n")
  cat("\n")
  cat("Longitudinal sub-model fixed effects:", deparse(object$Call$long.formula))
  longfixed <- object$coefficients$fixed$longitudinal
  colnames(longfixed) <- ""
  print(longfixed)
  cat("\n")
  cat("Time-to-event sub-model fixed effects:", deparse(object$Call$surv.formula))
  cat("\n")
  if (as.logical(as.character(object$Call$strat)) == TRUE) {
    cat("Strat: TRUE\n")
  } else {
    cat("Strat: FALSE\n")
  }
  survfixed <- object$coefficients$fixed$survival
  if (is.null(survfixed)) {
    cat("\n", "No survival baseline covariates specified", "\n")
  } else {
    print(survfixed)
    cat("\n")
  }
  cat("\n")
  cat("Latent association:")
  lat <- data.frame(object$coefficients$latent)
  names(lat) <- ""
  print(lat)
  cat("\n")
  cat("Variance components:\n")
  v1 <- 0.5
  if (variance) {
    v1 <- 1
  }
  D.diag <- diag(object$rand_cov$D)^v1
  names(D.diag) <- NULL
  sige <- as.numeric((object$sigma.e)^v1)
  names(sige) <- ""
  if (is.null(object$rand_cov$A) == FALSE) {
    A.diag <- diag(object$rand_cov$A)^v1
    names(A.diag) <- NULL
    type <- c("Individual level", rep("", length(D.diag) - 1), "Study level",
              rep("", length(A.diag) - 1), "Residual")
    paranames <- c(rownames(object$rand_cov$D), rownames(object$rand_cov$A),
                   "")
    values <- round(c(D.diag, A.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type, Name = paranames, Value = values))
    print(varoutput)
  } else {
    type <- c("Individual level", rep("", length(D.diag) - 1), "Residual")
    paranames <- c(rownames(object$rand_cov$D), "")
    values <- round(c(D.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type, Name = paranames, Value = values))
    print(varoutput)
  }
  if (!variance) {
    cat(" Note: the above are standard deviations\n")
  }
  cat("\n")
  if (object$convergence == TRUE) {
    cat("Convergence at iteration:", object$numIter, "\n")
  } else {
    cat("Convergence not achieved\n")
  }
  cat("\n")
  cat("Number of studies:", object$numstudies)
  cat("\n\n")
  cat("Number of individuals per study:\n")
  print(object$n.bystudy)
  cat("\n")
  cat("Number of longitudinal observations:\n")
  print(object$nobs)
  class(object) <- c("summary.jointmeta1", class(object))
  invisible(object)
}

