#' Print function for \code{jointmeta1} objects
#'
#' A function to print a \code{\link{jointmeta1.object}}.
#'
#' @param x a jointmeta1.object, the result of fitting a jointmeta1 model
#' @param ... additional arguments; currently none are used.
#'
#' @return An object inheriting from class \code{print.jointmeta1} with all
#'   components included in \code{x} (see \code{\link{jointmeta1.object}}).
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
#'     #print the fitted multi-study joint model
#'     print(onestagefit)
#'
print.jointmeta1 <- function(x, ...) {
  if (!inherits(x,"jointmeta1")) {
    stop("Variable x should be of class jointmeta1")
  }
  cat("\nCall:\n", paste(deparse(x$Call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Random effects joint meta model\n")
  cat(" Data:", deparse(x$Call$data), "\n")
  cat("\n")
  cat("Longitudinal sub-model fixed effects:", deparse(x$Call$long.formula))
  longfixed <- x$coefficients$fixed$longitudinal
  colnames(longfixed) <- ""
  print(longfixed)
  cat("\n")
  cat("Time-to-event sub-model fixed effects:", deparse(x$Call$surv.formula))
  cat("\n")
  if (as.logical(as.character(x$Call$strat)) == TRUE) {
    cat("Strat: TRUE\n")
  } else {
    cat("Strat: FALSE\n")
  }
  survfixed <- x$coefficients$fixed$survival
  if (is.null(survfixed)) {
    cat("\n", "No survival baseline covariates specified", "\n")
  } else {
    print(survfixed)
    cat("\n")
  }
  cat("\n")
  cat("Latent association:")
  lat <- data.frame(x$coefficients$latent)
  names(lat) <- ""
  print(lat)
  cat("\n")
  cat("Variance components:\n")
  v1 <- 1
  D.diag <- diag(x$rand_cov$D)^v1
  names(D.diag) <- NULL
  sige <- as.numeric((x$sigma.e)^v1)
  names(sige) <- ""
  if (is.null(x$rand_cov$A) == FALSE) {
    A.diag <- diag(x$rand_cov$A)^v1
    names(A.diag) <- NULL
    type <- c("Individual level", rep("", length(D.diag) - 1), "Study level",
              rep("", length(A.diag) - 1), "Residual")
    paranames <- c(rownames(x$rand_cov$D), rownames(x$rand_cov$A),
                   "")
    values <- round(c(D.diag, A.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type, Name = paranames, Value = values))
    print(varoutput)
  } else {
    type <- c("Individual level", rep("", length(D.diag) - 1), "Residual")
    paranames <- c(rownames(x$rand_cov$D), "")
    values <- round(c(D.diag, sige), 7)
    varoutput <- data.frame(cbind(Type = type, Name = paranames, Value = values))
    print(varoutput)
  }
  cat("\n")
  cat("Number of studies:", x$numstudies)
  cat("\n\n")
  cat("Number of individuals per study:\n")
  print(x$n.bystudy)
  cat("\n")
  cat("Number of longitudinal observations:\n")
  print(x$nobs)
  class(x) <- c("summary.jointmeta1", class(x))
  invisible(x)
}
