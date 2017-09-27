#' Extract confidence intervals
#'
#' \code{confint} returns the bootstrapped confidence intervals from a
#' \code{jointmeta1SE} object, which holds the results of bootstrapping the fit
#' from a \code{jointmeta1} function fit.
#'
#' @param object A jointmeta1SE object - the result of running the bootstrapping
#'   function \code{jointmetaSE} on a jointmeta1 object
#' @param parm A vector indicating what parameters to return confidence
#'   intervals for.  This should either be a vector of character strings, where
#'   any parameters matching the supplied parameters have their estimates and
#'   confidence intervals returned, or a numeric or integer vector of values
#'   indicating what rows of the results of the bootstapping procedure to print
#'   out
#' @param level A numerical value greater than 0 and less than 1 that indicates
#'   the requireed level of the confidence interval.  The default is \code{level
#'   = 0.95} giving 95\% confidence intervals.
#' @param ... additional arguments; currently none are used.
#'
#' @return Returns the name of variables, the part of the joint model they
#'   relate to (sub-model, variance parameter...), their estimate and their 95\%
#'   confidence interval
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
#'    #note: for demonstration purposes only - max.it restricted to 3
#'    #model would need more iterations to truely converge
#'    onestagefit<-jointmeta1(data = jointdat2, long.formula = Y ~ 1 + time +
#'                            + treat + study, long.rand.ind = c('int', 'time'),
#'                            long.rand.stud = c('treat'),
#'                            sharingstrct = 'randprop',
#'                            surv.formula = Surv(survtime, cens) ~ treat,
#'                            study.name = 'study', strat = TRUE, max.it=3)
#'
#'     \dontrun{
#'         #calculate the SE
#'         onestagefitSE <- jointmetaSE(fitted = onestagefit, n.boot = 200)
#'
#'         #extract confidence intervals
#'         confint(onestagefitSE)
#'     }
#'
confint.jointmeta1SE <- function(object, parm = NULL, level = 0.95, ...) {
  if (class(object) != "jointmeta1SE") {
    stop("fittedSE should be of class jointmeta1SE")
  }
  if (class(level) != "numeric") {
    stop("Please specify level as a numeric between 0 and 1")
  }
  if (level < 0 || level > 1) {
    stop("Please specify level as a numeric between 0 and 1")
  }
  out <- object$results
  if (level != 0.95) {
    boots <- object$bootstraps
    ci1 <- c()
    ci2 <- c()
    alpha2 <- ((1 - level)/2) * 100
    for (i in 1:nrow(out)) {
      if (nrow(boots) < 100) {
        ci1[i] <- 0
        ci2[i] <- 0
      } else {
        ci1[i] <- sort(as.numeric(boots[, i]))[alpha2/100 * nrow(boots)]
        ci2[i] <- sort(as.numeric(boots[, i]))[(100 - alpha2)/100 *
                                                 nrow(boots)]
      }
    }
    names(out)[(ncol(out) - 1):ncol(out)] <- paste((level * 100), c("%Lower",
                                                                    "%Upper"), sep = "")
    out[, 5] <- round(ci1, 4)
    out[, 6] <- round(ci2, 4)
  }
  if (is.null(parm)) {
    out<-out[, c(1:3, 5, 6)]
  } else {
    component <- as.character(out[, 1])
    temp <- component[1]
    for (i in 2:length(component)) {
      if (component[i] == "") {
        component[i] <- temp
      } else {
        temp <- component[i]
      }
    }
    out[, 1] <- component
    if (class(parm) == "character") {
      selection <- which(as.character(object$results[, 2]) %in% parm)
    } else if (class(parm) %in% c("numeric", "integer")) {
      selection <- parm
    } else {
      stop("Please supply parm as vector of character strings of parameter names,
           or numbers of their row location in supplied object")
    }
    out <- out[selection, c(1:3, 5, 6)]
    out[which(!((1:nrow(out)) %in% match(unique(out[, 1]), out[, 1]))),
        1] <- ""
    }
  return(out)
}
