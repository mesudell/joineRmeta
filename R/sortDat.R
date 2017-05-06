#' Sort data by order of increasing survival time
#'
#' This function is used internally in \code{jointmeta1}.  It takes longitudinal
#' and survival datasets, and orders them by increasing survival time.
#'
#' @param longdat longitudinal dataset with factor variables and interaction
#'   terms expanded out
#' @param survdat survival dataset with factor variables and interaction terms
#'   expanded out
#' @param longdat2 longitudinal dataset with original variables
#' @param survdat2 survival dataset with original variables
#'
#' @return a list containing four ordered versions of the inputted datasets. In
#'   the output \code{long.s} is the the ordered version of \code{longdat},
#'   \code{surv.s} is the ordered version of \code{survdat}, \code{long.s2} is
#'   the ordered version of longdat2, and \code{surv.s2} is the ordered version
#'   of survdat2.
#'
#' @keywords internal
sortDat <- function(longdat, survdat, longdat2, survdat2) {
  nameslong <- names(longdat)
  namessurv <- names(survdat)
  nameslong2 <- names(longdat2)
  namessurv2 <- names(survdat2)
  longid <- longdat[, 1]
  nn <- diff(match(unique(longid), longid))
  nn[length(nn) + 1] <- length(longid) - sum(nn)
  svec <- rep(survdat[, 2], nn)
  orderedsvec <- order(svec)
  sort.long <- longdat[orderedsvec, ]
  sort.long2 <- longdat2[orderedsvec, ]
  os <- order(survdat[, 2])
  sort.surv <- survdat[os, ]
  sort.surv2 <- survdat2[os, ]
  sort.long <- data.frame(sort.long)
  names(sort.long) <- nameslong
  sort.surv <- data.frame(sort.surv)
  names(sort.surv) <- namessurv
  sort.long2 <- data.frame(sort.long2)
  names(sort.long2) <- nameslong2
  sort.surv2 <- data.frame(sort.surv2)
  names(sort.surv2) <- namessurv2
  output <- list(long.s = sort.long, surv.s = sort.surv, long.s2 = sort.long2,
                 surv.s2 = sort.surv2)
}
