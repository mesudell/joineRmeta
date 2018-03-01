#' Function for longitudinal starting values
#'
#' Internal function to estimate the starting values for the EM algorithm for
#' the longitudinal sub-model
#'
#' @param long.formula.orig the original longitudinal formula as supplied to the
#'   function call
#' @param longdat2 the longitudinal dataset, without factors expanded into dummy
#'   variables, ordered by increasing survival time
#' @inheritParams jointmeta1
#' @inheritParams EMalgRandprop
#'
#'
#' @return A list of results from the separate longitudinal fit is returned. The
#'   elements of this list are: \describe{
#'
#'   \item{\code{beta1}}{a data frame containing the estimates of the fixed
#'   effects from the longitudinal sub-model}
#'
#'   \item{\code{sigma.e}}{the estimate of the variance for the measurement
#'   error variance for the longitudinal model}
#'
#'   \item{\code{D}}{a data frame containing the estimate of the covariance
#'   matrix for the individual level random effects}
#'
#'   \item{\code{log.like.long}}{the log-likelihood for the separate
#'   longitudinal model}
#'
#'   \item{\code{randstart.ind}}{the conditional modes of the individual level
#'   random effects given the data and the parameter estimates from the separate
#'   longitudinal model. The data frame has number of columns equal to the
#'   number of individual level random effects, and the number of rows equal to
#'   the number of individuals in the dataset}
#'
#'   \item{\code{randstart.ind.cov}}{a list of the conditional covariance
#'   matrices for each individual for the individual level random effects. The
#'   list is of length equal to the number of individuals in the dataset. Each
#'   element of the list is a matrix of dimensions equal to the number of study
#'   level random effects in the model.}
#'
#'   \item{\code{A}}{a data frame containing the estimate of the covariance
#'   matrix for the study level random effects.  This is present in the results
#'   only if study level random effects are included in the model.}
#'
#'   \item{\code{randstart.stud}}{the conditional modes of the study level
#'   random effects given the data and the parameter estimates from the separate
#'   longitudinal model.  The data frame has number of columns equal to the
#'   number of study level random effects, and number of rows equal to the
#'   number of studies in the dataset.  This is present in the results only if
#'   study level random effects are included in the model.}
#'
#'   \item{\code{randstart.stud.cov}}{a list of the conditional covariance
#'   matrices for each study for the study level random effects. The list is of
#'   length equal to the number of studies in the dataset. Each element of the
#'   list is a matrix of dimensions equal to the number of study level random
#'   effects in the model.  This is present in the results only if study level
#'   random effects are included in the model.}
#'
#'   \item{\code{modelfit}}{the initial longitudinal model fit, fitted using the
#'   \code{\link[lme4]{lmer}} function}
#'
#'   }
#' @keywords internal
#' @import lme4
longst <- function(longdat, long.formula.orig, long.rand.ind, long.rand.stud = NULL,
                   longdat2, id.name, study.name, studies) {
  if ("(Intercept)" %in% long.rand.ind) {
    remainingrandind <- long.rand.ind[-which((long.rand.ind %in% "(Intercept)") ==
                                               TRUE)]
    if (length(remainingrandind) > 0) {
      remainingrandind <- unlist(lapply(1:length(remainingrandind),
                                        function(u) {
                                          names(longdat[5:ncol(longdat)])[match(remainingrandind[u],
                                                                                names(longdat[5:ncol(longdat)]))]
                                        }))
      rf <- paste("(", paste("1", paste(remainingrandind, collapse = " + "),
                             sep = " + "), "|", id.name, ")", collapse = "", sep = "")
    } else {
      rf <- paste("(1|", id.name, ")", collapse = "", sep = "")
    }
    num.rand.ind <- 1 + length(remainingrandind)
  }
  if ("noint" %in% long.rand.ind) {
    remainingrandind <- long.rand.ind[-which((long.rand.ind %in% "noint") ==
                                               TRUE)]
    if (length(remainingrandind) > 0) {
      remainingrandind <- names(longdat[5:ncol(longdat)])[which(names(longdat[5:ncol(longdat)]) %in%
                                                                  remainingrandind)]
      rf <- paste("(", paste("-1", paste(remainingrandind, collapse = " + "),
                             sep = " + "), "|", id.name, ")", collapse = "", sep = "")
    }
    if (length(remainingrandind) == 0) {
      stop("Please specify at least one individual
           level random effect in long.rand.ind")
    }
    num.rand.ind <- length(remainingrandind)
    }
  if (is.null(long.rand.stud) == FALSE) {
    if (study.name %in% long.rand.stud) {
      remainingrandstudy <- long.rand.stud[-which((long.rand.stud %in%
                                                     study.name) == TRUE)]
      if (length(remainingrandstudy) > 0) {
        remainingrandstudy <- names(longdat[5:ncol(longdat)])[which(names(longdat[5:ncol(longdat)]) %in%
                                                                      remainingrandstudy)]
        rf <- paste(rf, " + (", paste("1", paste(remainingrandstudy,
                                                 collapse = " + "), sep = " + "), "|", study.name, ")",
                    collapse = "", sep = "")
      }
      if (length(remainingrandstudy) == 0) {
        rf <- paste(rf, " + (1|", study.name, ")", collapse = "",
                    sep = "")
      }
      num.rand.study <- 1 + length(remainingrandstudy)
    }
    if ((study.name %in% long.rand.stud) == FALSE) {
      cols <- c()
      for (i in 1:length(long.rand.stud)) {
        cols <- c(cols, grep(long.rand.stud[i], names(longdat[5:ncol(longdat)])))
      }
      long.rand.stud2 <- names(longdat[5:ncol(longdat)])[cols]
      rf <- paste(rf, " + (", paste("-1", paste(long.rand.stud2,
                                                collapse = " + "), sep = " + "), "|", study.name, ")",
                  collapse = "", sep = "")
      num.rand.study <- length(long.rand.stud)
    }
  }
  long.formula.orig2 <- as.character(long.formula.orig)
  long.formula.orig2[3] <- gsub("\\(Intercept\\)", "1", long.formula.orig2[3])
  long.start <- lmer(as.formula(paste(long.formula.orig2[2], "~", long.formula.orig2[3],
                                      " + ", rf, collapse = "")), data = data.frame(longdat), na.action = na.omit)
  randout <- ranef(long.start, condVar = TRUE)
  if (is.null(long.rand.stud) == FALSE) {
    r <- num.rand.study
    if (r > 1) {
      A <- as.data.frame(VarCorr(long.start)[which(names(VarCorr(long.start)) %in%
                                                     study.name)])
    }
    if (r == 1) {
      A <- as.data.frame(VarCorr(long.start)[which(names(VarCorr(long.start)) %in%
                                                     study.name)])
    }
    randstart.stud <- ranef(long.start)[[2]]
    rownames(randstart.stud) <- studies[as.numeric(rownames(randstart.stud))]
    randstart.stud.cov <- lapply(1:length(studies), function(u) {
      attr(randout[[study.name]], "postVar")[1:r, 1:r, u]
    })
    names(randstart.stud.cov) <- rownames(randstart.stud)
    if ((study.name %in% long.rand.stud) == TRUE) {
      remainingrandstudy <- long.rand.stud[-which((long.rand.stud %in%
                                                     study.name) == TRUE)]
      if (length(remainingrandstudy) > 0) {
        rownames(A) <- c("(Intercept)", remainingrandstudy)
        colnames(A) <- c("(Intercept)", remainingrandstudy)
        colnames(randstart.stud) <- c("(Intercept)", remainingrandstudy)
        randstart.stud.cov <- lapply(randstart.stud.cov, function(u) {
          rownames(u) <- c("(Intercept)", remainingrandstudy)
          colnames(u) <- c("(Intercept)", remainingrandstudy)
          u
        })
      } else {
        rownames(A) <- "(Intercept)"
        colnames(A) <- "(Intercept)"
        colnames(randstart.stud) <- "(Intercept)"
        randstart.stud.cov <- lapply(randstart.stud.cov, function(u) {
          names(u) <- "(Intercept)"
          u
        })
      }
    }
    if ((study.name %in% long.rand.stud) == FALSE) {
      rownames(A) <- long.rand.stud
      colnames(A) <- long.rand.stud
      colnames(randstart.stud) <- long.rand.stud
      if (r > 1) {
        randstart.stud.cov <- lapply(randstart.stud.cov, function(u) {
          rownames(u) <- long.rand.stud
          colnames(u) <- long.rand.stud
          u
        })
      } else {
        randstart.stud.cov <- lapply(randstart.stud.cov, function(u) {
          names(u) <- long.rand.stud
          u
        })
      }
    }
  }
  q <- num.rand.ind
  randstart.ind <- randout[[1]]
  randstart.ind.cov <- lapply(1:length(unique(longdat[[id.name]])), function(u) {
    attr(randout[[id.name]], "postVar")[1:q, 1:q, u]
  })
  names(randstart.ind.cov) <- rownames(randstart.ind)
  if (q > 1) {
    D <- as.data.frame(VarCorr(long.start)[which(names(VarCorr(long.start)) %in%
                                                   id.name)])
    if ("noint" %in% long.rand.ind) {
      long.rand.ind.temp <- long.rand.ind[-which(long.rand.ind %in%
                                                   "noint")]
    } else {
      long.rand.ind.temp <- long.rand.ind
    }
    rownames(D) <- long.rand.ind.temp
    colnames(D) <- long.rand.ind.temp
    randstart.ind.cov <- lapply(randstart.ind.cov, function(u) {
      rownames(u) <- long.rand.ind
      colnames(u) <- long.rand.ind.temp
      u
    })
  } else if (q == 1) {
    D <- as.data.frame(VarCorr(long.start)[which(names(VarCorr(long.start)) %in%
                                                   id.name)])
    if ("noint" %in% long.rand.ind) {
      long.rand.ind.temp <- long.rand.ind[-which(long.rand.ind %in%
                                                   "noint")]
    } else {
      long.rand.ind.temp <- long.rand.ind
    }
    rownames(D) <- long.rand.ind.temp
    colnames(D) <- long.rand.ind.temp
    randstart.ind.cov <- lapply(randstart.ind.cov, function(u) {
      names(u) <- long.rand.ind.temp
      u
    })
  }
  sigma.e <- as.numeric(as.data.frame(VarCorr(long.start))[nrow(as.data.frame(VarCorr(long.start))),
                                                           5])^2
  ll <- as.numeric(logLik(long.start))
  beta1 <- fixef(long.start)
  out <- list(beta1 = data.frame(beta1), sigma.e = sigma.e, D = D, log.like.long = ll,
              randstart.ind = randstart.ind, randstart.ind.cov = randstart.ind.cov,
              modelfit = long.start)
  if (is.null(long.rand.stud) == FALSE) {
    out$A <- A
    out$randstart.stud <- randstart.stud
    out$randstart.stud.cov <- randstart.stud.cov
  }
  return(out)
}
