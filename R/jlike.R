#' Function to calculate joint likelihood used in the jointmeta1 function
#'
#' @param q the number of individual level random effects
#' @param likeests a list of values required to calculated the log-likelihood
#'   for the fitted joint model.  This list has the following elements:
#'   \describe{
#'
#'   \item{\code{beta1}}{a data frame containing the estimates of the
#'   coefficients of the fixed effect parameters of the longitudinal sub-model.}
#'
#'   \item{\code{beta2}}{a data frame containing the estimates of the
#'   coefficients of the fixed effect parameters of the survival sub-model.}
#'
#'   \item{\code{sigma.e}}{the estimate of the variance for the measurement
#'   errors in the joint model.}
#'
#'   \item{\code{D}}{the estimated covariance matrix for the individual level
#'   random effects}
#'
#'   \item{\code{A}}{the estimated covariance matrix for the study level random
#'   effects.  This is only present in the output if study level random effects
#'   were specified in the function call to \code{jointmeta1}.}
#'
#'   \item{\code{random2}}{a list of matrices containing the conditional modes
#'   of the individual level random effects given the supplied data and the
#'   estimated parameters of the joint model. The list is of length equal to the
#'   number of studies in the dataset, and each element of the list has number
#'   of rows equal to the number of individuals in the study, and number of
#'   columns equal to the number of specified individual level random effects.}
#'
#'   \item{\code{random3}}{a matrix containing the conditional modes of the
#'   study level random effects given the supplied data and the estimated
#'   parameters of the joint model.  The matrix has number of rows equal to the
#'   number of studies, and number of columns equal to the number of specified
#'   study level random effects.}
#'
#'   \item{\code{long.rand.ind.form}}{a character string giving the formulation
#'   of the individual level random effects.}
#'
#'   \item{\code{long.rand.stud.form}}{a character string giving the formulation
#'   of the study level random effects if included in the model.}
#'
#'   \item{\code{conv}}{a logical value indicating whether convergence of the EM
#'   algorithm was achieved or not.}
#'
#'   \item{\code{iters}}{the number of iterations completed by the EM algorithm}
#'
#'   \item{\code{n.bystudy}}{the number of individuals present in each study in
#'   the data supplied to the function.}
#'
#'   \item{\code{haz}}{the estimated baseline hazard.  If \code{strat = TRUE} in
#'   the function call to \code{jointmeta1} then this is a list of length equal
#'   to the number of studies in the supplied dataset, each element of the list
#'   being the baseline hazard for the corresponding study. Otherwise there is a
#'   common baseline across all studies in the dataset and this is one vector.}
#'
#'   \item{\code{rs}}{a counter to indicate the last how many unique event times
#'   had occured by the individual's survival time - this is for use during
#'   further calculation in the joint model EM algorithm.  If a stratified
#'   baseline this is a list of numerical vectors, whereas if the baseline is
#'   not stratified this is a single numeric vector.}
#'
#'   \item{\code{sf}}{the unique event times observed in the dataset. If a
#'   stratified baseline this is a list of numerical vectors, whereas if the
#'   baseline is not stratified this is a single numeric vector.}
#'
#'   }
#' @param randstart.ind a list of the conditional modes of the individual level
#'   random effects in each study given the data and the estimates of the
#'   separate longitudinal model parameters
#' @param randstart.ind.cov a list of the conditional covariance matrices for
#'   each individual for the individual level random effects given the data and
#'   the estimates of the separate longitudinal model parameters
#' @param r the number of study level random effects (if included in the model)
#' @param randstart.stud a data frame containing the conditional modes of the
#'   study level random effects given the data and the estimates of the separate
#'   longitudinal model parameters.  This is only present if study level random
#'   effects were specified in the \code{jointmeta1} function call.
#' @param randstart.stud.cov a list of conditional covariance matrices for each
#'   study for the study level random effects given the data and the estimates
#'   of the separate longitudinal model parameters. This is only present if
#'   study level random effects were specified in the \code{jointmeta1} function
#'   call.
#' @inheritParams jointmeta1
#' @inheritParams EMalgRandprop
#'
#' @return A list containing three elements: \describe{
#'
#'   \item{\code{log.like}}{the overall log-likelihood for the fitted joint
#'   model. }
#'
#'   \item{\code{longlog.like}}{the portion of the log-likelihood attributable
#'   to the longitudinal sub-model.}
#'
#'   \item{\code{survlog.like}}{the portion of the log-likelihood attributable
#'   to the survival sub-model.}
#'
#'   }
#'
#' @importFrom Matrix bdiag
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#'
#' @keywords internal
jlike <- function(data, longdat, survdat, q, likeests, lgpt, studies, p1,
                  p2, long.rand.ind, randstart.ind, randstart.ind.cov, r = NULL, long.rand.stud = NULL,
                  randstart.stud = NULL, randstart.stud.cov = NULL, strat, study.name,
                  id.name) {
  numstudies <- length(studies)
  ids.surv <- survdat[, 1]
  ids.bystudy <- lapply(1:numstudies, function(u) {
    ids.surv[which(survdat[, 4] == studies[u])]
  })
  names(ids.bystudy) <- studies
  id.long <- longdat[, 1]
  id.long.bystudy <- lapply(1:numstudies, function(u) {
    id.long[which(id.long %in% ids.bystudy[[u]])]
  })
  Y <- longdat[, 2]
  Y.bystudy <- lapply(1:numstudies, function(u) {
    Y[id.long %in% ids.bystudy[[u]]]
  })
  tt <- longdat[, 3]
  tt.bystudy <- lapply(1:numstudies, function(u) {
    tt[id.long %in% ids.bystudy[[u]]]
  })
  X1 <- as.matrix(longdat[, 5:ncol(longdat)])
  X1.bystudy <- lapply(1:numstudies, function(u) {
    X1[id.long %in% ids.bystudy[[u]], ]
  })
  n <- nrow(survdat)
  n.bystudy <- sapply(ids.bystudy, length)
  s <- survdat[, 2]
  s.bystudy <- lapply(1:numstudies, function(u) {
    s[ids.surv %in% ids.bystudy[[u]]]
  })
  cen <- survdat[, 3]
  cen.bystudy <- lapply(1:numstudies, function(u) {
    cen[ids.surv %in% ids.bystudy[[u]]]
  })
  nn <- lapply(1:numstudies, function(u) {
    nn <- diff(match(unique(id.long.bystudy[[u]]), id.long.bystudy[[u]]))
    nn <- c(nn, length(id.long.bystudy[[u]]) - sum(nn))
  })
  X2 <- rep(0, numstudies)
  if (p2 > 0) {
    X2 <- as.matrix(survdat[, 5:dim(survdat)[2]])
    X2.bystudy <- lapply(1:numstudies, function(u) {
      as.matrix(X2[ids.surv %in% ids.bystudy[[u]], ])
    })
  } else {
    beta2x <- matrix(0, n, 1)
    beta2x.bystudy <- lapply(1:numstudies, function(u) {
      as.matrix(rep(0, n.bystudy[[u]]))
    })
  }
  beta1 <- likeests$beta1[, 1]
  sigma.e <- likeests$sigma.e
  beta2 <- likeests$beta2[, 1]
  haz <- likeests$haz
  sf <- likeests$sf
  rs <- likeests$rs
  if (strat) {
    rs <- rs[match(names(ids.bystudy), names(rs))]
    sf <- sf[match(names(ids.bystudy), names(sf))]
  }
  N <- sum(do.call(sum, nn))
  time.long <- colnames(longdat)[3]
  D <- likeests$D
  g <- gauss.quad.prob(lgpt, "normal", sigma = sqrt(0.5))
  ab <- g$nodes
  w <- g$weights * sqrt(pi)
  resid <- lapply(1:numstudies, function(u) {
    Y.bystudy[[u]] - X1.bystudy[[u]] %*% beta1
  })
  if ("(Intercept)" %in% long.rand.ind) {
    long.rand.ind2 <- long.rand.ind
    long.rand.ind2[which(long.rand.ind2 == "(Intercept)")] <- "1"
    long.rand.ind.form <- paste(long.rand.ind2, collapse = " + ")
  }
  if ("noint" %in% long.rand.ind) {
    long.rand.ind2 <- long.rand.ind[-which(long.rand.ind == "noint")]
    long.rand.ind.form <- paste("-1", paste(long.rand.ind2, collapse = " + "),
                                sep = " + ")
  }
  Z2.form <- as.formula(paste("~", long.rand.ind.form, sep = ""))
  Z2.frame <- model.frame(Z2.form, data = longdat)
  tZ2 <- model.matrix(Z2.form, Z2.frame)
  Z2 <- t(tZ2)
  tZ2.bystudy <- lapply(1:numstudies, function(u) {
    tZ2[which(id.long %in% ids.bystudy[[u]]), ]
  })
  Z2.bystudy <- lapply(1:numstudies, function(u) {
    Z2[, which(id.long %in% ids.bystudy[[u]])]
  })
  s1.2 <- rep(1:(q - 1), (q - 1):1)
  s2.2 <- sequence((q - 1):1) + rep(1:(q - 1), (q - 1):1)
  Z2.form.surv <- as.formula(paste(gsub(time.long, names(survdat)[2], Z2.form),collapse=""))
  survdat$survdatorder <- 1:nrow(survdat)
  tempsurv <- survdat[, which(!(names(survdat) %in% names(longdat)))]
  tempsurv <- cbind(survdat[, 1], tempsurv)
  names(tempsurv)[1] <- names(survdat)[1]
  tempdat <- merge(tempsurv, longdat[match(unique(longdat[, which(names(longdat) %in%
                                                                    id.name)]), longdat[, which(names(longdat) %in% id.name)]), c(1,
                                                                                                                                  5:ncol(longdat))], by = id.name)
  tempdat <- tempdat[order(tempdat$survdatorder), ]
  tempdat <- tempdat[, -which(names(tempdat) %in% "survdatorder")]
  survdat <- survdat[, -which(names(survdat) %in% "survdatorder")]
  Z2.frame.surv <- model.frame(Z2.form.surv, data = tempdat)
  Z2.surv <- t(model.matrix(Z2.form.surv, Z2.frame.surv))
  Z2.surv.bystudy <- lapply(1:numstudies, function(u) {
    Z2.surv[, which(ids.surv %in% ids.bystudy[[u]])]
  })
  Z2.event.bystudy <- lapply(1:numstudies, function(u) {
    lapply(1:n.bystudy[[u]], function(v) {
      if (strat) {
        rstemp <- rs[[u]][which(names(rs[[u]]) %in% ids.bystudy[[u]][[v]])]
      } else {
        rstemp <- rs[which(names(rs) %in% ids.bystudy[[u]][[v]])]
      }
      if (rstemp > 0) {
        if (q > 1) {
          temp <- matrix(rep(tZ2.bystudy[[u]][v, ], rstemp), ncol = q,
                         byrow = TRUE)
        } else {
          temp <- matrix(rep(tZ2.bystudy[[u]][v], rstemp), ncol = q,
                         byrow = TRUE)
        }
        colnames(temp) <- colnames(tZ2)
        if (time.long %in% colnames(temp)) {
          if (strat) {
            temp[, colnames(temp) == time.long] <- sf[[u]][1:rstemp]
          } else {
            temp[, colnames(temp) == time.long] <- sf[1:rstemp]
          }
        }
        t(temp)
      } else {
        matrix(rep(0, q))
      }
    })
  })
  cnn.bystudy <- lapply(1:numstudies, function(u) {
    c(0, cumsum(nn[[u]]))
  })
  if (p2 == 0) {
    beta2x <- lapply(1:numstudies, function(u) {
      matrix(0, n.bystudy[[u]], 1)
    })
  } else {
    beta2x <- lapply(1:numstudies, function(u) {
      X2.bystudy[[u]] %*% beta2[1:p2]
    })
  }
  sigma.ei <- lapply(1:numstudies, function(u) {
    sigma.e * diag(max(nn[[u]]))
  })
  cov2 <- lapply(1:numstudies, function(u) {
    D %*% Z2.bystudy[[u]]
  })
  randstart.ind.bystudy <- likeests$random2
  l1 <- 0
  l2 <- 0
  if (is.null(r)) {
    gmat <- matrix(0, lgpt^(q), q)
    gmat[, 1] <- rep(ab, each = lgpt^(q - 1))
    if (q > 1) {
      gmat[, 2] <- rep(ab, lgpt)
      w <- as.vector(w %x% w)
      if (q > 2) {
        gmat[, 3] <- rep(ab, each = lgpt)
        w <- as.vector(w %x% g$weights * sqrt(pi))
      }
    }
    newu.all <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[u], function(v) {
        eig <- eigen(D)
        if (length(which(eig$values < 0)) == 0) {
          if (q > 1) {
            B2 <- eig$vectors %*% diag(sqrt(eig$values))
          } else {
            B2 <- eig$vectors %*% matrix(sqrt(eig$values))
          }
          gmat %*% B2 + rep(randstart.ind.bystudy[[u]][v, ], each = nrow(gmat))
        } else {
          warning(paste("Singular matrix encountered during pseudo adaptive
                        procedure, likelihood estimation"))
          gmat + rep(randstart.ind.bystudy[[u]][v, ], each = nrow(gmat))
        }

      })
      })
    lr <- 0
  } else {
    A <- likeests$A
    gmat <- matrix(0, lgpt^(q + r), q + r)
    counter <- rep(1, q + r)
    counter1 <- 1
    while (counter[1] <= lgpt) {
      gmat[counter1, ] <- ab[counter]
      counter1 <- counter1 + 1
      breakout <- FALSE
      countvar <- (q + r)
      while (breakout == FALSE) {
        while (countvar > 0) {
          if (counter[countvar] < lgpt) {
            counter[countvar] <- counter[countvar] + 1
            breakout <- TRUE
            countvar <- 0
          } else {
            counter[countvar] <- counter[countvar] + 1
            for (ii in (q + r):2) {
              if (counter[ii] > lgpt) {
                counter[ii] <- 1
                counter[ii - 1] <- counter[ii - 1] + 1
              }
            }
            breakout <- TRUE
            countvar <- 0
          }
        }
      }
    }
    for (count in 2:(r + q)) {
      w <- as.vector(w %x% g$weights * sqrt(pi))
    }
    lr <- 0.5 * r * sum(numstudies) * log(pi)
    if (study.name %in% long.rand.stud) {
      long.rand.stud2 <- long.rand.stud
      long.rand.stud2[which(long.rand.stud2 == study.name)] <- "1"
      long.rand.stud.form <- paste(long.rand.stud2, collapse = " + ")
    } else {
      long.rand.stud2 <- long.rand.stud
      long.rand.stud.form <- paste("-1", paste(long.rand.stud2, collapse = " + "),
                                   sep = " + ")
    }
    Z3.form <- as.formula(paste("~", long.rand.stud.form, sep = ""))
    Z3.frame <- model.frame(Z3.form, data = longdat)
    tZ3 <- model.matrix(Z3.form, Z3.frame)
    Z3 <- t(tZ3)
    tZ3.bystudy <- lapply(1:numstudies, function(u) {
      tZ3[which(id.long %in% ids.bystudy[[u]]), ]
    })
    Z3.bystudy <- lapply(1:numstudies, function(u) {
      Z3[, which(id.long %in% ids.bystudy[[u]])]
    })
    s1.3 <- rep(1:(r - 1), (r - 1):1)
    s2.3 <- sequence((r - 1):1) + rep(1:(r - 1), (r - 1):1)
    Z3.form.surv <- as.formula(paste(gsub(time.long, names(survdat)[2], Z3.form),collapse=""))
    Z3.frame.surv <- model.frame(Z3.form.surv, data = tempdat)
    Z3.surv <- t(model.matrix(Z3.form.surv, Z3.frame.surv))
    Z3.surv.bystudy <- lapply(1:numstudies, function(u) {
      Z3.surv[, which(ids.surv %in% ids.bystudy[[u]])]
    })
    Z3.event.bystudy <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[[u]], function(v) {
        if (strat) {
          rstemp <- rs[[u]][which(names(rs[[u]]) %in% ids.bystudy[[u]][[v]])]
        } else {
          rstemp <- rs[which(names(rs) %in% ids.bystudy[[u]][[v]])]
        }
        if (rstemp > 0) {
          if (r > 1) {
            temp <- matrix(rep(tZ3.bystudy[[u]][v, ], rstemp),
                           ncol = r, byrow = TRUE)
          } else {
            temp <- matrix(rep(tZ3.bystudy[[u]][v], rstemp), ncol = r,
                           byrow = TRUE)
          }
          colnames(temp) <- colnames(tZ3)
          if (time.long %in% colnames(temp)) {
            if (strat) {
              temp[, colnames(temp) == time.long] <- sf[[u]][1:rstemp]
            } else {
              temp[, colnames(temp) == time.long] <- sf[1:rstemp]
            }
          }
          t(temp)
        } else {
          matrix(rep(0, r))
        }
      })
    })
    cov3 <- lapply(1:numstudies, function(u) {
      A %*% Z3.bystudy[[u]]
    })
    randstart.stud <- likeests$random3
    eig2 <- eigen(D)
    eig3 <- eigen(A)
    if (length(which(eig2$values < 0)) == 0 && length(which(eig3$values <
                                                            0)) == 0) {
      if (q > 1) {
        B2 <- eig2$vectors %*% diag(sqrt(eig2$values))
      } else {
        B2 <- eig2$vectors %*% matrix(sqrt(eig2$values))
      }
      if (r > 1) {
        B3 <- eig3$vectors %*% diag(sqrt(eig3$values))
      } else {
        B3 <- eig3$vectors %*% matrix(sqrt(eig3$values))
      }
      B <- bdiag(B3, B2)
      newu.all <- lapply(1:numstudies, function(u) {
        lapply(1:n.bystudy[u], function(v) {
          gmat %*% B + rep(c(randstart.stud[u, ], randstart.ind.bystudy[[u]][v,
                                                                             ]), each = nrow(gmat))
        })
      })
    } else {
      warning(paste("Singular matrix encountered during pseudo adaptive
                    procedure, likelihood estimation"))
      newu.all <- lapply(1:numstudies, function(u) {
        lapply(1:n.bystudy[u], function(v) {
          gmat + rep(c(randstart.stud[u, ], randstart.ind.bystudy[[u]][v,
                                                                       ]), each = nrow(gmat))
        })
      })
    }



  }
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  counter <- 1
  for (k in 1:numstudies) {
    for (i in 1:n.bystudy[k]) {
      if (strat) {
        rstemp <- rs[[k]][which(names(rs[[k]]) %in% ids.bystudy[[k]][[i]])]
        haz.i <- haz[[k]][rstemp]
        haz.i.rs <- haz[[k]][1:rstemp]
      } else {
        rstemp <- rs[which(names(rs) %in% ids.bystudy[[k]][[i]])]
        haz.i <- haz[rstemp]
        haz.i.rs <- haz[1:rstemp]
      }
      newu <- newu.all[[k]][[i]]
      if (q > 1) {
        tZ2.i <- tZ2.bystudy[[k]][(cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                               1], ]
      } else {
        tZ2.i <- tZ2.bystudy[[k]][(cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                               1]]
      }
      U21.2 <- cov2[[k]][, (cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                        1]]
      Z3iU31 <- matrix(0, nrow = nn[[k]][[i]], ncol = nn[[k]][[i]])
      if (is.null(r) == FALSE) {
        U31.2 <- cov3[[k]][, (cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                          1]]
        if (r > 1) {
          tZ3.i <- tZ3.bystudy[[k]][(cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                                 1], ]
          Z3iU31 <- tZ3.i %*% U31.2
        } else {
          tZ3.i <- tZ3.bystudy[[k]][(cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                                 1]]
          Z3iU31 <- tcrossprod(tZ3.i, U31.2)
        }
        if (q == 1) {
          expZ2u2 <- exp(beta2[p2 + 1] * (newu[, (r + 1):(r + q)] *
                                            Z2.surv.bystudy[[k]][i]))
        } else {
          expZ2u2 <- exp(beta2[p2 + 1] * (newu[, (r + 1):(r + q)] %*%
                                            Z2.surv.bystudy[[k]][, i]))
        }
        if (r == 1) {
          expZ3u3 <- exp(beta2[p2 + 2] * (newu[, 1:r] * Z3.surv.bystudy[[k]][i]))
        } else {
          expZ3u3 <- exp(beta2[p2 + 2] * (newu[, 1:r] %*% Z3.surv.bystudy[[k]][,
                                                                               i]))
        }
        expZu <- expZ2u2 * expZ3u3
        expZu.event <- exp(beta2[p2 + 1] * (newu[, (r + 1):(r +
                                                              q)] %*% Z2.event.bystudy[[k]][[i]])) * exp(beta2[p2 +
                                                                                                                 2] * (newu[, 1:r] %*% Z3.event.bystudy[[k]][[i]]))
      } else {
        if (q == 1) {
          expZu <- exp(beta2[p2 + 1] * (newu[, 1:q] * Z2.surv.bystudy[[k]][i]))
        } else {
          expZu <- exp(beta2[p2 + 1] * (newu[, 1:q] %*% Z2.surv.bystudy[[k]][,
                                                                             i]))
        }
        expZu.event <- exp(beta2[p2 + 1] * (newu[, 1:q] %*% Z2.event.bystudy[[k]][[i]]))
      }
      if (q == 1) {
        U11.2 <- tcrossprod(tZ2.i, U21.2) + Z3iU31 + sigma.ei[[k]][1:nn[[k]][i],
                                                                   1:nn[[k]][i]]
      } else {
        U11.2 <- tZ2.i %*% U21.2 + Z3iU31 + sigma.ei[[k]][1:nn[[k]][i],
                                                          1:nn[[k]][i]]
      }
      if (rstemp > 0) {
        den <- sum(((haz.i * exp(beta2x[[k]][i, ]) * expZu)^cen.bystudy[[k]][[i]]) *
                     (exp(-((exp(beta2x[[k]][i, ]) * expZu.event) %*% haz.i.rs))) *
                     w)
      } else {
        den <- sum((exp(-((exp(beta2x[[k]][i, ]) * expZu.event) %*%
                            matrix(0)))) * w)
      }
      resid.i <- resid[[k]][(cnn.bystudy[[k]][i] + 1):cnn.bystudy[[k]][i +
                                                                         1]]
      l2 <- l2 + 0
      if (den > 0) {
        l2 <- l2 + log(den)
      }
      l1 <- l1 - nn[[k]][i] * 0.5 * log(2 * pi) - 0.5 * log(det(U11.2)) -
        0.5 * sum(resid.i * solve(U11.2, resid.i))
      setTxtProgressBar(pb, counter)
      counter <- counter + 1
    }
  }
  close(pb)
  ll <- l1 + l2 - 0.5 * q * sum(n.bystudy) * log(pi) - lr
  list(log.like = ll, longlog.like = l1, survlog.like = ll - l1)
}
