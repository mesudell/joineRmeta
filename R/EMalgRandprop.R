#' EM algorithm function used in jointmeta1
#'
#' Function to run EM algorithm during one stage model fit.  Used when the
#' jointmeta1 function is called.
#'
#' @param data the original \code{jointdata} as supplied to the
#'   \code{jointmeta1}] call
#' @param longdat the longitudinal data with factors and interaction terms
#'   expanded, ordered by increasing survival time
#' @param survdat the survival data with factors and interaction terms expanded,
#'   ordered by increasing survival time
#' @param id.name character string specifying the id variable in the dataset
#' @param time.long the name of the variable holding the longitudinal time
#'   covariate
#' @param paraests a list of the estimates present from the separate
#'   longitudinal and survival fits.  Same structure as \code{sepests} if
#'   requested in a \code{\link{jointmeta1.object}}
#' @param studies the names of the studies present in the supplied data
#' @param p1 the number of fixed effects included in the longitudinal sub-model
#' @param p2 the number of fixed effects included in the survival sub-model
#' @param q the number of individual level random effects
#' @param r the number of study level random effects, set to \code{NULL} if no
#'   study level random effects included in the model
#' @inheritParams jointmeta1
#' @inheritParams longst
#'
#' @return This function returns a list of the estimates of parameters and other
#'   information from the run of the EM algorithm.  The list has the following
#'   components: \describe{
#'
#'   \item{\code{beta1}}{a data frame containing the estimates of the fixed
#'   effect parameters from the longitudinal sub-model.}
#'
#'   \item{\code{beta2}}{a data frame containing the estimates of the fixed
#'   effect parameters from the survival sub-model.}
#'
#'   \item{\code{sigma.e}}{the estimate of the variance of the measurement
#'   errors.}
#'
#'   \item{\code{haz}}{the estimated baseline hazard.  If \code{strat = TRUE} in
#'   the function call to \code{jointmeta1} then this is a list of length equal
#'   to the number of studies in the supplied dataset, each element of the list
#'   being the baseline hazard for the corresponding study. Otherwise there is a
#'   common baseline across all studies in the dataset and this is one vector.}
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
#'   the data supplied to the function.} }
#'
#' @seealso \code{\link{jointmeta1}}, \code{\link{tojointdata}},
#'   \code{\link[joineR]{jointdata}},\code{\link[lme4]{lmer}},
#'   \code{\link[coxph]{survival}}
#'
#' @keywords internal
#'
#' @import survival statmod
#'

EMalgRandprop <- function(data, longdat, survdat, long.rand.ind, long.rand.stud = NULL,
                          id.name, study.name, gpt, max.it, tol, time.long, surv.formula, long.formula,
                          long.formula.orig, paraests, studies, p1, p2, strat, print.detail,
                          bootrun = FALSE, q, r = NULL) {
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
  X1 <- as.matrix(longdat[, 5:ncol(longdat)])
  X1.bystudy <- lapply(1:numstudies, function(u) {
    X1[id.long %in% ids.bystudy[[u]], ]
  })
  n <- nrow(survdat)
  n.bystudy <- sapply(ids.bystudy, length)
  names(n.bystudy) <- studies
  s <- survdat[, 2]
  s.bystudy <- lapply(1:numstudies, function(u) {
    s[ids.surv %in% ids.bystudy[[u]]]
  })
  names(s.bystudy) <- studies
  cen <- survdat[, 3]
  cen.bystudy <- lapply(1:numstudies, function(u) {
    cen[ids.surv %in% ids.bystudy[[u]]]
  })
  names(cen.bystudy) <- studies
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
  beta1 <- paraests$beta1[, 1]
  names(beta1) <- rownames(paraests$beta1)
  D <- as.matrix(paraests$D)
  sigma.e <- paraests$sigma.e
  if (is.null(long.rand.stud) == FALSE) {
    beta2 <- c(paraests$beta2, 0, 0)
  } else {
    beta2 <- c(paraests$beta2, 0)
  }
  haz <- paraests$haz
  sf <- paraests$sf
  rs <- paraests$rs
  nev <- paraests$nev
  if (strat) {
    rs <- rs[match(names(ids.bystudy), names(rs))]
    sf <- sf[match(names(ids.bystudy), names(sf))]
    haz <- haz[match(names(ids.bystudy), names(haz))]
    nev <- nev[match(names(ids.bystudy), names(nev))]
  }
  nn <- lapply(1:numstudies, function(u) {
    nn <- diff(match(unique(id.long.bystudy[[u]]), id.long.bystudy[[u]]))
    nn <- c(nn, length(id.long.bystudy[[u]]) - sum(nn))
  })
  N <- sum(do.call(sum, nn))
  g <- gauss.quad.prob(gpt, "normal", sigma = sqrt(0.5))
  ab <- g$nodes
  w2 <- g$weights * sqrt(pi)
  gmat2 <- matrix(0, gpt^q, q)
  gmat2[, 1] <- rep(ab, each = gpt^(q - 1))
  if (q > 1) {
    gmat2[, 2] <- rep(ab, gpt)
    w2 <- as.vector(w2 %x% w2)
    if (q > 2) {
      gmat2[, 3] <- rep(ab, each = gpt)
      w2 <- as.vector(w2 %x% g$weights * sqrt(pi))
    }
  }
  EU.2 <- lapply(1:numstudies, function(u) {
    matrix(0, n.bystudy[[u]], q)
  })
  EUU.2 <- lapply(1:numstudies, function(u) {
    matrix(0, n.bystudy[[u]], sum(1:q))
  })
  EexpU.2 <- lapply(1:numstudies, function(u) {
    if (strat) {
      matrix(0, n.bystudy[[u]], length(haz[[u]]))
    } else {
      matrix(0, n.bystudy[[u]], length(haz))
    }
  })
  if ("(Intercept)" %in% long.rand.ind) {
    long.rand.ind2 <- long.rand.ind
    long.rand.ind2[which(long.rand.ind2 == "(Intercept)")] <- "1"
    long.rand.ind.form <- paste(long.rand.ind2, collapse = "+")
  }
  if ("noint" %in% long.rand.ind) {
    long.rand.ind2 <- long.rand.ind[-which(long.rand.ind == "noint")]
    long.rand.ind.form <- paste("-1", paste(long.rand.ind2, collapse = "+"),
                                sep = "+")
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
  Z2.form.surv <- as.formula(gsub(time.long, names(survdat)[2], Z2.form))
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
  tZ2.surv <- model.matrix(Z2.form.surv, Z2.frame.surv)
  tZ2.surv.bystudy <- lapply(1:numstudies, function(u) {
    tZ2.surv[which(ids.surv %in% ids.bystudy[[u]]), ]
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
  Z2.event.sqr.bystudy <- lapply(1:numstudies, function(u) {
    lapply(1:n.bystudy[u], function(v) {
      if (strat) {
        rstemp <- rs[[u]][which(names(rs[[u]]) %in% ids.bystudy[[u]][[v]])]
      } else {
        rstemp <- rs[which(names(rs) %in% ids.bystudy[[u]][[v]])]
      }
      if (rstemp > 0) {
        temp <- (Z2.event.bystudy[[u]][[v]])^2
        if (q > 1) {
          if (rstemp == 1) {
            temp <- rbind(temp, matrix(Z2.event.bystudy[[u]][[v]][s1.2,
                                                                  ] * Z2.event.bystudy[[u]][[v]][s2.2, ]))
          } else {
            temp <- rbind(temp, Z2.event.bystudy[[u]][[v]][s1.2,
                                                           ] * Z2.event.bystudy[[u]][[v]][s2.2, ])
          }
        }
        temp
      } else {
        matrix(rep(0, sum(1:q)))
      }
    })
  })
  Z2b2event <- lapply(1:numstudies, function(u) {
    lapply(1:n.bystudy[u], function(v) {
      if (strat) {
        matrix(0, nrow = length(haz[[u]]), ncol = q)
      } else {
        matrix(0, nrow = length(haz), ncol = q)
      }
    })
  })
  Z2b2eventsqr <- lapply(1:numstudies, function(u) {
    lapply(1:n.bystudy[u], function(v) {
      if (strat) {
        matrix(0, nrow = length(haz[[u]]), ncol = sum(1:q))
      } else {
        matrix(0, nrow = length(haz), ncol = sum(1:q))
      }
    })
  })
  newu.2.all <- lapply(1:numstudies, function(u) {
    lapply(1:n.bystudy[u], function(v) {
      idtemp <- ids.bystudy[[u]][v]
      B2temp <- paraests$randstart.ind.cov[[which(names(paraests$randstart.ind.cov) %in%
                                                    idtemp)]]
      eig <- eigen(B2temp)

      cm.2 <- matrix(as.numeric(as.character(paraests$randstart.ind[which(rownames(paraests$randstart.ind) %in%
                                                                            idtemp), ])), gpt^q, q, TRUE)
      if (length(which(eig$values < 0)) == 0) {
        if (q > 1) {
          B2 <- eig$vectors %*% diag(sqrt(eig$values))
        } else {
          B2 <- eig$vectors %*% matrix(sqrt(eig$values))
        }
        gmat2 %*% B2 + cm.2
      } else {
        warning(paste("Singular matrix encountered during pseudo adaptive
                      procedure, individual level"))
        gmat2 + cm.2
      }
    })
    })
  randstart.ind <- paraests$randstart.ind
  randstart.ind.bystudy <- lapply(1:numstudies, function(u) {
    out <- as.matrix(randstart.ind[which(rownames(randstart.ind) %in%
                                           ids.bystudy[[u]]), ])
    rownames(out) <- rownames(randstart.ind)[which(rownames(randstart.ind) %in%
                                                     ids.bystudy[[u]])]
    out
  })
  randstart.ind.bystudy <- lapply(1:numstudies, function(u) {
    randstart.ind.bystudy[[u]][match(ids.bystudy[[u]], rownames(randstart.ind.bystudy[[u]])),
                               ]
  })
  if (is.null(long.rand.stud) == FALSE) {
    A <- as.matrix(paraests$A)
    w3 <- g$weights * sqrt(pi)
    gmat3 <- matrix(0, gpt^r, r)
    gmat3[, 1] <- rep(ab, each = gpt^(r - 1))
    if (r > 1) {
      gmat3[, 2] <- rep(ab, gpt)
      w3 <- as.vector(w3 %x% w3)
      if (r > 2) {
        gmat3[, 3] <- rep(ab, each = gpt)
        w3 <- as.vector(w3 %x% g$weights * sqrt(pi))
      }
    }
    EU.3 <- matrix(0, numstudies, r)
    EUU.3 <- matrix(0, numstudies, sum(1:r))
    EexpU.3 <- lapply(1:numstudies, function(u) {
      if (strat) {
        matrix(0, n.bystudy[[u]], length(haz[[u]]))
      } else {
        matrix(0, n.bystudy[[u]], length(haz))
      }
    })
    if (study.name %in% long.rand.stud) {
      long.rand.stud2 <- long.rand.stud
      long.rand.stud2[which(long.rand.stud2 == study.name)] <- "1"
      long.rand.stud.form <- paste(long.rand.stud2, collapse = "+")
    } else {
      long.rand.stud2 <- long.rand.stud
      long.rand.stud.form <- paste("-1", paste(long.rand.stud2, collapse = "+"),
                                   sep = "+")
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
    Z3.form.surv <- as.formula(gsub(time.long, names(survdat)[2], Z3.form))
    Z3.frame.surv <- model.frame(Z3.form.surv, data = tempdat)
    tempdat <- NULL
    tZ3.surv <- model.matrix(Z3.form.surv, Z3.frame.surv)
    Z3.surv <- t(tZ3.surv)
    tZ3.surv.bystudy <- lapply(1:numstudies, function(u) {
      tZ3.surv[which(ids.surv %in% ids.bystudy[[u]]), ]
    })
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
    Z3.event.sqr.bystudy <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[u], function(v) {
        if (strat) {
          rstemp <- rs[[u]][which(names(rs[[u]]) %in% ids.bystudy[[u]][[v]])]
        } else {
          rstemp <- rs[which(names(rs) %in% ids.bystudy[[u]][[v]])]
        }
        if (rstemp > 0) {
          temp <- (Z3.event.bystudy[[u]][[v]])^2
          if (r > 1) {
            if (rstemp == 1) {
              temp <- rbind(temp, matrix(Z3.event.bystudy[[u]][[v]][s1.3] *
                                           Z3.event.bystudy[[u]][[v]][s2.3, ]))
            } else {
              temp <- rbind(temp, Z3.event.bystudy[[u]][[v]][s1.3,
                                                             ] * Z3.event.bystudy[[u]][[v]][s2.3, ])
            }
          }
          temp
        } else {
          matrix(rep(0, sum(1:r)))
        }
      })
    })
    randselector1 <- rep(1:r, each = q)
    randselector2 <- rep(1:q, times = r)
    Z3Z2.event.sqr.bystudy <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[u], function(v) {
        if (strat) {
          rstemp <- rs[[u]][which(names(rs[[u]]) %in% ids.bystudy[[u]][[v]])]
        } else {
          rstemp <- rs[which(names(rs) %in% ids.bystudy[[u]][[v]])]
        }
        if (rstemp > 0) {
          Z3.event.bystudy[[u]][[v]][randselector1, ] * Z2.event.bystudy[[u]][[v]][randselector2,
                                                                                   ]
        } else {
          matrix(rep(0, length(randselector2)))
        }
      })
    })
    Z3b3event <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[u], function(v) {
        if (strat) {
          matrix(0, nrow = length(haz[[u]]), ncol = r)
        } else {
          matrix(0, nrow = length(haz), ncol = r)
        }
      })
    })
    Z3b3eventsqr <- lapply(1:numstudies, function(u) {
      lapply(1:n.bystudy[u], function(v) {
        if (strat) {
          matrix(0, nrow = length(haz[[u]]), ncol = sum(1:r))
        } else {
          matrix(0, nrow = length(haz), ncol = sum(1:r))
        }
      })
    })
    newu.3.all <- lapply(1:numstudies, function(u) {
      studtemp <- studies[u]
      B3temp <- paraests$randstart.stud.cov[[which(names(paraests$randstart.stud.cov) %in%
                                                     studtemp)]]
      eig <- eigen(B3temp)
      cm.3 <- matrix(as.numeric(as.character(paraests$randstart.stud[which(rownames(paraests$randstart.stud) %in%
                                                                             studtemp), ])), gpt^r, r, TRUE)
      if (length(which(eig$values < 0)) == 0) {
        if (r > 1) {
          B3 <- eig$vectors %*% diag(sqrt(eig$values))
        } else {
          B3 <- eig$vectors %*% matrix(sqrt(eig$values))
        }
        gmat3 %*% B3 + cm.3
      } else {
        warning(paste("Singular matrix encountered during pseudo adaptive
                      procedure, study level"))
        gmat3 + cm.3
      }
    })
    names(newu.3.all) <- studies
    randstart.stud <- paraests$randstart.stud
    survnames <- NULL
    if (p2 > 0) {
      survnames <- paste("T.", names(beta2)[1:(length(beta2) - 2)],
                         sep = "")
    }
    paranames <- c(paste("Y.", names(beta1), sep = ""), survnames,
                   "gamma_ind_0", "gamma_stud_0", "sigma.e", paste("D", paste(rep(1:q,
                                                                                  times = q), rep(1:q, each = q), sep = ""), sep = ""), paste("A",
                                                                                                                                              paste(rep(1:r, times = r), rep(1:r, each = r), sep = ""),
                                                                                                                                              sep = ""))
    } else {
      r <- NULL
      survnames <- NULL
      if (p2 > 0) {
        survnames <- paste("T.", names(beta2)[1:(length(beta2) - 1)],
                           sep = "")
      }
      paranames <- c(paste("Y.", names(beta1), sep = ""), survnames,
                     "gamma_ind_0", "sigma.e", paste("D", paste(rep(1:q, times = q),
                                                                rep(1:q, each = q), sep = ""), sep = ""))
    }
  conv <- FALSE
  for (it in 1:max.it) {
    if (p2 > 0) {
      beta2x <- lapply(1:numstudies, function(u) {
        X2.bystudy[[u]] %*% beta2[1:p2]
      })
    } else {
      beta2x <- lapply(1:numstudies, function(u) {
        matrix(rep(0, n.bystudy[u]))
      })
    }
    ebeta2x <- lapply(1:numstudies, function(u) {
      exp(beta2x[[u]])
    })
    for (k in 1:numstudies) {
      if (is.null(long.rand.stud) == FALSE) {
        if (it == 1) {
          b2est <- randstart.ind.bystudy[[k]]
        } else {
          b2est <- EU.2[[k]]
        }
        if (it == 1) {
          b3est <- as.numeric(randstart.stud[k, ])
        } else {
          b3est <- as.numeric(EU.3[k, ])
        }
        newu.3 <- newu.3.all[[k]]
        newu.3sqr <- newu.3^2
        if (r > 1) {
          newu.3sqr <- cbind(newu.3sqr, newu.3[, s1.3] * newu.3[,
                                                                s2.3])
        }
        egDUs.3 <- t(apply(exp(newu.3 %*% (Z3.surv.bystudy[[k]] *
                                             beta2[(p2 + 2)])), 1, function(u) {
                                               u^cen.bystudy[[k]]
                                             }))
        egDUsf.3 <- lapply(1:n.bystudy[k], function(v) {
          exp(newu.3 %*% (Z3.event.bystudy[[k]][[v]] * beta2[(p2 +
                                                                2)]))
        })
        if (q > 1) {
          egDUsf.2 <- lapply(1:n.bystudy[k], function(v) {
            exp(as.numeric(b2est[v, ]) %*% (Z2.event.bystudy[[k]][[v]] *
                                              beta2[(p2 + 1)]))
          })
        } else {
          egDUsf.2 <- lapply(1:n.bystudy[k], function(v) {
            exp(as.numeric(b2est[v]) %*% (Z2.event.bystudy[[k]][[v]] *
                                            beta2[(p2 + 1)]))
          })
        }
        ess.3 <- lapply(1:n.bystudy[k], function(v) {
          if (strat) {
            exp(-(ebeta2x[[k]][v, ] * egDUsf.3[[v]]) %*% ((haz[[k]][1:rs[[k]][which(names(rs[[k]]) %in%
                                                                                      ids.bystudy[[k]][[v]])]] * as.numeric(egDUsf.2[[v]]))))
          } else {
            exp(-(ebeta2x[[k]][v, ] * egDUsf.3[[v]]) %*% ((haz[1:rs[which(names(rs) %in%
                                                                            ids.bystudy[[k]][[v]])]] * as.numeric(egDUsf.2[[v]]))))
          }
        })
        f.3 <- rowSums(do.call(cbind, lapply(1:n.bystudy[k], function(v) {
          rowSums(egDUs.3[, v] * ess.3[[v]] * w3)
        })))
        den.3 <- sum(f.3)
        EU.3[k, 1:r] <- f.3 %*% newu.3/den.3
        EUU.3[k, 1:sum(1:r)] <- f.3 %*% newu.3sqr/den.3
      }
      for (i in 1:n.bystudy[k]) {
        if (strat) {
          rstemp <- rs[[k]][which(names(rs[[k]]) %in% ids.bystudy[[k]][[i]])]
        } else {
          rstemp <- rs[which(names(rs) %in% ids.bystudy[[k]][[i]])]
        }
        newu.2 <- newu.2.all[[k]][[i]]
        newu.2sqr <- newu.2^2
        if (q > 1) {
          newu.2sqr <- cbind(newu.2sqr, newu.2[, s1.2] * newu.2[,
                                                                s2.2])
        }
        egDUs.2.ind <- 1
        egDUs.3.const <- 1
        egDUs.3.ind <- 1
        egDUs.2.const <- 1
        if (cen.bystudy[[k]][i] == 1) {
          if (q > 1) {
            egDUs.2.ind <- exp(newu.2 %*% (tZ2.surv.bystudy[[k]][i,
                                                                 ] * beta2[(p2 + 1)]))
          } else {
            egDUs.2.ind <- exp(newu.2 %*% (tZ2.surv.bystudy[[k]][i] *
                                             beta2[(p2 + 1)]))
          }
        }
        egDUsf.2.ind <- exp(newu.2 %*% (Z2.event.bystudy[[k]][[i]] *
                                          beta2[(p2 + 1)]))
        if (is.null(long.rand.stud) == FALSE) {
          if (cen.bystudy[[k]][i] == 1) {
            if (r > 1) {
              egDUs.3.const <- exp(b3est %*% (tZ3.surv.bystudy[[k]][i,
                                                                    ] * beta2[(p2 + 2)]))
              egDUs.3.ind <- exp(newu.3 %*% (tZ3.surv.bystudy[[k]][i,
                                                                   ] * beta2[(p2 + 2)]))
            } else {
              egDUs.3.const <- exp(b3est %*% (tZ3.surv.bystudy[[k]][i] *
                                                beta2[(p2 + 2)]))
              egDUs.3.ind <- exp(newu.3 %*% (tZ3.surv.bystudy[[k]][i] *
                                               beta2[(p2 + 2)]))
            }
            if (q > 1) {
              egDUs.2.const <- exp(b2est[i, ] %*% (tZ2.surv.bystudy[[k]][i,
                                                                         ] * beta2[(p2 + 1)]))
            } else {
              egDUs.2.const <- exp(b2est[i] %*% (tZ2.surv.bystudy[[k]][i] *
                                                   beta2[(p2 + 1)]))
            }
          }
          if (q > 1) {
            egDUsf.2.const <- exp(b2est[i, ] %*% (Z2.event.bystudy[[k]][[i]] *
                                                    beta2[(p2 + 1)]))
          } else {
            egDUsf.2.const <- exp(b2est[i] %*% (Z2.event.bystudy[[k]][[i]] *
                                                  beta2[(p2 + 1)]))
          }
          egDUsf.3.const <- exp(b3est %*% (Z3.event.bystudy[[k]][[i]] *
                                             beta2[(p2 + 2)]))
          egDUsf.3.ind <- exp(newu.3 %*% (Z3.event.bystudy[[k]][[i]] *
                                            beta2[(p2 + 2)]))
          if (strat) {
            ess.3.ind <- exp(-(ebeta2x[[k]][i, ] * egDUsf.3.ind) %*%
                               t(haz[[k]][1:rstemp] * egDUsf.2.const))
          } else {
            ess.3.ind <- exp(-(ebeta2x[[k]][i, ] * egDUsf.3.ind) %*%
                               t(haz[1:rstemp] * egDUsf.2.const))
          }
          f.3.ind <- egDUs.3.ind * ess.3.ind * w3
          den.3.ind <- sum(f.3.ind)
          C.3 <- egDUsf.3.ind[, 1:rstemp]
          EexpU.3[[k]][i, 1:rstemp] <- f.3.ind[, 1] %*% C.3/den.3.ind
          Z3b3event[[k]][[i]][1:ncol(Z3.event.bystudy[[k]][[i]]),
                              ] <- (t(Z3.event.bystudy[[k]][[i]] * as.numeric(f.3.ind[,
                                                                                      1] %*% newu.3)))/den.3.ind
          Z3b3eventsqr[[k]][[i]][1:ncol(Z3.event.sqr.bystudy[[k]][[i]]),
                                 ] <- (as.numeric(f.3.ind[, 1] %*% newu.3sqr) * Z3.event.sqr.bystudy[[k]][[i]])/den.3.ind
        }
        if (strat) {
          if (is.null(long.rand.stud) == FALSE) {
            ess.2 <- exp(-(ebeta2x[[k]][i, ] * egDUsf.2.ind) %*%
                           t(haz[[k]][1:rstemp] * egDUsf.3.const))
          } else {
            ess.2 <- exp(-(ebeta2x[[k]][i, ] * egDUsf.2.ind) %*%
                           (haz[[k]][1:rstemp]))
          }
        } else {
          if (is.null(long.rand.stud) == FALSE) {
            ess.2 <- exp(-(ebeta2x[[k]][i, ] * egDUsf.2.ind) %*%
                           t(haz[1:rstemp] * egDUsf.3.const))
          } else {
            ess.2 <- exp(-(ebeta2x[[k]][i, ] * egDUsf.2.ind) %*%
                           haz[1:rstemp])
          }
        }
        f.2 <- egDUs.2.ind * ess.2 * w2
        den.2 <- sum(f.2)
        EU.2[[k]][i, 1:q] <- f.2[, 1] %*% newu.2/den.2
        EUU.2[[k]][i, 1:sum(1:q)] <- f.2[, 1] %*% newu.2sqr/den.2
        C.2 <- egDUsf.2.ind[, 1:rstemp]
        EexpU.2[[k]][i, 1:rstemp] <- f.2[, 1] %*% C.2/den.2
        Z2b2event[[k]][[i]][1:ncol(Z2.event.bystudy[[k]][[i]]),
                            ] <- (t(Z2.event.bystudy[[k]][[i]] * as.numeric(f.2[,
                                                                                1] %*% newu.2)))/den.2
        Z2b2eventsqr[[k]][[i]][1:ncol(Z2.event.sqr.bystudy[[k]][[i]]),
                               ] <- (as.numeric(f.2[, 1] %*% newu.2sqr) * Z2.event.sqr.bystudy[[k]][[i]])/den.2
      }
    }
    parac <- data.frame(c(beta1, beta2, sigma.e, D))
    if (is.null(long.rand.stud) == FALSE) {
      parac <- data.frame(c(beta1, beta2, sigma.e, D, A))
    }
    EexpU <- lapply(1:numstudies, function(u) {
      if (is.null(long.rand.stud) == FALSE) {
        EexpU.2[[u]] * EexpU.3[[u]]
      } else {
        EexpU.2[[u]]
      }
    })
    if (strat == FALSE) {
      haz <- nev/colSums(do.call(rbind, EexpU) * do.call(rbind, ebeta2x)[,
                                                                         1])
    } else {
      haz <- lapply(1:numstudies, function(u) {
        nev[[u]]/colSums(EexpU[[u]] * ebeta2x[[u]][, 1])
      })
    }
    EUmat.2 <- lapply(1:numstudies, function(u) {
      apply(EU.2[[u]], 2, rep, nn[[u]])
    })
    EUUmat.2 <- lapply(1:numstudies, function(u) {
      apply(EUU.2[[u]], 2, rep, nn[[u]])
    })
    Ut.2 <- lapply(1:numstudies, function(u) {
      rowSums(EUmat.2[[u]] * tZ2.bystudy[[u]])
    })
    if (is.null(long.rand.stud) == FALSE) {
      EUmat.3 <- lapply(1:numstudies, function(u) {
        matrix(rep(EU.3[u, ], sum(nn[[u]])), ncol = r)
      })
      EUUmat.3 <- lapply(1:numstudies, function(u) {
        matrix(rep(EUU.3[u, ], sum(nn[[u]])), ncol = sum(1:r))
      })
      Ut.3 <- lapply(1:numstudies, function(u) {
        rowSums(EUmat.3[[u]] * tZ3.bystudy[[u]])
      })
      beta1 <- as.numeric(solve(crossprod(do.call(rbind, X1.bystudy)),
                                crossprod(do.call(rbind, X1.bystudy), do.call(c, Y.bystudy) -
                                            (do.call(c, Ut.2) + do.call(c, Ut.3)))))
      if (q == 1) {
        if (r == 1) {
          sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                                X1.bystudy) %*% beta1) + (do.call(c, Z2.bystudy) *
                                                                                            do.call(c, lapply(1:numstudies, function(u) {
                                                                                              rep(EU.2[[u]], nn[[u]])
                                                                                            }))) + (do.call(c, Z3.bystudy) * do.call(c, lapply(1:numstudies,
                                                                                                                                               function(u) {
                                                                                                                                                 EU.3[rep(rep(u, n.bystudy[[u]]), nn[[u]])]
                                                                                                                                               })))))^2)/sum(N)
        } else {
          temp1 <- do.call(rbind, lapply(1:numstudies, function(u) {
            EU.3[rep(rep(u, n.bystudy[[u]]), nn[[u]]), ]
          }))
          temp2 <- do.call(cbind, Z3.bystudy)
          out2 <- unlist(lapply(1:nrow(temp1), function(u) {
            sum(temp1[u, ] * temp2[, u])
          }))
          sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                                X1.bystudy) %*% beta1) + (do.call(c, Z2.bystudy) *
                                                                                            do.call(c, lapply(1:numstudies, function(u) {
                                                                                              rep(EU.2[[u]], nn[[u]])
                                                                                            }))) + out2))^2)/sum(N)
        }
      } else {
        if (r == 1) {
          temp1 <- do.call(rbind, lapply(1:numstudies, function(u) {
            EU.2[[u]][rep(1:nrow(EU.2[[u]]), nn[[u]]), ]
          }))
          temp2 <- do.call(cbind, Z2.bystudy)
          out1 <- unlist(lapply(1:nrow(temp1), function(u) {
            sum(temp1[u, ] * temp2[, u])
          }))
          sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                                X1.bystudy) %*% beta1) + out1 + (do.call(c, Z3.bystudy) *
                                                                                                   do.call(c, lapply(1:numstudies, function(u) {
                                                                                                     EU.3[rep(rep(u, n.bystudy[[u]]), nn[[u]])]
                                                                                                   })))))^2)/sum(N)
        } else {
          temp1 <- do.call(rbind, lapply(1:numstudies, function(u) {
            EU.2[[u]][rep(1:nrow(EU.2[[u]]), nn[[u]]), ]
          }))
          temp2 <- do.call(cbind, Z2.bystudy)
          out1 <- unlist(lapply(1:nrow(temp1), function(u) {
            sum(temp1[u, ] * temp2[, u])
          }))
          temp1 <- do.call(rbind, lapply(1:numstudies, function(u) {
            EU.3[rep(rep(u, n.bystudy[[u]]), nn[[u]]), ]
          }))
          temp2 <- do.call(cbind, Z3.bystudy)
          out2 <- unlist(lapply(1:nrow(temp1), function(u) {
            sum(temp1[u, ] * temp2[, u])
          }))
          sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                                X1.bystudy) %*% beta1) + out1 + out2))^2)/sum(N)
        }
      }
      diag(A) <- colMeans(EUU.3)[1:r]
      if (r > 1) {
        A[lower.tri(A)] <- colMeans(EUU.3)[-(1:r)]
        A[upper.tri(A)] <- t(A)[upper.tri(A)]
      }
    } else {
      beta1 <- as.numeric(solve(crossprod(do.call(rbind, X1.bystudy)),
                                crossprod(do.call(rbind, X1.bystudy), do.call(c, Y.bystudy) -
                                            (do.call(c, Ut.2)))))
      if (q == 1) {
        sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                              X1.bystudy) %*% beta1) + (do.call(c, Z2.bystudy) * do.call(c,
                                                                                                                         lapply(1:numstudies, function(u) {
                                                                                                                           rep(EU.2[[u]], nn[[u]])
                                                                                                                         })))))^2)/sum(N)
      } else {
        temp1 <- do.call(rbind, lapply(1:numstudies, function(u) {
          EU.2[[u]][rep(1:nrow(EU.2[[u]]), nn[[u]]), ]
        }))
        temp2 <- do.call(cbind, Z2.bystudy)
        out <- unlist(lapply(1:nrow(temp1), function(u) {
          sum(temp1[u, ] * temp2[, u])
        }))
        sigma.e <- colSums((do.call(c, Y.bystudy) - ((do.call(rbind,
                                                              X1.bystudy) %*% beta1) + (out)))^2)/sum(N)


      }
    }
    diag(D) <- colMeans(do.call(rbind, EUU.2))[1:q]
    if (q > 1) {
      D[lower.tri(D)] <- colMeans(do.call(rbind, EUU.2))[-(1:q)]
      D[upper.tri(D)] <- t(D)[upper.tri(D)]
    }
    if (is.null(long.rand.stud) == FALSE) {
      EexpU <- lapply(1:numstudies, function(u) {
        lapply(1:n.bystudy[u], function(v) {
          EexpU.2[[u]][v, ] * EexpU.3[[u]][v, ]
        })
      })
      fd <- vector("numeric", p2 + q + r)
      sd <- matrix(0, p2 + q + r, p2 + q + r)
    } else {
      EexpU <- lapply(1:numstudies, function(u) {
        lapply(1:n.bystudy[u], function(v) {
          EexpU.2[[u]][v, ]
        })
      })
      fd <- vector("numeric", p2 + q)
      sd <- matrix(0, p2 + q, p2 + q)
    }
    if (q > 1) {
      fd[(p2 + 1):(p2 + q)] <- colSums(do.call(c, cen.bystudy) *
                                         (do.call(rbind, EU.2) * do.call(rbind, tZ2.surv.bystudy))) -
        colSums(do.call(rbind, ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies,
                                                                     function(u) {
                                                                       if (strat) {
                                                                         haztemp <- haz[[u]]
                                                                       } else {
                                                                         haztemp <- haz
                                                                       }
                                                                       do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                         colSums(Z2b2event[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                       }))
                                                                     })))
    } else {
      fd[(p2 + 1):(p2 + q)] <- sum(do.call(c, cen.bystudy) * (do.call(rbind,
                                                                      EU.2) * do.call(c, tZ2.surv.bystudy))) - colSums(do.call(rbind,
                                                                                                                               ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                                                                 if (strat) {
                                                                                                                                   haztemp <- haz[[u]]
                                                                                                                                 } else {
                                                                                                                                   haztemp <- haz
                                                                                                                                 }
                                                                                                                                 do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                   colSums(Z2b2event[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                                                                                 }))
                                                                                                                               })))
    }
    if (is.null(long.rand.stud) == FALSE) {
      if (r > 1) {
        fd[(p2 + q + 1):(p2 + q + r)] <- colSums(do.call(c, cen.bystudy) *
                                                   (EU.3[rep(1:nrow(EU.3), n.bystudy), ]) * do.call(rbind,
                                                                                                    tZ3.surv.bystudy)) - colSums(do.call(rbind, ebeta2x)[,
                                                                                                                                                         1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                                                                                           if (strat) {
                                                                                                                                                             haztemp <- haz[[u]]
                                                                                                                                                           } else {
                                                                                                                                                             haztemp <- haz
                                                                                                                                                           }
                                                                                                                                                           do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                                             colSums(Z3b3event[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                                                                                                           }))
                                                                                                                                                         })))
      } else {
        fd[(p2 + q + 1):(p2 + q + r)] <- sum(do.call(c, cen.bystudy) *
                                               (EU.3[rep(1:nrow(EU.3), n.bystudy), ]) * do.call(c, tZ3.surv.bystudy)) -
          colSums(do.call(rbind, ebeta2x)[, 1] * do.call(rbind,
                                                         lapply(1:numstudies, function(u) {
                                                           if (strat) {
                                                             haztemp <- haz[[u]]
                                                           } else {
                                                             haztemp <- haz
                                                           }
                                                           do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                             colSums(Z3b3event[[u]][[v]] * EexpU[[u]][[v]] *
                                                                       haztemp)
                                                           }))
                                                         })))
      }
    }
    if (q > 1) {
      if (q == 2) {
        sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)][upper.tri(sd[(p2 +
                                                                 1):(p2 + q), (p2 + 1):(p2 + q)])] <- -sum(do.call(rbind,
                                                                                                                   ebeta2x)[, 1] * 0.5 * do.call(rbind, lapply(1:numstudies,
                                                                                                                                                               function(u) {
                                                                                                                                                                 if (strat) {
                                                                                                                                                                   haztemp <- haz[[u]]
                                                                                                                                                                 } else {
                                                                                                                                                                   haztemp <- haz
                                                                                                                                                                 }
                                                                                                                                                                 do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                                                   sum(Z2b2eventsqr[[u]][[v]][, (q + 1):sum(1:q)] *
                                                                                                                                                                         EexpU[[u]][[v]] * haztemp)
                                                                                                                                                                 }))
                                                                                                                                                               })))
      } else {
        sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)][upper.tri(sd[(p2 +
                                                                 1):(p2 + q), (p2 + 1):(p2 + q)])] <- -colSums(do.call(rbind,
                                                                                                                       ebeta2x)[, 1] * 0.5 * do.call(rbind, lapply(1:numstudies,
                                                                                                                                                                   function(u) {
                                                                                                                                                                     if (strat) {
                                                                                                                                                                       haztemp <- haz[[u]]
                                                                                                                                                                     } else {
                                                                                                                                                                       haztemp <- haz
                                                                                                                                                                     }
                                                                                                                                                                     do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                                                       sum(Z2b2eventsqr[[u]][[v]][, (q + 1):sum(1:q)] *
                                                                                                                                                                             EexpU[[u]][[v]] * haztemp)
                                                                                                                                                                     }))
                                                                                                                                                                   })))
      }
    }
    if (is.null(long.rand.stud) == FALSE) {
      if (r > 1) {
        if (r == 2) {
          sd[(p2 + q + 1):(p2 + q + r), (p2 + q + 1):(p2 + q +
                                                        r)][upper.tri(sd[(p2 + q + 1):(p2 + q + r), (p2 + q +
                                                                                                       1):(p2 + q + r)])] <- -sum(do.call(rbind, ebeta2x)[,
                                                                                                                                                          1] * 0.5 * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                                                                                            if (strat) {
                                                                                                                                                              haztemp <- haz[[u]]
                                                                                                                                                            } else {
                                                                                                                                                              haztemp <- haz
                                                                                                                                                            }
                                                                                                                                                            do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                                              sum(Z3b3eventsqr[[u]][[v]][, (r + 1):sum(1:r)] *
                                                                                                                                                                    EexpU[[u]][[v]] * haztemp)
                                                                                                                                                            }))
                                                                                                                                                          })))
        } else {
          sd[(p2 + q + 1):(p2 + q + r), (p2 + q + 1):(p2 + q +
                                                        r)][upper.tri(sd[(p2 + q + 1):(p2 + q + r), (p2 + q +
                                                                                                       1):(p2 + q + r)])] <- -colSums(do.call(rbind, ebeta2x)[,
                                                                                                                                                              1] * 0.5 * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                                                                                                if (strat) {
                                                                                                                                                                  haztemp <- haz[[u]]
                                                                                                                                                                } else {
                                                                                                                                                                  haztemp <- haz
                                                                                                                                                                }
                                                                                                                                                                do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                                                  sum(Z3b3eventsqr[[u]][[v]][, (r + 1):sum(1:r)] *
                                                                                                                                                                        EexpU[[u]][[v]] * haztemp)
                                                                                                                                                                }))
                                                                                                                                                              })))
        }
      }
      if (q == 1 && r == 1) {
        sd[(p2 + 1):(p2 + q), ((p2 + q + 1):(p2 + q + r))] <- -colSums(do.call(rbind,
                                                                               ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                 if (strat) {
                                                                                   haztemp <- haz[[u]]
                                                                                 } else {
                                                                                   haztemp <- haz
                                                                                 }
                                                                                 do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                   sum(Z3b3event[[u]][[v]] * Z2b2event[[u]][[v]] * EexpU.2[[u]][v,
                                                                                                                                                ] * EexpU.3[[u]][v, ] * haztemp)
                                                                                 }))
                                                                               })))
      } else {
        sd[(p2 + 1):(p2 + q), ((p2 + q + 1):(p2 + q + r))] <- -colSums(do.call(rbind,
                                                                               ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                 if (strat) {
                                                                                   haztemp <- haz[[u]]
                                                                                 } else {
                                                                                   haztemp <- haz
                                                                                 }
                                                                                 do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                   colSums(Z3b3event[[u]][[v]][, randselector1] * Z2b2event[[u]][[v]][,
                                                                                                                                                      randselector2] * EexpU[[u]][[v]] * haztemp)
                                                                                 }))
                                                                               })))
      }
    }
    if (p2 > 0) {
      fd[1:p2] <- c(colSums((do.call(c, cen.bystudy) * do.call(rbind,
                                                               X2.bystudy)) - (do.call(rbind, X2.bystudy) * as.numeric(do.call(rbind,
                                                                                                                               ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                                                                                 if (strat) {
                                                                                                                                   haztemp <- haz[[u]]
                                                                                                                                 } else {
                                                                                                                                   haztemp <- haz
                                                                                                                                 }
                                                                                                                                 do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                   sum(EexpU[[u]][[v]] * haztemp)
                                                                                                                                 }))
                                                                                                                               }))))))
      sd[(1:p2), (p2 + 1):(p2 + q)] <- -(t(do.call(rbind, X2.bystudy)) %*%
                                           (do.call(rbind, ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies,
                                                                                                 function(u) {
                                                                                                   if (strat) {
                                                                                                     haztemp <- haz[[u]]
                                                                                                   } else {
                                                                                                     haztemp <- haz
                                                                                                   }
                                                                                                   do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                     colSums(Z2b2event[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                                                   }))
                                                                                                 }))))
      if (is.null(long.rand.stud) == FALSE) {
        sd[(1:p2), (p2 + q + 1):(p2 + q + r)] <- -(t(do.call(rbind,
                                                             X2.bystudy)) %*% (do.call(rbind, ebeta2x)[, 1] * do.call(rbind,
                                                                                                                      lapply(1:numstudies, function(u) {
                                                                                                                        if (strat) {
                                                                                                                          haztemp <- haz[[u]]
                                                                                                                        } else {
                                                                                                                          haztemp <- haz
                                                                                                                        }
                                                                                                                        do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                          colSums(Z3b3event[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                                                                        }))
                                                                                                                      }))))
      }
      sd <- sd + t(sd)
      for (i in 1:p2) {
        for (j in 1:p2) {
          sd[i, j] <- -(sum(do.call(rbind, X2.bystudy)[, i] * do.call(rbind,
                                                                      X2.bystudy)[, j] * do.call(rbind, ebeta2x)[, 1] * do.call(rbind,
                                                                                                                                lapply(1:numstudies, function(u) {
                                                                                                                                  if (strat) {
                                                                                                                                    haztemp <- haz[[u]]
                                                                                                                                  } else {
                                                                                                                                    haztemp <- haz
                                                                                                                                  }
                                                                                                                                  do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                                    sum(EexpU[[u]][[v]] * haztemp)
                                                                                                                                  }))
                                                                                                                                }))))
        }
      }
    }
    if (q == 1) {
      sd[(p2 + 1), (p2 + 1)] <- -(colSums(do.call(rbind, ebeta2x)[,
                                                                  1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                    if (strat) {
                                                                      haztemp <- haz[[u]]
                                                                    } else {
                                                                      haztemp <- haz
                                                                    }
                                                                    do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                      colSums(Z2b2eventsqr[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                    }))
                                                                  }))))[1:q]
    } else {
      diag(sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)]) <- -(colSums(do.call(rbind,
                                                                          ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                                            if (strat) {
                                                                              haztemp <- haz[[u]]
                                                                            } else {
                                                                              haztemp <- haz
                                                                            }
                                                                            do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                              colSums(Z2b2eventsqr[[u]][[v]] * EexpU[[u]][[v]] * haztemp)
                                                                            }))
                                                                          }))))[1:q]
    }
    if (is.null(long.rand.stud) == FALSE) {
      if (r == 1) {
        sd[(p2 + q + 1), (p2 + q + 1)] <- -(colSums(do.call(rbind,
                                                            ebeta2x)[, 1] * do.call(rbind, lapply(1:numstudies, function(u) {
                                                              if (strat) {
                                                                haztemp <- haz[[u]]
                                                              } else {
                                                                haztemp <- haz
                                                              }
                                                              do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                colSums(Z3b3eventsqr[[u]][[v]] * EexpU[[u]][[v]] *
                                                                          haztemp)
                                                              }))
                                                            }))))[1:r]
      } else {
        diag(sd[(p2 + q + 1):(p2 + q + r), (p2 + q + 1):(p2 + q +
                                                           r)]) <- -(colSums(do.call(rbind, ebeta2x)[, 1] * do.call(rbind,
                                                                                                                    lapply(1:numstudies, function(u) {
                                                                                                                      if (strat) {
                                                                                                                        haztemp <- haz[[u]]
                                                                                                                      } else {
                                                                                                                        haztemp <- haz
                                                                                                                      }
                                                                                                                      do.call(rbind, lapply(1:n.bystudy[u], function(v) {
                                                                                                                        colSums(Z3b3eventsqr[[u]][[v]] * EexpU[[u]][[v]] *
                                                                                                                                  haztemp)
                                                                                                                      }))
                                                                                                                    }))))[1:r]
      }
    }
    if (is.null(long.rand.stud) == FALSE) {
      if (r == 1 && q == 1) {
        fd <- fd
        sd <- sd
      } else {
        if (p2 > 0) {
          fd <- c(fd[1:p2], sum(fd[(p2 + 1):(p2 + q)]), sum(fd[(p2 +
                                                                  q + 1):(p2 + q + r)]))
          sd.p2 <- sd[1:p2, 1:p2]
          sd.q <- sum(sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)])
          sd.r <- sum(sd[(p2 + q + 1):(p2 + q + r), (p2 + q + 1):(p2 +
                                                                    q + r)])
          if (p2 > 1) {
            if (q > 1) {
              sd.p2q <- rowSums(sd[1:p2, (p2 + 1):(p2 + q)])
            } else {
              sd.p2q <- sd[1:p2, (p2 + 1):(p2 + q)]
            }
            if (r > 1) {
              sd.p2r <- rowSums(sd[1:p2, (p2 + q + 1):(p2 + q +
                                                         r)])
            } else {
              sd.p2r <- sd[1:p2, (p2 + q + 1):(p2 + q + r)]
            }
          } else {
            sd.p2q <- sum(sd[1:p2, (p2 + 1):(p2 + q)])
            sd.p2r <- sum(sd[1:p2, (p2 + q + 1):(p2 + q + r)])
          }
          sd.qr <- sum(sd[(p2 + 1):(p2 + q), (p2 + q + 1):(p2 +
                                                             q + r)])
          sd <- rbind(cbind(sd.p2, sd.p2q, sd.p2r), c(sd.p2q, sd.q,
                                                      sd.qr), c(sd.p2r, sd.qr, sd.r))
          rownames(sd) <- colnames(sd) <- NULL
        } else {
          fd <- c(sum(fd[(p2 + 1):(p2 + q)]), sum(fd[(p2 + q +
                                                        1):(p2 + q + r)]))
          sd.q <- sum(sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)])
          sd.r <- sum(sd[(p2 + q + 1):(p2 + q + r), (p2 + q + 1):(p2 +
                                                                    q + r)])
          sd.qr <- sum(sd[(p2 + 1):(p2 + q), (p2 + q + 1):(p2 +
                                                             q + r)])
          sd <- rbind(c(sd.q, sd.qr), c(sd.qr, sd.r))
          rownames(sd) <- colnames(sd) <- NULL
        }
      }
    } else {
      if (q == 1) {
        fd <- fd
        sd <- sd
      } else {
        if (p2 > 0) {
          fd <- c(fd[1:p2], sum(fd[(p2 + 1):(p2 + q)]))
          sd.p2 <- sd[1:p2, 1:p2]
          sd.q <- sum(sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)])
          if (p2 > 1) {
            sd.p2q <- rowSums(sd[1:p2, (p2 + 1):(p2 + q)])
          } else {
            sd.p2q <- sum(sd[1:p2, (p2 + 1):(p2 + q)])
          }
          sd <- rbind(cbind(sd.p2, sd.p2q), c(sd.p2q, sd.q))
          rownames(sd) <- colnames(sd) <- NULL
        } else {
          fd <- c(sum(fd[(p2 + 1):(p2 + q)]))
          sd <- sum(sd[(p2 + 1):(p2 + q), (p2 + 1):(p2 + q)])
          rownames(sd) <- colnames(sd) <- NULL
        }
      }
    }
    beta2 <- beta2 - solve(sd, fd)
    para <- data.frame(c(beta1, beta2, sigma.e, D))
    if (is.null(long.rand.stud) == FALSE) {
      para <- data.frame(c(beta1, beta2, sigma.e, D, A))
    }
    dd <- abs(parac - para)
    if (print.detail) {
      print(paste("Iteration: ", it, sep = ""))
      print("Current parameter estimates:")
      detail <- data.frame(para)
      rownames(detail) <- paranames
      colnames(detail) <- NULL
      print(detail)
    }
    if (max(dd) < tol) {
      conv <- TRUE
      break
    }
  }
  if (conv != TRUE && bootrun == FALSE) {
    print("Not converged")
  }
  if (strat) {
    names(haz) <- studies
  }
  results <- list(beta1 = data.frame(beta1), beta2 = data.frame(beta2),
                  sigma.e = sigma.e, D = D, haz = haz, random2 = EU.2, conv = conv,
                  iters = it, long.rand.ind.form = long.rand.ind.form, n.bystudy = n.bystudy)
  if (is.null(long.rand.stud) == FALSE) {
    results$A <- A
    results$random3 <- EU.3
    results$long.rand.stud.form <- long.rand.stud.form
  }
  results
}
