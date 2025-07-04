#' Data simulation function for the case with both individual level random
#' effects and study level random effects
#'
#' Internal function to simulate joint data with random effect at the both
#' individual level and the study level.  Used inside
#' \code{\link{simjointmeta}}.
#'
#' @param gamma are the association parameters.  If different association
#'   parameters are supplied for each study in the dataset, this is a list of
#'   vectors each of length equal to the total number of random effects. If the
#'   same association parameters are supplied for each study in the dataset then
#'   this is a vector of length equal to the number of random effects.  If
#'   separate association parameters are defined for different random effects
#'   (i.e. if \code{sepassoc = TRUE}) then the elements in each of these vectors
#'   are not necessarily identical. If \code{sepassoc = FALSE} then the
#'   association parameters within each vector are identical.  In each vector of
#'   association parameters the first \code{q} values are the association
#'   parameters for the individual level random effects, and the remaining
#'   values are the association parameters for the study level random effects.
#' @param q the number of individual level random effects
#' @param r the number of study level random effects
#' @inheritParams simjointmeta
#'
#'
#' @return This function returns a list with three named elements.  The first
#'   element is named \code{'longdat'}, the second \code{'survdat'}, the third
#'   \code{'percentevent'}.  Each of these elements is a list of length equal to
#'   the number of studies specified to simulate in the function call.
#'
#'   The element \code{'longdat'} is a list of the simulated longitudinal data
#'   sets.  Each longitudinal dataset contains the following variables:
#'   \describe{
#'
#'   \item{\code{id}}{a numeric id variable}
#'
#'   \item{\code{Y}}{the continuous longitudinal outcome}
#'
#'   \item{\code{time}}{the numeric longitudinal time variable}
#'
#'   \item{\code{study}}{a study membership variable}
#'
#'   \item{\code{intercept}}{an intercept term}
#'
#'   \item{\code{treat}}{a treatment assignment variable to one of two treatment
#'   groups}
#'
#'   \item{\code{ltime}}{a duplicate of the longitudinal time variable}
#'
#'   }
#'
#'   The element \code{'survdat'} is a list of the simulated survival data sets.
#'   Each survival dataset contains the following variables: \describe{
#'
#'   \item{\code{id}}{a numeric id variable}
#'
#'   \item{\code{survtime}}{the numeric survival times}
#'
#'   \item{\code{cens}}{the censoring indicator}
#'
#'   \item{\code{study}}{a study membership variable}
#'
#'   \item{\code{treat}}{a treatment assignment variable to one of two treatment
#'   groups}
#'
#'   }
#'
#'   The element \code{'percentevent'} is a list of the percentage of events
#'   over censorings seen in the simulated survival data.
#'
#' @keywords internal
#'
#' @importFrom stats runif rnorm rbinom
#'
simdat2randlevels <- function(k, n, rand_ind, rand_stud, sepassoc, ntms,
                              longmeasuretimes, beta1, beta2, gamma, sigb_ind, sigb_stud, vare, theta0,
                              theta1, censoring, censlam, truncation, trunctime, q, r) {
  treat <- lapply(1:k, function(u) {
    rbinom(n[u], 1, 0.5)
  })
  X2 <- lapply(1:k, function(u) {
    cbind(treat = treat[[u]])
  })
  id <- lapply(1:k, function(u) {
    if (u > 1) {
      cs <- sum(n[1:(u - 1)])
    } else {
      cs <- 0
    }
    cs + (1:n[u])
  })
  idl <- lapply(1:k, function(u) {
    rep(id[[u]], each = ntms)
  })
  treatl <- lapply(1:k, function(u) {
    rep(treat[[u]], each = ntms)
  })
  time <- lapply(1:k, function(u) {
    rep(longmeasuretimes, length = n[[u]] * ntms)
  })
  X1 <- lapply(1:k, function(u) {
    cbind(intercept = 1, treat = treatl[[u]], ltime = time[[u]])
  })
  U_ind <- lapply(1:k, function(u) {
    mvrnorm(n[u], mu = rep(0, q), Sigma = sigb_ind)
  })
  U_stud <- mvrnorm(k, mu = rep(0, r), Sigma = sigb_stud)
  U_ind_l <- lapply(1:k, function(u) {
    U_ind[[u]][rep(1:n[u], each = ntms), ]
  })
  U_stud_l <- lapply(1:k, function(u) {
    U_stud[rep(u, each = (ntms * n[u])), ]
  })
  if (rand_ind == "intslope") {
    Z2 <- lapply(1:k, function(u) {
      X1[[u]][, which(colnames(X1[[u]]) %in% c("intercept", "ltime"))]
    })
  } else {
    Z2 <- lapply(1:k, function(u) {
      X1[[u]][, which(colnames(X1[[u]]) %in% c("intercept"))]
    })
  }
  if (rand_stud == "inttreat") {
    Z3 <- lapply(1:k, function(u) {
      X1[[u]][, which(colnames(X1[[u]]) %in% c("intercept", "treat"))]
    })
  } else if (rand_stud == "int") {
    Z3 <- lapply(1:k, function(u) {
      X1[[u]][, which(colnames(X1[[u]]) %in% c("intercept"))]
    })
  } else {
    Z3 <- lapply(1:k, function(u) {
      X1[[u]][, which(colnames(X1[[u]]) %in% c("treat"))]
    })
  }
  Z2b2 <- lapply(1:k, function(u) {
    t((Z2[[u]]) * as.vector(U_ind_l[[u]]))
  })
  Z3b3 <- lapply(1:k, function(u) {
    t((Z3[[u]]) * as.vector(U_stud_l[[u]]))
  })
  Y <- lapply(1:k, function(u) {
    X1[[u]] %*% beta1 + colSums(Z2b2[[u]]) + colSums(Z3b3[[u]]) + sqrt(vare) *
      rnorm(n[u] * ntms)
  })
  if (rand_stud == "inttreat") {
    Z3b3.surv <- lapply(1:k, function(u) {
      t(cbind(intercept = 1, treat = treat[[u]])) * as.vector(U_stud[[u]])
    })
  } else if (rand_stud == "int") {
    Z3b3.surv <- lapply(1:k, function(u) {
      t(cbind(intercept = rep(1, n[u]))) * as.vector(U_stud[[u]])
    })
  } else {
    Z3b3.surv <- lapply(1:k, function(u) {
      t(cbind(treat = treat[[u]])) * as.vector(U_stud[[u]])
    })
  }
  beta2x <- lapply(1:k, function(u) {
    X2[[u]] %*% beta2
  })
  uu <- lapply(1:k, function(u) {
    runif(n[u])
  })
  if (rand_ind == "int") {
    if (inherits(gamma,"list")) {
      survtime <- lapply(1:k, function(u) {
        -log(uu[[u]])/exp(theta0[[u]] + beta2x[[u]] + (U_ind[[u]] *
                                                         gamma[[u]][1]) + colSums(Z3b3.surv[[u]] * gamma[[u]][(q +
                                                                                                                 1):(q + r)]))
      })
    } else {
      survtime <- lapply(1:k, function(u) {
        -log(uu[[u]])/exp(theta0[[u]] + beta2x[[u]] + (U_ind[[u]] *
                                                         gamma[1]) + colSums(Z3b3.surv[[u]] * gamma[(q + 1):(q +
                                                                                                               r)]))
      })
    }
  } else {
    if (inherits(gamma,"list")) {
      survtime <- lapply(1:k, function(u) {
        temp1 <- ((gamma[[u]][2] * U_ind[[u]][, 2]) + theta1[[u]])
        temp2 <- exp(theta0[[u]] + beta2x[[u]] + (gamma[[u]][1] *
                                                    U_ind[[u]][, 1]) + colSums(Z3b3.surv[[u]] * gamma[[u]][(q +
                                                                                                              1):(q + r)]))
        ii <- (temp1 < 0) & (uu[[u]] < exp(temp2/temp1))
        survtime <- rep(0, n[u])
        survtime[ii] <- Inf
        survtime[!ii] <- (1/temp1[!ii]) * log(1 + (((temp1[!ii]) *
                                                      (-log(uu[[u]][!ii])))/temp2[!ii]))
        survtime
      })
    } else {
      survtime <- lapply(1:k, function(u) {
        temp1 <- ((gamma[2] * U_ind[[u]][, 2]) + theta1[[u]])
        temp2 <- exp(theta0[[u]] + beta2x[[u]] + (gamma[1] * U_ind[[u]][,
                                                                        1]) + colSums(Z3b3.surv[[u]] * gamma[(q + 1):(q + r)]))
        ii <- (temp1 < 0) & (uu[[u]] < exp(temp2/temp1))
        survtime <- rep(0, n[u])
        survtime[ii] <- Inf
        survtime[!ii] <- (1/temp1[!ii]) * log(1 + (((temp1[!ii]) *
                                                      (-log(uu[[u]][!ii])))/temp2[!ii]))
        survtime
      })
    }
  }
  if (censoring) {
    censtime <- lapply(1:k, function(u) {
      -log(runif(n[u]))/censlam[u]
    })
  } else {
    censtime <- lapply(1:k, function(u) {
      rep(Inf, n[u])
    })
  }
  if (truncation) {
    censtime <- lapply(1:k, function(u) {
      pmin(censtime[[u]], trunctime)
    })
  }
  ii <- lapply(1:k, function(u) {
    censtime[[u]] < survtime[[u]]
  })
  survtime <- lapply(1:k, function(u) {
    temp <- survtime[[u]]
    temp <- unlist(lapply(1:n[u], function(v) {
      if (censtime[[u]][v] < survtime[[u]][v]) {
        temp[v] <- censtime[[u]][v]
      } else {
        temp[v] <- survtime[[u]][v]
      }
    }))
    temp
  })
  cens <- lapply(1:k, function(u) {
    temp <- rep(1, n[u])
    temp[ii[[u]]] <- 0
    temp
  })
  ls <- lapply(1:k, function(u) {
    rep(survtime[[u]], each = ntms)
  })
  Y <- lapply(1:k, function(u) {
    Y[[u]][ls[[u]] > time[[u]]]
  })
  X1 <- lapply(1:k, function(u) {
    X1[[u]][ls[[u]] > time[[u]], ]
  })
  idl <- lapply(1:k, function(u) {
    idl[[u]][ls[[u]] > time[[u]]]
  })
  time <- lapply(1:k, function(u) {
    time[[u]][ls[[u]] > time[[u]]]
  })
  percentevent <- lapply(1:k, function(u) {
    val <- 100 * sum(cens[[u]])/n[u]
    cat(val, "% experienced event\n")
    val
  })
  longdat <- lapply(1:k, function(u) {
    data.frame(id = idl[[u]], Y = Y[[u]], time = time[[u]], study = rep(u,
                                                                        length(idl[[u]])), X1[[u]])
  })
  survdat <- lapply(1:k, function(u) {
    data.frame(id = id[[u]], survtime = survtime[[u]], cens = cens[[u]],
               study = rep(u, length(id[[u]])), X2 = X2[[u]])
  })
  list(longdat = longdat, survdat = survdat, percentevent = percentevent)
}
