#' Bootstrapping function to obtain standard errors for jointmeta1 fit
#'
#' This function takes the results of a \code{jointmeta1} fit and bootstraps
#'     it to find the standard errors of the parameter estimates.
#'
#' @param fitted a \code{jointmeta1} object
#' @param n.boot the number of bootstraps to conduct.  Note that confidence
#'     intervals will only be calculated if \code{n.boot} is greater than 100.
#' @param gpt the number of quadrature points over which the integration with
#'     respect to the random effects will be performed.  Will default to
#'     \code{gpt = 3}.
#' @param max.it the maximum number of iterations that the EM algorithm will
#'     perform for each bootstrap fit. Will default to 350.
#' @param tol the tolerance level before convergence of the algorithm is
#'     considered to have occurred. Default value is tol = 0.001.
#' @param print.detail this argument determines the level of printing that is
#'     done during the bootstrapping. If \code{TRUE} then the parameter
#'     estimates from each bootstrap sample are output.  Otherwise a progress
#'     bar is printed to indicated the proportion of bootstraps currently
#'     completed.
#' @param overalleffects this argument indicates what if any overall effects
#'     will have their standard errors and confidence intervals calculated
#'     during the bootstrap procedure. An example of an overall effect would be
#'     the combined value of a treatment effect, and a treatment by study
#'     membership interaction.  The overall treatment effect (the sum of these
#'     two values) could be of interest in an investigation.  This argument is a
#'     list containing two elements, \code{long} and \code{surv}.  Each of these
#'     elements contains a list of vectors, each of which contains the names of
#'     the parameters that make up the required overall effects
#'
#' @return a list containing three elements: \describe{
#'
#'     \item{\code{results}}{a data frame containing the estimates, standard
#'     errors and 95% confidence intervals for the parameters from the model
#'     and any overall effects requested.}
#'
#'     \item{\code{covmat}}{the covariance matrix for the model parameters}
#'
#'     \item{\code{bootstraps}}{a data frame containing the results of each
#'     bootstrap}
#' }
#'
#' @details This function takes the results of a one stage joint model fit to
#'     data from multiple studies using the function \code{jointmeta1} and
#'     performs \code{n.boot} bootstraps to determine the standard errors of the
#'     parameter estimates, and their confidence intervals if
#'     \code{n.boot > 100}.
#'
#'     The parameter \code{overalleffects} is designed for use in cases where
#'     interaction terms are included in the model specification, for example
#'     a model fitted using \code{jointmeta1} which includes both \code{treat}
#'     and \code{treat:study} where \code{treat} is a binary treatment indicator
#'     variable and \code{study} is a study indicator. In this case it may be of
#'     interest to calculate the confidence interval for the value of
#'     \code{treat + treat:study} for a given study. This is done by calculating
#'     the value of the expression for each bootstrap, and calculating the
#'     standard errors for the expression in the same way as for the other
#'     parameters.  Any overall effects to be calculated for the longitudinal
#'     sub-model are supplied as a list named \code{long} in the list
#'     \code{overalleffects}, with each element of this list containing a vector
#'     of the character names of the fixed effects to be summed to form an
#'     overall effect.  Overall effects from the survival model are specified in
#'     a similar way to an element named \code{surv}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' jointdat<-tojointdata(longitudinal = simdat2$longitudinal,
#'                       survival = simdat2$survival, id = "id",
#'                       longoutcome = "Y",
#'                       timevarying = c("time","ltime"),
#'                       survtime = "survtime", cens = "cens",
#'                       time = "time")
#'
#' onestagefit4 <- jointmeta1(data = jointdat,
#'                            long.formula = Y ~ 1 + time + treat + study,
#'                            long.rand.ind = c("int", "time"),
#'                            long.rand.stud = c("treat"),
#'                            sharingstrct = "randprop",
#'                            surv.formula = Surv(survtime, cens) ~ treat,
#'                            study.name = "study", strat = TRUE)
#'
#' onestagefit4SE <- jointmetaSE(fitted = onestagefit4, n.boot = 200)}
#'
#'
jointmetaSE<-function (fitted, n.boot, gpt, max.it, tol, print.detail = FALSE,
                       overalleffects = NULL) {
  if(inherits(fitted,"jointmeta1")) {
    if("long.rand.stud" %in% names(fitted$Call)) {
      data <- fitted$data
      id <- fitted$data$subj.col
      time.long <- fitted$data$time.col
      q <- length(diag(fitted$rand_cov$D))
      r<-length(diag(fitted$rand_cov$A))
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          overallnames.long<-unlist(
            lapply(1:length(overalleffects$long),
                   function(u) {
                     paste(unlist(overalleffects$long[u]), collapse = " and ")
                   }))
          if(is.null(overalleffects$surv) == FALSE) {
            overallnames.surv<-unlist(
              lapply(1:length(overalleffects$surv),
                     function(u) {
                       paste(unlist(overalleffects$surv[u]), collapse = " and ")
                     }))
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           overallnames.long,
                           names(fitted$coefficients$fixed$survival),
                           overallnames.surv,
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""),
                           paste("b3_", 0:(r - 1), sep = ""),
                           "Residual")
          }else{
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           overallnames.long,
                           names(fitted$coefficients$fixed$survival),
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""),
                           paste("b3_", 0:(r - 1), sep = ""),
                           "Residual")
          }
        }else{
          if(is.null(overalleffects$surv) == FALSE) {
            overallnames.surv<-unlist(
              lapply(1:length(overalleffects$surv),
                     function(u) {
                       paste(unlist(overalleffects$surv[u]),
                             collapse = " and ")
                       }))
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           names(fitted$coefficients$fixed$survival),
                           overallnames.surv,
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""),
                           paste("b3_", 0:(r - 1), sep = ""),
                           "Residual")
          }
        }
      }else{
          paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                       names(fitted$coefficients$fixed$survival),
                       names(fitted$coefficients$latent),
                       paste("b2_", 0:(q - 1), sep = ""),
                       paste("b3_", 0:(r - 1), sep = ""),
                       "Residual")
      }
      compnames <- rep("", length(paranames))
      lb1 <- length(fitted$coefficients$fixed$longitudinal[, 1])
      lb2 <- length(fitted$coefficients$fixed$survival)
      lg <- length(fitted$coefficients$latent)
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          if(is.null(overalleffects$surv) == FALSE) {
            lb1.1<-length(overalleffects$long)
            lb2.1<-length(overalleffects$surv)
            compnames[1] <- "Longitudinal"
            compnames[lb1 + 1]<-"Longitudinal Overall"
            compnames[lb1 + lb1.1 + 1] <- "Survival"
            compnames[lb1 + lb1.1 + lb2 + 1] <- "Survival Overall"
            compnames[lb1 + lb1.1 + lb2 + lb2.1 + 1] <- "Association"
            compnames[lb1 + lb1.1 + lb2 + lb2.1 + lg + 1] <- "Variance"
          }else{
            lb1.1<-length(overalleffects$long)
            compnames[1] <- "Longitudinal"
            compnames[lb1 + 1]<-"Longitudinal Overall"
            compnames[lb1 + lb1.1 + 1] <- "Survival"
            compnames[lb1 + lb1.1 + lb2 + 1] <- "Association"
            compnames[lb1 + lb1.1 + lb2 + lg + 1] <- "Variance"
          }
        }else{
          lb2.1<-length(overalleffects$surv)
          compnames[1] <- "Longitudinal"
          compnames[lb1 + 1] <- "Survival"
          compnames[lb1 + lb2 + 1] <- "Survival Overall"
          compnames[lb1 + lb2 + lb2.1 + 1] <- "Association"
          compnames[lb1 + lb2 + lb2.1 + lg + 1] <- "Variance"
        }
      }else{
        compnames[1] <- "Longitudinal"
        compnames[lb1 + 1] <- "Survival"
        compnames[lb1 + lb2 + 1] <- "Association"
        compnames[lb1 + lb2 + lg + 1] <- "Variance"
      }
      if (missing(gpt)) {
        gpt <- 3
      }
      if (missing(max.it)) {
        max.it <- 500
      }
      if (missing(tol)) {
        tol <- 0.001
      }
      long.formula<-fitted$Call$long.formula
      long.rand.ind<-as.character(fitted$Call$long.rand.ind[-1])
      long.rand.stud<-as.character(fitted$Call$long.rand.stud[-1])
      sharingstrct<-fitted$Call$sharingstrct
      if(sharingstrct !=  "randprop") {
        stop("Currently only the randprop sharing structure is supported")
      }
      surv.formula<-fitted$Call$surv.formula
      study.name<-fitted$Call$study.name
      if(as.character(fitted$Call$strat) == "T"||
         as.character(fitted$Call$strat) == "TRUE") {
        strat<-T
      }else{strat<-F}
      if("longsep" %in% names(fitted$Call)) {
        longsep<-as.logical(as.character(fitted$Call$longsep))
      }else{
        longsep<-F
      }
      if("survsep" %in% names(fitted$Call)) {
        survsep<-as.logical(as.character(fitted$Call$survsep))
      }else{
        survsep<-F
      }
      data.surv <- cbind(fitted$data$survival, fitted$data$baseline)
      surv.frame <- model.frame(surv.formula, data = data.surv)
      if (dim(surv.frame)[2] == 1) {
        n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
          dim(as.matrix(fitted$coefficients$latent))[1] +
          dim(as.matrix(diag(fitted$rand_cov$D)))[1] +
          dim(as.matrix(diag(fitted$rand_cov$A)))[1] + 1
      }else {
        n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
          dim(as.matrix(fitted$coefficients$fixed$survival))[1] +
          dim(as.matrix(fitted$coefficients$latent))[1] +
          dim(as.matrix(diag(fitted$rand_cov$D)))[1] +
          dim(as.matrix(diag(fitted$rand_cov$A)))[1] + 1
      }
      if(is.null(overalleffects$long) == FALSE) {
        n.est<-n.est + length(overalleffects$long)
      }
      if(is.null(overalleffects$surv) == FALSE) {
        n.est<-n.est + length(overalleffects$surv)
      }
      out <- matrix(0, n.boot + 2, n.est)
      nsubj <- as.numeric(table(data$baseline[,
                                              which(names(data$baseline) %in%
                                                      study.name)]))
      names(nsubj)<-studies<-names(table(
        data$baseline[, which(names(data$baseline) %in% study.name)]))
      numstudies<-length(nsubj)
      ids.surv<-data$survival[, 1]
      ids.bystudy<-lapply(1:numstudies, function(u) {
        ids.surv[which(
          data$baseline[, which(names(data$baseline) %in% study.name)] ==
            studies[u])]
        })
      id.long<-data$longitudinal[, 1]
      id.long.bystudy<-lapply(1:numstudies, function(u) {
        id.long[which(id.long %in% ids.bystudy[[u]])]
        })
      data.bystudy<-lapply(1:numstudies, function(u) {
        idstemp<-data$subject[which(
          data$baseline[, which(names(data$baseline) == study.name)] ==
            names(nsubj)[[u]])]
        subset(data, idstemp)
      })
      for(i in 1:numstudies) {class(data.bystudy[[i]])<-"jointdata"}
      nn<-lapply(1:numstudies, function(u) {
        nn <- diff(match(unique(id.long.bystudy[[u]]), id.long.bystudy[[u]]))
        nn <- c(nn, length(id.long.bystudy[[u]]) - sum(nn))
      })
      if(print.detail == FALSE) {
        pb <- txtProgressBar(min = 0, max = n.boot, style = 3)
        counter<-1
      }
      for(i in 1:n.boot) {
        s.new<-lapply(1:numstudies, function(u) {
          dataout<-sample.jointdata(data.bystudy[[u]],
                                    size = nsubj[u], replace = TRUE)
          newids<-1:nsubj[u]
          if(u>1) {
            dataout$longitudinal[, 1]<-dataout$longitudinal[, 1] +
              sum(nsubj[1:(u-1)])
            dataout$survival[, 1]<-dataout$survival[, 1] + sum(nsubj[1:(u-1)])
            dataout$baseline[, 1]<-dataout$subject<-dataout$baseline[, 1] +
              sum(nsubj[1:(u-1)])
          }
          dataout
        })
        s.new.long<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {
                                     s.new[[u]]$longitudinal
                                     }))
        s.new.surv<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {
                                     s.new[[u]]$survival
                                     }))
        s.new.base<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {
                                     s.new[[u]]$baseline
                                     }))
        s.new.subject<-do.call(c,
                               lapply(1:numstudies,
                                      function(u) {
                                        s.new[[u]]$subject
                                        }))
        s.new2<-jointdata(longitudinal = s.new.long,
                          survival = s.new.surv,
                          baseline = s.new.base,
                          id.col = s.new[[1]]$subj.col,
                          time.col = s.new[[1]]$time.col)
        class(s.new2)<-"jointdata"
        fitb<-NULL
        tryCatch(suppressWarnings(fitb <- jointmeta1(data = s.new2, long.formula = long.formula,
                                    long.rand.ind = long.rand.ind,
                                    long.rand.stud = long.rand.stud,
                                    sharingstrct = sharingstrct,
                                    surv.formula = surv.formula, gpt = gpt,
                                    max.it = max.it, tol = tol,
                                    study.name = study.name, strat = strat,
                                    longsep = longsep, survsep = survsep,
                                    bootrun = TRUE)),
                 error = function(e) {
                   fitb<-NULL
                 })
        if(is.null(fitb) == FALSE) {
          b1 <- as.numeric(as.vector(as.matrix(
            fitb$coefficients$fixed$longitudinal[, 1])))
          b3 <- as.numeric(as.vector(as.matrix(fitb$coefficients$latent)))
          b4 <- c(as.numeric(as.vector(as.matrix(diag(fitb$rand_cov$D)))),
                  as.numeric(as.vector(as.matrix(diag(fitb$rand_cov$A)))))
          b5 <- as.numeric(as.vector(as.matrix(fitb$sigma.e)))
          if(is.null(overalleffects$long) == FALSE) {
            overall.long<-overalleffects$long
            b1.overall<-rep(NA, length(overall.long))
            for(count in 1:length(overall.long)) {
              b1.overall[count]<-sum(b1[which(
                rownames(fitb$coefficients$fixed$longitudinal) %in%
                  unlist(overall.long[count]))])
            }
            b1<-c(b1, b1.overall)
          }
          if (dim(surv.frame)[2] != 1) {
            b2 <- as.numeric(as.vector(as.matrix(
              fitb$coefficients$fixed$survival)))
            if(is.null(overalleffects$surv) == FALSE) {
              overall.surv<-overalleffects$surv
              b2.overall<-rep(NA, length(overall.surv))
              for(count in 1:length(overall.surv)) {
                b2.overall[count]<-sum(b2[which(names(
                  fitted$coefficients$fixed$survival) %in%
                    unlist(overall.surv[count]))])
              }
              b2<-c(b2, b2.overall)
            }
            out[i, ] <- c(b1, b2, b3, b4, b5)
            ests <- out[i, ]
            if (print.detail) {
              detail <- data.frame(iteration = i, t(ests))
              names(detail) <- c("Iteration", paranames)
              print(detail)
            }
          }else {
            out[i, ] <- c(b1, b3, b4, b5)
            ests <- out[i, ]
            if (print.detail) {
              detail <- data.frame(iteration = i, t(ests))
              names(detail) <- c("Iteration", paranames)
              print(detail)
            }
          }
        }
        if(print.detail == FALSE) {
          setTxtProgressBar(pb, counter)
          counter<-counter + 1
        }
      }
      if(print.detail == FALSE) {close(pb)}
      i<-which(out[, 1] == 0)
      out1 <- out[-i, ]
      se <- 0
      ci1 <- 0
      ci2 <- 0
      if (n.boot == 1) {
        out <- matrix(out, nrow = 1)
      }
      for (i in 1:length(out[1, ])) {
        se[i] <- sqrt(var(as.numeric(out[, i])))
        if (n.boot < 100) {
          ci1[i] <- 0
          ci2[i] <- 0
        }
        else {
          ci1[i] <- sort(as.numeric(out[, i]))[2.5/100 * n.boot]
          ci2[i] <- sort(as.numeric(out[, i]))[97.5/100 * n.boot]
        }
      }
      covmat<-matrix(0, nrow = ncol(out), ncol = ncol(out))
      diag(covmat)<-se^2
      means<-colMeans(out)
      for(i in 1:(ncol(out)-1)) {
        for(j in (i + 1):ncol(out)) {
          covmat[i, j]<-(1/(nrow(out)-1))*
            sum((out[, i]-means[i])*(out[, j]-means[j]))
        }
      }
      covmat<-covmat + t(covmat)-diag(diag(covmat))
      colnames(covmat)<-rownames(covmat)<-paranames
      if (dim(surv.frame)[2] != 1) {
        if(is.null(overalleffects) == FALSE) {
          if(is.null(overalleffects$long) == FALSE) {
            overall.long<-overalleffects$long
            b1.overall.out<-rep(NA, length(overall.long))
            for(count in 1:length(overall.long)) {
              b1.overall.out[count]<-sum(as.numeric(
                as.vector(as.matrix(
                  fitted$coefficients$fixed$longitudinal)))[
                    which(rownames(
                      fitted$coefficients$fixed$longitudinal) %in%
                        unlist(overall.long[count]))])
            }
            if(is.null(overalleffects$surv) == FALSE) {
              overall.surv<-overalleffects$surv
              b2.overall.out<-rep(NA, length(overall.surv))
              for(count in 1:length(overall.surv)) {
                b2.overall.out[count]<-sum(as.numeric(
                  as.vector(as.matrix(
                    fitted$coefficients$fixed$survival)))[
                      which(names(fitted$coefficients$fixed$survival) %in%
                              unlist(overall.surv[count]))])
              }
              b1 <- data.frame(cbind(compnames,
                                     paranames,
                                     round(c(
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$longitudinal))),
                                       b1.overall.out,
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$survival))),
                                       b2.overall.out,
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$latent))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$D)))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$A)))),
                                       as.numeric(
                                         as.vector(as.matrix(fitted$sigma.e)))),
                                       4),
                                     round(cbind(se), 4),
                                     round(ci1, 4), round(ci2, 4)))
            }else{
              b1 <- data.frame(cbind(compnames,
                                     paranames,
                                     round(c(
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$longitudinal))),
                                       b1.overall.out,
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$survival))),
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$latent))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$D)))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$A)))),
                                       as.numeric(
                                         as.vector(as.matrix(fitted$sigma.e)))),
                                       4), round(cbind(se), 4),
                                     round(ci1, 4), round(ci2, 4)))
            }
          }else{
            if(is.null(overalleffects$surv) == FALSE) {
              b1 <- data.frame(cbind(compnames,
                                     paranames,
                                     round(c(
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$longitudinal))),
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$
                                             fixed$survival))),
                                       b2.overall.out,
                                       as.numeric(
                                         as.vector(as.matrix(
                                           fitted$coefficients$latent))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$D)))),
                                       as.numeric(
                                         as.vector(
                                           as.matrix(diag(fitted$rand_cov$A)))),
                                       as.numeric(
                                         as.vector(as.matrix(fitted$sigma.e)))),
                                       4), round(cbind(se), 4),
                                     round(ci1, 4), round(ci2, 4)))
            }
          }
        }else{
          b1 <- data.frame(cbind(compnames,
                                 paranames,
                                 round(c(
                                   as.numeric(
                                     as.vector(as.matrix(
                                       fitted$coefficients$
                                         fixed$longitudinal))),
                                   as.numeric(
                                     as.vector(as.matrix(
                                       fitted$coefficients$fixed$survival))),
                                   as.numeric(
                                     as.vector(as.matrix(
                                       fitted$coefficients$latent))),
                                   as.numeric(
                                     as.vector(
                                       as.matrix(diag(fitted$rand_cov$D)))),
                                   as.numeric(
                                     as.vector(
                                       as.matrix(diag(fitted$rand_cov$A)))),
                                   as.numeric(
                                     as.vector(as.matrix(fitted$sigma.e)))),4),
                                 round(cbind(se), 4),
                                 round(ci1, 4), round(ci2, 4)))
        }
      } else {
        if(is.null(overalleffects) == FALSE) {
          if(is.null(overalleffects$long) == FALSE) {
            overall.long<-overalleffects$long
            b1.overall.out<-rep(NA, length(overall.long))
            for(count in 1:length(overall.long)) {
              b1.overall.out[count]<-sum(
                as.numeric(
                  as.vector(
                    as.matrix(
                      fitted$coefficients$fixed$longitudinal)))[
                        which(rownames(
                          fitted$coefficients$fixed$longitudinal) %in%
                            unlist(overall.long[count]))])
            }

            b1 <- data.frame(cbind(compnames,
                                   paranames,
                                   round(c(
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$
                                           fixed$longitudinal))),
                                     b1.overall.out,
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$latent))),
                                     as.numeric(
                                       as.vector(
                                         as.matrix(diag(fitted$rand_cov$D)))),
                                     as.numeric(
                                       as.vector(
                                         as.matrix(diag(fitted$rand_cov$A)))),
                                     as.numeric(
                                       as.vector(as.matrix(fitted$sigma.e)))),
                                     4), round(cbind(se), 4),
                                   round(ci1, 4), round(ci2,4)))
          }
        }else{
          b1 <- data.frame(cbind(compnames,
                                 paranames,
                                 round(c(
                                   as.numeric(
                                     as.vector(as.matrix(
                                       fitted$coefficients$
                                         fixed$longitudinal))),
                                   as.numeric(
                                     as.vector(as.matrix(
                                       fitted$coefficients$latent))),
                                   as.numeric(
                                     as.vector(
                                       as.matrix(diag(fitted$rand_cov$D)))),
                                   as.numeric(
                                     as.vector(
                                       as.matrix(diag(fitted$rand_cov$A)))),
                                   as.numeric(
                                     as.vector(as.matrix(fitted$sigma.e)))),
                                   4), round(cbind(se), 4),
                                 round(ci1, 4), round(ci2, 4)))
        }
      }
      names(b1)[1:6] <- c("Component", "Parameter", "Estimate",
                          "SE", "95%Lower", "95%Upper")
      output<-list(results = b1, covmat = covmat, bootstraps = out)
      class(output)<-"jointmeta1SE"
      return(output)
    }else{
      data <- fitted$data
      id <- fitted$data$subj.col
      time.long <- fitted$data$time.col
      q <- length(diag(fitted$rand_cov$D))
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          overallnames.long<-unlist(
            lapply(1:length(overalleffects$long),
                   function(u) {
                     paste(unlist(overalleffects$long[u]), collapse = " and ")
                     }))
          if(is.null(overalleffects$surv) == FALSE) {
            overallnames.surv<-unlist(
              lapply(1:length(overalleffects$surv),
                     function(u) {
                       paste(unlist(overalleffects$surv[u]),
                             collapse = " and ")
                       }))
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           overallnames.long,
                           names(fitted$coefficients$fixed$survival),
                           overallnames.surv,
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""), "Residual")
          }else{
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           overallnames.long,
                           names(fitted$coefficients$fixed$survival),
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""), "Residual")
          }
        }else{
          if(is.null(overalleffects$surv) == FALSE) {
            overallnames.surv<-unlist(
              lapply(1:length(overalleffects$surv),
                     function(u) {
                       paste(unlist(overalleffects$surv[u]),
                             collapse = " and ")
                       }))
            paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                           names(fitted$coefficients$fixed$survival),
                           overallnames.surv,
                           names(fitted$coefficients$latent),
                           paste("b2_", 0:(q - 1), sep = ""), "Residual")
          }
        }
      }else{
        paranames <- c(row.names(fitted$coefficients$fixed$longitudinal),
                       names(fitted$coefficients$fixed$survival),
                       names(fitted$coefficients$latent),
                       paste("b2_", 0:(q - 1), sep = ""), "Residual")
      }
      compnames <- rep("", length(paranames))
      lb1 <- length(fitted$coefficients$fixed$longitudinal[, 1])
      lb2 <- length(fitted$coefficients$fixed$survival)
      lg <- length(fitted$coefficients$latent)
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          if(is.null(overalleffects$surv) == FALSE) {
            lb1.1<-length(overalleffects$long)
            lb2.1<-length(overalleffects$surv)
            compnames[1] <- "Longitudinal"
            compnames[lb1 + 1]<-"Longitudinal Overall"
            compnames[lb1 + lb1.1 + 1] <- "Survival"
            compnames[lb1 + lb1.1 + lb2 + 1] <- "Survival Overall"
            compnames[lb1 + lb1.1 + lb2 + lb2.1 + 1] <- "Association"
            compnames[lb1 + lb1.1 + lb2 + lb2.1 + lg + 1] <- "Variance"
          }else{
            lb1.1<-length(overalleffects$long)
            compnames[1] <- "Longitudinal"
            compnames[lb1 + 1]<-"Longitudinal Overall"
            compnames[lb1 + lb1.1 + 1] <- "Survival"
            compnames[lb1 + lb1.1 + lb2 + 1] <- "Association"
            compnames[lb1 + lb1.1 + lb2 + lg + 1] <- "Variance"
          }
        }else{
          lb2.1<-length(overalleffects$surv)
          compnames[1] <- "Longitudinal"
          compnames[lb1 + 1] <- "Survival"
          compnames[lb1 + lb2 + 1] <- "Survival Overall"
          compnames[lb1 + lb2 + lb2.1 + 1] <- "Association"
          compnames[lb1 + lb2 + lb2.1 + lg + 1] <- "Variance"
        }
      }else{
        compnames[1] <- "Longitudinal"
        compnames[lb1 + 1] <- "Survival"
        compnames[lb1 + lb2 + 1] <- "Association"
        compnames[lb1 + lb2 + lg + 1] <- "Variance"
      }
      if (missing(gpt)) {
        gpt <- 3
      }
      if (missing(max.it)) {
        max.it <- 350
      }
      if (missing(tol)) {
        tol <- 0.001
      }
      long.formula<-fitted$Call$long.formula
      long.rand.ind<-as.character(fitted$Call$long.rand.ind[-1])
      sharingstrct<-fitted$Call$sharingstrct
      surv.formula<-fitted$Call$surv.formula
      study.name<-fitted$Call$study.name
      if(sharingstrct != "randprop") {
        stop("Currently only the randprop sharing structure is supported")
      }
      if(as.character(fitted$Call$strat) == "T"||
         as.character(fitted$Call$strat) == "TRUE") {
        strat<-T
      }else{strat<-F}
      if("longsep" %in% names(fitted$Call)) {
        longsep<-as.logical(as.character(fitted$Call$longsep))
      }else{
        longsep<-F
      }
      if("survsep" %in% names(fitted$Call)) {
        survsep<-as.logical(as.character(fitted$Call$survsep))
      }else{
        survsep<-F
      }
      data.surv <- cbind(fitted$data$survival, fitted$data$baseline)
      surv.frame <- model.frame(surv.formula, data = data.surv)
      if (dim(surv.frame)[2] == 1) {
        n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
          dim(as.matrix(fitted$coefficients$latent))[1] +
          dim(as.matrix(diag(fitted$rand_cov$D)))[1] + 1
      }else {
        n.est <- dim(as.matrix(fitted$coefficients$fixed$longitudinal))[1] +
          dim(as.matrix(fitted$coefficients$fixed$survival))[1] +
          dim(as.matrix(fitted$coefficients$latent))[1] +
          dim(as.matrix(diag(fitted$rand_cov$D)))[1] + 1
      }
      if(is.null(overalleffects$long) == FALSE) {
        n.est<-n.est + length(overalleffects$long)
      }
      if(is.null(overalleffects$surv) == FALSE) {
        n.est<-n.est + length(overalleffects$surv)
      }
      out <- matrix(0, n.boot + 2, n.est)
      nsubj <- as.numeric(
        table(data$baseline[, which(names(data$baseline) %in% study.name)]))
      names(nsubj)<-studies<-names(
        table(data$baseline[, which(names(data$baseline) %in% study.name)]))
      numstudies<-length(nsubj)
      ids.surv<-data$survival[, 1]
      ids.bystudy<-lapply(1:numstudies,
                          function(u) {
                            ids.surv[
                              which(
                                data$baseline[, which(names(data$baseline) %in%
                                                        study.name)] ==
                                  studies[u])]
                            })
      id.long<-data$longitudinal[, 1]
      id.long.bystudy<-lapply(1:numstudies,
                              function(u) {
                                id.long[which(id.long %in% ids.bystudy[[u]])]
                              })
      data.bystudy<-lapply(1:numstudies, function(u) {
        idstemp<-data$subject[
          which(data$baseline[, which(names(data$baseline) == study.name)] ==
                  names(nsubj)[[u]])]
        subset(data, idstemp)
      })
      for(i in 1:numstudies) {class(data.bystudy[[i]])<-"jointdata"}
      nn<-lapply(1:numstudies, function(u) {
        nn <- diff(match(unique(id.long.bystudy[[u]]), id.long.bystudy[[u]]))
        nn <- c(nn, length(id.long.bystudy[[u]]) - sum(nn))
      })
      if(print.detail == FALSE) {
        pb <- txtProgressBar(min = 0, max = n.boot, style = 3)
        counter<-1
      }
      for(i in 1:n.boot) {
        s.new<-lapply(1:numstudies, function(u) {
          dataout<-sample.jointdata(data.bystudy[[u]],
                                    size = nsubj[u], replace = T)
          newids<-1:nsubj[u]
          if(u>1) {
            dataout$longitudinal[, 1]<-dataout$longitudinal[, 1] +
              sum(nsubj[1:(u-1)])
            dataout$survival[, 1]<-dataout$survival[, 1] + sum(nsubj[1:(u-1)])
            dataout$baseline[, 1]<-dataout$subject<-dataout$baseline[, 1] +
              sum(nsubj[1:(u-1)])
          }
          dataout
        })
        s.new.long<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {s.new[[u]]$longitudinal}))
        s.new.surv<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {s.new[[u]]$survival}))
        s.new.base<-do.call(rbind,
                            lapply(1:numstudies,
                                   function(u) {s.new[[u]]$baseline}))
        s.new.subject<-do.call(c,
                               lapply(1:numstudies,
                                      function(u) {s.new[[u]]$subject}))
        s.new2<-jointdata(longitudinal = s.new.long,
                          survival = s.new.surv, baseline = s.new.base,
                          id.col = s.new[[1]]$subj.col,
                          time.col = s.new[[1]]$time.col)
        class(s.new2)<-"jointdata"
        fitb<-NULL
        tryCatch(suppressWarnings(fitb <- jointmeta1(data = s.new2, long.formula = long.formula,
                                    long.rand.ind = long.rand.ind,
                                    sharingstrct = sharingstrct,
                                    surv.formula = surv.formula, gpt = gpt,
                                    max.it = max.it, tol = tol,
                                    study.name = study.name, strat = strat,
                                    longsep = longsep, survsep = survsep,
                                    print.detail = FALSE,
                                    bootrun = TRUE)),
                 error = function(e) {
                   fitb<-NULL
                                    })
        if(is.null(fitb) == FALSE) {
          b1 <- as.numeric(as.vector(
            as.matrix(fitb$coefficients$fixed$longitudinal[, 1])))
          b3 <- as.numeric(as.vector(as.matrix(fitb$coefficients$latent)))
          b4 <- as.numeric(as.vector(as.matrix(diag(fitb$rand_cov$D))))
          b5 <- as.numeric(as.vector(as.matrix(fitb$sigma.e)))
          if(is.null(overalleffects$long) == FALSE) {
            overall.long<-overalleffects$long
            b1.overall<-rep(NA, length(overall.long))
            for(count in 1:length(overall.long)) {
              b1.overall[count]<-sum(
                b1[which(rownames(
                  fitb$coefficients$fixed$longitudinal) %in%
                    unlist(overall.long[count]))])
            }
            b1<-c(b1, b1.overall)
          }
          if (dim(surv.frame)[2] != 1) {
            b2 <- as.numeric(
              as.vector(as.matrix(fitb$coefficients$fixed$survival)))
            if(is.null(overalleffects$surv) == FALSE) {
              overall.surv<-overalleffects$surv
              b2.overall<-rep(NA, length(overall.surv))
              for(count in 1:length(overall.surv)) {
                b2.overall[count]<-sum(
                  b2[which(names(
                    fitted$coefficients$fixed$survival) %in%
                      unlist(overall.surv[count]))])
              }
              b2<-c(b2, b2.overall)
            }
            out[i, ] <- c(b1, b2, b3, b4, b5)
            ests <- out[i, ]
            if (print.detail) {
              detail <- data.frame(iteration = i, t(ests))
              names(detail) <- c("Iteration", paranames)
              print(detail)
            }
          }else {
            out[i, ] <- c(b1, b3, b4, b5)
            ests <- out[i, ]
            if (print.detail) {
              detail <- data.frame(iteration = i, t(ests))
              names(detail) <- c("Iteration", paranames)
              print(detail)
            }
          }
        }
        if(print.detail == FALSE) {
          setTxtProgressBar(pb, counter)
          counter<-counter + 1
        }
      }
      if(print.detail == FALSE) {close(pb)}
      i <- 1
      while (out[i, 1] != 0) i = i + 1
      out <- out[1:(i - 1), ]
      se <- 0
      ci1 <- 0
      ci2 <- 0
      if (n.boot == 1) {
        out <- matrix(out, nrow = 1)
      }
      for (i in 1:length(out[1, ])) {
        se[i] <- sqrt(var(as.numeric(out[, i])))
        if (n.boot < 100) {
          ci1[i] <- 0
          ci2[i] <- 0
        }
        else {
          ci1[i] <- sort(as.numeric(out[, i]))[2.5/100 * n.boot]
          ci2[i] <- sort(as.numeric(out[, i]))[97.5/100 * n.boot]
        }
      }
    }
    covmat<-matrix(0, nrow = ncol(out), ncol = ncol(out))
    diag(covmat)<-se^2
    means<-colMeans(out)
    for(i in 1:(ncol(out)-1)) {
      for(j in (i + 1):ncol(out)) {
        covmat[i, j]<-(1/(nrow(out)-1))*
          sum((out[, i]-means[i])*(out[, j]-means[j]))
      }
    }
    covmat<-covmat + t(covmat)-diag(diag(covmat))
    colnames(covmat)<-rownames(covmat)<-paranames
    if (dim(surv.frame)[2] != 1) {
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          overall.long<-overalleffects$long
          b1.overall.out<-rep(NA, length(overall.long))
          for(count in 1:length(overall.long)) {
            b1.overall.out[count]<-sum(
              as.numeric(as.vector(as.matrix(
                fitted$coefficients$fixed$longitudinal)))[
                  which(rownames(fitted$coefficients$fixed$longitudinal) %in%
                          unlist(overall.long[count]))])
          }
          if(is.null(overalleffects$surv) == FALSE) {
            overall.surv<-overalleffects$surv
            b2.overall.out<-rep(NA, length(overall.surv))
            for(count in 1:length(overall.surv)) {
              b2.overall.out[count]<-sum(
                as.numeric(as.vector(as.matrix(
                  fitted$coefficients$fixed$survival)))[
                    which(names(fitted$coefficients$fixed$survival) %in%
                            unlist(overall.surv[count]))])
            }
            b1 <- data.frame(cbind(compnames,
                                   paranames,
                                   round(c(
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$
                                           fixed$longitudinal))),
                                     b1.overall.out,
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$fixed$survival))),
                                     b2.overall.out,
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$latent))),
                                     as.numeric(
                                       as.vector(as.matrix(diag(fitted$rand_cov$D)))),
                                     as.numeric(
                                       as.vector(as.matrix(fitted$sigma.e)))),
                                     4), round(cbind(se), 4),
                                   round(ci1, 4), round(ci2, 4)))
          }else{
            b1 <- data.frame(cbind(compnames,
                                   paranames,
                                   round(c(
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$
                                           fixed$longitudinal))),
                                     b1.overall.out,
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$fixed$survival))),
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$latent))),
                                     as.numeric(
                                       as.vector(as.matrix(diag(fitted$rand_cov$D)))),
                                     as.numeric(
                                       as.vector(as.matrix(fitted$sigma.e)))),
                                     4), round(cbind(se), 4),
                                   round(ci1, 4), round(ci2, 4)))
          }
        }else{
          if(is.null(overalleffects$surv) == FALSE) {
            b1 <- data.frame(cbind(compnames,
                                   paranames,
                                   round(c(
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$
                                           fixed$longitudinal))),
                                     as.numeric(
                                       as.vector(as.matrix(
                                         fitted$coefficients$fixed$survival))),
                                     b2.overall.out,
                                     as.numeric(as.vector(
                                       as.matrix(fitted$coefficients$latent))),
                                     as.numeric(as.vector(
                                       as.matrix(diag(fitted$rand_cov$D)))),
                                     as.numeric(as.vector(
                                       as.matrix(fitted$sigma.e)))), 4),
                                   round(cbind(se), 4),
                                   round(ci1, 4), round(ci2, 4)))
          }
        }
      }else{
        b1 <- data.frame(cbind(compnames,
                               paranames,
                               round(c(
                                 as.numeric(as.vector(
                                   as.matrix(
                                     fitted$coefficients$fixed$longitudinal))),
                                 as.numeric(as.vector(
                                   as.matrix(
                                     fitted$coefficients$fixed$survival))),
                                 as.numeric(as.vector(
                                   as.matrix(fitted$coefficients$latent))),
                                 as.numeric(as.vector(
                                   as.matrix(diag(fitted$rand_cov$D)))),
                                 as.numeric(as.vector(
                                   as.matrix(fitted$sigma.e)))), 4),
                               round(cbind(se), 4),
                               round(ci1, 4), round(ci2, 4)))
      }
    } else {
      if(is.null(overalleffects) == FALSE) {
        if(is.null(overalleffects$long) == FALSE) {
          overall.long<-overalleffects$long
          b1.overall.out<-rep(NA, length(overall.long))
          for(count in 1:length(overall.long)) {
            b1.overall.out[count]<-sum(
              as.numeric(
                as.vector(
                  as.matrix(
                    fitted$coefficients$fixed$longitudinal)))[
                      which(rownames(fitted$coefficients$fixed$longitudinal) %in%
                              unlist(overall.long[count]))])
          }
          b1 <- data.frame(cbind(compnames,
                                 paranames,
                                 round(c(
                                   as.numeric(as.vector(
                                     as.matrix(
                                       fitted$coefficients$
                                         fixed$longitudinal))),
                                   b1.overall.out,
                                   as.numeric(as.vector(
                                     as.matrix(fitted$coefficients$latent))),
                                   as.numeric(as.vector(
                                     as.matrix(diag(fitted$rand_cov$D)))),
                                   as.numeric(as.vector(
                                     as.matrix(fitted$sigma.e)))), 4),
                                 round(cbind(se), 4),
                                 round(ci1, 4), round(ci2, 4)))
        }
      }else{
        b1 <- data.frame(cbind(compnames,
                               paranames,
                               round(c(
                                 as.numeric(as.vector(
                                   as.matrix(
                                     fitted$coefficients$fixed$longitudinal))),
                                 as.numeric(as.vector(
                                   as.matrix(fitted$coefficients$latent))),
                                 as.numeric(as.vector(
                                   as.matrix(diag(fitted$rand_cov$D)))),
                                 as.numeric(as.vector(
                                   as.matrix(fitted$sigma.e)))),4),
                               round(cbind(se), 4),
                               round(ci1, 4),
                               round(ci2,4)))
      }
    }
    names(b1)[1:6] <- c("Component", "Parameter", "Estimate",
                        "SE", "95%Lower", "95%Upper")
    output<-list(results = b1, covmat = covmat, bootstraps = out)
    class(output)<-"jointmeta1SE"
    return(output)
  }else{
    stop("Supplied fitted should be of class jointmeta1")
  }
}
