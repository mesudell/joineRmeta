#' Function to pool joint model fits in two stage MA
#'
#' This function takes joint model fits from either \code{\link[joineR]{joint}}
#' or \code{\link[JM]{jointModel}} and pools the information from the fits in
#' the second stage of a two stage meta-analysis (MA).
#'
#' @param fits a list of joint modelling fits.  These fits should all be of the
#'   same type (i.e. all fitted using \code{\link[joineR]{joint}} or all fitted
#'   using \code{\link[JM]{jointModel}}), with the same model specification.
#' @param SE a list of the results from \code{\link[joineR]{jointSE}}.  Only to
#'   be supplied if the model fits supplied in \code{fits} are all fitted using
#'   the \code{\link[joineR]{joineR}} package.
#' @param longpar a vector of character strings of parameters from the
#'   longitudinal sub-model for which meta-analyses should be performed
#' @param survpar a vector of character strings of parameters from the survival
#'   sub-model for which meta-analyses should be performed
#' @param assoc a logical indicating whether a meta-analysis should be performed
#'   for the association parameter(s)
#' @param studynames a vector of character strings containing the names for the
#'   studies present in the dataset that the joint models were fitted to. These
#'   character strings if supplied are used to label the meta-analyses performed
#'   by the function
#'
#' @return This function returns a list of results for the two stage MA.  These
#'   results are split by the type of parameter being pooled.  If the names of
#'   longitudinal parameters were supplied to \code{longpar} then an element
#'   named \code{longMA} will be present in the results.  If the names of
#'   survival parameters were supplied to \code{survpar} then if the supplied
#'   joint model fits were fitted using the \code{joint} function from the
#'   \code{joineR} package, an element named \code{survMA.direct} will be
#'   present in the results.  If the supplied joint model fits were fitted using
#'   the \code{jointModel} function from the \code{JM} package, two elements
#'   named \code{survMA.direct} and \code{survMA.overall} will be present.  If
#'   \code{assoc = TRUE} then an element labelled \code{assocMA} will be present
#'   in the results.
#'
#'   Each element of each of these components of the results (\code{longMA},
#'   \code{survMA.direct}, \code{assocMA}...) is of class \code{metagen}, and is
#'   the result of using the \code{\link[meta]{metagen}} function on the results
#'   of joint models fitted to multiple studies in the dataset. This method
#'   pools the supplied information in fixed and random MA using inverse
#'   variance weighting.  Forest plots can be produced for these results simply
#'   by applying the function \code{\link[meta]{forest}} to the objects of class
#'   \code{metagen} / \code{meta} supplied in the results.
#'
#' @details This function is designed to take study specific joint model fits,
#'   and pool them in a meta-analysis.  This is the second stage in a two stage
#'   IPD meta-analysis.  The function can handle fits using the
#'   \code{\link[joineR]{joint}} function or the \code{\link[JM]{jointModel}},
#'   but not a mixture of the two.  This is due to the differing types of
#'   association structures available in each package (for more details see
#'   \code{\link[joineR]{joineR}}, and \code{\link[JM]{JM}}).
#'
#'   The fixed effects in the survival sub-models estimated in both packages
#'   have slightly different interpretations.  The joint model fits modelled
#'   using the \code{joineR} package link the sub-models using shared zero mean
#'   random effects (see Henderson et al (2000)).  However the joint model fits
#'   modelled using the \code{JM} package link the sub-models using sharing
#'   structures that involve both the fixed and random effects.  If a parameter
#'   specified in survpar is also present in the fixed effects of the
#'   longitudinal sub-model, a direct effect of the parameter on the risk of an
#'   event can be extracted from the survival sub-model, as well as the overall
#'   effect resulting from the sum of fixed effect in the survival sub-model,
#'   and the presence of the parameter in the longitudinal sub-model, present in
#'   the sharing structure of the joint model. As such, if a parameter specified
#'   in \code{survpar} is also present as a fixed effect in the longitudinal
#'   sub-model, and the fixed and random effects make up the sharing structure
#'   linking the sub-models, the overall parameter effect is found by
#'   \eqn{\beta_2 + (\alpha * \beta_1)}, where \eqn{\alpha} is the association
#'   parameter, \eqn{\beta_2} is the coefficient for the parameter in question
#'   from the survival sub-model, and \eqn{\beta_1} is the coefficient for the
#'   parameter in question from the longitudinal sub-model.  For more
#'   information about overall effects versus direct effects see Ibrahim et al
#'   (2010), Rizopoulos (2012) and Gould et al (2015).  Because both a direct
#'   and an overall effect of the survival parameters can be extracted from the
#'   model, both are present in the results if the joint models supplied in the
#'   fits are fitted using the \code{JM} package.
#'
#' @export
#' @import JM joineR meta msm
#'
#' @seealso \code{\link[joineR]{joint}}, \code{\link[JM]{jointModel}},
#'   \code{\link[joineR]{jointSE}}, \code{\link[meta]{metagen}}
#'
#' @references Ibrahim et al (2010) Basic Concepts and Methods for Joint Models
#'   of Longitudinal and Survival Data. JOURNAL OF CLINICAL ONCOLOGY 28 (10):
#'   2796-2801
#'
#'   Rizopoulos (2012) Joint Models for Longitudinal and Time-to-Event Data With
#'   Applications in R. Chapman and Hall/CRC Biostatistics Series
#'
#'   Henderson et al (2000) Joint modelling of longitudinal measurements and
#'   event time data. Biostatistics, 1,4, pp. 465–480
#'
#'   Gould et al (2015) Joint modeling of survival and longitudinal non-survival
#'   data: current methods and issues. Report of the DIA Bayesian joint modeling
#'   working group.  Statistics in Medicine 34(14): 2181–2195.
#'   doi:10.1002/sim.6141.
#'
#' @examples
#'
#'     #change data to jointdata format
#'     jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                           survival = simdat$survival, id = 'id',
#'                           longoutcome = 'Y',
#'                           timevarying = c('time','ltime'),
#'                           survtime = 'survtime', cens = 'cens',
#'                           time = 'time')
#'
#'     #ensure variables are correctly formatted
#'     jointdat$baseline$study <- as.factor(jointdat$baseline$study)
#'     jointdat$baseline$treat <- as.factor(jointdat$baseline$treat)
#'
#'     #subset the data by study
#'     for(i in 1:length(unique(jointdat$baseline$study))){
#'          idstemp<-jointdat$baseline$id[which(jointdat$baseline$study %in%
#'                                          unique(jointdat$baseline$study)[i])]
#'          temp<-subset(jointdat,idstemp)
#'          class(temp)<-'jointdata'
#'          assign(paste('jointdat',unique(jointdat$baseline$study)[i],
#'                       sep='.'),temp)
#'       }
#'
#'     #############################
#'     ### Example using joineR fits
#'
#'     #use the joineR package to fit study specific joint models
#'     joineRfit1<-joint(data = jointdat.1, long.formula = Y ~ 1 + time + treat,
#'                   surv.formula = Surv(survtime, cens) ~ treat,
#'                   model = 'intslope')
#'
#'     joineRfit2<-joint(data = jointdat.2, long.formula = Y ~ 1 + time + treat,
#'                   surv.formula = Surv(survtime, cens) ~ treat,
#'                   model = 'intslope')
#'
#'     joineRfit3<-joint(data = jointdat.3, long.formula = Y ~ 1 + time + treat,
#'                   surv.formula = Surv(survtime, cens) ~ treat,
#'                   model = 'intslope')
#'
#'     joineRfit4<-joint(data = jointdat.4, long.formula = Y ~ 1 + time + treat,
#'                   surv.formula = Surv(survtime, cens) ~ treat,
#'                   model = 'intslope')
#'
#'     joineRfit5<-joint(data = jointdat.5, long.formula = Y ~ 1 + time + treat,
#'                   surv.formula = Surv(survtime, cens) ~ treat,
#'                   model = 'intslope')
#'
#'     joineRfit1SE<-jointSE(fitted = joineRfit1, n.boot = 200)
#'     joineRfit2SE<-jointSE(fitted = joineRfit2, n.boot = 200)
#'     joineRfit3SE<-jointSE(fitted = joineRfit3, n.boot = 200)
#'     joineRfit4SE<-jointSE(fitted = joineRfit4, n.boot = 200)
#'     joineRfit5SE<-jointSE(fitted = joineRfit5, n.boot = 200)
#'
#'     joineRfits<-list(joineRfit1, joineRfit1SE,
#'     joineRfit2, joineRfit2SE,
#'     joineRfit3, joineRfit3SE,
#'     joineRfit4, joineRfit4SE,
#'     joineRfit5, joineRfit5SE)
#'
#'     names(joineRfits)<-c('joineRfit1', 'joineRfit1SE',
#'                          'joineRfit2', 'joineRfit2SE',
#'                          'joineRfit3', 'joineRfit3SE',
#'                          'joineRfit4', 'joineRfit4SE',
#'                          'joineRfit5', 'joineRfit5SE')
#'
#'     #perform the second stage of the two stage MA
#'     MAjoineRfits<-jointmeta2(fits = joineRmodels, SE = joineRmodelsSE,
#'                   longpar = c('time', 'treat1'), survpar = 'treat1',
#'                   assoc = TRUE, studynames = c('Study 1','Study 2',
#'                   'Study 3', 'Study 4', 'Study 5'))
#'
#'     #produce forest plots
#'     library(meta)
#'     forest(MAjoineRfits$longMA$treat1)
#'
#'     ###########################
#'     ### Example using JM fits
#'
#'     longdat.1<-merge(jointdat.1$longitudinal,jointdat.1$baseline,by='id')
#'     longdat.2<-merge(jointdat.2$longitudinal,jointdat.2$baseline,by='id')
#'     longdat.3<-merge(jointdat.3$longitudinal,jointdat.3$baseline,by='id')
#'     longdat.4<-merge(jointdat.4$longitudinal,jointdat.4$baseline,by='id')
#'     longdat.5<-merge(jointdat.5$longitudinal,jointdat.5$baseline,by='id')
#'
#'     survdat.1<-merge(jointdat.1$survival,jointdat.1$baseline,by='id')
#'     survdat.2<-merge(jointdat.2$survival,jointdat.2$baseline,by='id')
#'     survdat.3<-merge(jointdat.3$survival,jointdat.3$baseline,by='id')
#'     survdat.4<-merge(jointdat.4$survival,jointdat.4$baseline,by='id')
#'     survdat.5<-merge(jointdat.5$survival,jointdat.5$baseline,by='id')
#'
#'
#'     #fit with current value sharing structure
#'     library('nlme', 'JM')
#'     lmefit1<-lme(fixed = Y ~ 1 + time + treat, data = longdat.1,
#'                    random = ~ 1 + time|id, control=lmeControl(opt = 'optim'))
#'     survfit1<-coxph(Surv(survtime, cens) ~ treat, data = survdat.1, x = TRUE)
#'     JMfit1<-jointModel(lmeObject = lmefit1, survObject = survfit1,
#'                    timeVar = 'time',
#'                    parameterization = 'value', method = 'spline-PH-aGH')
#'
#'     lmefit2<-lme(fixed = Y ~ 1 + time + treat, data = longdat.2,
#'                    random = ~ 1 + time|id, control=lmeControl(opt = 'optim'))
#'     survfit2<-coxph(Surv(survtime, cens) ~ treat, data = survdat.2, x = TRUE)
#'     JMfit2<-jointModel(lmeObject = lmefit2, survObject = survfit2,
#'                    timeVar = 'time',
#'                    parameterization = 'value', method = 'spline-PH-aGH')
#'
#'     lmefit3<-lme(fixed = Y ~ 1 + time + treat, data = longdat.3,
#'                    random = ~ 1 + time|id, control=lmeControl(opt = 'optim'))
#'     survfit3<-coxph(Surv(survtime, cens) ~ treat, data = survdat.3, x = TRUE)
#'     JMfit3<-jointModel(lmeObject = lmefit3, survObject = survfit3,
#'                    timeVar = 'time',
#'                    parameterization = 'value', method = 'spline-PH-aGH')
#'
#'     lmefit4<-lme(fixed = Y ~ 1 + time + treat, data = longdat.4,
#'                    random = ~ 1 + time|id, control=lmeControl(opt = 'optim'))
#'     survfit4<-coxph(Surv(survtime, cens) ~ treat, data = survdat.4, x = TRUE)
#'     JMfit4<-jointModel(lmeObject = lmefit4, survObject = survfit4,
#'                    timeVar = 'time',
#'                    parameterization = 'value', method = 'spline-PH-aGH')
#'
#'     lmefit5<-lme(fixed = Y ~ 1 + time + treat, data = longdat.5,
#'                    random = ~ 1 + time|id, control=lmeControl(opt = 'optim'))
#'     survfit5<-coxph(Surv(survtime, cens) ~ treat, data = survdat.5, x = TRUE)
#'     JMfit5<-jointModel(lmeObject = lmefit5, survObject = survfit5,
#'                    timeVar = 'time',
#'                    parameterization = 'value', method = 'spline-PH-aGH')
#'
#'     JMfits<-list(JMfit1,
#'                  JMfit2,
#'                  JMfit3,
#'                  JMfit4,
#'                  JMfit5)
#'
#'     names(JMfits)<-c('JMfit1',
#'                      'JMfit2',
#'                      'JMfit3',
#'                      'JMfit4',
#'                      'JMfit5')
#'
#'     #perform the MA of the supplied fits
#'     MAJMfits<-jointmeta2(fits = JMfits, longpar = c('time', 'treat1'),
#'                   survpar = 'treat1', assoc = TRUE,
#'                   studynames = c('Study 1','Study 2', 'Study 3', 'Study 4',
#'                   'Study 5'))
#'
#'     #produce forest plots
#'     library(meta)
#'     forest(MAJMfits$longMA$treat1)
#'
#'
#'
jointmeta2 <- function(fits, SE = NULL, longpar = NULL, survpar = NULL,
                       assoc = TRUE, studynames = NULL) {
  if (class(fits) != "list") {
    stop("fits should be supplied a list of joint model fits")
  }
  if (is.null(longpar) == FALSE) {
    if (class(longpar) != "character") {
      stop("Please supply the longitudinal and survival
           parameters to meta-analyse as character strings")
    }
    }
  if (is.null(survpar) == FALSE) {
    if (class(survpar) != "character") {
      stop("Please supply the longitudinal and survival parameters
           to meta-analyse as character strings")
    }
    }
  jointclass <- unlist(lapply(fits, class))
  if (length(unique(jointclass)) > 1) {
    stop("Some of the joint modelling fits are different classes -
         consider subgrouping")
  }
  if (is.null(studynames) == FALSE) {
    if (class(studynames) != "character") {
      stop("Studynames should be supplied as a character vector")
    }
  } else {
    studynames <- as.character(1:length(fits))
  }
  jointclass <- unique(jointclass)
  if (jointclass == "joint") {
    if (is.null(SE)) {
      stop("No SE suplied to the function")
    }

    if (length(SE) != length(fits)) {
      stop("Supplied fits and SE differ in length")
    }
    if (length(unique(unlist(lapply(fits, function(u) {
      as.character(u$call$model)
    })))) > 1) {
      stop("Some of the joint model fits have differing random
           effects structures")
    }
    if (length(unique(unlist(lapply(fits, function(u) {
      if (is.null(u$call$sepassoc)) {
        FALSE
      } else {
        as.character(u$call$sepassoc)
      }
    })))) > 1) {
      stop("Some of the joint model fits have differing association structures")
    }
    for (i in 1:length(fits)) {
      if (is.null(longpar) == FALSE) {
        if (length(which(!(longpar %in% rownames(fits[[i]]$coefficients$fixed$longitudinal)))) >
            0) {
          stop("Specified longpar covariate not found
               in some longitudinal formulas of some joint fits")
        }
        }
      if (is.null(survpar) == FALSE) {
        if (length(which(!(survpar %in% names(fits[[i]]$coefficients$fixed$survival)))) >
            0) {
          stop("Specified survpar covariate not found
               in some survival formulas of some joint fits")
        }
        }
      }
    if (is.null(longpar) == FALSE) {
      longest <- longpres <- data.frame(matrix(NA, ncol = length(longpar),
                                               nrow = length(fits)))
      names(longest) <- names(longpres) <- longpar
    }
    if (is.null(survpar) == FALSE) {
      survest.direct <- survpres.direct <- data.frame(matrix(NA,
                                                             ncol = length(survpar), nrow = length(fits)))
      names(survest.direct) <- names(survpres.direct) <- survpar
    }
    if (assoc) {
      assocest <- assocpres <- data.frame(matrix(NA, ncol = length(fits[[1]]$coefficients$latent),
                                                 nrow = length(fits)))
      names(assocest) <- names(assocpres) <- names(fits[[1]]$coefficients$latent)
    }
    samplesizes <- rep(NA, length(fits))
    results <- list()
    for (i in 1:length(fits)) {
      tempfit <- fits[[i]]
      samplesizes[i] <- nrow(fits[[i]]$data$survival)
      if (is.null(longpar) == FALSE) {
        longest[i, ] <- tempfit$coefficients$fixed$longitudinal[unlist(lapply(1:length(longpar),
                                                                              function(u) {
                                                                                grep(longpar[u], rownames(tempfit$coefficients$fixed$longitudinal))
                                                                              })), ]
        longSE <- SE[[i]][which(SE[[i]][, 1] == "Longitudinal"):(which(SE[[i]][,
                                                                               1] == "Survival") - 1), ]
        longSE <- longSE[unlist(lapply(1:length(longpar), function(u) {
          grep(longpar[u], longSE[, 2])
        })), ]
        longpres[i, ] <- as.numeric(as.character(longSE[, 4]))
        longMA <- lapply(1:ncol(longest), function(u) {
          metagen(TE = longest[, u], seTE = longpres[, u], studlab = studynames,
                  sm = "MD", comb.fixed = TRUE, comb.random = TRUE)
        })
        names(longMA) <- longpar
        results$longMA <- longMA
      }
      if (is.null(survpar) == FALSE) {
        survest.direct[i, ] <- tempfit$coefficients$fixed$survival[unlist(lapply(1:length(survpar),
                                                                                 function(u) {
                                                                                   grep(survpar[u], names(tempfit$coefficients$fixed$survival))
                                                                                 }))]
        survSE <- SE[[i]][which(SE[[i]][, 1] == "Survival"):(which(SE[[i]][,
                                                                           1] == "Association") - 1), ]
        survSE <- survSE[unlist(lapply(1:length(survpar), function(u) {
          grep(survpar[u], survSE[, 2])
        })), ]
        survpres.direct[i, ] <- as.numeric(as.character(survSE[,
                                                               4]))
        survMA.direct <- lapply(1:ncol(survest.direct), function(u) {
          metagen(TE = survest.direct[, u], seTE = survpres.direct[,
                                                                   u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                  comb.random = TRUE)
        })
        names(survMA.direct) <- survpar
        results$survMA.direct <- survMA.direct
      }
      if (assoc) {
        assocest[i, ] <- tempfit$coefficients$latent
        assocSE <- SE[[i]][which(SE[[i]][, 1] == "Association"):(which(SE[[i]][,
                                                                               1] == "Variance") - 1), ]
        assocpres[i, ] <- as.numeric(as.character(assocSE[, 4]))
        assocMA <- lapply(1:ncol(assocest), function(u) {
          metagen(TE = assocest[, u], seTE = assocpres[, u], studlab = studynames,
                  sm = "MD", comb.fixed = TRUE, comb.random = TRUE)
        })
        names(assocMA) <- names(assocest)
        results$assocMA <- assocMA
      }
    }
    return(results)
      } else if (jointclass == "jointModel") {
        if (length(unique(unlist(lapply(fits, function(u) {
          as.character(u$parameterization)
        })))) > 1) {
          stop("Some of the joint model fits have differing association structures")
        }
        if (length(unique(unlist(lapply(fits, function(u) {
          temp <- unlist(strsplit(u$method, "-"))[[2]]
        })))) > 1) {
          stop("Different methods used (PH/AFT) for survival sub-model fit")
        }
        parametrization <- lapply(fits, function(u) {
          u$call$parameterization
        })
        uniqueparameterization <- unique(unlist(parametrization))
        randvar <- lapply(fits, function(u) {
          colnames(ranef(u))
        })
        uniquerandvar <- unique(unlist(randvar))
        if (uniqueparameterization %in% c("value", "both")) {
          for (i in 1:length(fits)) {
            for (j in 1:length(uniquerandvar)) {
              if (!(uniquerandvar[j] %in% randvar[[i]])) {
                separated <- unlist(strsplit(uniquerandvar[j], "\\:"))
                perms <- apply(permutations(n = length(separated),
                                            r = length(separated), v = separated), 1, paste,
                               collapse = ":")
                contained <- FALSE
                for (k in 1:length(perms)) {
                  if (perms[k] %in% randvar[[i]]) {
                    contained <- TRUE
                  }
                }
                if (contained == FALSE) {
                  stop("Random effect specification differs
                       between supplied model fits")
                }
                }
              }
            }
          }
        fixedvar <- lapply(fits, function(u) {
          names(fixef(u))
        })
        uniquefixedvar <- unique(unlist(fixedvar))
        if (uniqueparameterization %in% c("value", "both")) {
          for (i in 1:length(fits)) {
            for (j in 1:length(uniquefixedvar)) {
              if (!(uniquefixedvar[j] %in% fixedvar[[i]])) {
                separated <- unlist(strsplit(uniquefixedvar[j], "\\:"))
                perms <- apply(permutations(n = length(separated),
                                            r = length(separated), v = separated), 1, paste,
                               collapse = ":")
                contained <- FALSE
                for (k in 1:length(perms)) {
                  if (perms[k] %in% fixedvar[[i]]) {
                    contained <- TRUE
                  }
                }
                if (contained == FALSE) {
                  stop("Fixed effect specification differs
                       between supplied model fits")
                }
                }
              }
            }
          }
        if (uniqueparameterization %in% c("slope", "both")) {
          dform.fixed <- lapply(fits, function(u) {
            u$derivForm$fixed
          })
          dform.indFixed <- lapply(fits, function(u) {
            names(u$coefficients$betas)[u$derivForm$indFixed]
          })
          dform.random <- lapply(fits, function(u) {
            u$derivForm$random
          })
          dform.indRandom <- lapply(fits, function(u) {
            rownames(u$coefficients$D)[u$derivForm$indRandom]
          })
          uniquedform.indfixed <- unique(unlist(dform.indFixed))
          uniquedform.indrandom <- unique(unlist(dform.indRandom))
          for (i in 1:length(fits)) {
            for (j in 1:length(uniquedform.indfixed)) {
              if (!(uniquedform.indfixed[j] %in% dform.indFixed[[i]])) {
                separated <- unlist(strsplit(uniquedform.indfixed[j],
                                             "\\:"))
                perms <- apply(permutations(n = length(separated),
                                            r = length(separated), v = separated), 1, paste,
                               collapse = ":")
                contained <- FALSE
                for (k in 1:length(perms)) {
                  if (perms[k] %in% dform.indFixed[[i]]) {
                    contained <- TRUE
                  }
                }
                if (contained == FALSE) {
                  stop("DerivForm fixed effect specification differs
                       between supplied model fits")
                }
                }
              }
            }
          for (i in 1:length(fits)) {
            for (j in 1:length(uniquedform.indrandom)) {
              if (!(uniquedform.indrandom[j] %in% dform.indRandom[[i]])) {
                separated <- unlist(strsplit(uniquedform.indrandom[j],
                                             "\\:"))
                perms <- apply(permutations(n = length(separated),
                                            r = length(separated), v = separated), 1, paste,
                               collapse = ":")
                contained <- FALSE
                for (k in 1:length(perms)) {
                  if (perms[k] %in% dform.indRandom[[i]]) {
                    contained <- TRUE
                  }
                }
                if (contained == FALSE) {
                  stop("DerivForm random effect specification differs between
                       supplied model fits")
                }
                }
              }
            }
          }
        for (i in 1:length(fits)) {
          if (is.null(longpar) == FALSE) {
            if (length(which(!(longpar %in% names(fits[[i]]$coefficients$betas)))) >
                0) {
              stop("Specified longpar covariate not found in some
                   longitudinal formulas of some joint fits")
            }
            }
          if (is.null(survpar) == FALSE) {
            if (length(which(!(survpar %in% names(fits[[i]]$coefficients$gammas)))) >
                0) {
              stop("Specified survpar covariate not found in some
                   survival formulas of some joint fits")
            }
            }
          }
        if (is.null(longpar) == FALSE) {
          longest <- longpres <- data.frame(matrix(NA, ncol = length(longpar),
                                                   nrow = length(fits)))
          names(longest) <- names(longpres) <- longpar
        }
        if (is.null(survpar) == FALSE) {
          survest.direct <- survpres.direct <- data.frame(matrix(NA,
                                                                 ncol = length(survpar), nrow = length(fits)))
          names(survest.direct) <- names(survpres.direct) <- survpar
          if (uniqueparameterization %in% c("value", "both")) {
            survest.indirect.val <- data.frame(matrix(NA, ncol = length(survpar),
                                                      nrow = length(fits)))
            names(survest.indirect.val) <- survpar
          }
          if (uniqueparameterization %in% c("slope", "both")) {
            survest.indirect.slope <- data.frame(matrix(NA, ncol = length(survpar),
                                                        nrow = length(fits)))
            names(survest.indirect.slope) <- survpar
          }
          survest.overall <- survpres.overall <- data.frame(matrix(NA,
                                                                   ncol = length(survpar), nrow = length(fits)))
          names(survest.overall) <- names(survpres.overall) <- survpar
        }
        if (assoc) {
          numassoc <- 1
          namesassoc <- c("Assoct")
          if (uniqueparameterization == "both") {
            numassoc <- 2
            namesassoc <- c(namesassoc, "Assoct.s")
          }
          assocest <- assocpres <- data.frame(matrix(NA, ncol = numassoc,
                                                     nrow = length(fits)))
          names(assocest) <- names(assocpres) <- namesassoc
        }
        samplesizes <- rep(NA, length(fits))
        results <- list()
        for (i in 1:length(fits)) {
          samplesizes[i] <- fits[[i]]$n
          confints.long <- confint(fits[[i]])[grepl("Y.", rownames(confint(fits[[i]]))),
                                              ]
          confints.surv <- confint(fits[[i]])[grepl("T.", rownames(confint(fits[[i]]))),
                                              ]
          assocest[i, ] <- confints.surv[grepl("Assoct", rownames(confints.surv)),
                                         2]
          assocpres[i, ] <- summary(fits[[i]])$"CoefTable-Event"[, 2][grepl("Assoct",
                                                                            rownames(summary(fits[[i]])$"CoefTable-Event"))]
          if (is.null(longpar) == FALSE) {
            longest[i, ] <- confints.long[match(paste("Y.", longpar,
                                                      sep = ""), rownames(confints.long)), 2]
            longpres[i, ] <- summary(fits[[i]])$"CoefTable-Long"[,
                                                                 2][match(longpar, names(summary(fits[[i]])$"CoefTable-Long"[,
                                                                                                                             2]))]
            longMA <- lapply(1:ncol(longest), function(u) {
              metagen(TE = longest[, u], seTE = longpres[, u], studlab = studynames,
                      sm = "MD", comb.fixed = TRUE, comb.random = TRUE)
            })
            names(longMA) <- longpar
            results$longMA <- longMA
          }
          if (is.null(survpar) == FALSE) {
            survest.direct[i, ] <- confints.surv[match(paste("T.",
                                                             survpar, sep = ""), rownames(confints.surv)), 2]
            survpres.direct[i, ] <- summary(fits[[i]])$"CoefTable-Event"[,
                                                                         2][match(survpar, names(summary(fits[[i]])$"CoefTable-Event"[,
                                                                                                                                      2]))]
            if (uniqueparameterization %in% c("value", "both")) {
              indirect.val.id <- match(survpar, names(fits[[i]]$coefficients$betas))
              for (j in 1:length(indirect.val.id)) {
                if (is.na(indirect.val.id[j])) {
                  survest.indirect.val[i, j] <- 0
                } else {
                  survest.indirect.val[i, j] <- assocest[i, 1] * fits[[i]]$coefficients$betas[indirect.val.id[j]]
                }
              }
            }
            if (uniqueparameterization %in% c("slope", "both")) {
              slopetempnames <- slopetempnameslong <- names(fits[[i]]$coefficients$betas)[fits[[i]]$derivForm$indFixed]
              for (j in 1:length(slopetempnames)) {
                if (slopetempnames[[j]] == fits[[i]]$timeVar) {
                  slopetempnames[[j]] <- "1"
                } else if (grepl(paste(fits[[i]]$timeVar, ":", sep = ""),
                                 slopetempnames[j])) {
                  slopetempnames[j] <- gsub(paste(fits[[i]]$timeVar,
                                                  ":", sep = ""), "\\1", slopetempnames[j])
                } else if (grepl(paste(":", fits[[i]]$timeVar, sep = ""),
                                 slopetempnames[j])) {
                  slopetempnames[j] <- gsub(paste(":", fits[[i]]$timeVar,
                                                  sep = ""), "\\1", slopetempnames[j])
                }
              }
              for (j in 1:length(survpar)) {
                if (survpar[j] %in% slopetempnames) {
                  loc <- fits[[i]]$derivForm$indFixed[which(slopetempnames ==
                                                              survpar[j])]
                  if (uniqueparameterization == "both") {
                    survest.indirect.slope[i, j] <- assocest[i, 2] *
                      fits[[i]]$coefficients$betas[loc]
                  } else {
                    survest.indirect.slope[i, j] <- assocest[i, 1] *
                      fits[[i]]$coefficients$betas[loc]
                  }
                } else {
                  survest.indirect.slope[i, j] <- 0
                }
              }
            }
            if (uniqueparameterization == "value") {
              survest.overall[i, ] <- survest.direct[i, ] + survest.indirect.val[i,
                                                                                 ]
              for (j in 1:length(survpar)) {
                if (is.na(indirect.val.id[[j]]) == FALSE) {
                  namesselector <- c(paste("T.", survpar[j], sep = ""),
                                     "T.alpha", paste("Y.", names(fits[[i]]$coefficients$betas)[indirect.val.id[j]],
                                                      sep = ""))
                  selector <- match(namesselector, rownames(vcov(fits[[i]])))
                  covmat <- vcov(fits[[i]])[selector, selector]
                  survpres.overall[i, j] <- deltamethod(~x1 + (x2 *
                                                                 x3), c(survest.direct[i, j], assocest[i, 1], fits[[i]]$coefficients$betas[indirect.val.id[j]]),
                                                        covmat)
                } else {
                  survpres.overall[i, j] <- survpres.direct[i, j]
                }
              }
              survMA.direct <- lapply(1:ncol(survest.direct), function(u) {
                metagen(TE = survest.direct[, u], seTE = survpres.direct[,
                                                                         u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.direct) <- survpar
              survMA.overall <- lapply(1:ncol(survest.overall), function(u) {
                metagen(TE = survest.overall[, u], seTE = survpres.overall[,
                                                                           u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.overall) <- survpar
              results$survMA.direct <- survMA.direct
              results$survMA.overall <- survMA.overall
            } else if (uniqueparameterization == "slope") {
              survest.overall[i, ] <- survest.direct[i, ] + survest.indirect.slope[i,
                                                                                   ]
              for (j in 1:length(survpar)) {
                if (survpar[j] %in% slopetempnames) {
                  loc <- which(slopetempnames == survpar[j])
                  namesselector <- c(paste("T.", survpar[j], sep = ""),
                                     "T.alpha.s", paste("Y.", names(fits[[i]]$coefficients$betas)[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                          slopetempnameslong[loc])], sep = ""))
                  selector <- match(namesselector, rownames(vcov(fits[[i]])))
                  covmat <- vcov(fits[[i]])[selector, selector]
                  survpres.overall[i, j] <- deltamethod(~x1 + (x2 *
                                                                 x3), c(survest.direct[i, j], assocest[i, 1], fits[[i]]$coefficients$betas[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                                                                   slopetempnameslong[loc])]), covmat)
                } else {
                  survpres.overall[i, j] <- survpres.direct[i, j]
                }
              }
              survMA.direct <- lapply(1:ncol(survest.direct), function(u) {
                metagen(TE = survest.direct[, u], seTE = survpres.direct[,
                                                                         u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.direct) <- survpar
              survMA.overall <- lapply(1:ncol(survest.overall), function(u) {
                metagen(TE = survest.overall[, u], seTE = survpres.overall[,
                                                                           u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.overall) <- survpar
              results$survMA.direct <- survMA.direct
              results$survMA.overall <- survMA.overall
            } else if (uniqueparameterization == "both") {
              survest.overall[i, ] <- survest.direct[i, ] + survest.indirect.val[i,
                                                                                 ] + survest.indirect.slope[i, ]
              for (j in 1:length(survpar)) {
                if (is.na(indirect.val.id[[j]]) == FALSE) {
                  if (survpar[j] %in% slopetempnames) {
                    loc <- which(slopetempnames == survpar[j])
                    namesselector <- c(paste("T.", survpar[j], sep = ""),
                                       "T.alpha", paste("Y.", names(fits[[i]]$coefficients$betas)[indirect.val.id[j]],
                                                        sep = ""), "T.alpha.s", paste("Y.", names(fits[[i]]$coefficients$betas)[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                                                        slopetempnameslong[loc])], sep = ""))
                    selector <- match(namesselector, rownames(vcov(fits[[i]])))
                    covmat <- vcov(fits[[i]])[selector, selector]
                    survpres.overall[i, j] <- deltamethod(~x1 + (x2 *
                                                                   x3) + (x4 * x5), c(survest.direct[i, j], assocest[i,
                                                                                                                     1], fits[[i]]$coefficients$betas[indirect.val.id[j]],
                                                                                      assocest[i, 2], fits[[i]]$coefficients$betas[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                                                           slopetempnameslong[loc])]), covmat)

                  } else {
                    namesselector <- c(paste("T.", survpar[j], sep = ""),
                                       "T.alpha", paste("Y.", names(fits[[i]]$coefficients$betas)[indirect.val.id[j]],
                                                        sep = ""))
                    selector <- match(namesselector, rownames(vcov(fits[[i]])))
                    covmat <- vcov(fits[[i]])[selector, selector]
                    survpres.overall[i, j] <- deltamethod(~x1 + (x2 *
                                                                   x3), c(survest.direct[i, j], assocest[i, 1],
                                                                          fits[[i]]$coefficients$betas[indirect.val.id[j]]),
                                                          covmat)
                  }
                } else {
                  if (survpar[j] %in% slopetempnames) {
                    loc <- which(slopetempnames == survpar[j])
                    namesselector <- c(paste("T.", survpar[j], sep = ""),
                                       "T.alpha.s", paste("Y.", names(fits[[i]]$coefficients$betas)[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                            slopetempnameslong[loc])], sep = ""))
                    selector <- match(namesselector, rownames(vcov(fits[[i]])))
                    covmat <- vcov(fits[[i]])[selector, selector]
                    survpres.overall[i, j] <- deltamethod(~x1 + (x2 *
                                                                   x3), c(survest.direct[i, j], assocest[i, 2],
                                                                          fits[[i]]$coefficients$betas[which(names(fits[[i]]$coefficients$betas) %in%
                                                                                                               slopetempnameslong[loc])]), covmat)
                  } else {
                    survpres.overall[i, j] <- survpres.direct[i, j]
                  }
                }
              }
              survMA.direct <- lapply(1:ncol(survest.direct), function(u) {
                metagen(TE = survest.direct[, u], seTE = survpres.direct[,
                                                                         u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.direct) <- survpar
              survMA.overall <- lapply(1:ncol(survest.overall), function(u) {
                metagen(TE = survest.overall[, u], seTE = survpres.overall[,
                                                                           u], studlab = studynames, sm = "MD", comb.fixed = TRUE,
                        comb.random = TRUE)
              })
              names(survMA.overall) <- survpar
              results$survMA.direct <- survMA.direct
              results$survMA.overall <- survMA.overall
            } else {
              stop("parametrization of joint model not value,
               slope or both and so not recognised")
            }
          }
        }
        if (assoc) {
          assocMA <- lapply(1:ncol(assocest), function(u) {
            metagen(TE = assocest[, u], seTE = assocpres[, u], studlab = studynames,
                    sm = "MD", comb.fixed = TRUE, comb.random = TRUE)
          })
          names(assocMA) <- names(assocest)
          results$assocMA <- assocMA
        }
        class(results) <- "jointmeta2"
        return(results)
          } else {
            stop("Function currently supports only joineR or JM model fits")
          }
        }
