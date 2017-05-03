#' Function for survival starting values
#'
#' Internal function to estimate the starting values for the EM algorithm for
#'     the survival sub-model
#'
#' @param survdat2 the survival data with original variables (factors and
#' interaction terms not expanded), ordered by increasing survival time
#' @inheritParams jointmeta1
#' @inheritParams EMalgRandprop
#'
#' @return A list of the results from the initial survival fit is returned.
#'     This list contains the following elements: \describe{
#'
#'     \item{\code{beta2}}{a vector of the estimated coefficients for fixed
#'     effects included in the survival sub-model. If there are no fixed effects
#'     included in the survival sub-model then this returns \code{NULL}}.
#'
#'     \item{\code{haz}}{the estimate of the baseline hazard estimated from the
#'     separate survival model.  If \code{strat = TRUE} then this is a list of
#'     length equal to the number of studies in the dataset, each element of
#'     which is a vector of length equal to the number of events in each study.
#'     If \code{strat = FALSE} then this is a vector of length equal to the
#'     number of events in the entire dataset.}
#'
#'     \item{\code{rs}}{the number of events that have occured by the individual
#'     in question's survival time. If \code{strat = TRUE} then this is a list
#'     of length equal to the number of studies in the dataset, each element of
#'     which is a vector of length equal to the number of individuals in each
#'     study. If \code{strat = FALSE} then this is a vector of length equal to
#'     the number of individuals in the entire dataset.}
#'
#'     \item{\code{sf}}{the survival times where at least one event was
#'     observed. If \code{strat = TRUE} then this is a list of
#'     length equal to the number of studies in the dataset, each element of
#'     which is a vector of length equal to the number of events in each study.
#'     If \code{strat = FALSE} then this is a vector of length equal to the
#'     number of events in the entire dataset.}
#'
#'     \item{\code{nev}}{the number of events that occur at each unique event
#'     time.  If \code{strat = TRUE} then this is a list of
#'     length equal to the number of studies in the dataset, each element of
#'     which is a vector of length equal to the number of events in each study.
#'     If \code{strat = FALSE} then this is a vector of length equal to the
#'     number of events in the entire dataset.}
#'
#'     \item{\code{log.like.surv}}{the value of the log-likelihood from the
#'     separate survival analysis}
#'
#'     \item{\code{modelfit}}{the initial survival model fit, fitted using
#'     the \code{\link[survival]{coxph}} function from the
#'     \code{\link[survival]{survival}} package.}
#'
#'     }
#'
#'
#' @keywords internal
#' @import survival
survst <- function(survdat, surv.formula, survdat2, strat, study.name = NULL) {
  survdat2 <- survdat2[order(survdat2[, 2]), ]
  n <- length(survdat[, 2])
  s <- survdat[, 2]
  cen <- survdat[, 3]
  study<-survdat[, 4]
  ids.surv<-survdat[, 1]
  p2 <- dim(survdat)[2] - 4
  studies<-unique(study)
  numstudies<-length(studies)
  if(strat == FALSE) {
    surv.start <- coxph(surv.formula, data = survdat, x = TRUE)
    surv.start.f <- survfit(surv.start)
    sf <- surv.start.f$time[surv.start.f$n.event != 0]
    nf <- length(sf)
    nev <- surv.start.f$n.event[surv.start.f$n.event != 0]
    if (p2 > 0) {
      haz <- coxph.detail(surv.start)$hazard
    }else {
      haz <- surv.start.f$n.event/surv.start.f$n.risk
      haz <- haz[surv.start.f$n.event > 0]
    }
    rs <- rep(1:nf, c(diff(match(sf, s)), n + 1 - match(sf, s)[nf]))
    if(cen[1] == 0) {
      rs <- c(rep(0, which(cen == 1)[1]-1), rs)
    }
    names(rs) <- survdat[,1]
    beta2 <- coef(surv.start)
    ll <- surv.start$loglik - sum(cen)
    list(beta2 = beta2, haz = haz, rs = rs, sf = sf, nev = nev,
         log.like.surv = ll, modelfit = surv.start)
  }else{
    ids.bystudy <- lapply(1:numstudies, function(u) {
      ids.surv[which(survdat[, 4] == studies[u])]
      })
    names(ids.bystudy) <- studies
    n.bystudy <- unlist(lapply(1:numstudies, function(u) {
      length(ids.bystudy[[u]])
      }))
    cen.bystudy <- lapply(1:numstudies, function(u) {
      cen[which(ids.surv%in%ids.bystudy[[u]])]
      })
    names(cen.bystudy) <- studies
    s.bystudy <- lapply(1:numstudies, function(u) {
      s[which(ids.surv%in%ids.bystudy[[u]])]
      })
    names(s.bystudy) <- studies
    surv.formula.strat <- as.character(surv.formula)
    surv.formula.strat[3] <- paste(surv.formula.strat[3],
                                   "+strata(", study.name, ")",
                                   sep = "", collapse = "")
    surv.formula.strat <- as.formula(paste(surv.formula.strat[2],
                                           surv.formula.strat[1],
                                           surv.formula.strat[3], sep = ""))
    surv.start <- coxph(surv.formula.strat, data = survdat, x = TRUE)
    surv.start.f <- survfit(surv.start)
    data <- data.frame(cbind(basehaz(surv.start),
                             surv.start.f$n.event, surv.start.f$n.risk))
    data.bystudy <- lapply(1:numstudies, function(u) {
      data[which(data[, 3] == studies[u]), ]
      })
    names(data.bystudy) <- studies
    sf.bystudy <- lapply(1:numstudies, function(u) {
      data.bystudy[[u]][which(data.bystudy[[u]][, 4] != 0), 2]
      })
    names(sf.bystudy) <- studies
    nf <- unlist(lapply(1:numstudies, function(u) {length(sf.bystudy[[u]])}))
    nev.bystudy <- lapply(1:numstudies, function(u) {
      data.bystudy[[u]][which(data.bystudy[[u]][, 4] != 0), 4]
      })
    names(nev.bystudy) <- studies
    if (p2 > 0) {
      haz <- lapply(1:numstudies, function(u) {
        data.bystudy[[u]][which(data.bystudy[[u]][, 4] != 0), 1]
        })
    }else {
      haz <- lapply(1:numstudies, function(u) {
        datatemp <- data.bystudy[[u]]
        temp <- datatemp[, 4]/datatemp[, 5]
        temp[which(temp!=0)]
      })
    }
    names(haz) <- studies
    rs<-lapply(1:numstudies, function(u) {
      eventslocs <- which(data.bystudy[[u]]$surv.start.f.n.event!=0)
      eventcounter <- 0
      loccounter <- 1
      temp <- c()
      for(i in 1:length(eventslocs)) {
        temp <- c(temp,rep(eventcounter, data.bystudy[[u]]$surv.start.f.n.risk[
          (loccounter)]-data.bystudy[[u]]$surv.start.f.n.risk[
            (eventslocs[[i]])]))
        eventcounter <- eventcounter + 1
        loccounter <- eventslocs[i]
      }
      temp <- c(temp, rep(length(eventslocs),length(ids.bystudy[[u]])-length(temp)))
      names(temp) <- ids.bystudy[[u]]
      temp
    })

    names(rs) <- studies
    beta2 <- coef(surv.start)
    ll <- surv.start$loglik - sum(as.numeric(do.call(c, cen.bystudy)))
    list(beta2 = beta2, haz = haz, rs = rs, sf = sf.bystudy, nev = nev.bystudy,
         log.like.surv = ll, modelfit = surv.start)
  }
}
