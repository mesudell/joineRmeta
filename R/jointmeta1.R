#' One stage joint meta function
#'
#' Function to allow a one stage joint model (data from all studies analysed in
#' one model)  to be fitted to data from multiple studies.  The function allows
#' one longitudinal and one time-to-event outcome, and can accommodate baseline
#' hazard stratified or not stratified by study, as well as random effects at
#' the individual level and the study level.  Currently only zero mean random
#' effects only proportional association supported - see Wulfsohn and Tsiatis
#' 1997
#'
#' @param data an object of class jointdata containing the variables named in
#'   the model formulae
#' @param long.formula a formula object with the response varaible, and the
#'   covariates to include in the longitudinal sub-model
#' @param long.rand.ind a vector of character strings to indicate what variables
#'   to assign individual level random effects to.  A maximum of three
#'   individual level random effects can be assigned.  To assign a random
#'   intercept include 'int' in the vector.  To not include an individual level
#'   random intercept include 'noint' in the vector.  For example to fit a model
#'   with individual level random intercept and random slope set
#'   \code{long.rand.ind = c('int', 'time')}, where \code{'time'} is the
#'   longitudinal time variable in the \code{data}.
#' @param long.rand.stud a vector of character strings to indicate what
#'   variables to assign study level random effects to.  If no study level
#'   random effects then this either not specified in function call or set to
#'   \code{NULL}.  If a study level random intercept is required, include the
#'   name of the study membership variable for example \code{long.rand.stud =
#'   'study'}.
#' @param sharingstrct currently must be set to \code{'randprop'}.  This gives a
#'   model that shares the zero mean random effects (at both individual and
#'   study level if specified) between the sub-models.  Separate association
#'   parameters are calculated for the linear combination of random effects at
#'   each level.  There are plans to expand to more sharing structures in the
#'   future.
#' @param surv.formula a formula object with the survival time, censoring
#'   indicator and the covariates to include in the survival sub-model.  The
#'   response must be a survival object as returned by the
#'   \code{\link[survival]{Surv}} function.
#' @param gpt the number of quadrature points across which the integration with
#'   respect to the random effects will be performed.  If random effects are
#'   specified at both the individual and the study level, the same number of
#'   quadrature points is used in both cases.  Defaults to \code{gpt = 5}.
#' @param lgpt the number of quadrature points which the log-likelihood is
#'   evaluated over following a model fit.  This defaults to \code{lgpt = 7}.
#' @param max.it the maximum number of iterations of the EM algorithm that the
#'   function will perform.  Defaults to \code{max.it = 350} although more
#'   iterations could be required for large complex datasets.
#' @param tol the tolerance level used to determine convergence in the EM
#'   algorithm.  Defaults to \code{tol = 0.001}.
#' @param study.name a character string denoting the name of the variable in the
#'   baseline dataset in \code{data} holding study membership, for example
#'   \code{study.name = 'study'}.
#' @param strat logical value: if \code{TRUE} then the survival sub-model is
#'   calculated with a baseline stratified by study.  Otherwise baseline is
#'   unstratified
#' @param longsep logical value: if \code{TRUE} then parameter estimates, model
#'   fit and the log-likelihood from a separate linear mixed model analysis of
#'   the longitudinal data are returned (see the \code{\link[lme4]{lmer}}
#'   function).  The separate longitudinal model fit has the same specification
#'   as the longitudinal sub-model of the joint model.
#' @param survsep logical value: if \code{TRUE} then parameter estimates, model
#'   fit and log-likelihood from a separate analysis of the survival data using
#'   the Cox Proportional Hazards model are returned (see
#'   \code{\link[survival]{coxph}} function for more details).  This survival
#'   fit has the same specification (apart from the association structure) as
#'   the survival sub-model in the joint model.
#' @param bootrun logical value: if \code{TRUE} then the log-likelihood for the
#'   model is not calculated.  This option is available so that when
#'   bootstrapping to obtain standard errors, as the log-likelihood is not
#'   needed, it is not calculated, thus speeding up the bootstrapping process.
#' @param print.detail logical value: if \code{TRUE} then details of the
#'   parameter estimates at each iteration of the EM algorithm are printed to
#'   the console.
#'
#' @section Details: The \code{jointmeta1} function fits a one stage joint model
#'   to survival and longitudinal data from multiple studies. This model is an
#'   extension of the model proposed by Wulfsohn and Tsiatis (1997).  The model
#'   must contain at least one individual level random effect (specified using
#'   the \code{long.rand.ind} argument).  The model can also contain study level
#'   random effects (specified using the  \code{long.rand.stud} argument), which
#'   can differ from the individual level random effects. The maximum number of
#'   random effects that can be specified at each level is three. Note that the
#'   fitting and bootstrapping time increases as the number of included random
#'   effects increases.  The model can also include a baseline hazard stratified
#'   by study, or can utilise a common baseline across the studies in the
#'   dataset.  Interaction terms can be specified in either the longitudinal or
#'   the survival sub-model.
#'
#'   The longitudinal sub-model is a mixed effects model.  If both individual
#'   level and study level random effects are included in the function call,
#'   then the sub-model has the following format:
#'
#'   \deqn{Y_{kij} = X_{1kij}\beta_{1} + Z^{(2)}_{1kij}b^{(2)}_{ki} +
#'   Z^{(3)}_{1kij}b^{(3)}_{k} + \epsilon_{kij}}
#'
#'   Otherwise, if only individual level random effects are included in the
#'   function call, then the longitudinal sub-model has the following format:
#'
#'   \deqn{Y_{kij} = X_{1kij}\beta_{1} + Z^{(2)}_{1kij}b^{(2)}_{ki} +
#'   \epsilon_{kij}}
#'
#'   In the above equation, \eqn{Y} represents the longitudinal outcome and
#'   \eqn{X_1} represents the design matrix for the longitudinal fixed effects.
#'   The subscript 1 is used to distinguish between items from the longitudinal
#'   sub-model and items from the survival sub-model (which contain a subscript
#'   2).  The design matrices for random effects are represented using \eqn{Z},
#'   fixed effect coefficients are represented by \eqn{\beta}, random effects by
#'   \eqn{b} and the measurement error by \eqn{\epsilon}.  Study membership is
#'   represented by the subscript \eqn{k} whilst individuals are identified by
#'   \eqn{i} and time points at which they are measured by \eqn{j}.  The
#'   longitudinal outcome is assumed continuous.
#'
#'   Currently this function only supports one linking structure between the
#'   sub-models, namely a random effects only proportional sharing structure. In
#'   this structure, the zero mean random effects from the longitudinal
#'   sub-model are inserted into the survival sub-model, with a common
#'   association parameter for each level of random effects.  Therefore the
#'   survival sub-model (for a case without baseline stratified by study) takes
#'   the following format:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(X_{2ki}\beta_{2} +
#'   \alpha^{(2)}(Z^{(2)}_{1ki}b^{(2)}_{ki}) +
#'   \alpha^{(3)}(Z^{(3)}_{1ki}b^{(3)}_{k})) }
#'
#'   Otherwise, if only individual level random effects are included in the
#'   function call, this reduces to:
#'
#'   \deqn{\lambda_{ki}(t) = \lambda_{0}(t)exp(X_{2ki}\beta_{2} +
#'   \alpha^{(2)}(Z^{(2)}_{1ki}b^{(2)}_{ki}) }
#'
#'   In the above equation, \eqn{\lambda_{ki}(t)} represents the survival time
#'   of the individual \eqn{i} in study \eqn{k}, and \eqn{\lambda_{0}(t)}
#'   represents the baseline hazard.  If a stratified baseline hazard were
#'   specified this would be replaced by \eqn{\lambda_{0k}(t)}.  The design
#'   matrix for the fixed effects in the survival sub-model is represented by
#'   \eqn{X_{2ki}}, with fixed effect coefficients represented by
#'   \eqn{\beta_{2}}.  Association parameters quantifying the link between the
#'   sub-models are represented by \eqn{\alpha} terms.
#'
#'   The model is fitted using an EM algorithm, starting values for which are
#'   extracted from initial separate longitudinal and survival fits.  Pseudo
#'   adaptive Gauss - Hermite quadrature is used to evaluate functions of the
#'   random effects in the EM algorithm, see Rizopoulos 2012.
#'
#'
#' @return An object of class jointmeta1 See \code{\link{jointmeta1.object}}
#'
#' @export
#'
#' @import survival stats
#'
#' @references Wulfsohn, M.S. and A.A. Tsiatis, A Joint Model for Survival and
#'   Longitudinal Data Measured with Error. 1997, International Biometric
#'   Society. p. 330
#'
#'   Rizopoulos, D. (2012) Fast fitting of joint models for longitudinal and
#'   event time data using a pseudo-adaptive Gaussian quadrature rule.
#'   Computational Statistics & Data Analysis 56 (3) p.491-501
#'
#'
#'
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

jointmeta1 <- function(data, long.formula, long.rand.ind, long.rand.stud = NULL,
                       sharingstrct = c("randprop", "randsep", "value", "slope", "valandslope"),
                       surv.formula, gpt, lgpt, max.it, tol, study.name, strat = F, longsep = F,
                       survsep = F, bootrun = F, print.detail = F) {
  if (class(data) != "jointdata") {
    stop("Data should be supplied in jointdata format -
         run tojointdata function if not in jointfdataformat")
  }
  if (sharingstrct != "randprop") {
    stop("Currently only randprop sharing structure supported")
  }
  Call <- match.call()
  id.name <- data$subj.col
  time.long <- data$time.col
  long.formula <- as.formula(long.formula)
  long.formula.orig <- long.formula
  surv.formula <- as.formula(surv.formula)
  if (missing(gpt)) {
    gpt <- 5
  }
  if (missing(lgpt)) {
    lgpt <- 7
  }
  if (missing(max.it)) {
    max.it <- 350
  }
  if (missing(tol)) {
    tol <- 0.001
  }
  if (missing(bootrun)) {
    bootrun <- FALSE
  }
  if (missing(sharingstrct)) {
    stop("No sharing structure specified")
  }
  if ((sharingstrct %in% c("randprop", "randsep", "value", "slope", "valandslope")) ==
      FALSE) {
    stop("Invalid sharing structure specified")
  }
  if (sharingstrct != "randprop") {
    stop("Currently jointmeta only supports randprop sharing structures")
  }
  if (missing(long.rand.ind) == TRUE) {
    stop("Please specify at least one random effect
         at the individual level in long.rand.ind")
  }
  if (length(long.rand.ind) == 0) {
    stop("Please specify at least one random effect
         at the individual level in long.rand.ind")
  }
  if (length(which(("noint" == long.rand.ind) == F)) == 0) {
    stop("Please specify at least one random effect
         at the individual level in long.rand.ind")
  }
  if (("int" %in% long.rand.ind) == TRUE) {
    if (("noint" %in% long.rand.ind) == TRUE) {
      stop("Both the option for no random intercept (noint)
           and random intercept (int) specified in long.rand.ind")
    }
    }
  if (("int" %in% long.rand.ind) == TRUE) {
    long.rand.ind[which((long.rand.ind %in% "int") == TRUE)] <- "(Intercept)"
    if (which(long.rand.ind %in% "(Intercept)") != 1) {
      long.rand.ind <- long.rand.ind[-which(long.rand.ind %in% "(Intercept)")]
      long.rand.ind <- c("(Intercept)", long.rand.ind)
    }
  }
  if (missing(study.name)) {
    stop("Please supply name of study indicator variable to
         \"study.name\" in the function call")
  }
  if (is.null(long.rand.stud) == F) {
    if (study.name %in% long.rand.stud) {
      if (which(long.rand.stud %in% study.name) != 1) {
        long.rand.stud <- long.rand.stud[-which(long.rand.stud %in%
                                                  study.name)]
        long.rand.stud <- c(study.name, long.rand.stud)
      }
    }
  }
  studies <- as.character(unique(data$baseline[[study.name]]))
  numstudies <- length(studies)
  if (any(sapply(data$baseline, "class") == "factor")) {
    data$baseline <- droplevels(data$baseline)
  }
  longdat2 <- merge(data$longitudinal, data$baseline, by = id.name, sort = FALSE)
  long.frame <- model.frame(long.formula, data = longdat2, na.action = na.pass)
  long.cov <- model.matrix(long.formula, long.frame)
  long.terms <- terms(long.formula, data = longdat2)
  long.names <- colnames(long.cov)
  rll <- !is.na(data$longitudinal[[names(long.frame[1])]])
  for (i in 1:length(rll)) {
    if (length(which(is.na(long.cov[i, ]))) > 0) {
      rll[i] <- FALSE
    }
  }
  q <- 0
  for (count in 1:length(long.rand.ind)) {
    if (long.rand.ind[count] != "noint") {
      q <- q + 1
      if (length(which(grepl(long.rand.ind[count], colnames(long.cov)) ==
                       TRUE)) == 0) {
        if (grepl(".", long.rand.ind[count])) {
          temp <- unlist(strsplit(long.rand.ind[count], "."))
          combs <- expand.grid(1:length(temp), 1:length(temp))
          present <- FALSE
          for (i in 1:nrow(combs)) {
            if (!(combs[i, 1] == combs[i, 2])) {
              if (length(which(grepl(paste(temp[combs[i, 1]], temp[combs[i,
                                                                         2]], sep = "."), colnames(long.cov))) == TRUE) >
                  0) {
                present <- TRUE
                long.rand.ind[count] <- paste(temp[combs[i, 1]],
                                              temp[combs[i, 2]], sep = ".")
              }
            }
          }
        }
        if (!present) {
          stop("Individual level random effects included
               in model with no corresponding fixed effect")
        }
        }
      }
    }
  if (q > 3) {
    stop("Model only supports maximum of three individual level random effects")
  }
  if (is.null(long.rand.stud) == FALSE) {
    r <- 0
    for (count in 1:length(long.rand.stud)) {
      if (long.rand.stud[count] != study.name) {
        r <- r + 1
        if (length(which(grepl(long.rand.stud[count], colnames(long.cov)) ==
                         TRUE)) == 0) {
          if (grepl(".", long.rand.stud[count])) {
            temp <- unlist(strsplit(long.rand.stud[count], "."))
            combs <- expand.grid(1:length(temp), 1:length(temp))
            present <- FALSE
            for (i in 1:nrow(combs)) {
              if (!(combs[i, 1] == combs[i, 2])) {
                if (length(which(grepl(paste(temp[combs[i, 1]],
                                             temp[combs[i, 2]], sep = "."), colnames(long.cov))) ==
                           TRUE) > 0) {
                  present <- TRUE
                  long.rand.stud[count] <- paste(temp[combs[i,
                                                            1]], temp[combs[i, 2]], sep = ".")
                }
              }
            }
          }
          if (!present) {
            stop("Study level random effects included
                 in model with no corresponding fixed effect")
          }
          }
        } else {
          r <- r + 1
      }
    }
    if (r > 3) {
      stop("Model only supports maximum of three study level random effects")
    }
    } else {
      r <- NULL
    }
  longdat <- cbind(data$longitudinal[[id.name]][rll], long.frame[, 1][rll],
                   data$longitudinal[[time.long]][rll], longdat2[[study.name]][rll],
                   long.cov[rll, ])
  longdat <- as.data.frame(longdat)
  missingids <- unique(data$longitudinal[[id.name]][!rll])
  names(longdat) <- c(id.name, names(long.frame)[1], time.long, study.name,
                      long.names)
  long.formula <- as.formula(paste(as.character(long.formula)[2], "~",
                                   paste(names(longdat)[5:ncol(longdat)], collapse = " + "), sep = ""))
  p1 <- length(5:ncol(longdat))
  notinteractionterms <- names(longdat[, 5:ncol(longdat)])[!(grepl(":",
                                                                   names(longdat[, 5:ncol(longdat)])))]
  for (count in 1:length(long.rand.ind)) {
    if (length(grep(paste("^", long.rand.ind[count], "$", sep = ""),
                    notinteractionterms)) > 0) {
      long.rand.ind[count] <- notinteractionterms[grep(paste("^",
                                                             long.rand.ind[count], "$", sep = ""), notinteractionterms)]
    } else if (length(grep(paste("^", long.rand.ind[count], "$", sep = ""),
                           notinteractionterms)) == 0) {
      if (long.rand.ind[count] %in% colnames(data$baseline)) {
        if (class(data$baseline[, which(colnames(data$baseline) ==
                                        long.rand.ind[count])]) == "factor") {
          formtemp <- as.formula(paste("~", colnames(data$baseline)[which(colnames(data$baseline) ==
                                                                            long.rand.ind[count])]))
          matrixtemp <- model.matrix(formtemp, data$baseline)
          long.rand.ind[count] <- colnames(matrixtemp)[2:ncol(matrixtemp)]
        }
      } else if (long.rand.ind[count] %in% colnames(data$longitudinal)) {
        if (class(data$longitudinal[, which(colnames(data$longitudinal) ==
                                            long.rand.ind[count])]) == "factor") {
          formtemp <- as.formula(paste("~", colnames(data$longitudinal)[which(colnames(data$longitudinal) ==
                                                                                long.rand.ind[count])]))
          matrixtemp <- model.matrix(formtemp, data$longitudinal)
          long.rand.ind[count] <- colnames(matrixtemp)[2:ncol(matrixtemp)]
        }
      }
    }
  }
  q <- length(long.rand.ind)
  if (q > 3) {
    stop("Model only supports maximum of three individual level random effects")
  }
  if (is.null(long.rand.stud) == FALSE) {
    for (count in 1:length(long.rand.stud)) {
      if (long.rand.stud[count] != study.name) {
        if (length(grep(paste("^", long.rand.stud[count], "$",
                              sep = ""), notinteractionterms)) > 0) {
          long.rand.stud[count] <- notinteractionterms[grep(paste("^",
                                                                  long.rand.stud[count], "$", sep = ""), notinteractionterms)]
        } else if (length(grep(paste("^", long.rand.stud[count],
                                     "$", sep = ""), notinteractionterms)) == 0) {
          if (long.rand.stud[count] %in% colnames(data$baseline)) {
            if (class(data$baseline[, which(colnames(data$baseline) ==
                                            long.rand.stud[count])]) == "factor") {
              formtemp <- as.formula(paste("~", colnames(data$baseline)[which(colnames(data$baseline) ==
                                                                                long.rand.stud[count])]))
              matrixtemp <- model.matrix(formtemp, data$baseline)
              long.rand.stud[count] <- colnames(matrixtemp)[2:ncol(matrixtemp)]
            }
          } else if (long.rand.stud[count] %in% colnames(data$longitudinal)) {
            if (class(data$longitudinal[, which(colnames(data$longitudinal) ==
                                                long.rand.stud[count])]) == "factor") {
              formtemp <- as.formula(paste("~", colnames(data$longitudinal)[which(colnames(data$longitudinal) ==
                                                                                    long.rand.stud[count])]))
              matrixtemp <- model.matrix(formtemp, data$longitudinal)
              long.rand.stud[count] <- colnames(matrixtemp)[2:ncol(matrixtemp)]
            }
          }
        }
      }
    }
    r <- length(long.rand.stud)
    if (r > 3) {
      stop("Model only supports maximum of three study level random effects")
    }
  }
  surv.frame <- model.frame(surv.formula, data = cbind(data$survival,
                                                       data$baseline))
  srv <- model.extract(surv.frame, "response")
  surv.terms <- terms(surv.formula, data = cbind(data$survival, data$baseline))
  attr(surv.terms, "intercept") <- 1
  surv.cov <- model.matrix(surv.terms, data = cbind(data$survival, data$baseline))
  namestemp <- colnames(surv.cov)
  surv.cov <- as.matrix(surv.cov[, -1])
  colnames(surv.cov) <- namestemp[-1]
  rss <- as.integer(row.names(surv.cov))
  survdat <- cbind(data$survival[[id.name]][rss], srv[rss, 1], srv[rss,
                                                                   2], data$baseline[[study.name]][rss], surv.cov)
  survdat <- as.data.frame(survdat)
  names(survdat) <- c(id.name, surv.formula[2][[1]][[2]], surv.formula[2][[1]][[3]],
                      study.name, colnames(surv.cov))
  if (dim(survdat)[2] > 4) {
    survdat[, 5:dim(survdat)[2]] <- scale(survdat[, 5:dim(survdat)[2]],
                                          scale = FALSE)
  }
  survdat2 <- data.frame(data$survival[[id.name]][rss], srv[rss, 1],
                         srv[rss, 2], data$baseline[[study.name]][rss], surv.frame[, -1])
  if (ncol(survdat) > 4) {
    surv.formula <- as.formula(paste(as.character(surv.formula)[2],
                                     "~", paste(names(survdat)[5:ncol(survdat)], collapse = " + "),
                                     sep = ""))
    names(survdat2) <- c(id.name, surv.formula[2][[1]][[2]], surv.formula[2][[1]][[3]],
                         study.name, colnames(surv.frame)[2:ncol(surv.frame)])
  } else {
    surv.formula <- as.formula(paste(as.character(surv.formula)[2],
                                     "~ 1", sep = ""))
    names(survdat2) <- c(id.name, surv.formula[2][[1]][[2]], surv.formula[2][[1]][[3]],
                         study.name, colnames(surv.cov))
  }
  survdat[, 4] <- survdat2[, 4]
  if (ncol(survdat) > 4) {
    p2 <- length(5:ncol(survdat))
  } else {
    p2 <- 0
  }
  rll2 <- rep(TRUE, nrow(survdat2))
  for (i in 1:length(rll2)) {
    if (length(which(is.na(survdat2[i, ]))) > 0) {
      rll2[i] <- FALSE
    }
  }
  if (length(which(rll2 == FALSE)) > 0) {
    missingids <- c(missingids, survdat2[!rll2, 1])
  }
  if (length(missingids) > 0) {
    survdat <- survdat[!(survdat[, 1] %in% missingids), ]
    survdat2 <- survdat2[!(survdat[, 1] %in% missingids), ]
    longdat2 <- longdat2[!(longdat2[, 1] %in% missingids), ]
  }
  sorted <- sortDat(longdat, survdat, longdat2, survdat2)
  longdat <- as.data.frame(sorted$long.s)
  survdat <- as.data.frame(sorted$surv.s)
  longdat2 <- as.data.frame(sorted$long.s2)
  survdat2 <- as.data.frame(sorted$surv.s2)
  if (is.null(long.rand.stud)) {
    ldaests <- longst(longdat = longdat, long.formula.orig = long.formula,
                      long.rand.ind = long.rand.ind, longdat2 = longdat2, id.name = id.name,
                      study.name = study.name, studies = studies)
  } else {
    ldaests <- longst(longdat = longdat, long.formula.orig = long.formula,
                      long.rand.ind = long.rand.ind, long.rand.stud = long.rand.stud,
                      longdat2 = longdat2, id.name = id.name, study.name = study.name,
                      studies = studies)

  }
  if (strat) {
    survests <- survst(survdat = survdat, surv.formula = surv.formula,
                       survdat2 = survdat2, strat = strat, study.name = study.name)
  } else {
    survests <- survst(survdat = survdat, surv.formula = surv.formula,
                       survdat2 = survdat2, strat = strat, study.name = study.name)
  }
  sep.ll <- ldaests$log.like + survests$log.like[2]
  sep.loglik <- list(seplhood = sep.ll, sepy = ldaests$log.like, sepn = survests$log.like[2])
  paraests <- c(ldaests, survests)
  if (sharingstrct == "randprop") {
    if (bootrun == FALSE) {
      message("Running EM algorithm...")
    }
    jointfit <- EMalgRandprop(data = data, longdat = longdat, survdat = survdat,
                              long.rand.ind = long.rand.ind, long.rand.stud = long.rand.stud,
                              id.name = id.name, study.name = study.name, gpt = gpt, max.it = max.it,
                              tol = tol, time.long = time.long, surv.formula = surv.formula,
                              long.formula = long.formula, long.formula.orig = long.formula.orig,
                              paraests = paraests, studies = studies, p1 = p1, p2 = p2, strat = strat,
                              print.detail = print.detail, bootrun = bootrun, q = q, r = r)
    likeests <- c(jointfit, list(rs = survests$rs, sf = survests$sf))
    beta1 <- jointfit$beta1
    rownames(beta1) <- rownames(paraests$beta1)
    if (p2 > 0) {
      beta2 <- jointfit$beta2[1:p2, ]
      names(beta2) <- names(paraests$beta2)
    } else {
      beta2 <- NULL
    }
    fixed <- list(longitudinal = beta1, survival = beta2)
    D <- jointfit$D
    random_ind <- jointfit$random2
    ids.bystudy <- lapply(1:numstudies, function(u) {
      survdat[which(survdat[, 4] == studies[u]), 1]
    })
    random_ind <- lapply(1:numstudies, function(u) {
      randtemp <- random_ind[[u]]
      colnames(randtemp) <- paste("b2_", 0:(ncol(randtemp) - 1),
                                  sep = "")
      rownames(randtemp) <- ids.bystudy[[u]]
      randtemp
    })
    random <- list(random_ind = random_ind)
    if ("(Intercept)" %in% long.rand.ind) {
      long.rand.ind2 <- long.rand.ind
      long.rand.ind2[which(long.rand.ind2 == "(Intercept)")] <- "1"
      long.rand.ind.form <- paste(long.rand.ind2, collapse = " + ")
    }
    if ("noint" %in% long.rand.ind) {
      long.rand.ind2 <- long.rand.ind[-which(long.rand.ind == "noint")]
      long.rand.ind.form <- paste("-1", long.rand.ind2, sep = " + ")
    }
    n.bystudy <- jointfit$n.bystudy
    if (is.null(long.rand.stud) == FALSE) {
      A <- jointfit$A
      latent <- jointfit$beta2[(p2 + 1):(p2 + 2), ]
      names(latent) <- c(paste("gamma_ind_", 0, sep = ""), paste("gamma_stud_",
                                                                 0, sep = ""))
      random_stud <- jointfit$random3
      colnames(random_stud) <- paste("b3_", 0:(ncol(random_stud) -
                                                 1), sep = "")
      rownames(random_stud) <- studies
      random$random_stud <- random_stud
      randstart.stud.l <- paraests$randstart.stud
      randstart.stud.cov.l <- paraests$randstart.stud.cov
      if (study.name %in% long.rand.stud) {
        long.rand.stud2 <- long.rand.stud
        long.rand.stud2[which(long.rand.stud2 == study.name)] <- "1"
        long.rand.stud.form <- paste(long.rand.stud2, collapse = " + ")
      } else {
        long.rand.stud.form <- paste("-1", paste(long.rand.stud,
                                                 collapse = " + "), sep = " + ")
      }
    } else {
      latent <- jointfit$beta2[(p2 + 1), ]
      names(latent) <- paste("gamma_ind_", 0, sep = "")
      randstart.stud.l <- NULL
      randstart.stud.cov.l <- NULL
    }
    coefficients <- list(fixed = fixed, random = random, latent = latent)
    if (bootrun == FALSE) {
      message("Calculating log-likelihood...")
      jointll <- jlike(data = data, longdat = longdat, survdat = survdat,
                       q = q, likeests = likeests, lgpt = lgpt, studies = studies,
                       p1 = p1, p2 = p2, long.rand.ind = long.rand.ind, randstart.ind = paraests$randstart.ind,
                       randstart.ind.cov = paraests$randstart.ind.cov, r = r,
                       long.rand.stud = long.rand.stud, randstart.stud = randstart.stud.l,
                       randstart.stud.cov = randstart.stud.cov.l, strat = strat,
                       study.name = study.name, id.name = id.name)
      numpara <- p1 + p2 + (q^2) + 2
      if (!is.null(long.rand.stud)) {
        numpara <- numpara + (r^2) + 1
      }
      AIC <- (2 * numpara) - (2 * jointll$log.like)
      loglik <- list(jointlhood = jointll$log.like, jointy = jointll$longlog.like,
                     jointn = jointll$survlog.like)

    } else {
      loglik <- "Not Calculated"
      AIC <- "Not Calculated"
    }
    sepests <- list(longests = sep(ldaests, longsep), survests = sep(survests,
                                                                     survsep))
    formulae <- list(lformula = long.formula, sformula = surv.formula,
                     rand_ind_formula = as.formula(paste("~", long.rand.ind.form,
                                                         sep = "")))

    rand_cov <- list(D = jointfit$D)
    if (is.null(long.rand.stud) == FALSE) {
      formulae$rand_stud_formula <- as.formula(paste("~", long.rand.stud.form,
                                                     sep = ""))
      rand_cov$A <- jointfit$A
    }
    nobs <- table(longdat[[study.name]])
    names(nobs) <- studies
    results <- list(coefficients = coefficients, sigma.e = jointfit$sigma.e,
                    rand_cov = rand_cov, hazard = jointfit$haz, loglik = loglik,
                    numIter = jointfit$iters, convergence = jointfit$conv, sharingstrct = sharingstrct,
                    sepests = sepests, sep.loglik = sep.loglik, data = data, Call = Call,
                    numstudies = numstudies, n.bystudy = n.bystudy,
                    missingids = missingids, nobs = nobs, AIC = AIC)
    class(results) <- "jointmeta1"
    results
  }
}
