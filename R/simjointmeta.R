#' Simulation of multi-study joint data
#'
#' Function to allow the simulation of a correlated single continuous
#'     longitudinal outcome and a single survival outcome for data from multiple
#'     studies.  The longitudinal sub-model contains a fixed intercept, time
#'     (slope) term and a binary treatment assignment covariate, whilst the
#'     survival sub-model contains only a binary treatment assignment covariate.
#'
#' @param k the number of studies to be simulated
#' @param n a vector of length equal to k denoting the number of individuals to
#'     simulate per study
#' @param sepassoc a logical taking value \code{FALSE} if proportional
#'     association is required, \code{TRUE} if a separate association parameter
#'     is required for each random effect shared between the sub-models
#' @param ntms the maximum possible number of longitudinal measurements - should
#'     equal the length of the supplied \code{longmeasuretimes}
#' @param longmeasuretimes a vector giving the exact times of the longitudinal
#'     measurement times.  If this is not specified in the function call then
#'     the measurement times of the longitudinal outcome are set to start at 0
#'     then take integer values up to and including \code{ntms - 1}.
#' @param beta1 a vector of the fixed effects for the longitudinal sub-model.
#'     Here the first element gives the coefficient for a fixed or population
#'     intercept, the second gives the coefficient for the binary treatment
#'     assignment covariate and the third element gives the covariate for the
#'     time (slope) covariate
#' @param beta2 the coefficient for the binary treatment assignment covariate
#' @param rand_ind a character string specifying the individual level random
#'     effects structure.  If \code{rand_ind = "intslope"} then there is an
#'     individual specific random intercept and random time (slope) term
#'     included in the model.  If \code{rand_ind = "int"} then the model
#'     includes only a individual specific random intercept.
#' @param rand_stud a character string specifying the study level random effects
#'     structure.  If this is set to \code{NULL} or not specified in the
#'     function call then no study level random effects are included in the
#'     model that the data is simulated from.  There are three options if data
#'     is to be simulated with random effects at the study level.  If a study
#'     level random intercept only is to be included, then set
#'     \code{rand_stud = "int"}.  Else if a study level random treatment
#'     assignment term only is to be included then set
#'     \code{rand_stud = "treat"}.  Finally if both a study level random
#'     intercept and a study level random treatment effect is to be included,
#'     then set \code{rand_stud = "inttreat"}.
#' @param gamma_ind parameter specifying the level of association between the
#'     longitudinal and survival outcomes attributable to the individual
#'     deviation from the population longitudinal trajectory. If different
#'     association parameters are required for each
#'     study then a list of length equal to the number of studies should be
#'     supplied to \code{gamma_ind}.  If \code{sepassoc = TRUE} then
#'     \code{gamma_ind} should be either a vector of values of length equal to
#'     the number of individual level random effects, or a list of vectors each
#'     of length equal to the number of individual level random effects.
#'     However if \code{sepassoc = FALSE} then \code{gamma_ind} should be
#'     supplied as a single value, or a list of single values.
#' @param gamma_stud parameter specifying the level of association between the
#'     longitudinal and survival outcomes attributable to the study level
#'     deviation from the overall population longitudinal trajectory. If
#'     different association parameters are required for each
#'     study then a list of length equal to the number of studies should be
#'     supplied to \code{gamma_stud}.  If \code{sepassoc = TRUE} then
#'     \code{gamma_stud} should be either a vector of values of length equal to
#'     the number of study level random effects, or a list of vectors each
#'     of length equal to the number of study level random effects.
#'     However if \code{sepassoc = FALSE} then \code{gamma_stud} should be
#'     supplied as a single value, or a list of single values.  This parameter
#'     should only be present if \code{rand_stud} is specified in the function
#'     call.
#' @param sigb_ind the covariance matrix for the individual level random
#'     effects. This should have number of rows and columns equal to the
#'     number of individual level random effects.
#' @param sigb_stud the covariance matrix for the study level random
#'     effects. This should have number of rows and columns equal to the
#'     number of study level random effects.  This should only be specified if
#'     \code{rand_stud} is specified in the function call.
#' @param vare the variance of the measurement error term
#' @param theta0 parameter defining the distribution of the survival times.
#'     A separate parameter can be defined per study or a common parameter
#'     across all studies.  See Bender et al 2005 for advice on approximating
#'     appropriate values for \code{theta0} and \code{theta1} the using extreme
#'     value distribution.
#' @param theta1 parameter defining the distribution of the survival times.
#'     A separate parameter can be defined per study or a common parameter
#'     across all studies.  See Bender et al 2005 for advice on approximating
#'     appropriate values for \code{theta0} and \code{theta1} the using extreme
#'     value distribution.
#' @param censoring a logical indicating whether the simulated survival times
#'     should be censored or not
#' @param censlam the lambda parameter controlling the simulated exponentially
#'     distributed censoring times.  This can either be supplied as one value
#'     for all studies simulated, or a vector of length equal to the number of
#'     studies in the dataset.
#' @param truncation a logical value to specify whether the simulated survival
#'     times should be truncated at a specified time or not.
#' @param trunctime if \code{truncation = TRUE} then the survival times will be
#'     truncated at the specified \code{trunctime}
#'
#' @return This function returns a list with three named elements.  The first
#'     element is named \code{"longdat"}, the second \code{"survdat"}, the third
#'     \code{"percentevent"}.  Each of these elements is a list of length equal
#'     to the number of studies specified to simulate in the function call.
#'
#'     The element \code{"longdat"} is a list of the simulated longitudinal data
#'     sets.  Each longitudinal dataset contains the following variables:
#'     \describe{
#'
#'     \item{\code{id}}{a numeric id variable}
#'
#'     \item{\code{Y}}{the continuous longitudinal outcome}
#'
#'     \item{\code{time}}{the numeric longitudinal time variable}
#'
#'     \item{\code{study}}{a study membership variable}
#'
#'     \item{\code{intercept}}{an intercept term}
#'
#'     \item{\code{treat}}{a treatment assignment variable to one of two
#'         treatment groups}
#'
#'     \item{\code{ltime}}{a duplicate of the longitudinal time variable}
#'
#'     }
#'
#'     The element \code{"survdat"} is a list of the simulated survival data
#'     sets.  Each survival dataset contains the following variables:
#'     \describe{
#'
#'     \item{\code{id}}{a numeric id variable}
#'
#'     \item{\code{survtime}}{the numeric survival times}
#'
#'     \item{\code{cens}}{the censoring indicator}
#'
#'     \item{\code{study}}{a study membership variable}
#'
#'     \item{\code{treat}}{a treatment assignment variable to one of two
#'         treatment groups}
#'
#'     }
#'
#'     The element \code{"percentevent"} is a list of the percentage of events
#'     over censorings seen in the simulated survival data.
#'
#' @details  This function allows the simulation of a single continuous
#'     longitudinal and a single survival outcome which are potentially
#'     correlated.  The model simulates data under a joint model with a zero
#'     mean random effects only sharing structure.  The longitudinal sub-model
#'     is adjusted by a fixed or population intercept, time (slope) term and a
#'     binary treatment assignment covariate.  The survival sub-model is
#'     adjusted by only the fixed or population binary treatment assignment
#'     covariate.
#'
#'     Random effects can be specified at either just the individual level,
#'     or at both the individual and study level.  For the options for the
#'     random effects see the above parameter definitions.
#'
#'     The parameters controlling the distributions for the survival times and
#'     the censoring times can be identical across the studies, or separate
#'     values can be supplied for each study.  Similarly the association
#'     parameters can be identical across studies, or unique to each study.
#'
#'     The simulated longitudinal information is capped at each individual's
#'     survival time.  If \code{truncation= TRUE} then the survival times are
#'     truncated at the specified \code{trunctime}.
#'
#'     For description of the methodology of simulating this data see Bender et
#'     al 2005, and Austin 2012.
#'
#'     Note that this function does not return data in a \code{jointdata}
#'     format.  Function \code{\link{tojointdata}} can help to reformat this
#'     data into a \code{jointdata} format.
#'
#' @seealso \code{\link{tojointdata}}
#'
#' @export
#' @import MASS
#'
#' @references Bender et al (2005) Generating survival times to simulate Cox
#'     proportional hazards models. Statistics in Medicine 24:1713–1723
#'
#'     Austin (2012) Generating survival times to simulate Cox proportional
#'     hazards models with time-varying covariates. Statistics in Medicine
#'     31: 3946–3958
#'
#' @examples
#'  #simulated data without study level variation specified
#'  exampledat1<-simjointmeta(k = 5, n = rep(500, 5), sepassoc = FALSE,
#'               ntms = 5, longmeasuretimes = c(0, 1, 2, 3, 4),
#'               beta1 = c(1, 2, 3), beta2 = 1, rand_ind = "intslope",
#'               rand_stud = NULL, gamma_ind = 1,
#'               sigb_ind = matrix(c(1,0.5,0.5,1.5),nrow=2), vare = 0.01,
#'               theta0 = -3, theta1 = 1, censoring = TRUE, censlam = exp(-3),
#'               truncation = FALSE, trunctime = max(longmeasuretimes))
#'
#'  #simulated data with different parameters for each study for the
#'  #association parameters, censoring distribution parameters and survival time
#'  #parameters
#'  gamma_ind_set<-list(c(0.5, 1), c(0.4, 0.9), c(0.6, 1.1), c(0.5, 0.9),
#'                      c(0.4, 1.1))
#'  gamma_stud_set<-list(c(0.6, 1.1), c(0.5, 1), c(0.5, 0.9), c(0.4, 1.1),
#'                      c(0.4, 0.9))
#'  censlamset<-c(exp(-3), exp(-2.9), exp(-3.1), exp(-3), exp(-3.05))
#'  theta0set<-c(-3, -2.9, -3, -2.9, -3.1)
#'  theta1set<-c(1, 0.9, 1.1, 1, 0.9)
#'
#'  exampledat2<-simjointmeta(k = 5, n = rep(500, 5), sepassoc = TRUE, ntms = 5,
#'                            longmeasuretimes = c(0, 1, 2, 3, 4),
#'                            beta1 = c(1, 2, 3), beta2 = 1,
#'                            rand_ind = "intslope", rand_stud = "inttreat",
#'                            gamma_ind = gamma_ind_set,
#'                            gamma_stud = gamma_stud_set,
#'                            sigb_ind = matrix(c(1, 0.5, 0.5, 1.5), nrow = 2),
#'                            sigb_stud = matrix(c(1, 0.5, 0.5, 1.5), nrow = 2),
#'                            vare = 0.01, theta0 = theta0set,
#'                            theta1 = theta1set, censoring = TRUE,
#'                            censlam = censlamset, truncation = FALSE,
#'                            trunctime = max(longmeasuretimes))
#'
#'
simjointmeta <- function(k = 5, n = rep(500, 5), sepassoc = FALSE, ntms = 5,
                         longmeasuretimes = c(0, 1, 2, 3, 4),
                         beta1 = c(1, 1, 1), beta2 = 1,
                         rand_ind = c("intslope", "int"),
                         rand_stud = c("int", "inttreat", "treat", NULL),
                         gamma_ind = 1, gamma_stud = NULL, sigb_ind,
                         sigb_stud = NULL, vare = 0.01, theta0 = -3, theta1 = 1,
                         censoring = TRUE, censlam = exp(-3),
                         truncation = FALSE,
                         trunctime = max(longmeasuretimes)) {
  if(!is.null(gamma_stud) && !is.null(rand_stud) && !is.null(sigb_stud)) {
    if(length(n) != k) {
      stop("Number of studies differs between k and length of n")
    }
    rand_ind <- match.arg(rand_ind)
    if(rand_ind != "intslope" && rand_ind != "int") {
      stop(paste("Unknown individual level random effects specification:",
                 rand_ind))
    }
    rand_stud <- match.arg(rand_stud)
    if(rand_stud != "int" && rand_stud != "inttreat" && rand_stud != "treat") {
      stop(paste("Unknown study level random effects specification:",
                 rand_stud))
    }
    if(missing(sigb_ind)) {
      stop("Missing argument: sigb_ind")
    }
    if(missing(sigb_stud)) {
      stop("Missing argument: sigb_stud")
    }
    if(missing(ntms)) {
      stop("Missing argument: ntms")
    }
    samegamma <- TRUE
    if(class(gamma_ind) == "list" && class(gamma_stud) == "list") {
      samegamma <- FALSE
    }else if((class(gamma_ind) == "list" && class(gamma_stud) != "list") ||
             (class(gamma_ind) != "list" && class(gamma_stud) == "list")) {
      stop("one but not both of gamma_ind and gamma_stud supplied as
           varying between study - specify as both varying or both constant")
    }
    if(!(length(theta0) %in% c(1, k))) {
      stop("Supply either one theta0 or one per study")
    }
    if(length(theta0) == 1) {
      theta0 <- rep(theta0, k)
    }
    if(!(length(theta1) %in% c(1, k))) {
      stop("Supply either one theta1 or one per study")
    }
    if(length(theta1) == 1) {
      theta1 <- rep(theta1, k)
    }
    if(!(length(censlam) %in% c(1, k))) {
      stop("Supply either one censlam or one per study")
    }
    if(length(censlam) == 1) {
      censlam <- rep(censlam, k)
    }
    if(rand_ind == "intslope") {
      q <- 2
    }else if(rand_ind == "int") {
      q <- 1
    }
    if(rand_stud == "inttreat") {
      r <- 2
    }else{
      r <- 1
    }
    lat = q+r
    if(!sepassoc) {
      lat = 2
      if(samegamma) {
        if((length(gamma_ind)+length(gamma_stud)) != lat) {
          cat("Warning: Number of association parameters
              do not match model choice\n")
        }
        gamma = c(rep(gamma_ind, q), rep(gamma_stud, r))
      }else{
        for(i in 1:k) {
          if((length(gamma_ind[[i]])+length(gamma_stud[[i]])) != lat) {
            cat("Warning: Number of association parameters
                do not match model choice\n")
          }
        }
        gamma = lapply(1:k, function(u) {c(rep(gamma_ind[[u]], q),
                                           rep(gamma_stud[[u]], r))})
      }
    }else{
      if(samegamma) {
        if((length(gamma_ind)+length(gamma_stud)) != lat) {
          cat("Warning: Number of association parameters
              do not match model choice\n")
        }
        gamma = c(gamma_ind, gamma_stud)
      }else{
        for(i in 1:k) {
          if((length(gamma_ind[[i]])+length(gamma_stud[[i]])) != lat) {
            cat("Warning: Number of association parameters
                do not match model choice\n")
          }
        }
        gamma = lapply(1:k, function(u) {c(gamma_ind[[u]], gamma_stud[[u]])})
      }
    }
    if(length(sigb_ind) != q^2) {
      cat("Warning: Dimension of individual level covariance matrix
          does not match chosen rand_ind\n")
      if(length(sigb_ind)>q^2) {
        sigb_ind = sigb_ind[1:q, 1:q]
      }
      else{
        sigb_ind = diag(q)*sigb_ind[1]
      }
    }
    if(length(sigb_stud) != r^2) {
      cat("Warning: Dimension of individual level covariance matrix
          does not match chosen rand_ind\n")
      if(length(sigb_stud)>r^2) {
        sigb_stud = sigb_stud[1:r, 1:r]
      }
      else{
        sigb_stud = diag(r)*sigb_stud[1]
      }
    }
    if(q == 1) {
      if(sigb_ind<0) {
        stop("Variance must be positive")
      }
    }else{
      if(!isSymmetric(sigb_ind)) {
        stop("Individual level Covariance matrix is not symmetric")
      }
      if(any(eigen(sigb_ind)$values<0) || (det(sigb_ind)<=0)) {
        stop("Individual level Covariance matrix must be
             positive semi-definite")
      }
    }
    if(r == 1) {
      if(sigb_stud<0) {
        stop("Variance must be positive")
      }
    }else{
      if(!isSymmetric(sigb_stud)) {
        stop("Study level Covariance matrix is not symmetric")
      }
      if(any(eigen(sigb_stud)$values<0) || (det(sigb_stud)<=0)) {
        stop("Study level Covariance matrix must be positive semi-definite")
      }
    }
    if(missing(longmeasuretimes)) {
      longmeasuretimes <- 0:(ntms-1)
    }
    sim <- simdat2randlevels(k = k, n = n, rand_ind = rand_ind,
                             rand_stud = rand_stud, sepassoc = sepassoc,
                             ntms = ntms, longmeasuretimes = longmeasuretimes,
                             beta1 = beta1, beta2 = beta2, gamma = gamma,
                             sigb_ind = sigb_ind, sigb_stud = sigb_stud,
                             vare = vare, theta0 = theta0, theta1 = theta1,
                             censoring = censoring, censlam = censlam,
                             truncation = truncation, trunctime = trunctime,
                             q = q, r = r)
    list(longitudinal = sim$longdat, survival = sim$survdat,
         percentevent = sim$percentevent)
  }else{
    if(!is.null(gamma_stud) || !is.null(sigb_stud) || !is.null(rand_stud)) {
      stop("Some but not all of gamma_stud, rand_stud and sigb_stud supplied")
    }else{
      if(length(n) != k) {
        stop("Number of studies differes between k and length of n")
      }
      rand_ind <-  match.arg(rand_ind)
      if(rand_ind != "intslope" && rand_ind != "int") {
        stop(paste("Unknown individual level random effects specification:",
                   rand_ind))
      }
      if(missing(sigb_ind)) {
        stop("Missing argument: sigb_ind")
      }
      samegamma <- TRUE
      if(class(gamma_ind) == "list") {
        samegamma <- FALSE
      }
      if(!(length(theta0) %in% c(1, k))) {
        stop("Supply either one theta0 or one per study")
      }
      if(length(theta0) == 1) {
        theta0 <- rep(theta0, k)
      }
      if(!(length(theta1) %in% c(1, k))) {
        stop("Supply either one theta1 or one per study")
      }
      if(length(theta1) == 1) {
        theta1 <- rep(theta1, k)
      }
      if(!(length(censlam) %in% c(1, k))) {
        stop("Supply either one censlam or one per study")
      }
      if(length(censlam) == 1) {
        censlam <- rep(censlam, k)
      }
      if(rand_ind == "intslope") {
        q <- 2
      }else if(rand_ind == "int") {
        q <- 1
      }
      if(missing(ntms)) {
        stop("Missing argument: ntms")
      }
      lat = q
      if(!sepassoc) {
        lat = 1
        if(samegamma) {
          if(length(gamma_ind) != lat) {
            cat("Warning: Number of association parameters
                do not match model choice\n")
          }
          gamma = rep(gamma_ind, q)
        }else{
          for(i in 1:k) {
            if(length(gamma_ind[[i]]) != lat) {
              cat("Warning: Number of association parameters
                  do not match model choice\n")
            }
          }
          gamma = lapply(1:k, function(u) {rep(gamma_ind[[u]], q)})
        }
      }else{
        if(samegamma) {
          if(length(gamma_ind) != lat) {
            cat("Warning: Number of association parameters
                do not match model choice\n")
          }
          gamma = gamma_ind
        }else{
          for(i in 1:k) {
            if(length(gamma_ind[[i]]) != lat) {
              cat("Warning: Number of association parameters
                  do not match model choice\n")
            }
          }
          gamma = lapply(1:k, function(u) {gamma_ind[[u]]})
        }
      }

      if(length(sigb_ind) != q^2) {
        cat("Warning: Dimension of individual level covariance matrix
            does not match chosen rand_ind\n")
        if(length(sigb_ind)>q^2) {
          sigb_ind = sigb_ind[1:q, 1:q]
        }
        else{
          sigb_ind = diag(q)*sigb_ind[1]
        }
      }
      if(q == 1) {
        if(sigb_ind<0) {
          stop("Variance must be positive")
        }
      }else{
        if(!isSymmetric(sigb_ind)) {
          stop("Individual level Covariance matrix is not symmetric")
        }
        if(any(eigen(sigb_ind)$values<0) || (det(sigb_ind)<=0)) {
          stop("Individual level Covariance matrix must
               be positive semi-definite")
        }
      }
      if(missing(longmeasuretimes)) {
        longmeasuretimes <- 0:(ntms-1)
      }
      sim <- simdat1randlevel(k = k, n = n, rand_ind = rand_ind,
                              sepassoc = sepassoc, ntms = ntms,
                              longmeasuretimes = longmeasuretimes,
                              beta1 = beta1, beta2 = beta2, gamma = gamma,
                              sigb_ind = sigb_ind, vare = vare,
                              theta0 = theta0, theta1 = theta1,
                              censoring = censoring, censlam = censlam,
                              truncation = truncation, trunctime = trunctime,
                              q = q)
      list(longitudinal = sim$longdat, survival = sim$survdat,
           percentevent = sim$percentevent)
    }
  }
}
