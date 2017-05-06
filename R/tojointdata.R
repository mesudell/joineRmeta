#' Function to change multi-study data into jointdata format
#'
#' This function is designed to take data in various formats from multiple
#' studies and output the data in a \code{jointdata} format.
#'
#' @param dataset a dataset or list of datasets.  The datasets can be of either
#'   long or wide format holding both the longitudinal and survival information
#'   (and any baseline information).  Either the parameter \code{dataset} should
#'   be specified, or the parameters \code{longitudinal} and \code{survival}
#'   (and \code{baseline} if available) should be supplied.
#' @param longitudinal a dataset or list of datasets in long format containing
#'   the longitudinal outcome and any time varying covariates.  It can also
#'   contain baseline information.  This should not be supplied if
#'   \code{dataset} is supplied to the function call.
#' @param survival a dataset or list of datasets in wide format (one row per
#'   individual) containing the survival information (survival time and
#'   censoring variable).  It can also contain baseline information. This should
#'   not be supplied if \code{dataset} is supplied to the function call.
#' @param baseline a dataset or list of datasets in wide format (one row per
#'   individual) containing any baseline information.  This variable does not
#'   have to be supplied if there is no baseline information, or if it is
#'   already contained in the longitudinal or survival datasets.  This should
#'   not be supplied if \code{dataset} is supplied to the function call.
#' @param id a character string holding the name of the id variable.  This
#'   should be present in all datasets supplied to this function.
#' @param longoutcome a character string holding the name of the longitudinal
#'   outcome.
#' @param timevarying a vector of character strings indicating the names of the
#'   time varying covariates in the dataset
#' @param survtime a character string denoting the name of the surivival time
#'   variable in the dataset
#' @param cens a character string denoting the name of the censoring variable in
#'   the dataset
#' @param time a character string to label the time variable if the data is
#'   transformed from wide to long format, or the name of the variable holding
#'   the time variable if the data is supplied in long format
#' @param longtimes if wide data, the labels that denote the time varying
#'   variables to allow the longitudinal data to be returned in long data
#'   format.
#'
#' @return A jointdata object, see \code{\link[joineR]{jointdata}}.
#'
#' @details The data supplied to the \code{jointmeta1} function has to be in a
#'   \code{\link[joineR]{jointdata}} format.  However it is conceivable that
#'   data supplied from multiple studies could come in a range of formats.
#'
#'   This function can handle a range of formats in which data from multiple
#'   studies may be supplied.  This are discussed below.  We refer to wide data
#'   as data that contains both time varying and time stationary data, with one
#'   line per individual with variables measured over time recorded in multiple
#'   columns.  We refer to long data as data with multiple lines per individual,
#'   with time varying covariates differing between rows, but time stationary
#'   covariates identical between rows.  Survival data are considered data with
#'   survival outcome (survival time and censoring variable) with or without
#'   baseline data.  Longitiudinal datasets are considered long format datasets
#'   containing time varying data potentially with baseline data. Baseline
#'   datasets are considered wide format datasets containing non-time varying
#'   data measured at baseline. This function can take data in the following
#'   formats and output a \code{jointdata} object. \describe{
#'
#'   \item{\code{One wide dataset}}{One dataset in wide format (one row per
#'   individual) supplied to the parameter \code{dataset} in the function call
#'   with \code{longitudinal}, \code{survival} and \code{baseline} set to
#'   \code{NULL}, or left unspecified in the function call.  This dataset would
#'   contain all the data from all studies, and the survival, longitudinal and
#'   any baseline information would all be present in the same dataset.}
#'
#'   \item{\code{One long dataset}}{One dataset in long format (multiple rows
#'   for each individual) supplied to the parameter \code{dataset} in the
#'   function call with \code{longitudinal}, \code{survival} and \code{baseline}
#'   set to \code{NULL}, or left unspecified in the function call.  This dataset
#'   would contain all the data from all studies, and the survival, longitudinal
#'   and any baseline information would all be present in the same dataset.}
#'
#'   \item{\code{A list of study specific wide datasets}}{One dataset for each
#'   study, each in wide format (one row per individual), supplied to the
#'   parameter \code{dataset} in the function call with \code{longitudinal},
#'   \code{survival} and \code{baseline} set to \code{NULL}, or left unspecified
#'   in the function call.  The data from each study would contain the survival,
#'   longitudinal and any baseline information.}
#'
#'   \item{\code{A list of study specific long datasets}}{One dataset for each
#'   study, each in long format (multiple rows per individual), supplied to the
#'   parameter \code{dataset} in the function call with \code{longitudinal},
#'   \code{survival} and \code{baseline} set to \code{NULL} or left unspecified
#'   in the function call.  The data from each study would contain the survival,
#'   longitudinal and any baseline information.}
#'
#'   \item{\code{One longitudinal and one survival dataset with or without an
#'   additional baseline dataset}}{In this case all the longitudinal and time
#'   varying data for all studies is supplied in a single dataset in long format
#'   to the parameter \code{longitudinal}.  All the survival data for all
#'   studies is supplied in a single dataset in wide format to the parameter
#'   \code{survival}.  Baseline data can be present in these two datasets, or
#'   can also be supplied as a dataset to the parameter \code{baseline}. If
#'   \code{longitudinal} and \code{survival} are specified, then parameter
#'   \code{dataset} should be set to \code{NULL} or left unspecified in the
#'   function call. The parameter \code{baseline} is optional, but should only
#'   be specified if parameter \code{dataset} is \code{NULL} or unspecified.}
#'
#'   \item{\code{A list of longitudinal and a list of survival datasets with or
#'   without a list of baseline datasets}}{In this case the longitudinal and
#'   time varying data for each study is supplied as one element of a list of
#'   long format datasets to the parameter \code{longitudinal}.  The survival
#'   data for each study is supplied as one element of a list of wide format
#'   datasets to the parameter \code{survival}. Baseline data can be present in
#'   these two sets of datasets, or can be supplied as an additional list of
#'   datasets one for each study to the parameter \code{baseline}. If
#'   \code{longitudinal} and \code{survival} are specified (\code{baseline} is
#'   optional), then parameter \code{dataset} should be set to \code{NULL}, or
#'   left unspecified in the function call.}
#'
#'   }
#'
#'   The specified id variable should be present in all datasets supplied to the
#'   function.  Variables containing the same information should be identically
#'   named in each supplied dataset, for example if a variable \code{'age'} is
#'   present in one dataset, denoting age of individual at baseline,
#'   corresponding variables in other datasets also supplying age at baseline
#'   should also be named \code{'age'}.  Similarly, different variables should
#'   not share the same name across different datasets, for example there should
#'   not be a variable named \code{'age'} in the longitudinal dataset denoting
#'   individual's age at last longitudinal measurement along with a variable
#'   \code{'age'} in the baseline dataset that denotes age of the individual at
#'   baseline.  Before supplying data to this function, names of variables in
#'   each dataset should be checked to confirm that common variables share the
#'   same name, and differing variables are appropriately distinguished from
#'   each other.
#'
#' @export
#' @import joineR
#'
#' @seealso \code{\link[joineR]{jointdata}}, \code{\link{jointmeta1}}
#'
#' @examples
#'
#'    #simdat is a simulated dataset available in the joineRmeta package
#'    #it is supplied as a list of longitudinal and a list of survival datasets,
#'    #each list is of length equal to the number of studies in the entire
#'    #dataset.
#'    jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                          survival = simdat$survival, id = 'id',
#'                          longoutcome = 'Y', timevarying = c('time','ltime'),
#'                          survtime = 'survtime', cens = 'cens', time = 'time')
#'
tojointdata <- function(dataset = NULL, longitudinal = NULL, survival = NULL,
                        baseline = NULL, id, longoutcome, timevarying = NULL, survtime, cens,
                        time = NULL, longtimes = NULL) {
  if (missing(id)) {
    stop("Name of id variable has not been supplied")
  }
  if (length(id) > 1) {
    stop("More than one name supplied for the id variable")
  }
  if (missing(longoutcome)) {
    stop("At least one longitudinal variable name must
         be supplied to longoutcome")
  }
  if (missing(survtime)) {
    stop("No variable name supplied for survtime")
  }
  if (missing(cens)) {
    stop("No variable name supplied for cens")
  }
  if (is.null(dataset) == FALSE) {
    if (!is.null(longitudinal) || !is.null(survival) || !is.null(baseline)) {
      stop("More than one of dataset, and combination of \n
           longitudinal, survival and baseline datasets supplied")
    }
    if (class(dataset) == "list") {
      dataset <- Reduce(function(x, y) merge(x, y, all = TRUE), dataset)
    }
    idcol <- which(names(dataset) %in% id)
    datatype <- "long"
    if (nrow(dataset) == length(unique(dataset[, idcol]))) {
      datatype <- "wide"
    }
    if (datatype == "wide") {
      survcols <- which(names(dataset) %in% c(id, survtime, cens))
      survivaltemp <- dataset[, survcols]
      survivaltemp <- survivaltemp[, c(which(names(survivaltemp) %in%
                                               id), which(names(survivaltemp) %in% survtime), which(names(survivaltemp) %in%
                                                                                                      cens))]
      if (is.null(timevarying) == FALSE) {
        timevariable <- lapply(c(longoutcome, timevarying), function(u) {
          names(dataset)[grepl(u, names(dataset))]
        })
      } else {
        timevariable <- lapply(longoutcome, function(u) {
          names(dataset)[grepl(u, names(dataset))]
        })
      }
      if (is.null(longtimes)) {
        longtimes <- 1:length(timevariable[[1]])
      }
      longtemp <- reshape(dataset, varying = timevariable, v.names = c(longoutcome,
                                                                       timevarying), times = longtimes, direction = "long")
      if (is.null(timevarying)) {
        longtemp <- longtemp[, which(names(longtemp) %in% c(id,
                                                            longoutcome, "time"))]
      } else {
        longtemp <- longtemp[, which(names(longtemp) %in% c(id,
                                                            longoutcome, "time", timevarying))]
      }
      if (is.null(time)) {
        time <- "time"
      }
      if (length(which(!(names(longtemp) %in% c(id, longoutcome,
                                                time)))) > 0) {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id),
                                 which(names(longtemp) %in% longoutcome), which(names(longtemp) %in%
                                                                                  time), which(!(names(longtemp) %in% c(id, longoutcome,
                                                                                                                        time))))]
      } else {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id),
                                 which(names(longtemp) %in% longoutcome))]
      }
      droprows <- c()
      for (i in 1:nrow(longtemp)) {
        if (length(which(is.na(longtemp[i, ]))) == length(c(longoutcome,
                                                            timevarying))) {
          droprows <- c(droprows, i)
        }
      }
      if (length(droprows) > 0) {
        longtemp <- longtemp[-droprows, ]
      }
      rownames(longtemp) <- 1:nrow(longtemp)
      rows <- 1:nrow(survivaltemp)
      survivalremove <- c()
      if (length(rows[is.na(survivaltemp[, 2])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    2])])
      }
      if (length(rows[is.na(survivaltemp[, 3])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    3])])
      }
      survivalremove <- unique(survivalremove)
      if (length(survivalremove) > 0) {
        survivaltemp <- survivaltemp[-survivalremove, ]
      }
      rownames(survivaltemp) <- 1:nrow(survivaltemp)
      baselinecols <- which(!(names(dataset) %in% c(names(survivaltemp),
                                                    unlist(timevariable))))
      if (length(baselinecols) > 0) {
        baselinetemp <- dataset[, c(idcol, baselinecols)]
      } else {
        baselinetemp <- NULL
      }
      baselinetemp.surv <- survival[, which(!(names(survival) %in%
                                                c(survtime, cens)))]
    } else {
      if (is.null(time)) {
        stop("No variable specified as holding the
             time variable for the longitudinal data")
      }
      if (is.null(timevarying)) {
        stop("No time varying covariates specified -
             there should at least be time of longitudinal measurement")
      }
      longtemp <- dataset[, which(names(dataset) %in% c(id, longoutcome,
                                                        timevarying))]
      if (length(which(!(names(longtemp) %in% c(id, longoutcome,
                                                time)))) > 0) {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id),
                                 which(names(longtemp) %in% longoutcome), which(names(longtemp) %in%
                                                                                  time), which(!(names(longtemp) %in% c(id, longoutcome,
                                                                                                                        time))))]
      } else {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id),
                                 which(names(longtemp) %in% longoutcome))]
      }
      droprows <- c()
      for (i in 1:nrow(longtemp)) {
        if (length(which(is.na(longtemp[i, ]))) == length(c(longoutcome,
                                                            timevarying))) {
          droprows <- c(droprows, i)
        }
      }
      if (length(droprows) > 0) {
        longtemp <- longtemp[-droprows, ]
      }
      rownames(longtemp) <- 1:nrow(longtemp)
      survivalremove <- c()
      if (length(rows[is.na(survivaltemp[, 2])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    2])])
      }
      if (length(rows[is.na(survivaltemp[, 3])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    3])])
      }
      survivalremove <- unique(survivalremove)
      if (length(survivalremove) > 0) {
        survivaltemp <- survivaltemp[-survivalremove, ]
      }
      rownames(survivaltemp) <- 1:nrow(survivaltemp)
      survivaltemp <- dataset[match(unique(dataset[, idcol]), dataset[,
                                                                      idcol]), which(names(dataset) %in% c(id, survtime, cens))]
      survivaltemp <- survivaltemp[, c(which(names(survivaltemp) %in%
                                               id), which(names(survivaltemp) %in% survtime), which(names(survivaltemp) %in%
                                                                                                      cens))]
      baselinecols <- which(!(names(dataset) %in% c(names(survivaltemp),
                                                    names(longtemp))))
      if (length(baselinecols) > 0) {
        baselinetemp <- dataset[match(unique(dataset[, idcol]),
                                      dataset[, idcol]), c(idcol, baselinecols)]
      } else {
        baselinetemp <- NULL
      }
      }

    } else {
      if (is.null(time)) {
        stop("No variable specified as holding the
             time variable for the longitudinal data")
      }
      if (is.null(timevarying)) {
        stop("No time varying covariates specified -
             there should at least be time of longitudinal measurement")
      }
      if (is.null(longitudinal) || is.null(survival)) {
        stop("One of longitudinal or survival lists of datasets is missing")
      }
      if (class(longitudinal) == "list") {
        longitudinal <- Reduce(function(x, y) merge(x, y, all = TRUE),
                               longitudinal)
      }
      if (class(survival) == "list") {
        survival <- Reduce(function(x, y) merge(x, y, all = TRUE),
                           survival)
      }
      idcol.long <- which(names(longitudinal) %in% id)
      idcol.surv <- which(names(survival) %in% id)
      longtemp <- longitudinal[, which(names(longitudinal) %in% c(id,
                                                                  longoutcome, timevarying))]
      if (length(which(!(names(longtemp) %in% c(id, longoutcome, time)))) >
          0) {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id), which(names(longtemp) %in%
                                                                         longoutcome), which(names(longtemp) %in% time), which(!(names(longtemp) %in%
                                                                                                                                   c(id, longoutcome, time))))]
      } else {
        longtemp <- longtemp[, c(which(names(longtemp) %in% id), which(names(longtemp) %in%
                                                                         longoutcome))]
      }
      droprows <- c()
      for (i in 1:nrow(longtemp)) {
        if (length(which(is.na(longtemp[i, ]))) == length(c(longoutcome,
                                                            timevarying))) {
          droprows <- c(droprows, i)
        }
      }
      if (length(droprows) > 0) {
        longtemp <- longtemp[-droprows, ]
      }
      rownames(longtemp) <- 1:nrow(longtemp)
      baselinetemp.long <- longitudinal[match(unique(longitudinal[, idcol.long]),
                                              longitudinal[, idcol.long]), which(!(names(longitudinal) %in%
                                                                                     c(longoutcome, timevarying)))]
      survivaltemp <- survival[, which(names(survival) %in% c(id, survtime,
                                                              cens))]
      survivaltemp <- survivaltemp[, c(which(names(survivaltemp) %in%
                                               id), which(names(survivaltemp) %in% survtime), which(names(survivaltemp) %in%
                                                                                                      cens))]
      rows <- 1:nrow(survivaltemp)
      survivalremove <- c()
      if (length(rows[is.na(survivaltemp[, 2])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    2])])
      }
      if (length(rows[is.na(survivaltemp[, 3])]) > 0) {
        survivalremove <- c(survivalremove, rows[is.na(survivaltemp[,
                                                                    3])])
      }
      survivalremove <- unique(survivalremove)
      if (length(survivalremove) > 0) {
        survivaltemp <- survivaltemp[-survivalremove, ]
      }
      baselinetemp.surv <- survival[, which(!(names(survival) %in% c(survtime,
                                                                     cens)))]
      baselinetemp <- merge(baselinetemp.long, baselinetemp.surv)
      if (is.null(baseline) == FALSE) {
        if (class(baseline) == "list") {
          baseline <- Reduce(function(x, y) merge(x, y, all = TRUE),
                             baseline)
        }
        baselinetemp <- merge(baselinetemp, baseline)
      }
      }
  longidunique <- unique(longtemp[, 1])
  survidunique <- unique(survivaltemp[, 1])
  baseidunique <- unique(baselinetemp[, 1])
  commonid <- intersect(intersect(longidunique, survidunique), baseidunique)
  longtemp <- longtemp[which(longtemp[, 1] %in% commonid), ]
  survivaltemp <- survivaltemp[which(survivaltemp[, 1] %in% commonid),
                               ]
  baselinetemp <- baselinetemp[which(baselinetemp[, 1] %in% commonid),
                               ]
  if (any(sapply(baselinetemp, "class") == "factor")) {
    baselinetemp <- droplevels(baselinetemp)
  }
  if (any(sapply(longtemp, "class") == "factor")) {
    longtemp <- droplevels(longtemp)
  }
  if (any(sapply(survivaltemp, "class") == "factor")) {
    survivaltemp <- droplevels(survivaltemp)
  }
  longtemp <- longtemp[order(longtemp[, which(names(longtemp) %in% id)]),
                       ]
  survivaltemp <- survivaltemp[order(survivaltemp[, which(names(survivaltemp) %in%
                                                            id)]), ]
  baselinetemp <- baselinetemp[order(baselinetemp[, which(names(baselinetemp) %in%
                                                            id)]), ]
  out <- jointdata(longitudinal = longtemp, survival = survivaltemp,
                   baseline = baselinetemp, id.col = id, time.col = time)
  class(out) <- "jointdata"
  out
}
