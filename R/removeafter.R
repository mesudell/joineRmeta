#' Code to remove longitudinal information recorded after survival outcome
#'
#' This function is designed to remove any longitudinal information recorded
#' after the survival time for each individual.  If the survival event is not
#' terminal, it is possible that longitudinal information is available in the
#' data after the survival time, but it should not contribute to the joint
#' analysis.  This function takes and returns a \code{jointdata} object.
#'
#' @param data a jointdata object (see \code{\link[joineR]{jointdata}})
#' @param longitudinal a character string denoting name of the variable holding
#'   the longitudinal outcome of interest.
#' @param survival a character string denoting the name of the variable holding
#'   the survival time for the event of interest.
#' @param id a character string denoting the name of the variable holding the id
#'   variable for the data.
#' @param time a character string denoting the name of the variable holding the
#'   longitudinal time variable.
#'
#' @return a jointdata object, see \code{\link[joineR]{jointdata}}
#'
#' @export
#'
#' @details This function removes any longitudinal information recorded for an
#'   individual after their survival time. A joint data object should have the
#'   id column as the first column in each of the \code{survival},
#'   \code{longitudinal} and \code{baseline} datasets.  In the \code{survival}
#'   dataset the second column should be the survival time and the third should
#'   be the censoring variable.  In the \code{longitudinal} dataset, the second
#'   column should be longitudinal outcome, the third the longitudinal time
#'   variable and the remaining columns any other time varying covariates.  The
#'   baseline dataset should have columns 2 and onwards containing time
#'   stationary covariates such as treatment group assignment or study
#'   membership.
#'
#'   This function does not need to be run on the results of the multi-study
#'   data simulation function \code{\link{simjointmeta}}, because the
#'   longitudinal data simulated under this function is already capped at the
#'   individual's survival time.
#'
#'
#' @seealso \code{\link[joineR]{jointdata}}
#'
#' @examples
#'  \dontrun{
#'  #the dataset simdat3 in this package contains joint data where longitudinal
#'  #data exists after individual's survival times.
#'  str(simdat3)
#'
#'  #first this data needs to be changed to a jointdata object
#'  jointdat3<-tojointdata(longitudinal = simdat3$longitudinal,
#'                   survival = simdat3$survival, id = 'id', longoutcome = 'Y',
#'                   timevarying = c('time','ltime'), survtime = 'survtime',
#'                   cens = 'cens',time = 'time')
#'
#'  #then additional data recorded after the survival time can be removed
#'  jointdat3.1<-removeafter(data = jointdat3, longitudinal = 'Y',
#'                   survival = 'survtime', id = 'id', time = 'time')
#'
#'  #we can compare the two datasets to see the removed data
#'  str(jointdat3)
#'  str(jointdat3.1)
#'  }
#'
removeafter <- function(data, longitudinal, survival, id, time) {
  if (class(data) != "jointdata") {
    stop("data should be supplied in jointdata format")
  }
  if (class(survival) != "character") {
    stop("survival should be the character name
         of the time-to-event outcome of interest")
  }
  if (class(longitudinal) != "character") {
    stop("longitudinal should be the character name
         of the longitudinal outcome of interest")
  }
  if (class(id) != "character") {
    stop("id should be the character name
         of the id variable of the data")
  }
  if (class(time) != "character") {
    stop("time should be the character name
         of the longitudinal time variable of the data")
  }
  if (!(survival %in% names(data$survival))) {
    stop("Supplied survival variable is not in survival data
         in supplied dataset")
  }
  if (!(longitudinal %in% names(data$longitudinal))) {
    stop("Supplied longitudinal variable is not in longitudinal data
         in supplied dataset")
  }
  if (!(time %in% names(data$longitudinal))) {
    stop("Supplied longitudinal time variable is not in longitudinal data
         in supplied dataset")
  }
  if (!(id %in% names(data$longitudinal)) || !(id %in% names(data$survival))) {
    stop("Supplied id variable is not missing from at least one of
         the longitudinal data and the survival data in supplied dataset")
  }
  idcol.long <- which(names(data$longitudinal) %in% id)
  idcol.surv <- which(names(data$survival) %in% id)
  survcol <- which(names(data$survival) %in% survival)
  longcol <- which(names(data$longitudinal) %in% longitudinal)
  longtimecol <- which(names(data$longitudinal) %in% time)
  uniqueids <- data$survival[, idcol.surv]
  templongdat <- data$longitudinal[1, ]
  pb <- txtProgressBar(min = 0, max = length(uniqueids), style = 3)
  counter <- 1
  for (i in 1:length(uniqueids)) {
    tempdat <- data$longitudinal[which(data$longitudinal[, idcol.long] %in%
                                         uniqueids[i]), ]
    survtemp <- data$survival[, survcol][which(data$survival[, idcol.surv] %in%
                                                 uniqueids[i])]
    tempdat <- tempdat[which(tempdat[, longtimecol] <= survtemp), ]
    templongdat <- rbind(templongdat, tempdat)
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  close(pb)
  templongdat <- templongdat[-1, ]
  if (length(which(is.na(templongdat[, longcol])) == TRUE) > 0) {
    templongdat <- templongdat[-which(is.na(templongdat[, longcol])),
                               ]
  }
  rownames(templongdat) <- 1:nrow(templongdat)
  data$longitudinal <- templongdat
  uniqueids <- unique(data$longitudinal[, idcol.long])
  data <- subset(data, uniqueids)
  return(data)
}
