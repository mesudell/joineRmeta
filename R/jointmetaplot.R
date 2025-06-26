#' Produce plots of longitudinal and survival outcomes
#'
#' This function can produce plots for each study in the dataset of the
#' longitudinal trajectory pannelled by event type with or without a smoother,
#' and kaplan-meier plots for each study which plot the survival probability
#' against time.
#'
#' @param dataset a \code{jointdata} object
#' @param study the name of the variable holding study membership in the
#'   supplied dataset
#' @param longoutcome the name of the variable holding the longitudinal outcome
#'   in the supplied dataset
#' @param longtime the name of the variable holding the longitudinal time
#'   variable in the supplied dataset
#' @param survtime the name of the variable holding the survival time variable
#'   in the supplied dataset
#' @param cens the name of the variable holding the censoring variable in the
#'   supplied dataset
#' @param id the name of the variable holding the id in the supplied dataset
#' @param smoother a logical indicating whether or not a smoother should be
#'   displayed on the longitudinal plot
#' @param studynames a vector of character strings giving the names to label the
#'   study plots by - the first element of this vector will be the label for the
#'   plots for the first study in the dataset for both the longitudinal and the
#'   survival plots
#' @param type option to select what type of plots should be returned.  If just
#'   plots of the longitudinal trajectories are required then \code{type =
#'   'Longitudinal'}.  Else if just plots of the survival probabilities are
#'   required then \code{type = 'Survival'}.  Finally if both survival and
#'   longitudinal plots are required then this should be set to \code{type =
#'   'Both'}
#' @param eventby an optional character string giving a grouping variable that
#'   the graph of survival probability by time will be split by.
#' @param eventconfint a logical value indicating whether the survival plot
#'   should contain confidence intervals or not.  Defaults to \code{FALSE}.
#'
#' @return Returns an object of class \code{'jointplots'}.  This contains an
#'   element labelled \code{'longplots'} if \code{type} in the function call is
#'   set to one of \code{'Longitudinal'} or \code{'Both'}, and an element
#'   labelled \code{'eventplots'} if \code{type} in the function call is set to
#'   one of \code{'Survival'} or \code{'Both'}.  The element \code{'longplots'}
#'   is a list of ggplot2 objects plotting the longitudinal trajectories for
#'   each study, and is of length equal to the number of studies in the supplied
#'   dataset.  The element \code{'eventplots'} is a list of ggplot2 objects
#'   plotting the survival probabilities for each study and is of length equal
#'   to the number of studies in the supplied dataset.
#'
#'   To plot a particular graph, it can be called by position from the relevent
#'   element of the returned \code{'jointplots'} in the same way that an element
#'   in a particular position is called from a list, or it can be called by name
#'   if \code{study.names} supplied to the function call.
#'
#'   This function supplies separate plots for each study in the dataset.  To
#'   arrange these plots into one grid, use the function
#'   \code{\link{jointmetaplotall}}.
#'
#' @export
#' @import ggplot2 survival
#'
#' @seealso \code{\link[joineR]{jointdata}}, \code{\link[ggplot2]{ggplot}},
#'   \code{\link{jointmetaplotall}}
#'
#' @examples
#'     #change data to jointdata format
#'     jointdat<-tojointdata(longitudinal = simdat$longitudinal,
#'                           survival = simdat$survival, id = 'id',
#'                           longoutcome = 'Y', timevarying = c('time','ltime'),
#'                           survtime = 'survtime', cens = 'cens',
#'                           time = 'time')
#'
#'     #ensure variables are correctly formatted
#'     jointdat$baseline$study <- as.factor(jointdat$baseline$study)
#'     jointdat$baseline$treat <- as.factor(jointdat$baseline$treat)
#'
#'     #produce plots
#'     sepplots<-jointmetaplot(dataset = jointdat, study = 'study',
#'                         longoutcome = 'Y', longtime = 'time',
#'                         survtime = 'survtime', cens = 'cens', id = 'id',
#'                         smoother = TRUE, studynames = c('Study 1', 'Study 2',
#'                         'Study 3', 'Study 4', 'Study 5'), type = 'Both')

jointmetaplot <- function(dataset, study, longoutcome, longtime, survtime,
                          cens, id, smoother = FALSE, studynames = NULL, type = c("Longitudinal",
                                                                                  "Event", "Both"), eventby = NULL, eventconfint = FALSE) {
  if (!inherits(dataset,"jointdata")) {
    stop("Please run tojointdata function before attempting to plot data \n
         this will give data in jointdata format")
  }
  if (missing(study)) {
    stop("Variable holding study membership has not been supplied to study")
  }
  studycol <- which(names(dataset$baseline) %in% study)
  if (missing(longoutcome)) {
    stop("No variable specified as the longitudinal outcome")
  }
  if (missing(cens)) {
    stop("No variable specified as the event indicators to cens")
  }
  if (missing(id)) {
    stop("No variable specified as holding the id variable")
  }
  if (is.null(studynames) == FALSE) {
    if (!inherits(studynames, "character")) {
      stop("studynames should be supplied as a character string")
    }
    if (length(studynames) != length(unique(dataset$baseline[, studycol]))) {
      stop("Different number of studynames supplied to
           number of studies in dataset")
    }
    }
  if (missing(type)) {
    stop("No type of graph requested")
  }
  if (!(type %in% c("Longitudinal", "Event", "Both"))) {
    stop("type should be supplied as one of \"Longitudinal\",
         \"Event\", \"Both\")")
  }
  studies <- as.character(unique(dataset$baseline[, studycol]))
  numstudies <- length(studies)
  ids.bystudy <- lapply(1:numstudies, function(u) {
    dataset$baseline[which(dataset$baseline[, studycol] == studies[u]),
                     which(names(dataset$baseline) %in% id)]
  })
  dataset.bystudy <- lapply(1:numstudies, function(u) {
    out <- subset(dataset, ids.bystudy[[u]])
    class(out) <- "jointdata"
    out
  })
  plots <- list()
  plotnames <- c()
  for (i in 1:numstudies) {
    if (is.null(studynames)) {
      plotnametemp <- paste("studyplot.", studies[[i]], sep = "")
    } else {
      plotnametemp <- paste("studyplot.", studynames[[i]], sep = "")
    }
    plotnames <- c(plotnames, plotnametemp)
  }
  if (type %in% c("Longitudinal", "Both")) {
    longtimemin <- floor(min(dataset$longitudinal[, which(names(dataset$longitudinal) %in%
                                                            longtime)]))
    longtimemax <- ceiling(max(dataset$longitudinal[, which(names(dataset$longitudinal) %in%
                                                              longtime)]))
    longoutmin <- floor(min(dataset$longitudinal[, which(names(dataset$longitudinal) %in%
                                                           longoutcome)]))
    longoutmax <- ceiling(max(dataset$longitudinal[, which(names(dataset$longitudinal) %in%
                                                             longoutcome)]))

    longplots <- lapply(1:numstudies, function(u) {
      datatemp <- merge(dataset.bystudy[[u]]$longitudinal, dataset.bystudy[[u]]$survival)
      if (is.null(dataset.bystudy[[u]]$baseline) == FALSE) {
        datatemp <- merge(datatemp, dataset.bystudy[[u]]$baseline)
      }
      datatemp[, which(names(datatemp) %in% dataset.bystudy[[u]]$longtime)] <- as.numeric(datatemp[,
                                                                                                   which(names(datatemp) %in% dataset.bystudy[[u]]$longtime)])
      facetformula <- as.formula(paste(". ~ ", cens, sep = ""))
      p <- ggplot(data = datatemp, aes(y = .data[[longoutcome]], x = .data[[longtime]],
                                              group = id)) + geom_line() + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
        facet_grid(facetformula) + xlab(longtime) + ylab(longoutcome) +
        xlim(longtimemin, longtimemax) + ylim(longoutmin, longoutmax)
      if (is.null(studynames)) {
        p <- p + ggtitle(studies[[u]])
      } else {
        p <- p + ggtitle(studynames[u])
      }
      if (smoother) {
        p <- p + geom_smooth(aes(group = 1), method = "loess",
                             col = "red", se = FALSE)
      }
      p
    })
    names(longplots) <- plotnames
    plots$longplots <- longplots
  }
  if (type %in% c("Event", "Both")) {
    if (is.null(eventby) == FALSE) {
      if (!(eventby %in% names(dataset$baseline))) {
        stop("Supplied eventby is not in the dataset")
      }
    }
    if (missing(survtime)) {
      stop("No variable specified as the survival time")
    }
    if (is.null(eventby)) {
      eventform <- as.formula(paste("Surv(", survtime, ",", cens,
                                    ")~1"))
    } else {
      eventform <- as.formula(paste("Surv(", survtime, ",", cens,
                                    ")~", eventby))
    }
    survtimemax <- ceiling(max(dataset$survival[, which(names(dataset$survival) %in%
                                                          survtime)]))
    eventplots <- lapply(1:numstudies, function(u) {
      datatemp <- dataset.bystudy[[u]]$survival
      if (is.null(dataset.bystudy[[u]]$baseline) == FALSE) {
        datatemp <- merge(datatemp, dataset.bystudy[[u]]$baseline)
      }
      kmplot <- survfit(eventform, data = datatemp, conf.type = "log-log")
      if (is.null(eventby) == FALSE) {
        tempsum <- summary(kmplot)
        kmplotdata <- do.call(data.frame, lapply(c("time", "surv",
                                                   "strata", "upper", "lower"), function(x) {
                                                     tempsum[x]
                                                   }))

        tempdata <- lapply(1:length(unique(kmplotdata$strata)),
                           function(x) {
                             temp <- kmplotdata[kmplotdata$strata == unique(kmplotdata$strata)[x],
                                                ]
                             for (j in 1:nrow(temp)) {
                               if (j == 1) {
                                 kmplottemp <- temp[1, ]
                               } else {
                                 if (temp$surv[j] == temp$surv[j - 1]) {
                                   kmplottemp <- rbind(kmplottemp, temp[j, ])
                                 } else {
                                   temprow <- temp[(j - 1), ]
                                   temprow$time <- temp$time[j]
                                   kmplottemp <- rbind(kmplottemp, temprow)
                                   kmplottemp <- rbind(kmplottemp, temp[j, ])
                                 }
                               }
                             }
                             kmplottemp
                           })
        kmplotdata <- do.call(rbind, tempdata)

        if (eventconfint) {
          groups <- lapply(1:length(unique(kmplotdata$strata)),
                           function(x) {
                             kmplotdata[kmplotdata$strata == unique(kmplotdata$strata)[x],
                                        ]
                           })
          groups1 <- lapply(1:length(unique(kmplotdata$strata)),
                            function(x) {
                              tempdata <- groups[[x]]
                              temprow <- tempdata[1, ]
                              remove <- c()
                              for (i in 1:nrow(tempdata)) {
                                if (is.na(tempdata$upper[i]) || is.na(tempdata$lower[[i]])) {
                                  remove <- c(remove, i)
                                }
                              }
                              if (length(remove) > 0) {
                                tempdata <- tempdata[-remove, ]
                              }
                              temprow$time <- 0
                              temprow$surv <- 1
                              temprow$upper <- 1
                              temprow$lower <- 1
                              tempdata <- rbind(temprow, tempdata)
                              tempdata
                            })
          for (i in 1:length(groups1)) {
            groups1[[i]]$upper[which(groups1[[i]]$surv %in% (unique(groups1[[i]]$surv)[(length(unique(groups1[[i]]$surv)) -
                                                                                          1):length(unique(groups1[[i]]$surv))]))[-1]] <- NA
            groups1[[i]]$lower[which(groups1[[i]]$surv %in% (unique(groups1[[i]]$surv)[(length(unique(groups1[[i]]$surv)) -
                                                                                          1):length(unique(groups1[[i]]$surv))]))[-1]] <- NA
          }

          kmplotdata <- do.call(rbind, groups1)
        }
        p <- ggplot(data = kmplotdata, aes(x = time, y = surv,
                                           group = strata, colour = strata)) + geom_line() + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
          xlab("Time") + ylab("Survival") + scale_colour_brewer(palette = "Set1",
                                                                name = "Groups") + xlim(0, survtimemax) + ylim(0, 1)

        if (is.null(studynames)) {
          p <- p + ggtitle(studies[[u]])
        } else {
          p <- p + ggtitle(studynames[u])
        }
        if (eventconfint) {
          p <- p + geom_line(aes(y = .data$upper), linetype = 2, na.rm = TRUE)
          p <- p + geom_line(aes(y = .data$lower), linetype = 2, na.rm = TRUE)
        }
        p
      } else {
        time <- kmplot$time
        surv <- kmplot$surv
        kmplotdata <- data.frame(time = time, surv = surv)
        if (eventconfint) {
          kmplotdata$upper <- kmplot$upper
          kmplotdata$lower <- kmplot$lower
        }
        for (j in 1:nrow(kmplotdata)) {
          if (j == 1) {
            kmplottemp <- kmplotdata[1, ]
          } else {
            if (kmplotdata$surv[j] == kmplotdata$surv[j - 1]) {
              kmplottemp <- rbind(kmplottemp, kmplotdata[j, ])
            } else {
              temprow <- kmplotdata[(j - 1), ]
              temprow$time <- kmplotdata$time[j]
              kmplottemp <- rbind(kmplottemp, temprow)
              kmplottemp <- rbind(kmplottemp, kmplotdata[j, ])
            }
          }
        }
        kmplotdata <- kmplottemp
        if (eventconfint) {
          kmplotdata$upper[which(kmplotdata$surv %in% (unique(kmplotdata$surv)[(length(unique(kmplotdata$surv)) -
                                                                                  1):length(unique(kmplotdata$surv))]))[-1]] <- NA

          kmplotdata$lower[which(kmplotdata$surv %in% (unique(kmplotdata$surv)[(length(unique(kmplotdata$surv)) -
                                                                                  1):length(unique(kmplotdata$surv))]))[-1]] <- NA
        }

        p <- ggplot(data = kmplotdata, aes(x = time, y = surv)) +
          geom_line() + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                           legend.position = "bottom") + xlab("Time") + ylab("Survival") +
          xlim(0, survtimemax) + ylim(0, 1)
        if (is.null(studynames)) {
          p <- p + ggtitle(studies[[u]])
        } else {
          p <- p + ggtitle(studynames[u])
        }
        if (eventconfint) {
          p <- p + geom_line(aes(y = .data$upper), linetype = 2, na.rm = TRUE)
          p <- p + geom_line(aes(y = .data$lower), linetype = 2, na.rm = TRUE)
        }
        p
      }
    })
    names(eventplots) <- plotnames
    plots$eventplots <- eventplots
  }
  class(plots) <- "jointplots"
  return(plots)
}
