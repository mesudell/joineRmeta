#' Arrange study plots into a grid
#'
#' This function is designed to take the output from \code{\link{jointmetaplot}}
#' and output the study plots of each type arranged into a grid.
#'
#' @param plotlist the output from running the \code{\link{jointmetaplot}}
#'   function.
#' @param ncol the number of columns of the grid to arrange the plots in.  This
#'   must be supplied to the function
#' @param nrow the number of rows of the grid to arrange the plot in.  This is
#'   an optional parameter, which if not supplied is calculated in the function
#'   based on the number of supplied plots and the specified value of
#'   \code{ncol}.
#' @param top a character string to act as the title for the plots
#' @param type option to select what type of plots should be returned.  If just
#'   the grid of the longitudinal trajectories are required then \code{type =
#'   'Longitudinal'}.  Else if just the grid of the survival probabilities
#'   graphs are required then \code{type = 'Survival'}.  Finally if grids of
#'   both survival and longitudinal plots are required then this should be set
#'   to \code{type = 'Both'}.  If both, then the same title as supplied to
#'   \code{top} will be used, similarly for \code{ncol} and \code{nrow}.
#'
#' @return An object of class \code{'jointplotsall'} is returned.  If in the
#'   function call \code{type = 'Longitudinal'} or \code{type = 'Both'} then the
#'   element in the returned object names \code{'longall'} is the arranged grid
#'   of longitudinal trajectory plots from each study in the dataset.  If
#'   \code{type = 'Survival'} or \code{type = 'Both'} then the element in the
#'   returned object labelled \code{'eventsall'} is the arranged grid of the
#'   survival probability plots from each study in the dataset.  The arranged
#'   grids can either be printed by name, or by extracting them as you would an
#'   element from a list.
#'
#' @export
#' @import ggplot2 grid gridExtra
#' @importFrom grDevices recordPlot
#'
#'
#' @seealso \code{\link{jointmetaplot}}
#'
#' @examples
#'     \dontrun{
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
#'     #note that inclusion of a smoother sometime results in error messages
#'     #see ggplot2 for error message interpretation
#'     sepplots<-jointmetaplot(dataset = jointdat, study = 'study',
#'                         longoutcome = 'Y', longtime = 'time',
#'                         survtime = 'survtime', cens = 'cens', id = 'id',
#'                         smoother = TRUE, studynames = c('Study 1', 'Study 2',
#'                         'Study 3', 'Study 4', 'Study 5'), type = 'Both')
#'
#'     allplot2<-jointmetaplotall(plotlist = sepplots, ncol =2,
#'              top = 'All studies', type = 'Both')
#'    }
#'
jointmetaplotall <- function(plotlist, ncol, nrow = NULL, top = NULL, type = c("Longitudinal",
                                                                               "Event", "Both")) {
  if (class(plotlist) != "jointplots") {
    stop("plotlist should be a list of plots as produced by jointmetaplot,
         of class jointplots")
  }
  if (class(ncol) != "numeric") {
    stop("ncol should be a numeric value")
  }
  if (!(floor(ncol) == ncol)) {
    stop("ncol should be a whole number")
  }
  if (is.null(top)) {
    top <- "All Studies"
  } else {
    if (class(top) != "character") {
      stop("If supplied, top should be a character string")
    }
  }
  if (missing(type)) {
    stop("No type of graph requested")
  }
  if (!(type %in% c("Longitudinal", "Event", "Both"))) {
    stop("type should be supplied as one of \" Longitudinal \",
         \"Event\", \"Both\")")
  }
  if (is.null(nrow) == FALSE) {
    if (class(nrow) != "numeric") {
      stop("nrow should be a numeric value")
    }
    if (!(floor(nrow) == nrow)) {
      stop("nrow should be a whole number")
    }
    if ((ncol * nrow) < length(plotlist)) {
      stop("Supplied ncol and nrow don't
           leave enough spaces for all the plots in plotlist")
    }
    } else {
      if (is.null(plotlist$longplots) == FALSE) {
        nrow <- ceiling(length(plotlist$longplots)/ncol)
      } else {
        nrow <- ceiling(length(plotlist$eventplots)/ncol)
      }
  }
  plots <- list()
  if (type %in% c("Longitudinal", "Both")) {
    if (is.null(plotlist$longplots)) {
      stop("No longitudinal plots in supplied plot list")
    }
    do.call(grid.arrange, c(plotlist$longplots, list(ncol = ncol, nrow = nrow,
                                                     top = textGrob(top, gp = gpar(fontsize = 20, font = 2)))))
    plots$longall <- recordPlot()
  }
  if (type %in% c("Event", "Both")) {
    if (is.null(plotlist$eventplots)) {
      stop("No time-to-event plots in supplied plot list")
    }
    eventplots <- plotlist$eventplots
    tmp <- ggplot_gtable(ggplot_build(eventplots[[1]]))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) {
      eventlegend <- tmp$grobs[[leg]]
      eventplots <- lapply(1:length(eventplots), function(u) {
        tempplot <- eventplots[[u]]
        tempplot <- tempplot + theme(legend.position = "none")
      })
      blankplot <- ggplot(data.frame()) + geom_blank() + theme_bw() +
        theme(panel.border = element_blank())
      if (length(eventplots) < (ncol * nrow)) {
        for (i in (length(eventplots) + 1):(ncol * nrow)) {
          eventplots[[i]] <- blankplot
        }
      }
      eventplots[[length(eventplots) + 1]] <- eventlegend
      layoutset <- rbind(do.call(rbind, split(1:(ncol * nrow), ceiling(seq_along(1:(ncol *
                                                                                      nrow))/ncol))), rep((ncol * nrow) + 1, ncol))
      heights <- c(rep(2.5, nrow), 0.2)
      nrow <- nrow + 1
    } else {
      blankplot <- ggplot(data.frame()) + geom_blank() + theme_bw() +
        theme(panel.border = element_blank())
      if (length(eventplots) < (ncol * nrow)) {
        for (i in (length(eventplots) + 1):(ncol * nrow)) {
          eventplots[[i]] <- blankplot
        }
      }
      layoutset <- do.call(rbind, split(1:(ncol * nrow), ceiling(seq_along(1:(ncol *
                                                                                nrow))/ncol)))
      heights <- rep(2.5, nrow)
    }
    do.call(grid.arrange, c(eventplots, list(ncol = ncol, nrow = nrow,
                                             layout_matrix = layoutset, top = textGrob(top, gp = gpar(fontsize = 20,
                                                                                                      font = 2)), heights = heights)))
    plots$eventsall <- recordPlot()
  }
  class(plots) <- "jointplotsall"
  plots
}
