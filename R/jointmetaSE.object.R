#' A \code{jointmeta1SE} object
#'
#' An object returned by the \code{jointmeta1SE} function, inheriting from class
#' \code{jointmeta1SE} representing the results of bootstrapping a fit from the
#' \code{jointmeta1} function.  Objects of this class have methods for the
#' \code{\link{print}} and \code{\link{vcov}} functions.
#'
#' @author Maria Sudell (\email{mesudell@@liverpool.ac.uk})
#' @seealso \code{\link{jointmeta1}}, \code{\link{jointmeta1SE}}.
#' @return A list with the following components. \describe{
#'
#'   \item{\code{results}}{a data frame containing the estimates, standard
#'   errors and 95\% confidence intervals for the parameters from the model and
#'   any overall effects requested.}
#'
#'   \item{\code{covmat}}{the covariance matrix for the model parameters}
#'
#'   \item{\code{bootstraps}}{a data frame containing the results of each
#'   bootstrap}
#'
#'   }
#'
#'

"jointmeta1SE.object" <- NULL
