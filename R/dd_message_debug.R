#' messageDebug
#'
#' Extension of message method, instead yolked to a flag defining level of verbosity.
#'
#' @param fittingObject core fitting object
#' @param msg (char) message
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
messageDebug <- function(fittingObject, msg) {
  if (fittingObject[[ "verbose"  ]] == TRUE) message(msg)
}
