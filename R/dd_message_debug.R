#' message_debug
#'
#' Extension of message method, instead yolked to a flag defining level of verbosity.
#'
#' @param fittingObject core fitting object
#' @param msg (char) message
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
message_debug <- function(fittingObject, msg) {
  if (fittingObject[[ "verbose"  ]] == TRUE) message(msg)
}
