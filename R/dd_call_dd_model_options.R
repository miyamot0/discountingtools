#' dd_model_options
#'
#' This call builds out the model(s) included in subsequent regressions.
#'
#' @param fittingObject core dd fitting object
#' @param plan (char vector) This vector contains a list of possible model candidates.
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_model_options <- function(fittingObject, plan) {
  message_debug(fittingObject, "Setting Model Options")

  fittingObject[[ "models" ]] = plan

  fittingObject
}
