#' dd_metric_options
#'
#' This call builds out the metrics used to conduct cross-model comparisons (e.g., ED50, MBAUC)
#'
#' @param fittingObject core dd fitting object
#' @param metrics (char vector) This vector contains a list of possible cross-model metrics.
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_metric_options <- function(fittingObject, metrics) {
  message_debug(fittingObject, "Setting Cross Model Metrics")
  fittingObject[[ "metrics" ]] = metrics

  fittingObject
}
