#' fit_dd_curves
#'
#' This is the entry point for users. It constructs a core fitting object that is passed through the program, with branching options based on those specified by the user.
#'
#' @param data (dataframe) assigned data
#' @param settings (named list) mappings
#' @param maxValue (num) A parameter
#' @param verbose (bool) output level (default FALSE)
#' @param strategy (char) fit to individual ids (default) or group
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom rlang enexpr
#' @export
fit_dd_curves <- function(data, settings, maxValue, strategy = "ind", verbose = FALSE) {

  fittingObject = list()                             # Primary object
  fittingObject[[ "settings" ]] = enexpr(settings)   # Settings
  fittingObject[[ "data"     ]] = data               # Stored data
  fittingObject[[ "models"   ]] = character(0)       # Model selections
  fittingObject[[ "strategy" ]] = strategy           # Analytical strategy
  fittingObject[[ "metrics"  ]] = character(0)       # Cross-model Metrics
  fittingObject[[ "results"  ]] = list()             # Result frame
  fittingObject[[ "maxValue" ]] = maxValue           # Max level (A)
  fittingObject[[ "verbose"  ]] = verbose            # Output level

  class(fittingObject) <- c("discountingtools")

  if (is.null(fittingObject$settings[["Delays"]]))     stop('No Delays aesthetic specified')
  if (is.null(fittingObject$settings[["Values"]]))     stop('No Values aesthetic specified')
  if (is.null(fittingObject$settings[["Individual"]])) stop('No Individual aesthetic specified')

  fittingObject
}