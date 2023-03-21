#' fit_dd_curves
#'
#' This is the entry point for users. It constructs a core fitting object that is passed through the program, with branching options based on those specified by the user.
#'
#' @param data (dataframe) assigned data
#' @param settings (named list) mappings
#' @param maxValue (num) A parameter
#' @param verbose (bool) output level (default FALSE)
#' @param strategy (char) fit to individual ids (default) or group
#' @param plan (char vector) This vector contains a list of possible model candidates.
#' @param metrics (char vector) This vector contains a list of possible cross-model metrics.
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom rlang enexpr
#' @export
fit_dd_curves <- function(data, settings, maxValue = NULL, strategy = "ind", verbose = FALSE, plan = NULL, metrics = c('lned50', 'mbauc', 'logmbauc')) {

  cached_settings = enexpr(settings)

  fittingObject = list()                             # Primary object
  fittingObject[[ "settings" ]] = cached_settings    # Settings
  fittingObject[[ "data"     ]] = data               # Stored data
  fittingObject[[ "strategy" ]] = strategy           # Analytical strategy
  fittingObject[[ "metrics"  ]] = character(0)       # Cross-model Metrics
  fittingObject[[ "results"  ]] = list()             # Result frame
  fittingObject[[ "maxValue" ]] = maxValue           # Max level (A)
  fittingObject[[ "verbose"  ]] = verbose            # Output level
  fittingObject[[ "models"   ]] = plan
  fittingObject[[ "metrics"  ]] = metrics

  class(fittingObject) <- c("discountingtools")

  if (!("Delays" %in% names(cached_settings)))     stop('No Delays aesthetic specified')
  if (!("Values" %in% names(cached_settings)))     stop('No Values aesthetic specified')
  if (!("Individual" %in% names(cached_settings))) stop('No Individual aesthetic specified')

  if (is.null(fittingObject[["models"]]))          stop('No models specified')
  if (is.null(fittingObject[["maxValue"]]))        stop('No maximum value specified')

  if (!(strategy %in% c('ind', 'group')))          stop('Only `ind` or `group` strategies supported')

  fittingObject
}
