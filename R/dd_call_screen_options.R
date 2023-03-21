#' dd_screen_options
#'
#' This call applies screening criteria to a data dataset. Specifically, it can be used to apply criteria (no filtering) or apply criteria and filter based on one or more criteria (e.g., JB1, JB2)
#'
#' @param fittingObject core fitting object
#' @param screen (bool) set screen TRUE or FALSE (i.e. NULL)
#' @param JB1Flag (num) bounce constant per authors (set at initial defaults)
#' @param JB2Flag (num) extremity change constant per authors (set at initial defaults)
#' @param filterPassing (char vector) which JB criteria to retain in dataset, e.g. c("JB1", "JB2")
#'
#' @export
dd_screen <- function(fittingObject, screen = TRUE, JB1Flag = 0.2, JB2Flag = 0.1, filterPassing = NULL) {
  message_debug(fittingObject, "Setting Screening Options")
  fittingObject[[ "screen"  ]] = screen
  fittingObject[[ "JB1Flag" ]] = JB1Flag
  fittingObject[[ "JB2Flag" ]] = JB2Flag

  # TODO: validate passing
  #

  if (!is.logical(screen))                         stop('screen must be a boolean')
  if (!is.numeric(JB1Flag))                        stop('JB1Flag must be numeric')
  if (!is.numeric(JB2Flag))                        stop('JB2Flag must be numeric')

  if (screen == TRUE) {
    if (!is.null(filterPassing)) {
      if (!(filterPassing %in% c('JB1', 'JB2')))   stop('Only `JB1` or `JB2` screening supported')

      fittingObject[[ "filterPassing" ]] = filterPassing
    }
  }

  fittingObject
}
