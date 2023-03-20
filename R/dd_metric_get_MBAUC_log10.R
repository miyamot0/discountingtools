#' Scoring for the most probable model area, in log10 space
#'
#' In this set of methods, the area beneath the fitted model is calculated and divided by the maximum area using numerical integration methods.  All delays are calculated in the log base 10 scale.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return area beneath the fitted model, in log10 space
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
get_MBAUC_log10 <- function(fittingObject, id) {
  probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]

  message_debug(fittingObject, paste0("Cross Model Metric (Log10 MBAUC)[", probableModel, "]: ", id))

  if (probableModel == "noise")          fittingObject = dd_mbauc_log10_noise(fittingObject, id)
  if (probableModel == "mazur")          fittingObject = dd_mbauc_log10_mazur(fittingObject, id)
  if (probableModel == "exponential")    fittingObject = dd_mbauc_log10_exponential(fittingObject, id)
  if (probableModel == "laibson")        fittingObject = dd_mbauc_log10_laibson(fittingObject, id)
  if (probableModel == "greenmyerson")   fittingObject = dd_mbauc_log10_greenmyerson(fittingObject, id)
  if (probableModel == "rachlin")        fittingObject = dd_mbauc_log10_rachlin(fittingObject, id)
  if (probableModel == "ebertprelec")    fittingObject = dd_mbauc_log10_ebertprelec(fittingObject, id)
  if (probableModel == "bleichrodt")     fittingObject = dd_mbauc_log10_bleichrodt(fittingObject, id)
  if (probableModel == "rodriguezlogue") fittingObject = dd_mbauc_log10_rodriguezlogue(fittingObject, id)

  fittingObject
}
