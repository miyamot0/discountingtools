#' Scoring for the most probable model area
#'
#' In this set of methods, the area beneath the fitted model is calculated and divided by the maximum area using numerical integration methods.  All delays are calculated in the normal scale.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return area beneath the fitted model
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
get_MBAUC <- function(fittingObject, id) {
  probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]

  message_debug(fittingObject, paste0("Cross Model Metric (MBAUC)[", probableModel, "]: ", id))

  if (probableModel == "noise")          fittingObject = dd_mbauc_noise(fittingObject, id)
  if (probableModel == "mazur")          fittingObject = dd_mbauc_mazur(fittingObject, id)
  if (probableModel == "exponential")    fittingObject = dd_mbauc_exponential(fittingObject, id)
  if (probableModel == "laibson")        fittingObject = dd_mbauc_laibson(fittingObject, id)
  if (probableModel == "greenmyerson")   fittingObject = dd_mbauc_greenmyerson(fittingObject, id)
  if (probableModel == "rachlin")        fittingObject = dd_mbauc_rachlin(fittingObject, id)
  if (probableModel == "ebertprelec")    fittingObject = dd_mbauc_ebertprelec(fittingObject, id)
  if (probableModel == "bleichrodt")     fittingObject = dd_mbauc_bleichrodt(fittingObject, id)
  if (probableModel == "rodriguezlogue") fittingObject = dd_mbauc_rodriguezlogue(fittingObject, id)

  fittingObject
}
