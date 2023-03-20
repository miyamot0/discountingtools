#' #' get_ED50
#' #'
#' #' Scoring for the log ED50
#' #'
#' #' Methods for routing computations of the ln(ED50). When possible, exact solutions are provided and models without a straightforward exact solution have ed50 calculated numerically using a bisection search.
#' #'
#' #' @param fittingObject core dd fitting object
#' #' @param id id tag
#' #'
#' #' @return natural logarithm of the Effective Delay 50%
#' #' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' get_ED50 <- function(fittingObject, id) {
#'   probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]
#'
#'   message_debug(fittingObject, paste0("Cross Model Metric (ED50)[", probableModel, "]: ", id))
#'
#'   if (probableModel == "noise")          fittingObject = dd_ed50_noise(fittingObject, id)
#'   if (probableModel == "mazur")          fittingObject = dd_ed50_mazur(fittingObject, id)
#'   if (probableModel == "exponential")    fittingObject = dd_ed50_exponential(fittingObject, id)
#'   if (probableModel == "laibson")        fittingObject = dd_ed50_laibson(fittingObject, id)
#'   if (probableModel == "greenmyerson")   fittingObject = dd_ed50_greenmyerson(fittingObject, id)
#'   if (probableModel == "rachlin")        fittingObject = dd_ed50_rachlin(fittingObject, id)
#'   if (probableModel == "ebertprelec")    fittingObject = dd_ed50_ebertprelec(fittingObject, id)
#'   if (probableModel == "bleichrodt")     fittingObject = dd_ed50_bleichrodt(fittingObject, id)
#'   if (probableModel == "rodriguezlogue") fittingObject = dd_ed50_rodriguezlogue(fittingObject, id)
#'
#'   fittingObject
#' }
