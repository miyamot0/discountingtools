
#' dd_probableModel
#'
#' Return model probabilities
#'
#' @param fittingObject
#' @param id
#'
#' @return
#' @export
dd_probableModel <- function(fittingObject, id) {

  modelComparison     = list(
    Models            = fittingObject$models,
    ProbableModel     = NA,
    ProbableModelProb = NA,
    BFSum             = 0.0,
    BFs               = list(),
    Probs             = list(),
    BICs              = list()
  )

  paramList <- c("Noise.BIC")

  for (model in as.character(fittingObject$models)) {

    #print(model)

    if (model == "mazur")          paramList = c(paramList, "Mazur.BIC")
    if (model == "exponential")    paramList = c(paramList, "Exponential.BIC")
    if (model == "laibson")        paramList = c(paramList, "Laibson.BIC")
    if (model == "greenmyerson")   paramList = c(paramList, "GreenMyerson.BIC")
    if (model == "rachlin")        paramList = c(paramList, "Rachlin.BIC")
    if (model == "ebertprelec")    paramList = c(paramList, "EbertPrelec.BIC")
    if (model == "bleichrodt")     paramList = c(paramList, "Bleichrodt.BIC")
    if (model == "rodriguezlogue") paramList = c(paramList, "RodriguezLogue.BIC")
  }

  print(paramList)

  # TODO: Check if Noise.BIC == Inf

  fittingObject
}
