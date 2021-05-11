
#' dd_probableModel
#'
#' Return model probabilities
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_probableModel <- function(fittingObject, id) {

  modelComparison     = list(
    BFs               = list(),
    BFSum             = 0.0,
    Probs             = list(),
    ProbableModel     = NA,
    ProbableModelProb = NA
  )

  currentResults = fittingObject$results[[as.character(id)]]

  # Perfect fit for noise model, hacky workaround
  if (currentResults$noise$BIC == Inf) {
    for (model in as.character(fittingObject$models)) {
      modelComparison$BFs[[ model ]] = NULL
      modelComparison$BFSum = NULL
    }

    for (model in as.character(fittingObject$models))
      modelComparison$Probs[[ model ]] = 0

    modelComparison$ProbableModel     = "noise"
    modelComparison$ProbableModelProb = 1

  } else {
    for (model in as.character(fittingObject$models)) {
      modelComparison$BFs[[ model ]] = exp(-.5*(currentResults[[model]]$BIC - currentResults$noise$BIC))
      modelComparison$BFSum = modelComparison$BFSum + modelComparison$BFs[[ model ]]
    }

    for (model in as.character(fittingObject$models))
      modelComparison$Probs[[ model ]] = modelComparison$BFs[[ model ]] / modelComparison$BFSum
  }

  sortedProbs = sort(unlist(modelComparison[["Probs"]]), decreasing = TRUE)

  mostProbModel <- names(sortedProbs)[1]

  fittingObject$rotation[[as.character(id)]] = list(
    ProbableModel       = mostProbModel,
    ProbableModel.BF    = modelComparison$BFs[[mostProbModel]],
    ProbableModel.Prob  = modelComparison$Probs[[mostProbModel]]
  )

  fittingObject
}
