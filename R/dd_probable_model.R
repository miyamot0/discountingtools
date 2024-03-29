#' dd_probable_model
#'
#' This method is used to perform approximate Bayesian model selection using extracted Bayes Factors from calculated BIC values.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_probable_model <- function(fittingObject, id) {

  modelComparison     = list(
    BFs               = list(),
    BFSum             = 0.0,
    Probs             = list(),
    ProbableModel     = NA,
    ProbableModelProb = NA
  )

  currentResults = fittingObject$results[[as.character(id)]]

  # Perfect fit for noise model, hacky workaround
  if (is.infinite(currentResults$noise$BIC)) {
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

  for (metric in fittingObject[["metrics"]]) {
    if (metric == "lned50")   fittingObject$ed50[[as.character(id)]]       = currentResults[[mostProbModel]]$ED50
    if (metric == "mbauc")    fittingObject$mbauc[[as.character(id)]]      = currentResults[[mostProbModel]]$MBAUC
    if (metric == "logmbauc") fittingObject$mbauclog10[[as.character(id)]] = currentResults[[mostProbModel]]$Log10MBAUC
  }

  fittingObject
}
