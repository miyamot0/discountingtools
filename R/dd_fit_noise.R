
#' dd_fit_noise
#'
#' This fits an intercept only model to the data. Its trash, but its a testable alternative that inferring usefulness from an R2 value
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_noise <- function(fittingObject, id) {

  modelResults = list(
    Model      = "noise",
    Intercept  = NA,
    RMSE       = NA,
    BIC        = NA,
    AIC        = NA,
    ED50       = NA,
    MBAUC      = NA,
    Log10MBAUC = NA
  )

  modelFitNoise <- NULL

  currentData = fittingObject$data[
    which(fittingObject$data[,
      as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]
  currentData$ddY = currentData[,as.character(fittingObject$settings['Values'])]
  currentData$ddY = currentData$ddY / as.numeric(fittingObject[[ "maxValue" ]])

  try(modelFitNoise <- stats::lm(ddY ~ 1, currentData), silent = TRUE)

  if (!is.null(modelFitNoise)) {
    modelResults[[ "Intercept"   ]]  = modelFitNoise$coefficients[["(Intercept)"]]
    modelResults[[ "RMSE"        ]]  = summary(modelFitNoise)[["sigma"]]
    modelResults[[ "ED50"        ]]  = NA
    modelResults[[ "MBAUC"       ]]  = modelFitNoise$coefficients[["(Intercept)"]]
    modelResults[[ "Log10MBAUC"  ]]  = modelFitNoise$coefficients[["(Intercept)"]]
    modelResults[[ "BIC"         ]]  = ifelse(summary(modelFitNoise)[["sigma"]] == 0,
                                              Inf,
                                              stats::BIC(modelFitNoise))
    modelResults[[ "AIC"         ]]  = ifelse(summary(modelFitNoise)[["sigma"]] == 0,
                                              Inf,
                                              stats::AIC(modelFitNoise))
  }

  fittingObject$results[[as.character(id)]][["noise"]] = modelResults

  fittingObject
}
