##
## Copyright 2021 Shawn Gilroy
##
## This file is part of discountingtools.
##
## discountingtools is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, version 2.
##
## discountingtools is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with discountingtools. If not, see <http://www.gnu.org/licenses/gpl-2.0.html>.

library(rlang)
library(dplyr)

#' fitDDCurves
#'
#' Core fitting object
#'
#' @param data assigned data
#' @param settings mappings
#' @param maxValue A parameter
#' @param verbose output level (default FALSE)
#'
#' @return
#' @export
fitDDCurves <- function(data, settings, maxValue, verbose = FALSE) {

  fittingObject = list()                             # Primary object
  fittingObject[[ "settings" ]] = enexpr(settings)   # Aesthetics
  fittingObject[[ "data"     ]] = data               # Stored data
  fittingObject[[ "models"   ]] = character(0)       # Model selections
  fittingObject[[ "strategy" ]] = "ind"              # Analytical strategy
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

#' dd_modelOptions
#'
#' method to specify which models to include as candidates
#'
#' @param fittingObject core dd fitting object
#' @param plan vector of model candidates
#'
#' @return
#' @export
dd_modelOptions <- function(fittingObject, plan) {
  messageDebug(fittingObject, "Setting Model Options")

  fittingObject[[ "models" ]] = plan

  fittingObject
}

#' dd_metricOptions
#'
#' method to indicate metrics of interest in the analysis
#'
#' @param fittingObject core dd fitting object
#' @param metrics vector specifying metrics
#'
#' @return
#' @export
dd_metricOptions <- function(fittingObject, metrics) {
  messageDebug(fittingObject, "Setting Cross Model Metrics")
  fittingObject[[ "metrics" ]] = metrics

  fittingObject
}

#' dd_analyze
#'
#' this is the business end of the analytical process
#'
#' @param fittingObject core dd fitting object
#'
#' @return
#' @export
dd_analyze <- function(fittingObject, modelSelection = TRUE) {
  messageDebug(fittingObject, "Beginning Model Fitting(s)")

  fittingObject[[ "ModelSelection" ]] = modelSelection

  # Add in noise model as a comparator
  if (!("noise" %in% fittingObject[["models"]]) & modelSelection == TRUE)
    fittingObject[["models"]] = c("noise", fittingObject[["models"]])

  # loop through individual id's
  for (id in unique(fittingObject$data[[as.character(fittingObject$settings['Individual'])]])) {
    messageDebug(fittingObject, paste("Fitting:", id))

    fittingObject$results[[as.character(id)]] = list()

    for (model in fittingObject[["models"]]) {
      messageDebug(fittingObject, paste("Fitting:", id, "Rotation:", model))

      if (model == "noise")          fittingObject = dd_fit_noise(          fittingObject, id)
      if (model == "mazur")          fittingObject = dd_fit_mazur(          fittingObject, id)
      if (model == "exponential")    fittingObject = dd_fit_exponential(    fittingObject, id)
      if (model == "laibson")        fittingObject = dd_fit_laibson(        fittingObject, id)
      if (model == "greenmyerson")   fittingObject = dd_fit_greenmyerson(   fittingObject, id)
      if (model == "rachlin")        fittingObject = dd_fit_rachlin(        fittingObject, id)
      if (model == "ebertprelec")    fittingObject = dd_fit_ebertprelec(    fittingObject, id)
      if (model == "bleichrodt")     fittingObject = dd_fit_bleichrodt(     fittingObject, id)
      if (model == "rodriguezlogue") fittingObject = dd_fit_rodriguezlogue( fittingObject, id)

    }

    if (modelSelection) {
      fittingObject = dd_probableModel(fittingObject, id)

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   fittingObject = getED50(fittingObject, id)
        if (metric == "mbauc")    fittingObject = getMBAUC(fittingObject, id)
        if (metric == "logmbauc") fittingObject = getMBAUCLog10(fittingObject, id)
      }
    }
  }

  fittingObject
}

#' dd_probableModel
#'
#' Return model probabilities
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
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

  fittingObject
}

#' Scoring for the log ED50
#'
#' Methods for routing computations of the ln(ED50). When possible, exact solutions are provided and models without a straightforward exact solution have ed50 calculated numerically using a bisection search.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#'
#' @return natural logarithm of the Effective Delay 50%
getED50 <- function(fittingObject, id) {
  probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]

  messageDebug(fittingObject, paste0("Cross Model Metric (ED50)[", probableModel, "]: ", id))

  if (probableModel == "noise")          fittingObject = dd_ed50_noise(fittingObject, id)
  if (probableModel == "mazur")          fittingObject = dd_ed50_mazur(fittingObject, id)
  if (probableModel == "exponential")    fittingObject = dd_ed50_exponential(fittingObject, id)
  if (probableModel == "laibson")        fittingObject = dd_ed50_laibson(fittingObject, id)
  if (probableModel == "greenmyerson")   fittingObject = dd_ed50_greenmyerson(fittingObject, id)
  if (probableModel == "rachlin")        fittingObject = dd_ed50_rachlin(fittingObject, id)
  if (probableModel == "ebertprelec")    fittingObject = dd_ed50_ebertprelec(fittingObject, id)
  if (probableModel == "bleichrodt")     fittingObject = dd_ed50_bleichrodt(fittingObject, id)
  if (probableModel == "rodriguezlogue") fittingObject = dd_ed50_rodriguezlogue(fittingObject, id)

  fittingObject
}

#' Scoring for the most probable model area
#'
#' In this set of methods, the area beneath the fitted model is calculated and divided by the maximum area using numerical integration methods.  All delays are calculated in the normal scale.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#'
#' @return area beneath the fitted model
getMBAUC <- function(fittingObject, id) {
  probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]

  messageDebug(fittingObject, paste0("Cross Model Metric (MBAUC)[", probableModel, "]: ", id))

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

#' Scoring for the most probable model area, in log10 space
#'
#' In this set of methods, the area beneath the fitted model is calculated and divided by the maximum area using numerical integration methods.  All delays are calculated in the log base 10 scale.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#'
#' @return area beneath the fitted model, in log10 space
getMBAUCLog10 <- function(fittingObject, id) {
  probableModel = fittingObject$rotation[[as.character(id)]][["ProbableModel"]]

  messageDebug(fittingObject, paste0("Cross Model Metric (Log10 MBAUC)[", probableModel, "]: ", id))

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

#' Generalized residual call
#'
#' General, shared method for coordinating nls.lm fitting calls. Routes a supplied "valueFunction" with observed data and supplied parameters.
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return residual value of referenced function
residualFunction <- function(params, x, value, valueFunction, jacobianFunction)
{
  value - do.call("valueFunction", c(list(x = x), as.list(params)))
}

#' Generalized Jacobian call
#'
#' General, shared method for constructing the Jacobian matrix. Routes a supplied "jacobianFunction" with pre-computed derivatives to construct matrix with observed data and supplied parameters.
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return difference value for jacobian
jacobianMatrix <- function(params, x, value, valueFunction, jacobianFunction)
{
  -do.call("jacobianFunction", c(list(x = x), as.list(params)))
}
