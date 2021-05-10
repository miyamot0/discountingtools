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
#'
#' @return
#' @export
fitDDCurves <- function(data, settings, maxValue) {

  fittingObject = list()                             # Primary object
  fittingObject[[ "settings" ]] = enexpr(settings)   #
  fittingObject[[ "data"     ]] = data               # Stored data
  fittingObject[[ "models"   ]] = character(0)       #
  fittingObject[[ "strategy" ]] = "ind"              #
  fittingObject[[ "metrics"  ]] = character(0)       #
  fittingObject[[ "results"  ]] = list()             #
  fittingObject[[ "maxValue" ]] = maxValue           #

  class(fittingObject) <- c("discountingtools")

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
dd_analyze <- function(fittingObject) {

  # Add in noise model as a comparator
  if (!("noise" %in% fittingObject[["models"]]))
    fittingObject[["models"]] = c("noise", fittingObject[["models"]])

  # loop through individual id's
  for (id in unique(fittingObject$data[[as.character(fittingObject$settings['Individual'])]])) {

    # add in keyed list of results
    fittingObject$results[[as.character(id)]] = list()

    for (model in fittingObject[["models"]]) {

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
  }

  fittingObject
}
