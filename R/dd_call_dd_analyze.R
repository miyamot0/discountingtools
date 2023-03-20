
#' dd_analyze
#'
#' This call is the workhorse of the program. Based on the settings applied, this method applies all relevant methods and calculations to the supplied data.
#'
#' @param fittingObject core dd fitting object
#' @param modelSelection (bool) this flag determines whether or not a model selection procedure will be applied in the results frame.
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_analyze <- function(fittingObject, modelSelection = FALSE) {

  if (length(fittingObject[["models"]]) < 1) {
    stop("No model options were selected")
  }

  ## TODO: checks for screening call
  if ("screen" %in% names(fittingObject) | "filterPassing" %in% names(fittingObject)) {
    message_debug(fittingObject, "Beginning JB Screening")

    fittingObject = johnsonBickelScreen(fittingObject)
  }

  ## TODO: Check for filter passing
  if ("filterPassing" %in% names(fittingObject)) {
    message_debug(fittingObject, "Filtering based on JB Screen")

    if ("JB1" %in% fittingObject[[ "filterPassing" ]])
      fittingObject$data = fittingObject$data[fittingObject$data$JB1 == TRUE, ]

    if ("JB2" %in% fittingObject[[ "filterPassing" ]])
      fittingObject$data = fittingObject$data[fittingObject$data$JB2 == TRUE, ]
  }

  ## TODO: Check for analytical strategy
  if (fittingObject[[ "strategy" ]] == "group") {
    if (is.null(fittingObject$settings[["Group"]])) stop('No Group aesthetic specified')

    message_debug(fittingObject, "Casting Individual Ids to Group Ids")

    vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])
    newGrpIds = match(fittingObject$data[,as.character(fittingObject$settings['Group'])], vecGroups)
    fittingObject$data[,as.character(fittingObject$settings['Individual'])] <- newGrpIds
  }

  message_debug(fittingObject, "Beginning Model Fitting(s)")

  fittingObject[[ "ModelSelection" ]] = modelSelection

  # Add in noise model as a comparator
  if (!("noise" %in% fittingObject[["models"]]) & modelSelection == TRUE)
    fittingObject[["models"]] = c("noise", fittingObject[["models"]])

  # loop through individual id's
  for (id in unique(fittingObject$data[[as.character(fittingObject$settings['Individual'])]])) {
    message_debug(fittingObject, paste("Fitting:", id))

    fittingObject$results[[as.character(id)]] = list()

    for (model in fittingObject[["models"]]) {
      message_debug(fittingObject, paste("Fitting:", id, "Rotation:", model))

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

    if (modelSelection == TRUE) {
      fittingObject = dd_probable_model(fittingObject, id)

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   fittingObject = get_ED50(fittingObject, id)
        if (metric == "mbauc")    fittingObject = get_MBAUC(fittingObject, id)
        if (metric == "logmbauc") fittingObject = get_MBAUC_log10(fittingObject, id)
      }
    }
  }

  fittingObject
}
