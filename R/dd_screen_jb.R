#' Perform Johnson & Bickel Screen
#'
#' This function applies the Johnson & Bickel screening criteria to included data series. The result of this procedure is a TRUE/FALSE response to one of two screening criteria.
#'
#' @param fittingObject core fitting object
#'
#' @return A data frame of model screenings
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return
#' @export
johnsonBickelScreen <- function(fittingObject) {

  listOfIds = unique(fittingObject$data[[as.character(fittingObject$settings['Individual'])]])

  for (id in listOfIds) {
    messageDebug(fittingObject, paste("JB Screen: ", id))

    currentData = fittingObject$data[
      which(fittingObject$data[,
                               as.character(fittingObject$settings['Individual'])] == id),]

    fittingObject$data[
      which(fittingObject$data[,
                               as.character(fittingObject$settings['Individual'])] == id), "JB1"] = TRUE

    fittingObject$data[
      which(fittingObject$data[,
                               as.character(fittingObject$settings['Individual'])] == id), "JB2"] = TRUE

    currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]
    currentData$ddY = currentData[,as.character(fittingObject$settings['Values'])]
    currentData$ddY = currentData$ddY / as.numeric(fittingObject[[ "maxValue" ]])

    currentData = currentData[order(currentData$ddX), ]

    for (index in 2:length(currentData$ddX)) {
      prev = currentData[index - 1, "ddY"]
      curr = currentData[index,     "ddY"]

      if ((curr - prev) > as.numeric(fittingObject[[ "JB1Flag" ]])) {
        messageDebug(fittingObject, paste("JB Screen: ", id, "[Fail JB1]"))

        fittingObject$data[
          which(fittingObject$data[,
                                   as.character(fittingObject$settings['Individual'])] == id), "JB1"] = FALSE
      }
    }

    prev <- currentData[1,                       "ddY"]
    curr <- currentData[length(currentData$ddX), "ddY"]

    if ((prev - curr) < as.numeric(fittingObject[[ "JB2Flag" ]])) {
      messageDebug(fittingObject, paste("JB Screen: ", id, "[Fail JB2]"))

      fittingObject$data[
        which(fittingObject$data[,
                                 as.character(fittingObject$settings['Individual'])] == id), "JB2"] = FALSE
    }
  }

  fittingObject
}
