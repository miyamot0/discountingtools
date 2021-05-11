
#' dd_fit_rodriguezlogue
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_rodriguezlogue <- function(fittingObject, id) {

  modelResults = list(
    Model      = "rodriguezlogue",
    Lnk        = NA,
    Beta       = NA,
    RMSE       = NA,
    BIC        = NA,
    AIC        = NA
  )

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]
  currentData$ddY = currentData[,as.character(fittingObject$settings['Values'])]
  currentData$ddY = currentData$ddY / as.numeric(fittingObject[[ "maxValue" ]])

  startParams = dd_start_rodriguezlogue(currentData)

  modelFitRodriguezLogue <- NULL

  try(modelFitRodriguezLogue <- nls.lm(par              = startParams,
                                       fn               = residualFunction,
                                       jac              = jacobianMatrix,
                                       valueFunction    = RodriguezLogueDiscountFunc,
                                       jacobianFunction = RodriguezLogueDiscountGradient,
                                       x                = currentData$ddX,
                                       value            = currentData$ddY,
                                       control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitRodriguezLogue)) {

    modelResults[[ "Lnk"    ]] = modelFitRodriguezLogue$par[["lnk"]]
    modelResults[[ "Beta"   ]] = modelFitRodriguezLogue$par[["beta"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitRodriguezLogue$deviance/length(modelFitRodriguezLogue$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitRodriguezLogue))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitRodriguezLogue))
    modelResults[[ "Status" ]] = paste("Code:", modelFitRodriguezLogue$info,
                                       "- Message:", modelFitRodriguezLogue$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][["rodriguezlogue"]] = modelResults

  fittingObject
}

#' dd_start_rodriguezlogue
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @return
#' @export
dd_start_rodriguezlogue <- function(currentData) {

  startlnK <- seq(-12, 12, 1)
  startBeta <- seq(-12, 12, 1)

  presort <- expand.grid(startlnK  = startlnK,
                         startBeta = startBeta)
  presort$sumSquares <- NA

  getSS <- function(presort, index, Y, X) {
    projections <- (1 + X * exp(presort[index,]$startlnK))^(-exp(presort[index,]$startBeta) / exp(presort[index,]$startlnK))
    sqResidual <- (Y - projections)^2
    sum(sqResidual)
  }

  for (j in 1:nrow(presort)) {
    presort[j, ]$sumSquares <- getSS(presort, j, currentData$ddY, currentData$ddX)
  }

  presort <- presort[order(presort$sumSquares),]

  ini.par <- c(lnk  = presort[1,]$startlnK,
               beta = presort[1,]$startBeta)

  ini.par
}

#' dd_ed50_rodriguezlogue
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_ed50_rodriguezlogue <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["rodriguezlogue"]][["Lnk"]]
  b   = fittingObject$results[[as.character(id)]][["rodriguezlogue"]][["Beta"]]

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  lowDelay <- 0
  highDelay <- max(currentData$ddX)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst  <- integrandRodriguezLogue(  lowDelay, lnk, b)
    midEst  <- integrandRodriguezLogue( (lowDelay+highDelay)/2, lnk, b)
    highEst <- integrandRodriguezLogue(  highDelay, lnk, b)

    if (lowEst > 0.5 && midEst > 0.5) {
      lowDelay <- (lowDelay+highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      lowDelay <- lowDelay
      highDelay <- (lowDelay+highDelay)/2

    }
  }

  fittingObject$ed50[[as.character(id)]] = log((lowDelay+highDelay)/2)

  fittingObject
}
