
#' dd_fit_greenmyerson
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_greenmyerson <- function(fittingObject, id) {

  modelResults = list(
    Model      = "greenmyerson",
    Lnk        = NA,
    S          = NA,
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

  startParams = dd_start_greenmyerson(currentData)

  modelFitGreenMyerson <- NULL

  try(modelFitGreenMyerson <- nls.lm(par              = startParams,
                                     fn               = residualFunction,
                                     jac              = jacobianMatrix,
                                     valueFunction    = myersonHyperboloidDiscountFunc,
                                     jacobianFunction = myersonHyperboloidDiscountGradient,
                                     x                = currentData$ddX,
                                     value            = currentData$ddY,
                                     control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitGreenMyerson)) {

    modelResults[[ "Lnk"   ]]  = modelFitGreenMyerson$par[["lnk"]]
    modelResults[[ "S"  ]]     = modelFitGreenMyerson$par[["s"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitGreenMyerson$deviance/length(modelFitGreenMyerson$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitGreenMyerson))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitGreenMyerson))
    modelResults[[ "Status" ]] = paste("Code:", modelFitGreenMyerson$info,
                                       "- Message:", modelFitGreenMyerson$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][["greenmyerson"]] = modelResults

  fittingObject
}

#' dd_start_laibson
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @return
#' @export
dd_start_greenmyerson <- function(currentData) {

  startlnK <- seq(-12, 12, 1)
  starts <- seq(.01, 10, 0.01)

  lengthLnK <- length(startlnK)
  lengthX    = nrow(currentData)
  lengthS <- length(starts)

  SSlnK <- rep(startlnK, lengthS)
  SSs <- sort(rep(starts, lengthLnK))

  sumSquares <- rep(NA, lengthS * lengthLnK)

  SY <- rep(currentData$ddY, lengthS * lengthLnK)

  SlnK <- rep(sort(rep(startlnK,lengthX)), lengthS)
  Ss <- sort(rep(starts, lengthX * lengthLnK))

  projection <- (1 + exp(SlnK) * currentData$ddX)^(-Ss)
  sqResidual <- (SY - projection)^2

  for (j in 1:(lengthS*lengthLnK)) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort <- data.frame(SSlnK, SSs, sumSquares)
  sorted  <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk = sorted$SSlnK[1],
               s   = sorted$SSs[1])

  ini.par
}

#' dd_ed50_greenmyerson
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_ed50_greenmyerson <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["greenmyerson"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["greenmyerson"]][["S"]]

  fittingObject$ed50[[as.character(id)]] = log( (2^(1/s)-1)/exp(lnk))

  fittingObject
}
