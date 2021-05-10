
#' dd_fit_ebertprelec
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_ebertprelec <- function(fittingObject, id) {

  modelResults = list(
    Model      = "ebertprelec",
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

  startParams = dd_start_ebertprelec(currentData)

  modelFitEbertPrelec <- NULL

  try(modelFitEbertPrelec <- nls.lm(par              = startParams,
                                    fn               = residualFunction,
                                    jac              = jacobianMatrix,
                                    valueFunction    = ebertPrelecDiscountFunc,
                                    jacobianFunction = ebertPrelecDiscountGradient,
                                    x                = currentData$ddX,
                                    value            = currentData$ddY,
                                    control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitEbertPrelec)) {

    modelResults[[ "Lnk"   ]]  = modelFitEbertPrelec$par[["lnk"]]
    modelResults[[ "S"  ]]     = modelFitEbertPrelec$par[["s"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitEbertPrelec$deviance/length(modelFitEbertPrelec$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitEbertPrelec))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitEbertPrelec))
    modelResults[[ "Status" ]] = paste("Code:", modelFitEbertPrelec$info,
                                       "- Message:", modelFitEbertPrelec$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][[
    (length(fittingObject$results[[as.character(id)]]) + 1)]] = modelResults

  fittingObject
}

#' dd_start_ebertprelec
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @return
#' @export
dd_start_ebertprelec <- function(currentData) {

  startlnK <- seq(-12, 12, 0.1)
  starts   <- seq(.01, 1, 0.01)

  lengthLnK <- length(startlnK)
  lengthX    = nrow(currentData)
  lengthS <- length(starts)

  SSlnK <- rep(startlnK, lengthS)
  SSs <- sort(rep(starts, lengthLnK))

  sumSquares <- rep(NA, lengthS * lengthLnK)

  SY <- rep(currentData$ddY, lengthS * lengthLnK)

  SlnK <- rep(sort(rep(startlnK,lengthX)), lengthS)
  Ss <- sort(rep(starts, lengthX * lengthLnK))

  projection <- exp(-(exp(SlnK) * currentData$ddX)^Ss)
  sqResidual <- (SY - projection)^2

  for (j in 1:(lengthS*lengthLnK)) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort <- data.frame(SSlnK, SSs, sumSquares)
  sorted  <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk = sorted$SSlnK[1],
               s   = sorted$SSs[1])

  ini.par
}
