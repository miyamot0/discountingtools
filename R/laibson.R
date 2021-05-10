
#' dd_fit_laibson
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_laibson <- function(fittingObject, id) {

  modelResults = list(
    Model      = "laibson",
    Beta       = NA,
    Delta      = NA,
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

  startParams = dd_start_laibson(currentData)

  modelFitLaibson <- NULL

  try(modelFitLaibson <- nls.lm(par              = startParams,
                                fn               = residualFunction,
                                jac              = jacobianMatrix,
                                valueFunction    = betaDeltaDiscountFunc,
                                jacobianFunction = betaDeltaDiscountGradient,
                                x                = currentData$ddX,
                                value            = currentData$ddY,
                                upper            = c(beta = 1, delta = 1),
                                lower            = c(beta = 0, delta = 0),
                                control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitLaibson)) {

    modelResults[[ "Beta"   ]] = modelFitLaibson$par[["beta"]]
    modelResults[[ "Delta"  ]] = modelFitLaibson$par[["delta"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitLaibson$deviance/length(modelFitLaibson$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitLaibson))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitLaibson))
    modelResults[[ "Status" ]] = paste("Code:", modelFitLaibson$info,
                                       "- Message:", modelFitLaibson$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][[
    (length(fittingObject$results[[as.character(id)]]) + 1)]] = modelResults

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
dd_start_laibson <- function(currentData) {

  startbeta   <- seq(0, 1, 0.1)
  startdelta  <- seq(0, 1, 0.01)

  lengthX    = nrow(currentData)
  lengthBeta  <- length(startbeta)
  lengthDelta <- length(startdelta)

  startBeta   <- rep(sort(rep(startbeta, lengthX)), lengthDelta)
  startdelta  <- sort(rep(startdelta, lengthX * lengthBeta))

  sumSquares  <- rep(NA, lengthBeta * lengthDelta)

  SY          <- rep(currentData$ddY, lengthBeta * lengthDelta)
  projection  <- startBeta * startdelta^currentData$ddX

  sqResidual  <- (SY - projection)^2

  SSbeta      <- rep(startbeta, lengthDelta)
  SSdelta     <- sort(rep(startdelta, lengthBeta))

  for (j in 1:(lengthBeta * lengthDelta)) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort <- data.frame(SSbeta,
                        SSdelta,
                        sumSquares)

  sorted  <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(beta  = sorted$SSbeta[1],
               delta = sorted$SSdelta[1])

  ini.par
}
