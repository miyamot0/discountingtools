
#' dd_fit_exponential
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_exponential <- function(fittingObject, id) {

  modelResults = list(
    Model      = "exponential",
    Lnk        = NA,
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

  startParams = dd_start_exponential(currentData)

  modelFitExponential <- NULL

  try(modelFitExponential <- nls.lm(par              = startParams,
                                    fn               = residualFunction,
                                    jac              = jacobianMatrix,
                                    valueFunction    = exponentialDiscountFunc,
                                    jacobianFunction = exponentialDiscountGradient,
                                    x                = currentData$ddX,
                                    value            = currentData$ddY,
                                    control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitExponential)) {

    modelResults[[ "Lnk"    ]] = modelFitExponential$par[["lnk"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitExponential$deviance/length(modelFitExponential$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitExponential))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitExponential))
    modelResults[[ "Status" ]] = paste("Code:", modelFitExponential$info,
                                          "- Message:", modelFitExponential$message,
                                          sep = " ")
  }

  fittingObject$results[[as.character(id)]][["exponential"]] = modelResults

  fittingObject
}

#' dd_start_exponential
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#' @param increment step size for span
#'
#' @return
#' @export
dd_start_exponential <- function(currentData, increment = 1) {

  startlnK  <- seq(-15,  15, increment)
  lengthLnK  = length(startlnK)
  lengthX    = nrow(currentData)
  MlnK       = sort(rep(startlnK, lengthX))

  sumSquares = rep(NA,lengthLnK)
  MX         = rep(currentData$ddX, lengthLnK)
  MY         = rep(currentData$ddY, lengthLnK)

  # Projections
  projection = exp(-exp(MlnK) * MX)
  sqResidual = (MY - projection)^2

  for (j in 1:lengthLnK) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort    = data.frame(startlnK, sumSquares)
  sorted     = presort[order(presort[ ,"sumSquares"]), ]
  ini.par    = c(lnk = sorted$startlnK[1])

  ini.par
}

#' dd_ed50_exponential
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_ed50_exponential <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["exponential"]][["Lnk"]]

  fittingObject$ed50[[as.character(id)]] = log(log(2)/exp(lnk))

  fittingObject
}

#' dd_mbauc_exponential
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_mbauc_exponential <- function(fittingObject, id) {

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  maxX        = max(currentData$ddX)
  minX        = min(currentData$ddX)
  maximumArea = maxX - minX

  lnk = fittingObject$results[[as.character(id)]][["exponential"]][["Lnk"]]

  fittingObject$mbauc[[as.character(id)]] = stats::integrate(integrandExp,
                                                             lower = minX,
                                                             upper = maxX,
                                                             lnK = lnk)$value/maximumArea

  fittingObject
}

#' Exponential Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandExp <- function(x, lnK) { exp(-exp(lnK)*x) }

#' Exponential Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandExpLog <- function(x, lnK) { exp(-exp(lnK)*(10^x)) }
