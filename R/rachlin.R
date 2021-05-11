
#' dd_fit_rachlin
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_rachlin <- function(fittingObject, id) {

  modelResults = list(
    Model      = "rachlin",
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

  startParams = dd_start_rachlin(currentData)

  modelFitRachlin <- NULL

  try(modelFitRachlin <- nls.lm(par              = startParams,
                                fn               = residualFunction,
                                jac              = jacobianMatrix,
                                valueFunction    = rachlinHyperboloidDiscountFunc,
                                jacobianFunction = rachlinHyperboloidDiscountGradient,
                                x                = currentData$ddX,
                                value            = currentData$ddY,
                                control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitRachlin)) {

    modelResults[[ "Lnk"   ]]  = modelFitRachlin$par[["lnk"]]
    modelResults[[ "S"  ]]     = modelFitRachlin$par[["s"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitRachlin$deviance/length(modelFitRachlin$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitRachlin))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitRachlin))
    modelResults[[ "Status" ]] = paste("Code:", modelFitRachlin$info,
                                       "- Message:", modelFitRachlin$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][["rachlin"]] = modelResults

  fittingObject
}

#' dd_start_rachlin
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @return
#' @export
dd_start_rachlin <- function(currentData) {

  startlnK <- seq(-12, 12, 1)
  starts   <- seq(.01, 10, .01)

  lengthLnK <- length(startlnK)
  lengthX    = nrow(currentData)
  lengthS <- length(starts)

  SSlnK <- rep(startlnK, lengthS)
  SSs <- sort(rep(starts, lengthLnK))

  sumSquares <- rep(NA, lengthS * lengthLnK)

  SY <- rep(currentData$ddY, lengthS * lengthLnK)

  SlnK <- rep(sort(rep(startlnK,lengthX)), lengthS)
  Ss <- sort(rep(starts, lengthX * lengthLnK))

  projection <- (1 + exp(SlnK) * (currentData$ddX^Ss))^(-1)
  sqResidual <- (SY - projection)^2

  for (j in 1:(lengthS*lengthLnK)) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort <- data.frame(SSlnK, SSs, sumSquares)
  sorted  <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk = sorted$SSlnK[1],
               s   = sorted$SSs[1])

  ini.par
}

#' dd_ed50_rachlin
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_ed50_rachlin <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["rachlin"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["rachlin"]][["S"]]

  fittingObject$ed50[[as.character(id)]] = log( (1/(exp(lnk)))^(1/s))

  fittingObject
}

#' dd_mbauc_rachlin
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_mbauc_rachlin <- function(fittingObject, id) {

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  maxX        = max(currentData$ddX)
  minX        = min(currentData$ddX)
  maximumArea = maxX - minX

  lnk = fittingObject$results[[as.character(id)]][["rachlin"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["rachlin"]][["S"]]

  fittingObject$mbauc[[as.character(id)]] = stats::integrate(integrandRachlin,
                                                             lower = minX,
                                                             upper = maxX,
                                                             lnK   = lnk,
                                                             s     = s)$value/maximumArea

  fittingObject
}

#' dd_mbauc_log10_rachlin
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_mbauc_log10_rachlin <- function(fittingObject, id) {

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  maxX        = log10(max(currentData$ddX))
  minX        = log10(min(currentData$ddX))
  maximumArea = maxX - minX

  lnk = fittingObject$results[[as.character(id)]][["rachlin"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["rachlin"]][["S"]]

  fittingObject$mbauclog10[[as.character(id)]] = stats::integrate(integrandRachlinLog,
                                                                  lower = minX,
                                                                  upper = maxX,
                                                                  lnK   = lnk,
                                                                  s     = s)$value/maximumArea

  fittingObject
}

#' Rachlin Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
rachlinHyperboloidDiscountFunc <- function(x, lnk, s)
{
  func <- (1+exp(lnk)*(x^s))^(-1)
  eval(func)
}

#' Rachlin Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
rachlinHyperboloidDiscountGradient <- function(x, lnk, s)
{
  func <- expression((1+exp(lnk)*x)^(-s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Rachlin Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandRachlin <- function(x, lnK, s) { (1+exp(lnK)*(x^s))^(-1) }

#' Rachlin Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandRachlinLog <- function(x, lnK, s) { (1+exp(lnK)*((10^x)^s))^(-1) }
