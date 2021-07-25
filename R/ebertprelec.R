
#' dd_fit_ebertprelec
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
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

  fittingObject$results[[as.character(id)]][["ebertprelec"]] = modelResults

  fittingObject
}

#' dd_start_ebertprelec
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @return
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

#' dd_ed50_ebertprelec
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
dd_ed50_ebertprelec <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["ebertprelec"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["ebertprelec"]][["S"]]

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  lowDelay <- 0
  highDelay <- max(currentData$ddX)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst  <- integrandEbertPrelec(  lowDelay, lnk, s)
    midEst  <- integrandEbertPrelec( (lowDelay+highDelay)/2, lnk, s)
    highEst <- integrandEbertPrelec(  highDelay, lnk, s)

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

#' dd_mbauc_ebertprelec
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
dd_mbauc_ebertprelec <- function(fittingObject, id) {

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  maxX        = max(currentData$ddX)
  minX        = min(currentData$ddX)
  maximumArea = maxX - minX

  lnk = fittingObject$results[[as.character(id)]][["ebertprelec"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["ebertprelec"]][["S"]]

  fittingObject$mbauc[[as.character(id)]] = stats::integrate(integrandEbertPrelec,
                                                             lower = minX,
                                                             upper = maxX,
                                                             lnK   = lnk,
                                                             s     = s)$value/maximumArea

  fittingObject
}

#' dd_mbauc_log10_ebertprelec
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
dd_mbauc_log10_ebertprelec <- function(fittingObject, id) {

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  maxX        = log10(max(currentData$ddX))
  minX        = log10(min(currentData$ddX))
  maximumArea = maxX - minX

  lnk = fittingObject$results[[as.character(id)]][["ebertprelec"]][["Lnk"]]
  s   = fittingObject$results[[as.character(id)]][["ebertprelec"]][["S"]]

  fittingObject$mbauclog10[[as.character(id)]] = stats::integrate(integrandEbertPrelecLog,
                                                                  lower = minX,
                                                                  upper = maxX,
                                                                  lnK   = lnk,
                                                                  s     = s)$value/maximumArea

  fittingObject
}

#' Ebert & Prelec Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
#' @export
ebertPrelecDiscountFunc <- function(x, lnk, s)
{
  func <- exp(-(exp(lnk)*x)^s)
  eval(func)
}

#' Ebert & Prelec Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
ebertPrelecDiscountGradient <- function(x, lnk, s)
{
  func <- expression(exp(-(exp(lnk)*x)^s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Ebert & Prelec's Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandEbertPrelec <- function(x, lnK, s) {  exp(-(exp(lnK)*x)^s) }

#' Ebert & Prelec's ep Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
integrandEbertPrelecLog <- function(x, lnK, s) {  exp(-(exp(lnK)*(10^x))^s) }
