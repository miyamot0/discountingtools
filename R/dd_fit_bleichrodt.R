
#' dd_fit_bleichrodt
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom minpack.lm nls.lm nls.lm.control
dd_fit_bleichrodt <- function(fittingObject, id) {

  modelResults = list(
    Model      = "bleichrodt",
    Lnk        = NA,
    S          = NA,
    Beta       = NA,
    RMSE       = NA,
    BIC        = NA,
    AIC        = NA,
    ED50       = NA,
    MBAUC      = NA,
    Log10MBAUC = NA
  )

  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]
  currentData$ddY = currentData[,as.character(fittingObject$settings['Values'])]
  currentData$ddY = currentData$ddY / as.numeric(fittingObject[[ "maxValue" ]])

  startParams = dd_start_bleichrodt(currentData)

  modelFitBleichrodt <- NULL

  try(modelFitBleichrodt <- nls.lm(par              = startParams,
                                   fn               = residualFunction,
                                   jac              = jacobianMatrix,
                                   valueFunction    = dd_discount_func_bleichrodt_crdi,
                                   jacobianFunction = dd_discount_grad_bleichrodt_crdi,
                                   x                = currentData$ddX,
                                   value            = currentData$ddY,
                                   upper            = c(beta = 1,
                                                        lnk = Inf,
                                                        s = Inf),
                                   lower            = c(beta = 0,
                                                        lnk = -Inf,
                                                        s = -Inf),
                                   control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitBleichrodt)) {

    modelResults[[ "Lnk"    ]] = modelFitBleichrodt$par[["lnk"]]
    modelResults[[ "S"      ]] = modelFitBleichrodt$par[["s"]]
    modelResults[[ "Beta"   ]] = modelFitBleichrodt$par[["beta"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitBleichrodt$deviance/length(modelFitBleichrodt$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitBleichrodt))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitBleichrodt))
    modelResults[[ "ED50"        ]] = dd_ed50_bleichrodt(
      Lnk = modelFitBleichrodt$par[["lnk"]],
      s   = modelFitBleichrodt$par[["s"]],
      b   = modelFitBleichrodt$par[["beta"]],
      currentData
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_bleichrodt(
      A = 1,
      Lnk = modelFitBleichrodt$par[["lnk"]],
      s   = modelFitBleichrodt$par[["s"]],
      b   = modelFitBleichrodt$par[["beta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_bleichrodt(
      A = 1,
      Lnk = modelFitBleichrodt$par[["lnk"]],
      s   = modelFitBleichrodt$par[["s"]],
      b   = modelFitBleichrodt$par[["beta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Status" ]] = paste("Code:", modelFitBleichrodt$info,
                                       "- Message:", modelFitBleichrodt$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][["bleichrodt"]] = modelResults

  fittingObject
}

#' dd_start_bleichrodt
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_start_bleichrodt <- function(currentData) {

  startlnK  <- seq(-12,  12, 1)
  starts    <- seq(0.01, 1,  0.1)
  startBeta <- seq(0.1,  1,  0.1)

  # new Pre sort
  presort <- expand.grid(startlnK  = startlnK,
                         starts    = starts,
                         startBeta = startBeta)

  presort$sumSquares <- NA

  # clean, merge, or move
  getSS <- function(presort, index, Y, X) {
    projections <- presort[index,]$startBeta*exp(-(exp(presort[index,]$startlnK)*X^presort[index,]$starts))
    sqResidual  <- (Y - projections)^2
    sum(sqResidual)
  }

  for (j in 1:nrow(presort)) {
    presort[j, ]$sumSquares <- getSS(presort, j, currentData$ddY, currentData$ddX)
  }

  presort <- presort[order(presort$sumSquares),]

  ini.par <- c(beta = presort[1,]$startBeta,
               lnk = presort[1,]$startlnK,
               s = presort[1,]$starts)

  ini.par
}

#' dd_ed50_bleichrodt
#'
#' @param Lnk parameter
#' @param s parameter
#' @param b parameter
#' @param currentData currentData
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_ed50_bleichrodt <- function(Lnk, s, b, currentData) {
  lowDelay <- 0
  highDelay <- max(currentData$ddX)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst  <- integrand_bleichrodt_crdi(  lowDelay, Lnk, s, b)
    midEst  <- integrand_bleichrodt_crdi( (lowDelay + highDelay)/2, Lnk, s, b)
    highEst <- integrand_bleichrodt_crdi(  highDelay, Lnk, s, b)

    if (lowEst > 0.5 && midEst > 0.5) {
      lowDelay <- (lowDelay + highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      lowDelay <- lowDelay
      highDelay <- (lowDelay + highDelay)/2

    }
  }

  return(log((lowDelay + highDelay)/2))
}

#' dd_mbauc_bleichrodt
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param b parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_bleichrodt <- function(A, Lnk, s, b, startDelay, endDelay) {
  maxX        = endDelay
  minX        = startDelay
  maximumArea = (maxX - minX) * A

  area = stats::integrate(integrand_bleichrodt_crdi,
                          lower = minX,
                          upper = maxX,
                          lnK   = Lnk,
                          s     = s,
                          beta  = b)$value/maximumArea

  return(area)
}

#' dd_mbauc_log10_bleichrodt
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param b parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_bleichrodt <- function(A, Lnk, s, b, startDelay, endDelay) {
  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = (maxX - minX) * A

  area = stats::integrate(integrand_bleichrodt_crdi_log10,
                          lower = minX,
                          upper = maxX,
                          lnK   = Lnk,
                          s     = s,
                          beta  = b)$value/maximumArea

  return(area)
}

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_bleichrodt_crdi <- function(x, lnk, s, beta)
{
  func <- beta * exp(-exp(lnk)*x^s)
  eval(func)
}

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_bleichrodt_crdi <- function(x, lnk, s, beta)
{
  func <- beta * exp(-exp(lnk)*x^s)
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")),
    eval(stats::deriv(func, "beta")))
}

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
integrand_bleichrodt_crdi <- function(x, lnK, s, beta) {  beta * exp(-exp(lnK)*x^s) }

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
integrand_bleichrodt_crdi_log10 <- function(x, lnK, s, beta) {  beta * exp(-exp(lnK)*(10^x)^s) }
