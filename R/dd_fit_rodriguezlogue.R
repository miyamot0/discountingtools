
#' dd_fit_rodriguezlogue
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_rodriguezlogue <- function(fittingObject, id) {

  modelResults = list(
    Model      = "rodriguezlogue",
    Lnk        = NA,
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

  startParams = dd_start_rodriguezlogue(currentData)

  modelFitRodriguezLogue <- NULL

  try(modelFitRodriguezLogue <- nls.lm(par              = startParams,
                                       fn               = residualFunction,
                                       jac              = jacobianMatrix,
                                       valueFunction    = dd_discount_func_rodriguezlogue,
                                       jacobianFunction = dd_discount_grad_rodriguezlogue,
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
    modelResults[[ "ED50"        ]] = dd_ed50_rodriguezlogue(
      Lnk = modelFitRodriguezLogue$par[["lnk"]],
      b   = modelFitRodriguezLogue$par[["beta"]]
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_rodriguezlogue(
      A = 1,
      Lnk = modelFitRodriguezLogue$par[["lnk"]],
      b   = modelFitRodriguezLogue$par[["beta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_rodriguezlogue(
      A = 1,
      Lnk = modelFitRodriguezLogue$par[["lnk"]],
      b   = modelFitRodriguezLogue$par[["beta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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
#' @param Lnk parameter
#' @param b parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_ed50_rodriguezlogue <- function(Lnk, b) {
  currentData = fittingObject$data[
    which(fittingObject$data[,
                             as.character(fittingObject$settings['Individual'])] == id),]

  currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]

  lowDelay <- 0
  highDelay <- max(currentData$ddX)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst  <- dd_integrand_rodriguezlogue(  lowDelay, Lnk, b)
    midEst  <- dd_integrand_rodriguezlogue( (lowDelay + highDelay)/2, Lnk, b)
    highEst <- dd_integrand_rodriguezlogue(  highDelay, Lnk, b)

    if (lowEst > 0.5 && midEst > 0.5) {
      lowDelay <- (lowDelay + highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      lowDelay <- lowDelay
      highDelay <- (lowDelay + highDelay) / 2

    }
  }

  return(log((lowDelay + highDelay)/2))
}

#' dd_mbauc_rodriguezlogue
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param b parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_rodriguezlogue <- function(A, Lnk, b, startDelay, endDelay) {
  maxX        = endDelay
  minX        = startDelay
  maximumArea = (maxX - minX) * A

  fittingObject$mbauc[[as.character(id)]] = stats::integrate(dd_integrand_rodriguezlogue,
                                                             lower = minX,
                                                             upper = maxX,
                                                             lnK   = Lnk,
                                                             beta  = b)$value/maximumArea

  fittingObject
}

#' dd_mbauc_log10_rodriguezlogue
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param b parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_rodriguezlogue <- function(A, Lnk, b, startDelay, endDelay) {
  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = (maxX - minX) * A

  fittingObject$mbauclog10[[as.character(id)]] = stats::integrate(dd_integrand_rodriguezlogue_log10,
                                                                  lower = minX,
                                                                  upper = maxX,
                                                                  lnK   = Lnk,
                                                                  beta  = b)$value/maximumArea

  fittingObject
}

#' Rodriguez & Logue Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_rodriguezlogue <- function(x, lnk, beta)
{
  func <- (1 + x * exp(lnk))^(-exp(beta) / exp(lnk))
  eval(func)
}

#' Rodriguez & Logue Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_rodriguezlogue <- function(x, lnk, beta)
{
  func <- expression((1 + x * exp(lnk))^(-exp(beta) / exp(lnk)))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "beta")))
}

#' Rodriguez & Logue Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_rodriguezlogue <- function(x, lnK, beta) { (1 + x * exp(lnK))^(-exp(beta) / exp(lnK)) }

#' Rodriguez & Logue Integrand helper
#'
#' This integrand helper is a projection of the integrand (log10) with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_rodriguezlogue_log10 <- function(x, lnK, beta) { (1 + (10^x) * exp(lnK))^(-exp(beta) / exp(lnK)) }
