
#' dd_fit_ebertprelec
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom minpack.lm nls.lm nls.lm.control
dd_fit_ebertprelec <- function(fittingObject, id) {

  modelResults = list(
    Model      = "ebertprelec",
    Lnk        = NA,
    S          = NA,
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

  startParams = dd_start_ebertprelec(currentData)

  modelFitEbertPrelec <- NULL

  try(modelFitEbertPrelec <- nls.lm(par              = startParams,
                                    fn               = residualFunction,
                                    jac              = jacobianMatrix,
                                    valueFunction    = dd_discount_func_ebertprelec,
                                    jacobianFunction = dd_grad_func_ebertprelec,
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
    modelResults[[ "ED50"        ]] = dd_ed50_ebertprelec(
      Lnk = modelFitEbertPrelec$par[["lnk"]],
      s   = modelFitEbertPrelec$par[["s"]],
      currentData
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_ebertprelec(
      A = 1,
      Lnk = modelFitEbertPrelec$par[["lnk"]],
      s   = modelFitEbertPrelec$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_ebertprelec(
      A = 1,
      Lnk = modelFitEbertPrelec$par[["lnk"]],
      s   = modelFitEbertPrelec$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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
#' @param Lnk parameter
#' @param s parameter
#' @param currentData currentData
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_ed50_ebertprelec <- function(Lnk, s, currentData) {
  lowDelay <- 0
  highDelay <- max(currentData$ddX)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst  <- integrand_ebertprelec(  lowDelay, Lnk, s)
    midEst  <- integrand_ebertprelec( (lowDelay + highDelay)/2, Lnk, s)
    highEst <- integrand_ebertprelec(  highDelay, Lnk, s)

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

#' dd_mbauc_ebertprelec
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_ebertprelec <- function(A, Lnk, s, startDelay, endDelay) {
  maxX        = endDelay
  minX        = startDelay
  maximumArea = (maxX - minX) * A

  area = stats::integrate(integrand_ebertprelec,
                          lower = minX,
                          upper = maxX,
                          lnK   = Lnk,
                          s     = s)$value/maximumArea

  return(area)
}

#' dd_mbauc_log10_ebertprelec
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_ebertprelec <- function(A, Lnk, s, startDelay, endDelay) {
  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = (maxX - minX) * A

  area = stats::integrate(integrand_ebertprelec_log10,
                          lower = minX,
                          upper = maxX,
                          lnK   = Lnk,
                          s     = s)$value/maximumArea

  return(area)
}

#' Ebert & Prelec Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_ebertprelec <- function(x, lnk, s)
{
  func <- exp(-(exp(lnk)*x)^s)
  eval(func)
}

#' Ebert & Prelec Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_grad_func_ebertprelec <- function(x, lnk, s)
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
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
integrand_ebertprelec <- function(x, lnK, s) {  exp(-(exp(lnK)*x)^s) }

#' Ebert & Prelec's ep Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
integrand_ebertprelec_log10 <- function(x, lnK, s) {  exp(-(exp(lnK)*(10^x))^s) }
