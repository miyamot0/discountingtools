
#' dd_fit_greenmyerson
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_greenmyerson <- function(fittingObject, id) {

  modelResults = list(
    Model      = "greenmyerson",
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

  startParams = dd_start_greenmyerson(currentData)

  modelFitGreenMyerson <- NULL

  try(modelFitGreenMyerson <- nls.lm(par              = startParams,
                                     fn               = residualFunction,
                                     jac              = jacobianMatrix,
                                     valueFunction    = dd_discount_func_greenmyerson,
                                     jacobianFunction = dd_discount_grad_greenmyerson,
                                     x                = currentData$ddX,
                                     value            = currentData$ddY,
                                     control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitGreenMyerson)) {

    modelResults[[ "Lnk"         ]] = modelFitGreenMyerson$par[["lnk"]]
    modelResults[[ "S"           ]] = modelFitGreenMyerson$par[["s"]]
    modelResults[[ "RMSE"        ]] = sqrt(modelFitGreenMyerson$deviance/length(modelFitGreenMyerson$fvec))
    modelResults[[ "BIC"         ]] = stats::BIC(logLik.nls.lm(modelFitGreenMyerson))
    modelResults[[ "AIC"         ]] = stats::AIC(logLik.nls.lm(modelFitGreenMyerson))
    modelResults[[ "ED50"        ]] = dd_ed50_greenmyerson(
      Lnk = modelFitGreenMyerson$par[["lnk"]],
      s   = modelFitGreenMyerson$par[["s"]]
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_greenmyerson(
      A = 1,
      Lnk = modelFitGreenMyerson$par[["lnk"]],
      s   = modelFitGreenMyerson$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_greenmyerson(
      A = 1,
      Lnk = modelFitGreenMyerson$par[["lnk"]],
      s   = modelFitGreenMyerson$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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
#' @param Lnk parameter
#' @param s parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_ed50_greenmyerson <- function(Lnk, s) {
  return(log( (2^(1/s) - 1)/exp(Lnk)))
}

#' dd_mbauc_rachlin
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_mbauc_greenmyerson <- function(A, Lnk, s, startDelay, endDelay) {
  mgFinal = (A * ((exp(Lnk) * endDelay + 1)^(1 - s))) / (exp(Lnk) * (1 - s))
  mgInitial = (A * ((exp(Lnk) * startDelay + 1)^(1 - s))) / (exp(Lnk) * (1 - s))

  return((mgFinal - mgInitial) / ((endDelay - startDelay) * A))
}

#' dd_mbauc_log10_greenmyerson
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_greenmyerson <- function(A, Lnk, s, startDelay, endDelay) {

  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = maxX - minX

  area = stats::integrate(dd_integrand_myersongreen_log10,
                          lower = minX,
                          upper = maxX,
                          lnK   = Lnk,
                          s     = s)$value/maximumArea

  return(area)
}

#' Green & Myerson Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_greenmyerson <- function(x, lnk, s)
{
  func <- (1 + exp(lnk)*x)^(-s)
  eval(func)
}

#' Green & Myerson Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_greenmyerson <- function(x, lnk, s)
{
  func <- expression((1 + exp(lnk)*x)^(-s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Green & Myerson Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_myersongreen_log10 <- function(x, lnK, s) { (1 + exp(lnK)*(10^x))^(-s) }
