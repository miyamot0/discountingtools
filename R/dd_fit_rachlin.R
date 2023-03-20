
#' dd_fit_rachlin
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_rachlin <- function(fittingObject, id) {

  modelResults = list(
    Model      = "rachlin",
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

  startParams = dd_start_rachlin(currentData)

  modelFitRachlin <- NULL

  try(modelFitRachlin <- nls.lm(par              = startParams,
                                fn               = residualFunction,
                                jac              = jacobianMatrix,
                                valueFunction    = dd_discount_func_rachlin,
                                jacobianFunction = dd_discount_grad_rachlin,
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
    modelResults[[ "ED50"        ]] = dd_ed50_rachlin(
      Lnk = modelFitRachlin$par[["lnk"]],
      s   = modelFitRachlin$par[["s"]]
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_rachlin(
      A = 1,
      Lnk = modelFitRachlin$par[["lnk"]],
      s   = modelFitRachlin$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_rachlin(
      A = 1,
      Lnk = modelFitRachlin$par[["lnk"]],
      s   = modelFitRachlin$par[["s"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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
#' @param Lnk parameter
#' @param s parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_ed50_rachlin <- function(Lnk, s) {
  return(log( (1/(exp(Lnk)))^(1/s)))
}

#' dd_mbauc_laibson
#'
#' @param A maximum value
#' @param Lnk parameter value
#' @param s parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_mbauc_rachlin <- function(A, Lnk, s, startDelay, endDelay) {
  rachFinal   <- A * endDelay * gauss_2F1((1.0), (1.0/s), (1 + (1.0/s)), (-exp(Lnk) * (endDelay )^s))
  rachInitial <- A * startDelay * gauss_2F1((1.0), (1.0/s), (1 + (1.0/s)), (-exp(Lnk) * (startDelay )^s))

  return((rachFinal - rachInitial) / ((endDelay - startDelay) * A))
}

#' dd_mbauc_log10_rachlin
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_rachlin <- function(A, Lnk, s, startDelay, endDelay) {

  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = (maxX - minX) * A

  fittingObject$mbauclog10[[as.character(id)]] = stats::integrate(dd_integrand_rachlin_log10,
                                                                  lower = minX,
                                                                  upper = maxX,
                                                                  lnK   = Lnk,
                                                                  s     = s)$value/maximumArea

  fittingObject
}

#' Rachlin Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_rachlin <- function(x, lnk, s)
{
  func <- (1 + exp(lnk)*(x^s))^(-1)
  eval(func)
}

#' Rachlin Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_rachlin <- function(x, lnk, s)
{
  func <- expression((1 + exp(lnk)*x)^(-s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Rachlin Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_rachlin_log10 <- function(x, lnK, s) { (1 + exp(lnK)*((10^x)^s))^(-1) }
