
#' dd_fit_laibson
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_laibson <- function(fittingObject, id) {

  modelResults = list(
    Model      = "laibson",
    Beta       = NA,
    Delta      = NA,
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

  startParams = dd_start_laibson(currentData)

  modelFitLaibson <- NULL

  try(modelFitLaibson <- nls.lm(par              = startParams,
                                fn               = residualFunction,
                                jac              = jacobianMatrix,
                                valueFunction    = dd_discount_func_laibson,
                                jacobianFunction = dd_discount_grad_laibson,
                                x                = currentData$ddX,
                                value            = currentData$ddY,
                                upper            = c(beta = 1, delta = 1),
                                lower            = c(beta = 0, delta = 0),
                                control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitLaibson)) {

    modelResults[[ "Beta"        ]] = modelFitLaibson$par[["beta"]]
    modelResults[[ "Delta"       ]] = modelFitLaibson$par[["delta"]]
    modelResults[[ "RMSE"        ]] = sqrt(modelFitLaibson$deviance/length(modelFitLaibson$fvec))
    modelResults[[ "BIC"         ]] = stats::BIC(logLik.nls.lm(modelFitLaibson))
    modelResults[[ "AIC"         ]] = stats::AIC(logLik.nls.lm(modelFitLaibson))
    modelResults[[ "ED50"        ]] = dd_ed50_laibson(
      b = modelFitLaibson$par[["beta"]],
      d = modelFitLaibson$par[["delta"]]
    )
    modelResults[[ "MBAUC"       ]] = dd_mbauc_laibson(
      A = 1,
      b = modelFitLaibson$par[["beta"]],
      d = modelFitLaibson$par[["delta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_laibson(
      A = 1,
      b = modelFitLaibson$par[["beta"]],
      d = modelFitLaibson$par[["delta"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Status" ]] = paste("Code:", modelFitLaibson$info,
                                       "- Message:", modelFitLaibson$message,
                                       sep = " ")
  }

  fittingObject$results[[as.character(id)]][["laibson"]] = modelResults

  fittingObject
}

#' dd_start_laibson
#'
#' Extract starting parameters
#'
#' @param currentData  current data set
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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

#' dd_ed50_laibson
#'
#' @param b beta param
#' @param d delta param
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_ed50_laibson <- function(b, d) {
  return(log(log( (1/(2*b)),base = d)))
}

#' dd_mbauc_laibson
#'
#' @param A maximum value
#' @param b parameter value
#' @param d parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_mbauc_laibson <- function(A, b, d, startDelay, endDelay) {
  bdFinal = (-A * b * exp(-(1 - d) * endDelay)) / (1 - d)
  bdInitial = (-A * b * exp(-(1 - d) * startDelay)) / (1 - d)

  return((bdFinal - bdInitial) / ((endDelay - startDelay) * A))
}

#' dd_mbauc_log10_laibson
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_laibson <- function(A, b, d, startDelay, endDelay) {
  maxX        = log10(endDelay)
  minX        = log10(startDelay)
  maximumArea = (maxX - minX) * A

  fittingObject$mbauclog10[[as.character(id)]] = stats::integrate(integrandBetaDeltaLog,
                                                                  lower = minX,
                                                                  upper = maxX,
                                                                  beta  = b,
                                                                  delta = d)$value/maximumArea

  fittingObject
}

#' Beta Delta Value Function
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_laibson <- function(x, beta, delta)
{
  func <- beta*delta^x
  eval(func)
}

#' Beta Delta Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_laibson <- function(x, beta, delta)
{
  func <- expression(beta*delta^x)
  c(eval(stats::D(func, "delta")),
    eval(stats::D(func, "beta")))
}

#' Beta Delta Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_laibson_log10 <- function(x, beta, delta) { beta*delta^(10^x) }
