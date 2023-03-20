
#' dd_fit_exponential
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_exponential <- function(fittingObject, id) {

  modelResults = list(
    Model      = "exponential",
    Lnk        = NA,
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

  startParams = dd_start_exponential(currentData)

  modelFitExponential <- NULL

  try(modelFitExponential <- nls.lm(par              = startParams,
                                   fn               = residualFunction,
                                   jac              = jacobianMatrix,
                                   valueFunction    = dd_discount_func_exponential,
                                   jacobianFunction = dd_discount_grad_exponential,
                                   x                = currentData$ddX,
                                   value            = currentData$ddY,
                                   control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitExponential)) {

    modelResults[[ "Lnk"         ]] = modelFitExponential$par[["lnk"]]
    modelResults[[ "RMSE"        ]] = sqrt(modelFitExponential$deviance/length(modelFitExponential$fvec))
    modelResults[[ "BIC"         ]] = stats::BIC(logLik.nls.lm(modelFitExponential))
    modelResults[[ "AIC"         ]] = stats::AIC(logLik.nls.lm(modelFitExponential))
    modelResults[[ "ED50"        ]] = dd_ed50_exponential(modelFitExponential$par[["lnk"]])
    modelResults[[ "MBAUC"       ]] = dd_mbauc_exponential(
      A = 1,
      Lnk = modelFitExponential$par[["lnk"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_exponential(
      A = 1,
      Lnk = modelFitExponential$par[["lnk"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Status"      ]] = paste("Code:", modelFitExponential$info,
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
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
#' @param Lnk log transformed rate parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_ed50_exponential <- function(Lnk) {
  return(log(log(2)/exp(Lnk)))
}

#' dd_mbauc_exponential
#'
#' @param A maximum value of good
#' @param Lnk log transformed rate parameter
#' @param startDelay start delay
#' @param endDelay end delay
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_mbauc_exponential <- function(A, Lnk, startDelay, endDelay) {
  expFinal = (-A * exp(-exp(Lnk) * endDelay)) / exp(Lnk)
  expInitial = (-A * exp(-exp(Lnk) * startDelay)) / exp(Lnk)

  return((expFinal - expInitial) / ((endDelay - startDelay) * A))
}

#' dd_mbauc_log10_exponential
#'
#' @param A maximum value of good
#' @param Lnk log transformed rate parameter
#' @param startDelay start delay
#' @param endDelay end delay
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_exponential <- function(A, Lnk, startDelay, endDelay) {

  maximumArea = (endDelay - startDelay) * A

  area = stats::integrate(dd_integrand_exponential_log10,
                          lower = startDelay,
                          upper = endDelay,
                          lnK = Lnk)$value/maximumArea

  return(area)
}

#' Exponential discounting function
#'
#' @param x observation at point n (X)
#' @param Lnk fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_exponential <- function(x, Lnk)
{
  func <- exp(-exp(Lnk)*x)
  eval(func)
}

#' Exponential Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_exponential <- function(x, lnk)
{
  func <- expression(exp(-exp(lnk)*x))
  c(eval(stats::D(func, "lnk")))
}

#' Exponential Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_exponential_log10 <- function(x, lnK) { exp(-exp(lnK)*(10^x)) }
