
#' dd_fit_mazur
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_fit_mazur <- function(fittingObject, id) {

  modelResults = list(
    Model      = "mazur",
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

  startParams = dd_start_mazur(currentData)

  modelFitHyperbolic <- NULL

  try(modelFitHyperbolic <- nls.lm(par              = startParams,
                                   fn               = residualFunction,
                                   jac              = jacobianMatrix,
                                   valueFunction    = dd_discount_func_mazur,
                                   jacobianFunction = dd_discount_grad_mazur,
                                   x                = currentData$ddX,
                                   value            = currentData$ddY,
                                   control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitHyperbolic)) {

    modelResults[[ "Lnk"         ]] = modelFitHyperbolic$par[["lnk"]]
    modelResults[[ "RMSE"        ]] = sqrt(modelFitHyperbolic$deviance/length(modelFitHyperbolic$fvec))
    modelResults[[ "BIC"         ]] = stats::BIC(logLik.nls.lm(modelFitHyperbolic))
    modelResults[[ "AIC"         ]] = stats::AIC(logLik.nls.lm(modelFitHyperbolic))
    modelResults[[ "ED50"        ]] = dd_ed50_mazur(modelFitHyperbolic$par[["lnk"]])
    modelResults[[ "MBAUC"       ]] = dd_mbauc_mazur(
      A = 1,
      Lnk = modelFitHyperbolic$par[["lnk"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Log10MBAUC"  ]] = dd_mbauc_log10_mazur(
      A = 1,
      Lnk = modelFitHyperbolic$par[["lnk"]],
      startDelay = min(currentData$ddX),
      endDelay = max(currentData$ddX)
    )
    modelResults[[ "Status"      ]] = paste("Code:", modelFitHyperbolic$info,
                                          "- Message:", modelFitHyperbolic$message,
                                          sep = " ")
  }

  fittingObject$results[[as.character(id)]][["mazur"]] = modelResults

  fittingObject
}

#' dd_start_mazur
#'
#' Extract starting parameters
#'
#' @param currentData current data set
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_start_mazur <- function(currentData) {

  startlnK   = seq(-12,  12, 1)
  lengthLnK  = length(startlnK)
  lengthX    = nrow(currentData)
  MlnK       = sort(rep(startlnK, lengthX))

  sumSquares = rep(NA,lengthLnK)
  MX         = rep(currentData$ddX, lengthLnK)
  MY         = rep(currentData$ddY, lengthLnK)

  projection = (1 + exp(MlnK)*MX)^(-1)
  sqResidual = (MY - projection)^2

  for (j in 1:lengthLnK) sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  presort    = data.frame(startlnK, sumSquares)
  sorted     = presort[order(presort[ ,"sumSquares"]), ]
  ini.par    = c(lnk = sorted$startlnK[1])

  ini.par
}

#' dd_ed50_mazur
#'
#' @param Lnk log transformed rate parameter
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_ed50_mazur <- function(Lnk) {
  return(log(1/(exp(Lnk))))
}

#' dd_mbauc_mazur
#'
#' @param A maximum value
#' @param Lnk logged parameter value
#' @param startDelay time point
#' @param endDelay time point
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_mbauc_mazur <- function(A, Lnk, startDelay, endDelay) {
  hypFinal = (A * log((exp(Lnk) * endDelay) + 1)) / exp(Lnk)
  hypInitial = (A * log((exp(Lnk) * startDelay) + 1)) / exp(Lnk)

  return((hypFinal - hypInitial) / ((endDelay - startDelay) * A))
}

#' dd_mbauc_log10_mazur
#'
#' @param A maximum value of good
#' @param Lnk log transformed rate parameter
#' @param startDelay start delay
#' @param endDelay end delay
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_mbauc_log10_mazur <- function(A, Lnk, startDelay, endDelay) {

  maximumArea = (endDelay - startDelay) * A

  area = stats::integrate(dd_integrand_mazur_log10,
                          lower = startDelay,
                          upper = endDelay,
                          lnK = Lnk)$value/maximumArea

  return(area)
}

#' Hyperbolic Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export
dd_discount_func_mazur <- function(x, lnk)
{
  func <- (1 + exp(lnk)*x)^(-1)
  eval(func)
}

#' Hyperbolic Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#'
#' @return projected, subjective value
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_discount_grad_mazur <- function(x, lnk)
{
  func <- expression((1 + exp(lnk)*x)^(-1))
  c(eval(stats::D(func, "lnk")))
}

#' Hyperbolic Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#'
#' @return Numerical Integration Projection
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
dd_integrand_mazur_log10 <- function(x, lnK) { (1 + exp(lnK)*(10^x))^(-1) }
