
#' dd_fit_mazur
#'
#' This fits a hyperbolic model to the data.
#'
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_fit_mazur <- function(fittingObject, id) {

  modelResults = list(
    Model      = "mazur",
    Lnk        = NA,
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

  startParams = dd_start_mazur(currentData)

  modelFitHyperbolic <- NULL

  try(modelFitHyperbolic <- nls.lm(par              = startParams,
                                   fn               = residualFunction,
                                   jac              = jacobianMatrix,
                                   valueFunction    = hyperbolicDiscountFunc,
                                   jacobianFunction = hyperbolicDiscountGradient,
                                   x                = currentData$ddX,
                                   value            = currentData$ddY,
                                   control          = nls.lm.control(maxiter = 1000)),
      silent = TRUE)

  if (!is.null(modelFitHyperbolic)) {

    modelResults[[ "Lnk"    ]] = modelFitHyperbolic$par[["lnk"]]
    modelResults[[ "RMSE"   ]] = sqrt(modelFitHyperbolic$deviance/length(modelFitHyperbolic$fvec))
    modelResults[[ "BIC"    ]] = stats::BIC(logLik.nls.lm(modelFitHyperbolic))
    modelResults[[ "AIC"    ]] = stats::AIC(logLik.nls.lm(modelFitHyperbolic))
    modelResults[[ "Status" ]] = paste("Code:", modelFitHyperbolic$info,
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
#' @return
#' @export
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
#' @param fittingObject core dd fitting object
#' @param id id tag
#'
#' @return
#' @export
dd_ed50_mazur <- function(fittingObject, id) {

  lnk = fittingObject$results[[as.character(id)]][["mazur"]][["Lnk"]]

  fittingObject$ed50[[as.character(id)]] = log(1/(exp(lnk)))

  fittingObject
}
