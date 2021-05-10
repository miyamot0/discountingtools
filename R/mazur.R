
#' Title
#'
#' @param fittingObject
#' @param id id tag
#'
#' @return
#' @export
dd_fit_mazur <- function(fittingObject, id) {

  modelResults = list(
    Model      = "mazur",
    Intercept  = NA,
    RMSE       = NA,
    BIC        = NA,
    AIC        = NA
  )

  # modelFitNoise <- NULL
  #
  # try(modelFitNoise <- stats::lm(Y ~ 1, fittingObject$data),
  #     silent = TRUE)
  #
  # if (!is.character(modelFitNoise)) {
  #
  #   modelResults[[ "Intercept" ]] = modelFitNoise$coefficients[["(Intercept)"]]
  #   modelResults[[ "RMSE"      ]] = summary(modelFitNoise)[["sigma"]]
  #   modelResults[[ "BIC"       ]] = ifelse(summary(modelFitNoise)[["sigma"]] == 0, Inf, stats::BIC(modelFitNoise))
  #   modelResults[[ "AIC"       ]] = ifelse(summary(modelFitNoise)[["sigma"]] == 0, Inf, stats::AIC(modelFitNoise))
  #
  # }
  #
  # fittingObject$results[[(length(fittingObject[["results"]]) + 1)]] = modelResults

  fittingObject
}

#' Title
#'
#' @return
#' @export
dd_start_mazur <- function() {
  # Starts
  startlnK <- seq(-12,  12, 1)
  lengthLnK <- length(startlnK)

  # Params expand
  MlnK <- sort(rep(startlnK, lengthX))

  # Init SS vector
  sumSquares <- rep(NA,lengthLnK)

  # Observed Data
  MX <- rep(dat$X, lengthLnK)
  MY <- rep(dat$Y, lengthLnK)

  # Projections
  projection <- (1 + exp(MlnK)*MX)^(-1)
  sqResidual <- (MY - projection)^2

  for (j in 1:lengthLnK){
    sumSquares[j] <- sum(sqResidual[(j - 1) * lengthX + 1:lengthX])

  }

  # Sort starting estimates
  presort <- data.frame(startlnK, sumSquares)
  sorted  <- presort[order(presort[ ,"sumSquares"]), ]
  ini.par <- c(lnk=sorted$startlnK[1])

  ini.par
}
