

#' Scoring for the most probable model area, in log10 space
#'
#' In this set of methods, the area beneath the fitted model is calculated and divided by the maximum area using numerical integration methods.  All delays are calculated in the log base 10 scale.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return area beneath the fitted model, in log10 space
getModelAUCLog10Scaled <- function(dat, results) {

  maximumArea <- log10(max(dat$X)) - log10(min(dat$X))

  returnValue <- NaN

  if (results[["probable.model"]] == "Hyperbolic") {
    returnValue <- stats::integrate(integrandHypLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["Hyperbolic.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Exponential") {
    returnValue <- stats::integrate(integrandExpLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["Exponential.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Laibson") {
    returnValue <- stats::integrate(integrandBetaDeltaLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    beta = results[["Laibson.beta"]],
                                    delta = results[["Laibson.delta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "GreenMyerson") {
    returnValue <- stats::integrate(integrandMyersonLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["GreenMyerson.lnk"]],
                                    s = results[["GreenMyerson.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Rachlin") {
    returnValue <- stats::integrate(integrandRachlinLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["Rachlin.lnk"]],
                                    s = results[["Rachlin.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "EbertPrelec") {
    returnValue <- stats::integrate(integrandEbertPrelecLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["EbertPrelec.lnk"]],
                                    s = results[["EbertPrelec.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    returnValue <- stats::integrate(integrandBleichrodtCRDILog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["BleichrodtCRDI.lnk"]],
                                    s = results[["BleichrodtCRDI.s"]],
                                    beta = results[["BleichrodtCRDI.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "RodriguezLogue") {
    returnValue <- stats::integrate(integrandRodriguezLogueLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["RodriguezLogue.lnk"]],
                                    beta = results[["RodriguezLogue.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Noise") {
    returnValue <- results[["Noise.mean"]]
  }

  returnValue
}
