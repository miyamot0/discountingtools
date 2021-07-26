#' minpack.lm logLik hack
#'
#' This function constructs a class, derived from an nls.lm object, similar to that of the logLik function in nls. This allows for native calls of the AIC and BIC functions from stats, using nls.lm fit objects.
#'
#' @param fit nls.lm fitted model
#' @param REML determine whether or not to use ML (FALSE by default)
#' @param ... inherit other args as necessary
#' @author Katharine Mullen <kate@@few.vu.nl>
#' @return provide a logLik class for AIC/BIC
logLik.nls.lm <- function(fit, REML = FALSE, ...)
{
  logLikelihood <- -length(fit$fvec) * (log(2 * pi) + 1 - log(length(fit$fvec)) + log(sum(fit$fvec^2)))/2

  attr(logLikelihood, "df") <- 1L + length(stats::coef(fit))
  attr(logLikelihood, "nobs") <- attr(logLikelihood, "nall") <- length(fit$fvec)

  class(logLikelihood) <- "logLik"

  logLikelihood
}

#' summary.discountingtools
#'
#' Override summary output
#'
#' @param fittingObject
#'
#' @return
#' @export summary.discountingtools
#' @export
summary.discountingtools <- function(fittingObject) {

  localCopy <- fittingObject$results

  buildColNames = c("ID")

  if (!is.null(fittingObject$settings[["Group"]])) {
    buildColNames = c(buildColNames, "Group")
  }

  for (m in fittingObject$models) {
    if (m == "noise") {
      buildColNames = c(buildColNames,
                        "Noise.Intercept",
                        "Noise.RMSE",
                        "Noise.BIC",
                        "Noise.AIC")
    } else if (m == "mazur") {
      buildColNames = c(buildColNames,
                        "Mazur.Lnk",
                        "Mazur.RMSE",
                        "Mazur.BIC",
                        "Mazur.AIC",
                        "Mazur.Status")
    } else if (m == "exponential") {
      buildColNames = c(buildColNames,
                        "Exponential.Lnk",
                        "Exponential.RMSE",
                        "Exponential.BIC",
                        "Exponential.AIC",
                        "Exponential.Status")
    } else if (m == "laibson") {
      buildColNames = c(buildColNames,
                        "Laibson.Beta",
                        "Laibson.Delta",
                        "Laibson.RMSE",
                        "Laibson.BIC",
                        "Laibson.AIC",
                        "Laibson.Status")
    } else if (m == "greenmyerson") {
      buildColNames = c(buildColNames,
                        "GreenMyerson.Lnk",
                        "GreenMyerson.S",
                        "GreenMyerson.RMSE",
                        "GreenMyerson.BIC",
                        "GreenMyerson.AIC",
                        "GreenMyerson.Status")
    } else if (m == "rachlin") {
      buildColNames = c(buildColNames,
                        "Rachlin.Lnk",
                        "Rachlin.S",
                        "Rachlin.RMSE",
                        "Rachlin.BIC",
                        "Rachlin.AIC",
                        "Rachlin.Status")
    } else if (m == "ebertprelec") {
      buildColNames = c(buildColNames,
                        "EbertPrelec.Lnk",
                        "EbertPrelec.S",
                        "EbertPrelec.RMSE",
                        "EbertPrelec.BIC",
                        "EbertPrelec.AIC",
                        "EbertPrelec.Status")
    } else if (m == "bleichrodt") {
      buildColNames = c(buildColNames,
                        "Bleichrodt.Lnk",
                        "Bleichrodt.S",
                        "Bleichrodt.Beta",
                        "Bleichrodt.RMSE",
                        "Bleichrodt.BIC",
                        "Bleichrodt.AIC",
                        "Bleichrodt.Status")
    } else if (m == "rodriguezlogue") {
      buildColNames = c(buildColNames,
                        "RodriguezLogue.Lnk",
                        "RodriguezLogue.Beta",
                        "RodriguezLogue.RMSE",
                        "RodriguezLogue.BIC",
                        "RodriguezLogue.AIC",
                        "RodriguezLogue.Status")

    }
  }

  if (fittingObject$ModelSelection == TRUE) buildColNames = append(buildColNames, c("ProbableModel",
                                                                                    "ProbableModel.BF",
                                                                                    "ProbableModel.Prob"))

  for (metric in fittingObject[["metrics"]]) {
    if (metric == "lned50")   buildColNames = append(buildColNames, c("LnED50"))
    if (metric == "mbauc")    buildColNames = append(buildColNames, c("MBAUC"))
    if (metric == "logmbauc") buildColNames = append(buildColNames, c("Log10MBAUC"))
  }

  nRows    = length(names(localCopy))
  resFrame = data.frame(matrix(ncol = length(buildColNames),
                               nrow = nRows))

  colnames(resFrame) <- buildColNames

  resFrame$ID <- names(localCopy)

  for (name in names(localCopy)) {
    index = which(names(localCopy) == name)

    for (res in localCopy[[name]]) {

      if (res$Model == "noise") {
        resFrame[index, c("Noise.Intercept",
                          "Noise.RMSE",
                          "Noise.BIC",
                          "Noise.AIC")] = as.data.frame(res)[, c("Intercept",
                                                                 "RMSE",
                                                                 "BIC",
                                                                 "AIC")]

      } else if (res$Model == "mazur") {
        resFrame[index, c("Mazur.Lnk",
                          "Mazur.RMSE",
                          "Mazur.BIC",
                          "Mazur.AIC",
                          "Mazur.Status")] = as.data.frame(res)[, c("Lnk",
                                                                    "RMSE",
                                                                    "BIC",
                                                                    "AIC",
                                                                    "Status")]
      } else if (res$Model == "exponential") {
        resFrame[index, c("Exponential.Lnk",
                          "Exponential.RMSE",
                          "Exponential.BIC",
                          "Exponential.AIC",
                          "Exponential.Status")] = as.data.frame(res)[, c("Lnk",
                                                                          "RMSE",
                                                                          "BIC",
                                                                          "AIC",
                                                                          "Status")]
      } else if (res$Model == "laibson") {
        resFrame[index, c("Laibson.Beta",
                          "Laibson.Delta",
                          "Laibson.RMSE",
                          "Laibson.BIC",
                          "Laibson.AIC",
                          "Laibson.Status")] = as.data.frame(res)[, c("Beta",
                                                                      "Delta",
                                                                      "RMSE",
                                                                      "BIC",
                                                                      "AIC",
                                                                      "Status")]
      } else if (res$Model == "greenmyerson") {
        resFrame[index, c("GreenMyerson.Lnk",
                          "GreenMyerson.S",
                          "GreenMyerson.RMSE",
                          "GreenMyerson.BIC",
                          "GreenMyerson.AIC",
                          "GreenMyerson.Status")] = as.data.frame(res)[, c("Lnk",
                                                                           "S",
                                                                           "RMSE",
                                                                           "BIC",
                                                                           "AIC",
                                                                           "Status")]
      } else if (res$Model == "rachlin") {
        resFrame[index, c("Rachlin.Lnk",
                          "Rachlin.S",
                          "Rachlin.RMSE",
                          "Rachlin.BIC",
                          "Rachlin.AIC",
                          "Rachlin.Status")] = as.data.frame(res)[, c("Lnk",
                                                                      "S",
                                                                      "RMSE",
                                                                      "BIC",
                                                                      "AIC",
                                                                      "Status")]
      } else if (res$Model == "ebertprelec") {
        resFrame[index, c("EbertPrelec.Lnk",
                          "EbertPrelec.S",
                          "EbertPrelec.RMSE",
                          "EbertPrelec.BIC",
                          "EbertPrelec.AIC",
                          "EbertPrelec.Status")] = as.data.frame(res)[, c("Lnk",
                                                                          "S",
                                                                          "RMSE",
                                                                          "BIC",
                                                                          "AIC",
                                                                          "Status")]
      } else if (res$Model == "bleichrodt") {
        resFrame[index, c("Bleichrodt.Lnk",
                          "Bleichrodt.S",
                          "Bleichrodt.Beta",
                          "Bleichrodt.RMSE",
                          "Bleichrodt.BIC",
                          "Bleichrodt.AIC",
                          "Bleichrodt.Status")] = as.data.frame(res)[, c("Lnk",
                                                                         "S",
                                                                         "Beta",
                                                                         "RMSE",
                                                                         "BIC",
                                                                         "AIC",
                                                                         "Status")]
      } else if (res$Model == "rodriguezlogue") {
        resFrame[index, c("RodriguezLogue.Lnk",
                          "RodriguezLogue.Beta",
                          "RodriguezLogue.RMSE",
                          "RodriguezLogue.BIC",
                          "RodriguezLogue.AIC",
                          "RodriguezLogue.Status")] = as.data.frame(res)[, c("Lnk",
                                                                             "Beta",
                                                                             "RMSE",
                                                                             "BIC",
                                                                             "AIC",
                                                                             "Status")]
      }
    }

    if (fittingObject$ModelSelection == TRUE) {
      resFrame[index, "ProbableModel"]      = fittingObject$rotation[[name]]$ProbableModel
      resFrame[index, "ProbableModel.BF"]   = fittingObject$rotation[[name]]$ProbableModel.BF
      resFrame[index, "ProbableModel.Prob"] = fittingObject$rotation[[name]]$ProbableModel.Prob
    }

    for (metric in fittingObject[["metrics"]]) {
      if (metric == "lned50")   resFrame[index, "LnED50"]     = fittingObject$ed50[[name]]
      if (metric == "mbauc")    resFrame[index, "MBAUC"]      = fittingObject$mbauc[[name]]
      if (metric == "logmbauc") resFrame[index, "Log10MBAUC"] = fittingObject$mbauclog10[[name]]
    }

    if (!is.null(fittingObject$settings[["Group"]])) {
        resFrame[index, "Group"] = unique(results$data[
          which(results$data[,as.character(results$settings['Individual'])] == name),
          as.character(results$settings['Group'])])
    }
  }

  resFrame
}

#' plot.discountingtools
#'
#' Override plot output
#'
#' @param fittingObject core frame
#' @param which plottype
#' @param position0 position legend
#' @param ylab0 y axis label
#' @param xlab0 x axis label
#' @param logAxis axis designation
#' @param yMin y axis lower limit
#'
#' @return
#' @export plot.discountingtools
#' @export
plot.discountingtools <- function(fittingObject, which = "ind", position0 = "bottomleft", ylab0 = "Subjective Value", xlab0 = "Delay", logAxis = "x", yMin = 0.01) {

  if (which == "ind")        plotIndividualRainbow(fittingObject, position0, ylab0, xlab0, logAxis, yMin)
  if (which == "group")      plotGroupRainbow(fittingObject,      position0, ylab0, xlab0, logAxis, yMin)
  if (which == "ED50")       plotRainbowCross(fittingObject, metric = "LnED50")
  if (which == "MBAUC")      plotRainbowCross(fittingObject, metric = "MBAUC")
  if (which == "Log10MBAUC") plotRainbowCross(fittingObject, metric = "Log10MBAUC")
}

#' plotIndividualRainbow
#'
#' @param fittingObject core frame
#' @param position0 position legend
#' @param ylab0 y axis label
#' @param xlab0 x axis label
#' @param logAxis axis designation
#' @param yMin y axis lower limit
#'
#' @return
plotIndividualRainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin) {

  preDraw = TRUE
  yLimits = c(0, fittingObject$maxValue)

  vecModels = fittingObject$models
  vecColors = rainbow(length(vecModels), alpha = 1)

  legendBuildModel = NA
  legendBuildColor = NA

  for (id in names(fittingObject$results)) {

    ogData = subset(fittingObject$data, ids == id)

    # Hack: Check if even multiple models

    if (is.null(fittingObject$rotation)) {
      model = names(results$results[[id]])
    } else {
      model  = fittingObject$rotation[[id]]$ProbableModel
    }

    result = fittingObject$results[[id]][[model]]

    xs = seq(min(ogData[,as.character(fittingObject$settings['Delays'])]),
             max(ogData[,as.character(fittingObject$settings['Delays'])]), length.out = 2000)

    if (model == "noise")          yhat = rep(result$Intercept, length(xs))

    if (model == "bleichrodt")     yhat = BleichrodtCRDIDiscountFunc(xs,     result$Lnk,  result$S, result$Beta)
    if (model == "ebertprelec")    yhat = ebertPrelecDiscountFunc(xs,        result$Lnk,  result$S)
    if (model == "exponential")    yhat = exponentialDiscountFunc(xs,        result$Lnk)
    if (model == "greenmyerson")   yhat = myersonHyperboloidDiscountFunc(xs, result$Lnk,  result$S)
    if (model == "laibson")        yhat = betaDeltaDiscountFunc(xs,          result$Beta, result$Delta)
    if (model == "mazur")          yhat = hyperbolicDiscountFunc(xs,         result$Lnk)
    if (model == "rachlin")        yhat = rachlinHyperboloidDiscountFunc(xs, result$Lnk,  result$S)
    if (model == "rodriguezlogue") yhat = RodriguezLogueDiscountFunc(xs,     result$Lnk,  result$Beta)

    if (length(vecColors) == 1) {
      col = vecColors
    } else {
      col = vecColors[match(model, vecModels)]
    }

    if (!(model %in% legendBuildModel)) {
      if (is.na(legendBuildModel)) {
        legendBuildModel = c(model)
        legendBuildColor = c(col)
      } else {
        legendBuildModel = c(legendBuildModel, model)
        legendBuildColor = c(legendBuildColor, col)
      }
    }

    if (grepl("y", logAxis) == TRUE) {
      yhat    = yhat[yhat >= 0]
      yLimits = c(yMin, fittingObject$maxValue)
    }

    if (preDraw) {
      plot(xs, yhat * fittingObject$maxValue,
           type = "l",
           ylim = yLimits,
           log  = logAxis,
           main = "Summary Fits",
           col  = col,
           ylab = ylab0,
           xlab = xlab0)

      preDraw = FALSE
    } else {
      lines(xs, yhat * fittingObject$maxValue,
            col  = col)
    }
  }

  legend(position0,
         legend = legendBuildModel,
         col    = legendBuildColor,
         lty    = 1)

}

#' plotGroupRainbow
#'
#' @param fittingObject core frame
#' @param position0 position legend
#' @param ylab0 y axis label
#' @param xlab0 x axis label
#' @param logAxis axis designation
#' @param yMin y axis lower limit
#'
#' @return
plotGroupRainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin) {

  if (is.null(fittingObject$settings[["Group"]])) stop('No Group aesthetic specified')

  preDraw = TRUE
  yLimits = c(0, fittingObject$maxValue)

  vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])
  vecColors = rainbow(length(vecGroups), alpha = 1)

  #unique(results$data[,as.character(results$settings['Group'])])

  for (id in names(fittingObject$results)) {

    ogData = subset(fittingObject$data, ids == id)

    model  = fittingObject$rotation[[id]]$ProbableModel
    result = fittingObject$results[[id]][[model]]

    xs = seq(min(ogData[,as.character(fittingObject$settings['Delays'])]),
             max(ogData[,as.character(fittingObject$settings['Delays'])]), length.out = 2000)

    if (model == "noise")          yhat = rep(result$Intercept, length(xs))

    if (model == "bleichrodt")     yhat = BleichrodtCRDIDiscountFunc(xs,     result$Lnk,  result$S, result$Beta)
    if (model == "ebertprelec")    yhat = ebertPrelecDiscountFunc(xs,        result$Lnk,  result$S)
    if (model == "exponential")    yhat = exponentialDiscountFunc(xs,        result$Lnk)
    if (model == "greenmyerson")   yhat = myersonHyperboloidDiscountFunc(xs, result$Lnk,  result$S)
    if (model == "laibson")        yhat = betaDeltaDiscountFunc(xs,          result$Beta, result$Delta)
    if (model == "mazur")          yhat = hyperbolicDiscountFunc(xs,         result$Lnk)
    if (model == "rachlin")        yhat = rachlinHyperboloidDiscountFunc(xs, result$Lnk,  result$S)
    if (model == "rodriguezlogue") yhat = RodriguezLogueDiscountFunc(xs,     result$Lnk,  result$Beta)

    col = vecColors[match(ogData[1, as.character(fittingObject$settings['Group'])], vecGroups)]

    if (grepl("y", logAxis) == TRUE) {
      yhat    = yhat[yhat >= 0]
      yLimits = c(yMin, fittingObject$maxValue)
    }

    if (preDraw) {
      plot(xs, yhat * fittingObject$maxValue,
           type = "l",
           ylim = yLimits,
           log  = logAxis,
           main = "Summary Fits",
           col  = col,
           ylab = ylab0,
           xlab = xlab0)

      preDraw = FALSE
    } else {
      lines(xs, yhat * fittingObject$maxValue,
            col  = col)
    }
  }

  legend(position0,
         legend = vecGroups,
         col    = vecColors,
         lty    = 1)
}

#' plotRainbowCross
#'
#' @param fittingObject core fitting object
#' @param metric ...
#'
#' @return
plotRainbowCross <- function(fittingObject, metric) {

  vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])
  vecColors = rainbow(length(vecGroups), alpha = 1)

  resultFrame = summary(fittingObject)

  print(histogram(as.formula(paste("~", metric)),
            data   = resultFrame,
            type   = "p",
            groups = Group,
            panel  = function(...)
              panel.superpose(...,
                              panel.groups = panel.histogram,
                              col          = vecColors,
                              alpha        = 0.5),
                              auto.key     = list(columns    = length(vecColors),
                                                  rectangles = FALSE,
                                                  col        = vecColors))
  )
}

#' messageDebug
#'
#' @param fittingObject core fitting object
#' @param msg message
#'
#' @return
messageDebug <- function(fittingObject, msg) {
  if (fittingObject[[ "verbose"  ]] == TRUE) message(msg)
}
