#' minpack.lm logLik hack
#'
#' This function constructs a class, derived from an nls.lm object, similar to that of the logLik function in nls. This allows for native calls of the AIC and BIC functions from stats, using nls.lm fit objects.
#'
#' @param fit nls.lm fitted model
#' @param REML determine whether or not to use ML (FALSE by default)
#' @param ... inherit other args as necessary
#'
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

#' Perform Johnson & Bickel Screen
#'
#' This function applies the Johnson & Bickel screening criteria to included data series. The result of this procedure is a TRUE/FALSE response to one of two screening criteria.
#'
#' @param fittingObject core fitting object
#'
#' @return A data frame of model screenings
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return
#' @export
johnsonBickelScreen <- function(fittingObject) {

  listOfIds = unique(fittingObject$data[[as.character(fittingObject$settings['Individual'])]])

  for (id in listOfIds) {
    messageDebug(fittingObject, paste("JB Screen: ", id))

    currentData = fittingObject$data[
      which(fittingObject$data[,
           as.character(fittingObject$settings['Individual'])] == id),]

    fittingObject$data[
      which(fittingObject$data[,
           as.character(fittingObject$settings['Individual'])] == id), "JB1"] = TRUE

    fittingObject$data[
      which(fittingObject$data[,
           as.character(fittingObject$settings['Individual'])] == id), "JB2"] = TRUE

    currentData$ddX = currentData[,as.character(fittingObject$settings['Delays'])]
    currentData$ddY = currentData[,as.character(fittingObject$settings['Values'])]
    currentData$ddY = currentData$ddY / as.numeric(fittingObject[[ "maxValue" ]])

    currentData = currentData[order(currentData$ddX), ]

    for (index in 2:length(currentData$ddX)) {
      prev = currentData[index - 1, "ddY"]
      curr = currentData[index,     "ddY"]

      if ((curr - prev) > as.numeric(fittingObject[[ "JB1Flag" ]])) {
        messageDebug(fittingObject, paste("JB Screen: ", id, "[Fail JB1]"))

        fittingObject$data[
          which(fittingObject$data[,
            as.character(fittingObject$settings['Individual'])] == id), "JB1"] = FALSE
      }
    }

    prev <- currentData[1,                       "ddY"]
    curr <- currentData[length(currentData$ddX), "ddY"]

    if ((prev - curr) < as.numeric(fittingObject[[ "JB2Flag" ]])) {
      messageDebug(fittingObject, paste("JB Screen: ", id, "[Fail JB2]"))

      fittingObject$data[
        which(fittingObject$data[,
          as.character(fittingObject$settings['Individual'])] == id), "JB2"] = FALSE
    }
  }

  fittingObject
}

#' summary.discountingtools
#'
#' Override summary output. Rather than display the core fitting object, a data frame block of results is provided to the user for easy interpretation and further analysis
#'
#' @param fittingObject core fitting object
#' @param detailed enable additional model metrics (default TRUE)
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export summary.discountingtools
#' @export
summary.discountingtools <- function(fittingObject, detailed = TRUE) {

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

  if (fittingObject$ModelSelection == TRUE)
    buildColNames = append(buildColNames, c("ProbableModel",
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

  if (detailed == FALSE) {
    resFrame = resFrame[,!grepl(".RMSE",   colnames(resFrame))]
    resFrame = resFrame[,!grepl(".AIC",    colnames(resFrame))]
    resFrame = resFrame[,!grepl(".Status", colnames(resFrame))]
    resFrame = resFrame[,!grepl(".BF",     colnames(resFrame))]
    resFrame = resFrame[,!grepl(".Prob",   colnames(resFrame))]
  }

  resFrame
}

#' plot.discountingtools
#'
#' This method overrides the base plot function to provide various plots relevant to the user.
#'
#' @param fittingObject core fitting object
#' @param which (char) type of plot to show, based on fits
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param logAxis (char) axis designation
#' @param yMin (num) y axis lower limit
#' @param id (num) participant number to focus
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export plot.discountingtools
#' @export
plot.discountingtools <- function(fittingObject, which = "ind", position0 = "bottomleft", ylab0 = "Subjective Value", xlab0 = "Delay", logAxis = "x", yMin = 0.01, id = NULL, plotit = TRUE) {

  if (plotit) {
    if (which == "ind" & is.null(id))        plotIndividualRainbow(fittingObject,     position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "ind" & !is.null(id))       plotIndividualDetailed(fittingObject,    position0, ylab0, xlab0, logAxis, yMin, id, plotit)
    if (which == "group")                    plotGroupRainbow(fittingObject,          position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "model")                    plotModelCharacterization(fittingObject, position0, ylab0, xlab0, plotit)

    if (which == "ED50")                     plotRainbowCross(fittingObject, metric = "LnED50",     plotit)
    if (which == "MBAUC")                    plotRainbowCross(fittingObject, metric = "MBAUC",      plotit)
    if (which == "Log10MBAUC")               plotRainbowCross(fittingObject, metric = "Log10MBAUC", plotit)
  } else {
    if (which == "ind" & is.null(id))        out = plotIndividualRainbow(fittingObject,     position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "ind" & !is.null(id))       out = plotIndividualDetailed(fittingObject,    position0, ylab0, xlab0, logAxis, yMin, id, plotit)
    if (which == "group")                    out = plotGroupRainbow(fittingObject,          position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "model")                    out = plotModelCharacterization(fittingObject, position0, ylab0, xlab0, plotit)

    if (which == "ED50")                     out = plotRainbowCross(fittingObject, metric = "LnED50",     plotit)
    if (which == "MBAUC")                    out = plotRainbowCross(fittingObject, metric = "MBAUC",      plotit)
    if (which == "Log10MBAUC")               out = plotRainbowCross(fittingObject, metric = "Log10MBAUC", plotit)

    return(out)
  }
}

#' plotIndividualRainbow
#'
#' This specific implementation shows cross-model fits, with series characterized by different models illustrated with different colors. A legend is also provided for convenience of interpretation.
#'
#' @param fittingObject core fitting object
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param logAxis (char) axis designation
#' @param yMin (num) y axis lower limit
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plotIndividualRainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, plotit) {

  if (plotit) {
    preDraw = TRUE
    yLimits = c(0, fittingObject$maxValue)

    vecModels = fittingObject$models
    vecColors = rainbow(length(vecModels), alpha = 1)

    preBuiltLegend   = FALSE
    legendBuildModel = NA
    legendBuildColor = NA

    for (id in names(fittingObject$results)) {

      ogData = subset(fittingObject$data, ids == id)

      # Hack: Check if even multiple models

      if (is.null(fittingObject$rotation)) {
        model = names(fittingObject$results[[id]])
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

      modelP = gsub("ebertprelec",    "ebert prelec",    model)
      modelP = gsub("greenmyerson",   "green myerson",   modelP)
      modelP = gsub("rodriguezlogue", "rodriguez logue", modelP)

      modelC = tools::toTitleCase(modelP)

      if (!(modelC %in% legendBuildModel)) {
        if (!preBuiltLegend) {
          legendBuildModel = c(modelC)
          legendBuildColor = c(col)

          preBuiltLegend   = TRUE
        } else {
          legendBuildModel = c(legendBuildModel, modelC)
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
           lty    = 1,
           bty    = "n")

  } else {
    outputframe = NULL

    for (id in names(fittingObject$results)) {

      ogData = subset(fittingObject$data, ids == id)

      # Hack: Check if even multiple models

      if (is.null(fittingObject$rotation)) {
        model = names(fittingObject$results[[id]])
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

      modelP = gsub("ebertprelec",    "ebert prelec",    model)
      modelP = gsub("greenmyerson",   "green myerson",   modelP)
      modelP = gsub("rodriguezlogue", "rodriguez logue", modelP)

      modelC = tools::toTitleCase(modelP)

      if (!(modelC %in% legendBuildModel)) {
        if (!preBuiltLegend) {
          legendBuildModel = c(modelC)
          legendBuildColor = c(col)

          preBuiltLegend   = TRUE
        } else {
          legendBuildModel = c(legendBuildModel, modelC)
          legendBuildColor = c(legendBuildColor, col)
        }
      }

      if (grepl("y", logAxis) == TRUE) {
        yhat    = yhat[yhat >= 0]
        yLimits = c(yMin, fittingObject$maxValue)
      }

      tempFrame = data.frame(
        ID    = rep(id, length(xs)),
        X     = xs,
        Y     = yhat * fittingObject$maxValue,
        Model = rep(modelC, length(xs))
      )

      if (is.null(outputframe)) {
        outputframe = tempFrame
      } else {
        outputframe = rbind(outputframe,
                            tempFrame)
      }
    }
  }

  if (!plotit) outputframe
}

#' plotIndividualDetailed
#'
#' This implementation of plot singles out a particular responder, providing the fits to the observed data as well as the probability that the "probable" model characterizes the data
#'
#' @param fittingObject core fitting object
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param logAxis (char) axis designation
#' @param yMin (num) y axis lower limit
#' @param id (num) participant id
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plotIndividualDetailed <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, id, plotit) {
  if (!(id %in% names(fittingObject$results))) stop('id not found in results')

  if (plotit) {
    if (grepl("y", logAxis) == TRUE) {
      yLimits = c(yMin, fittingObject$maxValue)
    } else {
      yLimits = c(0,    fittingObject$maxValue)
    }

    vecModels = fittingObject$models
    vecColors = rainbow(length(vecModels), alpha = 1)

    preBuiltLegend   = FALSE
    legendBuildModel = NA
    legendBuildColor = NA

    ogData = subset(fittingObject$data, ids == id)

    plot(ogData[,as.character(fittingObject$settings['Delays'])],
         ogData[,as.character(fittingObject$settings['Values'])],
         type = "p",
         ylim = yLimits,
         log  = logAxis,
         main = "Summary Fits",
         col  = "black",
         pch  = 19,
         ylab = ylab0,
         xlab = xlab0)

    for (model in names(fittingObject$results[[id]])) {
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
        probString = ""

        modelP = gsub("ebertprelec",    "ebert prelec",    model)
        modelP = gsub("greenmyerson",   "green myerson",   modelP)
        modelP = gsub("rodriguezlogue", "rodriguez logue", modelP)

        modelC = tools::toTitleCase(modelP)
        modelC = gsub(" ", "~", modelC)

        if (grepl(model, fittingObject$rotation[[id]]$ProbableModel) == TRUE) {
          probString = paste0("~(",
                              as.character(round(fittingObject$rotation[[id]]$ProbableModel.Prob, 3)),
                              ")")

          modelP = parse(text = paste0("bold(",modelC, probString,")"))
        } else {
          modelP = parse(text = paste0(modelC, probString))
        }

        if (!preBuiltLegend) {
          legendBuildModel = c(modelP)
          legendBuildColor = c(col)

          preBuiltLegend   = TRUE
        } else {
          legendBuildModel = c(legendBuildModel, modelP)
          legendBuildColor = c(legendBuildColor, col)
        }
      }

      if (grepl("y", logAxis) == TRUE) {
        yhat    = yhat[yhat >= 0]
      }

      lines(xs, yhat * fittingObject$maxValue,
            col  = col)
    }

    legend(position0,
           legend = legendBuildModel,
           col    = legendBuildColor,
           lty    = 1,
           bty    = "n")
  } else {
    outputframe = NULL
    legendBuildModel = NA

    ogData = subset(fittingObject$data, ids == id)

    for (model in names(fittingObject$results[[id]])) {
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

      if (!(model %in% legendBuildModel)) {
        probString = ""

        modelP = gsub("ebertprelec",    "ebert prelec",    model)
        modelP = gsub("greenmyerson",   "green myerson",   modelP)
        modelP = gsub("rodriguezlogue", "rodriguez logue", modelP)

        modelC = tools::toTitleCase(modelP)
        modelC = gsub(" ", "~", modelC)

        if (grepl(model, fittingObject$rotation[[id]]$ProbableModel) == TRUE) {
          probString = paste0("~(",
                              as.character(round(fittingObject$rotation[[id]]$ProbableModel.Prob, 3)),
                              ")")

          modelP = parse(text = paste0("bold(",modelC, probString,")"))
        } else {
          modelP = parse(text = paste0(modelC, probString))
        }
      }

      tempFrame = data.frame(
        ID    = rep(id, length(xs)),
        X     = xs,
        Y     = yhat * fittingObject$maxValue,
        Model = rep(modelC, length(xs))
      )

      if (is.null(outputframe)) {
        outputframe = tempFrame
      } else {
        outputframe = rbind(outputframe,
                            tempFrame)
      }
    }
  }

  if (!plotit) outputframe
}

#' plotGroupRainbow
#'
#' Convenience method for illustrating individual fits when characterized by some a priori grouping.
#'
#' @param fittingObject core fitting object
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param logAxis (char) axis designation
#' @param yMin (num) y axis lower limit
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plotGroupRainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, plotit) {

  if (is.null(fittingObject$settings[["Group"]])) stop('No Group aesthetic specified')

  if (plotit) {
    preDraw = TRUE
    yLimits = c(0, fittingObject$maxValue)

    vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])
    vecColors = rainbow(length(vecGroups), alpha = 1)

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
  } else {
    outputframe = NULL

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

      tempFrame = data.frame(
        ID    = rep(id, length(xs)),
        Group = rep(ogData[1, as.character(fittingObject$settings['Group'])], length(xs)),
        X     = xs,
        Y     = yhat * fittingObject$maxValue,
        Model = rep(model, length(xs))
      )

      if (is.null(outputframe)) {
        outputframe = tempFrame
      } else {
        outputframe = rbind(outputframe,
                            tempFrame)
      }
    }
  }

  if (!plotit) outputframe
}

#' plotRainbowCross
#'
#' @param fittingObject core fitting object
#' @param metric (char) the cross model metric to be displayed
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plotRainbowCross <- function(fittingObject, metric, plotit) {

  if (!("Group" %in% names(fittingObject$settings))) {
    vecGroups = "sample"

    vecColors = rainbow(length(vecGroups), alpha = 1)

    resultFrame = summary(fittingObject)

    if (plotit) {
      print(histogram(as.formula(paste("~", metric)),
                      data   = resultFrame,
                      type   = "p"))
    }
  } else {
    vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])

    vecColors = rainbow(length(vecGroups), alpha = 1)

    resultFrame = summary(fittingObject)

    if (plotit) {
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
                                          col        = vecColors)))
    }
  }

  if (!plotit) resultFrame
}

#' plotModelCharacterization
#'
#' @param fittingObject core fitting object
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plotModelCharacterization <- function(fittingObject, position0, ylab0, xlab0, plotit) {

  if (!("Group" %in% names(fittingObject$settings))) {
    resultFrame = summary(fittingObject)

    prePlot = table(resultFrame$ProbableModel)
    prePlotDf = data.frame(
      Counts = as.numeric(prePlot),
      Model  = attr(prePlot, "dimnames")[[1]]
    )

    prePlotDfFinal = prePlotDf

    if (plotit) {
      print(barchart(Counts ~ Model,
               data = prePlotDfFinal,
               #groups = year,
               main = "Model Characterization",
               #xlab = "Yield Value",
               #stack = TRUE,
               #auto.key = list(space = "right"),
               scales = list(x = list(rot = 45))))
    }
  } else {
    resultFrame = summary(fittingObject)

    prePlotDfFinal = NULL

    for (grp in unique(resultFrame$Group)) {
      subsetFrame = subset(resultFrame, Group == grp)

      prePlot = table(subsetFrame$ProbableModel)

      prePlotDf = data.frame(
        Counts = as.numeric(prePlot),
        Model  = attr(prePlot, "dimnames")[[1]],
        Group  = rep(grp, length(as.numeric(prePlot)))
      )

      if (is.null(prePlotDfFinal)) {
        prePlotDfFinal = prePlotDf
      } else {
        prePlotDfFinal = rbind(prePlotDfFinal,
                               prePlotDf)
      }
    }

    if (plotit) {
      print(barchart(Counts ~ Model | Group,
               data = prePlotDfFinal,
               groups = Group,
               main = "Model Characterization",
               #xlab = "Yield Value",
               stack = TRUE,
               #auto.key = list(space = "right"),
               scales = list(x = list(rot = 45))))
    }
  }

  if (!plotit) prePlotDfFinal
}

#' messageDebug
#'
#' Extension of message method, instead yolked to a flag defining level of verbosity.
#'
#' @param fittingObject core fitting object
#' @param msg (char) message
#'
#' @return
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
messageDebug <- function(fittingObject, msg) {
  if (fittingObject[[ "verbose"  ]] == TRUE) message(msg)
}
