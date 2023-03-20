#' plot_individual_rainbow
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
plot_individual_rainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, plotit) {

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
      if (model == "bleichrodt")     yhat = dd_discount_func_bleichrodt_crdi(xs, result$Lnk,  result$S, result$Beta)
      if (model == "ebertprelec")    yhat = dd_discount_func_ebertprelec(xs,     result$Lnk,  result$S)
      if (model == "exponential")    yhat = dd_discount_func_exponential(xs,     result$Lnk)
      if (model == "greenmyerson")   yhat = dd_discount_func_greenmyerson(xs,    result$Lnk,  result$S)
      if (model == "laibson")        yhat = dd_discount_func_laibson(xs,         result$Beta, result$Delta)
      if (model == "mazur")          yhat = dd_discount_func_mazur(xs,           result$Lnk)
      if (model == "rachlin")        yhat = dd_discount_func_rachlin(xs,         result$Lnk,  result$S)
      if (model == "rodriguezlogue") yhat = dd_discount_func_rodriguezlogue(xs,  result$Lnk,  result$Beta)

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
      if (model == "bleichrodt")     yhat = dd_discount_func_bleichrodt_crdi(xs, result$Lnk,  result$S, result$Beta)
      if (model == "ebertprelec")    yhat = dd_discount_func_ebertprelec(xs,     result$Lnk,  result$S)
      if (model == "exponential")    yhat = dd_discount_func_exponential(xs,     result$Lnk)
      if (model == "greenmyerson")   yhat = dd_discount_func_greenmyerson(xs,    result$Lnk,  result$S)
      if (model == "laibson")        yhat = dd_discount_func_laibson(xs,         result$Beta, result$Delta)
      if (model == "mazur")          yhat = dd_discount_func_mazur(xs,           result$Lnk)
      if (model == "rachlin")        yhat = dd_discount_func_rachlin(xs,         result$Lnk,  result$S)
      if (model == "rodriguezlogue") yhat = dd_discount_func_rodriguezlogue(xs,  result$Lnk,  result$Beta)

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
