#' plot_individual_detailed
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom grDevices rainbow
#' @importFrom graphics lines legend
plot_individual_detailed <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, id, plotit) {
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

      if (model == "bleichrodt")     yhat = dd_discount_func_bleichrodt_crdi(xs, result$Lnk,  result$S, result$Beta)
      if (model == "ebertprelec")    yhat = dd_discount_func_ebertprelec(xs,     result$Lnk,  result$S)
      if (model == "exponential")    yhat = dd_discount_func_exponential(xs,     result$Lnk)
      if (model == "greenmyerson")   yhat = dd_discount_func_greenmyerson(xs,    result$Lnk,  result$S)
      if (model == "laibson")        yhat = dd_discount_func_laibson(xs,         result$Beta, result$Delta)
      if (model == "mazur")          yhat = dd_discount_func_mazur(xs,           result$Lnk)
      if (model == "rachlin")        yhat = dd_discount_func_rachlin(xs,         result$Lnk,  result$S)
      if (model == "rodriguezlogue") yhat = dd_discount_func_rodriguezlogue(xs,  result$Lnk,  result$Beta)

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
