#' plot_group_rainbow
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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom grDevices rainbow
#' @importFrom graphics lines legend
plot_group_rainbow <- function(fittingObject, position0, ylab0, xlab0, logAxis, yMin, plotit) {

  if (is.null(fittingObject$settings[["Group"]])) stop('No Group aesthetic specified')

  results = summary(fittingObject)

  if (!("ProbableModel" %in% names(results)) & length(fittingObject$models) > 1) {
    stop('Cannot plot individual fits without selecting a single model or using model selection')
  }

  if (plotit) {
    preDraw = TRUE
    yLimits = c(0, fittingObject$maxValue)

    vecGroups = unique(fittingObject$data[,as.character(fittingObject$settings['Group'])])
    vecColors = rainbow(length(vecGroups), alpha = 1)

    for (id in names(fittingObject$results)) {
      ogData = fittingObject$data[
        which(fittingObject$data[,as.character(fittingObject$settings['Individual'])] == id),]

      model = fittingObject$models[1]

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

      ogData = fittingObject$data[
        which(fittingObject$data[,as.character(fittingObject$settings['Individual'])] == id),]

      model = fittingObject$models[1]

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
