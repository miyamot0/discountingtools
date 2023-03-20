#' plotModelCharacterization
#'
#' @param fittingObject core fitting object
#' @param position0 (char) position of legend
#' @param ylab0 (char) y axis label
#' @param xlab0 (char) x axis label
#' @param plotit (logical) bool of whether or not to print visual or output plotting frame
#'
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
                     main = "Model Characterization",
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
                     stack = TRUE,
                     scales = list(x = list(rot = 45))))
    }
  }

  if (!plotit) prePlotDfFinal
}
