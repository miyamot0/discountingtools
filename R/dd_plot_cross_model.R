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
