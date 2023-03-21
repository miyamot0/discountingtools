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
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export plot.discountingtools
#' @export
plot.discountingtools <- function(fittingObject, which = "ind", position0 = "bottomleft", ylab0 = "Subjective Value", xlab0 = "Delay", logAxis = "x", yMin = 0.01, id = NULL, plotit = TRUE) {

  if (plotit) {
    if (which == "ind" & is.null(id))        plot_individual_rainbow(fittingObject,     position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "ind" & !is.null(id))       plot_individual_detailed(fittingObject,    position0, ylab0, xlab0, logAxis, yMin, id, plotit)
    if (which == "group")                    plot_group_rainbow(fittingObject,          position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "model")                    plot_model_characterization(fittingObject, position0, ylab0, xlab0, plotit)

    if (which == "ED50")                     plot_cross_rainbow(fittingObject, metric = "ProbableModel.LnED50",     plotit)
    if (which == "MBAUC")                    plot_cross_rainbow(fittingObject, metric = "ProbableModel.MBAUC",      plotit)
    if (which == "Log10MBAUC")               plot_cross_rainbow(fittingObject, metric = "ProbableModel.Log10MBAUC", plotit)
  } else {
    if (which == "ind" & is.null(id))        out = plot_individual_rainbow(fittingObject,     position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "ind" & !is.null(id))       out = plot_individual_detailed(fittingObject,    position0, ylab0, xlab0, logAxis, yMin, id, plotit)
    if (which == "group")                    out = plot_group_rainbow(fittingObject,          position0, ylab0, xlab0, logAxis, yMin, plotit)
    if (which == "model")                    out = plot_model_characterization(fittingObject, position0, ylab0, xlab0, plotit)

    if (which == "ED50")                     out = plot_cross_rainbow(fittingObject, metric = "LnED50",     plotit)
    if (which == "MBAUC")                    out = plot_cross_rainbow(fittingObject, metric = "MBAUC",      plotit)
    if (which == "Log10MBAUC")               out = plot_cross_rainbow(fittingObject, metric = "Log10MBAUC", plotit)

    return(out)
  }
}
