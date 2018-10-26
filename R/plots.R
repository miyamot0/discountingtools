#' Display of all fitted series with ED50 metric
#'
#' This method constructs a figure that displays all fitted models as well as the probability that they are the "true" model.  The ED50 metric is also provided for the most probable model.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @param lineWidth Line width
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom reshape melt
#' @importFrom ggplot2 geom_point geom_line aes annotation_logticks element_blank element_line element_rect element_text expand_limits guide_legend guides labs scale_x_continuous scale_colour_manual theme theme_bw xlab ylab
#' @importFrom scales trans_breaks trans_format math_format
#' @return display figure
displayED50Figure <- function(dat, results, lineWidth = 1) {
  samuelsonK <- NA
  ainslieK <- NA
  betaConstant <- NA
  deltaConstant <- NA
  myerK <- NA
  myerS <- NA
  rachK <- NA
  rachS <- NA
  epK <- NA
  epS <- NA
  crdiK <- NA
  crdiS <- NA
  crdiBeta <- NA
  ghK <- NA
  ghBeta <- NA

  digitPrecision <- 3

  endDelay <- max(dat$X)
  delaySeries = 1:(endDelay+1)
  expSeries  = rep(NA,endDelay+1)
  hypSeries  = rep(NA,endDelay+1)
  quaSeries  = rep(NA,endDelay+1)
  myerSeries = rep(NA,endDelay+1)
  rachSeries = rep(NA,endDelay+1)
  epSeries = rep(NA,endDelay+1)
  crdiSeries = rep(NA,endDelay+1)
  ghSeries = rep(NA,endDelay+1)

  legend = c(paste("Noise(",
                   round(results[["Noise.prob"]], digitPrecision),
                   ")",
                   sep = ""))
  colors = c("red")

  if ("Exponential.lnk" %in% names(results)) {
    samuelsonK <- results[["Exponential.lnk"]]
    legend = c(legend, paste("Exponential(",
                             round(results[["Exponential.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "cadetblue1")
  }

  if ("Hyperbolic.lnk" %in% names(results)) {
    ainslieK <- results[["Hyperbolic.lnk"]]
    legend = c(legend, paste("Hyperbolic(",
                             round(results[["Hyperbolic.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "chartreuse")
  }

  if ("Laibson.beta" %in% names(results)) {
    betaConstant <- results[["Laibson.beta"]]
    deltaConstant <- results[["Laibson.delta"]]
    legend = c(legend, paste("Laibson(",
                             round(results[["Laibson.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "darkslategray")
  }

  if ("GreenMyerson.lnk" %in% names(results)) {
    myerK <- results[["GreenMyerson.lnk"]]
    myerS <- results[["GreenMyerson.s"]]
    legend = c(legend, paste("GreenMyerson(",
                             round(results[["GreenMyerson.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "purple")
  }

  if ("Rachlin.lnk" %in% names(results)) {
    rachK <- results[["Rachlin.lnk"]]
    rachS <- results[["Rachlin.s"]]
    legend = c(legend, paste("Rachlin(",
                             round(results[["Rachlin.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "coral3")
  }

  if ("EbertPrelec.lnk" %in% names(results)) {
    epK <- results[["EbertPrelec.lnk"]]
    epS <- results[["EbertPrelec.s"]]
    legend = c(legend, paste("EbertPrelec(",
                             round(results[["EbertPrelec.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "dodgerblue")
  }

  if ("BleichrodtCRDI.lnk" %in% names(results)) {
    crdiK <- results[["BleichrodtCRDI.lnk"]]
    crdiS <- results[["BleichrodtCRDI.s"]]
    crdiBeta <- results[["BleichrodtCRDI.beta"]]
    legend = c(legend, paste("Bleichrodt CRDI(",
                             round(results[["BleichrodtCRDI.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "springgreen4")
  }

  if ("RodriguezLogue.lnk" %in% names(results)) {
    ghK <- results[["RodriguezLogue.lnk"]]
    ghBeta <- results[["RodriguezLogue.beta"]]
    legend = c(legend, paste("RodriguezLogue(",
                             round(results[["RodriguezLogue.prob"]], digitPrecision),
                             ")",
                             sep = ""))
    colors = c(colors, "burlywood")
  }

  for (delay in delaySeries)
  {
    delaySeries[delay] = delay-1

    if(!is.na(samuelsonK))
    {
      expSeries[delay] = (1.0 * exp(-(exp(samuelsonK))*delay))
    }

    if(!is.na(ainslieK))
    {
      hypSeries[delay] = 1.0 * (1+exp(ainslieK)*delay)^(-1)
    }

    if(!is.na(betaConstant))
    {
      quaSeries[delay] = 1.0 * ((betaConstant)*(deltaConstant)^delay)
    }

    if(!is.na(myerK))
    {
      myerSeries[delay] = 1.0 * (1+exp(myerK)*delay)^(-myerS)
    }

    if(!is.na(rachK))
    {
      rachSeries[delay] = 1.0 * (1 + exp(rachK)*(delay^rachS))^(-1)
    }

    if(!is.na(epK))
    {
      epSeries[delay] = 1.0 * exp(-(exp(epK)*delay)^epS)
    }

    if(!is.na(crdiK))
    {
      crdiSeries[delay] = 1.0 * crdiBeta * exp(-exp(crdiK)*delay^crdiS)
    }

    if(!is.na(ghK))
    {
      ghSeries[delay] = 1.0 * (1 + delay * exp(ghK))^(-ghBeta / exp(ghK))
    }
  }

  totalFrame = data.frame(Delays = delaySeries)

  mData <- data.frame(X = dat$X, Y = dat$Y)

  totalFrame = data.frame(Delays = delaySeries)

  if(!is.na(samuelsonK))
  {
    totalFrame["Exponential"] <- expSeries
  }

  if(!is.na(ainslieK))
  {
    totalFrame["Hyperbolic"] <- hypSeries
  }

  if(!is.na(betaConstant))
  {
    totalFrame["QuasiHyperbolic"] <- quaSeries
  }

  if(!is.na(myerK))
  {
    totalFrame["HyperboloidM"] <- myerSeries
  }

  if(!is.na(rachK))
  {
    totalFrame["HyperboloidR"] <- rachSeries
  }

  if(!is.na(epK))
  {
    totalFrame["EbertPrelec"] <- epSeries
  }

  if(!is.na(crdiK))
  {
    totalFrame["BleichrodtCRDI"] <- crdiSeries
  }

  if(!is.na(ghK))
  {
    totalFrame["RodriguezLogue"] <- ghSeries
  }

  totalFrame$Noise <- results[["Noise.mean"]]

  totalFrame.melt <- reshape::melt(totalFrame, id = c("Delays"))

  logChart <- ggplot2::ggplot() +
    geom_line(data = totalFrame.melt,
              size = 0.65,
              aes(x = totalFrame.melt$Delays,
                  y = totalFrame.melt$value,
                  colour = totalFrame.melt$variable)) +
    scale_colour_manual(name = "Model\n(Probability)", labels = legend, values = colors) +
    geom_point(data = mData,
               aes(x = mData$X,
                   y = mData$Y),
               size = 2,
               shape = 21,
               show.legend = F) +
    expand_limits(y = 0) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    ylab("Value") +
    scale_x_continuous(trans = 'log',
                       limits = c(1 , endDelay),
                       breaks = trans_breaks('log', function(x) exp(x)),
                       labels = trans_format('log', math_format(.x))) +
    annotation_logticks(sides = "b") +
    xlab("ln(Delay)") +
    labs(title = paste("Participant ", dat$id[1], sep = ""),
         subtitle = paste("ln(", results[["probable.model"]],
                          " ED50) = ",
                          round(results[["probable.ED50"]], 5),
                          sep = "")) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "white",
                                      fill = FALSE,
                                      size = 0),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          text = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.title.align = 0.5,
          legend.position = "bottom",
          legend.key = element_rect(fill = "transparent",
                                    colour = "transparent")) +
    guides(col = guide_legend(ncol = 3))

  print(logChart)
}

#' Display most probable fitted series with AUC metric
#'
#' This method constructs a figure that displays the most probable model as well as the probability that it is the "true" model.  The model-based AUC metric is also provided for the most probable model.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @param lineWidth Line width
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return display figure
displayAUCFigure <- function(dat, results, lineWidth = 1) {

  samuelsonK <- NA
  ainslieK <- NA
  betaConstant <- NA
  deltaConstant <- NA
  myerK <- NA
  myerS <- NA
  rachK <- NA
  rachS <- NA
  epK <- NA
  epS <- NA
  crdiK <- NA
  crdiS <- NA
  crdiBeta <- NA
  ghK <- NA
  ghBeta <- NA

  endDelay <- max(dat$X)
  delaySeries = 1:(endDelay+1)
  expSeries  = rep(NA,endDelay+1)
  hypSeries  = rep(NA,endDelay+1)
  quaSeries  = rep(NA,endDelay+1)
  myerSeries = rep(NA,endDelay+1)
  rachSeries = rep(NA,endDelay+1)
  epSeries = rep(NA,endDelay+1)
  crdiSeries = rep(NA,endDelay+1)
  ghSeries = rep(NA,endDelay+1)

  legend = c("Empirical: ")
  colors = c("black", "black")

  if ("Exponential.lnk" %in% names(results)) {
    samuelsonK <- results[["Exponential.lnk"]]
  }

  if ("Hyperbolic.lnk" %in% names(results)) {
    ainslieK <- results[["Hyperbolic.lnk"]]
  }

  if ("Laibson.beta" %in% names(results)) {
    betaConstant <- results[["Laibson.beta"]]
    deltaConstant <- results[["Laibson.delta"]]
  }

  if ("GreenMyerson.lnk" %in% names(results)) {
    myerK <- results[["GreenMyerson.lnk"]]
    myerS <- results[["GreenMyerson.s"]]
  }

  if ("Rachlin.lnk" %in% names(results)) {
    rachK <- results[["Rachlin.lnk"]]
    rachS <- results[["Rachlin.s"]]
  }

  if ("EbertPrelec.lnk" %in% names(results)) {
    epK <- results[["EbertPrelec.lnk"]]
    epS <- results[["EbertPrelec.s"]]
  }

  if ("BleichrodtCRDI.lnk" %in% names(results)) {
    crdiK <- results[["BleichrodtCRDI.lnk"]]
    crdiS <- results[["BleichrodtCRDI.s"]]
    crdiBeta <- results[["BleichrodtCRDI.beta"]]
  }

  if ("RodriguezLogue.lnk" %in% names(results)) {
    ghK <- results[["RodriguezLogue.lnk"]]
    ghBeta <- results[["RodriguezLogue.beta"]]
  }

  for (delay in delaySeries)
  {
    delaySeries[delay] = delay-1

    if(!is.na(samuelsonK))
    {
      expSeries[delay] = (1.0 * exp(-(exp(samuelsonK))*delay))
    }

    if(!is.na(ainslieK))
    {
      hypSeries[delay] = 1.0 * (1+exp(ainslieK)*delay)^(-1)
    }

    if(!is.na(betaConstant))
    {
      quaSeries[delay] = 1.0 * ((betaConstant)*(deltaConstant)^delay)
    }

    if(!is.na(myerK))
    {
      myerSeries[delay] = 1.0 * (1+exp(myerK)*delay)^(-myerS)
    }

    if(!is.na(rachK))
    {
      rachSeries[delay] = 1.0 * (1 + exp(rachK)*(delay^rachS))^(-1)
    }

    if(!is.na(epK))
    {
      epSeries[delay] = 1.0 * exp(-(exp(epK)*delay)^epS)
    }

    if(!is.na(crdiK))
    {
      crdiSeries[delay] = 1.0 * crdiBeta * exp(-exp(crdiK)*delay^crdiS)
    }

    if(!is.na(ghK))
    {
      ghSeries[delay] = 1.0 * (1 + delay * exp(ghK))^(-ghBeta / exp(ghK))
    }
  }

  mData <- data.frame(X = dat$X, Y = dat$Y)

  lineColor <- "black"

  if (results[["probable.model"]] == "Hyperbolic") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = hypSeries)
    legend = c(legend, paste("Hyperbolic(",
                             round(results[["Hyperbolic.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Exponential") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = expSeries)
    legend = c(legend, paste("Exponential(",
                             round(results[["Exponential.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Laibson") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = quaSeries)
    legend = c(legend, paste("Laibson(",
                             round(results[["Laibson.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "GreenMyerson") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = myerSeries)
    legend = c(legend, paste("GreenMyerson(",
                             round(results[["GreenMyerson.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Rachlin") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = rachSeries)
    legend = c(legend, paste("Rachlin(",
                             round(results[["Rachlin.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "EbertPrelec") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = epSeries)
    legend = c(legend, paste("EbertPrelec(",
                             round(results[["EbertPrelec.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = crdiSeries)
    legend = c(legend, paste("BleichrodtCRDI(",
                             round(results[["BleichrodtCRDI.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "RodriguezLogue") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = ghSeries)
    legend = c(legend, paste("RodriguezLogue(",
                             round(results[["RodriguezLogue.prob"]], 5),
                             ")",
                             sep = ""))

  } else {
    totalFrame = data.frame(Delays = delaySeries)
    totalFrame$ModelArea <- results[["Noise.mean"]]

    legend = c(legend, paste("Noise(",
                             round(results[["Noise.prob"]], 5),
                             ")",
                             sep = ""))

  }

  graphics::plot(totalFrame$Delays, totalFrame$ModelArea, type = "l", ylim = c(0,1),
                 main = paste("Participant ", dat$id[1],
                              "\n", results[["probable.model"]], " Model Area = ", round(results[["probable.AUC"]], 5),
                              sep = ""),
                 xlab = "Delays",
                 ylab = "Value",
                 col = lineColor,
                 lwd = lineWidth)

  graphics::points(dat$X,
                   dat$Y,
                   type = "p",
                   cex = 2,
                   pch = 18)

  graphics::lines(dat$X,
                  dat$Y,
                  col = "black",
                  lty = 2,
                  lwd = 0.5)

  legend("topright",
         legend = legend,
         col = colors,
         lty = c(2, 1),
         lwd = 3,
         title = "Model (Probability)")
}

#' Display most probable fitted series with AUC metric (log10)
#'
#' This method constructs a figure that displays the most probable model as well as the probability that it is the "true" model.  The model-based AUC metric in log base 10 scale is also provided for the most probable model.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @param lineWidth Line width
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return display figure
displayLogAUCFigure <- function(dat, results, lineWidth = 1) {

  samuelsonK <- NA
  ainslieK <- NA
  betaConstant <- NA
  deltaConstant <- NA
  myerK <- NA
  myerS <- NA
  rachK <- NA
  rachS <- NA
  epK <- NA
  epS <- NA
  crdiK <- NA
  crdiS <- NA
  crdiBeta <- NA
  ghK <- NA
  ghBeta <- NA

  endDelay <- max(dat$X)
  delaySeries = 1:(endDelay+1)
  expSeries  = rep(NA,endDelay+1)
  hypSeries  = rep(NA,endDelay+1)
  quaSeries  = rep(NA,endDelay+1)
  myerSeries = rep(NA,endDelay+1)
  rachSeries = rep(NA,endDelay+1)
  epSeries = rep(NA,endDelay+1)
  crdiSeries = rep(NA,endDelay+1)
  ghSeries = rep(NA,endDelay+1)

  legend = c(paste("Noise: ", round(results[["Noise.prob"]], 5), sep = ""))
  colors = c("red")

  legend = c("Empirical:")
  colors = c("black", "black")

  if ("Exponential.lnk" %in% names(results)) {
    samuelsonK <- results[["Exponential.lnk"]]
  }

  if ("Hyperbolic.lnk" %in% names(results)) {
    ainslieK <- results[["Hyperbolic.lnk"]]
  }

  if ("Laibson.beta" %in% names(results)) {
    betaConstant <- results[["Laibson.beta"]]
    deltaConstant <- results[["Laibson.delta"]]
  }

  if ("GreenMyerson.lnk" %in% names(results)) {
    myerK <- results[["GreenMyerson.lnk"]]
    myerS <- results[["GreenMyerson.s"]]
  }

  if ("Rachlin.lnk" %in% names(results)) {
    rachK <- results[["Rachlin.lnk"]]
    rachS <- results[["Rachlin.s"]]
  }

  if ("EbertPrelec.lnk" %in% names(results)) {
    epK <- results[["EbertPrelec.lnk"]]
    epS <- results[["EbertPrelec.s"]]
  }

  if ("BleichrodtCRDI.lnk" %in% names(results)) {
    crdiK <- results[["BleichrodtCRDI.lnk"]]
    crdiS <- results[["BleichrodtCRDI.s"]]
    crdiBeta <- results[["BleichrodtCRDI.beta"]]
  }

  if ("RodriguezLogue.lnk" %in% names(results)) {
    ghK <- results[["RodriguezLogue.lnk"]]
    ghBeta <- results[["RodriguezLogue.beta"]]
  }

  for (delay in delaySeries)
  {
    delaySeries[delay] = delay-1

    if(!is.na(samuelsonK))
    {
      expSeries[delay] = (1.0 * exp(-(exp(samuelsonK))*delay))
    }

    if(!is.na(ainslieK))
    {
      hypSeries[delay] = 1.0 * (1+exp(ainslieK)*delay)^(-1)
    }

    if(!is.na(betaConstant))
    {
      quaSeries[delay] = 1.0 * ((betaConstant)*(deltaConstant)^delay)
    }

    if(!is.na(myerK))
    {
      myerSeries[delay] = 1.0 * (1+exp(myerK)*delay)^(-myerS)
    }

    if(!is.na(rachK))
    {
      rachSeries[delay] = 1.0 * (1 + exp(rachK)*(delay^rachS))^(-1)
    }

    if(!is.na(epK))
    {
      epSeries[delay] = 1.0 * exp(-(exp(epK)*delay)^epS)
    }

    if(!is.na(crdiK))
    {
      crdiSeries[delay] = 1.0 * crdiBeta * exp(-exp(crdiK)*delay^crdiS)
    }

    if(!is.na(ghK))
    {
      ghSeries[delay] = 1.0 * (1 + delay * exp(ghK))^(-ghBeta / exp(ghK))
    }
  }

  mData <- data.frame(X = dat$X, Y = dat$Y)

  lineColor <- "black"

  if (results[["probable.model"]] == "Hyperbolic") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = hypSeries)
    legend = c(legend, paste("Hyperbolic(",
                             round(results[["Hyperbolic.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Exponential") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = expSeries)
    legend = c(legend, paste("Exponential(",
                             round(results[["Exponential.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Laibson") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = quaSeries)
    legend = c(legend, paste("Laibson(",
                             round(results[["Laibson.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "GreenMyerson") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = myerSeries)
    legend = c(legend, paste("GreenMyerson(",
                             round(results[["GreenMyerson.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "Rachlin") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = rachSeries)
    legend = c(legend, paste("Rachlin(",
                             round(results[["Rachlin.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "EbertPrelec") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = epSeries)
    legend = c(legend, paste("EbertPrelec(",
                             round(results[["EbertPrelec.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = crdiSeries)
    legend = c(legend, paste("BleichrodtCRDI(",
                             round(results[["BleichrodtCRDI.prob"]], 5),
                             ")",
                             sep = ""))

  } else if (results[["probable.model"]] == "RodriguezLogue") {
    totalFrame = data.frame(Delays = delaySeries,
                            ModelArea = ghSeries)
    legend = c(legend, paste("RodriguezLogue(",
                             round(results[["RodriguezLogue.prob"]], 5),
                             ")",
                             sep = ""))

  } else {
    totalFrame = data.frame(Delays = delaySeries)
    totalFrame$ModelArea <- results[["Noise.mean"]]

    legend = c(legend, paste("Noise(",
                             round(results[["Noise.prob"]], 5),
                             ")",
                             sep = ""))

  }

  graphics::plot(totalFrame$Delays, totalFrame$ModelArea, type = "l", ylim = c(0,1),
                 main = paste("Participant ", dat$id[1],
                              "\n", results[["probable.model"]], " Scaled Model Area = ", round(results[["probable.Log10AUC"]], 5),
                              sep = ""),
                 xlab = "Delays",
                 ylab = "Value",
                 col = lineColor,
                 lwd = lineWidth,
                 log = "x")

  graphics::points(dat$X,
                   dat$Y,
                   type = "p",
                   cex = 2,
                   pch = 18)

  graphics::lines(dat$X,
                  dat$Y,
                  col = "black",
                  lty = 2,
                  lwd = 0.5)

  legend("topright",
         legend = legend,
         col = colors,
         lty = c(2, 1),
         lwd = 3,
         title = "Model (Probability)")
}
