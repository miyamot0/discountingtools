#' summary.discountingtools
#'
#' Override summary output. Rather than display the core fitting object, a data frame block of results is provided to the user for easy interpretation and further analysis
#'
#' @param fittingObject core fitting object
#' @param detailed enable additional model metrics (default FALSE)
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export summary.discountingtools
#' @export
summary.discountingtools <- function(fittingObject, detailed = FALSE) {

  localCopy <- fittingObject$results

  buildColNames = c("ID", "Strategy")

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

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Noise.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Noise.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Noise.Log10MBAUC"))
      }
    } else if (m == "mazur") {
      buildColNames = c(buildColNames,
                        "Mazur.Lnk",
                        "Mazur.RMSE",
                        "Mazur.BIC",
                        "Mazur.AIC",
                        "Mazur.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Mazur.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Mazur.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Mazur.Log10MBAUC"))
      }
    } else if (m == "exponential") {
      buildColNames = c(buildColNames,
                        "Exponential.Lnk",
                        "Exponential.RMSE",
                        "Exponential.BIC",
                        "Exponential.AIC",
                        "Exponential.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Exponential.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Exponential.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Exponential.Log10MBAUC"))
      }
    } else if (m == "laibson") {
      buildColNames = c(buildColNames,
                        "Laibson.Beta",
                        "Laibson.Delta",
                        "Laibson.RMSE",
                        "Laibson.BIC",
                        "Laibson.AIC",
                        "Laibson.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Laibson.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Laibson.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Laibson.Log10MBAUC"))
      }
    } else if (m == "greenmyerson") {
      buildColNames = c(buildColNames,
                        "GreenMyerson.Lnk",
                        "GreenMyerson.S",
                        "GreenMyerson.RMSE",
                        "GreenMyerson.BIC",
                        "GreenMyerson.AIC",
                        "GreenMyerson.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("GreenMyerson.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("GreenMyerson.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("GreenMyerson.Log10MBAUC"))
      }
    } else if (m == "rachlin") {
      buildColNames = c(buildColNames,
                        "Rachlin.Lnk",
                        "Rachlin.S",
                        "Rachlin.RMSE",
                        "Rachlin.BIC",
                        "Rachlin.AIC",
                        "Rachlin.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Rachlin.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Rachlin.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Rachlin.Log10MBAUC"))
      }
    } else if (m == "ebertprelec") {
      buildColNames = c(buildColNames,
                        "EbertPrelec.Lnk",
                        "EbertPrelec.S",
                        "EbertPrelec.RMSE",
                        "EbertPrelec.BIC",
                        "EbertPrelec.AIC",
                        "EbertPrelec.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("EbertPrelec.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("EbertPrelec.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("EbertPrelec.Log10MBAUC"))
      }
    } else if (m == "bleichrodt") {
      buildColNames = c(buildColNames,
                        "Bleichrodt.Lnk",
                        "Bleichrodt.S",
                        "Bleichrodt.Beta",
                        "Bleichrodt.RMSE",
                        "Bleichrodt.BIC",
                        "Bleichrodt.AIC",
                        "Bleichrodt.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("Bleichrodt.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("Bleichrodt.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("Bleichrodt.Log10MBAUC"))
      }
    } else if (m == "rodriguezlogue") {
      buildColNames = c(buildColNames,
                        "RodriguezLogue.Lnk",
                        "RodriguezLogue.Beta",
                        "RodriguezLogue.RMSE",
                        "RodriguezLogue.BIC",
                        "RodriguezLogue.AIC",
                        "RodriguezLogue.Status")

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   buildColNames = append(buildColNames, c("RodriguezLogue.LnED50"))
        if (metric == "mbauc")    buildColNames = append(buildColNames, c("RodriguezLogue.MBAUC"))
        if (metric == "logmbauc") buildColNames = append(buildColNames, c("RodriguezLogue.Log10MBAUC"))
      }
    }
  }

  if (fittingObject$ModelSelection == TRUE) {
    buildColNames = append(buildColNames, c("ProbableModel",
                                            "ProbableModel.BF",
                                            "ProbableModel.Prob"))

    for (metric in fittingObject[["metrics"]]) {
      if (metric == "lned50")   buildColNames = append(buildColNames, c("ProbableModel.LnED50"))
      if (metric == "mbauc")    buildColNames = append(buildColNames, c("ProbableModel.MBAUC"))
      if (metric == "logmbauc") buildColNames = append(buildColNames, c("ProbableModel.Log10MBAUC"))
    }
  }

  nRows    = length(names(localCopy))
  resFrame = data.frame(matrix(ncol = length(buildColNames),
                               nrow = nRows))

  colnames(resFrame) <- buildColNames

  resFrame$ID <- names(localCopy)

  for (name in names(localCopy)) {
    index = which(names(localCopy) == name)

    resFrame[index, "Strategy"] = fittingObject[[ "strategy" ]]

    for (res in localCopy[[name]]) {

      if (res$Model == "noise") {
        resFrame[index, c("Noise.Intercept",
                          "Noise.RMSE",
                          "Noise.BIC",
                          "Noise.AIC",
                          "Noise.LnED50",
                          "Noise.MBAUC",
                          "Noise.Log10MBAUC")] = as.data.frame(res)[, c("Intercept",
                                                                        "RMSE",
                                                                        "BIC",
                                                                        "AIC",
                                                                        "ED50",
                                                                        "MBAUC",
                                                                        "Log10MBAUC")]

      } else if (res$Model == "mazur") {
        resFrame[index, c("Mazur.Lnk",
                          "Mazur.RMSE",
                          "Mazur.BIC",
                          "Mazur.AIC",
                          "Mazur.Status",
                          "Mazur.LnED50",
                          "Mazur.MBAUC",
                          "Mazur.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                    "RMSE",
                                                                    "BIC",
                                                                    "AIC",
                                                                    "Status",
                                                                    "ED50",
                                                                    "MBAUC",
                                                                    "Log10MBAUC")]
      } else if (res$Model == "exponential") {
        resFrame[index, c("Exponential.Lnk",
                          "Exponential.RMSE",
                          "Exponential.BIC",
                          "Exponential.AIC",
                          "Exponential.Status",
                          "Exponential.LnED50",
                          "Exponential.MBAUC",
                          "Exponential.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                              "RMSE",
                                                                              "BIC",
                                                                              "AIC",
                                                                              "Status",
                                                                              "ED50",
                                                                              "MBAUC",
                                                                              "Log10MBAUC")]
      } else if (res$Model == "laibson") {
        resFrame[index, c("Laibson.Beta",
                          "Laibson.Delta",
                          "Laibson.RMSE",
                          "Laibson.BIC",
                          "Laibson.AIC",
                          "Laibson.Status",
                          "Laibson.LnED50",
                          "Laibson.MBAUC",
                          "Laibson.Log10MBAUC")] = as.data.frame(res)[, c("Beta",
                                                                      "Delta",
                                                                      "RMSE",
                                                                      "BIC",
                                                                      "AIC",
                                                                      "Status",
                                                                      "ED50",
                                                                      "MBAUC",
                                                                      "Log10MBAUC")]
      } else if (res$Model == "greenmyerson") {
        resFrame[index, c("GreenMyerson.Lnk",
                          "GreenMyerson.S",
                          "GreenMyerson.RMSE",
                          "GreenMyerson.BIC",
                          "GreenMyerson.AIC",
                          "GreenMyerson.Status",
                          "GreenMyerson.LnED50",
                          "GreenMyerson.MBAUC",
                          "GreenMyerson.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                           "S",
                                                                           "RMSE",
                                                                           "BIC",
                                                                           "AIC",
                                                                           "Status",
                                                                           "ED50",
                                                                           "MBAUC",
                                                                           "Log10MBAUC")]
      } else if (res$Model == "rachlin") {
        resFrame[index, c("Rachlin.Lnk",
                          "Rachlin.S",
                          "Rachlin.RMSE",
                          "Rachlin.BIC",
                          "Rachlin.AIC",
                          "Rachlin.Status",
                          "Rachlin.LnED50",
                          "Rachlin.MBAUC",
                          "Rachlin.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                      "S",
                                                                      "RMSE",
                                                                      "BIC",
                                                                      "AIC",
                                                                      "Status",
                                                                      "ED50",
                                                                      "MBAUC",
                                                                      "Log10MBAUC")]
      } else if (res$Model == "ebertprelec") {
        resFrame[index, c("EbertPrelec.Lnk",
                          "EbertPrelec.S",
                          "EbertPrelec.RMSE",
                          "EbertPrelec.BIC",
                          "EbertPrelec.AIC",
                          "EbertPrelec.Status",
                          "EbertPrelec.LnED50",
                          "EbertPrelec.MBAUC",
                          "EbertPrelec.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                          "S",
                                                                          "RMSE",
                                                                          "BIC",
                                                                          "AIC",
                                                                          "Status",
                                                                          "ED50",
                                                                          "MBAUC",
                                                                          "Log10MBAUC")]
      } else if (res$Model == "bleichrodt") {
        resFrame[index, c("Bleichrodt.Lnk",
                          "Bleichrodt.S",
                          "Bleichrodt.Beta",
                          "Bleichrodt.RMSE",
                          "Bleichrodt.BIC",
                          "Bleichrodt.AIC",
                          "Bleichrodt.Status",
                          "Bleichrodt.LnED50",
                          "Bleichrodt.MBAUC",
                          "Bleichrodt.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                         "S",
                                                                         "Beta",
                                                                         "RMSE",
                                                                         "BIC",
                                                                         "AIC",
                                                                         "Status",
                                                                         "ED50",
                                                                         "MBAUC",
                                                                         "Log10MBAUC")]
      } else if (res$Model == "rodriguezlogue") {
        resFrame[index, c("RodriguezLogue.Lnk",
                          "RodriguezLogue.Beta",
                          "RodriguezLogue.RMSE",
                          "RodriguezLogue.BIC",
                          "RodriguezLogue.AIC",
                          "RodriguezLogue.Status",
                          "RodriguezLogue.LnED50",
                          "RodriguezLogue.MBAUC",
                          "RodriguezLogue.Log10MBAUC")] = as.data.frame(res)[, c("Lnk",
                                                                             "Beta",
                                                                             "RMSE",
                                                                             "BIC",
                                                                             "AIC",
                                                                             "Status",
                                                                             "ED50",
                                                                             "MBAUC",
                                                                             "Log10MBAUC")]
      }
    }

    if (fittingObject$ModelSelection == TRUE) {
      resFrame[index, "ProbableModel"]      = fittingObject$rotation[[name]]$ProbableModel
      resFrame[index, "ProbableModel.BF"]   = fittingObject$rotation[[name]]$ProbableModel.BF
      resFrame[index, "ProbableModel.Prob"] = fittingObject$rotation[[name]]$ProbableModel.Prob

      for (metric in fittingObject[["metrics"]]) {
        if (metric == "lned50")   resFrame[index, "ProbableModel.LnED50"]     = fittingObject$ed50[[name]]
        if (metric == "mbauc")    resFrame[index, "ProbableModel.MBAUC"]      = fittingObject$mbauc[[name]]
        if (metric == "logmbauc") resFrame[index, "ProbableModel.Log10MBAUC"] = fittingObject$mbauclog10[[name]]
      }
    }
  }

  if (detailed == FALSE) {
    resFrame = resFrame[,!grepl(".RMSE",   colnames(resFrame))]
    resFrame = resFrame[,!grepl(".AIC",    colnames(resFrame))]
    resFrame = resFrame[,!grepl(".BIC",    colnames(resFrame))]
    resFrame = resFrame[,!grepl(".Status", colnames(resFrame))]
    resFrame = resFrame[,!grepl(".BF",     colnames(resFrame))]
    resFrame = resFrame[,!grepl(".Prob",   colnames(resFrame))]
  }

  resFrame
}
