#' summary.discountingtools
#'
#' Override summary output. Rather than display the core fitting object, a data frame block of results is provided to the user for easy interpretation and further analysis
#'
#' @param fittingObject core fitting object
#' @param detailed enable additional model metrics (default TRUE)
#'
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @export summary.discountingtools
#' @export
summary.discountingtools <- function(fittingObject, detailed = TRUE) {

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

    resFrame[index, "Strategy"] = fittingObject[[ "strategy" ]]

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
