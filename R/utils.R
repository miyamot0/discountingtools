#' minpack.lm logLik hack
#'
#' This function constructs a class, derived from an nls.lm object, similar to that of the logLik function in nls. This allows for native calls of the AIC and BIC functions from stats, using nls.lm fit objects.
#'
#' @param fit nls.lm fitted model
#' @param REML determine whether or not to use ML (FALSE by default)
#' @param ... inherit other args as necessary
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

#' summary.discountingtools
#'
#' Override summary output
#'
#' @param fittingObject
#'
#' @return
#' @export summary.discountingtools
#' @export
summary.discountingtools <- function(fittingObject) {

  localCopy <- fittingObject$results

  buildColNames = c("ID")

  # TODO: build out additional models

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

  nRows    = length(names(localCopy))
  resFrame = data.frame(matrix(ncol = length(buildColNames),
                               nrow = nRows))

  colnames(resFrame) <- buildColNames

  resFrame$ID <- names(localCopy)

  ### Load results

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
  }

  resFrame
}
