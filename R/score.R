#' Perform Discounting Model Selection
#'
#' This function takes a data frame of temporal discounting values (X, Y) and performs approximate Bayesian model selection using the Bayesian Information Criterion and returns the log of the Effective Delay 50 and numerical integration area.
#'
#' @param dat data frame with X column and Y column (0 <= Y <= 1)
#' @param A OPTIONAL: Modify line sizes for figures
#' @param models vector of models to include in selection
#' @param detailed include additional output details in results
#' @param figures OPTIONAL: show figure of results
#' @param summarize OPTIONAL: Descriptive observations
#' @param lineSize OPTIONAL: Modify line sizes for figures
#' @return A data frame of model parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @importFrom stats lm nls
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @return data frame of fitted model parameters
discountingModelSelectionCall <- function(dat, A = NULL, models = c("noise"), detailed = FALSE, figures = FALSE, summarize = FALSE, lineSize = 1) {

  lengthX <- length(dat$X)

  # List to return
  returnList <- list()

  # Lists to carry
  bicList <- list()

  # Re-used references
  tempList <- NULL
  projection <- NULL
  sqResidual <- NULL
  sumSquares <- NULL
  presort <- NULL
  sorted <- NULL
  ini.par <- NULL

  modelFitNoise <- NULL

  try(modelFitNoise <- stats::lm(Y ~ 1, dat),
      silent=TRUE)

  if(!is.character(modelFitNoise)) {

    if (detailed == TRUE) {
      tempList <- list(Noise.mean = modelFitNoise$coefficients[["(Intercept)"]],
                       Noise.RMSE = summary(modelFitNoise)[["sigma"]],
                       Noise.BIC  = ifelse(summary(modelFitNoise)[["sigma"]] == 0, Inf, stats::BIC(modelFitNoise)),
                       Noise.AIC  = ifelse(summary(modelFitNoise)[["sigma"]] == 0, Inf, stats::AIC(modelFitNoise)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(Noise.BIC = tempList$Noise.BIC))

    } else {
      tempList <- list(Noise.mean = modelFitNoise$coefficients[["(Intercept)"]])

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(Noise.BIC  = ifelse(summary(modelFitNoise)[["sigma"]] == 0, Inf, stats::BIC(modelFitNoise))))
    }
  }

  if ('exponential' %in% models) {
    # Starts
    startlnK  <- seq(-12,  12, 1)
    lengthLnK <- length(startlnK)

    # Params expand
    MlnK <- sort(rep(startlnK, lengthX))

    # Init SS vector
    sumSquares<- rep(NA,lengthLnK)

    # Observed Data
    MX <- rep(dat$X, lengthLnK)
    MY <- rep(dat$Y, lengthLnK)

    # Projections
    projection<- exp(-exp(MlnK)*MX)
    sqResidual<- (MY-projection)^2

    for(j in 1:lengthLnK){
      sumSquares[j]  <-sum(sqResidual[(j-1)*lengthX +1:lengthX])

    }

    # Sort starting estimates
    presort   <-data.frame(startlnK,
                           sumSquares)
    sorted    <-presort[order(presort[ ,"sumSquares"]), ]
    ini.par   <-c(lnk=sorted$startlnK[1])

    modelFitExponential <- NULL

    try(modelFitExponential<-nls.lm(par = ini.par,
                                    fn = residualFunction,
                                    jac = jacobianMatrix,
                                    valueFunction = exponentialDiscountFunc,
                                    jacobianFunction = exponentialDiscountGradient,
                                    x = dat$X,
                                    value = dat$Y,
                                    control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitExponential)) {

      if (detailed == TRUE) {
        tempList <- list(Exponential.lnk = modelFitExponential$par[["lnk"]],
                         Exponential.RMSE = sqrt(modelFitExponential$deviance/length(modelFitExponential$fvec)),
                         Exponential.BIC  = stats::BIC(logLik.nls.lm(modelFitExponential)),
                         Exponential.AIC  = stats::AIC(logLik.nls.lm(modelFitExponential)),
                         Exponential.status = paste("Code:",
                                            modelFitExponential$info,
                                            "- Message:",
                                            modelFitExponential$message,
                                            sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Exponential.BIC = tempList$Exponential.BIC))

      } else {
        tempList <- list(Exponential.lnk = modelFitExponential$par[["lnk"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Exponential.BIC = stats::BIC(logLik.nls.lm(modelFitExponential))))

      }

    }
  }

  if ('hyperbolic' %in% models) {
    # Starts
    startlnK <- seq(-12,  12, 1)
    lengthLnK <- length(startlnK)

    # Params expand
    MlnK <- sort(rep(startlnK, lengthX))

    # Init SS vector
    sumSquares <- rep(NA,lengthLnK)

    # Observed Data
    MX <- rep(dat$X, lengthLnK)
    MY <- rep(dat$Y, lengthLnK)

    # Projections
    projection <- (1+exp(MlnK)*MX)^(-1)
    sqResidual <- (MY-projection)^2

    for(j in 1:lengthLnK){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

    }

    # Sort starting estimates
    presort <-data.frame(startlnK,
                         sumSquares)
    sorted  <-presort[order(presort[ ,"sumSquares"]), ]
    ini.par <-c(lnk=sorted$startlnK[1])

    modelFitHyperbolic <- NULL

    try(modelFitHyperbolic <- nls.lm(par = ini.par,
                                     fn = residualFunction,
                                     jac = jacobianMatrix,
                                     valueFunction = hyperbolicDiscountFunc,
                                     jacobianFunction = hyperbolicDiscountGradient,
                                     x = dat$X,
                                     value = dat$Y,
                                     control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitHyperbolic)) {

      if (detailed == TRUE) {
        tempList <- list(Hyperbolic.lnk  = modelFitHyperbolic$par[["lnk"]],
                         Hyperbolic.RMSE = sqrt(modelFitHyperbolic$deviance/length(modelFitHyperbolic$fvec)),
                         Hyperbolic.BIC  = stats::BIC(logLik.nls.lm(modelFitHyperbolic)),
                         Hyperbolic.AIC  = stats::AIC(logLik.nls.lm(modelFitHyperbolic)),
                         Hyperbolic.status = paste("Code:",
                                              modelFitHyperbolic$info,
                                              "- Message:",
                                              modelFitHyperbolic$message,
                                              sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Hyperbolic.BIC = tempList$Hyperbolic.BIC))

      } else {
        tempList <- list(Hyperbolic.lnk  = modelFitHyperbolic$par[["lnk"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Hyperbolic.BIC = stats::BIC(logLik.nls.lm(modelFitHyperbolic))))

      }
    }
  }

  if ('laibson' %in% models) {
    # Starts
    startbeta   <- seq(0, 1, 0.1)
    startdelta  <- seq(0, 1, 0.01)

    lengthBeta  <- length(startbeta)
    lengthDelta <- length(startdelta)

    # Params expand
    startBeta   <- rep(sort(rep(startbeta, lengthX)),lengthDelta)
    startdelta  <- sort(rep(startdelta,lengthX*lengthBeta))

    # Init SS vector
    sumSquares  <- rep(NA,lengthBeta*lengthDelta)

    # Observed Data
    SY          <- rep(dat$Y,lengthBeta*lengthDelta)

    # Projections
    projection  <- startBeta*startdelta^dat$X
    sqResidual  <- (SY-projection)^2

    SSbeta      <- rep(startbeta,lengthDelta)
    SSdelta     <- sort(rep(startdelta,lengthBeta))

    for(j in 1:(lengthBeta*lengthDelta)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    # Sort starting estimates
    presort<-data.frame(SSbeta,
                        SSdelta,
                        sumSquares)

    sorted  <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(beta  = sorted$SSbeta[1],
                 delta = sorted$SSdelta[1])

    modelFitBetaDelta <- NULL

    try(modelFitBetaDelta <- nls.lm(par = ini.par,
                                    fn = residualFunction,
                                    jac = jacobianMatrix,
                                    valueFunction = betaDeltaDiscountFunc,
                                    jacobianFunction = betaDeltaDiscountGradient,
                                    x = dat$X,
                                    value = dat$Y,
                                    upper = c(beta = 1, delta = 1),
                                    lower = c(beta = 0, delta = 0),
                                    control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitBetaDelta)) {

      if (detailed == TRUE) {
        tempList <- list(Laibson.beta  = modelFitBetaDelta$par[["beta"]],
                         Laibson.delta  = modelFitBetaDelta$par[["delta"]],
                         Laibson.RMSE = sqrt(modelFitBetaDelta$deviance/length(modelFitBetaDelta$fvec)),
                         Laibson.BIC  = stats::BIC(logLik.nls.lm(modelFitBetaDelta)),
                         Laibson.AIC  = stats::AIC(logLik.nls.lm(modelFitBetaDelta)),
                         Laibson.status = paste("Code:",
                                           modelFitBetaDelta$info,
                                           "- Message:",
                                           modelFitBetaDelta$message,
                                           sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Laibson.BIC = tempList$Laibson.BIC))

      } else {
        tempList <- list(Laibson.beta  = modelFitBetaDelta$par[["beta"]],
                         Laibson.delta  = modelFitBetaDelta$par[["delta"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Laibson.BIC = stats::BIC(logLik.nls.lm(modelFitBetaDelta))))

      }
    }
  }

  if ('greenmyerson' %in% models) {
    # Starts
    startlnK <- seq(-12, 12, 1)
    starts <- seq(.01, 10, 0.01)

    lengthLnK <- length(startlnK)
    lengthS <- length(starts)

    # Params expand
    SSlnK <- rep(startlnK,lengthS)
    SSs <- sort(rep(starts,lengthLnK))

    # Init SS vector
    sumSquares <- rep(NA,lengthS*lengthLnK)

    # Observed Data
    SY <- rep(dat$Y,lengthS*lengthLnK)

    # Projections
    SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss <- sort(rep(starts,lengthX*lengthLnK))

    projection <- (1+exp(SlnK)*dat$X)^(-Ss)
    sqResidual <- (SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

    }

    # Sort starting estimates
    presort <- data.frame(SSlnK,SSs,sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    modelFitMyerson <- NULL

    try(modelFitMyerson <- nls.lm(par = ini.par,
                                  fn = residualFunction,
                                  jac = jacobianMatrix,
                                  valueFunction = myersonHyperboloidDiscountFunc,
                                  jacobianFunction = myersonHyperboloidDiscountGradient,
                                  x = dat$X,
                                  value = dat$Y,
                                  control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitMyerson)) {

      if (detailed == TRUE) {
        tempList <- list(GreenMyerson.lnk  = modelFitMyerson$par[["lnk"]],
                         GreenMyerson.s  = modelFitMyerson$par[["s"]],
                         GreenMyerson.RMSE = sqrt(modelFitMyerson$deviance/length(modelFitMyerson$fvec)),
                         GreenMyerson.BIC  = stats::BIC(logLik.nls.lm(modelFitMyerson)),
                         GreenMyerson.AIC  = stats::AIC(logLik.nls.lm(modelFitMyerson)),
                         GreenMyerson.status = paste("Code:",
                                           modelFitMyerson$info,
                                           "- Message:",
                                           modelFitMyerson$message,
                                           sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(GreenMyerson.BIC = tempList$GreenMyerson.BIC))

      } else {
        tempList <- list(GreenMyerson.lnk  = modelFitMyerson$par[["lnk"]],
                         GreenMyerson.s  = modelFitMyerson$par[["s"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(GreenMyerson.BIC = stats::BIC(logLik.nls.lm(modelFitMyerson))))

      }
    }
  }

  if ('rachlin' %in% models) {
    # Starts
    startlnK<-seq(-12,12,1)
    starts<-seq(.01,10,.01)

    lengthLnK<-length(startlnK)
    lengthS<-length(starts)

    # Params expand
    SSlnK<-rep(startlnK,lengthS)
    SSs<-sort(rep(starts,lengthLnK))

    # Init SS vector
    sumSquares<-rep(NA,lengthS*lengthLnK)

    # Observed Data
    SY<-rep(dat$Y,lengthS*lengthLnK)

    # Projections
    SlnK<-rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss<-sort(rep(starts,lengthX*lengthLnK))

    projection<-(1+exp(SlnK)*(dat$X^Ss))^(-1)
    sqResidual<-(SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

    }

    # Sort starting estimates
    presort <- data.frame(SSlnK,
                          SSs,
                          sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    modelFitRachlin <- NULL

    try(modelFitRachlin <- nls.lm(par = ini.par,
                                  fn = residualFunction,
                                  jac = jacobianMatrix,
                                  valueFunction = rachlinHyperboloidDiscountFunc,
                                  jacobianFunction = rachlinHyperboloidDiscountGradient,
                                  x = dat$X,
                                  value = dat$Y,
                                  control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitRachlin)) {

      if (detailed == TRUE) {
        tempList <- list(Rachlin.lnk  = modelFitRachlin$par[["lnk"]],
                         Rachlin.s  = modelFitRachlin$par[["s"]],
                         Rachlin.RMSE = sqrt(modelFitRachlin$deviance/length(modelFitRachlin$fvec)),
                         Rachlin.BIC  = stats::BIC(logLik.nls.lm(modelFitRachlin)),
                         Rachlin.AIC  = stats::AIC(logLik.nls.lm(modelFitRachlin)),
                         Rachlin.status = paste("Code:",
                                                modelFitRachlin$info,
                                                "- Message:",
                                                modelFitRachlin$message,
                                                sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Rachlin.BIC = tempList$Rachlin.BIC))

      } else {
        tempList <- list(Rachlin.lnk  = modelFitRachlin$par[["lnk"]],
                         Rachlin.s  = modelFitRachlin$par[["s"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(Rachlin.BIC = stats::BIC(logLik.nls.lm(modelFitRachlin))))

      }
    }
  }

  if ('ebertprelec' %in% models) {
    # Starts
    startlnK <- seq(-12, 12, 0.1)
    starts <- seq(.01, 1, 0.01)

    lengthLnK <- length(startlnK)
    lengthS <- length(starts)

    # Params expand
    SSlnK <- rep(startlnK,lengthS)
    SSs <- sort(rep(starts,lengthLnK))

    # Init SS vector
    sumSquares <- rep(NA,lengthS*lengthLnK)

    # Observed Data
    SY <- rep(dat$Y,lengthS*lengthLnK)

    # Projections
    SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss <- sort(rep(starts,lengthX*lengthLnK))

    projection <- exp(-(exp(SlnK)*dat$X)^Ss)
    sqResidual <- (SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

    }

    # Sort starting estimates
    presort <- data.frame(SSlnK,
                          SSs,
                          sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    modelFitep <- NULL

    try(modelFitep <- nls.lm(par = ini.par,
                             fn = residualFunction,
                             jac = jacobianMatrix,
                             valueFunction = ebertPrelecDiscountFunc,
                             jacobianFunction = ebertPrelecDiscountGradient,
                             x = dat$X,
                             value = dat$Y,
                             control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitep)) {

      if (detailed == TRUE) {
        tempList <- list(EbertPrelec.lnk  = modelFitep$par[["lnk"]],
                         EbertPrelec.s  = modelFitep$par[["s"]],
                         EbertPrelec.RMSE = sqrt(modelFitep$deviance/length(modelFitep$fvec)),
                         EbertPrelec.BIC  = stats::BIC(logLik.nls.lm(modelFitep)),
                         EbertPrelec.AIC  = stats::AIC(logLik.nls.lm(modelFitep)),
                         EbertPrelec.status = paste("Code:",
                                           modelFitep$info,
                                           "- Message:",
                                           modelFitep$message,
                                           sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(EbertPrelec.BIC = tempList$EbertPrelec.BIC))

      } else {
        tempList <- list(EbertPrelec.lnk  = modelFitep$par[["lnk"]],
                         EbertPrelec.s  = modelFitep$par[["s"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(EbertPrelec.BIC = stats::BIC(logLik.nls.lm(modelFitep))))

      }
    }
  }

  if ('crdi' %in% models) {
    # Starts
    startlnK <- seq(-12, 12, 1)
    starts <- seq(.01, 1, 0.1)
    startBeta <- seq(0.1, 1, 0.1)

    # new Pre sort
    presort <- expand.grid(startlnK = startlnK, starts = starts, startBeta = startBeta)
    presort$sumSquares <- NA

    # clean, merge, or move
    getSS <- function(presort, index, Y, X) {
      projections <- presort[index,]$startBeta*exp(-(exp(presort[index,]$startlnK)*X^presort[index,]$starts))
      sqResidual <- (Y - projections)^2
      sum(sqResidual)
    }

    for (j in 1:nrow(presort)) {
      presort[j, ]$sumSquares <-getSS(presort, j, dat$Y, dat$X)
    }

    presort <- presort[order(presort$sumSquares),]

    ini.par <- c(beta = presort[1,]$startBeta, lnk = presort[1,]$startlnK, s = presort[1,]$starts)

    modelFitCRDI <- NULL

    try(modelFitCRDI<-nls.lm(par = ini.par,
                             fn = residualFunction,
                             jac = jacobianMatrix,
                             valueFunction = BleichrodtCRDIDiscountFunc,
                             jacobianFunction = BleichrodtCRDIDiscountGradient,
                             x = dat$X,
                             value = dat$Y,
                             upper = c(beta = 1, lnk = Inf, s = Inf),
                             lower = c(beta = 0, lnk = -Inf, s = -Inf),
                             control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitCRDI)) {

      if (detailed == TRUE) {
        tempList <- list(BleichrodtCRDI.lnk  = modelFitCRDI$par[["lnk"]],
                         BleichrodtCRDI.s    = modelFitCRDI$par[["s"]],
                         BleichrodtCRDI.beta = modelFitCRDI$par[["beta"]],
                         BleichrodtCRDI.RMSE = sqrt(modelFitCRDI$deviance/length(modelFitCRDI$fvec)),
                         BleichrodtCRDI.BIC  = stats::BIC(logLik.nls.lm(modelFitCRDI)),
                         BleichrodtCRDI.AIC  = stats::AIC(logLik.nls.lm(modelFitCRDI)),
                         BleichrodtCRDI.status = paste("Code:",
                                                       modelFitCRDI$info,
                                                       "- Message:",
                                                       modelFitCRDI$message,
                                                       sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(BleichrodtCRDI.BIC = tempList$BleichrodtCRDI.BIC))

      } else {
        tempList <- list(BleichrodtCRDI.lnk  = modelFitCRDI$par[["lnk"]],
                         BleichrodtCRDI.s    = modelFitCRDI$par[["s"]],
                         BleichrodtCRDI.beta = modelFitCRDI$par[["beta"]])

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(BleichrodtCRDI.BIC = stats::BIC(logLik.nls.lm(modelFitCRDI))))

      }
    }
  }

  if ('genhyp' %in% models) {
    # Starts
    startlnK <- seq(-12, 12, 1)
    startBeta <- seq(-12, 12, 1)

    # new Pre sort
    presort <- expand.grid(startlnK = startlnK, startBeta = startBeta)
    presort$sumSquares <- NA

    # clean, merge, or move
    getSS <- function(presort, index, Y, X) {
      projections <- (1 + X * exp(presort[index,]$startlnK))^(-exp(presort[index,]$startBeta) / exp(presort[index,]$startlnK))
      sqResidual <- (Y - projections)^2
      sum(sqResidual)
    }

    for (j in 1:nrow(presort)) {
      presort[j, ]$sumSquares <-getSS(presort, j, dat$Y, dat$X)
    }

    presort <- presort[order(presort$sumSquares),]

    ini.par <- c(lnk = presort[1,]$startlnK, beta = presort[1,]$startBeta)

    modelFitGenHyp <- NULL

    try(modelFitGenHyp<-nls.lm(par = ini.par,
                               fn = residualFunction,
                               jac = jacobianMatrix,
                               valueFunction = RodriguezLogueDiscountFunc,
                               jacobianFunction = RodriguezLogueDiscountGradient,
                               x = dat$X,
                               value = dat$Y,
                               control = nls.lm.control(maxiter = 1000)), silent = TRUE)

    if (!is.character(modelFitGenHyp)) {

      if (detailed == TRUE) {
        tempList <- list(RodriguezLogue.lnk  = modelFitGenHyp$par[["lnk"]],
                         RodriguezLogue.beta = exp(modelFitGenHyp$par[["beta"]]),
                         RodriguezLogue.RMSE = sqrt(modelFitGenHyp$deviance/length(modelFitGenHyp$fvec)),
                         RodriguezLogue.BIC  = stats::BIC(logLik.nls.lm(modelFitGenHyp)),
                         RodriguezLogue.AIC  = stats::AIC(logLik.nls.lm(modelFitGenHyp)),
                         RodriguezLogue.status = paste("Code:",
                                                       modelFitGenHyp$info,
                                                       "- Message:",
                                                       modelFitGenHyp$message,
                                                       sep = " "))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(RodriguezLogue.BIC = tempList$RodriguezLogue.BIC))

      } else {
        tempList <- list(RodriguezLogue.lnk  = modelFitGenHyp$par[["lnk"]],
                         RodriguezLogue.beta = exp(modelFitGenHyp$par[["beta"]]))

        returnList <- c(returnList, tempList)

        bicList <- c(bicList, list(RodriguezLogue.BIC = stats::BIC(logLik.nls.lm(modelFitGenHyp))))

      }
    }
  }

  ### bfs
  bfList <- list()
  bfSum <- 0.0

  ### probs
  probList <- list()

  # Perfect fit, hacky workaround
  if (bicList[["Noise.BIC"]] == Inf) {
    for (i in 1:length(bicList)) {
      bfName = gsub("BIC", "BF", names(bicList)[i])
      bfList[[bfName]] = 0

      bfSum <- bfSum + bfList[[bfName]]
    }

    if (detailed == TRUE) {
      returnList <- c(returnList, bfList)
    }

    for (i in 1:length(bfList)) {
      probName = gsub("BF", "prob", names(bfList)[i])
      probList[[probName]] = 0
    }

    # Hack
    probList[["Noise.prob"]] <- 1

  } else {
    for (i in 1:length(bicList)) {
      bfName = gsub("BIC", "BF", names(bicList)[i])
      bfList[[bfName]] = exp(-.5*(bicList[[names(bicList)[i]]]-bicList[["Noise.BIC"]]))

      bfSum <- bfSum + bfList[[bfName]]

    }

    if (detailed == TRUE) {
      returnList <- c(returnList, bfList)
    }

    for (i in 1:length(bfList)) {
      probName = gsub("BF", "prob", names(bfList)[i])
      probList[[probName]] = bfList[[names(bfList)[i]]]/bfSum

    }
  }


  returnList <- c(returnList, probList)

  mostProb <- names(probList[which.max(probList)])
  probableModel = gsub(".prob", "", mostProb)
  returnList <- c(returnList, list(probable.model=probableModel))

  probableED50 <- getED50(dat, returnList)
  probableAUC <- getModelAUC(dat, returnList)
  probableAUCLog <- getModelAUCLog10Scaled(dat, returnList)
  returnList <- c(returnList, list(probable.ED50=probableED50,
                                   probable.AUC=probableAUC,
                                   probable.Log10AUC=probableAUCLog))

  ### Summarize here
  if (summarize == TRUE) {
    message("Model Comparison Summary:")

    bicList  <- bicList[order(unlist(bicList), decreasing=FALSE)]
    bfList   <- bfList[order(unlist(bfList), decreasing=FALSE)]
    probList <- probList[order(unlist(probList), decreasing=TRUE)]

    for (i in 1:length(probList)) {
      message(paste("Rank #: ", i, " = ", gsub(".prob", "", names(probList)[i])));
      message(paste(" Probability", " = ", round(as.numeric(probList[i]), 6)));
      message(paste(" BIC", " = ", round(as.numeric(bicList[i]), 6)));

      if (i == 1) {
        message(paste(" Most Probable ln(ED50)", " = ", round(probableED50, 6)));
        message(paste(" Most Probable Model AUC", " = ", round(probableAUC, 6)));
        message(paste(" Most Probable Model AUC (log)", " = ", round(probableAUCLog, 6)));

      }

      message("")
    }
  }

  ### Plotting here
  if (figures == "ed50") {
    suppressWarnings(displayED50Figure(dat, returnList, lineSize))

  } else if (figures == "auc") {
    suppressWarnings(displayAUCFigure(dat, returnList, lineSize))

  } else if (figures == "logauc") {
    suppressWarnings(displayLogAUCFigure(dat, returnList, lineSize))

  }

  returnFrame <- as.data.frame(returnList)
  returnFrame
}

#' Perform Discounting Model Selection
#'
#' This function takes a data frame of temporal discounting values (X, Y) and performs approximate Bayesian model selection using the Bayesian Information Criterion and returns the log of the Effective Delay 50 and numerical integration area.
#'
#' @param dat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param A OPTIONAL: The maximum value of the commodity (e.g., of 100 USD, out of 1, etc.)
#' @param models vector of models to include in selection. Includes "noise", "hyperbolic", "exponential", "laibson", "greenmyerson", "rachlin", "ebertprelec", "crdi", and "genhyp".  All models are included by specifying "all" or leaving the parameter at default ("all").
#' @param idCol string content identifying participants (e.g., idCol = "ids")
#' @param detailed include additional output details in results (TRUE/FALSE, FALSE by default)
#' @param figures OPTIONAL: show figure of results. Includes "ed50" for ed50 plotting, "auc" for numerical integration with normal delays, and "logauc" for numerical integration with log base 10 scaled delays.
#' @param summarize OPTIONAL: Descriptive observations (TRUE/FALSE, FALSE by default)
#' @param lineSize OPTIONAL: Modify line sizes for figures (numeric, 1 by default)
#' @return A data frame of model parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of fitted model parameters
#' @examples
#' dat <- data.frame(X=c(1,30,180,540,1080,2160),
#' Y=c(1,0.9,0.8,0.7,0.6,0.4),
#' ids=c(1,1,1,1,1,1))
#'
#' dat2<- data.frame(X=c(1,30,180,540,1080,2160),
#'                   Y=c(1,0.9,0.7,0.6,0.5, 0.4),
#'                   ids=c(2,2,2,2,2,2))
#'
#' dat <- rbind(dat, dat2)
#'
#' dat$Y <- dat$Y * 100
#'
#' results <- discountingModelSelection(dat, idCol = "ids", A = 100)
#'
#' @export
discountingModelSelection <- function(dat, A = NULL, models = c("all"), idCol = "id", detailed = FALSE, figures = FALSE, summarize = FALSE, lineSize = 1) {

  mModels <- c("noise", "hyperbolic", "exponential", "laibson", "greenmyerson", "rachlin", "ebertprelec", "crdi", "genhyp")

  if (!"all" %in% models) {

    if (!"noise" %in% models) {
      models <- c("noise", models)
    }

    if (length(intersect(models, mModels)) < 2) {
      stop("At least one model must be specified to perform a comparison")
    }
  } else {
    models <- mModels

  }

  if (!is.null(A) && is.numeric(A)) {
    dat$Y <- dat$Y / A
  }

  if (any(dat$Y > 1)) {
    stop("Y values exceed the range of 0<=Y<=1, Did you forget to assign an A (max value) argument?")
  }

  if (!idCol %in% colnames(dat)) {
    stop("Id column not found, please check naming")

  } else {
    colnames(dat)[colnames(dat) == idCol] <- 'id'

  }

  lengthReturn <- length(unique(dat$id))

  cl <- makeCluster(2)
  registerDoParallel(cores=detectCores())

  finalResult <- NA

  #for (id in unique(dat$id)) {
  finalResult <- foreach(i = unique(dat$id), .combine=rbind, .packages = c("discountingtools")) %dopar% {
#    indData
#
#    resultFrame <- data.frame(id = i,
#                              Equation = equation)

    localDat <- subset(dat, id == i)

    result <- data.frame(id = i)

    if (nrow(localDat) < 3) {
      message(paste("Dropping Series: ", i, ", Fewer than 3 data points", sep = ""));
    } else {
      message(paste("Calculating series id: ", i, sep = ""))

      screenRes <- johnsonBickelScreen(localDat)

      result <- cbind(id = i,
                      JB.C1 = screenRes$C1,
                      JB.C2 = screenRes$C2,
                      discountingModelSelectionCall(localDat, A, models, detailed, figures, summarize, lineSize))
    }

    result

    #


#    if (id == unique(dat$id)[1]) {
#
#      finalResult <- cbind(id = id,
#                           JB.C1 = screenRes$C1,
#                           JB.C2 = screenRes$C2,
#                           result)
#
#    } else {
#
#      finalResult <- rbind(finalResult,
#                           cbind(id = id,
#                                 JB.C1 = screenRes$C1,
#                                 JB.C2 = screenRes$C2,
#                                 result))
#
    #}
  }

  stopCluster(cl)

  finalResult
}

#' Perform Johnson & Bickel Screen
#'
#' This function applies the Johnson & Bickel screening criteria to included data series. The result of this procedure is a TRUE/FALSE response to one of two screening criteria. These are included by default in all model selection calls.
#'
#' @param dat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param idCol string content identifying participants (e.g., idCol = "ids")
#' @return A data frame of model screenings
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#' dat <- data.frame(X=c(1,30,180,540,1080,2160),
#' Y=c(1,0.9,0.8,0.7,0.6,0.4),
#' ids=c(1,1,1,1,1,1))
#'
#' dat2<- data.frame(X=c(1,30,180,540,1080,2160),
#'                   Y=c(1,0.9,0.7,0.6,0.5,0.4),
#'                   ids=c(2,2,2,2,2,2))
#'
#' dat <- rbind(dat, dat2)
#'
#' dat$Y <- dat$Y * 100
#' johnsonBickelScreen(dat, idCol = "ids")
#'
#' @export
johnsonBickelScreen <- function(dat, idCol = "id") {

  if (!idCol %in% colnames(dat)) {
    stop("Id column not found, please check naming")
  } else {
    colnames(dat)[colnames(dat) == idCol] <- 'id'
  }

  lengthReturn <- length(unique(dat$id))

  returnFrame <- data.frame(id = rep(NA, lengthReturn),
                            C1 = rep(NA, lengthReturn),
                            C2 = rep(NA, lengthReturn))

  mIndex <- 1

  for (i in unique(dat$id)) {
    subsetData <- dat[dat$id == i,]

    criteriaOne <- TRUE
    criteriaTwo <- TRUE

    subsetData <- subsetData[order(subsetData$X), ]

    for (index in 2:length(subsetData$X)) {
      prev = subsetData[index-1, ]$Y
      curr = subsetData[index, ]$Y

      if ((curr - prev) > 0.2) {
        criteriaOne = FALSE
      }
    }

    prev <- subsetData[1, ]$Y
    curr <- subsetData[length(subsetData$X), ]$Y

    if ((prev - curr) < 0.1) {
      criteriaTwo = FALSE
    }

    returnFrame[mIndex, ]$id <- i
    returnFrame[mIndex, ]$C1 <- criteriaOne
    returnFrame[mIndex, ]$C2 <- criteriaTwo

    mIndex <- mIndex + 1

  }

  returnFrame
}

#' Exponential Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesExponential(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesExponential <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 0.01
  }

  startlnK  <- seq(-15,  15, increment)
  lengthLnK <- length(startlnK)
  lengthX <- length(localDat$X)

  # Params expand
  MlnK <- sort(rep(startlnK, lengthX))

  # Init SS vector
  sumSquares<- rep(NA,lengthLnK)

  # Observed Data
  MX <- rep(localDat$X, lengthLnK)
  MY <- rep(localDat$Y/A, lengthLnK)

  # Projections
  projection<- exp(-exp(MlnK)*MX)
  sqResidual<- (MY-projection)^2

  for(j in 1:lengthLnK){
    sumSquares[j]  <-sum(sqResidual[(j-1)*lengthX +1:lengthX])

  }

  # Sort starting estimates
  presort   <-data.frame(startlnK, sumSquares)
  sorted    <-presort[order(presort[ ,"sumSquares"]), ]
  ini.par   <-c(lnk=sorted$startlnK[1])

  ini.par
}

#' Hyperbolic Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesHyperbolic(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesHyperbolic <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 0.01
  }

  startlnK  <- seq(-15,  15, increment)
  lengthLnK <- length(startlnK)
  lengthX <- length(localDat$X)

  # Params expand
  MlnK <- sort(rep(startlnK, lengthX))

  # Init SS vector
  sumSquares <- rep(NA,lengthLnK)

  # Observed Data
  MX <- rep(localDat$X, lengthLnK)
  MY <- rep(localDat$Y/A, lengthLnK)

  # Projections
  projection <- (1+exp(MlnK)*MX)^(-1)
  sqResidual <- (MY-projection)^2

  for(j in 1:lengthLnK){
    sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

  }

  # Sort starting estimates
  presort <-data.frame(startlnK,
                       sumSquares)
  sorted  <-presort[order(presort[ ,"sumSquares"]), ]

  ini.par <-c(lnk=sorted$startlnK[1])

  ini.par
}

#' Laibson Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesLaibson(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesLaibson <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 3
  }

  # Starts
  startbeta   <- seq(0, 1, 0.1/increment)
  startdelta  <- seq(0, 1, 0.01/increment)

  lengthBeta  <- length(startbeta)
  lengthDelta <- length(startdelta)
  lengthX <- length(localDat$X)

  # Params expand
  startBeta   <- rep(sort(rep(startbeta, lengthX)),lengthDelta)
  startdelta  <- sort(rep(startdelta,lengthX*lengthBeta))

  # Init SS vector
  sumSquares  <- rep(NA,lengthBeta*lengthDelta)

  # Observed Data
  SY          <- rep(localDat$Y/A,lengthBeta*lengthDelta)

  # Projections
  projection  <- startBeta*startdelta^localDat$X
  sqResidual  <- (SY-projection)^2

  SSbeta      <- rep(startbeta,lengthDelta)
  SSdelta     <- sort(rep(startdelta,lengthBeta))

  for(j in 1:(lengthBeta*lengthDelta)){
    sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
  }

  # Sort starting estimates
  presort<-data.frame(SSbeta,
                      SSdelta,
                      sumSquares)

  sorted  <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(beta  = sorted$SSbeta[1],
               delta = sorted$SSdelta[1])

  ini.par
}

#' Green/Myerson Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesGreenMyerson(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesGreenMyerson <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 3
  }

  # Starts
  startlnK <- seq(-12, 12, 1/increment)
  starts <- seq(.01, 10, 0.01/increment)

  lengthLnK <- length(startlnK)
  lengthS <- length(starts)
  lengthX <- length(localDat$X)

  # Params expand
  SSlnK <- rep(startlnK,lengthS)
  SSs <- sort(rep(starts,lengthLnK))

  # Init SS vector
  sumSquares <- rep(NA,lengthS*lengthLnK)

  # Observed Data
  SY <- rep(localDat$Y/A,lengthS*lengthLnK)

  # Projections
  SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
  Ss <- sort(rep(starts,lengthX*lengthLnK))

  projection <- (1+exp(SlnK)*localDat$X)^(-Ss)
  sqResidual <- (SY-projection)^2

  for(j in 1:(lengthS*lengthLnK)){
    sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

  }

  # Sort starting estimates
  presort <- data.frame(SSlnK,SSs,sumSquares)
  sorted <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk=sorted$SSlnK[1],
               s=sorted$SSs[1])

  ini.par
}

#' Rachlin Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesRachlin(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesRachlin <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 3
  }

  # Starts
  startlnK <- seq(-12, 12, 1/increment)
  starts <- seq(.01, 10, 0.01/increment)

  lengthLnK<-length(startlnK)
  lengthS<-length(starts)
  lengthX <- length(localDat$X)

  # Params expand
  SSlnK<-rep(startlnK,lengthS)
  SSs<-sort(rep(starts,lengthLnK))

  # Init SS vector
  sumSquares<-rep(NA,lengthS*lengthLnK)

  # Observed Data
  SY<-rep(localDat$Y/A,lengthS*lengthLnK)

  # Projections
  SlnK<-rep(sort(rep(startlnK,lengthX)),lengthS)
  Ss<-sort(rep(starts,lengthX*lengthLnK))

  projection<-(1+exp(SlnK)*(localDat$X^Ss))^(-1)
  sqResidual<-(SY-projection)^2

  for(j in 1:(lengthS*lengthLnK)){
    sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

  }

  # Sort starting estimates
  presort <- data.frame(SSlnK,
                        SSs,
                        sumSquares)
  sorted <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk=sorted$SSlnK[1],
               s=sorted$SSs[1])

  ini.par
}

#' Ebert/Prelec Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesEbertPrelec(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesEbertPrelec <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 3
  }

  # Starts
  startlnK <- seq(-12, 12, 1/increment)
  starts <- seq(.01, 10, 0.01/increment)

  lengthLnK <- length(startlnK)
  lengthS <- length(starts)
  lengthX <- length(localDat$X)

  # Params expand
  SSlnK <- rep(startlnK,lengthS)
  SSs <- sort(rep(starts,lengthLnK))

  # Init SS vector
  sumSquares <- rep(NA,lengthS*lengthLnK)

  # Observed Data
  SY <- rep(localDat$Y/A,lengthS*lengthLnK)

  # Projections
  SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
  Ss <- sort(rep(starts,lengthX*lengthLnK))

  projection <- exp(-(exp(SlnK)*localDat$X)^Ss)
  sqResidual <- (SY-projection)^2

  for(j in 1:(lengthS*lengthLnK)){
    sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])

  }

  # Sort starting estimates
  presort <- data.frame(SSlnK,
                        SSs,
                        sumSquares)
  sorted <- presort[order(presort[,"sumSquares"]),]
  ini.par <- c(lnk=sorted$SSlnK[1],
               s=sorted$SSs[1])

  ini.par
}

#' Bleichrodt et al. Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesBleichrodt(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesBleichrodt <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 2
  }

  startlnK <- seq(-12, 12, 1/increment)
  starts <- seq(.01, 1, 0.11/increment)
  startBeta <- seq(0.1, 1, 0.11/increment)

  lengthX <- length(localDat$X)

  # new Pre sort
  presort <- expand.grid(startlnK = startlnK, starts = starts, startBeta = startBeta)
  presort$sumSquares <- NA

  # clean, merge, or move
  getSS <- function(presort, index, Y, X) {
    projections <- presort[index,]$startBeta*exp(-(exp(presort[index,]$startlnK)*X^presort[index,]$starts))
    sqResidual <- (Y - projections)^2
    sum(sqResidual)
  }

  for (j in 1:nrow(presort)) {
    presort[j, ]$sumSquares <-getSS(presort, j, localDat$Y/A, localDat$X)
  }

  presort <- presort[order(presort$sumSquares),]

  ini.par <- c(beta = presort[1,]$startBeta, lnk = presort[1,]$startlnK, s = presort[1,]$starts)

  ini.par
}

#' Rodriguez/Logue Starting Values estimator
#'
#' This function applies a grid search for starting parameters
#'
#' @param A Maximum value of commodity
#' @param localDat data frame with X column and Y column (0 <= Y <= 1, NOTE capitalized letters)
#' @param overkill optional higher resolution search
#' @return A named vector of starting parameters
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
#' @return data frame of Screening Criteria
#' @examples
#'dat <- data.frame(X=c(1,30,180,540,1080,2160),
#'                  Y=c(1,0.9,0.8,0.7,0.6,0.4),
#'                  ids=c(1,1,1,1,1,1))
#'dat$Y <- dat$Y * 100
#'startValuesRodriguezLogue(A = 100, dat, overkill = TRUE)
#'
#' @export
startValuesRodriguezLogue <- function(A = 1, localDat, overkill = FALSE) {

  increment <- 1

  if (overkill == TRUE) {
    increment <- 3
  }

  # Starts
  startlnK <- seq(-12, 12, 1/increment)
  startBeta <- seq(-12, 12, 1/increment)

  # new Pre sort
  presort <- expand.grid(startlnK = startlnK, startBeta = startBeta)
  presort$sumSquares <- NA

  # clean, merge, or move
  getSS <- function(presort, index, Y, X) {
    projections <- (1 + X * exp(presort[index,]$startlnK))^(-exp(presort[index,]$startBeta) / exp(presort[index,]$startlnK))
    sqResidual <- (Y - projections)^2
    sum(sqResidual)
  }

  for (j in 1:nrow(presort)) {
    presort[j, ]$sumSquares <-getSS(presort, j, localDat$Y/A, localDat$X)
  }

  presort <- presort[order(presort$sumSquares),]

  ini.par <- c(lnk = presort[1,]$startlnK, beta = presort[1,]$startBeta)

  ini.par
}
