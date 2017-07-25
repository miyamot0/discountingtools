
#' Perform Discounting Model Selection
#'
#' This function takes a data frame of temporal discounting values (X, Y)
#' and performs approximate Bayesian model selection using the Bayesian
#' Information Criterion and returns the log of the Effective Delay
#' 50 and numerical integration area.
#'
#' Models: Exponential, Hyperbolic, BetaDelta, GreenMyerson, & Rachlin, Ebert & Prelec's Constant Sensitivity
#'
#' @param dat data frame with X column and Y column (0 <= Y <= 1)
#' @param models vector of models to include in selection
#' @param figures OPTIONAL: show figure of results
#' @param summarize OPTIONAL: Descriptive observations
#' @param lineSize OPTIONAL: Modify line sizes for figures
#' @return A data frame of model parameters
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @importFrom stats lm nls
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @return data frame of fitted model parameters
#' @examples
#' discountingModelSelection(data.frame(X=c(1,30,180,540,1080,2160),
#' Y=c(1.0,0.9,0.8,0.7,0.6,0.4)),
#' models=c("noise", "exponential", "hyperbolic", "bd", "mg", "rachlin", "cs",
#' Figures = FALSE))
#' @export
discountingModelSelection <- function(dat, models = c("noise"), figures = FALSE, summarize = FALSE, lineSize = 1) {

  if (!"noise" %in% models) {
    models <- c("noise", models)
  }

  mModels <- c("noise", "hyperbolic", "exponential", "bd", "mg", "rachlin", "cs")

  if (length(intersect(models, mModels)) < 2) {
    stop("At least one model must be specified to perform a comparison")
  }

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

  if(!is.character(try(lm(Y ~ 1, dat), silent=TRUE))) {
    modelFitNoise <- stats::lm(Y ~ 1, dat)

    tempList <- list(noise.mean = modelFitNoise$coefficients[["(Intercept)"]],
                     noise.RMSE = summary(modelFitNoise)[["sigma"]],
                     noise.BIC  = stats::BIC(modelFitNoise),
                     noise.AIC  = stats::AIC(modelFitNoise))

    returnList <- c(returnList, tempList)

    bicList <- c(bicList, list(noise.BIC = tempList$noise.BIC))
  }

  if ('exponential' %in% models) {
    startlnK  <- seq(-12,  12, 1)
    lengthLnK <- length(startlnK)

    sumSquares<- rep(NA,lengthLnK)

    MX <- rep(dat$X, lengthLnK)
    MY <- rep(dat$Y, lengthLnK)

    projection<- exp(-exp(lengthLnK)*MX)
    sqResidual<- (MY-projection)^2

    for(j in 1:lengthLnK){
      sumSquares[j]  <-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort   <-data.frame(startlnK, sumSquares)
    sorted    <-presort[order(presort[ ,"sumSquares"]), ]
    ini.par   <-c(lnk=sorted$startlnK[1])

    if (!is.character(try(nls.lm(par = ini.par,
                                 fn = residualFunction,
                                 jac = jacobianMatrix,
                                 valueFunction = exponentialDiscountFunc,
                                 jacobianFunction = exponentialDiscountGradient,
                                 x = dat$X,
                                 value = dat$Y,
                                 control = nls.lm.control(maxiter = 1000)), silent=TRUE))) {
      modelFitExponential<-nls.lm(par = ini.par,
                                  fn = residualFunction,
                                  jac = jacobianMatrix,
                                  valueFunction = exponentialDiscountFunc,
                                  jacobianFunction = exponentialDiscountGradient,
                                  x = dat$X,
                                  value = dat$Y,
                                  control = nls.lm.control(maxiter = 1000))

      tempList <- list(exp.lnk = modelFitExponential$par[["lnk"]],
                       exp.RMSE = summary(modelFitExponential)[["sigma"]],
                       exp.BIC  = stats::BIC(logLik.nls.lm(modelFitExponential)),
                       exp.AIC  = stats::AIC(logLik.nls.lm(modelFitExponential)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(exp.BIC = tempList$exp.BIC))
    }
  }

  if ('hyperbolic' %in% models) {
    startlnK <- seq(-12,  12, 1)
    lengthLnK <- length(startlnK)

    sumSquares <- rep(NA,lengthLnK)

    MlnK <- sort(rep(startlnK, lengthX))

    MX <- rep(dat$X, lengthLnK)
    MY <- rep(dat$Y, lengthLnK)

    projection <- (1+exp(MlnK)*MX)^(-1)
    sqResidual <- (MY-projection)^2

    for(j in 1:lengthLnK){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort <-data.frame(startlnK, sumSquares)
    sorted  <-presort[order(presort[ ,"sumSquares"]), ]
    ini.par <-c(lnk=sorted$startlnK[1])

    if (!is.character(try(nls.lm(par = ini.par,
                                 fn = residualFunction,
                                 jac = jacobianMatrix,
                                 valueFunction = hyperbolicDiscountFunc,
                                 jacobianFunction = hyperbolicDiscountGradient,
                                 x = dat$X,
                                 value = dat$Y,
                                 control = nls.lm.control(maxiter = 1000)), silent=FALSE))) {
      modelFitHyperbolic <- nls.lm(par = ini.par,
                                   fn = residualFunction,
                                   jac = jacobianMatrix,
                                   valueFunction = hyperbolicDiscountFunc,
                                   jacobianFunction = hyperbolicDiscountGradient,
                                   x = dat$X,
                                   value = dat$Y,
                                   control = nls.lm.control(maxiter = 1000))

      tempList <- list(Mazur.lnk  = modelFitHyperbolic$par[["lnk"]],
                       Mazur.RMSE = summary(modelFitHyperbolic)[["sigma"]],
                       Mazur.BIC  = stats::BIC(logLik.nls.lm(modelFitHyperbolic)),
                       Mazur.AIC  = stats::AIC(logLik.nls.lm(modelFitHyperbolic)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(Mazur.BIC = tempList$Mazur.BIC))
    }
  }

  if ('bd' %in% models) {
    startbeta   <- seq(0, 1, 0.1)
    startdelta  <- seq(0, 1, 0.01)

    lengthBeta  <- length(startbeta)
    lengthDelta <- length(startdelta)

    startBeta   <- rep(sort(rep(startbeta, lengthX)),lengthDelta)
    startdelta  <- sort(rep(startdelta,lengthX*lengthBeta))

    SY          <- rep(dat$Y,lengthBeta*lengthDelta)

    projection  <- startBeta*startdelta^dat$X
    sqResidual  <- (SY-projection)^2

    SSbeta      <- rep(startbeta,lengthDelta)
    SSdelta     <- sort(rep(startdelta,lengthBeta))
    sumSquares  <- rep(NA,lengthBeta*lengthDelta)

    for(j in 1:(lengthBeta*lengthDelta)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort<-data.frame(SSbeta,
                        SSdelta,
                        sumSquares)

    sorted  <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(beta  = sorted$SSbeta[1],
                 delta = sorted$SSdelta[1])

    if (!is.character(try(nls.lm(par = ini.par,
                                 fn = residualFunction,
                                 jac = jacobianMatrix,
                                 valueFunction = betaDeltaDiscountFunc,
                                 jacobianFunction = betaDeltaDiscountGradient,
                                 x = dat$X,
                                 value = dat$Y,
                                 upper = c(beta = 1, delta = 1),
                                 lower = c(beta = 0, delta = 0),
                                 control = nls.lm.control(maxiter = 1000)), silent=FALSE))) {
      modelFitBetaDelta<-nls.lm(par = ini.par,
                                fn = residualFunction,
                                jac = jacobianMatrix,
                                valueFunction = betaDeltaDiscountFunc,
                                jacobianFunction = betaDeltaDiscountGradient,
                                x = dat$X,
                                value = dat$Y,
                                upper = c(beta = 1, delta = 1),
                                lower = c(beta = 0, delta = 0),
                                control = nls.lm.control(maxiter = 1000))

      tempList <- list(BD.beta  = modelFitBetaDelta$par[["beta"]],
                       BD.delta  = modelFitBetaDelta$par[["delta"]],
                       BD.RMSE = summary(modelFitBetaDelta)[["sigma"]],
                       BD.BIC  = stats::BIC(logLik.nls.lm(modelFitBetaDelta)),
                       BD.AIC  = stats::AIC(logLik.nls.lm(modelFitBetaDelta)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(BD.BIC = tempList$BD.BIC))
    }
  }

  if ('mg' %in% models) {
    startlnK <- seq(-12, 12, 1)
    starts <- seq(.01, 10, 0.01)

    lengthLnK <- length(startlnK)
    lengthS <- length(starts)

    SSlnK <- rep(startlnK,lengthS)
    SSs <- sort(rep(starts,lengthLnK))

    sumSquares <- rep(NA,lengthS*lengthLnK)

    SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss <- sort(rep(starts,lengthX*lengthLnK))

    SY <- rep(dat$Y,lengthS*lengthLnK)

    projection <- (1+exp(SlnK)*dat$X)^(-Ss)
    sqResidual <- (SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort <- data.frame(SSlnK,SSs,sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    if (!is.character(try(nls.lm(par = ini.par,
                                 fn = residualFunction,
                                 jac = jacobianMatrix,
                                 valueFunction = myersonHyperboloidDiscountFunc,
                                 jacobianFunction = myersonHyperboloidDiscountGradient,
                                 x = dat$X,
                                 value = dat$Y,
                                 control = nls.lm.control(maxiter = 1000)), silent=FALSE))) {
      modelFitMyerson<-nls.lm(par = ini.par,
                              fn = residualFunction,
                              jac = jacobianMatrix,
                              valueFunction = myersonHyperboloidDiscountFunc,
                              jacobianFunction = myersonHyperboloidDiscountGradient,
                              x = dat$X,
                              value = dat$Y,
                              control = nls.lm.control(maxiter = 1000))

      tempList <- list(MG.lnk  = modelFitMyerson$par[["lnk"]],
                       MG.s  = modelFitMyerson$par[["s"]],
                       MG.RMSE = summary(modelFitMyerson)[["sigma"]],
                       MG.BIC  = stats::BIC(logLik.nls.lm(modelFitMyerson)),
                       MG.AIC  = stats::AIC(logLik.nls.lm(modelFitMyerson)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(MG.BIC = tempList$MG.BIC))
    }
  }

  if ('rachlin' %in% models) {
    startlnK<-seq(-12,12,1)
    starts<-seq(.01,10,.01)

    lengthLnK<-length(startlnK)
    lengthS<-length(starts)

    SSlnK<-rep(startlnK,lengthS)
    SSs<-sort(rep(starts,lengthLnK))

    sumSquares<-rep(NA,lengthS*lengthLnK)

    SlnK<-rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss<-sort(rep(starts,lengthX*lengthLnK))

    SY<-rep(dat$Y,lengthS*lengthLnK)

    projection<-(1+exp(SlnK)*(dat$X^Ss))^(-1)
    sqResidual<-(SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort <- data.frame(SSlnK,SSs,sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    if (!is.character(try(nls.lm(par = ini.par,
                                 fn = residualFunction,
                                 jac = jacobianMatrix,
                                 valueFunction = rachlinHyperboloidDiscountFunc,
                                 jacobianFunction = rachlinHyperboloidDiscountGradient,
                                 x = dat$X,
                                 value = dat$Y,
                                 control = nls.lm.control(maxiter = 1000)), silent=FALSE))) {
      modelFitRachlin<-nls.lm(par = ini.par,
                              fn = residualFunction,
                              jac = jacobianMatrix,
                              valueFunction = rachlinHyperboloidDiscountFunc,
                              jacobianFunction = rachlinHyperboloidDiscountGradient,
                              x = dat$X,
                              value = dat$Y,
                              control = nls.lm.control(maxiter = 1000))

      tempList <- list(Rachlin.lnk  = modelFitRachlin$par[["lnk"]],
                       Rachlin.s  = modelFitRachlin$par[["s"]],
                       Rachlin.RMSE = summary(modelFitRachlin)[["sigma"]],
                       Rachlin.BIC  = stats::BIC(logLik.nls.lm(modelFitRachlin)),
                       Rachlin.AIC  = stats::AIC(logLik.nls.lm(modelFitRachlin)))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(Rachlin.BIC = tempList$Rachlin.BIC))

    }
  }

  if ('cs' %in% models) {
    startlnK <- seq(-12, 12, 0.1)
    starts <- seq(.01, 1, 0.01)

    lengthLnK <- length(startlnK)
    lengthS <- length(starts)

    SSlnK <- rep(startlnK,lengthS)
    SSs <- sort(rep(starts,lengthLnK))

    sumSquares <- rep(NA,lengthS*lengthLnK)

    SlnK <- rep(sort(rep(startlnK,lengthX)),lengthS)
    Ss <- sort(rep(starts,lengthX*lengthLnK))

    SY <- rep(dat$Y,lengthS*lengthLnK)

    projection <- exp(-(exp(SlnK)*dat$X)^Ss)#(1+exp(SlnK)*dat$X)^(-Ss)
    sqResidual <- (SY-projection)^2

    for(j in 1:(lengthS*lengthLnK)){
      sumSquares[j]<-sum(sqResidual[(j-1)*lengthX +1:lengthX])
    }

    presort <- data.frame(SSlnK,SSs,sumSquares)
    sorted <- presort[order(presort[,"sumSquares"]),]
    ini.par <- c(lnk=sorted$SSlnK[1],
                 s=sorted$SSs[1])

    if (!is.character(try(nls(Y ~ exp(-(exp(lnk)*X)^s),
                              start = ini.par,
                              data = dat), silent=FALSE))) {
      modelFitCS<-nls(Y ~ exp(-(exp(lnk)*X)^s),
                      start = ini.par,
                      data = dat)

      tempList <- list(CS.lnk  = stats::coef(modelFitCS)[["lnk"]],
                       CS.s  = stats::coef(modelFitCS)[["s"]],
                       CS.RMSE = summary(modelFitCS)[["sigma"]],
                       CS.BIC  = stats::BIC(modelFitCS),
                       CS.AIC  = stats::AIC(modelFitCS))

      returnList <- c(returnList, tempList)

      bicList <- c(bicList, list(CS.BIC = tempList$CS.BIC))
    }
  }

  ### bfs
  bfList <- list()
  bfSum <- 0.0

  for (i in 1:length(bicList)) {
    bfName = gsub("BIC", "BF", names(bicList)[i])
    bfList[[bfName]] = exp(-.5*(bicList[[names(bicList)[i]]]-bicList[["noise.BIC"]]))

    bfSum <- bfSum + bfList[[bfName]]
  }

  returnList <- c(returnList, bfList)

  ### probs
  probList <- list()

  for (i in 1:length(bfList)) {
    probName = gsub("BF", "prob", names(bfList)[i])
    probList[[probName]] = bfList[[names(bfList)[i]]]/bfSum
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
      message("")
    }
  }

  ### Plotting here
  if (figures == "ed50") {
    message("Plotting ED50 Data...")
    displayED50Figure(dat, returnList, lineSize)
  } else if (figures == "auc") {
    message("Plotting Model Area Data...")
    displayAUCFigure(dat, returnList, lineSize)
  } else if (figures == "logauc") {
    displayLogAUCFigure(dat, returnList, lineSize)
  }

  returnFrame <- as.data.frame(returnList)
  returnFrame
}

#' Perform Johnson & Bickel Screen
#'
#' This function applies the Johnson & Bickel screening criteria to included data series
#'
#' @param dat data frame with X column and Y column (0 <= Y <= 1)
#' @param idCol id column
#' @return A data frame of model screenings
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return data frame of Screening Criteria
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
