#' Scoring for the log ED50
#'
#' Optionally, models without a straightforward exact solution
#' will be scored numerically using a bisection search
#'
#' Only Ebert & Prelec, 2007 requires this bisection search
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return natural logarithm of the Effective Delay 50%
getED50 <- function(dat, results) {
  returnValue <- NaN

  if (results[["probable.model"]] == "Hyperbolic") {
    returnValue <- log(1/(exp(results[["Hyperbolic.lnk"]])))
  } else if (results[["probable.model"]] == "Exponential") {
    returnValue <- log(log(2)/exp(results[["Exponential.lnk"]]))
  } else if (results[["probable.model"]] == "Laibson") {
    returnValue <- log(log( (1/(2*results[["Laibson.beta"]])),base=results[["Laibson.delta"]]))
  } else if (results[["probable.model"]] == "GreenMyerson") {
    returnValue <- log( (2^(1/results[["GreenMyerson.s"]])-1)/exp(results[["GreenMyerson.lnk"]]))
  } else if (results[["probable.model"]] == "Rachlin") {
    returnValue <- log( (1/(exp(results[["Rachlin.lnk"]])))^(1/results[["Rachlin.s"]]))
  } else if (results[["probable.model"]] == "EbertPrelec") {
    returnValue <- getED50ep(dat, results)
  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    returnValue <- getED50crdi(dat, results)
  } else if (results[["probable.model"]] == "RodriguezLogue") {
    returnValue <- getED50genhyp(dat, results)
  }

  returnValue
}

#' Numerically solve for ED50 value for Ebert & Prelec
#'
#' This method solves for ED50 for Ebert & Prelec using a point bisection procedure.
#' This procedure will continue for n (20 currently) by default until a value of 50%
#' is observed in the midpoint of two more moving delays.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return effective delay (value) for Ebert & Prelec ep
getED50ep <- function(dat, results) {
  lowDelay <- 0
  highDelay <- max(dat$X)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst <- integrandEbertPrelec(lowDelay, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])
    midEst <- integrandEbertPrelec((lowDelay+highDelay)/2, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])
    highEst <- integrandEbertPrelec(highDelay, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])

    if (lowEst > 0.5 && midEst > 0.5) {
      # Above 50% mark range
      lowDelay <- (lowDelay+highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      # Below 50% mark range
      lowDelay <- lowDelay
      highDelay <- (lowDelay+highDelay)/2

    }
  }

  # NECRO
  #for (i in seq(1, 20)) {
  #  lowEst <- integrandEbertPrelec(lowDelay, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])
  #  midEst <- integrandEbertPrelec((lowDelay+highDelay)/2, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])
  #  highEst <- integrandEbertPrelec(highDelay, results[["EbertPrelec.lnk"]], results[["EbertPrelec.s"]])

  #  if (lowEst > 0.5 && midEst > 0.5) {
  # Above 50% mark range
  #    lowDelay <- (lowDelay+highDelay)/2
  #    highDelay <- highDelay

  #  } else if (highEst < 0.5 && midEst < 0.5) {
  # Below 50% mark range
  #    lowDelay <- lowDelay
  #    highDelay <- (lowDelay+highDelay)/2

  #  }
  #}

  returnValue <- log((lowDelay+highDelay)/2)

  returnValue
}

#' Numerically solve for ED50 value for Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI)
#'
#' This method solves for ED50 for Bleichrodt CRDI using a point bisection procedure.
#' This procedure will continue for n (20 currently) by default until a value of 50%
#' is observed in the midpoint of two more moving delays.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return effective delay (value) for Ebert & Prelec ep
getED50crdi <- function(dat, results) {
  lowDelay <- 0
  highDelay <- max(dat$X)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst <- integrandBleichrodtCRDI(lowDelay, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])
    midEst <- integrandBleichrodtCRDI((lowDelay+highDelay)/2, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])
    highEst <- integrandBleichrodtCRDI(highDelay, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])

    if (lowEst > 0.5 && midEst > 0.5) {
      # Above 50% mark range
      lowDelay <- (lowDelay+highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      # Below 50% mark range
      lowDelay <- lowDelay
      highDelay <- (lowDelay+highDelay)/2

    }
  }

  # NECRO
  #for (i in seq(1, 20)) {
  #  lowEst <- integrandBleichrodtCRDI(lowDelay, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])
  #  midEst <- integrandBleichrodtCRDI((lowDelay+highDelay)/2, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])
  #  highEst <- integrandBleichrodtCRDI(highDelay, results[["BleichrodtCRDI.lnk"]], results[["BleichrodtCRDI.s"]], results[["BleichrodtCRDI.beta"]])

  #  if (lowEst > 0.5 && midEst > 0.5) {
  # Above 50% mark range
  #    lowDelay <- (lowDelay+highDelay)/2
  #    highDelay <- highDelay

  #  } else if (highEst < 0.5 && midEst < 0.5) {
  # Below 50% mark range
  #    lowDelay <- lowDelay
  #    highDelay <- (lowDelay+highDelay)/2

  #  }
  #}

  returnValue <- log((lowDelay+highDelay)/2)

  returnValue
}

#' Numerically solve for ED50 value for Rodriguez & Logue Generalized Hyperbolioid
#'
#' This method solves for ED50 for Rodriguez & Logue using a point bisection procedure.
#' This procedure will continue for n (20 currently) by default until a value of 50%
#' is observed in the midpoint of two more moving delays.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return effective delay (value) for Ebert & Prelec ep
getED50genhyp <- function(dat, results) {
  lowDelay <- 0
  highDelay <- max(dat$X)*10

  while ((highDelay - lowDelay) > 0.001) {
    lowEst <- integrandRodriguezLogue(lowDelay, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])
    midEst <- integrandRodriguezLogue((lowDelay+highDelay)/2, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])
    highEst<- integrandRodriguezLogue(highDelay, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])

    if (lowEst > 0.5 && midEst > 0.5) {
      # Above 50% mark range
      lowDelay <- (lowDelay+highDelay)/2
      highDelay <- highDelay

    } else if (highEst < 0.5 && midEst < 0.5) {
      # Below 50% mark range
      lowDelay <- lowDelay
      highDelay <- (lowDelay+highDelay)/2

    }
  }

  # NECRO
  #for (i in seq(1, 20)) {
  #  lowEst <- integrandRodriguezLogue(lowDelay, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])
  #  midEst <- integrandRodriguezLogue((lowDelay+highDelay)/2, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])
  #  highEst<- integrandRodriguezLogue(highDelay, results[["RodriguezLogue.lnk"]], results[["RodriguezLogue.beta"]])

  #  if (lowEst > 0.5 && midEst > 0.5) {
  # Above 50% mark range
  #    lowDelay <- (lowDelay+highDelay)/2
  #    highDelay <- highDelay

  #  } else if (highEst < 0.5 && midEst < 0.5) {
  # Below 50% mark range
  #    lowDelay <- lowDelay
  #    highDelay <- (lowDelay+highDelay)/2

  #  }
  #}

  returnValue <- log((lowDelay+highDelay)/2)

  returnValue
}
