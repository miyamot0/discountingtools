#' Generalized residual call
#'
#' General, shared method for coordinating with fitting functions
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return residual value of referenced function
residualFunction <- function(params, x, value, valueFunction, jacobianFunction)
{
  value - do.call("valueFunction", c(list(x = x), as.list(params)))
}

#' Generalized Jacobian call
#'
#' General shared method for constructing Jacobian
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return difference value for jacobian
jacobianMatrix <- function(params, x, value, valueFunction, jacobianFunction)
{
  -do.call("jacobianFunction", c(list(x = x), as.list(params)))
}

#' Exponential Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
exponentialDiscountFunc <- function(x, lnk)
{
  func <- exp(-exp(lnk)*x)
  eval(func)
}

#' Exponential Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
exponentialDiscountGradient <- function(x, lnk)
{
  func <- expression(exp(-exp(lnk)*x))
  c(eval(stats::D(func, "lnk")))
}

#' Hyperbolic Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
hyperbolicDiscountFunc <- function(x, lnk)
{
  func <- (1+exp(lnk)*x)^(-1)
  eval(func)
}

#' Hyperbolic Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
hyperbolicDiscountGradient <- function(x, lnk)
{
  func <- expression((1+exp(lnk)*x)^(-1))
  c(eval(stats::D(func, "lnk")))
}

#' Beta Delta Value Function
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
betaDeltaDiscountFunc <- function(x, beta, delta)
{
  func <- beta*delta^x
  eval(func)
}

#' Beta Delta Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
betaDeltaDiscountGradient <- function(x, beta, delta)
{
  func <- expression(beta*delta^x)
  c(eval(stats::D(func, "delta")),
    eval(stats::D(func, "beta")))
}

#' Green & Myerson Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
myersonHyperboloidDiscountFunc <- function(x, lnk, s)
{
  func <- (1+exp(lnk)*x)^(-s)
  eval(func)
}

#' Green & Myerson Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
myersonHyperboloidDiscountGradient <- function(x, lnk, s)
{
  func <- expression((1+exp(lnk)*x)^(-s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Rachlin Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
rachlinHyperboloidDiscountFunc <- function(x, lnk, s)
{
  func <- (1+exp(lnk)*(x^s))^(-1)
  eval(func)
}

#' Rachlin Gradient Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
rachlinHyperboloidDiscountGradient <- function(x, lnk, s)
{
  func <- expression((1+exp(lnk)*x)^(-s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Ebert & Prelec Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
ebertPrelecDiscountFunc <- function(x, lnk, s)
{
  func <- exp(-(exp(lnk)*x)^s)
  eval(func)
}

#' Ebert & Prelec Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
ebertPrelecDiscountGradient <- function(x, lnk, s)
{
  func <- expression(exp(-(exp(lnk)*x)^s))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")))
}

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
BleichrodtCRDIDiscountFunc <- function(x, lnk, s, beta)
{
  func <- beta * exp(-exp(lnk)*x^s)
  eval(func)
}

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
BleichrodtCRDIDiscountGradient <- function(x, lnk, s, beta)
{
  func <- beta * exp(-exp(lnk)*x^s)
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "s")),
    eval(stats::deriv(func, "beta")))
}

#' Rodriguez & Logue Value Function
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
RodriguezLogueDiscountFunc <- function(x, lnk, beta)
{
  func <- (1 + x * exp(lnk))^(-exp(beta) / exp(lnk))
  eval(func)
}

#' Rodriguez & Logue Helper for Nonlinear Fitting
#'
#' @param x observation at point n (X)
#' @param lnk fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return projected, subjective value
RodriguezLogueDiscountGradient <- function(x, lnk, beta)
{
  func <- expression((1 + x * exp(lnk))^(-exp(beta) / exp(lnk)))
  c(eval(stats::deriv(func, "lnk")),
    eval(stats::deriv(func, "beta")))
}

#' minpack.lm logLik hack
#'
#' @param fit nls.lm object
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

#' Exponential Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandExp <- function(x, lnK) { exp(-exp(lnK)*x) }

#' Exponential Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandExpLog <- function(x, lnK) { exp(-exp(lnK)*(10^x)) }

#' Hyperbolic Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandHyp <- function(x, lnK) { (1+exp(lnK)*x)^(-1) }

#' Hyperbolic Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandHypLog <- function(x, lnK) { (1+exp(lnK)*(10^x))^(-1) }

#' Beta Delta Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBetaDelta <- function(x, beta, delta) { beta*delta^x }

#' Beta Delta Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBetaDeltaLog <- function(x, beta, delta) { beta*delta^(10^x) }

#' Green & Myerson Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandMyerson <- function(x, lnK, s) { (1+exp(lnK)*x)^(-s) }

#' Green & Myerson Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandMyersonLog <- function(x, lnK, s) { (1+exp(lnK)*(10^x))^(-s) }

#' Rachlin Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRachlin <- function(x, lnK, s) { (1+exp(lnK)*(x^s))^(-1) }

#' Rachlin Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRachlinLog <- function(x, lnK, s) { (1+exp(lnK)*((10^x)^s))^(-1) }

#' Ebert & Prelec's Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandEbertPrelec <- function(x, lnK, s) {  exp(-(exp(lnK)*x)^s) }

#' Ebert & Prelec's ep Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandEbertPrelecLog <- function(x, lnK, s) {  exp(-(exp(lnK)*(10^x))^s) }

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBleichrodtCRDI <- function(x, lnK, s, beta) {  beta * exp(-exp(lnK)*x^s) }

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBleichrodtCRDILog <- function(x, lnK, s, beta) {  beta * exp(-exp(lnK)*(10^x)^s) }

#' Rodriguez & Logue Integrand helper
#'
#' This integrand helper is a projection of the integrand
#' with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRodriguezLogue <- function(x, lnK, beta) { (1 + x * exp(lnK))^(-beta / exp(lnK)) }

#' Rodriguez & Logue Integrand helper
#'
#' This integrand helper is a projection of the integrand (log10)
#' with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRodriguezLogueLog <- function(x, lnK, beta) { (1 + (10^x) * exp(lnK))^(-beta / exp(lnK)) }

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

#' Scoring for the most probable model area
#'
#' In this set of methods, the area beneath the fitted model is
#' calculated and divided by the maximum area using numerical integration
#' methods.  All delays are calculated in the normal scale.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return area beneath the fitted model
getModelAUC <- function(dat, results) {

  maximumArea <- max(dat$X) - min(dat$X)

  returnValue <- NaN

  if (results[["probable.model"]] == "Hyperbolic") {
    returnValue <- stats::integrate(integrandHyp,
                             lower = min(dat$X),
                             upper = max(dat$X),
                             lnK = results[["Hyperbolic.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Exponential") {
    returnValue <- stats::integrate(integrandExp,
                             lower = min(dat$X),
                             upper = max(dat$X),
                             lnK = results[["Exponential.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Laibson") {
    returnValue <- stats::integrate(integrandBetaDelta,
                             lower = min(dat$X),
                             upper = max(dat$X),
                             beta = results[["Laibson.beta"]],
                             delta = results[["Laibson.delta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "GreenMyerson") {
    returnValue <- stats::integrate(integrandMyerson,
                             lower = min(dat$X),
                             upper = max(dat$X),
                             lnK = results[["GreenMyerson.lnk"]],
                             s = results[["GreenMyerson.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Rachlin") {
    returnValue <- stats::integrate(integrandRachlin,
                             lower = min(dat$X),
                             upper = max(dat$X),
                             lnK = results[["Rachlin.lnk"]],
                             s = results[["Rachlin.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "EbertPrelec") {
    returnValue <- stats::integrate(integrandEbertPrelec,
                                    lower = min(dat$X),
                                    upper = max(dat$X),
                                    lnK = results[["EbertPrelec.lnk"]],
                                    s = results[["EbertPrelec.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    returnValue <- stats::integrate(integrandBleichrodtCRDI,
                                    lower = min(dat$X),
                                    upper = max(dat$X),
                                    lnK = results[["BleichrodtCRDI.lnk"]],
                                    s = results[["BleichrodtCRDI.s"]],
                                    beta = results[["BleichrodtCRDI.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "RodriguezLogue") {
    returnValue <- stats::integrate(integrandRodriguezLogue,
                                    lower = min(dat$X),
                                    upper = max(dat$X),
                                    lnK = results[["RodriguezLogue.lnk"]],
                                    beta = results[["RodriguezLogue.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Noise") {
    returnValue <- results[["Noise.mean"]]
  }

  returnValue
}

#' Scoring for the most probable model area, in log10 space
#'
#' In this set of methods, the area beneath the fitted model is
#' calculated and divided by the maximum area using numerical integration
#' methods.  All delays are calculated in the log base 10 scale.
#'
#' @param dat observed data
#' @param results Results of analyses for data series
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return area beneath the fitted model, in log10 space
getModelAUCLog10Scaled <- function(dat, results) {

  maximumArea <- log10(max(dat$X)) - log10(min(dat$X))

  returnValue <- NaN

  if (results[["probable.model"]] == "Hyperbolic") {
    returnValue <- stats::integrate(integrandHypLog,
                             lower = log10(min(dat$X)),
                             upper = log10(max(dat$X)),
                             lnK = results[["Hyperbolic.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Exponential") {
    returnValue <- stats::integrate(integrandExpLog,
                             lower = log10(min(dat$X)),
                             upper = log10(max(dat$X)),
                             lnK = results[["Exponential.lnk"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Laibson") {
    returnValue <- stats::integrate(integrandBetaDeltaLog,
                             lower = log10(min(dat$X)),
                             upper = log10(max(dat$X)),
                             beta = results[["Laibson.beta"]],
                             delta = results[["Laibson.delta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "GreenMyerson") {
    returnValue <- stats::integrate(integrandMyersonLog,
                             lower = log10(min(dat$X)),
                             upper = log10(max(dat$X)),
                             lnK = results[["GreenMyerson.lnk"]],
                             s = results[["GreenMyerson.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Rachlin") {
    returnValue <- stats::integrate(integrandRachlinLog,
                             lower = log10(min(dat$X)),
                             upper = log10(max(dat$X)),
                             lnK = results[["Rachlin.lnk"]],
                             s = results[["Rachlin.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "EbertPrelec") {
    returnValue <- stats::integrate(integrandEbertPrelecLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["EbertPrelec.lnk"]],
                                    s = results[["EbertPrelec.s"]])$value/maximumArea

  } else if (results[["probable.model"]] == "BleichrodtCRDI") {
    returnValue <- stats::integrate(integrandBleichrodtCRDILog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["BleichrodtCRDI.lnk"]],
                                    s = results[["BleichrodtCRDI.s"]],
                                    beta = results[["BleichrodtCRDI.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "RodriguezLogue") {
    returnValue <- stats::integrate(integrandRodriguezLogueLog,
                                    lower = log10(min(dat$X)),
                                    upper = log10(max(dat$X)),
                                    lnK = results[["RodriguezLogue.lnk"]],
                                    beta = results[["RodriguezLogue.beta"]])$value/maximumArea

  } else if (results[["probable.model"]] == "Noise") {
    returnValue <- results[["Noise.mean"]]
  }

  returnValue
}




