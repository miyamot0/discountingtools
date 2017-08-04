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
