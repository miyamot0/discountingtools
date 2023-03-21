#' Generalized residual call
#'
#' General, shared method for coordinating nls.lm fitting calls. Routes a supplied "valueFunction" with observed data and supplied parameters.
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#'
#' @return residual value of referenced function
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
residualFunction <- function(params, x, value, valueFunction, jacobianFunction)
{
  value - do.call("valueFunction", c(list(x = x), as.list(params)))
}

#' Generalized Jacobian call
#'
#' General, shared method for constructing the Jacobian matrix. Routes a supplied "jacobianFunction" with pre-computed derivatives to construct matrix with observed data and supplied parameters.
#'
#' @param params model parameters
#' @param x observation at point n (X)
#' @param value observation at point n (Y)
#' @param valueFunction function to get projected value
#' @param jacobianFunction function to create jacobian
#'
#' @return difference value for jacobian
#' @author Shawn Gilroy <sgilroy1@lsu.edu>
jacobianMatrix <- function(params, x, value, valueFunction, jacobianFunction)
{
  -do.call("jacobianFunction", c(list(x = x), as.list(params)))
}



#' minpack.lm logLik hack
#'
#' This function constructs a class, derived from an nls.lm object, similar to that of the logLik function in nls. This allows for native calls of the AIC and BIC functions from stats, using nls.lm fit objects.
#'
#' @param fit nls.lm fitted model
#' @param REML determine whether or not to use ML (FALSE by default)
#' @param ... inherit other args as necessary
#'
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

#' Workaround for varying bx for hypergeometric series
#'
#' Credit: Stéphane Laurent <laurent_step at yahoo.fr>
#' Source: https://stats.stackexchange.com/questions/33451/computation-of-hypergeometric-function-in-r
#' Licensed CC-BY-SA 3.0, as Per SA Guidelines
#'
#' @param a param
#' @param b param
#' @param c param
#' @param x param
#'
#' @author Stéphane Laurent <laurent_step at yahoo.fr>
#' @importFrom gsl hyperg_2F1
gauss_2F1 <- function(a, b, c, x){
  if (x >= 0 & x < 1) {
    hyperg_2F1(a, b, c, x)
  } else {
    hyperg_2F1(c - a, b, c, 1 - 1 / (1 - x)) / (1 - x)^b
  }
}
