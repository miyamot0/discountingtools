#' Exponential Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandExp <- function(x, lnK) { exp(-exp(lnK)*x) }

#' Exponential Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandExpLog <- function(x, lnK) { exp(-exp(lnK)*(10^x)) }

#' Hyperbolic Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandHyp <- function(x, lnK) { (1+exp(lnK)*x)^(-1) }

#' Hyperbolic Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandHypLog <- function(x, lnK) { (1+exp(lnK)*(10^x))^(-1) }

#' Beta Delta Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBetaDelta <- function(x, beta, delta) { beta*delta^x }

#' Beta Delta Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param beta fitted parameter
#' @param delta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandBetaDeltaLog <- function(x, beta, delta) { beta*delta^(10^x) }

#' Green & Myerson Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandMyerson <- function(x, lnK, s) { (1+exp(lnK)*x)^(-s) }

#' Green & Myerson Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandMyersonLog <- function(x, lnK, s) { (1+exp(lnK)*(10^x))^(-s) }

#' Rachlin Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRachlin <- function(x, lnK, s) { (1+exp(lnK)*(x^s))^(-1) }

#' Rachlin Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRachlinLog <- function(x, lnK, s) { (1+exp(lnK)*((10^x)^s))^(-1) }

#' Ebert & Prelec's Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandEbertPrelec <- function(x, lnK, s) {  exp(-(exp(lnK)*x)^s) }

#' Ebert & Prelec's ep Integrand helper (log10)
#'
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param s fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandEbertPrelecLog <- function(x, lnK, s) {  exp(-(exp(lnK)*(10^x))^s) }

#' Bleichrodt et al. Constant Relative Decreasing Impatience (CRDI) Integrand helper
#'
#' This integrand helper is a projection of the integrand with delays represented as normal
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
#' This integrand helper is a projection of the integrand with delays represented in the log base 10 scale
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
#' This integrand helper is a projection of the integrand with delays represented as normal
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRodriguezLogue <- function(x, lnK, beta) { (1 + x * exp(lnK))^(-beta / exp(lnK)) }

#' Rodriguez & Logue Integrand helper
#'
#' This integrand helper is a projection of the integrand (log10) with delays represented in the log base 10 scale
#'
#' @param x observation at point n (X)
#' @param lnK fitted parameter
#' @param beta fitted parameter
#' @author Shawn Gilroy <shawn.gilroy@temple.edu>
#' @return Numerical Integration Projection
integrandRodriguezLogueLog <- function(x, lnK, beta) { (1 + (10^x) * exp(lnK))^(-beta / exp(lnK)) }
