##
## Copyright 2021 Shawn Gilroy
##
## This file is part of discountingtools.
##
## discountingtools is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, version 2.
##
## discountingtools is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with discountingtools. If not, see <http://www.gnu.org/licenses/gpl-2.0.html>.

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

