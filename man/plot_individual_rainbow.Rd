% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dd_plot_ind.R
\name{plot_individual_rainbow}
\alias{plot_individual_rainbow}
\title{plot_individual_rainbow}
\usage{
plot_individual_rainbow(
  fittingObject,
  position0,
  ylab0,
  xlab0,
  logAxis,
  yMin,
  plotit
)
}
\arguments{
\item{fittingObject}{core fitting object}

\item{position0}{(char) position of legend}

\item{ylab0}{(char) y axis label}

\item{xlab0}{(char) x axis label}

\item{logAxis}{(char) axis designation}

\item{yMin}{(num) y axis lower limit}

\item{plotit}{(logical) bool of whether or not to print visual or output plotting frame}
}
\description{
This specific implementation shows cross-model fits, with series characterized by different models illustrated with different colors. A legend is also provided for convenience of interpretation.
}
\author{
Shawn Gilroy <sgilroy1@lsu.edu>
}
