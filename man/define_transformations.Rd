% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/running_inference.R
\name{define_transformations}
\alias{define_transformations}
\title{Define Transformation Function for MCMC Simulation}
\usage{
define_transformations(temp, c_R_D, SMC, decay, cov_SMC)
}
\arguments{
\item{temp}{Numeric vector or matrix of temperature values.}

\item{c_R_D}{Numeric vector representing case reduction from treatment or diagnostics.}

\item{SMC}{Numeric vector of SMC coverage.}

\item{decay}{Numeric vector for SMC efficacy decay over time.}

\item{cov_SMC}{Numeric vector of SMC target coverage.}
}
\value{
A transformation function compatible with \code{mcstate} model input.
}
\description{
Wraps the construction of a transformation function for model covariates (e.g., temperature, SMC coverage).
}
