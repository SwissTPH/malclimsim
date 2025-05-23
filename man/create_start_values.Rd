% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_setup.R
\name{create_start_values}
\alias{create_start_values}
\title{Create Starting Values for Parameters}
\usage{
create_start_values(
  params_to_estimate,
  control_params,
  min_max_start_values = NULL,
  random = TRUE,
  seed = 10,
  model,
  param_inputs
)
}
\arguments{
\item{params_to_estimate}{A character vector of parameter names that are to be estimated.}

\item{control_params}{A list of control parameters (must include \code{n_chains} to specify the number of chains).}

\item{min_max_start_values}{A named list with minimum and maximum start values for each parameter.}

\item{random}{A logical value. If \code{TRUE}, random values are drawn from a uniform distribution; if \code{FALSE}, values are evenly spaced between the lower and upper bounds.}

\item{seed}{An integer specifying the seed for random number generation (default is 10).}

\item{model}{A model object with a \code{param()} method that returns all parameter names.}

\item{param_inputs}{A list of parameter values passed to initialize the model, used to extract parameters.}
}
\value{
A matrix of starting values for each parameter and chain.
}
\description{
This function generates starting values for parameters based on the specified method. It can either draw random values from a uniform distribution or assign values from a specified interval.
}
