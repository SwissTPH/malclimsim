% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_simulation.R
\name{run_mcmc_simulations}
\alias{run_mcmc_simulations}
\title{Run Model Simulations for Each Parameter Set}
\usage{
run_mcmc_simulations(
  model,
  param_inputs,
  param_samples,
  start_date,
  end_date,
  prewarm_years = 2,
  days_per_year = 360
)
}
\arguments{
\item{model}{A function representing the epidemiological model.}

\item{param_samples}{A matrix of parameter sets sampled from the MCMC posterior.}

\item{start_date}{The simulation start date.}

\item{end_date}{The simulation end date.}

\item{prewarm_years}{The number of years to prewarm the model.}

\item{days_per_year}{The number of days in a simulation year.}
}
\value{
A list of data frames containing simulation results for each sampled parameter set.
}
\description{
This function runs the epidemiological model using different parameter sets sampled from the MCMC posterior.
}
\examples{
simulations <- run_mcmc_simulations(model, sampled_params, "2014-01-01", "2022-12-31", 2, 360)
}
