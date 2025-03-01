% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counterfactual_analysis.R
\name{calculate_eff_smc_confidence_intervals}
\alias{calculate_eff_smc_confidence_intervals}
\title{Calculate Confidence Intervals for eff_SMC}
\usage{
calculate_eff_smc_confidence_intervals(
  results,
  param_inputs,
  model,
  start_date,
  end_date,
  n_samples = 100
)
}
\arguments{
\item{results}{The MCMC results object containing parameter samples.}

\item{param_inputs}{A list of parameter values for the model.}

\item{model}{The model function to simulate from.}

\item{start_date}{The start date for the simulation.}

\item{end_date}{The end date for the simulation.}

\item{n_samples}{The number of parameter sets to sample from the MCMC results.}
}
\value{
A data frame with confidence intervals for eff_SMC effectiveness.
}
\description{
This function calculates confidence intervals for eff_SMC by sampling
MCMC results, extracting parameter sets, updating them, and simulating outcomes.
}
\examples{
ci <- calculate_eff_smc_confidence_intervals(results, param_inputs, model, "2021-01-01", "2021-12-31", 100)
}
