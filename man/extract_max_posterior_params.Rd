% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_simulation.R
\name{extract_max_posterior_params}
\alias{extract_max_posterior_params}
\title{Extract Parameters with Maximum Log Posterior}
\usage{
extract_max_posterior_params(results)
}
\arguments{
\item{results}{The MCMC results object containing parameter samples.}
}
\value{
A named vector of parameter values corresponding to the maximum log posterior.
}
\description{
This function finds and extracts the parameter set corresponding to the maximum log posterior value
from the MCMC results.
}
\examples{
# Assuming `results` contains the MCMC output
max_posterior_params <- extract_max_posterior_params(results)
}
