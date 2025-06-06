% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_simulation.R
\name{sample_params}
\alias{sample_params}
\title{Sample Parameters from MCMC Results}
\usage{
sample_params(results, n)
}
\arguments{
\item{results}{The MCMC results object containing parameter samples.}

\item{n}{The number of samples to extract.}
}
\value{
A list of named vectors, each containing one set of sampled parameters.
}
\description{
This function randomly selects \code{n} rows from the MCMC parameter output and extracts the
parameter values from each row.
}
\examples{
# Assuming `results` contains the MCMC output
sampled_parameters <- sample_params(results, 100)
}
