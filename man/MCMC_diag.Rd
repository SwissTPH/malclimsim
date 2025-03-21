% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_diagnostics.R
\name{MCMC_diag}
\alias{MCMC_diag}
\title{MCMC Diagnostics}
\usage{
MCMC_diag(
  results,
  params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"),
  thin = 1
)
}
\arguments{
\item{results}{A list containing the MCMC output, typically with at least
\code{results[[1]]} for chain data and \code{results[[2]]} for trace plot data.
Should include an element \code{n_chains} for the number of chains.}

\item{params}{A character vector specifying which diagnostics to display.
Options include: "trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance".}

\item{thin}{An integer specifying thinning interval for the chains.
Thinning reduces the number of samples by keeping every \code{thin}th sample.}
}
\value{
Diagnostic outputs based on the selected \code{params} are printed or plotted.
}
\description{
Provides diagnostic information for MCMC results, including trace plots,
chain convergence statistics, correlation structure, effective sample size,
autocorrelation function, posterior quantiles, and acceptance rates.
Designed for users to assess MCMC chain convergence and parameter mixing.
}
\examples{
# Assuming 'results' is a valid MCMC result list with required structure:
MCMC_diag(results, params = c("trace", "ess"))
}
