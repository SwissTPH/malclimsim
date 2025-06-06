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

MCMC_diag(
  results,
  params = c("trace", "gelman", "corr", "ess", "acf", "quantiles", "acceptance"),
  thin = 1
)
}
\arguments{
\item{results}{A fitted MCMC result object (e.g., from \code{mcstate::pmcmc()}). Must contain a \code{pars} matrix and \code{n_chains} value.}

\item{params}{A character vector specifying which diagnostics to compute/plot. Options include:
\itemize{
\item \code{"trace"}: Trace plots of each parameter.
\item \code{"gelman"}: Gelman-Rubin convergence diagnostics (requires \code{n_chains > 1}).
\item \code{"corr"}: Pairwise correlation plots (calls \code{plot_corr()}).
\item \code{"ess"}: Effective sample size.
\item \code{"acf"}: Autocorrelation plots.
\item \code{"quantiles"}: Posterior quantile summaries.
\item \code{"acceptance"}: Acceptance rate (based on rejection rate).
}}

\item{thin}{Integer; optional thinning interval. If greater than 1, chains are thinned accordingly. Default is 1 (no thinning).}
}
\value{
Diagnostic outputs based on the selected \code{params} are printed or plotted.

A list containing the selected diagnostics. Output elements are named according to \code{params}, and may include:
\itemize{
\item \code{$trace}: The \code{mcmc.list} object used for plotting.
\item \code{$gelman}: Output from \code{gelman.diag()}.
\item \code{$corr}: Placeholder message indicating correlation plot was created.
\item \code{$ess}: Effective sample sizes from \code{effectiveSize()}.
\item \code{$acf}: Placeholder message indicating ACF was plotted.
\item \code{$quantiles}: Matrix of posterior quantiles (2.5\%, 50\%, 97.5\%).
\item \code{$acceptance}: Acceptance rates per chain.
}
}
\description{
Provides diagnostic information for MCMC results, including trace plots,
chain convergence statistics, correlation structure, effective sample size,
autocorrelation function, posterior quantiles, and acceptance rates.
Designed for users to assess MCMC chain convergence and parameter mixing.

Computes and optionally plots a suite of MCMC diagnostics including trace plots, Gelman-Rubin statistics,
autocorrelation, effective sample size, and posterior quantiles.
}
\details{
This function wraps common diagnostics from the \code{coda} package and custom correlation plotting into a single call.
Useful for quickly assessing convergence and chain quality after MCMC sampling.
}
\examples{
# Assuming 'results' is a valid MCMC result list with required structure:
MCMC_diag(results, params = c("trace", "ess"))
\dontrun{
diagnostics <- MCMC_diag(results, params = c("trace", "gelman", "ess", "quantiles"))
}

}
