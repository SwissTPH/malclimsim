% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_diagnostics.R
\name{post_plot}
\alias{post_plot}
\title{Plot Posterior Distributions of Estimated Parameters}
\usage{
post_plot(
  results_list,
  params_to_estimate,
  dim_plot,
  show_true = TRUE,
  true_value = NULL,
  show_prior = FALSE,
  prior_n = 10000,
  title = "",
  run_labels = NULL,
  plot_type = "histogram"
)
}
\arguments{
\item{results_list}{A list containing MCMC results from multiple runs. Each element should contain posterior samples.}

\item{params_to_estimate}{A character vector of parameter names to include in the plot.}

\item{dim_plot}{A numeric vector of length 2 specifying the number of rows and columns for arranging plots.}

\item{show_true}{Logical; if \code{TRUE}, true parameter values will be indicated on the plots if provided.}

\item{true_value}{A named numeric vector of true values for parameters (optional). Must match \code{params_to_estimate} if used.}

\item{show_prior}{Logical; if \code{TRUE}, prior distributions will be overlaid on the plots.}

\item{prior_n}{Integer; the number of points to sample from the prior distribution.}

\item{title}{A character string specifying the title for the entire plot layout (optional).}

\item{run_labels}{A character vector specifying labels for each MCMC run (optional). Must match the length of \code{results_list}.}

\item{plot_type}{A character string specifying the type of plot: \code{"histogram"} or \code{"density"}.}
}
\value{
A combined plot object displaying posterior distributions for each parameter.
The plot includes optional lines for true values if \code{show_true = TRUE}, and prior distributions if \code{show_prior = TRUE}.
}
\description{
Generates density or histogram plots for each estimated parameter from MCMC output, with annotated quantiles
and optional true values if available. It helps users visualize the posterior distribution of model parameters.
}
