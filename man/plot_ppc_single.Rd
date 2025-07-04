% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_visualizations.R
\name{plot_ppc_single}
\alias{plot_ppc_single}
\title{Plot Single Time Series Comparison (Uncomplicated or Severe Cases)}
\usage{
plot_ppc_single(
  plot_data,
  ci_data = NULL,
  obs_data = NULL,
  obs_col = NULL,
  plot_title = "Time Series Comparison",
  xlim = NULL,
  ylim = NULL,
  severity = "Uncomplicated",
  scale_severe = 1,
  scale_severe_by_year = NULL
)
}
\arguments{
\item{plot_data}{Data frame with columns: date_ymd (Date), value (numeric), label}

\item{ci_data}{Optional data frame: date_ymd, lower, upper, label}

\item{obs_data}{Optional observed data (must have date_ymd)}

\item{obs_col}{Name of obs column in obs_data}

\item{plot_title}{Title}

\item{xlim}{x‐axis limits (Date vector)}

\item{ylim}{y‐axis limits (numeric vector)}

\item{severity}{Either "Uncomplicated" or "Severe"}

\item{scale_severe}{Scaling factor for severe cases (used if no per-year mapping is provided)}

\item{scale_severe_by_year}{Named vector of multipliers by year (optional)}
}
\value{
ggplot2 plot object
}
\description{
Displays time series with optional credible intervals, observed points, and
consistent legend formatting across scenarios (e.g., With/Without SMC, June/July 2023, etc.)
}
