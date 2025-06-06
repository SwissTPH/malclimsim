% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_visualizations.R
\name{plot_counterfactual_time_series}
\alias{plot_counterfactual_time_series}
\title{Plot Counterfactual Time Series Across Multiple Scenarios}
\usage{
plot_counterfactual_time_series(
  summaries_list,
  title = "Counterfactuals Over Time"
)
}
\arguments{
\item{summaries_list}{List of data frames returned by \code{summarize_simulation_ci()}.}

\item{title}{Title of the plot.}
}
\value{
A ggplot object.
}
\description{
Combines simulation summaries and generates a faceted time series plot.
}
