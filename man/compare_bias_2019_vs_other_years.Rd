% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted_model_assessment.R
\name{compare_bias_2019_vs_other_years}
\alias{compare_bias_2019_vs_other_years}
\title{Compare Mean Bias for 2019 vs Other Years}
\usage{
compare_bias_2019_vs_other_years(simulated_df, obs_cases)
}
\arguments{
\item{simulated_df}{A data frame containing the simulated data with a \code{date_ymd} column and metrics such as \code{inc_A}, \code{inc_C}, and \code{inc}.}

\item{obs_cases}{A data frame containing the observed data with a \code{date_ymd} column and metrics such as \code{inc_A}, \code{inc_C}, and \code{inc}.}
}
\value{
A list containing:
\itemize{
\item \code{smc_bias}: A data frame with total bias for 2019 (SMC).
\item \code{mean_non_smc_bias}: A data frame summarizing the mean bias for non-SMC years across all metrics.
}
}
\description{
This function compares the bias of simulated data against observed data for the year 2019 (SMC)
versus other years (non-SMC).
}
\examples{
simulated_df <- data.frame(
  date_ymd = as.Date(c("2018-01-01", "2019-01-01", "2019-02-01")),
  inc_A = c(5, 4, 6),
  inc_C = c(2, 3, 5),
  inc = c(7, 6, 11)
)
obs_cases <- data.frame(
  date_ymd = as.Date(c("2018-01-01", "2019-01-01", "2019-02-01")),
  inc_A = c(6, 5, 7),
  inc_C = c(3, 4, 6),
  inc = c(9, 9, 13)
)
compare_bias_2019_vs_other_years(simulated_df, obs_cases)

}
