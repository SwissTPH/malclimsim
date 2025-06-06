% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counterfactual_analysis.R
\name{calculate_cases_averted}
\alias{calculate_cases_averted}
\title{Calculate Cases Averted by SMC}
\usage{
calculate_cases_averted(with_smc_df, without_smc_df)
}
\arguments{
\item{with_smc_df}{Data frame containing the incidence data from the model simulated with SMC.}

\item{without_smc_df}{Data frame containing the incidence data from the model simulated without SMC.}
}
\value{
A list containing the number of cases averted and the percentage of cases averted.
}
\description{
This function calculates the number and percent of malaria cases averted due to SMC intervention by comparing
the simulated scenarios with and without SMC.
}
\examples{
cases_averted <- calculate_cases_averted(sim_with_smc, sim_without_smc)
}
