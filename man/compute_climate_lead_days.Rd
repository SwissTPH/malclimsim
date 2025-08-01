% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{compute_climate_lead_days}
\alias{compute_climate_lead_days}
\title{Compute how many days Climate must lead SMC}
\usage{
compute_climate_lead_days(smc_df, clim_df)
}
\arguments{
\item{smc_df}{Data.frame with a Date column \code{dates}.}

\item{clim_df}{Data.frame with a Date column \code{dates}.}
}
\value{
Integer >= 0: days by which Climate must precede SMC.
}
\description{
Given SMC and Climate data.frames (both with a Date column \code{dates}),
returns the integer number of days by which the Climate series must
start earlier than the SMC series start.
}
