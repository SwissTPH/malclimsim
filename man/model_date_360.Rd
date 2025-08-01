% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{model_date_360}
\alias{model_date_360}
\title{Convert model day index to calendar date assuming 360-day years}
\usage{
model_date_360(day_index, start_year)
}
\arguments{
\item{day_index}{Integer vector of day indices (starting from 1).
Represents the number of days since the start of the model.}

\item{start_year}{Integer. The starting calendar year for day_index = 1.}
}
\value{
A data frame with columns:
\itemize{
\item \code{year}: The calendar year.
\item \code{month}: The calendar month (1-12).
\item \code{day}: The day of the month (1-30).
}
}
\description{
Converts a given day index (e.g. from a model using 360-day years) into
a calendar-like date assuming each year has 12 months of 30 days.
}
