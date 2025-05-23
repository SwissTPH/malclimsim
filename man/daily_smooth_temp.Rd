% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate_processing_functions.R
\name{daily_smooth_temp}
\alias{daily_smooth_temp}
\title{#' Daily Smoothed Temperature Function
#'
#' This function takes a data frame with monthly temperature data and produces a daily smoothed
#' temperature series using a smoothing spline. The resulting data frame includes predicted temperatures
#' for every day in the range from the earliest to the latest date in the original data.
#'
#' @param temp_df A data frame containing two columns:
#'   \describe{
#'     \item{Date}{A character or numeric column representing the year and month (e.g., "2024-11").}
#'     \item{Temperature}{A numeric column representing the temperature for the corresponding month.}
#'   }
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{Date}{A Date column representing daily dates.}
#'     \item{Temperature}{A numeric column with the smoothed daily temperature predictions.}
#'   }
#'
#' @examples
#' # Example usage:
#' temp_data <- data.frame(
#'   Date = c("2024-01", "2024-02", "2024-03"),
#'   Temperature = c(30.5, 28.0, 25.3)
#' )
#' daily_smooth_temp(temp_data)
#'
#' @importFrom stats smooth.spline predict
#' @export}
\usage{
daily_smooth_temp(temp_df)
}
\description{
#' Daily Smoothed Temperature Function
#'
#' This function takes a data frame with monthly temperature data and produces a daily smoothed
#' temperature series using a smoothing spline. The resulting data frame includes predicted temperatures
#' for every day in the range from the earliest to the latest date in the original data.
#'
#' @param temp_df A data frame containing two columns:
#'   \describe{
#'     \item{Date}{A character or numeric column representing the year and month (e.g., "2024-11").}
#'     \item{Temperature}{A numeric column representing the temperature for the corresponding month.}
#'   }
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{Date}{A Date column representing daily dates.}
#'     \item{Temperature}{A numeric column with the smoothed daily temperature predictions.}
#'   }
#'
#' @examples
#' # Example usage:
#' temp_data <- data.frame(
#'   Date = c("2024-01", "2024-02", "2024-03"),
#'   Temperature = c(30.5, 28.0, 25.3)
#' )
#' daily_smooth_temp(temp_data)
#'
#' @importFrom stats smooth.spline predict
#' @export
}
