#' Print Maximum Log-Likelihood and Log-Posterior
#'
#' This function extracts and prints the maximum values of the log-likelihood and
#' log-posterior from the results of an MCMC run. Useful for quick evaluation of model fit.
#'
#' @param results A list containing the MCMC results from an inference run,
#' including the `coda_pars` matrix with columns for `log_likelihood` and `log_posterior`.
#'
#' @return NULL. The function prints the maximum values of log-likelihood and log-posterior.
#' @export
#'
#' @examples
#' # Assuming `results` is the result of an inference run with MCMC sampling:
#' max_ll_post(results)
max_ll_post <- function(results) {
  # Check if 'coda_pars' exists and contains required columns
  if (!"coda_pars" %in% names(results)) {
    stop("The results list must contain 'coda_pars'.")
  }

  if (!all(c("log_likelihood", "log_posterior") %in% colnames(results$coda_pars))) {
    stop("'coda_pars' must contain 'log_likelihood' and 'log_posterior' columns.")
  }

  # Extract MCMC posterior results
  coda_pars <- results$coda_pars

  # Find maximum values for log_likelihood and log_posterior
  max_log_likelihood <- max(coda_pars[ , "log_likelihood"], na.rm = TRUE)
  max_log_posterior <- max(coda_pars[ , "log_posterior"], na.rm = TRUE)

  # Print the results
  cat("Maximum Log-Likelihood:", max_log_likelihood, "\n")
  cat("Maximum Log-Posterior:", max_log_posterior, "\n")
}

#' Calculate Mean Absolute Error (MAE)
#'
#' Calculates the Mean Absolute Error between observed and simulated values.
#'
#' @param observed A numeric vector of observed values.
#' @param simulated A numeric vector of simulated values.
#' @return The MAE value.
#' @examples
#' mae_value <- calculate_mae(observed, simulated)
calculate_mae <- function(observed, simulated) {
  mae <- mean(abs(observed - simulated))
  return(mae)
}

#' Calculate Root Mean Squared Error (RMSE)
#'
#' Calculates the Root Mean Squared Error between observed and simulated values.
#'
#' @param observed A numeric vector of observed values.
#' @param simulated A numeric vector of simulated values.
#' @return The RMSE value.
#' @examples
#' rmse_value <- calculate_rmse(observed, simulated)
calculate_rmse <- function(observed, simulated) {
  rmse <- sqrt(mean((observed - simulated)^2))
  return(rmse)
}

#' Calculate Mean Absolute Percentage Error (MAPE)
#'
#' Calculates the Mean Absolute Percentage Error between observed and simulated values.
#'
#' @param observed A numeric vector of observed values.
#' @param simulated A numeric vector of simulated values.
#' @return The MAPE value as a percentage.
#' @examples
#' mape_value <- calculate_mape(observed, simulated)
calculate_mape <- function(observed, simulated) {
  mape <- mean(abs((observed - simulated) / observed)) * 100
  return(mape)
}

#' Calculate Bias
#'
#' Calculates the bias, which shows if the model tends to underpredict or overpredict.
#' Positive bias means overprediction, while negative bias means underprediction.
#'
#' @param observed A numeric vector of observed values.
#' @param simulated A numeric vector of simulated values.
#' @return The bias value.
#' @examples
#' bias_value <- calculate_bias(observed, simulated)
calculate_bias <- function(observed, simulated) {
  bias <- mean(simulated - observed)
  return(bias)
}

#' Plot Residuals Between Observed and Simulated Data
#'
#' Plots the residuals (difference between observed and simulated values) over time.
#'
#' @param dates A vector of dates corresponding to the observations.
#' @param observed A numeric vector of observed values.
#' @param simulated A numeric vector of simulated values.
#' @return A ggplot object showing residuals over time.
#' @examples
#' residual_plot <- plot_residuals(dates, observed, simulated)
plot_residuals <- function(dates, observed, simulated) {
  residuals <- observed - simulated
  library(ggplot2)
  p <- ggplot(data.frame(date = dates, residuals = residuals), aes(x = date, y = residuals)) +
    geom_line(color = "darkred", size = 1) +
    labs(title = "Residuals Over Time", x = "Date", y = "Residuals (Observed - Simulated)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

#' Assess Model Performance
#'
#' Calculates various error metrics (MAE, RMSE, MAPE, Bias) and plots residuals between
#' observed and simulated values.
#'
#' @param observed_df A data frame containing observed values (columns: `inc_A`, `inc_C`, `inc`).
#' @param simulated_df A data frame containing simulated values (columns: `inc_A`, `inc_C`, `inc`).
#' @param date_column The name of the column representing dates in `observed_df` and `simulated_df`.
#' @return A list of error metrics for each group (`inc_A`, `inc_C`, `inc`) and residual plots.
#' @examples
#' error_metrics <- assess_model_performance(observed_df, simulated_df, date_column = "date_ymd")
assess_model_performance <- function(observed_df, simulated_df, date_column = "date_ymd") {

  # Merging observed and simulated data based on the date column
  combined_data <- merge(observed_df, simulated_df, by = date_column, suffixes = c("_obs", "_sim"))

  # Error metrics for each group (inc_A, inc_C, inc)
  groups <- c("inc_A", "inc_C", "inc")

  error_metrics <- lapply(groups, function(group) {
    observed <- combined_data[[paste0(group, "_obs")]]
    simulated <- combined_data[[paste0(group, "_sim")]]

    mae <- calculate_mae(observed, simulated)
    rmse <- calculate_rmse(observed, simulated)
    mape <- calculate_mape(observed, simulated)
    bias <- calculate_bias(observed, simulated)

    list(
      group = group,
      MAE = mae,
      RMSE = rmse,
      MAPE = mape,
      Bias = bias
    )
  })

  # Create residual plots for each group
  residual_plots <- lapply(groups, function(group) {
    observed <- combined_data[[paste0(group, "_obs")]]
    simulated <- combined_data[[paste0(group, "_sim")]]
    plot_residuals(combined_data[[date_column]], observed, simulated)
  })

  # Return a list with metrics and plots
  return(list(
    error_metrics = error_metrics,
    residual_plots = residual_plots
  ))
}

#' Calculate Model Errors
#'
#' This function calculates various error metrics (MAE, RMSE, MAPE, R²) for observed and simulated data.
#'
#' @param observed_df Data frame containing the observed data.
#' @param simulated_df Data frame containing the simulated data.
#' @param date_column The name of the date column in both data frames.
#' @return A data frame containing error metrics for each group.
#' @export
#' @examples
#' error_metrics <- calculate_model_errors(obs_cases, simulated_df, date_column = "date_ymd")
calculate_model_errors <- function(observed_df, simulated_df, date_column) {

  # Ensure date compatibility and merge observed and simulated data
  observed_df[[date_column]] <- as.Date(observed_df[[date_column]])
  simulated_df[[date_column]] <- as.Date(simulated_df[[date_column]])
  combined_data <- merge(observed_df, simulated_df, by = date_column, suffixes = c("_obs", "_sim"))

  # Calculate residuals for each group: (observed - simulated)
  combined_data$residual_inc_A <- combined_data$inc_A_obs - combined_data$inc_A_sim
  combined_data$residual_inc_C <- combined_data$inc_C_obs - combined_data$inc_C_sim
  combined_data$residual_inc <- combined_data$inc_obs - combined_data$inc_sim

  # Calculate error metrics (MAE, RMSE, MAPE, R²) for each group
  mae_inc_A <- mean(abs(combined_data$residual_inc_A), na.rm = TRUE)
  mae_inc_C <- mean(abs(combined_data$residual_inc_C), na.rm = TRUE)
  mae_inc <- mean(abs(combined_data$residual_inc), na.rm = TRUE)

  rmse_inc_A <- sqrt(mean(combined_data$residual_inc_A^2, na.rm = TRUE))
  rmse_inc_C <- sqrt(mean(combined_data$residual_inc_C^2, na.rm = TRUE))
  rmse_inc <- sqrt(mean(combined_data$residual_inc^2, na.rm = TRUE))

  mape_inc_A <- mean(abs(combined_data$residual_inc_A / combined_data$inc_A_obs), na.rm = TRUE) * 100
  mape_inc_C <- mean(abs(combined_data$residual_inc_C / combined_data$inc_C_obs), na.rm = TRUE) * 100
  mape_inc <- mean(abs(combined_data$residual_inc / combined_data$inc_obs), na.rm = TRUE) * 100

  r2_inc_A <- 1 - (sum(combined_data$residual_inc_A^2, na.rm = TRUE) /
                     sum((combined_data$inc_A_obs - mean(combined_data$inc_A_obs, na.rm = TRUE))^2, na.rm = TRUE))

  r2_inc_C <- 1 - (sum(combined_data$residual_inc_C^2, na.rm = TRUE) /
                     sum((combined_data$inc_C_obs - mean(combined_data$inc_C_obs, na.rm = TRUE))^2, na.rm = TRUE))

  r2_inc <- 1 - (sum(combined_data$residual_inc^2, na.rm = TRUE) /
                   sum((combined_data$inc_obs - mean(combined_data$inc_obs, na.rm = TRUE))^2, na.rm = TRUE))

  # Compile error metrics into a data frame
  error_metrics <- data.frame(
    Group = c(">=5 years old (inc_A)", "<5 years old (inc_C)", "Total (inc)"),
    MAE = c(mae_inc_A, mae_inc_C, mae_inc),
    RMSE = c(rmse_inc_A, rmse_inc_C, rmse_inc),
    MAPE = c(mape_inc_A, mape_inc_C, mape_inc),
    R2 = c(r2_inc_A, r2_inc_C, r2_inc)
  )

  return(error_metrics)
}




