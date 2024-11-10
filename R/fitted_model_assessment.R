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

#' Plot Residuals of Observed vs Simulated Data
#'
#' This function creates residual plots for observed and simulated data, allowing the user to select specific groups.
#'
#' @param observed_df Data frame containing the observed data.
#' @param simulated_df Data frame containing the simulated data.
#' @param date_column The name of the date column in both data frames.
#' @param groups A character vector specifying which groups to include in the plot (`inc_A`, `inc_C`, `inc`).
#' @return A ggplot object displaying the residual plots.
#' @export
#' @examples
#' residual_plot <- plot_residuals(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C"))
plot_residuals <- function(observed_df, simulated_df, date_column, groups = c("inc_A", "inc_C", "inc")) {

  # Ensure date compatibility and merge observed and simulated data
  observed_df[[date_column]] <- as.Date(observed_df[[date_column]])
  simulated_df[[date_column]] <- as.Date(simulated_df[[date_column]])
  combined_data <- merge(observed_df, simulated_df, by = date_column, suffixes = c("_obs", "_sim"))

  # Calculate residuals for each group: (observed - simulated)
  combined_data$residual_inc_A <- combined_data$inc_A_obs - combined_data$inc_A_sim
  combined_data$residual_inc_C <- combined_data$inc_C_obs - combined_data$inc_C_sim
  combined_data$residual_inc <- combined_data$inc_obs - combined_data$inc_sim

  # Prepare data for plotting residuals for each selected group
  residuals_long <- data.frame(
    date_ymd = combined_data[[date_column]],
    inc_A = combined_data$residual_inc_A,
    inc_C = combined_data$residual_inc_C,
    inc = combined_data$residual_inc
  ) %>%
    tidyr::pivot_longer(cols = c("inc_A", "inc_C", "inc"),
                        names_to = "Group",
                        values_to = "Residual") %>%
    dplyr::filter(Group %in% groups)  # Filter selected groups

  # Create facet plot for residuals
  residual_plot <- ggplot(residuals_long, aes(x = date_ymd, y = Residual)) +
    geom_line(color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Group, scales = "free_y", ncol = 1,
               labeller = labeller(Group = c(
                 "inc_A" = ">=5 years old (inc_A)",
                 "inc_C" = "<5 years old (inc_C)",
                 "inc" = "Total incidence (inc)"
               ))) +
    labs(
      title = "Residuals: Observed vs Simulated Incidence",
      x = "Date",
      y = "Residuals (Observed - Simulated)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14, face = "bold")
    )

  return(residual_plot)
}



#' Calculate Error Metrics for Model Assessment
#'
#' This function calculates various error metrics (MAE, RMSE, MAPE, R^2) between observed and simulated data.
#'
#' @param observed_df Data frame containing the observed data.
#' @param simulated_df Data frame containing the simulated data.
#' @param date_column The name of the date column in both data frames.
#' @param groups A character vector specifying which groups to calculate error metrics for (e.g., `c("inc_A", "inc_C", "inc")`).
#' @return A data frame containing the calculated error metrics for each group.
#' @export
#' @examples
#' error_metrics <- calculate_error_metrics(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C", "inc"))
calculate_error_metrics <- function(observed_df, simulated_df, date_column, groups) {
  # Check if date_column exists in both data frames
  if (!(date_column %in% names(observed_df)) || !(date_column %in% names(simulated_df))) {
    stop("The specified date_column does not exist in one or both data frames.")
  }

  # Ensure date compatibility and merge observed and simulated data
  observed_df[[date_column]] <- as.Date(observed_df[[date_column]])
  simulated_df[[date_column]] <- as.Date(simulated_df[[date_column]])
  combined_data <- merge(observed_df, simulated_df, by = date_column, suffixes = c("_obs", "_sim"))

  error_metrics_list <- lapply(groups, function(group) {
    residuals <- combined_data[[paste0(group, "_obs")]] - combined_data[[paste0(group, "_sim")]]
    mae <- mean(abs(residuals), na.rm = TRUE)
    rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
    mape <- mean(abs(residuals / combined_data[[paste0(group, "_obs")]]), na.rm = TRUE) * 100
    r2 <- 1 - (sum(residuals^2, na.rm = TRUE) / sum((combined_data[[paste0(group, "_obs")]] - mean(combined_data[[paste0(group, "_obs")]], na.rm = TRUE))^2, na.rm = TRUE))

    data.frame(Group = group, MAE = mae, RMSE = rmse, MAPE = mape, R2 = r2)
  })

  # Combine error metrics into a data frame
  error_metrics_df <- do.call(rbind, error_metrics_list)
  return(error_metrics_df)
}

#' Assess Model Performance
#'
#' This function assesses the performance of a model by calculating error metrics and plotting residuals.
#'
#' @param observed_df Data frame containing the observed data.
#' @param simulated_df Data frame containing the simulated data.
#' @param date_column The name of the date column in both data frames.
#' @param groups A character vector specifying which groups to assess (e.g., `c("inc_A", "inc_C", "inc")`).
#' @return A list containing a data frame of error metrics and a ggplot object of residual plots.
#' @export
#' @examples
#' performance_results <- assess_model_performance(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C", "inc"))
assess_model_performance <- function(observed_df, simulated_df, date_column, groups) {
  # Calculate error metrics
  error_metrics <- calculate_error_metrics(observed_df, simulated_df, date_column, groups)

  # Plot residuals
  residual_plot <- plot_residuals(observed_df, simulated_df, date_column, groups)

  return(list(error_metrics = error_metrics, residual_plot = residual_plot))
}



