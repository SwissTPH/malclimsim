#' Extract and Print Maximum Log-Likelihood and Log-Posterior
#'
#' This function extracts and optionally prints the maximum values of the log-likelihood
#' and log-posterior from the results of an MCMC run. It is useful for quickly evaluating
#' the fit of a model.
#'
#' @param results A list containing the MCMC results from an inference run,
#' including the `coda_pars` matrix with columns for `log_likelihood` and `log_posterior`.
#' @param print_ll Boolean. If `TRUE`, the maximum log-likelihood and log-posterior values
#' will be printed to the console (default: `TRUE`).
#'
#' @return A list with two named elements:
#' \describe{
#'   \item{max_log_likelihood}{The maximum value of the log-likelihood.}
#'   \item{max_log_posterior}{The maximum value of the log-posterior.}
#' }
#'
#' @examples
#' # Simulate results with coda_pars
#' results <- list(
#'   coda_pars = data.frame(
#'     log_likelihood = rnorm(1000, mean = -100, sd = 10),
#'     log_posterior = rnorm(1000, mean = -110, sd = 12)
#'   )
#' )
#'
#' # Extract and print maximum log-likelihood and log-posterior
#' max_ll_post(results)
#'
#' # Suppress printing and only return the values
#' max_ll_post(results, print_ll = FALSE)
#'
#' @export
max_ll_post <- function(results, print_ll = TRUE) {
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
  if(print_ll){
    cat("Maximum Log-Likelihood:", max_log_likelihood, "\n")
    cat("Maximum Log-Posterior:", max_log_posterior, "\n")
  }
  return(list(max_log_likelihood = max_log_likelihood,
              max_log_posterior = max_log_posterior))
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
#'
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
#' This function calculates various error metrics (MAE, RMSE, MAPE, R^2, Bias) between observed and simulated data.
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
    bias <- sum(residuals, na.rm = TRUE) # Sum of residuals

    data.frame(Group = group, MAE = mae, RMSE = rmse, MAPE = mape, R2 = r2, Bias = bias)
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
#' # Example using minimal mock data
#' obs_cases <- data.frame(
#'   date_ymd = as.Date("2022-01-01") + 0:6,
#'   inc_A = rpois(7, lambda = 10),
#'   inc_C = rpois(7, lambda = 12),
#'   inc = rpois(7, lambda = 15)
#' )
#'
#' simulated_df <- obs_cases
#'
#' assess_model_performance(obs_cases, simulated_df, date_column = "date_ymd", groups = c("inc_A", "inc_C", "inc"))
assess_model_performance <- function(observed_df, simulated_df, date_column, groups) {
  # Calculate error metrics
  error_metrics <- calculate_error_metrics(observed_df, simulated_df, date_column, groups)

  # Plot residuals
  residual_plot <- plot_residuals(observed_df, simulated_df, date_column, groups)

  return(list(error_metrics = error_metrics, residual_plot = residual_plot))
}

#' Calculate Bias Per Year
#'
#' This function calculates the bias between simulated and observed data for each year.
#'
#' @param simulated_df A data frame containing the simulated data with a `date_ymd` column and metrics such as `inc_A`, `inc_C`, and `inc`.
#' @param obs_cases A data frame containing the observed data with a `date_ymd` column and metrics such as `inc_A`, `inc_C`, and `inc`.
#'
#' @return A data frame summarizing the total bias for each metric (`inc_A`, `inc_C`, `inc`) per year.
#' The output includes the columns:
#' - `year`: The year.
#' - `total_bias_inc_A`: Total bias for `inc_A` in the year.
#' - `total_bias_inc_C`: Total bias for `inc_C` in the year.
#' - `total_bias_inc`: Total bias for `inc` in the year.
#'
#' @examples
#' simulated_df <- data.frame(
#'   date_ymd = as.Date(c("2018-01-01", "2018-01-02", "2019-01-01")),
#'   inc_A = c(5, 3, 4),
#'   inc_C = c(2, 1, 3),
#'   inc = c(7, 4, 7)
#' )
#' obs_cases <- data.frame(
#'   date_ymd = as.Date(c("2018-01-01", "2018-01-02", "2019-01-01")),
#'   inc_A = c(6, 4, 5),
#'   inc_C = c(3, 2, 4),
#'   inc = c(9, 6, 9)
#' )
#' calculate_bias_per_year(simulated_df, obs_cases)
#'
#' @export
calculate_bias_per_year <- function(simulated_df, obs_cases) {
  # Merge the simulated and observed data
  combined_df <- merge(simulated_df, obs_cases, by = "date_ymd", suffixes = c("_sim", "_obs"))

  # Extract the year from the date
  combined_df$year <- as.numeric(format(as.Date(combined_df$date_ymd), "%Y"))

  # Calculate bias for each metric
  combined_df$bias_inc_A <- combined_df$inc_A_obs - combined_df$inc_A_sim
  combined_df$bias_inc_C <- combined_df$inc_C_obs - combined_df$inc_C_sim
  combined_df$bias_inc <- combined_df$inc_obs - combined_df$inc_sim

  # Group by year and calculate the mean bias per year
  bias_per_year <- combined_df %>%
    group_by(year) %>%
    summarise(
      total_bias_inc_A = sum(bias_inc_A, na.rm = TRUE),
      total_bias_inc_C = sum(bias_inc_C, na.rm = TRUE),
      total_bias_inc = sum(bias_inc, na.rm = TRUE)
    )

  return(bias_per_year)
}

#' Compare Mean Bias for 2019 vs Other Years
#'
#' This function compares the bias of simulated data against observed data for the year 2019 (SMC)
#' versus other years (non-SMC).
#'
#' @param simulated_df A data frame containing the simulated data with a `date_ymd` column and metrics such as `inc_A`, `inc_C`, and `inc`.
#' @param obs_cases A data frame containing the observed data with a `date_ymd` column and metrics such as `inc_A`, `inc_C`, and `inc`.
#'
#' @return A list containing:
#' - `smc_bias`: A data frame with total bias for 2019 (SMC).
#' - `mean_non_smc_bias`: A data frame summarizing the mean bias for non-SMC years across all metrics.
#'
#' @examples
#' simulated_df <- data.frame(
#'   date_ymd = as.Date(c("2018-01-01", "2019-01-01", "2019-02-01")),
#'   inc_A = c(5, 4, 6),
#'   inc_C = c(2, 3, 5),
#'   inc = c(7, 6, 11)
#' )
#' obs_cases <- data.frame(
#'   date_ymd = as.Date(c("2018-01-01", "2019-01-01", "2019-02-01")),
#'   inc_A = c(6, 5, 7),
#'   inc_C = c(3, 4, 6),
#'   inc = c(9, 9, 13)
#' )
#' compare_bias_2019_vs_other_years(simulated_df, obs_cases)
#'
#' @export
compare_bias_2019_vs_other_years <- function(simulated_df, obs_cases) {
  # Get the bias per year
  bias_per_year <- calculate_bias_per_year(simulated_df, obs_cases)

  # Separate the biases into 2019 (SMC) and non-2019 (non-SMC)
  smc_bias <- bias_per_year %>% filter(year == 2019)
  non_smc_bias <- bias_per_year %>% filter(year != 2019)

  # Calculate the mean bias for non-SMC years
  mean_non_smc_bias <- non_smc_bias %>%
    summarise(
      mean_bias_inc_A = mean(total_bias_inc_A, na.rm = TRUE),
      mean_bias_inc_C = mean(total_bias_inc_C, na.rm = TRUE),
      mean_bias_inc = mean(total_bias_inc, na.rm = TRUE)
    )

  # Return a list comparing SMC (2019) vs non-SMC (other years)
  comparison <- list(
    smc_bias = smc_bias,
    mean_non_smc_bias = mean_non_smc_bias
  )

  return(comparison)
}
