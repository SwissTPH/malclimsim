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
