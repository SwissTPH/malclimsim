#' Set eff_SMC Parameter to Zero
#'
#' This function takes a list of model parameters and sets the `eff_SMC` parameter to zero to simulate
#' a scenario without Seasonal Malaria Chemoprevention (SMC).
#'
#' @param param_inputs A list of parameter values for the model.
#' @return A list with the `eff_SMC` parameter set to zero.
#' @export
#' @examples
#' no_smc_params <- set_eff_smc_to_zero(param_inputs)
set_eff_smc_to_zero <- function(param_inputs) {
  # Set the SMC effectiveness parameter to zero
  param_inputs$eff_SMC <- 0
  return(param_inputs)
}

#' Simulate Model with SMC
#'
#' This function simulates the model using the parameters that maximize the log posterior, which include SMC.
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @return A data frame containing the simulation output with SMC.
#' @export
simulate_with_smc <- function(results, start_date, end_date, model,
                              mu_transform_C = NULL, mu_transform_A = NULL,
                              covariate_matrix = NULL) {
  # Extract the best-fit parameters
  max_posterior_params <- extract_max_posterior_params(results)

  # Update the parameters to set eff_SMC to zero
  param_inputs <- results$param_inputs
  updated_params <- update_param_list(param_inputs, max_posterior_params)
  no_smc_params <- set_eff_smc_to_zero(updated_params)

  # Run the model simulation without SMC
  simulation_output <- data_sim(model, no_smc_params, start_date, end_date,
                                month = TRUE, round = FALSE, save = FALSE,
                                month_unequal_days = FALSE,
                                mu_transform_C = mu_transform_C,
                                mu_transform_A = mu_transform_A,
                                covariate_matrix = covariate_matrix)
  return(sim_data)
}

#' Simulate Model Without SMC
#'
#' This function simulates the model with `eff_SMC` set to zero, representing a scenario without SMC intervention.
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @return A data frame containing the simulation output without SMC.
#' @export
simulate_without_smc <- function(results, start_date, end_date, model,
                                 mu_transform_C = NULL, mu_transform_A = NULL,
                                 covariate_matrix = NULL) {
  # Extract the best-fit parameters
  max_posterior_params <- extract_max_posterior_params(results)

  # Update the parameters to set eff_SMC to zero
  param_inputs <- results$param_inputs
  updated_params <- update_param_list(param_inputs, max_posterior_params)
  no_smc_params <- set_eff_smc_to_zero(updated_params)

  # Run the model simulation without SMC
  simulation_output <- data_sim(model, no_smc_params, start_date, end_date,
                                month = TRUE, round = FALSE, save = FALSE,
                                month_unequal_days = FALSE,
                                mu_transform_C = mu_transform_C,
                                mu_transform_A = mu_transform_A,
                                covariate_matrix = covariate_matrix)

  return(simulation_output)
}

#' Calculate Cases Averted by SMC
#'
#' This function calculates the number and percent of malaria cases averted due to SMC intervention by comparing
#' the simulated scenarios with and without SMC.
#'
#' @param with_smc_df Data frame containing the incidence data from the model simulated with SMC.
#' @param without_smc_df Data frame containing the incidence data from the model simulated without SMC.
#' @return A list containing the number of cases averted and the percentage of cases averted.
#' @export
#' @examples
#' cases_averted <- calculate_cases_averted(sim_with_smc, sim_without_smc)
calculate_cases_averted <- function(with_smc_df, without_smc_df) {

  # Ensure both data frames have the same dates for proper comparison
  if (!all(with_smc_df$date_ymd == without_smc_df$date_ymd)) {
    stop("Dates in with_smc_df and without_smc_df do not match.")
  }

  # Calculate the total number of cases for each scenario
  total_cases_with_smc <- sum(with_smc_df$inc, na.rm = TRUE)
  total_cases_without_smc <- sum(without_smc_df$inc, na.rm = TRUE)

  # Calculate the number and percentage of cases averted
  cases_averted <- total_cases_without_smc - total_cases_with_smc
  percent_averted <- (cases_averted / total_cases_without_smc) * 100

  return(list(
    cases_averted = cases_averted,
    percent_averted = percent_averted
  ))
}

#' Set SMC Coverage Level
#'
#' This function takes a list of model parameters and sets the `cov_SMC` parameter to a user-defined level.
#'
#' @param param_inputs A list of parameter values for the model.
#' @param coverage_level A numeric value (between 0 and 1) representing the desired coverage level. Default is 1 (100% coverage).
#' @return A list with the `cov_SMC` parameter updated to the specified coverage level.
#' @export
#' @examples
#' full_coverage_params <- set_smc_coverage(param_inputs, coverage_level = 1)
set_smc_coverage <- function(param_inputs, coverage_level = 1) {
  # Update all values of cov_SMC to the desired level
  param_inputs$cov_SMC <- rep(coverage_level, length(param_inputs$cov_SMC))
  return(param_inputs)
}


#' Simulate Model with Full SMC Coverage
#'
#' This function simulates the model with `cov_SMC` set to a user-defined level (default is 100% coverage).
#'
#' @param results The MCMC results object containing the parameter samples.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param model The model function to simulate from.
#' @param coverage_level A numeric value between 0 and 1 for the desired SMC coverage. Default is 1 (100% coverage).
#' @return A data frame containing the simulation output with full or specified SMC coverage.
#' @export
#' @examples
#' sim_full_coverage <- simulate_with_full_coverage(results, start_date, end_date, model, coverage_level = 1)
simulate_with_full_coverage <- function(results, start_date, end_date, model, coverage_level = 1) {
  # Extract the best-fit parameters
  max_posterior_params <- extract_max_posterior_params(results)

  # Update the parameters with the extracted values
  param_inputs <- results$param_inputs
  updated_params <- update_param_list(param_inputs, max_posterior_params)

  # Set the SMC coverage level to full (or specified level)
  updated_params <- set_smc_coverage(updated_params, coverage_level)

  # Run the model simulation with the updated coverage level
  simulation_output <- data_sim(model, updated_params, start_date, end_date,
                                month = TRUE, round = FALSE, save = FALSE,
                                month_unequal_days = FALSE)

  return(simulation_output)
}


#' Calculate Confidence Intervals for eff_SMC
#'
#' This function calculates confidence intervals for eff_SMC by sampling
#' MCMC results, extracting parameter sets, updating them, and simulating outcomes.
#'
#' @param results The MCMC results object containing parameter samples.
#' @param param_inputs A list of parameter values for the model.
#' @param model The model function to simulate from.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param n_samples The number of parameter sets to sample from the MCMC results.
#' @return A data frame with confidence intervals for eff_SMC effectiveness.
#' @export
#' @examples
#' ci <- calculate_eff_smc_confidence_intervals(results, param_inputs, model, "2021-01-01", "2021-12-31", 100)
calculate_eff_smc_confidence_intervals <- function(results, param_inputs, model, start_date, end_date, n_samples = 100) {
  # Sample parameters from MCMC results
  sampled_params <- sample_params(results, n_samples)

  # Simulate outcomes with and without SMC for each parameter set
  simulations_with_smc <- simulate_models(model, param_inputs, sampled_params, start_date, end_date)
  simulations_without_smc <- simulate_models(model, param_inputs,
                                             lapply(sampled_params, set_eff_smc_to_zero),
                                             start_date, end_date)

  # Calculate cases averted for each simulation
  cases_averted_list <- mapply(calculate_cases_averted,
                               simulations_with_smc,
                               simulations_without_smc,
                               SIMPLIFY = FALSE)

  # Extract effectiveness estimates
  eff_smc_estimates <- sapply(cases_averted_list, function(x) x$percent_averted)

  # Calculate confidence intervals
  ci <- quantile(eff_smc_estimates, probs = c(0.025, 0.975))

  return(data.frame(
    lower_ci = ci[1],
    upper_ci = ci[2],
    median = median(eff_smc_estimates)
  ))
}


#' Calculate Outcome Across Matched Simulation Lists
#'
#' @param o1 List of data frames (e.g., with SMC).
#' @param o2 List of data frames (e.g., no SMC).
#' @param outcome_fn Function that takes two data frames and returns a scalar outcome.
#'
#' @return Numeric vector of outcomes (e.g., effectiveness or cases averted) for each parameter set.
#' @export
calculate_estimate <- function(o1, o2, outcome_fn) {
  if (!(is.null(o1) || is.null(o2))) {
    if (length(o1) != length(o2)) stop("Lists must be of equal length.")
  }

  estimates <- vapply(seq_along(o1), function(i) {
    result <- tryCatch({
      val <- outcome_fn(o1[[i]], o2[[i]])
      if (is.na(val) || is.nan(val)) {
        warning(sprintf("NA or NaN in outcome for simulation %d", i))
      }
      val
    }, error = function(e) {
      warning(sprintf("Error in outcome_fn for simulation %d: %s", i, e$message))
      NA_real_
    })
    result
  }, numeric(1))

  return(estimates)
}


#' Plot Histogram of Estimated Effects with Confidence Interval and Trimmed X-Axis
#'
#' @param estimates Numeric vector of estimates (e.g., cases averted).
#' @param x_label Label for the x-axis.
#' @param title Title for the plot.
#' @param bins Number of histogram bins (default = 30).
#' @param ci_level Confidence interval level (default = 0.95).
#'
#' @return A ggplot object showing histogram, density, and CI lines.
#' @export
plot_estimate_distribution <- function(estimates,
                                       x_label = "Cases Averted",
                                       title = "Posterior Distribution of Estimated Cases Averted",
                                       bins = 30,
                                       ci_level = 0.95) {
  df <- data.frame(estimate = estimates)

  # Central estimates
  median_est <- median(estimates)
  lower_ci <- quantile(estimates, probs = (1 - ci_level) / 2, na.rm = TRUE)
  upper_ci <- quantile(estimates, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)

  # Axis trimming (0.25% to 99.75%)
  x_min <- quantile(estimates, 0.0025, na.rm = TRUE)
  x_max <- quantile(estimates, 0.9975, na.rm = TRUE)

  ggplot(df, aes(x = estimate)) +
    geom_histogram(aes(y = ..density..), bins = bins, fill = "#2C3E50", alpha = 0.85, color = "white") +
    geom_density(color = "#E74C3C", size = 1.3, adjust = 1.2) +
    geom_vline(xintercept = median_est, linetype = "dashed", color = "#2980B9", size = 1.1) +
    geom_vline(xintercept = lower_ci, linetype = "dotted", color = "#27AE60", size = 1) +
    geom_vline(xintercept = upper_ci, linetype = "dotted", color = "#27AE60", size = 1) +
    annotate("text", x = median_est, y = Inf,
             label = paste("Median =", round(median_est, 1)),
             vjust = -0.7, hjust = -0.1, size = 5, fontface = "bold", color = "#2980B9") +
    annotate("text", x = lower_ci, y = Inf,
             label = paste0(round(100 * (1 - ci_level) / 2), "% = ", round(lower_ci, 1)),
             vjust = -1.5, hjust = 1, size = 4.5, color = "#27AE60") +
    annotate("text", x = upper_ci, y = Inf,
             label = paste0(round(100 * (1 + ci_level) / 2), "% = ", round(upper_ci, 1)),
             vjust = -1.5, hjust = 0, size = 4.5, color = "#27AE60") +
    coord_cartesian(xlim = c(x_min, x_max)) +
    labs(
      title = title,
      x = x_label,
      y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 12)
    )
}

#' Evaluate Multiple SMC Scenarios
#'
#' Runs model simulations for a list of SMC coverage patterns and summarizes outcomes.
#'
#' @param patterns A named list of 12-element binary vectors representing monthly SMC coverage patterns.
#' @param smc_day_of_month Day of month for simulated SMC coverage start (default = 1).
#' @param model A compiled model object used for simulations.
#' @param param_inputs A named list of baseline parameter values.
#' @param param_samples A matrix of sampled parameter sets (rows = samples).
#' @param start_date Start date for the simulation.
#' @param end_date End date for the simulation.
#' @param avg_cov Average SMC coverage to apply during active months.
#' @param years A vector of years to apply SMC coverage.
#' @param exclude_years Years to exclude from summarizing coverage (default = 2023).
#' @param mu_transform_C Optional transformation function for mu_C.
#' @param mu_transform_A Optional transformation function for mu_A.
#' @param outcome_fn A function defining the outcome to summarize (default: total cases).
#' @param o1 Baseline simulations for comparison.
#' @param ci_level Confidence level for credible intervals (default = 0.95).
#' @param out_dir Directory to save plots (optional).
#' @param month Logical. Whether to summarize and simulate using weekly data. Default is FALSE.
#' @param apply_decay Logical. Whether to apply decay to SMC coverage. Default is TRUE.
#' @param use_SMC_as_covariate Logical. Whether or not SMC is included in the model or as a covariate in the observation model.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{outputs}{A named list of lists with estimates, plots, and summaries for each scenario.}
#'   \item{summaries}{A named list of time series summaries for each scenario.}
#'   \item{estimates}{A named list of scalar outcome estimates for each scenario.}
#' }
#' @export
evaluate_multiple_scenarios <- function(patterns,
                                        smc_day_of_month = 1,
                                        model,
                                        param_inputs,
                                        param_samples,
                                        start_date,
                                        end_date,
                                        avg_cov,
                                        years,
                                        exclude_years = 2023,
                                        mu_transform_C = NULL,
                                        mu_transform_A = NULL,
                                        outcome_fn = function(y1, y0) sum(y1$inc_C_transformed),
                                        o1 = NULL,
                                        ci_level = 0.95,
                                        out_dir = NULL,
                                        month = FALSE,
                                        apply_decay = TRUE,
                                        use_SMC_as_covariate = FALSE,
                                        noise = FALSE) {

  outputs <- list()
  summaries <- list()
  estimates <- list()

  for (label in names(patterns)) {
    month_pattern <- patterns[[label]]
    n_years <- (lubridate::year(end_date) + 1) - lubridate::year(start_date)
    months_active <- matrix(rep(month_pattern, n_years + 1), nrow = n_years + 1, byrow = TRUE)

    smc_schedule <- gen_smc_schedule(start_date, end_date, years = years,
                                     months_active = months_active, coverage = avg_cov,
                                     smc_day_of_month = smc_day_of_month)

    if(use_SMC_as_covariate){
      if (apply_decay) {
        smc_schedule$cov <- smc_schedule$cov * smc_schedule$decay
      }

      if (month) {
        smc_summary <- calculate_monthly_metrics(smc_schedule, exclude_years = exclude_years)
      } else {
        smc_summary <- calculate_weekly_metrics(smc_schedule, exclude_years = exclude_years)
      }



      covariate_matrix <- data.frame(
        date_ymd = as.Date(smc_summary$date_ymd),
        cov_SMC = smc_summary$cov
      )

      # Custom daily rates
      r_df <- get_population_scaling(n = nrow(covariate_matrix), month = month,
                                     growth_rate_C = 1.000071,
                                     growth_rate_A = 1.000092)
      # Apply population growth
      covariate_matrix$r_C <- r_df$r_C

    }else{
      n_weeks <- length(seq(min(smc_schedule$dates), max(smc_schedule$dates), by = "week"))
      r_df <- get_population_scaling(n = n_weeks, month = month,
                                     growth_rate_C = 1.000071,
                                     growth_rate_A = 1.000092)

      covariate_matrix <- data.frame(date_ymd = obs_cases$date_ymd,
                                     r_C = r_df$r_C)


      param_inputs$SMC <- smc_schedule$SMC
      param_inputs$decay <- smc_schedule$decay
      param_inputs$cov_SMC <- smc_schedule$cov

      if (!apply_decay) {
        param_inputs$decay <- rep(1, length((param_inputs$decay)))
      }
    }

    sims <- run_simulations_from_samples(
      model = model,
      param_inputs = param_inputs,
      param_samples = param_samples,
      start_date = start_date,
      end_date = end_date,
      prewarm_years = 2,
      mu_transform_C = mu_transform_C,
      mu_transform_A = mu_transform_A,
      covariate_matrix = covariate_matrix,
      month = month,
      noise = noise
    )

    est <- calculate_estimate(o1 = o1, o2 = sims, outcome_fn = outcome_fn)
    estimates[[label]] <- est

    plot <- plot_estimate_distribution(est, x_label = "", title = label, ci_level = ci_level)
    if (!is.null(out_dir)) {
      save_plot_dynamic(plot, paste0("hist_", gsub(" ", "_", label)), out_dir)
    }

    sim_summary <- summarize_simulation_ci(sims, variables = "inc_C_transformed") %>%
      dplyr::mutate(scenario = label)

    outputs[[label]] <- list(
      estimate = est,
      plot = plot,
      summary = sim_summary
    )

    summaries[[label]] <- sim_summary
  }

  return(list(
    outputs = outputs,
    summaries = summaries,
    estimates = estimates
  ))
}
