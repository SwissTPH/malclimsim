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
#' @param patterns A named list of 12-element binary vectors or lists of such vectors (per year), indicating months when SMC is applied.
#' @param smc_day_of_month Integer specifying the day of the month to begin SMC each round (default is 1).
#' @param model A compiled model object to simulate from.
#' @param param_inputs A named list of baseline model parameters.
#' @param param_samples A matrix or data frame of sampled parameter sets (one row per sample).
#' @param start_date Simulation start date (character or Date).
#' @param end_date Simulation end date (character or Date).
#' @param avg_cov Average SMC coverage to apply during active months.
#' @param years A vector of years to apply SMC coverage (e.g., 2018:2023).
#' @param exclude_years Years to exclude when computing summaries (default is 2023).
#' @param mu_transform_C Optional transformation function for compartment C incidence.
#' @param mu_transform_A Optional transformation function for compartment A incidence.
#' @param outcome_fn Function computing the scalar outcome of interest from two simulation outputs (default: sum of transformed incidence).
#' @param o1 Optional baseline simulation for comparison.
#' @param ci_level Confidence level for summary intervals (default is 0.95).
#' @param out_dir Output directory to save plots. If NULL, plots are not saved.
#' @param month Logical. Whether the model and summaries are in monthly (TRUE) or weekly (FALSE) time (default is FALSE).
#' @param apply_decay Logical. Whether to apply time-decay to SMC coverage values (default is TRUE).
#' @param use_SMC_as_covariate Logical. Whether SMC is used as a covariate in the observation model instead of being built into the transmission model (default is FALSE).
#' @param noise Logical. Whether to add stochastic noise to simulated incidence (default is FALSE).
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
  cases_df1 <- list()
  cases_df2 <- list()

  for (label in names(patterns)) {
    month_pattern <- patterns[[label]]

    # --- NEW BLOCK: unify to an n_years x 12 matrix ---
    if (is.list(month_pattern) && all(names(month_pattern) %in% as.character(years))) {
      # ensure rows are in the same order as `years`
      year_strs     <- as.character(years)
      mat           <- do.call(rbind, month_pattern[year_strs])
      rownames(mat) <- year_strs
      months_active <- mat
    } else if (is.numeric(month_pattern) && length(month_pattern) == 12) {
      # legacy single-pattern: repeat it for each year
      months_active <- matrix(
        rep(month_pattern, length(years)),
        nrow   = length(years),
        byrow  = TRUE
      )
    } else {
      stop("Each element of `patterns` must be either a 12-vector or a named list-of-12-vectors.")
    }

    smc_schedule <- gen_smc_schedule(start_date, end_date, years = years, months_30_days = month,
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
      n_weeks <- length(seq(min(smc_schedule$dates), max(smc_schedule$dates), by = "week")) - 1
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

    if(!is.null(mu_transform_C)){
      sim_summary <- summarize_simulation_ci(sims, variables = "inc_C_transformed") %>%
        dplyr::mutate(scenario = label)
    }else{
      sim_summary <- summarize_simulation_ci(sims, variables = "inc_C") %>%
        dplyr::mutate(scenario = label)
    }


    outputs[[label]] <- list(
      estimate = est,
      plot = plot,
      summary = sim_summary
    )

    summaries[[label]] <- sim_summary
    cases_df1[[label]] <- o1
    cases_df2[[label]] <- sims
  }

  return(list(
    outputs = outputs,
    summaries = summaries,
    estimates = estimates,
    o1 = cases_df1,
    o2 = cases_df2
  ))
}

#' Plot Cases Averted Across SMC Scenarios (with 95% Credible Intervals)
#'
#' Generates a bar plot of estimated cases averted for each SMC pattern/scenario,
#' including 95% credible interval error bars. Supports either:
#' * A list of numeric vectors (posterior samples), or
#' * A list of lists with `mean`, `lower`, and `upper`.
#'
#' @param estimates A named list where each element is either:
#'   - A numeric vector of posterior draws, or
#'   - A list with elements `mean`, `lower`, and `upper`.
#' @param title Character string for the main plot title.
#' @param x_lab Character string for the x-axis label.
#' @param y_lab Character string for the y-axis label.
#' @param horizontal Logical; if TRUE, produces a horizontal bar plot with patterns on the y-axis.
#' @param by_year Logical; if TRUE, divides all values by `n_years` to yield per-year estimates.
#' @param n_years Numeric; the number of years used for scaling if `by_year = TRUE`. Must be positive.
#' @param bar_colors Named character vector specifying fill colors for each SMC pattern. The names should match the names in `estimates`.
#'
#' @return A `ggplot2` bar plot object showing mean cases averted with 95% CI.
#' @export
plot_cases_averted_barplot <- function(estimates,
                                       title      = "Estimated Cases Averted for Different SMC Timings",
                                       x_lab      = "SMC Timing",
                                       y_lab      = "Mean Cases Averted (95% CI)",
                                       horizontal = FALSE,
                                       by_year    = FALSE,
                                       n_years    = 1,
                                       bar_colors = c(
                                         "4 rounds (July start)" = "#0072B2",  # Blue
                                         "4 rounds (June start)" = "#E69F00",  # Orange
                                         "5 rounds (July start)" = "#009E73",  # Green
                                         "5 rounds (June start)" = "#D55E00"   # Red
                                       )) {

  # 1. Convert input into tidy data frame
  est_df <- purrr::imap_dfr(estimates, function(est, label) {
    if (is.numeric(est) && !is.list(est)) {
      quant <- quantile(est, probs = c(0.025, 0.975), na.rm = TRUE)
      tibble::tibble(
        pattern       = label,
        cases_averted = mean(est, na.rm = TRUE),
        lower         = quant[1],
        upper         = quant[2]
      )
    } else if (is.list(est) && all(c("mean", "lower", "upper") %in% names(est))) {
      tibble::tibble(
        pattern       = label,
        cases_averted = est$mean,
        lower         = est$lower,
        upper         = est$upper
      )
    } else {
      stop("Each element of 'estimates' must be a numeric vector or a list with 'mean', 'lower', and 'upper'.")
    }
  })

  # 2. Optionally scale to per-year values
  if (by_year) {
    if (!is.numeric(n_years) || n_years <= 0) {
      stop("`n_years` must be a positive number when `by_year = TRUE`.")
    }
    est_df <- est_df %>%
      mutate(
        cases_averted = cases_averted / n_years,
        lower         = lower         / n_years,
        upper         = upper         / n_years
      )
    if (y_lab == "Mean Cases Averted (95% CI)") {
      y_lab <- "Mean Yearly Cases Averted (95% CI)"
    }
  }

  # 3. Sort and add colors
  est_df <- est_df %>%
    #arrange(desc(cases_averted)) %>%
    mutate(
      pattern = factor(pattern, levels = unique(pattern)),
      fill_color = bar_colors[as.character(pattern)]
    )

  # 4. Plot
  p <- ggplot(est_df, aes(x = pattern, y = cases_averted, fill = pattern)) +
    geom_col(alpha = 0.9) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.2, color = "black"
    ) +
    scale_fill_manual(values = bar_colors, drop = FALSE) +
    labs(
      title = title,
      x     = x_lab,
      y     = y_lab,
      fill  = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text  = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.position = "none"
    )

  # 5. Optional horizontal layout
  if (horizontal) {
    p <- p +
      coord_flip() +
      theme(axis.text.y = element_text(hjust = 1))
  } else {
    p <- p +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(p)
}
