#' Create Observation Function Configuration
#'
#' Constructs a configuration list specifying how to compare observed and modeled data.
#' This configuration is passed to the observation likelihood generator function
#' (`generate_incidence_comparison()`) and determines which outcomes are compared and how.
#'
#' @param use_monthly Logical. If `TRUE`, observation and model outputs are matched at monthly resolution;
#' otherwise, weekly data is assumed.
#' @param age_group Character. One of `"u5"`, `"o5"`, `"sep_ages"`, `"total"`, or `"all_ages"`.
#' This determines which age group(s) are included in the likelihood.
#' @param include_prev Logical. If `TRUE`, the observation function includes a prevalence likelihood term.
#' @param use_SMC_as_covariate Logical. If `TRUE`, the incidence is adjusted using a regression effect (`beta_1`) on SMC coverage.
#' @param log_link Logical. If `TRUE`, a log-link is used to transform the predicted incidence when adjusting by SMC coverage.
#' If `FALSE`, a linear form is used.
#' @param include_pop_growth Logical. If `TRUE`, model-predicted incidence is rescaled to account for
#' population growth using the `r_C` and/or `r_A` outputs from the model.
#'
#' @return A named list containing configuration options to pass to the observation likelihood function.
#'
#' @examples
#' obs_config <- make_obs_config(
#'   use_monthly = TRUE,
#'   age_group = "u5",
#'   include_prev = FALSE,
#'   use_SMC_as_covariate = TRUE,
#'   log_link = TRUE,
#'   include_pop_growth = TRUE
#' )
#'
#' @export
make_obs_config <- function(use_monthly = TRUE,
                            age_group = "u5",
                            include_prev = FALSE,
                            use_SMC_as_covariate = FALSE,
                            log_link = FALSE,
                            include_pop_growth = FALSE) {
  list(
    time = ifelse(use_monthly, "month", "week"),
    age_group = age_group,
    include_prev = include_prev,
    use_SMC_as_covariate = use_SMC_as_covariate,
    log_link = log_link,
    include_pop_growth = include_pop_growth
  )
}

#' Get Model Output Field Names for Observation Function
#'
#' Returns a named list of model state variables used for incidence comparisons,
#' based on the time scale and age group configuration.
#'
#' @param time Character string. Either "month" or "week".
#' @param age_group Character. One of: "total", "sep_ages", "u5", "o5", "all_ages".
#'
#' @return A named list with one or more of: mu, mu_C, mu_A, mu_C1, mu_C2, etc.
#' These correspond to expected incidence variables in the model output.
#' @export
get_field_mapping <- function(time, age_group) {
  prefix <- switch(time,
                   month = "month_inc",
                   week  = "wk_inc",
                   stop("Invalid time unit: must be 'month' or 'week'")
  )

  # Return mapping based on age group
  switch(age_group,
         total = list(
           mu = paste0(prefix, "_total")
         ),

         sep_ages = list(
           mu_C = paste0(prefix, "_C"),   # Children
           mu_A = paste0(prefix, "_A")    # Adults
         ),

         all_ages = list(
           mu_C1 = paste0(prefix, "_C1"), # Younger children
           mu_C2 = paste0(prefix, "_C2"), # Older children
           mu_A  = paste0(prefix, "_A")   # Adults
         ),

         u5 = list(
           mu_C = paste0(prefix, "_C")    # Under-5 children
         ),

         o5 = list(
           mu_A = paste0(prefix, "_A")    # Over-5 population
         ),

         stop("Invalid age group: must be one of 'total', 'sep_ages', 'u5', 'o5', or 'all_ages'")
  )
}

generate_incidence_comparison <- function(month,
                                          age_for_inf,
                                          incidence_df,
                                          include_prev,
                                          use_SMC_as_covariate = FALSE,
                                          log_link = TRUE,
                                          include_pop_growth = FALSE,
                                          r_C = NULL,
                                          r_A = NULL) {
  time <- ifelse(month, "month", "week")
  mapping <- get_field_mapping(time, age_for_inf)
  coverage_divisor <- if (month) 30 else 7

  # --- Helper: Unadjusted likelihood ---
  ll_inc_nb <- function(observed, predicted, size, obs_weight = 1, r_scale = 1) {
    if (is.na(observed)) return(0)
    #cat(paste0(round(predicted, 0), ","))
    dnbinom(x = observed, mu = predicted, size = size, log = TRUE)
  }

  ll_inc_nb <- function(observed, predicted, size, obs_weight = 1) {
    if (is.na(observed)) return(0)
    obs_weight * dnbinom(x = observed, mu = predicted, size = size, log = TRUE)
  }

  # --- Helper: SMC-adjusted incidence likelihood ---
  if (log_link) {
    ll_inc_nb_beta_adjusted <- function(inc_observed, inc_predicted, beta, coverage,
                                        size, treatment, obs_weight = 1, divisor = 30, r_scale = 1) {
      if (is.na(inc_observed)) return(0)
      coverage <- ifelse(is.null(coverage), 0, coverage / divisor)
      inc_predicted <- pmax(inc_predicted, 1e-6)
      mu_adjusted <- exp(log(inc_predicted) + beta * coverage) * r_scale
      size <- ifelse(size == 0, 1e-3, size)
      obs_weight * dnbinom(x = inc_observed, mu = mu_adjusted, size = size, log = TRUE)
    }
  } else {
    ll_inc_nb_beta_adjusted <- function(inc_observed, inc_predicted, beta, coverage,
                                        size, treatment, obs_weight = 1, divisor = 30, r_scale = 1) {
      if (is.na(inc_observed)) return(0)
      coverage <- ifelse(is.null(coverage), 0, coverage / divisor)
      inc_predicted <- pmax(inc_predicted, 1e-6)
      mu_adjusted <- inc_predicted * (1 - beta * coverage) * r_scale
      size <- ifelse(size == 0, 1e-3, size)
      obs_weight * dnbinom(x = inc_observed, mu = mu_adjusted, size = size, log = TRUE)
    }
  }

  # --- Helper: Prevalence likelihood ---
  ll_prev_beta <- function(observed, predicted, kappa) {
    obs_clamped <- pmin(pmax(observed, 1e-6), 1 - 1e-6)
    alpha <- predicted * kappa
    beta <- (1 - predicted) * kappa
    if (is.na(obs_clamped)) return(0)
    dbeta(obs_clamped, shape1 = alpha, shape2 = beta, log = TRUE)
  }

  # --- Main log-likelihood function ---
  return(function(state, observed, pars) {
    total_ll <- 0

    # --- Children ---
    if (!is.null(mapping$mu_C)) {
      mu_C <- state[mapping$mu_C, , drop = TRUE]
      r_C <- if (include_pop_growth) state["r_C", , drop = TRUE] else 1
      size_1 <- ifelse(is.null(pars$size_1), 1e-3, max(pars$size_1, 1e-3))

      if (use_SMC_as_covariate && !is.null(pars$beta_1)) {
        beta <- if (log_link) pars$beta_1 else pars$beta_2
        ll_C <- ll_inc_nb_beta_adjusted(
          inc_observed = observed$inc_C,
          inc_predicted = mu_C,
          beta = beta,
          coverage = observed$cov_SMC,
          size = size_1,
          treatment = observed$treatment,
          obs_weight = observed$obs_weight,
          divisor = coverage_divisor,
          r_scale = r_C
        )
      } else {
        mu_C_adj <- mu_C * r_C
        ll_C <- ll_inc_nb(observed = observed$inc_C,
                          predicted = mu_C_adj,
                          size = size_1,
                          obs_weight = observed$obs_weight)
      }

      total_ll <- total_ll + ll_C
    }

    # --- Adults ---
    if (!is.null(mapping$mu_A)) {
      mu_A <- state[mapping$mu_A, , drop = TRUE]
      r_A <- if (include_pop_growth) state["r_A", , drop = TRUE] else 1
      mu_A_adj <- mu_A * r_A
      size_2 <- ifelse(is.null(pars$size_2), 1e-3, max(pars$size_2, 1e-3))
      ll_A <- ll_inc_nb(observed = observed$inc_A,
                        predicted = mu_A_adj,
                        size = size_2,
                        obs_weight = observed$obs_weight)
      total_ll <- total_ll + ll_A
    }

    return(total_ll)
  })
}



#' Print Configuration for Observation Function
#'
#' Utility to confirm current likelihood settings (time unit, age group, etc.).
print_observation_structure <- function(month, age_for_inf, include_prev, use_beta_adjustment = FALSE) {
  cat("Observation Function Configuration:\n")
  cat("- Time unit:", ifelse(month, "Monthly", "Weekly"), "\n")
  cat("- Age group:", age_for_inf, "\n")
  cat("- Includes prevalence:", include_prev, "\n")
  cat("- Uses beta_1 adjustment:", use_beta_adjustment, "\n")
}

#' Get Observation Function Based on User Configuration
#'
#' This helper function allows users to specify key options such as time resolution,
#' age group, inclusion of prevalence, and whether to model SMC as a covariate using `beta_1`.
#'
#' @param month Logical. TRUE for monthly data, FALSE for weekly.
#' @param age_for_inf Character. One of: "total", "sep_ages", "u5", "o5", "all_ages".
#' @param incidence_df A data frame of observed incidence and prevalence values.
#' @param include_prev Logical. If TRUE, includes prevalence likelihood.
#' @param use_SMC_as_covariate Logical. If TRUE, uses `beta_1` covariate effect on incidence.
#'
#' @return A likelihood function that takes (state, observed, pars) as arguments.
#' @export
get_observation_function <- function(month = TRUE,
                                     age_for_inf = "u5",
                                     incidence_df,
                                     include_prev = TRUE,
                                     use_SMC_as_covariate = FALSE) {

  # Print configuration for clarity
  cat("Configuring Observation Function:\n")
  cat("- Time resolution:", ifelse(month, "Monthly", "Weekly"), "\n")
  cat("- Age group:", age_for_inf, "\n")
  cat("- Include prevalence in likelihood:", include_prev, "\n")
  cat("- Use SMC as covariate (via beta_1):", use_SMC_as_covariate, "\n\n")

  # Create and return the observation function
  obs_fn <- generate_incidence_comparison(
    month = month,
    age_for_inf = age_for_inf,
    incidence_df = incidence_df,
    include_prev = include_prev,
    use_SMC_as_covariate = use_SMC_as_covariate,
    log_link = TRUE
  )

  return(obs_fn)
}
