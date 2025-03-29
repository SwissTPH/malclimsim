#' Create Observation Function Configuration
#'
#' Constructs a configuration list specifying how to compare observed and modeled data.
#'
#' @param use_monthly Logical. TRUE for monthly data, FALSE for weekly.
#' @param age_group Character. One of "u5", "o5", "sep_ages", "total", "all_ages".
#' @param include_prev Logical. Include prevalence likelihood.
#' @param use_SMC_as_covariate Logical. Use beta_1 as covariate for SMC coverage.
#'
#' @return A named list suitable for passing to comparison functions.
#' @export
make_obs_config <- function(use_monthly = TRUE,
                            age_group = "u5",
                            include_prev = TRUE,
                            use_SMC_as_covariate = FALSE) {
  list(
    time = ifelse(use_monthly, "month", "week"),
    age_group = age_group,
    include_prev = include_prev,
    use_SMC_as_covariate = use_SMC_as_covariate
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


#' Generalized Likelihood Generator for Malaria Incidence & Prevalence
#'
#' Builds a log-likelihood function tailored to the configuration specified.
#' Includes internal log-likelihood helpers for safe use with mcstate/callr.
#'
#' @param month Logical. Use monthly (TRUE) or weekly (FALSE) data.
#' @param age_for_inf Age group string: "u5", "o5", "total", etc.
#' @param incidence_df Data frame of observed incidence/prevalence.
#' @param include_prev Logical. If TRUE, include prevalence components.
#' @param use_SMC_as_covariate Logical. If TRUE, apply beta_1 to adjust incidence by SMC coverage.
#'
#' @return A function: function(state, observed, pars) -> log-likelihood.
generate_incidence_comparison <- function(month,
                                          age_for_inf,
                                          incidence_df,
                                          include_prev,
                                          use_SMC_as_covariate = FALSE) {

  time <- ifelse(month, "month", "week")
  mapping <- get_field_mapping(time, age_for_inf)

  # --- Inlined helper functions ---

  ll_inc_nb <- function(observed, predicted, size) {
    if (is.na(observed)) return(0)
    dnbinom(x = observed, mu = predicted, size = size, log = TRUE)
  }

  ll_inc_nb_beta_adjusted <- function(inc_observed, inc_predicted, beta_1, coverage,
                                      size, treatment, obs_weight = 1) {
    if (is.na(inc_observed)) return(0)
    coverage <- ifelse(is.null(coverage), 0, coverage / 30)
    inc_predicted <- pmax(inc_predicted, 1e-6)
    mu_adjusted <- exp(log(inc_predicted) + beta_1 * coverage)
    size <- ifelse(size == 0, 1e-3, size)
    obs_weight * dnbinom(x = inc_observed, mu = mu_adjusted, size = size, log = TRUE)
  }

  ll_prev_beta <- function(observed, predicted, kappa) {
    obs_clamped <- pmin(pmax(observed, 1e-6), 1 - 1e-6)
    alpha <- predicted * kappa
    beta <- (1 - predicted) * kappa
    if (is.na(obs_clamped)) return(0)
    dbeta(obs_clamped, shape1 = alpha, shape2 = beta, log = TRUE)
  }

  # --- Likelihood function returned ---
  return(function(state, observed, pars) {
    total_ll <- 0  # Initialize total log-likelihood

    # --- Incidence: Children ---
    if (!is.null(mapping$mu_C)) {
      mu_C <- state[mapping$mu_C, , drop = TRUE]
      size_1 <- ifelse(is.null(pars$size_1), 1e-3, max(pars$size_1, 1e-3))

      if (use_SMC_as_covariate && !is.null(pars$beta_1)) {
        cov <- observed$cov_SMC
        ll_C <- ll_inc_nb_beta_adjusted(inc_observed = observed$inc_C, inc_predicted = mu_C,
                                        beta_1 = pars$beta_1, coverage = cov,
                                        size = size_1, treatment = observed$treatment,
                                        obs_weight = observed$obs_weight)
      } else {
        ll_C <- ll_inc_nb(observed$inc_C, mu_C, size_1)
      }

      total_ll <- total_ll + ll_C
    }

    # --- Incidence: Adults ---
    if (!is.null(mapping$mu_A)) {
      mu_A <- state[mapping$mu_A, , drop = TRUE]
      size_2 <- ifelse(is.null(pars$size_2), 1e-3, max(pars$size_2, 1e-3))
      ll_A <- ll_inc_nb(observed$inc_A, mu_A, size_2)
      total_ll <- total_ll + ll_A
    }

    # --- Prevalence ---
    if (include_prev) {
      if (!is.null(observed$prev_C) && !is.null(pars$kappa_C)) {
        mu_prev_C <- state["prev_C_1", , drop = TRUE]
        ll_prev_C <- ll_prev_beta(observed$prev_C, mu_prev_C, pars$kappa_C)
        total_ll <- total_ll + 12 * ll_prev_C
      }

      if (!is.null(observed$prev_A) && !is.null(pars$kappa_A)) {
        mu_prev_A <- state["prev_A_1", , drop = TRUE]
        ll_prev_A <- ll_prev_beta(observed$prev_A, mu_prev_A, pars$kappa_A)
        total_ll <- total_ll + 12 * ll_prev_A
      }
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
    use_SMC_as_covariate = use_SMC_as_covariate
  )

  return(obs_fn)
}
