#' Generate a Comparison Function for Malaria Incidence
#'
#' This function generates a likelihood comparison function based on the specified
#' time period (monthly or weekly) and age group for malaria incidence.
#'
#' @param month A logical value indicating whether to use monthly data (TRUE) or weekly data (FALSE).
#' @param age_for_inf A string specifying the age group for incidence:
#'        "total", "sep_ages", "all_ages", "u5", or "o5".
#' @param incidence_df A data frame containing observed malaria incidence data.
#' @param include_prev A logical value indicating if prevalence should be included in the likelihood function
#'
#' @return A function for comparing modeled and observed incidence based on a specified distribution (e.g., Negative Binomial).
#' @export
#'
#' @examples
#' # Example usage for monthly total incidence comparison:
#' comparison_fn <- generate_incidence_comparison(month = TRUE, age_for_inf = "total",
#' incidence_df = incidence_data, include_prev = include_prev)
#' # Then use the returned function in MCMC or other model fitting tasks
generate_incidence_comparison <- function(month, age_for_inf, incidence_df, include_prev) {
  # Helper function to calculate rolling average
  calculate_rolling_average <- function(rainfall, D) {
    D <- round(D)  # Ensure D is an integer
    rolling_avg <- zoo::rollapply(rainfall, width = D, FUN = mean, align = "right", fill = NA)
    # Handle edge cases
    rolling_avg[is.na(rolling_avg)] <- mean(rainfall[1:D], na.rm = TRUE)
    return(rolling_avg)
  }

  # Set up the comparison function based on conditions
  if(month && age_for_inf == 'total') {
    return(function(state, observed, pars = c("size_1")) {
      incidence_observed_total <- observed$inc
      mu_total <- state["month_inc_total", , drop = TRUE]  # Cumulative monthly incidence
      size_1 <- pars$size_1
      if (size_1 == 0) { size_1 <- 1e30 }
      return(dnbinom(x = incidence_observed_total, mu = mu_total, size = size_1, log = TRUE))
    })
  }


  # if(month && age_for_inf == 'sep_ages' && include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1", "size_2")) {
  #     # Incidence Log-Likelihood
  #     incidence_observed_C <- observed$inc_C
  #     incidence_observed_A <- observed$inc_A
  #     mu_C <- state["month_inc_C", , drop = TRUE]
  #     mu_A <- state["month_inc_A", , drop = TRUE]
  #     size_1 <- pars$size_1
  #     size_2 <- pars$size_2
  #     if (size_1 == 0) { size_1 <- 1e5 }
  #     if (size_2 == 0) { size_2 <- 1e5 }
  #     if(is.na(incidence_observed_C)){ll_C <- 0}else{ # if no observed value, does not contribute to likelihood
  #       ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE)
  #     }
  #     if(is.na(incidence_observed_A)){ll_A <- 0}else{
  #       ll_A <- dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE)
  #     }
  #     #print(paste("size_1:", round(size_1, 2), "size_2:",
  #     #            round(size_2, 2), "mu_C:", round(mu_C, 2), "mu_A:",
  #     #            round(mu_A, 2), "obs_C:", incidence_observed_C,
  #     #            "obs_A:", incidence_observed_A))
  #     return(ll_C + ll_A)
  #   })
  # }

  if (month && age_for_inf == 'sep_ages' && include_prev == TRUE) {
    return(function(state, observed, pars = c("size_1", "size_2", "kappa_C", "kappa_A")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["month_inc_C", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
      size_2 <- ifelse(pars$size_2 == 0, 1e-3, pars$size_2)

      ll_C <- ifelse(is.na(incidence_observed_C), 0,
                     dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))
      ll_A <- ifelse(is.na(incidence_observed_A), 0,
                     dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))

      # Prevalence Likelihood for Children
      prev_observed_C <- pmin(pmax(observed$prev_C, 1e-6), 1 - 1e-6)
      prev_model_C <- state["prev_C_1", , drop = TRUE]
      mu_prev_C <- prev_model_C
      kappa_C <- pars$kappa_C
      alpha_C <- mu_prev_C * kappa_C
      beta_C <- (1 - mu_prev_C) * kappa_C
      ll_prev_C <- ifelse(is.na(prev_observed_C), 0,
                          dbeta(prev_observed_C, shape1 = alpha_C, shape2 = beta_C, log = TRUE))

      # Prevalence Likelihood for Adults
      prev_observed_A <- pmin(pmax(observed$prev_A, 1e-6), 1 - 1e-6)
      prev_model_A <- state["prev_A_1", , drop = TRUE]
      mu_prev_A <- prev_model_A
      kappa_A <- pars$kappa_A
      alpha_A <- mu_prev_A * kappa_A
      beta_A <- (1 - mu_prev_A) * kappa_A
      ll_prev_A <- ifelse(is.na(prev_observed_A), 0,
                          dbeta(prev_observed_A, shape1 = alpha_A, shape2 = beta_A, log = TRUE))

      # Weight for prevalence likelihood
      prev_weight <- 12  # Adjust based on data frequency

      # Combine all likelihood components
      total_ll <- ll_C + ll_A + prev_weight * (ll_prev_C + ll_prev_A)

      # Diagnostics
      #cat("ll_C:", ll_C, " ll_A:", ll_A, " ll_prev_C:", ll_prev_C, " ll_prev_A:", ll_prev_A, " Total:", total_ll, "\n")

      return(total_ll)
    })
  }

  if (month && age_for_inf == 'sep_ages' && include_prev == FALSE) {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["month_inc_C", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
      size_2 <- ifelse(pars$size_2 == 0, 1e-3, pars$size_2)

      ll_C <- ifelse(is.na(incidence_observed_C), 0,
                     dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))
      ll_A <- ifelse(is.na(incidence_observed_A), 0,
                     dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))

      # Predicted prevalence from model
      prev_model_C <- state["prev_C_1", , drop = TRUE]
      prev_model_A <- state["prev_A_1", , drop = TRUE]

      # Apply penalty if prev_C is smaller than prev_A
      penalty <- ifelse(prev_model_C < prev_model_A, -1e6 * abs(prev_model_C - prev_model_A), 0)

      # Combine all likelihood components
      #total_ll <- ll_C + ll_A + penalty

      total_ll <- ifelse(observed$Treatment == 1, ll_C + ll_A + penalty,
             12 * (ll_C + ll_A + penalty))
      return(total_ll)
    })
  }


  if(month && age_for_inf == 'all_ages') {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      incidence_observed_C1 <- observed$inc_C1
      incidence_observed_C2 <- observed$inc_C2
      incidence_observed_A <- observed$inc_A
      mu_C1 <- state["month_inc_C1", , drop = TRUE]
      mu_C2 <- state["month_inc_C2", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_1 <- pars$size_1
      size_2 <- pars$size_2
      if (size_1 == 0) { size_1 <- 1e30 }
      if (size_2 == 0) { size_2 <- 1e30 }
      return(dnbinom(x = incidence_observed_C1, mu = mu_C1, size = size_1, log = TRUE) +
               dnbinom(x = incidence_observed_C2, mu = mu_C2, size = size_1, log = TRUE) +
               dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))
    })
  }

  # if(month && age_for_inf == 'u5') {
  #   return(function(state, observed, pars = c("size_1")) {
  #     if (is.na(observed$inc_C)) {
  #       return(NULL)
  #     }
  #     incidence_observed_C <- observed$inc_C # this is the observed data
  #     mu_C <- state["month_inc_C", , drop = TRUE] # this is "x"
  #     size_1 <- pars$size_1
  #     if (size_1 == 0) { size_1 <- 1e30 }
  #     return(dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))
  #   })
  # }

  # if(month && age_for_inf == 'u5') {
  #   return(function(state, observed, pars = c("size")) {
  #     incidence_observed_C <- observed$inc_C # this is the observed data
  #     mu_C <- state["month_inc_C", , drop = TRUE] # this is "x"
  #     size <- pars$size
  #     if (size == 0) { size <- 1e30 }
  #     ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size, log = TRUE)
  #     if(is.na(ll_C)){ll_C <- 0}
  #     return(ll_C)
  #   })
  # }

  if(month && age_for_inf == 'u5' & include_prev == TRUE) {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      mu_C <- state["month_inc_C", , drop = TRUE]
      size_1 <- pars$size_1
      if (size_1 == 0) { size_1 <- 1e5 }
      if(is.na(incidence_observed_C)){ll_C <- 0}else{ # if no observed value, does not contribute to likelihood
        ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE)
      }
      # Prevalence Penalty
      prev_observed_C <- observed$prev_C
      prev_model_C <- state["prev_C_2", , drop = TRUE]
      if (is.na(prev_observed_C)){penalty_prev_C <- 0}else{ # if no observed value, does not contribute to penalty
        #penalty_prev_C <- (prev_observed_C - prev_model_C)^2
        if(abs(prev_observed_C - prev_model_C) > 0.1) {
          penalty_prev_C <- 10000  # Large penalty for being outside range
        } else {
          penalty_prev_C <- 0      # No penalty if within range
        }

      }
      penalty_prev <- penalty_prev_C
      # **Determine an adaptive penalty weight**
      #lambda <- mean(abs(c(ll_C, ll_A))) / penalty_prev  # Avoid division by zero
      return(ll_C - penalty_prev)  # Use weighted penalty
    })
  }

  # if(month && age_for_inf == 'u5' & include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1", "size_2")) {
  #     # Incidence Log-Likelihood
  #     incidence_observed_C <- observed$inc_C
  #     mu_C <- state["month_inc_C", , drop = TRUE]
  #     size_1 <- pars$size_1
  #     if (size_1 == 0) { size_1 <- 1e5 }
  #     if(is.na(incidence_observed_C)){ll_C <- 0}else{ # if no observed value, does not contribute to likelihood
  #       ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE)
  #     }
  #     return(ll_C)  # Use weighted penalty
  #   })
  # }

  # if (month && age_for_inf == 'u5' && include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1")) {
  #     # Incidence Log-Likelihood
  #     incidence_observed_C <- observed$inc_C
  #     mu_C <- state["month_inc_C", , drop = TRUE]
  #     size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
  #
  #     ll_C <- ifelse(is.na(incidence_observed_C), 0,
  #                    dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))
  #
  #     return(ll_C)
  #   })
  # }

  # if (month && age_for_inf == 'u5' && include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1", "size_2")) {
  #
  #     # Incidence Log-Likelihood
  #     incidence_observed_C <- observed$inc_C
  #     mu_C <- max(state["month_inc_C", , drop = TRUE], 1e-6)
  #     size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
  #
  #     #print(mu_C)
  #
  #     ll_C <- ifelse(is.na(incidence_observed_C), 0,
  #                    dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))
  #
  #     # Predicted prevalence from model
  #     prev_model_C <- state["prev_C_1", , drop = TRUE]
  #
  #     # Combine all likelihood components
  #     #total_ll <- ll_C + ll_A + penalty
  #
  #     #total_ll <- ifelse(observed$Treatment == 1, ll_C,
  #     #                   12 * ll_C)
  #     total_ll <- ll_C
  #     return(total_ll)
  #   })
  # }

  # if (month && age_for_inf == 'u5' && include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1", "size_2", "beta_1")) {
  #
  #     # Retrieve observed incidence
  #     incidence_observed_C <- observed$inc_C
  #
  #     # Retrieve model-predicted incidence
  #     mu_C <- max(state["month_inc_C", , drop = TRUE], 1e-6)
  #
  #     # Retrieve SMC coverage for the current month
  #     coverage <- observed$cov_SMC / 30
  #
  #     # Apply the multiplicative reductive effect of SMC
  #     mu_C_adjusted <- mu_C * (1 - pars$beta_1 * coverage)
  #
  #     # Ensure size parameter is positive
  #     size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
  #
  #     # Calculate the negative binomial log-likelihood
  #     ll_C <- ifelse(is.na(incidence_observed_C), 0,
  #                    dnbinom(x = incidence_observed_C, mu = mu_C_adjusted, size = size_1, log = TRUE))
  #     if(observed$Treatment == 0){ll_C = 12 * ll_C}
  #
  #     # Return the total log-likelihood
  #     return(ll_C)
  #   })
  # }

  if (month && age_for_inf == 'u5' && include_prev == FALSE) {
    return(function(state, observed, pars = c("size_1", "size_2", "beta_1")) {

      # Retrieve observed incidence
      incidence_observed_C <- observed$inc_C

      # Retrieve model-predicted incidence (latent)
      mu_C <- max(state["month_inc_C", , drop = TRUE], 1e-6)  # Ensure positivity


      # Ensure size parameter is positive (to avoid numerical issues)
      size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)

      # Compute the NB likelihood using the adjusted mean
      ll_C <- ifelse(is.na(incidence_observed_C), 0,
                     dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))

      if(observed$Treatment == 0){ll_C = 1 * ll_C}

      return(ll_C)
    })
  }

  # if (month && age_for_inf == 'u5' && include_prev == FALSE) {
  #   return(function(state, observed, pars = c("size_1", "size_2", "beta_1")) {
  #
  #     # Retrieve observed incidence
  #     incidence_observed_C <- observed$inc_C
  #
  #     # Retrieve model-predicted incidence (latent)
  #     mu_C <- max(state["month_inc_C", , drop = TRUE], 1e-6)  # Ensure positivity
  #
  #     # Retrieve SMC coverage for the current month
  #     coverage <- ifelse(is.null(observed$cov_SMC), 0, observed$cov_SMC / 30)
  #
  #     # Apply log link transformation to ensure multiplicative effect
  #     log_mu_C_adjusted <- log(mu_C) + pars$beta_1 * coverage
  #     mu_C_adjusted <- exp(log_mu_C_adjusted)  # Ensures incidence is always positive
  #
  #     # Ensure size parameter is positive (to avoid numerical issues)
  #     size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)
  #
  #     # Compute the NB likelihood using the adjusted mean
  #     ll_C <- ifelse(is.na(incidence_observed_C), 0,
  #                    dnbinom(x = incidence_observed_C, mu = mu_C_adjusted, size = size_1, log = TRUE))
  #
  #     if(observed$Treatment == 0){ll_C = 8 * ll_C}
  #
  #     return(ll_C)
  #   })
  # }


  if (month && age_for_inf == 'u5' && include_prev == TRUE) {
    return(function(state, observed, pars = c("size_1", "kappa_C")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      mu_C <- state["month_inc_C", , drop = TRUE]
      size_1 <- ifelse(pars$size_1 == 0, 1e-3, pars$size_1)

      ll_C <- ifelse(is.na(incidence_observed_C), 0,
                     dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))

      # Prevalence Likelihood (Beta Distribution)
      prev_observed_C <- pmin(pmax(observed$prev_C, 1e-6), 1 - 1e-6)
      prev_model_C <- state["prev_C_2", , drop = TRUE]
      mu_prev_C <- prev_model_C
      kappa_C <- pars$kappa_C
      alpha_C <- mu_prev_C * kappa_C
      beta_C <- (1 - mu_prev_C) * kappa_C

      ll_prev_C <- ifelse(is.na(prev_observed_C), 0,
                          dbeta(prev_observed_C, shape1 = alpha_C, shape2 = beta_C, log = TRUE))

      # Prevalence weight (similar to sep_ages)
      prev_weight <- 12  # Adjust as necessary

      # Combine Likelihoods
      total_ll <- ll_C + prev_weight * ll_prev_C

      # Diagnostics for Debugging
      #cat("ll_C:", ll_C, " ll_prev_C:", ll_prev_C, " Total:", total_ll, "\n")

      return(total_ll)
    })
  }



  # if(month && age_for_inf == 'u5') {
  #   return(function(state, observed, pars) {
  #     incidence_observed_C <- observed$inc_C # this is the observed data
  #     mu_C <- state["month_inc_C", , drop = TRUE] # this is "x"
  #     return(dpois(x = incidence_observed_C, lambda = mu_C, log = TRUE))
  #   })
  # }

  if (month && age_for_inf == 'o5' && include_prev == FALSE) {
    return(function(state, observed, pars = c("size_2")) {
      # Incidence Log-Likelihood (No Prevalence)
      incidence_observed_A <- observed$inc_A
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_2 <- ifelse(pars$size_2 == 0, 1e-3, pars$size_2)

      ll_A <- ifelse(is.na(incidence_observed_A), 0,
                     dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))

      return(ll_A)
    })
  }

  if (month && age_for_inf == 'o5' && include_prev == TRUE) {
    return(function(state, observed, pars = c("size_2", "kappa_A")) {
      # Incidence Log-Likelihood
      incidence_observed_A <- observed$inc_A
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_2 <- ifelse(pars$size_2 == 0, 1e-3, pars$size_2)

      ll_A <- ifelse(is.na(incidence_observed_A), 0,
                     dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))

      # Prevalence Likelihood (Beta Distribution)
      prev_observed_A <- pmin(pmax(observed$prev_A, 1e-6), 1 - 1e-6)
      prev_model_A <- state["prev_A_1", , drop = TRUE]
      mu_prev_A <- prev_model_A
      kappa_A <- pars$kappa_A
      alpha_A <- mu_prev_A * kappa_A
      beta_A <- (1 - mu_prev_A) * kappa_A

      ll_prev_A <- ifelse(is.na(prev_observed_A), 0,
                          dbeta(prev_observed_A, shape1 = alpha_A, shape2 = beta_A, log = TRUE))

      # Prevalence Weight (align with sep_ages)
      prev_weight <- 12  # Adjust as needed

      # Combine Likelihoods
      total_ll <- ll_A + prev_weight * ll_prev_A

      # Diagnostics for Debugging
      #cat("ll_A:", ll_A, " ll_prev_A:", ll_prev_A, " Total:", total_ll, "\n")

      return(total_ll)
    })
  }


  # For weekly data
  if(!month && age_for_inf == 'total') {
    return(function(state, observed, pars = c("size_1")) {
      incidence_observed_total <- observed$inc
      mu_total <- state["wk_inc_total", , drop = TRUE]
      size_1 <- pars$size_1
      if (size_1 == 0) { size_1 <- 1e30 }
      return(dnbinom(x = incidence_observed_total, mu = mu_total, size = size_1, log = TRUE))
    })
  }

  if(!month && age_for_inf == 'sep_ages') {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["wk_inc_C", , drop = TRUE]
      mu_A <- state["wk_inc_A", , drop = TRUE]
      size_1 <- pars$size_1
      size_2 <- pars$size_2
      if (size_1 == 0) { size_1 <- 1e30 }
      if (size_2 == 0) { size_2 <- 1e30 }
      return(dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE) +
               dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))
    })
  }

  # If no condition is met, return NULL or a default message
  stop("Invalid combination of 'month' and 'age_for_inf'.")
}
