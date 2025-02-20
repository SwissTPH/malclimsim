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


  if(month && age_for_inf == 'sep_ages' && include_prev == FALSE) {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["month_inc_C", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_1 <- pars$size_1
      size_2 <- pars$size_2
      if (size_1 == 0) { size_1 <- 1e5 }
      if (size_2 == 0) { size_2 <- 1e5 }
      if(is.na(incidence_observed_C)){ll_C <- 0}else{ # if no observed value, does not contribute to likelihood
        ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE)
      }
      if(is.na(incidence_observed_A)){ll_A <- 0}else{
        ll_A <- dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE)
      }
      #print(paste("size_1:", round(size_1, 2), "size_2:",
      #            round(size_2, 2), "mu_C:", round(mu_C, 2), "mu_A:",
      #            round(mu_A, 2), "obs_C:", incidence_observed_C,
      #            "obs_A:", incidence_observed_A))
      return(ll_C + ll_A)
    })
  }

  if (month && age_for_inf == 'sep_ages' & include_prev == TRUE) {
    return(function(state, observed, pars = c("size_1", "size_2", "kappa_C", "kappa_A")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["month_inc_C", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_1 <- pars$size_1
      size_2 <- pars$size_2
      if (size_1 == 0) { size_1 <- 1e30 }
      if (size_2 == 0) { size_2 <- 1e30 }

      # Compute incidence likelihoods
      ll_C <- ifelse(is.na(incidence_observed_C), 0,
                     dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE))

      ll_A <- ifelse(is.na(incidence_observed_A), 0,
                     dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE))

      # Prevalence Likelihood for Children (C)
      prev_observed_C <- observed$prev_C
      prev_model_C <- state["prev_C_1", , drop = TRUE]
      mu_prev_C <- observed$prev_C  # Center for C
      mu_prev_C <- prev_model_C  # Use the model's predicted prevalence
      kappa_C <- pars$kappa_C  # Precision parameter for Beta (specific to Children)

      if (is.na(prev_observed_C)) {
        ll_prev_C <- 0  # No contribution if missing
      } else {
        alpha_C <- mu_prev_C * kappa_C
        beta_C <- (1 - mu_prev_C) * kappa_C
        ll_prev_C <- dbeta(prev_observed_C, shape1 = alpha_C, shape2 = beta_C, log = TRUE)
      }

      # Prevalence Likelihood for Adults (A)
      prev_observed_A <- observed$prev_A
      prev_model_A <- state["prev_A_1", , drop = TRUE]
      mu_prev_A <- observed$prev_A  # Center for A
      kappa_A <- pars$kappa_A  # Precision parameter for Beta (specific to Adults)

      if (is.na(prev_observed_A)) {
        ll_prev_A <- 0  # No contribution if missing
      } else {
        alpha_A <- mu_prev_A * kappa_A
        beta_A <- (1 - mu_prev_A) * kappa_A
        ll_prev_A <- dbeta(prev_observed_A, shape1 = alpha_A, shape2 = beta_A, log = TRUE)
      }

      # Combine all likelihood components
      #return(ll_C + ll_A + 12 * (ll_prev_C + ll_prev_A))  # Adding prevalence likelihoods for C and A
      return(ll_C + ll_A + (ll_prev_C))  # Adding prevalence likelihoods for C and A
    })
  }

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
      cat("ll_C:", ll_C, " ll_A:", ll_A, " ll_prev_C:", ll_prev_C, " ll_prev_A:", ll_prev_A, " Total:", total_ll, "\n")

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

  if(month && age_for_inf == 'u5' & include_prev == FALSE) {
    return(function(state, observed, pars = c("size_1", "size_2")) {
      # Incidence Log-Likelihood
      incidence_observed_C <- observed$inc_C
      mu_C <- state["month_inc_C", , drop = TRUE]
      size_1 <- pars$size_1
      if (size_1 == 0) { size_1 <- 1e5 }
      if(is.na(incidence_observed_C)){ll_C <- 0}else{ # if no observed value, does not contribute to likelihood
        ll_C <- dnbinom(x = incidence_observed_C, mu = mu_C, size = size_1, log = TRUE)
      }
      return(ll_C)  # Use weighted penalty
    })
  }


  # if(month && age_for_inf == 'u5') {
  #   return(function(state, observed, pars) {
  #     incidence_observed_C <- observed$inc_C # this is the observed data
  #     mu_C <- state["month_inc_C", , drop = TRUE] # this is "x"
  #     return(dpois(x = incidence_observed_C, lambda = mu_C, log = TRUE))
  #   })
  # }

  if(month && age_for_inf == 'o5') {
    return(function(state, observed, pars = c("size_2")) {
      incidence_observed_A <- observed$inc_A
      mu_A <- state["month_inc_A", , drop = TRUE]
      size_2 <- pars$size_2
      if (size_2 == 0) { size_2 <- 1e30 }
      ll_A <- dnbinom(x = incidence_observed_A, mu = mu_A, size = size_2, log = TRUE)
      if(is.na(ll_A)){ll_A <- 0}
      return(ll_A)
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
