#' Generate a Comparison Function for Malaria Incidence
#'
#' This function generates a likelihood comparison function based on the specified
#' time period (monthly or weekly) and age group for malaria incidence.
#'
#' @param month A logical value indicating whether to use monthly data (TRUE) or weekly data (FALSE).
#' @param age_for_inf A string specifying the age group for incidence:
#'        "total", "sep_ages", "all_ages", "u5", or "o5".
#' @param incidence_df A data frame containing observed malaria incidence data.
#'
#' @return A function for comparing modeled and observed incidence based on a specified distribution (e.g., Negative Binomial).
#' @export
#'
#' @examples
#' # Example usage for monthly total incidence comparison:
#' comparison_fn <- generate_incidence_comparison(month = TRUE, age_for_inf = "total", incidence_df = incidence_data)
#' # Then use the returned function in MCMC or other model fitting tasks
generate_incidence_comparison <- function(month, age_for_inf, incidence_df) {
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
    return(function(state, observed, pars = c("size")) {
      incidence_observed_total <- observed$inc
      mu_total <- state["month_inc_total", , drop = TRUE]  # Cumulative monthly incidence
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(dnbinom(x = incidence_observed_total, mu = mu_total, size = size, log = TRUE))
    })
  }


  if(month && age_for_inf == 'sep_ages') {
    return(function(state, observed, pars = c("size")) {
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["month_inc_C", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(0.12 * dnbinom(x = incidence_observed_C, mu = mu_C, size = size, log = TRUE) +
               0.88 * dnbinom(x = incidence_observed_A, mu = mu_A, size = size, log = TRUE))
    })
  }

  if(month && age_for_inf == 'all_ages') {
    return(function(state, observed, pars = c("size")) {
      incidence_observed_C1 <- observed$inc_C1
      incidence_observed_C2 <- observed$inc_C2
      incidence_observed_A <- observed$inc_A
      mu_C1 <- state["month_inc_C1", , drop = TRUE]
      mu_C2 <- state["month_inc_C2", , drop = TRUE]
      mu_A <- state["month_inc_A", , drop = TRUE]
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(dnbinom(x = incidence_observed_C1, mu = mu_C1, size = size, log = TRUE) +
               dnbinom(x = incidence_observed_C2, mu = mu_C2, size = size, log = TRUE) +
               dnbinom(x = incidence_observed_A, mu = mu_A, size = size, log = TRUE))
    })
  }

  if(month && age_for_inf == 'u5') {
    return(function(state, observed, pars = c("size")) {
      if (is.na(observed$inc_C)) {
        return(NULL)
      }
      incidence_observed_C <- observed$inc_C # this is the observed data
      mu_C <- state["month_inc_C", , drop = TRUE] # this is "x"
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(dnbinom(x = incidence_observed_C, mu = mu_C, size = size, log = TRUE))
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
    return(function(state, observed, pars = c("size")) {
      incidence_observed_A <- observed$inc_A
      mu_A <- state["month_inc_A", , drop = TRUE]
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(dnbinom(x = incidence_observed_A, mu = mu_A, size = size, log = TRUE))
    })
  }

  # For weekly data
  if(!month && age_for_inf == 'total') {
    return(function(state, observed, pars = c("size")) {
      incidence_observed_total <- observed$inc
      mu_total <- state["wk_inc_total", , drop = TRUE]
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(dnbinom(x = incidence_observed_total, mu = mu_total, size = size, log = TRUE))
    })
  }

  if(!month && age_for_inf == 'sep_ages') {
    return(function(state, observed, pars = c("size")) {
      incidence_observed_C <- observed$inc_C
      incidence_observed_A <- observed$inc_A
      mu_C <- state["wk_inc_C", , drop = TRUE]
      mu_A <- state["wk_inc_A", , drop = TRUE]
      size <- pars$size
      if (size == 0) { size <- 1e30 }
      return(0.12 * dnbinom(x = incidence_observed_C, mu = mu_C, size = size, log = TRUE) +
               0.88 * dnbinom(x = incidence_observed_A, mu = mu_A, size = size, log = TRUE))
    })
  }

  # If no condition is met, return NULL or a default message
  stop("Invalid combination of 'month' and 'age_for_inf'.")
}
