

#' Loading an Odin Model
#'
#' @param name name of the model written in Odin DSL to be loaded. This is a path to an R file,
#'
#' @return the loaded model to be used for simulation, inference, etc
#' @export
#'
#' @examples
load_model <- function(name){
  file <- system.file("models", paste0(name, ".R"), package = "malclimsim")
  return(odin.dust::odin_dust(file))
}

# Helper function to model the decay of SMC (Seasonal Malaria Chemoprevention) efficacy over time.
# The decay follows a sigmoid curve where efficacy reduces as time (x) increases.
# Args:
#   - x: A numeric value representing the number of days since the last SMC treatment.
#   - const: The decay constant controlling the steepness of the sigmoid curve (default = -0.1806).
# Returns:
#   - A numeric value representing the remaining efficacy of SMC treatment at day 'x'.
decay_SMC <- function(x, const = -0.1806){
  return(0.90 / (1 + exp(const * (32.97 - x))))  # Sigmoid function for efficacy decay
}


calc_decay_arr <- function(SMC, decay_func = decay_SMC, const = -0.1806) {
  # Initialize decay array with zeros, length same as the input SMC
  decay_arr <- numeric(length(SMC))

  # Iterate over each point where SMC is administered
  i <- 1
  while (i <= length(SMC)) {
    if (SMC[i] > 0) {  # If SMC was administered
      # Find the index of the next SMC administration
      next_SMC <- which(SMC[(i + 1):length(SMC)] > 0)

      if (length(next_SMC) == 0) {
        # No subsequent SMC - apply decay until the end of the array
        end <- length(SMC)
      } else {
        # There is a subsequent SMC - apply decay until just before next SMC
        end <- i + next_SMC[1]
      }

      # Calculate decay values from index `i` to `end` (inclusive)
      decay_values <- decay_func(seq(0, end - i), const = const)

      # Ensure we do not exceed the boundaries of `decay_arr`
      decay_arr[i:end] <- decay_values[1:(end - i + 1)]
    }

    # Move to the next index
    i <- i + 1
  }

  # Return a decay array that has exactly the same length as `SMC`
  decay_arr <- decay_arr[1:length(SMC)]

  return(decay_arr)
}


# Function to run the malaria model simulation
# Arguments:
#   odin_mod: the model to simulate
#   pars: parameters for the model
#   time_start: initial time step
#   n_particles: number of particles (repetitions for stochastic simulation)
#   sim_time: total simulation time (days)
sim_mod <- function(odin_mod, pars, time_start, n_particles, sim_time){
  # Initialize the model with the provided parameters and particles
  model <- odin_mod$new(pars = pars, time = time_start, n_particles = n_particles)

  # Initialize an array to store simulation output for each time step
  x <- array(NA, dim = c(model$info()$len, n_particles, sim_time + 1))

  # Loop through each time step and run the model
  for (t in 0:sim_time) {
    x[ , , (t + 1)] <- model$run(t)  # Run the model at time t and store the result
  }

  # Return the simulation result array and the model object
  return(list(x, model))
}

#' Simulate Data from Model
#'
#' This function runs the simulation using the specified model and parameters,
#' and returns either monthly or weekly aggregated incidence data for adults (`inc_A`)
#' and children (`inc_C`), as well as total incidence (`inc`).
#'
#' @param model The model to be simulated.
#' @param param_inputs A list of parameters required for the simulation.
#' @param start_date The start date for the simulation, as a `Date` object or character string.
#' @param end_date The end date for the simulation, as a `Date` object or character string.
#' @param month Logical; if `TRUE`, aggregate data monthly; if `FALSE`, aggregate weekly.
#' @param round Logical; if `TRUE`, round the results to whole numbers.
#' @param save Logical; if `TRUE`, save the results to a file.
#' @param file Character; file name for saving the results, if `save` is `TRUE`.
#' @param month_unequal_days Logical; if `TRUE`, adjust aggregation for unequal days in months.
#' @return A data frame containing the incidence data aggregated by month or week.
#' @export
#' @examples
#' # Running the simulation and getting monthly data
#' inc_data <- data_sim(model, param_inputs, start_date = "2021-01-01", end_date = "2021-12-31", month = TRUE)
data_sim <- function(model, param_inputs, start_date, end_date,
                     month = FALSE, round = TRUE, save = TRUE, file = "",
                     month_unequal_days = FALSE, return_EIR = FALSE){
  # Calculate the number of days in the SMC schedule
  n_days <- calculate_360_day_difference(start_date, end_date) + 1

  # Run the simulation using the model and parameters
  results <- sim_mod(model, pars = c(param_inputs), time_start = 0,
                     n_particles = 1, sim_time = n_days)

  # Extract the simulation output and the model
  x <- results[[1]]
  mod <- results[[2]]

  # If monthly aggregation is selected
  if(month){
    # Create monthly indices for aggregation
    month_ind <- seq(1, n_days, by = 30)

    # Adjust for months with unequal days, if specified
    if(month_unequal_days){
      month_ind <- which(param_inputs$day_count == 0)
    }

    # Extract monthly incidence for adults (inc_A) and children (inc_C)
    inc_A <- x[mod$info()$index$month_inc_A,,][month_ind]
    inc_C <- x[mod$info()$index$month_inc_C,,][month_ind]

    # Get the corresponding monthly dates
    month <- date_to_months(start_date = as.Date(start_date), end_date = as.Date(end_date))

    # Ensure the lengths are equal
    n_months <- min(length(month), length(inc_A), length(inc_C))
    month <- month[1:n_months]
    inc_A <- inc_A[1:n_months]
    inc_C <- inc_C[1:n_months]

    # Create a dataframe for monthly incidence data
    month_no <- 0:(n_months - 1)

    if(return_EIR){
      EIR <- x[mod$info()$index$EIR2,,][month_ind]
      EIR <- EIR[1:n_months]
      inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C,
                           inc = inc_A + inc_C, EIR = EIR)
    } else{inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C)}

  } else {
    # If weekly aggregation is selected
    # Create weekly indices for aggregation
    wk_ind <- seq(1, n_days, by = 7)

    # Extract weekly incidence for adults (inc_A) and children (inc_C)
    inc_A <- x[mod$info()$index$wk_inc_A,,][wk_ind]
    inc_C <- x[mod$info()$index$wk_inc_C,,][wk_ind]

    # Get the corresponding weekly dates
    week <- date_to_weeks(start_date = as.Date(start_date), end_date = as.Date(end_date))

    # Ensure the lengths are equal
    n_weeks <- min(length(week), length(inc_A), length(inc_C))
    week <- week[1:n_weeks]
    inc_A <- inc_A[1:n_weeks]
    inc_C <- inc_C[1:n_weeks]

    # Create a dataframe for weekly incidence data
    week_no <- 0:(n_weeks - 1)

    if(return_EIR){
      EIR <- x[mod$info()$index$EIR2,,][wk_ind]
      EIR <- EIR[1:n_weeks]
      inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C, EIR = EIR)
    } else{inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C)}
  }

  # Round the results if specified
  if(round){
    inc_df[3:5] <- round(inc_df[3:5])
  }

  # Save the dataframe to a file if specified
  if(save){
    saveRDS(inc_df, paste0(dir, file))
  }

  # Return the resulting incidence dataframe
  return(inc_df)
}


data_sim_for_inference <- function(model, param_inputs,
                                   dates, noise = FALSE, month = FALSE,
                                   month_unequal_days = FALSE){

  # incidence dataframe in specific format from data_sim
  if(month){
    incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                             end_date = dates[2], month = TRUE, save = FALSE,
                             month_unequal_days = month_unequal_days)
  }else{incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                                 end_date = dates[2], month = FALSE, save = FALSE,
                                 month_unequal_days = month_unequal_days)}

  if(noise){
    set.seed(seed)
    #incidence_df$inc <- rpois(nrow(incidence_df), lambda = incidence_df$inc)
    incidence_df$inc <- rnbinom(nrow(incidence_df), mu = incidence_df$inc, size = 100)
  }
  return(incidence_df)
}

#' Extract Parameters with Maximum Log Posterior
#'
#' This function finds and extracts the parameter set corresponding to the maximum log posterior value
#' from the MCMC results.
#'
#' @param results The MCMC results object containing parameter samples.
#' @return A named vector of parameter values corresponding to the maximum log posterior.
#' @export
#' @examples
#' # Assuming `results` contains the MCMC output
#' max_posterior_params <- extract_max_posterior_params(results)
extract_max_posterior_params <- function(results) {
  coda_pars <- results$coda_pars

  # Find the index of the row with the maximum log posterior
  max_posterior_index <- which.max(coda_pars[ , "log_posterior"])

  # Extract the row with the maximum log posterior
  max_posterior_row <- coda_pars[max_posterior_index, ]

  # Remove log_prior, log_likelihood, and log_posterior columns
  param_values <- max_posterior_row[!names(max_posterior_row) %in%
                                      c("log_prior", "log_likelihood", "log_posterior")]

  return(param_values)
}

#' Update Parameters in List with Values from Extracted Vector
#'
#' This function updates values in a list of parameters (e.g., `param_inputs`)
#' with values from a named vector of parameters obtained from an MCMC run.
#'
#' @param param_inputs A list of parameters where each element is a named parameter.
#' @param param_values A named vector of parameters, typically from the maximum
#' log posterior, to update values in `param_inputs`.
#' @return An updated list of parameters with replaced values where matches exist.
#' @export
#' @examples
#' # Assuming `param_inputs` is a list of parameters and `params_at_max_posterior`
#' # is a named vector with values to update
#' updated_params <- update_param_list(param_inputs, params_at_max_posterior)
update_param_list <- function(param_inputs, param_values) {
  # Loop over each element in the named vector 'param_values'
  for (param_name in names(param_values)) {
    # Check if the parameter exists in the list
    if (param_name %in% names(param_inputs)) {
      # Update the list's parameter with the value from 'param_values'
      param_inputs[[param_name]] <- param_values[[param_name]]
    }
  }
  return(param_inputs)
}

#' Simulate Model Using Parameters with Maximum Log Posterior
#'
#' This function takes the MCMC results, extracts the parameter values corresponding to the maximum
#' log posterior, updates the parameter list, and runs a simulation from the model. Additionally,
#' the function includes a pre-warm period, where the simulation starts earlier than the specified
#' `start_date` to allow the model to stabilize before the output period.
#'
#' @param results The MCMC results object containing parameter samples, including `log_posterior`.
#' @param start_date The start date for the desired simulation output, as a `Date` object or character string.
#' @param end_date The end date for the simulation, as a `Date` object or character string.
#' @param model The model function to simulate from. This function should accept parameters, start/end dates,
#'              and additional arguments for the simulation.
#' @param prewarm_years Integer, the number of years to simulate before `start_date` as a pre-warm period
#'                      (default is 2 years).
#' @param days_per_year Integer, the number of days in a year for the simulation (default is 360 days).
#'
#' @return A data frame containing the simulation results, filtered to include only the period
#'         from `start_date` to `end_date`.
#' @export
#'
#' @examples
#' # Assuming `results` contains the MCMC output, and `data_sim` is the simulation function
#' simulation_output <- simulate_with_max_posterior_params(
#'   results = results,
#'   start_date = "2021-01-01",
#'   end_date = "2021-12-31",
#'   model = data_sim,
#'   prewarm_years = 3,
#'   days_per_year = 360
#' )
simulate_with_max_posterior_params <- function(results, start_date, end_date, model, prewarm_years = 2, days_per_year = 360) {

  param_inputs <- results$param_inputs

  # Setup for allowing model to run some prior to inference (comparing to observations)
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = 360,
                                                 years_to_extend = prewarm_years)

  # Extract parameters with maximum log posterior
  max_posterior_params <- extract_max_posterior_params(results)

  # Update parameter list with the extracted parameters
  updated_params <- update_param_list(param_inputs_ext, max_posterior_params)

  # Extend the start date backward by the pre-warm period
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Run the model simulation for the pre-warm period + desired period
  extended_simulation_output <- data_sim(
    model = model,
    param_inputs = updated_params,
    start_date = prewarm_start_date,
    end_date = end_date,
    month = TRUE,
    round = FALSE,
    save = FALSE,
    month_unequal_days = FALSE
  )

  # Filter the simulation results to include only the desired period
  simulation_output <- extended_simulation_output[extended_simulation_output$date_ymd >= as.Date(start_date), ]

  return(simulation_output)
}



#' Sample Parameters from MCMC Results
#'
#' This function randomly selects `n` rows from the MCMC parameter output and extracts the
#' parameter values from each row.
#'
#' @param results The MCMC results object containing parameter samples.
#' @param n The number of samples to extract.
#' @return A list of named vectors, each containing one set of sampled parameters.
#' @examples
#' # Assuming `results` contains the MCMC output
#' sampled_parameters <- sample_params(results, 100)
sample_params <- function(results, n) {
  # Randomly sample `n` rows
  sampled_indices <- sample(1:nrow(results$coda_pars), size = n, replace = TRUE)
  sampled_params <- lapply(sampled_indices, function(idx) {
    # Extract row and remove columns for log_prior, log_likelihood, and log_posterior
    params <- results$coda_pars[idx, ]
    params <- params[!names(params) %in% c("log_prior", "log_likelihood", "log_posterior")]
    return(params)
  })
  return(sampled_params)
}

#' Simulate Models Using Sampled Parameters
#'
#' This function runs `n` model simulations, each using a different set of parameters sampled
#' from the MCMC results.
#'
#' @param model The model function to simulate from.
#' @param sampled_params A list of parameter sets to use for each simulation.
#' @param start_date The start date for the simulation.
#' @param end_date The end date for the simulation.
#' @param param_inputs The initial list of parameters to update for each simulation.
#' @return A list of data frames, each containing the simulation output for a different parameter set.
#' @examples
#' # Assuming `sampled_params` is a list of sampled parameter sets
#' simulations <- simulate_models(model = data_sim, param_inputs = param_inputs,
#'                                sampled_params = sampled_params,
#'                                start_date = "2021-01-01", end_date = "2021-12-31")
simulate_models <- function(model, param_inputs, sampled_params, start_date, end_date) {
  simulations <- lapply(sampled_params, function(params) {
    # Make a copy of param_inputs for each simulation to prevent overwriting
    current_params <- param_inputs

    # Update the copied parameters with the sampled parameters
    updated_params <- update_param_list(current_params, params)

    # Run the simulation with the updated parameters
    sim_data <- data_sim(model, updated_params, start_date = start_date, end_date = end_date,
                         month = TRUE, round = FALSE, save = FALSE, month_unequal_days = FALSE)
    return(sim_data)
  })
  return(simulations)
}

#' Calculate Quantiles for Incidence from Multiple Simulations
#'
#' This function calculates the 1st and 99th quantiles of the simulated incidence data across multiple simulations,
#' for each group (`inc_A`, `inc_C`, and `inc`).
#'
#' @param simulations A list of data frames, each containing the simulation output for different parameter sets.
#' @return A data frame containing the quantiles for each group (`inc_A`, `inc_C`, and `inc`) for each time point.
#' @examples
#' # Assuming `simulations` is a list of data frames with simulated incidence data
#' incidence_quantiles <- calculate_incidence_quantiles(simulations)
calculate_incidence_quantiles <- function(simulations) {
  # Extract incidence data for each group (inc_A, inc_C, inc) across all simulations
  inc_A_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc_A))
  inc_C_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc_C))
  inc_data <- do.call(cbind, lapply(simulations, function(sim) sim$inc))

  # Calculate quantiles along rows (for each time point) for each group
  quantiles_inc_A <- apply(inc_A_data, 1, function(x) {
    c(quantile(x, probs = 0.01), quantile(x, probs = 0.99))
  })
  quantiles_inc_C <- apply(inc_C_data, 1, function(x) {
    c(quantile(x, probs = 0.01), quantile(x, probs = 0.99))
  })
  quantiles_inc <- apply(inc_data, 1, function(x) {
    c(quantile(x, probs = 0.01), quantile(x, probs = 0.99))
  })

  # Format the output into a data frame
  quantiles_df <- data.frame(
    date_ymd = simulations[[1]]$date_ymd,
    inc_A_q01 = quantiles_inc_A[1, ],
    inc_A_q99 = quantiles_inc_A[2, ],
    inc_C_q01 = quantiles_inc_C[1, ],
    inc_C_q99 = quantiles_inc_C[2, ],
    inc_q01 = quantiles_inc[1, ],
    inc_q99 = quantiles_inc[2, ]
  )

  return(quantiles_df)
}

