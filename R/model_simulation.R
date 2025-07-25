#' Loading an Odin Model
#'
#' @param name name of the model written in Odin DSL to be loaded. This is a path to an R file,
#' @return the loaded model to be used for simulation, inference, etc
#' @export
#'
#' @examples
#' load_model("model_new_R_with_FOI")
load_model <- function(name){
  file <- system.file("models", paste0(name, ".R"), package = "malclimsim")
  return(odin.dust::odin_dust(file))
}


#' Import an Odin Model into the malclimsim Package
#'
#' @param model_path Full file path to the external Odin model (an R file) to be imported.
#' @param model_name Desired name (without ".R") for the model when stored inside the package.
#'
#' @return None. The function copies the model file into the package's models directory.
#' @export
#'
#' @examples
#' import_model("path/to/external_model.R", "my_new_model")
import_model <- function(model_path, model_name){
  # Construct the path to the "models" directory within the installed package
  package_models_path <- paste0(find.package("malclimsim"), "/models/")

  # Construct the destination path for the new model file within the package
  new_model_new_path <- paste0(package_models_path, model_name, ".R")

  # Copy the external model file to the destination path inside the package's models folder
  file.copy(from = model_path, to = new_model_new_path, overwrite = TRUE)
}


#' Compute SMC Efficacy Decay Over Time
#'
#' Models the decay of Seasonal Malaria Chemoprevention (SMC) efficacy as a sigmoid function
#' of time since last administration.
#'
#' @param x Numeric. Number of days since the last SMC treatment.
#' @param const Numeric. The decay constant that controls the steepness of the sigmoid curve. Default is -0.1806.
#'
#' @return Numeric. Remaining efficacy of SMC treatment at day \code{x}.
#'
#' @examples
#' decay_SMC(0)      # Should return near maximum efficacy
#' decay_SMC(60)     # Should return reduced efficacy
#'
#' @export
decay_SMC <- function(x, const = -0.1806){
  return(0.90 / (1 + exp(const * (32.97 - x))))
}

#' Calculate Normalized SMC Decay Array Over Time
#'
#' Generates a vector of normalized decayed SMC efficacies for a time series where SMC is applied intermittently.
#' Each decay curve is scaled so that its maximum value is 1.
#'
#' @param SMC Numeric vector. Time series indicating the presence (positive value) or absence (zero) of SMC at each time point.
#' @param decay_func Function. Function defining the decay curve. Default is \code{decay_SMC}.
#' @param const Numeric. Decay constant to pass to \code{decay_func}. Default is -0.1806.
#'
#' @return Numeric vector of the same length as \code{SMC}, with decay values normalized to a max of 1 per SMC round.
#'
#' @examples
#' smc_series <- c(0, 1, 0, 0, 0, 1, 0, 0)
#' calc_decay_arr(smc_series)
#'
#' @export
calc_decay_arr <- function(SMC, decay_func = decay_SMC, const = -0.1806) {
  decay_arr <- numeric(length(SMC))
  i <- 1
  while (i <= length(SMC)) {
    if (SMC[i] > 0) {
      next_SMC <- which(SMC[(i + 1):length(SMC)] > 0)
      end <- if (length(next_SMC) == 0) length(SMC) else i + next_SMC[1]
      decay_values <- decay_func(seq(0, end - i), const = const)
      #decay_values <- decay_values / max(decay_values, na.rm = TRUE)  # Normalize so max is 1
      decay_arr[i:end] <- decay_values[1:(end - i + 1)]
    }
    i <- i + 1
  }
  return(decay_arr)
}

#' Simulate the Malaria Model Over Time
#'
#' Runs a stochastic malaria transmission model using an `odin`-generated model object, storing state outputs across multiple time steps and particles.
#'
#' @param odin_mod An `odin` model object (usually created with `odin::odin()`), representing the compiled malaria model.
#' @param pars A named list of parameters to pass to the model.
#' @param time_start Numeric value specifying the starting time (e.g., 0).
#' @param n_particles Integer; number of particles used in the stochastic simulation (i.e., number of replicate simulations).
#' @param sim_time Integer; total number of time steps (days) to simulate.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{x}{A 3D array of simulation results with dimensions `[state, particle, time]`.}
#'   \item{model}{The `odin` model object after simulation (can be used to inspect internals).}
#' }
#'
#' @details
#' The model is initialized once using the `odin_mod$new()` constructor, and then `model$run(t)` is called at each time step `t` to advance the simulation.
#'
#' @examples
#' \dontrun{
#' mod <- odin::odin("malaria_model.R")
#' results <- sim_mod(mod, pars = list(beta = 0.2), time_start = 0, n_particles = 10, sim_time = 100)
#' }
#'
#' @export
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

#' Simulate Incidence Data Over Time (Monthly or Weekly)
#'
#' Runs a deterministic dust model simulation with an optional pre‑warm period,
#' then aggregates and returns incidence (and optional transformed incidence)
#' at either monthly or weekly resolution.  Weekly output may use a 360‑day or
#' 365‑day year for the prewarm and simulation period.
#'
#' @param model A dust model object (as returned by \code{load_model}) to simulate.
#' @param param_inputs Named list of model parameters, including any time‑varying
#'   vectors (e.g. \code{temp}, \code{c_R_D}, \code{cov_SMC}).
#' @param start_date Character or \code{Date} giving the first day of the analysis.
#' @param end_date   Character or \code{Date} giving the last day of the analysis.
#' @param prewarm_years Integer number of years (360‑ or 365‑day) to run before
#'   \code{start_date} as warm‑up.
#' @param month Logical; if \code{TRUE}, returns monthly aggregates (30‑day or
#'   unequal month lengths if \code{month_unequal_days = TRUE}); if \code{FALSE},
#'   returns weekly aggregates.
#' @param round Logical; if \code{TRUE}, round all incidence counts to integers.
#' @param save Logical; if \code{TRUE}, save the returned data frame to disk as
#'   an RDS file.
#' @param file Character; path/filename to use when saving (if \code{save = TRUE}).
#' @param month_unequal_days Logical; when \code{month = TRUE}, if \code{TRUE}
#'   use actual month boundaries from \code{param_inputs\$day_count}, otherwise
#'   use fixed 30‑day months.
#' @param return_EIR Logical; if \code{TRUE}, include monthly EIR values.
#' @param return_compartments Logical; currently unused.
#' @param mu_transform_A Optional function to transform adult incidence after
#'   simulation.  Should accept \code{(inc_df, param_inputs)} and return a vector.
#' @param mu_transform_C Optional function to transform child incidence after
#'   simulation.  Should accept \code{(inc_df, param_inputs)} and return a vector.
#' @param covariate_matrix Optional data frame of external covariates (must have
#'   a date column in the same name/location as \code{date_ymd} in the output);
#'   will be joined by \code{date_ymd}.
#' @param noise Logical; if \code{TRUE}, add negative‑binomial noise to each
#'   weekly or monthly incidence draw, using \code{size} as dispersion.
#' @param size Optional numeric; dispersion parameter for negative‑binomial noise.
#'
#' @return A data frame with columns:
#'   \item{date_ymd}{Date of each week (Sunday) or month (start day).}
#'   \item{week_no, month_no}{Index of each period, starting at 0.}
#'   \item{inc_A, inc_C, inc}{Adult, child, and total incidence.}
#'   \item{EIR_monthly}{Monthly EIR, if \code{return_EIR = TRUE}.}
#'   \item{inc_A_transformed, inc_C_transformed}{Transformed incidence, if
#'     \code{mu_transform_A} or \code{mu_transform_C} supplied.}
#'   plus any merged covariates from \code{covariate_matrix}.
#'
#' @examples
#' \dontrun{
#' sim_df <- data_sim(
#'   model             = model_1,
#'   param_inputs      = param_inputs_wk,
#'   start_date        = "2018-01-01",
#'   end_date          = "2023-12-31",
#'   prewarm_years     = 5,
#'   month             = FALSE,
#'   round             = FALSE,
#'   save              = FALSE,
#'   mu_transform_C    = mu_transform_C_wk,
#'   covariate_matrix  = covariate_matrix_wk,
#'   noise             = TRUE,
#'   size              = 10,
#'   year_360_days     = FALSE
#' )
#' }
#' @export
data_sim <- function(model, param_inputs, start_date, end_date,
                     prewarm_years = 2, month = FALSE, round = TRUE, save = TRUE, file = "",
                     month_unequal_days = FALSE, return_EIR = FALSE, return_compartments = FALSE,
                     mu_transform_A = NULL, mu_transform_C = NULL, covariate_matrix = NULL,
                     noise = FALSE, size = NULL) {

  # --- Determine days_per_year & total days ---
  if (month) {
    days_per_year <- 360
    # prewarm start in model calendar
    prewarm_start <- paste0(year(start_date) - prewarm_years, "-", format.Date(start_date, "%m-%d"))
    n_days    <- calculate_360_day_difference(prewarm_start, end_date) + 1
  } else {
    days_per_year <- if (month_unequal_days) 365 else 360
    prewarm_start <- as.Date(start_date) - prewarm_years * days_per_year
    n_days    <- as.integer(as.Date(end_date) - as.Date(prewarm_start)) + 1
  }

  # --- Extend time-varying inputs ---
  param_inputs <- extend_time_varying_inputs(
    param_inputs,
    days_per_year = days_per_year,
    years_to_extend = prewarm_years
  )

  # --- Run model ---
  results <- sim_mod(model, pars = param_inputs, time_start = 0,
                     n_particles = 1, sim_time = n_days)
  x   <- results[[1]]
  mod <- results[[2]]$info()$index

  # --- Aggregate ---
  if (month) {
    month_ind <- seq(1, n_days, by = 30)
    if (month_unequal_days) month_ind <- which(param_inputs$day_count == 0)

    inc_A <- x[mod$month_inc_A,,][month_ind]
    inc_C <- x[mod$month_inc_C,,][month_ind]
    dates <- date_to_months(prewarm_start, end_date)

    if(noise){
      inc_C <- rnbinom(length(inc_C), size = size, mu = inc_C)
      inc_A <- rnbinom(length(inc_A), size = size, mu = inc_A)
    }

    idx <- seq_len(min(length(dates), length(inc_A), length(inc_C)))
    inc_df <- data.frame(
      date_ymd = dates[idx], month_no = 0:(length(idx)-1),
      inc_A = inc_A[idx], inc_C = inc_C[idx],
      inc   = inc_A[idx] + inc_C[idx]
    )
    if (return_EIR) {
      EIR <- x[mod$EIR_monthly,,][month_ind][idx]
      inc_df$EIR_monthly <- EIR
    }
  } else {
    week_ind <- seq(1, n_days, by = 7)

    inc_A <- x[mod$wk_inc_A,,][week_ind]
    inc_C <- x[mod$wk_inc_C,,][week_ind]

    if(noise){
      inc_C <- rnbinom(length(inc_C), size = size, mu = inc_C)
      inc_A <- rnbinom(length(inc_A), size = size, mu = inc_A)
    }

    # use real calendar Sundays if using 365
    dates <- date_to_weeks(prewarm_start, end_date)
    #if (use_365_for_weeks) {
    #  dates <- date_to_weeks(prewarm_start, end_date)
    #} else {
    #  dates <- date_to_weeks(prewarm_start, n_days)  # your 360-day helper
    #}

    idx <- seq_len(min(length(dates), length(inc_A), length(inc_C)))
    inc_df <- data.frame(
      date_ymd = dates[idx], week_no = 0:(length(idx)-1),
      inc_A = inc_A[idx], inc_C = inc_C[idx],
      inc   = inc_A[idx] + inc_C[idx]
    )
    if (return_EIR) {
      EIR <- x[mod$EIR_monthly,,][week_ind][idx]
      inc_df$EIR_monthly <- EIR
    }
  }

  #inc_df$anom <- param_inputs$c_R_D[seq(1, length(param_inputs$c_R_D), by = 7)]
  #inc_df$temp <- param_inputs$temp[seq(1, length(param_inputs$temp), by = 7)]

  # --- Covariate merge & transform ---
  if (!is.null(covariate_matrix)) {
    inc_df <- dplyr::left_join(
      inc_df, covariate_matrix,
      by = c("date_ymd" = names(covariate_matrix)[1])
    )
  }
  if (!is.null(mu_transform_A)) inc_df$inc_A_transformed <- mu_transform_A(inc_df, param_inputs)
  if (!is.null(mu_transform_C)) inc_df$inc_C_transformed <- mu_transform_C(inc_df, param_inputs)

  # --- Round, trim, save ---
  if (round) inc_df[c("inc_A","inc_C","inc")] <- round(inc_df[c("inc_A","inc_C","inc")])
  inc_df <- inc_df[inc_df$date_ymd >= start_date & inc_df$date_ymd <= end_date, ]
  if (save) saveRDS(inc_df, file)

  return(inc_df)
}


#' Simulate Data for Inference
#'
#' This function generates a simulated dataset for inference based on a specified model and parameter values.
#' It supports options for monthly aggregation, adjusting for unequal days in months, and adding stochastic noise
#' to the incidence data.
#'
#' @param model A compiled model object compatible with `data_sim()`. Typically created using `odin.dust::odin_dust()`.
#' @param param_inputs A named list of parameter values required by the model.
#' @param dates A character vector of length two specifying the start and end dates (e.g., \code{c("2022-01-01", "2022-12-31")}).
#' @param noise Logical. If \code{TRUE}, random noise is added using a negative binomial distribution. Defaults to \code{FALSE}.
#' @param month Logical. If \code{TRUE}, the output is aggregated by month. Defaults to \code{FALSE}.
#' @param month_unequal_days Logical. If \code{TRUE}, accounts for unequal month lengths when aggregating. Defaults to \code{FALSE}.
#'
#' @return A data frame of simulated incidence data, in the same format as the output of \code{data_sim()}.
#' If \code{noise = TRUE}, the \code{inc} column will include added random variation.
#'
#' @details
#' - When \code{month = TRUE}, output is aggregated by month. If \code{month_unequal_days = TRUE}, aggregation is adjusted
#'   to reflect actual month lengths.
#' - If \code{noise = TRUE}, the function adds variability to the incidence estimates using a negative binomial distribution
#'   with size = 100.
#'
#' @examples
#' # Example (requires a model created using odin.dust, not run here):
#' \dontrun{
#' model <- odin.dust::odin_dust("path/to/model.R")
#' param_inputs <- list(beta = 0.3, gamma = 0.1)
#' dates <- c("2022-01-01", "2022-12-31")
#'
#' # Simulate data without noise
#' sim_data <- data_sim_for_inference(model, param_inputs, dates, noise = FALSE)
#'
#' # Simulate data with noise and monthly aggregation
#' sim_data_monthly <- data_sim_for_inference(model, param_inputs, dates,
#'                                            noise = TRUE, month = TRUE)
#' }
#'
#' @seealso \code{\link{data_sim}} for the underlying simulation function.
#' @export
data_sim_for_inference <- function(model, param_inputs,
                                   dates, noise = FALSE, month = FALSE,
                                   month_unequal_days = FALSE){

  # incidence dataframe in specific format from data_sim
  incidence_df <- data_sim(model, param_inputs, start_date = dates[1],
                           end_date = dates[2], month = month, save = FALSE,
                           month_unequal_days = month_unequal_days)

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
#' from MCMC output. It assumes that the MCMC results object includes a `coda_pars` matrix or data frame
#' with columns for `log_prior`, `log_likelihood`, and `log_posterior`, along with parameter names.
#'
#' @param results A list containing MCMC output, with a `coda_pars` element (a matrix or data frame)
#' that includes a column named `"log_posterior"` and columns for sampled parameters.
#'
#' @return A named numeric vector of parameter values corresponding to the maximum log posterior.
#'
#' @examples
#' # Create dummy MCMC output for illustration
#' set.seed(123)
#' results <- list(
#'   coda_pars = data.frame(
#'     param1 = rnorm(100),
#'     param2 = rnorm(100),
#'     log_prior = rnorm(100),
#'     log_likelihood = rnorm(100),
#'     log_posterior = rnorm(100)
#'   )
#' )
#'
#' # Extract parameters corresponding to the maximum log posterior
#' extract_max_posterior_params(results)
#'
#' @export
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

#' Update a Parameter List with New Values
#'
#' Replaces entries in a parameter list using a named vector of new values. Typically used to inject
#' sampled or optimized parameter values (e.g., from an MCMC posterior) into a base parameter list.
#'
#' @param param_inputs A named list containing parameter values, where each name corresponds to a model parameter.
#' @param param_values A named vector or list of parameter values to use for updating. Each name in this object
#'   should correspond to a name in \code{param_inputs}.
#'
#' @return A modified version of \code{param_inputs} where parameters matching names in \code{param_values}
#'   are replaced by the corresponding values.
#'
#' @examples
#' # Original parameter list
#' param_inputs <- list(alpha = 1, beta = 2, gamma = 3)
#'
#' # New parameter values, e.g., from posterior sample
#' params_at_max_posterior <- c(alpha = 0.9, gamma = 2.7)
#'
#' # Update the parameter list
#' updated_params <- update_param_list(param_inputs, params_at_max_posterior)
#' print(updated_params)
#'
#' @export
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
#' Runs a model simulation using the parameter set corresponding to the highest posterior density
#' (i.e., the maximum log posterior value) from the MCMC output. The function includes a prewarm period
#' prior to the specified `start_date` to stabilize the system before collecting simulation results.
#'
#' @param results A list containing MCMC results. Must include `log_posterior`, `param_inputs`,
#'                and samples of estimated parameters.
#' @param start_date Simulation start date (as `Date` or character string) for the output period.
#' @param end_date Simulation end date (as `Date` or character string).
#' @param model The model function to simulate. Should be compatible with `data_sim()`.
#' @param prewarm_years Number of years to simulate prior to `start_date` to allow system stabilization (default is 2).
#' @param days_per_year Number of days in a simulation year (default is 360).
#' @param mu_transform_A Optional transformation function for adult incidence (default is `NULL`).
#' @param mu_transform_C Optional transformation function for child incidence (default is `NULL`).
#' @param covariate_matrix Optional covariate matrix used in incidence transformation (default is `NULL`).
#'
#' @return A data frame containing model simulation results from `start_date` to `end_date`.
#'
#' @examples
#' \dontrun{
#' # Run simulation using max-posterior parameter set
#' sim_df <- simulate_with_max_posterior_params(
#'   results = results,
#'   start_date = "2021-01-01",
#'   end_date = "2021-12-31",
#'   model = data_sim
#' )
#' }
#'
#' @export
simulate_with_max_posterior_params <- function(results, start_date, end_date, model,
                                               prewarm_years = 2, days_per_year = 360,
                                               mu_transform_A = NULL,
                                               mu_transform_C = NULL,
                                               covariate_matrix = NULL) {
  param_inputs <- results$param_inputs
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = 360, years_to_extend = prewarm_years)
  max_posterior_params <- extract_max_posterior_params(results)
  updated_params <- update_param_list(param_inputs_ext, max_posterior_params)
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  simulation_output <- data_sim(
    model = model,
    param_inputs = updated_params,
    start_date = prewarm_start_date,
    end_date = end_date,
    month = TRUE,
    round = FALSE,
    save = FALSE,
    month_unequal_days = FALSE,
    mu_transform_A = mu_transform_A,
    mu_transform_C = mu_transform_C,
    covariate_matrix = covariate_matrix
  )

  simulation_output <- simulation_output[simulation_output$date_ymd >= as.Date(start_date), ]
  return(simulation_output)
}

#' Simulate Models Using Sampled Parameters
#'
#' Runs the epidemiological model multiple times using different parameter sets sampled from MCMC results.
#' Each simulation is optionally prewarmed using a specified number of years prior to the start date.
#'
#' @param model The model function to simulate from.
#' @param param_inputs A list of baseline parameters to be updated with each sampled parameter set.
#' @param sampled_params A list of named parameter sets, typically obtained from posterior samples.
#' @param start_date Character string giving the start date for the simulation (e.g., "2021-01-01").
#' @param end_date Character string giving the end date for the simulation (e.g., "2021-12-31").
#' @param prewarm_years Integer specifying how many years before `start_date` to simulate for model prewarming (default is 2).
#' @param days_per_year Integer specifying the number of days per simulation year (default is 360).
#'
#' @return A list of data frames. Each data frame contains the simulation results for a different parameter set.
#'
#' @examples
#' \dontrun{
#' # Example setup with a dummy model and parameters
#' model <- your_model_function  # Replace with your actual model function
#' param_inputs <- return_default_inputs()  # Replace with your actual input generator
#'
#' # Simulate posterior parameter samples
#' posterior_samples <- matrix(rnorm(30), nrow = 10, ncol = 3)
#' colnames(posterior_samples) <- c("alpha", "beta", "gamma")
#'
#' # Convert to list of named parameter sets
#' sampled_params <- lapply(1:nrow(posterior_samples), function(i) {
#'   as.list(posterior_samples[i, ])
#' })
#'
#' # Run simulations
#' simulations <- simulate_models(
#'   model = model,
#'   param_inputs = param_inputs,
#'   sampled_params = sampled_params,
#'   start_date = "2021-01-01",
#'   end_date = "2021-12-31"
#' )
#' }
#'
#' @export
simulate_models <- function(model, param_inputs, sampled_params, start_date, end_date, prewarm_years = 2, days_per_year = 360) {
  # Extend parameter inputs to accommodate the prewarm period
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = days_per_year, years_to_extend = prewarm_years)

  # Calculate the prewarm start date
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Run simulations with prewarming
  simulations <- lapply(sampled_params, function(params) {
    # Make a copy of param_inputs for each simulation to prevent overwriting
    current_params <- param_inputs_ext

    # Update the copied parameters with the sampled parameters
    updated_params <- update_param_list(current_params, params)

    # Run the simulation for the prewarm period and the desired period
    extended_sim_data <- data_sim(
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
    sim_data <- extended_sim_data[extended_sim_data$date_ymd >= as.Date(start_date), ]

    return(sim_data)
  })

  return(simulations)
}

#' Simulate Compartments Over Time with Prewarm Period
#'
#' This function simulates the compartments of an epidemiological model over a specified time period,
#' including an optional prewarm period to allow the model to reach equilibrium before the main simulation period.
#'
#' @param model An epidemiological model object used for simulation.
#' @param param_inputs A named list or vector of model parameters to be used in the simulation.
#' @param start_date A character string or Date object specifying the start date of the main simulation (e.g., "2014-01-01").
#' @param end_date A character string or Date object specifying the end date of the simulation (e.g., "2022-12-31").
#' @param prewarm_years Integer specifying the number of years to prewarm the model before the main simulation (default: 2).
#' @param days_per_year Integer specifying the number of days in a model year (default: 360).
#'
#' @return A data frame containing the simulated values of each compartment (SC, EC, IC, etc.)
#' over the simulation period with a weekly time resolution.
#'
#' @details The function simulates the dynamics of susceptible, exposed, infected, treated, and recovered
#' compartments for both children and adults. It first runs a prewarm period (if specified) to allow the model to stabilize,
#' before simulating the main study period. The results are aggregated at weekly intervals.
#'
#' @export
compartments_sim <- function(model, param_inputs, start_date, end_date, prewarm_years = 2, days_per_year = 360) {
  # Extend parameter inputs to accommodate the prewarm period
  param_inputs_ext <- extend_time_varying_inputs(param_inputs, days_per_year = days_per_year, years_to_extend = prewarm_years)

  # Calculate the prewarm start date
  prewarm_start_date <- paste0(year(as.Date(start_date)) - prewarm_years, "-", format(as.Date(start_date), "%m-%d"))

  # Calculate total simulation days
  total_days <- calculate_360_day_difference(prewarm_start_date, end_date) - 1

  # Run the simulation
  results <- sim_mod(model, pars = c(param_inputs_ext), time_start = 0,
                     n_particles = 1, sim_time = total_days)

  x <- results[[1]]
  mod <- results[[2]]

  # Weekly indices for output (skip prewarm)
  full_wk_ind <- seq(1, total_days, by = 7)
  start_index <- calculate_360_day_difference(prewarm_start_date, start_date)
  wk_ind <- full_wk_ind[full_wk_ind >= start_index]

  # Helper function for safe extraction
  extract_compartment <- function(index_name) {
    if (index_name %in% names(mod$info()$index)) {
      values <- x[mod$info()$index[[index_name]], , ]
      if (length(values) >= max(wk_ind)) {
        values[wk_ind]
      } else {
        stop(paste("Mismatch in compartment size for:", index_name))
      }
    } else {
      # Return NA vector if not present in model
      rep(NA_real_, length(wk_ind))
    }
  }


  # Vector of compartments to extract
  compartment_names <- c(
    "SC", "EC", "IC", "TrC", "RC", "SA", "EA", "IA", "TrA", "RA",
    "P_C", "P_A", "EIR2", "SMC_effect_2", "prev_total_1", "prev_C_1", "prev_A_1",
    "prev_total_2", "prev_C_2", "prev_A_2", "mu_SE_C_2", "mu_SE_A_2",
    "X2", "X_I", "X_AP", "X_ASP", "rain_effect_2", "temp_effect_2"
  )

  # Extract all compartments in a list
  compartments <- setNames(lapply(compartment_names, extract_compartment), compartment_names)

  # Rename for clarity (optional)
  names(compartments)[names(compartments) == "P_C"] <- "PC"
  names(compartments)[names(compartments) == "P_A"] <- "PA"
  names(compartments)[names(compartments) == "EIR2"] <- "EIR"
  names(compartments)[names(compartments) == "SMC_effect_2"] <- "SMC_effect"
  names(compartments)[names(compartments) == "prev_total_1"] <- "prev_total_with_R"
  names(compartments)[names(compartments) == "prev_C_1"] <- "prev_C_with_R"
  names(compartments)[names(compartments) == "prev_A_1"] <- "prev_A_with_R"
  names(compartments)[names(compartments) == "prev_total_2"] <- "prev_total_no_R"
  names(compartments)[names(compartments) == "prev_C_2"] <- "prev_C_no_R"
  names(compartments)[names(compartments) == "prev_A_2"] <- "prev_A_no_R"
  names(compartments)[names(compartments) == "mu_SE_C_2"] <- "mu_SE_C"
  names(compartments)[names(compartments) == "mu_SE_A_2"] <- "mu_SE_A"
  names(compartments)[names(compartments) == "X2"] <- "X"
  names(compartments)[names(compartments) == "rain_effect_2"] <- "rain_effect"
  names(compartments)[names(compartments) == "temp_effect_2"] <- "temp_effect"

  # Derived total population
  compartments$P <- compartments$PC + compartments$PA

  # Weekly dates
  compart_dates <- seq(as.Date(start_date), as.Date(end_date), by = "7 days")

  # Trim all vectors to minimum shared length
  min_length <- min(length(compart_dates), sapply(compartments, length))
  compart_dates <- compart_dates[1:min_length]
  compartments <- lapply(compartments, function(x) x[1:min_length])

  # Combine into final data frame
  compart_df <- data.frame(date = compart_dates, compartments)

  return(compart_df)
}

#' Sample Parameter Sets from MCMC Results
#'
#' Randomly samples rows from a matrix or data frame containing MCMC posterior samples. Useful for generating
#' parameter sets to use in posterior predictive simulations or scenario analyses.
#'
#' @param mcmc_results A matrix or data frame where each row represents a parameter set from the posterior distribution.
#' @param num_samples Integer specifying how many parameter sets to sample.
#'
#' @return A matrix of dimension \code{num_samples} × \code{ncol(mcmc_results)}, containing randomly selected parameter sets.
#'
#' @examples
#' # Simulate MCMC posterior with 100 samples and 3 parameters
#' posterior_samples <- matrix(rnorm(300), nrow = 100, ncol = 3)
#' colnames(posterior_samples) <- c("alpha", "beta", "gamma")
#'
#' # Draw 10 random samples from the MCMC results
#' sampled <- sample_mcmc_steps(posterior_samples, 10)
#' print(sampled)
#'
#' @export
sample_mcmc_steps <- function(mcmc_results, num_samples) {
  if (num_samples > nrow(mcmc_results)) {
    stop("Number of samples requested exceeds available MCMC steps.")
  }

  sampled_indices <- sample(1:nrow(mcmc_results), num_samples, replace = FALSE)
  return(mcmc_results[sampled_indices, , drop = FALSE])
}

#' Run Model Simulations for Each Parameter Set
#'
#' This function runs the epidemiological model using different parameter sets sampled from the MCMC posterior.
#'
#' @param model A function representing the epidemiological model, typically compiled via `odin.dust`.
#' @param param_inputs A list of default model parameters to be updated with sampled values.
#' @param param_samples A matrix or data frame of parameter sets sampled from the MCMC posterior. Each row is a distinct sample.
#' @param start_date Character string or Date object indicating the simulation start date (e.g., "2014-01-01").
#' @param end_date Character string or Date object indicating the simulation end date (e.g., "2022-12-31").
#' @param prewarm_years Number of years for model prewarming to reach equilibrium before the start date (default: 2).
#' @param days_per_year Number of simulation time steps per year (default: 360).
#'
#' @return A list of data frames, each containing the simulated compartment outputs for one parameter set.
#'
#' @examples
#' \dontrun{
#' # Assume model is loaded via `load_model()` and priors are defined
#' model <- load_model("model_det_1")
#'
#' # Get default parameter inputs and samples from MCMC
#' param_inputs <- return_default_param_inputs()
#' param_samples <- matrix(
#'   c(0, 0, 2, 5, 3, 0, 1, 0.5, 0.8, 1),  # Example parameter set (one row)
#'   nrow = 1,
#'   byrow = TRUE
#' )
#' colnames(param_samples) <- c("lag_R", "lag_T", "alpha", "sigma_LT", "sigma_RT",
#'                               "R_opt", "k1", "eff_SMC", "s", "b")  # Adjust as needed
#'
#' # Run simulation
#' simulations <- run_mcmc_simulations(
#'   model,
#'   param_inputs,
#'   param_samples,
#'   start_date = "2014-01-01",
#'   end_date = "2022-12-31",
#'   prewarm_years = 2,
#'   days_per_year = 360
#' )
#'
#' # View first simulation output
#' head(simulations[[1]])
#' }
#'
#' @export
run_mcmc_simulations <- function(model, param_inputs, param_samples, start_date, end_date, prewarm_years = 2, days_per_year = 360) {
  simulation_results <- lapply(1:nrow(param_samples), function(i) {
    param_sample <- as.list(param_samples[i, ])
    param_inputs <- update_param_list(param_inputs, param_sample)
    compartments_sim(model, param_inputs, start_date, end_date, prewarm_years, days_per_year)
  })

  return(simulation_results)
}

#' Summarize Simulation Outputs with Medians and Confidence Intervals
#'
#' Computes the pointwise median and confidence intervals for selected variables across
#' multiple simulation runs. Useful for summarizing posterior predictive distributions
#' of model outputs.
#'
#' @param simulation_results A list of data frames, where each data frame contains the model
#'   output for a different parameter set. Each data frame must include a \code{date_ymd} column.
#' @param ci_level Numeric value between 0 and 1 specifying the confidence level for the intervals
#'   (e.g., \code{0.95} for 95% confidence intervals).
#' @param variables A character vector indicating which compartment names to summarize
#'   (e.g., \code{c("SC", "EC", "prev_total_with_R")}). If \code{NULL}, all variables except
#'   \code{date_ymd} and \code{simulation_id} will be included.
#'
#' @return A data frame with one row per date and summary statistics (median, lower, and upper
#'   bounds) for each selected variable. Columns are named using the pattern \code{<variable>_median},
#'   \code{<variable>_lower}, and \code{<variable>_upper}.
#' @export
summarize_simulations <- function(simulation_results, ci_level = 0.95, variables = NULL) {
  # Combine all simulation runs into a single data frame
  all_results <- bind_rows(simulation_results, .id = "simulation_id")

  # If variables are NULL, use all columns except 'date' and 'simulation_id'
  if (is.null(variables)) {
    variables <- setdiff(colnames(all_results), c("date_ymd", "simulation_id"))
  } else {
    # Ensure selected variables exist in the data
    missing_vars <- setdiff(variables, colnames(all_results))
    if (length(missing_vars) > 0) {
      stop(paste("The following variables are not found in the simulation results:", paste(missing_vars, collapse = ", ")))
    }
  }

  # Group by date and calculate median and CI for selected variables
  summary_stats <- all_results %>%
    group_by(date_ymd) %>%
    summarize(across(all_of(variables), list(
      median = ~ median(.x, na.rm = TRUE),
      lower = ~ quantile(.x, probs = (1 - ci_level) / 2, na.rm = TRUE),
      upper = ~ quantile(.x, probs = 1 - (1 - ci_level) / 2, na.rm = TRUE)
    ), .names = "{.col}_{.fn}")) %>%
    ungroup()

  return(summary_stats)
}

#' Run Simulations from Sampled Parameter Sets
#'
#' This function runs a simulation model multiple times, each time using a different
#' row from a sampled parameter matrix (e.g., from posterior draws). It is used to
#' generate predictive simulations or counterfactual scenarios based on
#' parameter uncertainty. Each simulation returns output from the `data_sim()` function,
#' which integrates optional transformations and covariates.
#'
#' @param model A simulation model function (e.g., a compartmental ODE model)
#'   that takes in a parameter list and returns state trajectories.
#' @param param_inputs A named list of baseline input parameters. These serve as the
#'   default values before updating with each row from `param_samples`.
#' @param param_samples A matrix or data frame of sampled parameter values,
#'   where each row corresponds to a unique parameter set to simulate.
#' @param start_date,end_date The simulation date range as `Date` or character objects
#'   coercible to `Date`. Typically, these define the time window for prediction after prewarming.
#' @param prewarm_years Number of years before `start_date` used for prewarming the model
#'   (e.g., to reach equilibrium or load past covariates). Defaults to 2.
#' @param mu_transform_C Optional function to transform the underlying `mu`
#'   incidence into a predicted count (e.g., logistic or log link for cluster C).
#' @param mu_transform_A Optional transformation function for cluster A (if used).
#' @param covariate_matrix Optional covariate matrix used for transforming incidence
#'   (e.g., for climate-driven models or spatial heterogeneity).
#' @param month Logical. If `TRUE`, the output will be aggregated or indexed monthly
#'   (rather than by week or day).
#' @param noise Logical. If `TRUE`, draws the simulated incidence from a negative
#'   binomial distribution using `size_1` as the dispersion parameter.
#'
#' @return A list of data frames, each containing a simulation run with the
#'   same structure as the output of `data_sim()`. The list has one element per row
#'   in `param_samples`.
#'
#' @examples
#' \dontrun{
#' sim_list <- run_simulations_from_samples(
#'   model = malaria_model,
#'   param_inputs = base_params,
#'   param_samples = posterior_draws,
#'   start_date = "2016-01-01",
#'   end_date = "2022-12-31",
#'   mu_transform_C = logit_transform,
#'   covariate_matrix = climate_covariates,
#'   noise = TRUE
#' )
#' }
#'
#' @export
run_simulations_from_samples <- function(model, param_inputs, param_samples,
                                         start_date, end_date,
                                         prewarm_years = 2,
                                         mu_transform_C = NULL,
                                         mu_transform_A = NULL,
                                         covariate_matrix = NULL,
                                         month = TRUE,
                                         noise = FALSE) {
  lapply(1:nrow(param_samples), function(i) {
    updated_inputs <- update_param_list(param_inputs, as.list(param_samples[i, ]))
    size = updated_inputs$size_1
    #compartments_sim(model, updated_inputs, start_date, end_date, prewarm_years, days_per_year)
    data_sim(model, updated_inputs, as.Date(start_date), as.Date(end_date), prewarm_years = 2, save = FALSE,
             mu_transform_C = mu_transform_C, mu_transform_A = mu_transform_A, month = month,
             covariate_matrix = covariate_matrix, noise = noise, size = size)
  })
}

#' Create Long Format Data for Simulated Posterior Median and CI
#'
#' @param simulations List of simulation data frames.
#' @param variables Character vector of variable names.
#' @param ci_level Confidence interval width (e.g., 0.95).
#' @return Data frame in long format.
#' @export
summarize_simulation_ci <- function(simulations, variables, ci_level = 0.95) {
  summarized <- summarize_simulations(simulations, ci_level = ci_level, variables = variables)
  summarized %>%
    pivot_longer(-date_ymd, names_to = "var_stat", values_to = "value") %>%
    separate(var_stat, into = c("variable", "stat"), sep = "_(?=[^_]+$)") %>%
    pivot_wider(names_from = stat, values_from = value)
}

#' Run the model once at the maximum-posterior parameter set (deterministic)
#'
#' @param model        A dust or odin model object.
#' @param param_inputs Named list of baseline inputs.
#' @param best_params  Named vector or list of best-fit parameters.
#' @param start_date   Date or string start of sim window.
#' @param end_date     Date or string end of sim window.
#' @param ...          Additional args passed to \code{data_sim()}.
#' @return A data.frame in long form with \code{date_ymd}, \code{variable}, \code{value}.
#' @export
simulate_best <- function(model, param_inputs, best_params,
                          start_date, end_date, ...) {
  inputs <- update_param_list(param_inputs, best_params)
  df_wide <- data_sim(
    model        = model,
    param_inputs = inputs,
    start_date   = start_date,
    end_date     = end_date,
    noise        = FALSE,
    ...
  )
  tidyr::pivot_longer(
    df_wide,
    cols      = -date_ymd,
    names_to  = "variable",
    values_to = "value"
  )
}

