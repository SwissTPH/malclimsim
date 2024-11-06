

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

# Function to calculate the decay of SMC efficacy over time for an entire schedule.
# This function computes the decay between each SMC treatment, based on the efficacy decay curve.
# Args:
#   - SMC: A numeric vector where each element represents an SMC event (1 = SMC administered, 0 = no SMC).
#   - decay_func: The decay function used to calculate efficacy (default = decay_SMC).
#   - const: Decay constant to be passed to the decay function (default = -0.1806).
# Returns:
#   - A numeric vector representing the decayed efficacy over time.
calc_decay_arr <- function(SMC, decay_func = decay_SMC, const = -0.1806){
  decay_arr <- array(NA, dim = length(SMC))  # Initialize decay array

  # Loop through each SMC administration point
  for (i in 1:length(SMC)) {
    if (SMC[i] > 0) {  # Check if SMC was administered
      next_SMC <- i + which(SMC[(i + 1):length(SMC)] > 0)  # Find the next SMC administration

      # If no subsequent SMC, calculate decay until the end
      if (length(next_SMC) == 0) {
        end <- length(SMC)
        decay_arr[i:end] <- decay_func(seq(0, (end - i)), const = const)
      } else {  # If there is a subsequent SMC, calculate decay until the next SMC
        end <- next_SMC[1] - 1
        decay_arr[i:end] <- decay_func(seq(0, (end - i)), const = const)
      }
    }
  }

  # Set any missing values in the decay array to 0 (no decay if no SMC was administered)
  decay_arr[is.na(decay_arr)] <- 0
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

# Function to simulate malaria cases over time, either by month or week
# Arguments:
#   model_SMC: the model for SMC interventions
#   param_inputs: model parameters (including SMC schedule and efficacy)
#   start_date, end_date: date range for the simulation
#   month: whether to aggregate results by month (TRUE) or by week (FALSE)
#   round: whether to round the results (default is TRUE)
#   save: whether to save the results to a file (default is TRUE)
#   file: file path to save results if save = TRUE
#   month_unequal_days: account for unequal days in months (default is FALSE)
data_sim <- function(model_SMC, param_inputs, start_date, end_date,
                     month = FALSE, round = TRUE, save = TRUE, file = "",
                     month_unequal_days = FALSE){
  # Calculate the number of days in the SMC schedule
  n_days <- length(param_inputs$SMC)

  # Run the simulation using the model and parameters
  results <- sim_mod(model_SMC, pars = c(param_inputs), time_start = 0,
                     n_particles = 1, sim_time = n_days)

  # Extract the simulation output and the model
  x <- results[[1]]
  mod <- results[[2]]

  # If monthly aggregation is selected
  if(month){
    # Create monthly indices for aggregation
    month_ind <- seq(1, length(param_inputs$SMC), by = 30)

    # Adjust for months with unequal days, if specified
    if(month_unequal_days){
      month_ind <- which(param_inputs$day_count == 0)
    }

    # Extract monthly incidence for adults (inc_A) and children (inc_C)
    inc_A <- x[mod$info()$index$month_inc_A,,][month_ind]
    inc_C <- x[mod$info()$index$month_inc_C,,][month_ind]

    # Get the corresponding monthly dates
    month <- date_to_months(start_date = as.Date(start_date), end_date = as.Date(end_date))

    # Create a dataframe for monthly incidence data
    month_no <- 0:(length(inc_C) - 1)
    inc_df <- data.frame(date_ymd = month, month_no, inc_A, inc_C, inc = inc_A + inc_C)

    # If weekly aggregation is selected
  } else {
    # Create weekly indices for aggregation
    wk_ind <- seq(1, n_days, by = 7)

    # Extract weekly incidence for adults (inc_A) and children (inc_C)
    inc_A <- x[mod$info()$index$wk_inc_A,,][wk_ind]
    inc_C <- x[mod$info()$index$wk_inc_C,,][wk_ind]

    # Get the corresponding weekly dates
    week <- date_to_weeks(start_date = as.Date(start_date), end_date = as.Date(end_date))

    # Create a dataframe for weekly incidence data
    week_no <- 0:(length(inc_C) - 1)
    inc_df <- data.frame(week, week_no, inc_A, inc_C, inc = inc_A + inc_C)
  }

  # Round the results if specified
  if(round){
    inc_df[3:5] <- round(inc_df[3:5])
  }

  # Save the dataframe to a file if specified
  if(save){
    saveRDS(inc_df, paste(dir, file, sep = ""))
  }

  # Return the resulting incidence dataframe
  return(inc_df)
}

